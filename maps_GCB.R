



#######################################
# maps
rm(list=ls())
gc()
gc()
"%ni%" <- Negate("%in%")
mse <- function(x, y){mean((x-y)^2)}

library(devtools)
install_github("cranedroesch/panelNNET", ref = "earlystop")
library(panelNNET)
library(doParallel)
library(doBy)
library(glmnet)
library(dplyr)
library(nnls)
library(ggplot2)
library(broom)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(tidyverse)
library(scales)
library(mgcv)

AWS <- grepl('ubuntu', getwd())
desktop <- grepl(':', getwd())
laptop <- grepl('/home/andrew', getwd())
if(AWS){
  setwd("/home/ubuntu/projdir")
  system("sudo mkdir /home/ubuntu/projdir/outdir")
  outdir <- "/home/ubuntu/projdir/outdir"
  registerDoParallel(detectCores())
}
if(desktop){
}
if(laptop){
  setwd("/home/andrew/Dropbox/USDA/ARC/data")
  outdir <- "/home/andrew/Dropbox/USDA/ARC/output"
  registerDoParallel(detectCores())
}
dat <- readRDS("panel_corn.Rds")
dat <- subset(dat, state %in% c("17", "19", "27", "18", "39", "26", "21", "55", "29"))
# irr/dry
if (irrdry == "irr"){
  dat <- dat[dat$prop_irr > .5,]
  registerDoParallel(detectCores())
} else {
  dat <- dat[dat$prop_irr < .5,]
  # registerDoParallel(64)
}

# make TAP variable
dat$TAP <- rowSums(dat[,grepl('precip', colnames(dat))])/1000
dat$TAP2 <- dat$TAP^2
# nonparametric
X <- dat[,grepl('year|SR|soil_|precip|sdsfia|wind|minat|maxat|minrh|maxrh|lat|lon|prop_irr',colnames(dat))]
X <- X[sapply(X, sd) > 0]
# parametric
dat$y <- dat$year - min(dat$year) + 1
dat$y2 <- dat$y^2

# Xp <- cbind(dat[,c("y", "y2", "TAP", "TAP2")])
Xp <- cbind(dat[,c("y", "y2", "TAP", "TAP2")], dat[,grepl("SR", colnames(dat))])
Xp <- Xp[sapply(Xp, sd) > 0]

# OLS baseline
mm <- model.matrix(as.formula(paste0("~yield+y+y2+TAP+TAP2+",
                                     paste(colnames(dat)[grepl("SR", colnames(dat))], collapse = "+"),
                                     "-1")),
                   data = dat)



# maps
par(mfrow = c(1,1))
gispath <- paste0(getwd(),"/shapefiles")

USA_shp <- readOGR(dsn = gispath, layer = "statesaea")
head(as.data.frame(USA_shp))
USA_shp <- subset(USA_shp, STATE_NAME %in% c('Illinois','Iowa', "Minnesota", "Michigan", "Missouri", "Kentucky", "Wisconsin", "Indiana", "Ohio"))
USA_shp <- spTransform(USA_shp, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))

#Load county shapefile
fips_shp <- readOGR(dsn = gispath, layer = "FIPS2010")
fips_shp <- subset(fips_shp, STATE_NAME %in% c('Illinois','Iowa', "Minnesota", "Michigan", "Missouri", "Kentucky", "Wisconsin", "Indiana", "Ohio"))
fips_shp <- spTransform(fips_shp, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "))

theme_map <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.ticks.length=unit(0.3, "lines"),
      axis.ticks.margin=unit(0.5, "lines"),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.background=element_rect(fill="white", colour=NA),
      legend.key=element_rect(colour="white"),
      legend.key.size=unit(1.2, "lines"),
      legend.position="right",
      legend.text=element_text(size=rel(0.8)),
      legend.title=element_text(size=rel(0.8), face="bold", hjust=0),
      panel.background=element_blank(),
      #      panel.background = element_rect(fill = NULL, colour = "white"),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.spacing=unit(0, "lines"),
      plot.background=element_blank(),
      plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"),
      plot.title=element_text(size=rel(1.2), vjust = 2),
      strip.background=element_rect(fill="white", colour="white"),
      strip.text.x=element_text(size=rel(.8)),
      strip.text.y=element_text(size=rel(0.8), angle=-90) 
    )   
}

USA_fort_premade <- tidy(USA_shp, region = "STATE_FIPS")

############################
# maps of predictive skill
PNN_oob <- load_obj(paste0(outdir,"/pnnoob_GCB.Rda"))
ens_weights<- readRDS(file = paste0(outdir,"/ensemble_weights_GCB.Rds"))
hist_pred <- data.frame(PNN = rowMeans(PNN_oob, na.rm=T),
                        PNN_weighted = apply(PNN_oob, 1, weighted.mean, w = ens_weights, na.rm=T)
)
hist_pred <- cbind(dat[with(dat, order(fips, year)),c("yield", "fips", "year")], hist_pred)

#SR
dis <- as.data.frame(demeanlist(mm, list(dat$fips)))
m <- glmnet(y = dis$yield, x = as.matrix(dis[,-1]), lambda = 0, intercept = FALSE)
XB <- mm[,-1] %*% coef(m)[-1,]
fe <- (dat$yield-dis$yield) - (mm[,-1]-as.matrix(dis[,-1])) %*% coef(m)[-1,]
tm <- data.frame(fips = dat$fips, fe = fe)
tm <- tm[!duplicated(tm),]
tm <- summaryBy(fe~fips, data = tm, keep.names = T)
om <- data.frame(XB, fips = dat$fips, yield = dat$yield, year = dat$year)
pred <- merge(om, tm)
pred$OLS <- with(pred, fe+XB)

#bagged SR
nfolds <- 96
mergetemplate <- dat[,c("fips", "year")]
mergetemplate <- mergetemplate[with(mergetemplate, order(fips, year)),]
SRoob <- foreach(i = 1:nfolds, .combine = cbind) %dopar% {
  set.seed(i, kind = "L'Ecuyer-CMRG")
  samp <- sample(unique(dat$year), replace = TRUE)
  bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
  oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
  mmis <- mm[bsamp,]
  dis <- as.data.frame(demeanlist(mmis, list(dat$fips[bsamp])))
  m <- glmnet(y = dis$yield, x = as.matrix(dis[,-1]), lambda = 0, intercept = FALSE)
  XB <- mm[dat$year %in% oosamp,-1] %*% coef(m)[-1,]
  fe <- (dat$yield[bsamp]-dis$yield) - (mmis[,-1]-as.matrix(dis[,-1])) %*% coef(m)[-1,]
  tm <- data.frame(fips = dat$fips[bsamp], fe = fe)
  tm <- tm[!duplicated(tm),]
  tm <- summaryBy(fe~fips, data = tm, keep.names = TRUE)
  om <- data.frame(XB, year = dat$year[dat$year %in% oosamp], fips = dat$fips[dat$year %in% oosamp], yield = dat$yield[dat$year %in% oosamp])
  pred <- merge(om, tm)
  pred$pred <- with(pred, fe+XB)
  with_holes <- full_join(mergetemplate, pred)
  with_holes <- with_holes[with(with_holes, order(fips, year)),]
  return(with_holes$pred)
}
tosweep <- dat$yield[order(dat$fips, dat$year)]
SR_unbagged <- sweep(SRoob, 1, tosweep, "-")^2 %>% rowMeans(na.rm=T) %>% sqrt

SRoob <- rowMeans(SRoob, na.rm = T)
mergetemplate$OLS_bagged <- SRoob
hist_pred <- full_join(hist_pred, mergetemplate)
hist_pred <- hist_pred[with(hist_pred, order(fips, year)),]

err_frame <- cbind(hist_pred[,c("fips", "year")], 
                   sweep(hist_pred[,4:6], 1, hist_pred$yield, "-") %>% abs,
                   SR_unbagged)

# gather
predframe <- gather(err_frame, type, err, -fips, -year)
predframe$id <- predframe$fips
predframe$err <- predframe$err^2
sumpredframe <- summaryBy(err~fips+type, data = predframe, keep.names = T)
sumpredframe$err <- sumpredframe$err^.5
sumpredframe$id <- sumpredframe$fips

# fortify
fips_shp <- gBuffer(fips_shp, byid=TRUE, width=0)
fips_fort <- broom::tidy(fips_shp, region = "FIPS")
fips_fort <- inner_join(fips_fort, sumpredframe)
fips_fort$type[fips_fort$type == "SR_unbagged"] <- "OLS"
fips_fort$type[fips_fort$type == "OLS_bagged"] <- "Bagged OLS"
fips_fort$type[fips_fort$type == "PNN"] <- "PNN"
fips_fort$type[fips_fort$type == "PNN_weighted"] <- "Weighted PNN"
fips_fort$type <- factor(fips_fort$type, levels = c("OLS", "Bagged OLS", "PNN", "Weighted PNN"))

mapres <- ggplot(USA_fort_premade) + 
  aes(long,lat,group=group) + 
  geom_polygon(aes(x=long,y=lat, group=group, fill=err), 
               data=fips_fort)+
  geom_path(color="black", lwd = .1) +
  coord_fixed(1.3) +
  facet_wrap( ~ type, nrow = 1) +
  scale_fill_gradient2(low = "darkgreen", mid = "yellow",
                       high = "red", midpoint = mean(fips_fort$err, na.rm = TRUE), 
                       space = "Lab", name = "E|error|",
                       na.value = "grey50", guide = "colourbar") +
  ggtitle("Expected prediction Error")+theme_map()
mapres
dev.copy2pdf(file = "/home/andrew/Dropbox/USDA/yield_forecasting/tex/map_GCB_RMSE.pdf", height = 3, width = 8)
system("evince /home/andrew/Dropbox/USDA/yield_forecasting/tex/map_GCB_RMSE.pdf")

# diff pred map
OLS <- sumpredframe[sumpredframe$type == "SR_unbagged",]
sumpredframe <- sumpredframe[sumpredframe$type %ni% c("SR_unbagged", "PNN_weighted"),]
OLS$OLSerr <- OLS$err
dif_frame <- inner_join(sumpredframe, OLS[,c("fips", "OLSerr")])
dif_frame$dif <- with(dif_frame, OLSerr - err)

fips_fort <- broom::tidy(fips_shp, region = "FIPS")
fips_fort <- inner_join(fips_fort, dif_frame)
fips_fort$type[fips_fort$type == "OLS_bagged"] <- "Bagged OLS"
fips_fort$type[fips_fort$type == "PNN"] <- "SNN"
fips_fort$type[fips_fort$type == "PNN_weighted"] <- "Weighted PNN"
fips_fort$type <- factor(fips_fort$type, levels = c("Bagged OLS", "SNN", "Weighted PNN"))

mapres <- ggplot(USA_fort_premade) +
  aes(long,lat,group=group) +
  geom_polygon(aes(x=long,y=lat, group=group, fill=dif),
               data=fips_fort)+
  geom_path(color="black", lwd = .1) +
  coord_fixed(1.3) +
  facet_wrap( ~ type, nrow = 1) +
  scale_fill_gradientn(colors = c("red", "lightgrey", "light blue", "darkgreen"),
                       values = rescale(c(min(fips_fort$dif, na.rm=T), 
                                          0, 
                                          mean(dif_frame$dif[dif_frame$type == "OLS_bagged"]), 
                                          max(fips_fort$dif, na.rm=T))),
                       # values = rescale(c(-100, 0, 25)),
                       limits=range(fips_fort$dif, na.rm=T),
                       space = "Lab", name = str_wrap("difference", width = 8),
                       na.value = "grey50", guide = "colourbar")+
  ggtitle("Prediction improvement over OLS")+theme_map()
mapres
dev.copy2pdf(file = "/home/andrew/Dropbox/USDA/yield_forecasting/tex/map_GCB_RMSE_diff.pdf", height = 3, width = 4)
system("evince /home/andrew/Dropbox/USDA/yield_forecasting/tex/map_GCB_RMSE_diff.pdf")


###################
# Projections
predfiles <- list.files(pattern = "_predframe", path = "/home/andrew/Dropbox/USDA/ARC/output/", full.names = T)
predfiles <- predfiles[!grepl("v2", predfiles)]
for (i in 1:13){
  predframe <- load_obj(predfiles[i])
  predframe <- gather(predframe, type, pred, -fips, -year)
  tag <- ((predfiles[i] %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[1]
  # historical data
  dtreg <- gam(yield~s(year, k=4), method = "REML", data = dat)
  summary(dtreg)
  dat$yield_dt <- predict(dtreg, newdata = data.frame(year = 2016)) + residuals(dtreg)
  sumHist <- summaryBy(yield_dt~fips, data = dat, keep.names = T, FUN = function(x){c(mean(x), sd(x))})
  colnames(sumHist)[2:3] <- c("yield_hist", "sd_yield_hist")
  
  # summarize by periods
  sumMid <- summaryBy(pred~fips+type, data = filter(predframe, year > 2040 & year<2070), keep.names = T, FUN = function(x){c(mean(x), sd(x))})
  sumLate <- summaryBy(pred~fips+type, data = filter(predframe, year > 2070), keep.names = T, FUN = function(x){c(mean(x), sd(x))})
  colnames(sumLate)[3:4] <- colnames(sumMid)[3:4] <- c("yield", "sd_yield")
  sumMid$period <- "2040-70"
  sumLate$period <- "2070-99"
  sumMid$RCP <- ifelse(grepl("45", sumMid$type), "RCP4.5", "RCP8.5")
  sumMid$type <- gsub("45|85", "", sumMid$type)
  sumLate$RCP <- ifelse(grepl("45", sumLate$type), "RCP4.5", "RCP8.5")
  sumLate$type <- gsub("45|85", "", sumLate$type)
  
  scen <- rbind(sumMid, sumLate)
  tomap <- full_join(scen, sumHist)
  tomap$id <- tomap$fips
  tomap$pctdiff <- (tomap$yield/tomap$yield_hist -1) *100
  tomap <- filter(tomap, type %in% c("OLS", "PNN"))
  
  fips_fort <- broom::tidy(fips_shp, region = "FIPS")
  fips_fort <- inner_join(fips_fort, tomap)
  
  toadd <- data.frame(fips = NA, type = NA, yield = NA, sd_yield = NA, period = NA, RCP = NA, yield_hist = NA, sd_yield_hist = NA,
                      id= "999", pctdiff = c(-100, 25))
  tomap <- rbind(tomap, toadd)
  library(stringr)
  # percent difference
  mapPct <- ggplot(USA_fort_premade) +
    aes(long,lat,group=group) +
    geom_polygon(aes(x=long,y=lat, group=group, fill=pctdiff),
                 data=fips_fort)+
    geom_path(color="black", lwd = .1) +
    coord_fixed(1.3) +
    facet_wrap( ~ RCP+period+type, ncol = 4,labeller = function (labels) {
      labels <- lapply(labels, as.character)
      list(do.call(paste, c(labels, list(sep = ", "))))
    }) +
    scale_fill_gradientn(colors = c("brown", "lightgrey", "darkgreen"),
                         values = rescale(c(min(tomap$pctdiff), 0, max(tomap$pctdiff))),
                         # values = rescale(c(-100, 0, 25)),
                         limits=range(tomap$pctdiff),
                         space = "Lab", name = str_wrap("percent change", width = 8),
                         na.value = "grey50", guide = "colourbar")+
    theme_map() +ggtitle(paste0("Expected yield change: ", tag))
  print(mapPct)
  dev.copy2pdf(file = paste0("/home/andrew/Dropbox/USDA/yield_forecasting/tex/map_yieldchange_", tag,".pdf"), height = 6, width = 8)
  
  OLS <- tomap[tomap$type == "OLS",]
  PNN <- tomap[tomap$type == "PNN",]
  colnames(OLS)[colnames(OLS) == "pctdiff"] <- "OLS"
  colnames(PNN)[colnames(PNN) == "pctdiff"] <- "PNN"
  tomap2 <- (full_join(OLS[,c("id", "OLS","period","RCP")], PNN[,c("id", "PNN","period","RCP")]))
  tomap2$diff <- with(tomap2, PNN - OLS)
  toadd <- data.frame(id= "999", OLS = NA, period = NA, RCP = NA, PNN = NA, diff = c(-10, 25))
  tomap2 <- rbind(tomap2, toadd)
  fips_fort <- broom::tidy(fips_shp, region = "FIPS")
  fips_fort <- inner_join(fips_fort, tomap2)
  
  mapdiff <- ggplot(USA_fort_premade) +
    aes(long,lat,group=group) +
    geom_polygon(aes(x=long,y=lat, group=group, fill=diff),
                 data=fips_fort)+
    geom_path(color="black", lwd = .1) +
    coord_fixed(1.3) +
    facet_wrap( ~ RCP+period, ncol = 4,labeller = function (labels) {
      labels <- lapply(labels, as.character)
      list(do.call(paste, c(labels, list(sep = ", "))))
    }) +
    scale_fill_gradient2(low="orange", mid="grey80", high="blue", #colors in the scale
                         midpoint=0,    
                         limits=range(tomap2$diff, na.rm=T),
                         space = "Lab", name = "yield (bu/ac)",
                         na.value = "grey50", guide = "colourbar")+
    theme_map() +ggtitle(paste0("Difference between SNN and OLS: ", tag,"\n"))
  print(mapdiff)
  dev.copy2pdf(file = paste0("/home/andrew/Dropbox/USDA/yield_forecasting/tex/map_yielddiff_", tag,".pdf"), height = 6, width = 8)
}
