
#######################
# predict
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
library(mgcv)
library(tidyverse)


AWS <- grepl('ubuntu', getwd())
desktop <- grepl(':', getwd())
laptop <- grepl('/home/andrew', getwd())
if(AWS){
  setwd("/home/ubuntu/projdir/")
  outdir <- "/home/ubuntu/projdir/outdir"
  registerDoParallel(32) # too big, even on R4xlarge!
}
if(desktop){
}
if(laptop){
  setwd("/home/andrew/Dropbox/USDA/ARC/data")
  outdir <- "/home/andrew/Dropbox/USDA/ARC/output"
  registerDoParallel(detectCores())
}
cr <- "corn"
if (cr == "corn"){dat <- readRDS("panel_corn.Rds")}
dat <- subset(dat, state %in% c("17", "19", "27", "18", "39", "26", "21", "55", "29"))
dat <- dat[dat$prop_irr < .5,]
dat$TAP <- rowSums(dat[,grepl('precip', colnames(dat))])/1000
dat$TAP2 <- dat$TAP^2
X <- dat[,grepl('year|SR|soil_|precip|sdsfia|wind|minat|maxat|minrh|maxrh|lat|lon|prop_irr',colnames(dat))]
X <- X[sapply(X, sd) > 0]
dat$y <- dat$year - min(dat$year) + 1
dat$y2 <- dat$y^2
Xp <- cbind(dat[,c("y", "y2", "TAP", "TAP2")], dat[,grepl("SR", colnames(dat))])
Xp <- Xp[sapply(Xp, sd) > 0]

mm <- model.matrix(as.formula(paste0("~yield+y+y2+TAP+TAP2+",
                                     paste(colnames(dat)[grepl("SR", colnames(dat))], collapse = "+"),
                                     "-1")), data = dat)
dis <- as.data.frame(demeanlist(mm, list(dat$fips)))
m <- glmnet(y = dis$yield, x = as.matrix(dis[,-1]), lambda = 0, intercept = FALSE)
fe <- (dat$yield-dis$yield) - (mm[,-1]-as.matrix(dis[,-1])) %*% coef(m)[-1,]
XB <- mm[,-1] %*% coef(m)[-1,]
tm <- data.frame(fips = dat$fips, fe = fe)
tm <- tm[!duplicated(tm),]
tm <- summaryBy(fe~fips, keep.names = T, data = tm)

fi <- list.files(path = outdir, pattern = "v6", full.names = TRUE)


fi <- fi[grepl("v6", fi)]
foldkey <- fi %>% strsplit("dry") %>% sapply(function(x){unlist(x)[2]}) %>% 
  gsub(pattern = "_", replacement = "") %>% 
  gsub(pattern = ".Rda", replacement = "") %>% 
  as.numeric
if(AWS){
  scenario_panels <- list.files(path = "/home/ubuntu/MACA/panels", pattern = "panel_GCB_rcp")  
}
if(laptop){
  scenario_panels <- list.files(path = outdir, pattern = "panel_GCB_rcp")  
}
scenario_panels <- scenario_panels[!grepl("corn.Rds", scenario_panels)] #lose the original hadley
mod <- scenario_panels %>% strsplit("corn") %>% sapply(function(x){unlist(x)[2]}) %>% unique

# faster predict function for PCA
ppc <- function(x, PC, ncomp){
  x <- x %>% sweep(2, PC$center, "-")%>% sweep(2, PC$scale, "/")
  MatMult(as.matrix(x), PC$rotation[,1:ncomp])
}

nfolds <- 96
# loop through scenarios and make predictions
for(M in mod){
  tag <- gsub(".Rds","", M)
  # load data
  if(AWS){
    scfi <- list.files(pattern = M, path = "/home/ubuntu/MACA/panels", full.names = T)
  }
  if(laptop){
    scfi <- list.files(pattern = M, path = outdir, full.names = T)
  }
  dat45 <- readRDS(scfi[grepl("rcp45", scfi)])
  dat85 <- readRDS(scfi[grepl("rcp85", scfi)])
  dat45$prop_irr <- dat85$prop_irr <- 0

  # fix issue with NANs in the SR variables having the odd NAN in them.
  havenan <- dat45[,grepl("SR", colnames(dat45))] %>% rowSums %>% is.nan %>% which
  for (i in havenan){
    x <- dat45[i, grep("SR", colnames(dat45))] %>% as.numeric
    g <- 1:length(x)
    mnan <- gam(x~s(g), method = "REML", family = tw)
    pfill <- predict(mnan, newdata = data.frame(g = g[is.nan(x)]))
    dat45[i, which(sapply(dat45[i,], is.nan))] <- pfill
  }
  havenan <- dat85[,grepl("SR", colnames(dat85))] %>% rowSums %>% is.nan %>% which
  for (i in havenan){
    x <- dat85[i, grep("SR", colnames(dat85))] %>% as.numeric
    g <- 1:length(x)
    mnan <- gam(x~s(g), method = "REML", family = tw)
    pfill <- predict(mnan, newdata = data.frame(g = g[is.nan(x)]))
    dat85[i, which(sapply(dat85[i,], is.nan))] <- pfill
  }

  # make TAP variable
  dat45$TAP <- rowSums(dat45[,grepl('precip', colnames(dat45))])/1000
  dat45$TAP2 <- dat45$TAP^2
  dat85$TAP <- rowSums(dat85[,grepl('precip', colnames(dat85))])/1000
  dat85$TAP2 <- dat85$TAP^2
  # nonparametric
  X45 <- dat45[,grepl('year|SR|soil_|precip|sdsfia|wind|minat|maxat|minrh|maxrh|lat|lon|prop_irr',colnames(dat45))]
  X85 <- dat85[,grepl('year|SR|soil_|precip|sdsfia|wind|minat|maxat|minrh|maxrh|lat|lon|prop_irr',colnames(dat85))]
  X45 <- X45[,colnames(X)]
  X85 <- X85[,colnames(X)]

  # parametric
  dat45$y <- dat85$y <- max(dat$y)
  dat45$y2 <- dat85$y2 <- max(dat$y2)
  Xp45 <- cbind(dat45[,c("y", "y2", "TAP", "TAP2")], dat45[,grepl("SR", colnames(dat45))])
  Xp85 <- cbind(dat85[,c("y", "y2", "TAP", "TAP2")], dat85[,grepl("SR", colnames(dat85))])
  Xp45 <- Xp45[,colnames(Xp)]
  Xp85 <- Xp85[,colnames(Xp)]

  # SR pred
  mm45 <- model.matrix(as.formula(paste0("~y+y2+TAP+TAP2+",
                                         paste(colnames(dat)[grepl("SR", colnames(dat))], collapse = "+"),
                                         "-1")), data = dat45)
  mm85 <- model.matrix(as.formula(paste0("~y+y2+TAP+TAP2+",
                                         paste(colnames(dat)[grepl("SR", colnames(dat))], collapse = "+"),
                                         "-1")), data = dat85)


  XB45 <- mm45 %*% coef(m)[-1,]
  XB85 <- mm85 %*% coef(m)[-1,]
  om <- data.frame(XB45, XB85, fips = dat45$fips, year = dat45$year)
  pred <- merge(om, tm)
  pred$pred45 <- with(pred, XB45 + fe)
  pred$pred85 <- with(pred, XB85 + fe)
  # truncate
  pred$pred85[pred$pred85<0] <- 0
  pred$pred45[pred$pred45<0] <- 0
  pred$pred85[pred$pred85>200] <- 200
  pred$pred45[pred$pred45>200] <- 200

  # generate predframe
  predframe <- data.frame(fips = pred$fips, year = pred$year,
                          OLS45 = pred$pred45,
                          OLS85 = pred$pred85)

  # SR bagged
  mergetemplate <- predframe[,c("fips", "year")]
  SRoob_bagged_45 <- foreach(i = 1:nfolds, .combine = cbind) %dopar% {
    set.seed(i, kind = "L'Ecuyer-CMRG")
    samp <- sample(unique(dat$year), replace = TRUE)
    bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
    oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
    mmis <- mm[bsamp,]
    dis <- as.data.frame(demeanlist(mmis, list(dat$fips[bsamp])))
    m <- glmnet(y = dis$yield, x = as.matrix(dis[,-1]), lambda = 0, intercept = FALSE)
    XB <- mm45 %*% coef(m)[-1,]
    fe <- (dat$yield[bsamp]-dis$yield) - (mmis[,-1]-as.matrix(dis[,-1])) %*% coef(m)[-1,]
    tm <- data.frame(fips = dat$fips[bsamp], fe = fe)
    tm <- tm[!duplicated(tm),]
    tm <- summaryBy(fe~fips, data = tm, keep.names = TRUE)
    om <- data.frame(XB, year = dat45$year, fips = dat45$fips)
    pred <- merge(om, tm)
    pred <- full_join(pred, mergetemplate, by = c("fips", "year"))
    np <- pred$XB+pred$fe
    return(np)
  }
  # Truncate
  SRoob_bagged_45[SRoob_bagged_45 < 0] <- 0
  SRoob_bagged_45[SRoob_bagged_45 > 200] <- 200
  predframe$OLS_bagged45 <- rowMeans(SRoob_bagged_45)
  saveRDS(SRoob_bagged_45, file = paste0(outdir, "/",tag, "_SR_bagframe_45.Rds"))


  SRoob_bagged_85 <- foreach(i = 1:nfolds, .combine = cbind) %dopar% {
    set.seed(i, kind = "L'Ecuyer-CMRG")
    samp <- sample(unique(dat$year), replace = TRUE)
    bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
    oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
    mmis <- mm[bsamp,]
    dis <- as.data.frame(demeanlist(mmis, list(dat$fips[bsamp])))
    m <- glmnet(y = dis$yield, x = as.matrix(dis[,-1]), lambda = 0, intercept = FALSE)
    XB <- mm85 %*% coef(m)[-1,]
    fe <- (dat$yield[bsamp]-dis$yield) - (mmis[,-1]-as.matrix(dis[,-1])) %*% coef(m)[-1,]
    tm <- data.frame(fips = dat$fips[bsamp], fe = fe)
    tm <- tm[!duplicated(tm),]
    tm <- summaryBy(fe~fips, data = tm, keep.names = TRUE)
    om <- data.frame(XB, year = dat85$year, fips = dat85$fips)
    pred <- merge(om, tm)
    pred <- full_join(pred, mergetemplate, by = c("fips", "year"))
    np <- pred$XB+pred$fe
    return(np)
  }
  # Truncate
  SRoob_bagged_85[SRoob_bagged_85 < 0] <- 0
  SRoob_bagged_85[SRoob_bagged_85 > 200] <- 200
  predframe$OLS_bagged85 <- rowMeans(SRoob_bagged_85)
  saveRDS(SRoob_bagged_85, file = paste0(outdir, "/",tag, "_SR_bagframe_85.Rds"))

  # pnn pred
  pnnpred <- foreach(i = 1:length(fi), .combine = cbind) %dopar% {
    print(i)
    set.seed(foldkey[i], kind = "L'Ecuyer-CMRG")
    samp <- sample(unique(dat$year), replace = TRUE)
    bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
    oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
    obj <- load_obj(fi[i])
    PC <- readRDS(paste0(outdir,"/PCA_",foldkey[i], ".Rds"))
    Xpc <- ppc(X[bsamp,], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95)) %>% as.data.frame
    X45pc <- ppc(X45, PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95)) %>% as.data.frame
    parlist <- obj$parlist
    pnn <- panelNNET(y = dat$yield[bsamp],
                     X = Xpc,
                     hidden_units = rep(100, 10),
                     fe_var = dat$fips[bsamp],
                     maxit = 0,
                     time_var = dat$year[bsamp],
                     param = Xp[bsamp,],
                     verbose = T,
                     activation = 'lrelu',
                     parlist = parlist
                     )
    newfips <- dat45$fips[dat45$fips %in% dat$fips[dat$year %in% samp]]
    newyear <- dat45$year[dat45$fips %in% dat$fips[dat$year %in% samp]]
    predoos <- predict(pnn, newX = as.matrix(X45pc[dat45$fips %in% dat$fips[dat$year %in% samp],]),
                       new.param = as.matrix(Xp45[dat45$fips %in% dat$fips[dat$year %in% samp],]),
                       fe.newX = newfips)
    tm <- data.frame(fips = newfips, year = newyear, pred = predoos)
    with_holes <- full_join(mergetemplate, tm)
    return(with_holes$pred)
  }

  pnnpred[pnnpred < 0] <- 0
  pnnpred[pnnpred > 200] <- 200
  predframe$PNN45 <- rowMeans(pnnpred)
  saveRDS(pnnpred, file = paste0(outdir, "/",tag, "_PNN_bagframe_45.Rds"))

  pnnpred <- foreach(i = 1:length(fi), .combine = cbind) %dopar% {
    print(i)
    set.seed(foldkey[i], kind = "L'Ecuyer-CMRG")
    samp <- sample(unique(dat$year), replace = TRUE)
    bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
    oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
    obj <- load_obj(fi[i])
    PC <- readRDS(paste0(outdir,"/PCA_",foldkey[i], ".Rds"))
    Xpc <- ppc(X[bsamp,], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95)) %>% as.data.frame
    X85pc <- ppc(X85, PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95)) %>% as.data.frame
    parlist <- obj$parlist
    pnn <- panelNNET(y = dat$yield[bsamp],
                     X = Xpc,
                     hidden_units = rep(100, 10),
                     fe_var = dat$fips[bsamp],
                     maxit = 0,
                     time_var = dat$year[bsamp],
                     param = Xp[bsamp,],
                     verbose = T,
                     activation = 'lrelu',
                     parlist = parlist
    )
    newfips <- dat85$fips[dat85$fips %in% dat$fips[dat$year %in% samp]]
    newyear <- dat85$year[dat85$fips %in% dat$fips[dat$year %in% samp]]
    predoos <- predict(pnn, newX = as.matrix(X85pc[dat85$fips %in% dat$fips[dat$year %in% samp],]),
                       new.param = as.matrix(Xp85[dat85$fips %in% dat$fips[dat$year %in% samp],]),
                       fe.newX = newfips)
    tm <- data.frame(fips = newfips, year = newyear, pred = predoos)
    with_holes <- full_join(mergetemplate, tm)
    return(with_holes$pred)
  }

  pnnpred[pnnpred < 0] <- 0
  pnnpred[pnnpred > 200] <- 200
  predframe$PNN85 <- rowMeans(pnnpred)
  saveRDS(pnnpred, file = paste0(outdir, "/",tag, "_PNN_bagframe_85.Rds"))


  print(tag)
  save(predframe, file = paste0(outdir, "/",tag, "_predframe.Rda"))
}

###########################
# characterize the models by hottness/dryness/wetness

# start with baseline data
Ep <- mean(dat$TAP)
meanmax <- mean(as.matrix(dat[,grepl("maxat", colnames(dat))]) )
meanmin <- mean(as.matrix(dat[,grepl("minat", colnames(dat))]) )
meanmaxrh <- mean(as.matrix(dat[,grepl("maxrh", colnames(dat))]) )
meanminrh <- mean(as.matrix(dat[,grepl("minrh", colnames(dat))]) )
wind <- mean(as.matrix(dat[,grepl("wind", colnames(dat))]) )
light <- mean(as.matrix(dat[,grepl("sdsfia", colnames(dat))]) )
dwp <- mean(dat[,grepl("precip", colnames(dat))] %>% as.matrix !=0)
prop_over_30 <- (dat[,paste0("SR",30:40)] %>% sum) /(dat[,paste0("SR",0:29)] %>% sum)
tosweep <- c(meanmax, meanmin, meanmaxrh, meanminrh, Ep, dwp, prop_over_30, wind, light)


# do same for other panels
if (laptop){
  fi <- list.files(path = outdir, pattern = "panel_GCB", full.names = T)
}
if(AWS){
  fi <- list.files(path = "/home/ubuntu/MACA/panels", pattern = "panel_GCB", full.names = T)
}

model_sumstats <- foreach(i = 1:length(fi), .combine = rbind, .errorhandling = "remove") %do% {
  x <- readRDS(fi[i])
  tag <- (fi[i] %>% strsplit("GCB_") %>% unlist)[2]
  RCP <- substr(tag, 4, 5)
  mod <- (tag %>% gsub(pattern = ".Rds", replacement = "") %>% strsplit("_") %>% unlist)[3]
  #TAP
  x$TAP <- rowSums(x[,grepl('precip', colnames(x))])/1000
  EP_40 <- mean(x$TAP[x$year %in% 2040:2069])
  EP_70 <- mean(x$TAP[x$year %in% 2070:2099])
  # dwp
  dwp40 <- mean(x[x$year %in% 2040:2069,grepl("precip", colnames(x))] %>% as.matrix !=0)
  dwp70 <- mean(x[x$year %in% 2070:2099,grepl("precip", colnames(x))] %>% as.matrix !=0)
  # maxat
  meanmax40 <- mean(as.matrix(x[x$year %in% 2040:2069,grepl("maxat", colnames(x))]) )
  meanmax70 <- mean(as.matrix(x[x$year %in% 2070:2099,grepl("maxat", colnames(x))]) )
  meanmin40 <- mean(as.matrix(x[x$year %in% 2040:2069,grepl("minat", colnames(x))]) )
  meanmin70 <- mean(as.matrix(x[x$year %in% 2070:2099,grepl("minat", colnames(x))]) )
  #SR
  x[,paste0("SR",0:40)][is.na(x[,paste0("SR",0:40)])] <- 0
  prop_over_30_40 <- (x[x$year %in% 2040:2069,paste0("SR",30:40)] %>% sum) /(x[x$year %in% 2040:2069,paste0("SR",0:29)] %>% sum)
  prop_over_30_70 <- (x[x$year %in% 2070:2099,paste0("SR",30:40)] %>% sum) /(x[x$year %in% 2070:2099,paste0("SR",0:29)] %>% sum)
  # RH
  meanmaxrh40 <- mean(as.matrix(x[x$year %in% 2040:2069,grepl("maxrh", colnames(x))]) )
  meanmaxrh70 <- mean(as.matrix(x[x$year %in% 2070:2099,grepl("maxrh", colnames(x))]) )
  meanminrh40 <- mean(as.matrix(x[x$year %in% 2040:2069,grepl("minrh", colnames(x))]) )
  meanminrh70 <- mean(as.matrix(x[x$year %in% 2070:2099,grepl("minrh", colnames(x))]) )
  # wind
  wind40 <- mean(as.matrix(x[x$year %in% 2040:2069,grepl("wind", colnames(x))]) )
  wind70 <- mean(as.matrix(x[x$year %in% 2070:2099,grepl("wind", colnames(x))]) )
  # light
  light40 <- mean(as.matrix(x[x$year %in% 2040:2069,grepl("sdsfia", colnames(x))]) )
  light70 <- mean(as.matrix(x[x$year %in% 2070:2099,grepl("sdsfia", colnames(x))]) )
  out <- data.frame(mod, RCP, p = c(40, 70),
                    meanmax = c(meanmax40, meanmax70),
                    meanmin = c(meanmin40, meanmin70),
                    meanmaxrh = c(meanmaxrh40, meanmaxrh70),
                    meanminrh = c(meanminrh40, meanminrh70),
                    EP = c(EP_40, EP_70),
                    dwp = c(dwp40, dwp70),
                    prop_over_30 = c(prop_over_30_40, prop_over_30_70),
                    wind = c(wind40, wind70),
                    light = c(light40, light70)
                    )
  print(out)
  return(out)
}
saveRDS(model_sumstats, file = paste0(outdir, "/model_sumstats_GCB.Rds"))
model_sumstats <- readRDS(file = paste0(outdir, "/model_sumstats_GCB.Rds"))
model_difstats <- model_sumstats
model_difstats[,4:ncol(model_difstats)] <- sweep(model_difstats[,4:ncol(model_difstats)], 2, tosweep, "-")


#########################
# make figure showing how PNN-OLS difference depends on severity of model

# collect means  of models/time periods/estimators
predfiles <- list.files(pattern = "_predframe", path = "/home/andrew/Dropbox/USDA/ARC/output/", full.names = T)
predfiles <- predfiles[!grepl("v2", predfiles)]
allpreds <- foreach(i = predfiles, .combine = rbind) %do%{
  pf <- load_obj(i)
  pf$tag <- ((i %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[1]
  pf
}
head(allpreds)
meanpred <- foreach(i = predfiles, .combine = "+") %do%{
  pf <- load_obj(i)
  toadd <- as.matrix(pf[,3:8]) /13
  if (i == predfiles[1]){idx <- pf[,1:2]}
  toadd
}
meanpred <- cbind(idx, meanpred)
meanpred$tag <- "unweighted ensemble"

preds <- rbind(allpreds, meanpred)
preds$tag <- factor(preds$tag)

preds_g <- gather(preds, estimator, pred, OLS45:PNN85)
saveRDS(preds_g, file = paste0(outdir, "/preds_gathered_GCB.Rds"))
preds_g <- readRDS(file = paste0(outdir, "/preds_gathered_GCB.Rds"))
preds_g_late <- summaryBy(pred~ tag + estimator, keep.names = T, data = preds_g[preds_g$year %in% 2070:2099,])
preds_g_early <- summaryBy(pred~ tag + estimator, keep.names = T, data = preds_g[preds_g$year %in% 2040:2069,])
preds_g_early$p <- 40
preds_g_late$p <- 70
preds <- rbind(preds_g_early, preds_g_late)
preds$RCP <- ifelse(grepl("45", preds$estimator), 45, 85)
preds <- preds[!grepl("weighted|bagged", preds$estimator),]
preds$estimator <- substr(preds$estimator, 1, 3)
preds <- spread(preds, key = estimator, value = pred)
preds$diff <- preds$PNN - preds$OLS
preds$RCP <- as.factor(preds$RCP)
colnames(preds)[colnames(preds) == "tag"] <- "mod"
model_sumstats <- readRDS(file = paste0(outdir, "/model_sumstats_GCB.Rds"))
model_difstats <- model_sumstats
model_difstats[,4:ncol(model_difstats)] <- sweep(model_difstats[,4:ncol(model_difstats)], 2, tosweep, "-")
model_difstats <- inner_join(model_difstats, preds)
model_difstats$p[model_difstats$p == 40] <- "2040-69"
model_difstats$p[model_difstats$p == 70] <- "2070-99"
model_difstats$p <- as.factor(model_difstats$p)
colnames(model_difstats)[colnames(model_difstats) == "p"] <- "period"
gdif <- gather(model_difstats, key = "var", value = "value", meanmax:light)
gdif$OLS <- gdif$PNN <- NULL
gdif$var <- factor(gdif$var, levels = c("meanmax", "meanmin", "prop_over_30","meanmaxrh", "meanminrh","light","EP","dwp","wind"))
levels(gdif$var) <- c("Mean Max Temp", "Mean Min Temp", "Proportion GDD 30+", "Mean Max RelH", "Mean Min RelH",
                      "Mean Sunlight (w/m^2)", "Mean Precip (m) ",
                      "Proportion Days With Precip", "Mean Wind Speed")

ggplot(gdif, aes(x=value, y=diff))+ geom_point(aes(shape = RCP, col = period)) + 
  facet_wrap(~var, scales = "free_x") + xlab("difference from 1979-2016 average") +
  ylab("Yield prediction difference between SNN and OLS")+
  geom_vline(xintercept = 0, lty = 2)+geom_hline(yintercept = 0, lty  =2)
dev.copy2pdf(file = "/home/andrew/Dropbox/USDA/yield_forecasting/tex/diff_by_var.pdf", height = 7, width = 8)
system("evince /home/andrew/Dropbox/USDA/yield_forecasting/tex/diff_by_var.pdf")
sdif <- spread(gdif, key = var, value = value)
sdif[,5:ncol(sdif)] %>% cor %>% Matrix %>% image

forcor <- model_difstats[,c(4:5, 10 , 6:7, 12,8, 9, 11)] %>% cor
colnames(forcor) <- rownames(forcor) <- c("Mean Max Temp", "Mean Min Temp", "Proportion GDD 30+", "Mean Max RelH", "Mean Min RelH",
                                                   "Mean Sunlight (w/m^2)", "Mean Precip (m) ",
                                                   "Proportion Days With Precip", "Mean Wind Speed")
forcor <- round(forcor,2)*lower.tri(forcor, diag = T)
forcor[forcor == 0] <- NA
library(xtable)
xtable(forcor)

##############################
# Make predictions for future periods
dtreg <- gam(yield~s(year, k=4), method = "REML", data = dat)#, family = Gamma(link = "log"))
plot(dtreg)
summary(dtreg)
dat$yield_dt <- (predict(dtreg, newdata = data.frame(year = 2016)) + residuals(dtreg))

mt <- data.frame(year = 1979:2016)
# mt <- data.frame(year = 1:nrow(dat))
Ni <- foreach(i = 1:96, .combine = cbind)%do% {
  set.seed(i, kind = "L'Ecuyer-CMRG")
  samp <- sample(unique(dat$year), replace = TRUE)
  # bsamp <- foreach(y = samp, .combine = c) %do% {which(dat$year == y)}
  d <- data.frame(table(samp))
  colnames(d) <- c("year", paste0("y",i))
  d$year <- as.numeric(as.character(d$year))
  full_join(mt, d)[,2]
}
Ni[is.na(Ni)] <- 0

predfiles <- list.files(pattern = "_predframe", path = "/home/andrew/Dropbox/USDA/ARC/output/", full.names = T)
predfiles <- predfiles[!grepl("v2", predfiles)]
allpreds <- foreach(i = predfiles, .combine = rbind) %do%{
  pf <- load_obj(i)
  pf$tag <- ((i %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[1]
  pf
}
head(allpreds)
p40 <- summaryBy(OLS45+OLS85+OLS_bagged45+OLS_bagged85+PNN45+PNN85~tag, 
                 data = allpreds[allpreds$year %in% 2040:2069,],
                 keep.names = T)
p70 <- summaryBy(OLS45+OLS85+OLS_bagged45+OLS_bagged85+PNN45+PNN85~tag, 
                 data = allpreds[allpreds$year %in% 2070:2099,],
                 keep.names = T)

s70 <- s40 <- p40
s70[,-1] <- s40[,-1] <- NA

bagfiles <- list.files(pattern = "_bagframe", path = "/home/andrew/Dropbox/USDA/ARC/output/", full.names = T)
# JAB estimates
for (i in 1:length(bagfiles)){
  b <- readRDS(bagfiles[i])
  tag <- ((bagfiles[i] %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[1]
  est <- ((bagfiles[i] %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[c(2,4)] %>% 
    paste0(collapse = "") %>% gsub(pattern = ".Rds", replacement = "") %>% 
    gsub(pattern = "SR", replacement = "OLS")
  bb <- b[pf$year %in% 2040:2069,] %>% colMeans
  jab <- foreach(j = 1:38, .combine= "+")%dopar%{
    (mean(bb[Ni[j,]==0]) - mean(bb))^2
  }*(37/38)
  U <- (exp(1)-1)*var(bb)/96*38
  sd40 <- sqrt(jab - U)

  bb <- b[pf$year %in% 2070:2099,] %>% colMeans
  jab <- foreach(j = 1:38, .combine= "+")%dopar%{
    (mean(bb[Ni[j,]==0]) - mean(bb))^2
  }*(37/38)
  U <- (exp(1)-1)*var(bb)/96*38
  sd70 <- sqrt(jab - U)
  s40[which(s40$tag == tag), grepl(est, colnames(s40))] <- sd40
  s70[which(s70$tag == tag), grepl(est, colnames(s70))] <- sd70
  print(s40)
  print(s70)
}
s70[,4:5] <- s70[,2:3]
s40[,4:5] <- s40[,2:3]
# normal bootstrap SE's for unbagged estimates
s70[,2:3] <- NA
s40[,2:3] <- NA

for (i in 1:length(bagfiles[grepl("SR", bagfiles)])){
  b <- readRDS(bagfiles[grepl("SR", bagfiles)][i])
  tag <- ((bagfiles[grepl("SR", bagfiles)][i] %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[1]
  est <- ((bagfiles[grepl("SR", bagfiles)][i] %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[c(2,4)] %>% 
    paste0(collapse = "") %>% gsub(pattern = ".Rds", replacement = "") %>% 
    gsub(pattern = "SR", replacement = "OLS")
  sd40 <- b[pf$year %in% 2040:2069,] %>% colMeans %>% sd
  sd70 <- b[pf$year %in% 2070:2099,] %>% colMeans %>% sd
  s40[which(s40$tag == tag), grepl(est, colnames(s40))] <- sd40
  s70[which(s70$tag == tag), grepl(est, colnames(s70))] <- sd70
  print(s40)
  print(s70)
}

colnames(p70) <- gsub("PNN", "SNN", colnames(p70))
colnames(s70) <- gsub("PNN", "SNN", colnames(s70))
colnames(p40) <- gsub("PNN", "SNN", colnames(p40))
colnames(s40) <- gsub("PNN", "SNN", colnames(s40))

plotfun <- function(RCP, period, x = FALSE, y = FALSE, leg = FALSE){
  par(xpd = T)
  plot(idx, p40$OLS45, ylim = c(0,175), xlim = c(0, max(idx)+3), pch = 15, cex = ptsize, 
       col = "blue", xlab = "", xaxt = "n", yaxt = "n", ylab= "")
  abline(h = mean(dat$yield_dt), lty = 2)
  dd <- get(paste0("p",period))
  ss <- get(paste0("s",period))
  mods <- c("OLS", "OLS_bagged", "SNN")
  cols <- c("blue", "darkgreen","red")
  for (mod in 1:3){
    for(i in 1:13){segments(x0 = idx[i]+mod/2, 
                            y0 = dd[i,paste0(mods[mod],RCP)] + (ss[i,paste0(mods[mod],RCP)])*1.96,
                            y1 = dd[i,paste0(mods[mod],RCP)] - (ss[i,paste0(mods[mod],RCP)])*1.96,
                            col = cols[mod], lwd = lsize)
    }    
  }
  mtext(side = 3, line = -1.5, paste0(period,"-",period+29, ", RCP ", RCP))
  if (x == TRUE){
    mtext(side = 1, cex=1, at=idx+1, p40$tag, xpd=TRUE, las = 2)
  }
  if (y == T){
    mtext(side = 2, cex=1, line = .25,at=seq(25, 150, 25), seq(25, 150, 25), xpd=TRUE)
    mtext(side = 2, line = 1.5, "Yield (bu/ac)")
  }
  if(leg == TRUE){
    legend(x = 1, y = 100, legend = gsub("_", " ", mods), lwd = 2, col = cols, bty = "n", y.intersp = 1.5)
  }
}

par(mfrow = c(2,2))
par(mar = c(0,0,0,0))
par(oma = c(9.8,4.1, .5, .5))
idx <- (1:13)*3-2
ptsize <- 0
lsize <- 2
plotfun(45, 40, y = T, leg = T)
plotfun(45, 70)
plotfun(85, 40, x = T, y = T)
plotfun(85, 70, x = T)
dev.copy2pdf(file = "/home/andrew/Dropbox/USDA/yield_forecasting/tex/yield_by_model_JAB.pdf", height = 8, width = 8)
system("evince /home/andrew/Dropbox/USDA/yield_forecasting/tex/yield_by_model_JAB.pdf")

#########
# same figure, but for interannual
pf <- load_obj(predfiles[1])
fo <- colnames(pf)[colnames(pf) %ni% c("fips", "year", "yield", "yield_dt")] %>% 
  paste0(collapse = "+") %>% 
  paste0("~year") %>% 
  as.formula

smooth_list <- foreach(i = 1:length(predfiles)) %do% {
  pf <- load_obj(predfiles[i])
  pf <- pf[with(pf, order(year, fips)),]
  d <- summaryBy(fo, data = pf, keep.names = T)
  fitframe <- foreach(j = 3:8, .combine = cbind) %do% {
    foj <- as.formula(paste0(colnames(pf)[j], "~s(year)"))
    m <- gam(foj, data = d, method = "REML")
    # plot(m)
    p <- predict(m, newdata = data.frame(year = 2006:2099), se.fit = T)
    out <- data.frame(mu = p$fit, sd = p$se.fit)
    colnames(out) <- paste0(colnames(pf)[j], "_",colnames(out))
    return(out)
  }
  return(data.frame(year = 2006:2099, fitframe))
}

smooth_dump <- foreach(i = 1:length(predfiles), .combine = rbind) %do% {
  x <- smooth_list[[i]]
  x$tag <- ((predfiles[i] %>% strsplit("//_") %>% unlist)[2] %>% strsplit("_") %>% unlist)[1]
  x
}
colnames(smooth_dump) <- gsub("PNN","SNN", colnames(smooth_dump))
fo <- as.formula(paste0(paste(colnames(smooth_dump)[2:13], collapse = "+"), "~tag"))
smooth_2040 <- summaryBy(fo, data = smooth_dump[smooth_dump$year %in% 2040:2069,], keep.names = T)
smooth_2070 <- summaryBy(fo, data = smooth_dump[smooth_dump$year %in% 2070:2099,], keep.names = T)
o <- order(smooth_2070$SNN85_mu, decreasing =T)
smooth_2040 <- smooth_2040[o,]
smooth_2070 <- smooth_2070[o,]

smooth_2070$SNN85_mu / smooth_2070$OLS85_mu
smooth_2070$SNN85_mu - smooth_2070$OLS85_mu
# percentage point difference
(1-58.1/159.1) - (1-32.5/159.1)

plotfun <- function(RCP, period, x = FALSE, y = FALSE, leg = FALSE){
  par(xpd = T)
  plot(idx, smooth_2040$OLS45_mu, ylim = c(20,175), xlim = c(0, max(idx)+3), pch = 15, cex = ptsize, 
       col = "blue", xlab = "", xaxt = "n", yaxt = "n", ylab= "")
  abline(h = mean(dat$yield_dt), lty = 2)
  dd <- get(paste0("smooth_",period))
  mods <- c("OLS", "OLS_bagged", "SNN")
  cols <- c("blue", "darkgreen","red")
  for (mod in 1:3){
    for(i in 1:13){segments(x0 = idx[i]+mod/2, 
                            y0 = dd[i,paste0(mods[mod],RCP,"_mu")] + dd[i,paste0(mods[mod],RCP,"_sd")]*1.96, 
                            y1 = dd[i,paste0(mods[mod],RCP,"_mu")] - dd[i,paste0(mods[mod],RCP,"_sd")]*1.96, 
                            col = cols[mod], lwd = lsize)
    }    
  }
  mtext(side = 3, line = -1.5, paste0(period,"-",period+29, ", RCP ", RCP))
  if (x == TRUE){
    mtext(side = 1, cex=1, at=idx+1, smooth_2040$tag, xpd=TRUE, las = 2)
  }
  if (y == T){
    mtext(side = 2, cex=1, line = .25,at=seq(25, 150, 25), seq(25, 150, 25), xpd=TRUE)
    mtext(side = 2, line = 1.5, "Yield (bu/ac)")
  }
  if(leg == TRUE){
    legend(x = 1, y = 100, legend = gsub("_", " ", mods), lwd = 2, col = cols, bty = "n", y.intersp = 1.5)
  }
}

par(mfrow = c(2,2))
par(mar = c(0,0,0,0))
par(oma = c(9.8,4.1, .5, .5))
idx <- (1:13)*3-2
ptsize <- 0
lsize <- 2
plotfun(45, 2040, y = T, leg = T)
plotfun(45, 2070)
plotfun(85, 2040, x = T, y = T)
plotfun(85, 2070, x = T)
dev.copy2pdf(file = "/home/andrew/Dropbox/USDA/yield_forecasting/tex/yield_by_model.pdf", height = 8, width = 8)
system("evince /home/andrew/Dropbox/USDA/yield_forecasting/tex/yield_by_model.pdf")




