



#######################################
# variable importance
rm(list=ls())
gc()
gc()
"%ni%" <- Negate("%in%")
mse <- function(x, y){mean((x-y)^2)}

library(devtools)
install_github("cranedroesch/panelNNET", ref = "earlystop")
library(panelNNET)
library(doParallel)
library(dplyr)

AWS <- grepl('ubuntu', getwd())
desktop <- grepl(':', getwd())
laptop <- grepl('/home/andrew', getwd())
if(AWS){
  setwd("/home/ubuntu/projdir")
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
dat <- dat[dat$prop_irr < .5,]

# make TAP variable
dat$TAP <- rowSums(dat[,grepl('precip', colnames(dat))])/1000
dat$TAP2 <- dat$TAP^2
# nonparametric
X <- dat[,grepl('year|SR|soil_|precip|sdsfia|wind|minat|maxat|minrh|maxrh|lat|lon|prop_irr',colnames(dat))]
X <- X[sapply(X, sd) > 0]
# parametric
dat$y <- dat$year - min(dat$year) + 1
dat$y2 <- dat$y^2

Xp <- cbind(dat[,c("y", "y2", "TAP", "TAP2")], dat[,grepl("SR", colnames(dat))])
Xp <- Xp[sapply(Xp, sd) > 0]

# faster predict function for PCA
ppc <- function(x, PC, ncomp){
  x <- x %>% sweep(2, PC$center, "-")%>% sweep(2, PC$scale, "/")
  MatMult(as.matrix(x), PC$rotation[,1:ncomp])
}


###############################
# Variable importance
weights0 <- readRDS(file = "ensemble_weights_GCB.Rds")

ensemble_mse <- 246.6626
nfolds <- 96
#
varImp <- function(pterms, npterms, tag, reps = 5){
  pterms <- colnames(Xp)[colnames(Xp) %in% pterms]
  npterms <- colnames(X)[colnames(X) %in% npterms]
  out <- foreach (j = 1:reps, .combine = c)%do%{
    set.seed(round(runif(1, 0, 1e6)), kind = "L'Ecuyer-CMRG")
    scramp <- Xp[sample(1:nrow(Xp)),pterms]
    scramn <- X[sample(1:nrow(X)),npterms]
    Xps <- Xp
    Xs <- X
    Xps[,pterms] <- scramp
    Xs[,npterms] <- scramn
    # PNN OOB
    fi <- list.files(path = outdir, pattern = "v6", full.names = TRUE)
    # Make a foldkey
    foldkey <- fi %>% strsplit("dry") %>% sapply(function(x){unlist(x)[2]}) %>% 
      gsub(pattern = "_", replacement = "") %>% 
      gsub(pattern = ".Rda", replacement = "") %>% 
      as.numeric
    mergetemplate <- dat[,c("fips", "year")]
    mergetemplate <- mergetemplate[with(mergetemplate, order(fips, year)),]
    pnnoob <- foreach(i = 1:length(fi), .combine = cbind) %dopar% {
      set.seed(foldkey[i], kind = "L'Ecuyer-CMRG")
      samp <- sample(unique(dat$year), replace = TRUE)
      bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
      oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
      obj <- load_obj(fi[i])
      PC <- readRDS(paste0(outdir,"/PCA_", foldkey[i],".Rds"))
      Xpc <- ppc(X[bsamp,], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95)) %>% as.data.frame
      Xs <- ppc(Xs[dat$year %in% oosamp  & dat$fips %in% dat$fips[dat$year %in% samp],], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
      parlist <- obj$parlist
      pnn <- panelNNET(y = dat$yield[bsamp],
                       X = Xpc,
                       hidden_units = rep(100, 10),
                       fe_var = dat$fips[bsamp],
                       maxit = 0,
                       time_var = dat$year[bsamp],
                       param = Xp[bsamp,],
                       activation = 'lrelu',
                       parlist = parlist
                       )
      newfips <- dat$fips[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]]
      newyear <- dat$year[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]]
      predoos <- predict(pnn, newX = as.matrix(Xs),
                         new.param = as.matrix(Xps[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp],]),
                         fe.newX = newfips)
      tm <- data.frame(fips = newfips, year = newyear, pred = predoos)
      with_holes <- full_join(mergetemplate, tm, by = c("fips", "year"))
      with_holes <- with_holes[with(with_holes, order(fips, year)),]  
      return(with_holes$pred)
    }
    # pnnoob_pred <- rowMeans(pnnoob, na.rm = T)
    pnnoob_pred <- apply(pnnoob, 1, weighted.mean, w = weights0, na.rm = T)
    mergetemplate$pnnoob <- pnnoob_pred
    pnn_result <- full_join(dat[,c("fips", "year", "yield")], mergetemplate, by = c("fips", "year"))
    perm_MSE <- with(pnn_result, mse(yield, pnnoob))
    print(perm_MSE)
    return(perm_MSE)
  }
  outdf <- data.frame(tag = tag, MSE_diff = mean(out - ensemble_mse), MSE_sd = sd(out - ensemble_mse))
  print(outdf)
  save(outdf, file = paste0(outdir, "/",tag,".Rda"))
}



dvars <- c("precip", "sdsfia", "wind", "minat", "maxat", "minrh", "maxrh")
days <- 60:304
stops <- days[(days %% 10)==5]
dvcombs <- foreach(v=dvars, .combine = c) %do% {
  foreach(s = stops) %do% {
    span <- c((s - 5):(s+5))
    paste0(v,"_jday_",span)
  }
}
all_vars_span <- foreach(s = stops) %do% {
  span <- c((s - 5):(s+5))
  colnames(dat)[grepl(paste0("_jday_",span, collapse = "|"), colnames(dat))]
}
all_spans_var <- foreach(v = dvars) %do% {
  colnames(dat)[grepl(v, colnames(dat))]
}

NPtodo <- c(dvcombs, all_vars_span, all_spans_var, list(
          c("lat", "lon"),
          "prop_irr",
          colnames(X)[grep("soil", colnames(X))],
          "year", #NP only
          "year", # both
          NULL, # P only
          colnames(X)[grep("SR", colnames(X))], # NP only
          colnames(X)[grep("SR", colnames(X))], # both
          NULL, # P SR
          NULL # P TAP
          ))
Ptodo <- c(vector("list", length(dvcombs)), vector("list", length(all_vars_span)), vector("list", length(all_spans_var)), list(
            NULL, #c("lat", "lon"),
            NULL, #"prop_irr",
            NULL,#colnames(X)[grep("soil", colnames(X))],
            NULL,#"year", #NP only
            c("y", "y2"), #"year", # both
            c("y", "y2"), #NULL, # P only
            NULL, #colnames(X)[grep("SR", colnames(X))], # NP only
            colnames(X)[grep("SR", colnames(X))], # both
            colnames(X)[grep("SR", colnames(X))], #NULL, # P SR
            c("TAP", "TAP2")#NULL # P TAP
))
dvtags <- expand.grid(stops, dvars)[,c(2,1)]
dvtags <- paste(dvtags[,1], dvtags[,2], sep = "_")
stoptags <- paste0("span_", stops)
vartags <- paste0("all_", dvars)
tags <- c(dvtags, stoptags, vartags, 
          "latlon", 
          "prop_irr",
          "soil",
          "year_NP",
          "year_both",
          "year_P",
          "SR_NP",
          "SR_both",
          "SR_P",
          "TAP"
)
tags <- paste0("permImpV3_GCB_", tags)


final_permimp <- foreach(k = 1:length(tags), .combine = rbind) %do%{
  PT <- proc.time()
  vi <- varImp(pterms = Ptodo[[k]],
         npterms = NPtodo[[k]],
         tag = tags[k],
         reps = 2)
  print(proc.time() - PT)
  return(vi)
}

varImp(pterms = NULL, npterms = colnames(X), tag = "permImpV3_GCB_ALLNP", reps = 2)
varImp(pterms = colnames(Xp), npterms = NULL, tag = "permImpV3_GCB_ALLParm", reps = 2)
varImp(pterms = colnames(Xp)[colnames(Xp) %ni% c("y", "y2")], npterms = NULL, tag = "permImpV3_GCB_ParmNotY", reps = 2)


###########################
# figures
library(panelNNET)
fi <- list.files("/home/andrew/Dropbox/USDA/ARC/output/", pattern = "permImpV3")

# time-varying terms
tv <- fi[grepl("precip|sdsfia|wind|minat|maxat|minrh|maxrh|span", fi) & !grepl("all", fi)]
dvars <- c("precip", "sdsfia", "wind", "minat", "maxat", "minrh", "maxrh", "span")
toplot <- foreach(i = tv, .combine = rbind) %do%{
  load_obj(paste0("/home/andrew/Dropbox/USDA/ARC/output/", i))
}
toplot$day <- (toplot$tag %>% as.character %>% strsplit(split = "_")) %>% sapply(function(x){x[[4]]}) %>% as.numeric
toplot <- toplot[with(toplot, order(day)),]
par(mar = c(5,4,.5,2)+.1)
plot(1, cex = 0, xlim = c(60, 300), ylim = range(sqrt(toplot$MSE_diff), na.rm=T), xlab = "day", ylab = "Change in RMSE")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
for (v in dvars){
  x <- toplot[grepl(v, toplot$tag),]
  x <- x[with(x, order(day)),]
  x$MSE_diff[x$MSE_diff < 0] <- 0
  if (v == "span"){lty = 2; col = "black"} else {lty = 1; col = which(v == dvars)}
  with(x, lines(day, sqrt(MSE_diff), col = col, lty = lty, lwd = 2, ylim = range(sqrt(toplot$MSE_diff), na.rm=T)))
}
legend(x = 60, y = 3.3, lty = c(rep(1, 7), 2), col = c(1:7, "black"), 
       legend <- c("precip", "sunlight", "wind", "min temp", "max temp", "min relh", "max relh", "All"),
       y.intersp = 2.2, bty = "n")
dev.copy2pdf(file = "/home/andrew/Dropbox/USDA/yield_forecasting/tex/varimp_tv_GCB.pdf", height = 5, width = 5)
system("evince /home/andrew/Dropbox/USDA/yield_forecasting/tex/varimp_tv_GCB.pdf")

# non-TV
ntv <- fi[fi%ni%tv]
toplot <- foreach(i = ntv, .combine = rbind) %do%{
  load_obj(paste0("/home/andrew/Dropbox/USDA/ARC/output/", i))
}
toplot$MSE_diff <- sqrt(toplot$MSE_diff)
toplot <- rbind(toplot[!grepl("_ALL|ParmNot", toplot$tag),], toplot[grepl("_ALL|ParmNot", toplot$tag),])
toplot$vname <- c(
  "All parametric terms except time",
  "All parametric terms",
  "All nonparametric terms",
  "Time (parametric)",
  "Time (nonparametric)",
  "Time (both)",
  "Total Annual Precip",
  "GDD (parametric)",
  "GDD (nonparametric)",
  "GDD (both)",
  "Soil",
  "Proportion Irrigated",
  "Lat/Lon",
  "Wind",
  "Sunlight",
  "Precip (daily)",
  "Min Rel Humidity",
  "Max Rel Humidity",
  "Min Temp",
  "Max Temp"
) %>% rev
toplot <- toplot[with(toplot, order(MSE_diff)),]
toplot <- toplot[toplot$tag != "permImpV3_GCB_ParmNotY",]
par(mar = c(5.1, 11, 2.1, 2.1))
dval <- 20
barplot(toplot$MSE_diff,
        names.arg = toplot$vname,
        horiz = T, las = 2,
        col = c("blue", "blue", "green", "green", "blue", rep("red", 8), "blue", rep("red",5)) %>% rev,
        density = c(-1, rep(dval,4), -1, rep(dval, 13)) %>% rev,
        xlab = "Change in RMSE")
dev.copy2pdf(file = "/home/andrew/Dropbox/USDA/yield_forecasting/tex/varimp_ntv_GCB.pdf", height = 5, width = 5)
system("evince /home/andrew/Dropbox/USDA/yield_forecasting/tex/varimp_ntv_GCB.pdf")
