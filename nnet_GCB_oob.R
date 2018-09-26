
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
cr <- "corn"
irrdry <- "dry"
# load data
if (cr == "corn"){dat <- readRDS("panel_corn.Rds")}
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

# faster predict function for PCA
ppc <- function(x, PC, ncomp){
  x <- x %>% sweep(2, PC$center, "-")%>% sweep(2, PC$scale, "/")
  MatMult(as.matrix(x), PC$rotation[,1:ncomp])
}

# OLS baseline
mm <- model.matrix(as.formula(paste0("~yield+y+y2+TAP+TAP2+",
                                     paste(colnames(dat)[grepl("SR", colnames(dat))], collapse = "+"),
                                     "-1")),
                   data = dat)
nfolds <- 96
SRoos <- foreach(i = 1:nfolds, .combine = c) %dopar% {
  set.seed(i, kind = "L'Ecuyer-CMRG")
  samp <- sample(unique(dat$year), replace = TRUE)
  bsamp <- foreach(y = samp, .combine = c) %do% {which(dat$year == y)}
  oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
  mmis <- mm[bsamp,]
  dis <- as.data.frame(demeanlist(mmis, list(dat$fips[bsamp])))
  m <- glmnet(y = dis$yield, x = as.matrix(dis[,-1]), lambda = 0, intercept = FALSE)
  XB <- mm[dat$year %in% oosamp,-1] %*% coef(m)[-1,]
  fe <- (dat$yield[bsamp]-dis$yield) - (mmis[,-1]-as.matrix(dis[,-1])) %*% coef(m)[-1,]
  tm <- data.frame(fips = dat$fips[bsamp], fe = fe)
  tm <- tm[!duplicated(tm),]
  om <- data.frame(XB, fips = dat$fips[dat$year %in% oosamp], yield = dat$yield[dat$year %in% oosamp])
  pred <- merge(om, tm)
  pred$pred <- with(pred, fe+XB)
  mse((pred$yield), (pred$pred))
}
mean(SRoos)

# out-of-bag predictions and errors
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
SRoob_pred <- rowMeans(SRoob, na.rm = T)
mergetemplate$SRoob_pred <- SRoob_pred
SR_result <- full_join(dat[,c("fips", "year", "yield")], mergetemplate)
head(SR_result)
with(SR_result, mse(SRoob_pred, yield))
#################################
# stacking to get weights

# Make a vector of the filenames of all of the parameters
fi <- list.files(path = outdir, pattern = "v13", full.names = TRUE)

# Make a foldkey
foldkey <- fi %>% strsplit("dry") %>% sapply(function(x){unlist(x)[2]}) %>% 
           gsub(pattern = "_", replacement = "") %>% 
           gsub(pattern = ".Rda", replacement = "") %>% 
           as.numeric

mergetemplate <- dat[,c("fips", "year")]
mergetemplate <- mergetemplate[with(mergetemplate, order(fips, year)),]

# PNN OOB
mergetemplate <- dat[,c("fips", "year")]
mergetemplate <- mergetemplate[with(mergetemplate, order(fips, year)),]
fi <- list.files(path = outdir, pattern = "optimv13", full.names = T)
stubs <- fi %>% strsplit(split = ".Rda") %>% unlist
foldkey <- stubs %>% strsplit(split = "_") %>% sapply(function(x){x[length(x)]}) %>% as.numeric
pnnoob <- foreach(i = 1:length(fi), .combine = cbind, .errorhandling = "remove") %dopar% {
  set.seed(foldkey[i], kind = "L'Ecuyer-CMRG")
  samp <- sample(unique(dat$year), replace = TRUE)
  bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
  oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
  PC <- readRDS(paste0(outdir,"/PCA_v13_", foldkey[i],".Rds"))
  Xpc <- ppc(X[bsamp,], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
  Xtest <- ppc(X[dat$year %in% oosamp  & dat$fips %in% dat$fips[dat$year %in% samp],], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
  obj <- load_obj(fi[i])
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
                   parlist = parlist)
  newfips <- dat$fips[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]]
  newyear <- dat$year[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]]
  predoos <- predict(pnn, newX = as.matrix(Xtest),
                     new.param = as.matrix(Xp[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp],]),
                     fe.newX = newfips)
  tm <- data.frame(fips = newfips, year = newyear, pred = predoos)
  with_holes <- full_join(mergetemplate, tm)
  with_holes <- with_holes[with(with_holes, order(fips, year)),]
  return(with_holes$pred)
}
save(pnnoob, file = paste0(outdir,"/pnnoob_GCB_v13.Rda"))
load(paste0(outdir,"/pnnoob_GCB_v13.Rda"))

MSE_result <- function(weights){
  pnnoob_pred <- apply(pnnoob, 1, weighted.mean, w = weights, na.rm = T)
  mergetemplate$pnnoob <- pnnoob_pred
  pnn_result <- full_join(dat[,c("fips", "year", "yield")], mergetemplate)
  with(pnn_result, mse(yield, pnnoob))
}
MSE_result(rep(1, 204))
MSE_result(c(weights0))

MSE_result(c(baseweights0))
MSE_result(log(stackweights))

unique(typetag)
mask <- grepl("opt", typetag)
MSE_result(c(weights0^0)*mask)

SR_result <- function(weights){
  pnnoob_pred <- apply(SRoob, 1, weighted.mean, w = weights, na.rm = T)
  mergetemplate$pnnoob <- pnnoob_pred
  pnn_result <- full_join(dat[,c("fips", "year", "yield")], mergetemplate)
  with(pnn_result, mse(yield, pnnoob))
}
SR_result((1:96)^0)


