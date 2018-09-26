
#######################################
# bayesian hyperparameter search on PCs, doing the PCA train/test split appropriately
# NO PARAMETRIC TERMS

rm(list=ls())
gc()
gc()
"%ni%" <- Negate("%in%")
mse <- function(x, y){mean((x-y)^2)}

library(devtools)
install_github("cranedroesch/panelNNET", ref = "fixing_nopar", force = F)
library(panelNNET)
library(doParallel)
library(doBy)
library(glmnet)
library(dplyr)
library(randomForest)

AWS <- grepl('ubuntu', getwd())
desktop <- grepl(':', getwd())
laptop <- grepl('/home/andrew', getwd())
if(AWS){
  setwd("/home/ubuntu/projdir")
  system("mkdir /home/ubuntu/projdir/outdir")
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

# Xp <- cbind(dat[,c("y", "y2", "TAP", "TAP2")])
Xp <- cbind(dat[,c("y", "y2", "TAP", "TAP2")], dat[,grepl("SR", colnames(dat))])
Xp <- Xp[sapply(Xp, sd) > 0]

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

# faster predict function for PCA
ppc <- function(x, PC, ncomp){
  x <- x %>% sweep(2, PC$center, "-")%>% sweep(2, PC$scale, "/")
  MatMult(as.matrix(x), PC$rotation[,1:ncomp])
}

fi <- list.files(path = outdir, pattern = "v8", full.names = T)
stubs <- fi %>% strsplit(split = ".Rda") %>% unlist
foldkey <- stubs %>% substr(nchar(stubs) -1, nchar(stubs)) %>% gsub(pattern = "_", replacement = "") %>% as.numeric
# foldkey <- NULL
# starting fit
registerDoParallel(detectCores())
msgs <- foreach(g = 1:nfolds, .errorhandling = "pass")%dopar% {
  set.seed(g, kind = "L'Ecuyer-CMRG")
  samp <- sample(unique(dat$year), replace = TRUE)
  bsamp <- foreach(y = samp, .combine = c) %do% {which(dat$year == y)}
  oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
  # data frame of combinations to try
  Dlam <- 2^seq(-8, 5, by = .1)
  Dgravity <- seq(1.00001, 1.3, by = .01)
  Dstart_LR <- 10^seq(-6, -1, by = .05)
  Dbatchsize <- seq(10, 500, by = 10)
  Dpp <- seq(0, 1, by = .01)
  DLRSR <- seq(1.1, 3, by = .01)
  Ddrop <- seq(.3, 2, by = .01)
  ns <- 6
  if (g %in% foldkey){
    print(paste0("restarting ", g))
    obj <- load_obj(fi[which(foldkey == g)])
    parlist<-bestparlist <- obj$parlist
    besterr <- obj$err
    besthp <- obj$hyperparm
    metapar <- obj$metapar
    PC <- readRDS(file = paste0(outdir, "/PCA_", g,".Rds"))
    Xpc <- ppc(X[bsamp,], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
    Xtest <- ppc(X[dat$year %in% oosamp  & dat$fips %in% dat$fips[dat$year %in% samp],], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
  } else { # if g is not in the foldkey, start from scratch
    #PCA
    if (paste0("PCA_", g,".Rds") %in% list.files(path = outdir)){
      PC <- readRDS(file = paste0(outdir, "/PCA_", g,".Rds"))
    } else {
      print(paste0("starting PCA ", g))
      PC <- prcomp(X[bsamp,], scale. = T, center = T, retx = FALSE)
      saveRDS(PC, file = paste0(outdir, "/PCA_", g,".Rds"))
    }
    Xpc <- ppc(X[bsamp,], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
    Xtest <- ppc(X[dat$year %in% oosamp  & dat$fips %in% dat$fips[dat$year %in% samp],], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
    #
    print(paste0("starting initial (over)fit ", g))
    if (g %in% foldkey){
      obj <- load_obj(fi[which(foldkey == g)])
      parlist<-bestparlist <- obj$parlist
      besterr <- obj$err
      besthp <- obj$hyperparm
    } else {
      pnn <- panelNNET(y = dat$yield[bsamp],
                       X = Xpc,
                       hidden_units = rep(100, 10),
                       fe_var = NULL,
                       maxit = 1000,
                       lam = .01,
                       time_var = dat$year[bsamp],
                       param = NULL,
                       verbose = F,
                       report_interval = 1,
                       gravity = 1.01,
                       convtol = 1e-3,
                       activation = 'lrelu',
                       start_LR = .01,
                       parlist = NULL,
                       OLStrick = TRUE,
                       OLStrick_interval = 25,
                       batchsize = 256,
                       maxstopcounter = 25,
                       LR_slowing_rate = 2,
                       parapen = NULL,
                       return_best = TRUE
      )
      bestparlist <- pnn$parlist
      predoos <- predict(pnn, newX = as.matrix(Xtest),
                         new.param = NULL,
                         fe.newX = NULL)
      besterr <- mse(predoos, dat$yield[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]])
      besthp <- NULL
    }
    metapar <- data.frame(lam = sample(Dlam, ns),
                          gravity = sample(Dgravity, ns),
                          start_LR = sample(Dstart_LR, ns),
                          batchsize = sample(Dbatchsize, ns),
                          LRSR = sample(DLRSR, ns),
                          drop = sample(Ddrop, ns),
                          err = NA,
                          improvement = NA)
    metapar$drop[metapar$drop>1] <- 1
    print(paste0("starting random search ", g))
    for (i in 1:nrow(metapar)){
      # print(metapar[i,])
      # print(paste(g, i))
      pnn <- tryCatch(panelNNET(y = dat$yield[bsamp],
                                X = Xpc,
                                hidden_units = rep(100, 10),
                                fe_var = NULL,
                                maxit = 200,
                                lam = metapar$lam[i],
                                time_var = dat$year[bsamp],
                                param = NULL,
                                verbose = F,
                                report_interval = 1,
                                gravity = metapar$gravity[i],
                                convtol = 1e-3,
                                activation = 'lrelu',
                                start_LR = metapar$start_LR[i],
                                parlist = bestparlist,
                                OLStrick = TRUE,
                                OLStrick_interval = 20,
                                batchsize = metapar$batchsize[i],
                                maxstopcounter = 25,
                                LR_slowing_rate = metapar$LRSR[i],
                                parapen = NULL, 
                                dropout_hidden = metapar$drop[i],
                                dropout_input = metapar$drop[i]^.321,
                                return_best = TRUE,
                                stop_early = list(check_every = 20,
                                                  max_ES_stopcounter = 5,
                                                  y_test = dat$yield[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]],
                                                  X_test = as.matrix(Xtest),
                                                  P_test = NULL,
                                                  fe_test = NULL)
      ), error = function(e){e})
      if (inherits(pnn, "error")){
        print(paste("problem at ", g))
        metapar$err[i] <- besterr
      } else {
        metapar$err[i] <- pnn$mse_test
      }
      # print(metapar$err[i])
      metapar$improvement[i] <- besterr - metapar$err[i]
      if (metapar$err[i] == besterr){
        # print("no progress down loss function")
      } else if (metapar$err[i] < besterr){
        # print("new low!")
        besterr <- metapar$err[i]
        bestparlist <- pnn$parlist
      } else if (metapar$err[i] > besterr){
        # print("got worse :-(")
      }
    }
  }
  # model and optim to figure out what works best
  idx <- Nruns <- 0
  print(paste0("starting bayesian search", g))
  while(idx<200 & (mean(metapar$improvement)>.001 | nrow(metapar)<99)){
    idx <- idx+1
    Nruns <- Nruns + 1
    prgrid <- expand.grid(Dlam %>% sample(4), 
                          Dgravity %>% sample(4),
                          Dstart_LR %>% sample(4),
                          Dbatchsize %>% sample(4),
                          DLRSR %>% sample(4), 
                          Ddrop %>% sample(4)
    ) %>% as.data.frame
    colnames(prgrid) <- colnames(metapar)[1:6]
    prgrid$drop[prgrid$drop>1] <- 1
    if (idx %% 3 %in% 1:2 & length(unique(metapar$improved))!=1){
      rf <- randomForest(y = (metapar$improvement), x = metapar[,1:6], importance = T)
      pr <- predict(rf, newdata = prgrid)
      # varImpPlot(rf)
      # print(rf)
      f <- function(x){predict(rf, newdata = x)}
      # if(idx%%3 == 1){
      # print("optim")
      o <- optim(par = prgrid[which.max(pr),1:6],
                 fn = f,
                 lower = c(0, 1.001, min(Dstart_LR), 2, min(DLRSR), .3),
                 upper = c(2^8, 2, 2, 500, 3, 1),
                 method = "L-BFGS-B",
                 control = list(fnscale = -1)
      )
      # print(o)
      newrow <- data.frame(t(o$par), err = NA)    
    } else {
      # print("random")
      newrow <- data.frame(prgrid[sample(1:nrow(prgrid), 1),], err = NA)
    }
    newrow$drop[newrow$drop>1] <- 1
    # print(newrow)
    # fit net with new metapar
    pnn <- tryCatch(panelNNET(y = dat$yield[bsamp],
                              X = Xpc,
                              hidden_units = rep(100, 10),
                              fe_var = NULL,
                              maxit = 200,
                              lam = newrow$lam,
                              time_var = dat$year[bsamp],
                              param = NULL,
                              verbose = F,
                              report_interval = 1,
                              gravity = newrow$gravity,
                              convtol = 1e-3,
                              activation = 'lrelu',
                              start_LR = newrow$start_LR,
                              parlist = bestparlist,
                              OLStrick = T,
                              OLStrick_interval = 20,
                              batchsize = round(newrow$batchsize),
                              maxstopcounter = 10,
                              LR_slowing_rate = newrow$LRSR,
                              dropout_hidden = newrow$drop,
                              dropout_input = newrow$drop^.321,
                              parapen = NULL,
                              stop_early = list(check_every = 20,
                                                max_ES_stopcounter = 4,
                                                y_test = dat$yield[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]],
                                                X_test = as.matrix(Xtest),
                                                P_test = NULL,
                                                fe_test = NULL)
    ), error = function(e){e})
    if (inherits(pnn, "error")){
      print(paste("problem at ", g))
      newrow$err <- besterr
    } else {
      newrow$err <- pnn$mse_test
    }
    newrow$improvement <- besterr - newrow$err
    if (newrow$err == besterr){
      # print("no progress down loss function")
    } else if (newrow$err < besterr){
      # print("new low!")
      besterr <- newrow$err
      bestparlist <- pnn$parlist
      besthp <- newrow
      tosave <- list(hyperparm = newrow, parlist = bestparlist, err = besterr, fold = g, metapar = metapar)
      save(tosave, file = paste0(outdir,"/tempbestparlist_GCB_nopar_optimv8_", cr, "_",irrdry,"_",g,".Rda"))
      idx <- 0
    } else if (newrow$err > besterr){
      # print("got worse :-(")
    }
    if (nrow(metapar) >100){
      metapar <- rbind(metapar[-1,], newrow)  
    } else {
      metapar <- rbind(metapar, newrow)  
    }
    if(idx%%10 == 0){
      writeLines(paste0(
        "*******************************************\n"
        , 'Fold = ',g, "\n"
        , "best error = ",besterr, "\n"
        , "SRbest = ", SRoos[g], "\n"
        , "Tries since best = ", idx, " \n"
        , "Number of runs = ", Nruns, "\n"
        , "Best config = "
      ))   
      print(besthp)
    }
  }
}

#######################
# oob error
fi <- list.files(path = outdir, pattern = "v8", full.names = TRUE)

# Make a foldkey
foldkey <- fi %>% strsplit("dry") %>% sapply(function(x){unlist(x)[2]}) %>% 
  gsub(pattern = "_", replacement = "") %>% 
  gsub(pattern = ".Rda", replacement = "") %>% 
  as.numeric

mergetemplate <- dat[,c("fips", "year")]
mergetemplate <- mergetemplate[with(mergetemplate, order(fips, year)),]

pnnoob_nopar <- foreach(i = 1:length(fi), .combine = cbind) %dopar% {
  set.seed(foldkey[i], kind = "L'Ecuyer-CMRG")
  samp <- sample(unique(dat$year), replace = TRUE)
  bsamp <- foreach(y = samp, .combine = c) %dopar% {which(dat$year == y)}
  oosamp <- unique(dat$year)[unique(dat$year) %ni% samp]
  PC <- readRDS(paste0(outdir,"/PCA_", foldkey[i],".Rds"))
  Xpc <- ppc(X[bsamp,], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
  Xtest <- ppc(X[dat$year %in% oosamp  & dat$fips %in% dat$fips[dat$year %in% samp],], PC, sum(cumsum((PC$sdev^2)/sum(PC$sdev^2))<.95))
  obj <- load_obj(fi[i])
  parlist <- obj$parlist
  pnn <- panelNNET(y = dat$yield[bsamp],
                   X = Xpc,
                   hidden_units = rep(100, 10),
                   fe_var = NULL,
                   maxit = 0,
                   time_var = dat$year[bsamp],
                   param = NULL,
                   verbose = T,
                   activation = 'lrelu',
                   parlist = parlist)
  newfips <- dat$fips[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]]
  newyear <- dat$year[dat$year %in% oosamp & dat$fips %in% dat$fips[dat$year %in% samp]]
  predoos <- predict(pnn, newX = as.matrix(Xtest),
                     new.param = NULL,
                     fe.newX = NULL)
  tm <- data.frame(fips = newfips, year = newyear, pred = predoos)
  with_holes <- full_join(mergetemplate, tm)
  with_holes <- with_holes[with(with_holes, order(fips, year)),]
  return(with_holes$pred)
}
save(pnnoob_nopar, file = paste0(outdir,"/pnnoob_nopar_GCB.Rda"))

MSE_result <- function(weights){
  pnnoob_pred <- apply(pnnoob_nopar, 1, weighted.mean, w = weights, na.rm = T)
  mergetemplate$pnnoob <- pnnoob_pred
  pnn_result <- full_join(dat[,c("fips", "year", "yield")], mergetemplate)
  with(pnn_result, mse(yield, pnnoob))
}
MSE_result(rep(1, 96))

