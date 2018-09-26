

rm(list=ls())
gc()
gc()
"%ni%" <- Negate("%in%")
hx <- function(x, a=6, b = 6){x[1:a, 1:b]}

library(MASS)
library(foreach)
library(doParallel)
registerDoParallel(detectCores())
library(rgdal)
library(rpart)
library(rgeos)
library(doBy)
library(dplyr)
library(sp)
library(doBy)

AWS <- grepl('ubuntu', getwd())
desktop <- grepl(':', getwd())
laptop <- grepl('/home/andrew', getwd())
if(AWS){
  path <- '/home/ubuntu/'
  setwd(path)
}
if(desktop){
  path <- ''
  setwd(path)
}
if(laptop){
  path <- '/home/andrew/Dropbox/USDA/ARC/data'
  setwd(path)
}

if(desktop){
  cores <- detectCores()
  cl <- makeCluster(cores)
  print(cores)
  registerDoParallel(cl)
} else {
  registerDoParallel(detectCores())
  print(detectCores())
}

gispath <- paste0(getwd(),"/shapefiles")
macadir <- "/home/ubuntu/MACA/rdata/"
#S&R variables
days.in.range <- function( t0, t1 , tMin, tMax, noGrids )  {
  n <-  length(tMin)
  t0 <-  rep(t0, n)
  t1 <-  rep(t1, n)
  t0[t0 < tMin] <-  tMin[t0 < tMin]
  t1[t1 > tMax] <-  tMax[t1 > tMax]
  u <- function(z, ind) (z[ind] - tMin[ind])/(tMax[ind] - tMin[ind])  
  outside <-  t0 > tMax | t1 < tMin
  inside <-  !outside
  time.at.range <- ( 2/pi )*( asin(u(t1,inside)) - asin(u(t0,inside)) ) 
  return( sum(time.at.range)/noGrids ) 
}

DIR.wrapper <- function( t0, t1 , tMin, tMax, noGrids ) { #Wrapper function with error handling for NA's.  any year with a single NA outputs all NAs.  Some sort of imputation is feasible, as is dropping.  Both would introduce bias.
  out <- tryCatch(days.in.range(t0, t1 , tMin, tMax, noGrids),  error = function(e) e)
  if (inherits(out, 'error')){out <- NA}
  return(out)
}

# starter data frame
fips_shp <- readOGR(dsn = gispath, layer = "ers_fipscrop")
fips_shp <- spTransform(fips_shp, CRS("+proj=longlat +datum=WGS84"))
centroids <- gCentroid(fips_shp, byid=TRUE)
fips_shp$lat <- centroids@coords[,"y"]
fips_shp$lon <- centroids@coords[,"x"]
fips_shp$FIPS <- as.character(fips_shp$FIPS)
fips_shp$FIPS[substr(fips_shp$FIPS,1,1) == "0"] <- substr(fips_shp$FIPS[substr(fips_shp$FIPS,1,1) == "0"],2,nchar(fips_shp$FIPS))
fips_shp$FIPS <- as.numeric(fips_shp$FIPS)
fips_shp@data <- renameCol(fips_shp@data, "FIPS", "fips")
dat <- with(fips_shp@data, data.frame(fips, lat, lon))
impsoil <- read.csv("FIPSSSURGONumeric_imputed.csv")[,-1] # lose the X
dat <- inner_join(dat, impsoil)
colnames(dat)[colnames(dat) %in% colnames(impsoil)][-1] <- paste0("soil_",colnames(dat)[colnames(dat) %in% colnames(impsoil)][-1])

# parse filenames to get a list of models
fi <- list.files(path = macadir, pattern = "wide.rds")
filedf <- foreach(i = fi, .combine = rbind) %dopar% {
  strsplit(i, split = "_") %>% unlist
} %>% as.data.frame
colnames(filedf) <- c("term", "rcp", "year", "model", "x")
filedf$x <- NULL
filedf <- filter(filedf, term != "huss" & model != "NorESM1-M")

modlist <- unique(filedf$model) %>% as.character
varlist <- c("tasmax","tasmin","pr", "rhsmax", "rhsmin","rsds","uas","vas")
newvarlist <- c("maxat","minat","precip","maxrh","minrh", "sdsfia","uwind","vwind")
rcplist <- c("rcp45","rcp85")
for (m in modlist){
  for (r in rcplist){
    for (cr in c("corn")){#}, "soy", "wheat")){

      # merge on the county-level MACA data
      maca <- foreach (y = 2006:2099, .combine = rbind) %dopar%{
        macay <- foreach(v = varlist, .combine = cbind) %do% {
          fn <- paste(v, r, y, m, "wide.rds", sep = "_")
          vdat <- readRDS(paste0(macadir, fn))
          # rename columns
          tochange <- (colnames(vdat[3]) %>% strsplit(split = "_jday_") %>% unlist)[1]
          colnames(vdat) <- gsub(tochange, newvarlist[which(varlist == v)], colnames(vdat))
          # check and see if there are 366 days
          if (ncol(vdat) == 368){
            vdat <- vdat[,!grepl("_jday_60", colnames(vdat))]
            # rename the others
            colnames(vdat)[colnames(vdat) %in% paste0(newvarlist[which(varlist == v)],"_jday_", 61:366)] <- paste0(newvarlist[which(varlist == v)],"_jday_", 60:365)
          }
          # rename columns
          tochange <- (colnames(vdat[3]) %>% strsplit(split = "_jday_") %>% unlist)[1]
          colnames(vdat) <- gsub(tochange, newvarlist[which(varlist == v)], colnames(vdat))
          #subset counties
          vdat <- vdat[substr(vdat$fips, 1, 2) %in% c("17", "19", "27", "18", "39", "26", "21", "55", "29"),]
          vdat <- vdat[with(vdat, order(year, fips)),]
          if (v == varlist[1]){idx <- vdat[,c("year", "fips")]} 
          vdat <- vdat[colnames(vdat) %ni% c("year", "fips")]
          return(vdat)
        }
        return((data.frame(idx, macay)))
      }    

      # kelvin to celsius
      maca[,grepl("maxat|minat", colnames(maca))] <- maca[,grepl("maxat|minat", colnames(maca))] - 272.15
      # deal with wind
      wind <- sqrt(maca[,grep("vwind", colnames(maca))]^2 + maca[,grep("uwind", colnames(maca))]^2)
      colnames(wind) <- substring(colnames(wind), 2, nchar(colnames(wind)))
      maca <- maca[,!grepl("wind", colnames(maca))]
      maca <- data.frame(maca, wind)
      # lose non-growing-season
      if (cr == "corn" | cr == "soy"){
        colnums <- as.numeric(gsub("\\D", "", colnames(maca)))
        regex <- "minat|maxat|precip|maxrh|minrh|spech|sdsfia|wind"
        colnums[!grepl(regex, colnames(maca))] <- NA
        torm <- which((colnums < 60 | colnums >= 305) & !is.na(colnums))
        maca <- maca[,1:ncol(maca) %ni% torm]    
      } else if (cr == "wheat"){ #shift the year
        colnums <- as.numeric(gsub("\\D", "", colnames(maca)))
        regex <- "minat|maxat|precip|maxrh|minrh|spech|sdsfia|wind"
        colnums[!grepl(regex, colnames(maca))] <- NA
        prevyear <- maca[,(colnums >196 & !is.na(colnums))|colnames(maca) %in% c("year", "fips")]
        prevyear$key <- 1:nrow(prevyear)
        prevyear$year <- prevyear$year + 1 # counts for following year
        prevyear <- prevyear[prevyear$year <= 2016,]
        thisyear <- maca[,colnums <=196 | is.na(colnums)]
        thisyear$key <- 1:nrow(thisyear)
        maca <- inner_join(thisyear, prevyear, by = "key")
      }
      # derived variables
      tgrid <- -1:40
      noGrids <- sum(grepl('minat', colnames(maca)))
      SRmat <- foreach(j = 1:nrow(maca), .combine = rbind) %dopar% {
        tmin <- as.numeric(maca[j,grepl('minat', colnames(maca))])
        tmax <- as.numeric(maca[j,grepl('maxat', colnames(maca))])
        SRvars <- foreach(i = tgrid, .combine = c) %do% {
          t0 <- ifelse(i == min(tgrid), -100, i)
          t1 <- ifelse(i == max(tgrid), 100, i+1)
          DIR.wrapper(i,i+1,tmin, tmax, noGrids)    
        }
        return(SRvars)
      }
      SRmat <- as.data.frame(SRmat)
      colnames(SRmat) <- paste0('SR',gsub("-","neg",tgrid))
      maca <- cbind(maca, SRmat)
      maca$fips <- maca$fips %>% as.character %>% as.numeric
      # join to yield and soil data
      maca <- inner_join(dat, maca)
      # Make a state variable
      maca$state <- substr(maca$fips, 1, nchar(maca$fips)-3)
      maca$state <- as.factor(maca$state)
      maca$fips <- as.factor(maca$fips)
      # save
      saveRDS(maca, file = paste0("panel_GCB_", r,"_", cr,"_", m, ".Rds"))
    }
  }
}


