

rm(list=ls())
gc()
gc()
"%ni%" <- Negate("%in%")


library(MASS)
library(foreach)
library(doParallel)
registerDoParallel(detectCores())
library(rgdal)
library(rgeos)
library(doBy)
library(dplyr)

AWS <- grepl('ubuntu', getwd())
desktop <- grepl(':', getwd())
laptop <- grepl('/home/andrew', getwd())
if(AWS){
  path <- '/home/ubuntu/'
  setwd(path)
  macadir <- "/home/ubuntu/impacts_scale/data/maca_county"
}
if(desktop){
  path <- ''
  setwd(path)
}
if(laptop){
  path <- '/home/andrew/Dropbox/USDA/ARC/data'
  setwd(path)
  macadir <- paste0(getwd(), "/maca_county")
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


nass_yields_files <- list.files(pattern = "augmented.csv")
for (cr in c("corn", "soy", "wheat")){
  dry <- read.csv(nass_yields_files[grepl(cr, nass_yields_files) & grepl("dry", nass_yields_files)])
  irr <- read.csv(nass_yields_files[grepl(cr, nass_yields_files) & grepl("irr", nass_yields_files)])
  irr$X <- dry$X <- NULL
  #add tiny noise to exact 50% entries
  set.seed(1)
  irr$prop_irr[irr$prop_irr == .5] <- .500001
  dry$prop_irr[dry$prop_irr == .5] <- .500001
  irr <- subset(irr, prop_irr > .5)
  dry <- subset(dry, prop_irr < .5)
  overlaps <- intersect(dry$fips, irr$fips)
  dat <- rbind(dry[dry$fips %in% overlaps,],irr[irr$fips %in% overlaps,])
  dat$irrdry_reported <- TRUE
  nrdat <- rbind(dry[dry$fips %ni% overlaps,],irr[irr$fips %ni% overlaps,])
  nrdat$irrdry_reported <- FALSE
  dat <- rbind(dat, nrdat)
  rm(nrdat)
  # merge on weather data
  fips_shp <- readOGR(dsn = gispath, layer = "ers_fipscrop")
  fips_shp <- spTransform(fips_shp, CRS("+proj=longlat +datum=WGS84"))
  centroids <- gCentroid(fips_shp, byid=TRUE)
  fips_shp$lat <- centroids@coords[,"y"]
  fips_shp$lon <- centroids@coords[,"x"]
  fips_shp$FIPS <- as.character(fips_shp$FIPS)
  fips_shp$FIPS[substr(fips_shp$FIPS,1,1) == "0"] <- substr(fips_shp$FIPS[substr(fips_shp$FIPS,1,1) == "0"],2,nchar(fips_shp$FIPS))
  fips_shp$FIPS <- as.numeric(fips_shp$FIPS)
  fips_shp@data <- renameCol(fips_shp@data, "FIPS", "fips")
  tm <- with(fips_shp@data, data.frame(fips, lat, lon))
  dat <- merge(dat, tm)
  # lose singletons
  singletons <- names(which(table(dat$fips) == 1))
  dat <- dat[dat$fips %ni% singletons,]
  # merge in soil data
  impsoil <- read.csv("FIPSSSURGONumeric_imputed.csv")[,-1] # lose the X
  dat <- left_join(dat, impsoil)
  colnames(dat)[colnames(dat) %in% colnames(impsoil)][-1] <- paste0("soil_",colnames(dat)[colnames(dat) %in% colnames(impsoil)][-1])
  # merge on the county-level MACA data
  fi <- list.files(path = macadir, full.names = T)
  maca <- NULL
  for (i in fi){
    macai <- readRDS(i)
    colnames(macai)[colnames(macai) == "statyear"] <- 'year'
    inp <- inner_join(dat, macai, by = c('year', 'fips'))
    maca <- rbind(maca, inp)
    print(i)
  }
  # kelvin to celsius
  maca[,grepl("maxat|minat", colnames(maca))] <- maca[,grepl("maxat|minat", colnames(maca))] - 272.15
  
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
  # specify the days in a way that makes it simpler to query them later
  alpha_part <- gsub("\\d", "", colnames(maca))
  num_part <- as.numeric(gsub("\\D", "", colnames(maca)))
  alpha_part[!grepl(regex, colnames(maca))] <- NA
  num_part[!grepl(regex, colnames(maca))] <- NA
  colnames(maca)[!is.na(num_part)] <- paste0(alpha_part[!is.na(num_part)], 
                                             "_jday_",
                                             num_part[!is.na(num_part)]
  )
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
  # Make a state variable
  maca$state <- substr(maca$fips, 1, nchar(maca$fips)-3)
  maca$state <- as.factor(maca$state)
  maca$fips <- as.factor(maca$fips)
  # save
  saveRDS(maca, file = paste0("panel_", cr, ".Rds"))
}


