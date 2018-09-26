



rm(list=ls())
gc()
gc()
"%ni%" <- Negate("%in%")
library(foreach)


library(MASS)
library(panelNNET)
library(foreach)
library(doParallel)
library(rgdal)
library(rgeos)
library(doBy)

AWS <- grepl('ubuntu', getwd())
desktop <- grepl(':', getwd())
laptop <- grepl('/home/andrew', getwd())

if(laptop){
  path <- '/home/andrew/Dropbox/USDA/ARC/data/'
  setwd(path)
}

if (desktop){
  gispath <- paste0(getwd(),"/GISdata")
} else {
  gispath <- paste0(getwd(),"/shapefiles")
}


#######################
#Build the dataset
# corn

# load the NASS data on irrigated an non-irrigated yield
corn_irrdry <- rbind(read.csv("corn_NASS_irrdry_7980.csv"), read.csv("corn_NASS_irrdry.csv"))
# lose redundant colums
colnames(corn_irrdry) <- tolower(colnames(corn_irrdry))
corn_irrdry <- corn_irrdry[,c('year', "state.ansi", "county.ansi", "data.item", "value")]
# lose rows without fips codes
corn_irrdry <- corn_irrdry[!is.na(corn_irrdry$county.ansi),]
# make fips code
cofips <- as.character(corn_irrdry$county.ansi)
cofips[nchar(cofips) == 1] <- paste0("0", cofips[nchar(cofips) == 1])
cofips[nchar(cofips) == 2] <- paste0("0", cofips[nchar(cofips) == 2])
corn_irrdry$fips <- as.numeric(paste0(corn_irrdry$state.ansi,cofips))
corn_irrdry$county.ansi <- corn_irrdry$state.ansi <- NULL
# break out irr and dry
corn_irr <- subset(corn_irrdry, !grepl("NON", data.item))
corn_dry <- subset(corn_irrdry, grepl("NON", data.item))
rm(corn_irrdry)
corn_dry$data.item <- corn_irr$data.item <- NULL
colnames(corn_dry)[2] <- colnames(corn_irr)[2] <- "yield"

# load the aggregated nass data
corn_all <- rbind(read.csv("corn_NASS_all_1979_1980.csv"),
                  read.csv("corn_NASS_all_1981_1999.csv"),
                  read.csv("corn_NASS_all_2000_2016.csv"))
colnames(corn_all) <- tolower(colnames(corn_all))
corn_all <- corn_all[,c('year', "state.ansi", "county.ansi", "data.item", "value")]
corn_all <- corn_all[!is.na(corn_all$county.ansi),]
# make fips code
cofips <- as.character(corn_all$county.ansi)
cofips[nchar(cofips) == 1] <- paste0("0", cofips[nchar(cofips) == 1])
cofips[nchar(cofips) == 2] <- paste0("0", cofips[nchar(cofips) == 2])
corn_all$fips <- as.numeric(paste0(corn_all$state.ansi,cofips))
corn_all$county.ansi <- corn_all$state.ansi <- NULL
corn_all$data.item <- NULL
colnames(corn_all)[2] <- "yield"
# lose data from counties that I've got in the irr/dry data
corn_all <- subset(corn_all, fips %ni% unique(c(corn_dry$fips, corn_irr$fips)))

# load Nass irrigation data
corn_irrdat <- read.csv("corn_NASS_irrdry_area.csv")
# lose redundant colums
colnames(corn_irrdat) <- tolower(colnames(corn_irrdat))
corn_irrdat <- corn_irrdat[,c('year', "state.ansi", "county.ansi", "data.item", "value")]
# make fips code
cofips <- as.character(corn_irrdat$county.ansi)
cofips[nchar(cofips) == 1] <- paste0("0", cofips[nchar(cofips) == 1])
cofips[nchar(cofips) == 2] <- paste0("0", cofips[nchar(cofips) == 2])
corn_irrdat$fips <- as.numeric(paste0(corn_irrdat$state.ansi,cofips))
corn_irrdat$county.ansi <- corn_irrdat$state.ansi <- NULL
# fix commas
corn_irrdat$value <- gsub(",", "", as.character(corn_irrdat$value))
# set the D's to NA
corn_irrdat$value[grepl("D", corn_irrdat$value)] <- NA
corn_irrdat$value <- as.numeric(corn_irrdat$value)
colnames(corn_irrdat)[3] <- "acres"
corn_irr_acres <- subset(corn_irrdat, grepl("IRRIGATED", data.item))
corn_all_acres <- subset(corn_irrdat, !grepl("IRRIGATED", data.item))
corn_all_acres$data.item <- corn_irr_acres$data.item <- NULL
colnames(corn_all_acres)[2] <- "all_acres"
colnames(corn_irr_acres)[2] <- "irr_acres"
corn_irrdat <- merge(corn_all_acres, corn_irr_acres, all = TRUE)
# lose counties with no data
corn_irrdat <- corn_irrdat[!is.na(corn_irrdat$all_acres),]
# lose counties with less than 500 acres in corn
corn_irrdat <- corn_irrdat[corn_irrdat$all_acres >500,]
#summarize n_acres by year.  Need this for bubble plots
corn_irrdat_tmp <- corn_irrdat
corn_irrdat_tmp$irr_acres[is.na(corn_irrdat_tmp$irr_acres)] <- 0
corn_acreage <- summaryBy(all_acres + irr_acres ~fips, data = corn_irrdat_tmp, FUN = mean, 
                           na.rm = TRUE, keep.names = TRUE)

# lose counties that I've got irr/dry yield data for
corn_irrdat <- corn_irrdat[corn_irrdat$fips %in% corn_all$fips,]
# now that this is done, set remaining irr NA's to zero
corn_irrdat$irr_acres[is.na(corn_irrdat$irr_acres)] <- 0
corn_irrdat$prop_irr <- with(corn_irrdat, irr_acres/all_acres)
#summarize prop_irriaged by year
corn_irrdat <- summaryBy(prop_irr~fips, data = corn_irrdat, FUN = mean, 
                         na.rm = TRUE, keep.names = TRUE)

# merge the data from *all* into irr and dry, with the prop_irr variable
corn_dry$prop_irr <- 0
corn_irr$prop_irr <- 1
corn_all <- merge(corn_all, corn_irrdat[,c("fips", "prop_irr")],
                  by = "fips")
corn_all <- corn_all[,c(2,3,1,4)] #reorder columns
corn_irr <- rbind(corn_irr, subset(corn_all, prop_irr>=.5))
corn_dry <- rbind(corn_dry, subset(corn_all, prop_irr<=.5))

# save the data
write.csv(corn_irr, file = "corn_NASS_yields_irr_augmented.csv")
write.csv(corn_dry, file = "corn_NASS_yields_dry_augmented.csv")


