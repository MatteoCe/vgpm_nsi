#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

args

if(args[1] == 1){

# No specific packages are necessary to calculate the VGPM. Here, the package geosphere was used in previous version of the
# script to calculate the daylength of specific days, later removed to calculate it using base R

# define max number of cores for parallelization. Here it is set to half of the available, but can me modified
n.core <- detectCores()/2

#library(geosphere)
library(tidyverse)
library(parallel)

filename <- args[2]

# import the 8day period tsv file with all the chl, par and sst values for each lat/lon combination
file_read <- read_delim(paste("processed/",filename,"_tot_date_clean_head", sep=''), col_names = TRUE, delim='\t', show_col_types = FALSE)

# define functions for calculation of Net Primary Production following the VGPM algorithm. These functions will allow to
# calculate the euphotic zone and the "maximum daily net primary production found within a given water column and expressed
# in units of mg carbon fixed per mg chlorophyll per hour". A thorough description of the VGPM model and its calculation
# can be found at http://sites.science.oregonstate.edu/ocean.productivity/vgpm.model.php and
# http://sites.science.oregonstate.edu/ocean.productivity/vgpm.code.php

euphotic <- function(x) {
	if (x <  1) {
        chl_tot <- 38 * x^0.425
        z_eu <- 200 * chl_tot^-.293
        if (z_eu <= 102) {
            y <- 568.2 * chl_tot^-.746
        } else {
            y <- z_eu
        }
    } else {
        chl_tot <- 40.2 * x^0.507
        z_eu <- 200 * chl_tot^-.293
        if (z_eu <= 102) {
            y <- 568.2 * chl_tot^-.746
        } else {
            y <- z_eu
        }
	}
}

Pb_opt <- function(x) {
	if (x < -10) {
        y <- 0
    } else if (x < -1) {
        y <- 1.13
	} else if (x > 28.5) {
		y <- 4
	} else {
		y <- 1.2956 +
		2.749e-1*x +
		6.17e-2*x^2 -
		2.05e-2*x^3 +
		2.462e-3*x^4 -
		1.348e-4*x^5 +
		3.4132e-6*x^6 -
		3.27e-8*x^7
	}
}

#length_of_period <- function(x) {
#	y <- geosphere::daylength(x, as.Date(t[[1]]))
#}

# Create a function to calculate the daylength of a particular date. This is necessary for the calculation of the Net
# Primary Productivity.
length_of_period_mod <- function(x) {
	pmin(pmax((sin(0.8333 * pi/180) + sin(x * pi/180) * sin(P))/(cos(x *pi/180) * cos(P)), -1), 1)
}

# Now, the grid value will have to be converted in lat/lon decimal format. Thi is perfomed following the conversion
# instruction at http://orca.science.oregonstate.edu/faq01.php. In the specific:

# The north west corner of the start of the hdf file is at +90 lat, -180 lon.
#
# To obtain the location of the center of any pixel:
# take the number of rows and columns you are away from the NW corner,
# multiply by the grid spacing to get the change in latitude and longitude,
# subtract the change in latitude from +90 lat,
# add the change in longitude to -180 lon;
# shift the latitude down (subtract) by 1/2 of a grid spacing
# and shift the longitude over (add) by 1/2 of a grid spacing 

file_read$longitude <- (-180+((file_read$longitude-0.5)*1/6))+(1/2*1/6)
file_read$latitude <- (90-((file_read$latitude-0.5)*1/6))-(1/2*1/6)

# The entire analyses can be performed on a geographical subset of the world, indicating maximum and minimum latitudes and
# longitudes- This is performed relying on the informations gathered from a before-hand created tsv file, named region_coord,
# that indicates the bounding box values. Ifor the examples of the script, region_coord has values indicating the entire
# world. Multiple regions can be indicated in the region_coord file, one for each line and named in the first column.

region_coord <- read.csv("data/region_coord", header=T, sep="\t")
reg.group = list()

# Inside each chosen bounding box, which is selection of the original imported file, the functions "euphotic", "Pb_opt" and
# the calculation fo the irradiance function (also called the "Volume function") are applied to chl, sst and par values
# respectively. To speed up the analyses, this is performed in parallel using mclapply from the parallel package.

# For speed purposes, the script has been modified in order to use only the first date of the 8day time span
# as date for the calculation of day length for each point, instead of doing a mean of the daylength of each date in
# the 8day.

for (reg in levels(factor(region_coord$region))) {
    
    region_coord_sel <- region_coord[region_coord$region==reg,]
    
    minlon <- region_coord_sel$minlon
    maxlon <- region_coord_sel$maxlon
    minlat <- region_coord_sel$minlat
    maxlat <- region_coord_sel$maxlat

    tab1 <- file_read[(file_read$longitude > minlon & file_read$longitude < maxlon),]
    tab <- tab1[(tab1$latitude > minlat & tab1$latitude < maxlat),]
    
	if (nrow(tab) == 0) {
		
		next
		
		} else {
			
			}
	
	tab$z_eu <- unlist(mclapply(tab$chl, euphotic, mc.cores=n.core))
	tab$pb_opt <- unlist(mclapply(tab$sst, Pb_opt, mc.cores=n.core))
	tab$irrFunc <- 0.66125 * tab$par / (tab$par + 4.1)

    # based on the type of files, "8day" or "month", average daylength for the specific timeframe is calculated
	# however, as mentioned before, this is currently inactive for speed purposes.
    if (args[3] == "month") {
                
        t0 <- as.POSIXct(paste(unique(tab$date),"-01 12:00:00",sep=''), tz="UTC")
                    
        nummon <- as.numeric(lubridate::days_in_month(t0))
                    
        t <- seq.POSIXt(t0, by="1 day", length.out=1*nummon)
                    
    } else if (args[3] == "8day") {
                    
        t0 <- unique(tab$date)
                    
        date1 <- as.POSIXct(paste(stringr::str_replace(t0, "_.*", ""), "12:00:00",sep=''), tz="UTC")
                    
        #date2 <- as.POSIXct(paste(stringr::str_replace(t0, ".*_", ""), "12:00:00",sep=''), tz="UTC")
        #            
        #daynum <- lubridate::interval(date1, date2)/lubridate::days(1)
        #            
        #t <- seq.POSIXt(date1, by="1 day", length.out=1*daynum)
                    
    }
	
	# this variable "P" is used by the "length_of_period_mod" function to calculate the daylength. It is based on the date
	# in each loop, and will be merged in the function in later versions of the script
	P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (as.numeric(format(as.Date(date1), "%j")) - 186)))))
	
	# Apply the daylength calculations to each latitude values
	tab$a <- unlist(mclapply(as.numeric(tab$latitude), length_of_period_mod, mc.cores=n.core))
	tab$dayL <- 24 - (24/pi) * acos(tab$a)
    
	# Finally, calculate the VGPM model, using all the necessary data for each lat/lon point
    tab$npp <- tab$pb_opt * tab$chl * tab$dayL * tab$irrFunc * tab$z_eu
	# simplify table keeping only lat,lon and npp
	tab.clean <- tab[,c("longitude", "latitude", "npp")]
	
	# store data for each region in a list
    reg.group[[reg]] <- tab.clean

}

# rowwise-bind all regions.
reg.tab.call = do.call(rbind, reg.group)

# finally, write the entire table as a single tsv file, including all regions.
write.table(reg.tab.call, file=paste("processed/ext_",filename,".tab", sep=''), row.names=FALSE, sep='\t', quote=FALSE)

} else if (args[1] == 2) {

library(tidyverse)
library(parallel)

# define max number of cores for parallelization. In this case, they are all except one, but can be reduced
n.core <- detectCores()-1

# As in the first part of the script, the argument passed to the script call defined the type of file we're dealing with,
# monthly of 8day data. This also checks how many observations correspond to the "half-time" concept. If monthly data, then
# the half-time will be 6 moths, whereas for the 8day data it will be approximately 23 ("weeks", 365 days divided by 8 days,
# divided by two)
if (args[2] == "month") {
	t0 <- 6 # half of possible number of observetions in a year (12 months / 2)

} else if (args[2] == "8day") {
	t0 <- 23 # half of possible number of observetions in a year (46 "weeks" / 2)

}

# get a list with each chosen year as a vector from the "filelist" file in the processed/ subdirectory.
# This is done because the seasonality index is an annual index, thus the index will be calculate for each year
years <- unique(stringr::str_replace(scan(file="processed/filelist"), "...$", ""))

# create list of colnames for the final table based on the number of years and the two chosen seasonality indices.
# Previous versions of the script also calculated the seasonality index by Lutz et al. 2007, although, again for speed,
# that part is removed from the analyses, and the l.inspp column will be left empty and removed.
years_col <- unlist(lapply(1:length(years), function(x) {
	prim <- stringr::str_replace(years[x], "$", "_b.isnpp")
	sec <- stringr::str_replace(years[x], "$", "_l.isnpp")
	y <- c(prim, sec)
	}))

# create empty list to store the tables of each year
final.list = list()

# initiate the loop for each year

for (y in years) {
	
	# get a list of files belonging to the that year
	files <- list.files(path = "processed/", pattern = y) 
	
	# load them singularly into a list with all of them
	files_list <- lapply(files, function(f) {x <- read_delim(file=paste("processed/",f,sep=''), col_names = TRUE, delim='\t', show_col_types = FALSE)})
	
	# substitute empty tibbles with NULL
	files_list.null <- lapply(files_list, function(x) { if (nrow(x) == 0) {} else {x}})
	
	# remove empty observations
	files_list <- files_list.null[!vapply(files_list.null, is.null, logical(1))]
	
	# create another list of files, assuring they include only lat/lon and npp values
	simp <- lapply(files_list, function(x) {x[,c("longitude", "latitude", "npp")]})
	
	# merge all tables by longitude and latitude values, thus having the npp values of each 8day of that year for all
	# combinations of lat/lon
	simp2 <- purrr::reduce(simp, dplyr::full_join, by = c("longitude", "latitude"))
	
	# create a new field, t, which counts the number of observations that sum up, after ordering by decreasing
	# values of npp, to half of the production of that year
	simp2$t <- apply(simp2, 1, function(x) {length(which(cumsum(sort(x, decreasing=TRUE)) <= sum(sort(x, decreasing=TRUE))/2))})
	
	# calcaulate the Nornalized Seasonality Index (NSI) intoruced by Brown et al. 2014
	simp2$b.isnpp <- (t0-simp2$t)/t0
	
	# simplify table by including only lat/lon and the NSI for the year of the loop
	simp3 <- simp2[,c("longitude", "latitude", "b.isnpp")]
	
	# add the year annotation to the column name of the NSI, this will discriminate different columns that correspond to the
	# NSI for each year investigated
	colnames(simp3) <- c("longitude", "latitude", paste(y,"b.isnpp", sep='_'))
	
	# add the table to the empty list list
	final.list[[y]] <- simp3
	
	# save memory by removing all newly created objects in each loop that are not useful (aòthoug not that necessary as
	# they would be overwritten at each iteration)
	rm(files, files_list, files_list.null, simp, simp2, simp3)
	gc()
	
}

# merge all data for each lat/lon
final.tab <- purrr::reduce(final.list, dplyr::full_join, by = c("longitude", "latitude"))

# create a simplyfied version of the previously created table including only the NSI values for each year
b.isnpp.tab <- tibble(cbind(final.tab[,1:2],select(final.tab, matches("*b.isnpp"))))

# calculate the average values of each row, thus the average NSI for that years period
b.isnpp.tab$b.isnpp <- rowMeans(select(b.isnpp.tab, matches("*b.isnpp")), na.rm=TRUE)

# some lat/lon combinations might not enough obervations to calculate the index, thus will be empty and need to be removed
b.isnpp.tab.clean <- b.isnpp.tab[!is.na(b.isnpp.tab$b.isnpp), c("longitude", "latitude", "b.isnpp")]

# write table with values for each year
write.table(final.tab, file="processed/final/final.tab", row.names=FALSE, sep ='\t', quote = FALSE)

# write table with final valu of the average NSI
write.table(b.isnpp.tab.clean, file="processed/final/b_isnpp.tab", row.names=FALSE, sep ='\t', quote = FALSE)

# call shell to remove all remaining, useless files
system("rm -f processed/*tab", intern=FALSE)

}
	
