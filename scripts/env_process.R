#!/usr/bin/Rscript

library(magrittr)
library(purrr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

args

# define functions for calculation of Net Primary Production following 
# the VGPM algorithm from http://orca.science.oregonstate.edu/faq01.php and
# http://sites.science.oregonstate.edu/ocean.productivity/vgpm.code.php

# calculation of euphotic zone depth
euphotic <- function(x) {

  if (is.na(x)) {

    NA

  } else if (x < 1) {

    chl_tot <- 38 * x^0.425

    z_eu <- 200 * chl_tot^-0.293

    if (z_eu <= 102) {

      568.2 * chl_tot^-0.746

    } else {

      z_eu

    }

  } else {

    chl_tot <- 40.2 * x^0.507

    z_eu <- 200 * chl_tot^-0.293

    if (z_eu <= 102) {

      568.2 * chl_tot^-0.746

    } else {

      z_eu

    }

  }

}

# maximum photosynthetic efficiency
pb_opt <- function(x) {

  if (is.na(x)) {

    NA

  } else if (x < -10) {

    0

  } else if (x < -1) {

    1.13

  } else if (x > 28.5) {

    4

  } else {

    1.2956 +
    2.749e-1*x +
    6.17e-2*x^2 -
    2.05e-2*x^3 +
    2.462e-3*x^4 -
    1.348e-4*x^5 +
    3.4132e-6*x^6 -
    3.27e-8*x^7

  }

}

# irradiance function based on photosynthetically active radiation values
irradiance <- function(x) {

  if (is.na(x)) {

    NA

  } else {

    0.66125 * x / (x + 4.1)

  }

}

# calculation of daylength
LatToDayLength <- function(lat, yDay) {
  # get lat in radians
  gamma <- lat / 180.0 * pi
  
  # convert date into an angle
  psi <- yDay / 365.0 * 2.0 * pi
  
  # calc solar declination
  # Kirk page 35
  solarDec <- (0.39637
               - 22.9133 * cos(psi)
               + 4.02543 * sin(psi)
               - 0.38720 * cos(2 * psi)
               + 0.05200 * sin(2 * psi)) * pi / 180.0
  
  r <- -tan(gamma) * tan(solarDec)
  if (r <= -1) {
    return(24.0)
  } else if (abs(r) < 1) {
    return(24.0 * acos(r) / pi)
  } else {
    return(0)
  }
}

# irradiation function
irrFunc <- function(par) {
  return(0.66125 * par / (par + 4.1))
}

# Normalized seasonality index 

nsi_calc <- function(npp) {

  if (all(is.na(npp))) {

    return(NA)

  } else {

    prod <- c(npp, recursive = TRUE, use.names = FALSE) %>% sort(., decreasing = TRUE)

    t <- length(which(cumsum(prod) <= sum(prod) / 2))

    (23 - t) / 23

  }

}

# Based on the first argument provided to the Rscript command, different 
# calculations will be performed on data.

# Daily copernicus data corresponding to a 8-day period will be checked, null 
# points "-32767" will be turned to NA, coordinates transitioned to pixel points
# and salinity values turned to a mean for this period. This is done to resemble
# the original time period (8-day) provided by MODIS.

# The final weekly data will be appended to a "year" file, which will include 
# all 8-day periods for a specific year. 

if (args[1] == 1) {

  name <- args[2]

  tab <- readr::read_delim(file = paste0("tmp/processed/week/", name), 
                           delim = "\t", 
                           na = "-32767", 
                           col_names = FALSE,
                           show_col_types = FALSE)

  tab[, 1] <- round((180 + tab[, 1]) * 12, digits = 0)

  tab[, 2] <- round((90 - tab[, 2]) * 12, digits = 0)

  tab <- tab %>% dplyr::select(-c(1,2)) %>% 
                 dplyr::mutate_all(., function(x) {

                    (x - 0.001525925472378731) * 0.00152592547237873}) %>%

                 dplyr::mutate(m = rowMeans(., na.rm = TRUE)) %>%
                 dplyr::select(m) %>%
                 dplyr::bind_cols(tab[,c(1,2)], .)

  tab[is.nan(tab$m),]$m <- NA

  tab$m <- round(tab$m, digits = 3)

  tab <- tab %>% dplyr::rename(lon = X1, lat = X2, {{name}} := m)

  new_name <- stringr::str_remove_all(name, "\\-.*")

  if (!file.exists(paste0("tmp/processed/year/", new_name))) {

    readr::write_delim(tab,
                       file = paste0("tmp/processed/year/", new_name),
                       delim = "\t",
                       na = "NA",
                       escape = "none")

  } else {

    final_tab <- readr::read_delim(file = paste0("tmp/processed/year/", new_name), 
                                   delim = "\t",
                                   na = "NA",
                                   col_names = TRUE,
                                   show_col_types = FALSE)

    final_tab <- dplyr::bind_cols(final_tab, tab[, 3])

    rm(tab)

    readr::write_delim(final_tab, 
                       file = paste0("tmp/processed/year/", new_name),
                       delim = "\t",
                       na = "NA",
                       escape = "none",
                       append = FALSE)

    rm(final_tab)

  }

# Range and mean values will be calculated for salinity copernicus data for each
# full year, including all 8.day periods for that year.

} else if (args[1] == 2) {

  previous_year <- args[2]

  tab <- readr::read_delim(file = paste0("tmp/processed/year/", previous_year), 
                           delim = "\t",
                           na = "NA",
                           col_names = TRUE,
                           show_col_types = FALSE)

  tab <- tab %>% dplyr::select(-c(1,2)) %>% 
                 dplyr::mutate(m = rowMeans(., na.rm = TRUE), 
                               max = do.call(pmax, c(., na.rm = TRUE)),
                               min = do.call(pmin, c(., na.rm = TRUE))) %>%
                 dplyr::mutate(r = max - min) %>%
                 dplyr::select(m, r) %>%
                 dplyr::bind_cols(tab[,c(1,2)], .)

  tab[is.nan(tab$m),]$m <- NA

  tab$m <- round(tab$m, digits = 3)
  tab$r <- round(tab$r, digits = 3)

  m <- paste0(previous_year, "_m")
  r <- paste0(previous_year, "_r")

  tab <- tab %>% dplyr::rename({{m}} := m, {{r}} := r)

  if (!file.exists("variables/sss_final")) {

    readr::write_delim(tab, 
                       file = "variables/sss_final",
                       delim = "\t",
                       na = "NA",
                       escape = "none",
                       append = FALSE)

  } else {

    final_tab <- readr::read_delim(file = "variables/sss_final", 
                                   delim = "\t",
                                   col_names = TRUE,
                                   show_col_types = FALSE)

    final_tab <- dplyr::bind_cols(final_tab, tab[, c(3, 4)])

    rm(tab, m, r)

    readr::write_delim(final_tab, 
                       file = "variables/sss_final",
                       delim = "\t",
                       na = "NA",
                       escape = "none",
                       append = FALSE)

  }

# The final option used by the env_download shell script, will calculate 
# long-term averages of all statistics for the entire year.

} else if (args[1] == 3) {

  year_range <- args[2]

  range_pattern <- seq(stringr::str_remove(year_range, "\\_.*"), 
                       stringr::str_remove(year_range, ".*\\_"), 1) %>% 
                   paste(., collapse = "|")

  parameters <- list.files("variables/") %>% stringr::str_remove("_.*")

  for (par in parameters) {

    tab <- readr::read_delim(file = paste0("variables/", par, "_final"), 
                             delim = "\t",
                             col_names = TRUE,
                             show_col_types = FALSE)

    if (par == "npp") {

      tab <- tab %>% dplyr::select(., dplyr::matches(range_pattern)) %>%
                     dplyr::mutate(mean = rowMeans(dplyr::select(., dplyr::matches("_m")), 
                                                   na.rm = TRUE)) %>%
                     dplyr::mutate(total = rowMeans(dplyr::select(., dplyr::matches("_tot")), 
                                                    na.rm = TRUE)) %>%
                     dplyr::mutate(nsi = rowMeans(dplyr::select(., dplyr::matches("_nsi")), 
                                                   na.rm = TRUE)) %>%
                     dplyr::select(mean, total, nsi) %>%
                     dplyr::bind_cols(tab[,c(1, 2)], .)

      tab$mean <- round(tab$mean, digits = 4)
      tab$nsi <- round(tab$nsi, digits = 4)
      tab$total <- round(tab$total, digits = 4)

    } else if (any(par %in% c("chl", "par", "k490"))) {

      tab <- tab %>% dplyr::select(., dplyr::matches(range_pattern)) %>%
                     dplyr::mutate(mean = rowMeans(dplyr::select(., dplyr::matches("_m")), 
                                                   na.rm = TRUE)) %>%
                     dplyr::select(mean) %>%
                     dplyr::bind_cols(tab[,c(1, 2)], .)

      tab$mean <- round(tab$mean, digits = 4)

    } else if (any(par %in% c("sss", "sst"))) {

      tab <- tab %>% dplyr::select(., dplyr::matches(range_pattern)) %>%
                     dplyr::mutate(mean = rowMeans(dplyr::select(., dplyr::matches("_m")), 
                                                   na.rm = TRUE)) %>%
                     dplyr::mutate(range = rowMeans(dplyr::select(., dplyr::matches("_r")), 
                                                    na.rm = TRUE)) %>%
                     dplyr::select(mean, range) %>%
                     dplyr::bind_cols(tab[,c(1, 2)], .)

      tab$mean <- round(tab$mean, digits = 4)
      tab$range <- round(tab$range, digits = 4)

    }

    # as the sss values from copernicus are converted by me to pixel points and 
    # already refer to the center of the pixel (read faq at 
    # http://orca.science.oregonstate.edu/faq01.php), there is no need to shift the
    # coordinates of half a grid spacing to obtain the center
    if (par == "sss") {

      half_grid_spacing <- 0

    } else {

      half_grid_spacing <- (1/2 * 1/12)

    }

    # convert pixels to center of coordinates
    tab[, 1] <- round((( - 180 + (tab[, 1] / 12)) + half_grid_spacing), digits = 6)

    tab[, 2] <- round(((90 - (tab[, 2] / 12)) - half_grid_spacing), digits = 6)

    readr::write_delim(tab, 
                       file = paste0("variables/", par, "_", year_range),
                       delim = "\t",
                       na = "NA",
                       escape = "none",
                       append = FALSE)

  }

# As per argument 1, load data from 8.day periods of MODIS data, turn no data
# values "-9999" to NA and append environmental data to year file.

# Also, Net Primary Productivity will be calculated for each 8.day period

} else if (args[1] == 4) {

  name <- args[2]

  year <- stringr::str_remove_all(name, "\\-.*")

  base_params <- purrr::map(purrr::set_names(c("sst", "par", "chl", "k490")), function(var) {

    week_var <- paste(name, var, sep = "_")

    tab <- readr::read_delim(file = paste0("tmp/processed/week/", name, "_", var),
                             delim = "\t", 
                             na = "-9999", 
                             col_names = FALSE,
                             show_col_types = FALSE) %>%
           dplyr::rename(lon = X1, lat = X2, {{name}} := X3)

    # variable precision is reduce only for calculation of mean and range values
    # but kept at maximum precision for the calculation of NPP
    tab_npp_ancillary <- tab

    tab[, 3] <- round(tab[, 3], digits = 4)

    if (!file.exists(paste0("tmp/processed/year/", year, "_", var))) {

      readr::write_delim(tab,
                         file = paste0("tmp/processed/year/", year, "_", var),
                         delim = "\t",
                         na = "NA",
                         escape = "none")

    } else {

      final_tab <- readr::read_delim(file = paste0("tmp/processed/year/", year, "_", var), 
                                     delim = "\t",
                                     na = "NA",
                                     col_names = TRUE,
                                     show_col_types = FALSE)

      final_tab <- dplyr::bind_cols(final_tab, tab[, 3])

      readr::write_delim(final_tab, 
                         file = paste0("tmp/processed/year/", year, "_", var),
                         delim = "\t",
                         na = "NA",
                         escape = "none",
                         append = FALSE)

    }

    rm(tab)

    dplyr::rename(tab_npp_ancillary, {{var}} := {{name}})

  })

  tab <- dplyr::full_join(base_params[[1]], base_params[[2]]) %>% 
         dplyr::full_join(., base_params[[3]])

  gc()

  # get daylength mean for the date range specified in name (second argument 
  # supplied to Rscript ) and for each specific latitude point

  daylength_tab <- purrr::map_dbl(unique(tab$lat), function(point) {

    days <- seq(format(as.Date(stringr::str_remove(name, "_.*")), "%j"),
                format(as.Date(stringr::str_remove(name, ".*_")), "%j"))

    sapply(days, function(day) {

      lat <- (90 - (point * 1/12)) - (1/2 * 1/12)

      LatToDayLength(lat, day)

    }) %>% mean()

  })

  daylength_tab <- dplyr::bind_cols(tibble::tibble(lat = unique(tab$lat)), tibble::tibble(dayL = daylength_tab))

  # now comes the calculation of npp for each week in (milligrams Carbon per square meter per day)
  tab <- tab %>% dplyr::mutate(pb_opt = purrr::map_dbl(sst, pb_opt), 
                               z_eu = purrr::map_dbl(chl, euphotic),
                               irr = purrr::map_dbl(par, irradiance)) %>%
                 dplyr::left_join(daylength_tab, by = "lat") %>%
                 dplyr::transmute(npp = pb_opt * chl * dayL * irr * z_eu) %>%
                 dplyr::bind_cols(tab[, c(1, 2)], .) %>%
                 dplyr::rename({{name}} := npp)

  tab[, 3] <- round(tab[, 3], digits = 4)

  if (!file.exists(paste0("tmp/processed/year/", year, "_npp"))) {

    readr::write_delim(tab,
                       file = paste0("tmp/processed/year/", year, "_npp"),
                       delim = "\t",
                       na = "NA",
                       escape = "none")

  } else {

    final_tab <- readr::read_delim(file = paste0("tmp/processed/year/", year, "_npp"), 
                                   delim = "\t",
                                   na = "NA",
                                   col_names = TRUE,
                                   show_col_types = FALSE)

    final_tab <- dplyr::bind_cols(final_tab, tab[, 3])

    readr::write_delim(final_tab, 
                       file = paste0("tmp/processed/year/", year, "_npp"),
                       delim = "\t",
                       na = "NA",
                       escape = "none",
                       append = FALSE)

  }

# Calculate Normalized Seasonality Index for yearly NPP data, and then calculate
# range and mean values for other parameters, except for NPP for which also the 
# total yearly NPP will be summed.

} else if (args[1] == 5) {

  year <- args[2]

  parameters <- list.files("tmp/processed/year/") %>% stringr::str_remove(".*_")

  for (par in parameters) {

    tab <- readr::read_delim(file = paste0("tmp/processed/year/", year, "_", par), 
                             delim = "\t",
                             na = "NA",
                             col_names = TRUE,
                             show_col_types = FALSE)

    if (par == "npp") {

      # first, calculate the Normalized seasonality index ...
      tab$nsi <- apply(tab[, -c(1,2)], 1, function(x) {

                   nsi_calc(x)

                 })

      # ... then apply mean and sum
      tab <- tab %>% dplyr::select(-c("lon", "lat", "nsi")) %>% 
                     dplyr::mutate(m = rowMeans(., na.rm = TRUE), 
                                   tot = rowSums(., na.rm = TRUE)) %>%
                     dplyr::select(m, tot) %>%
                     dplyr::bind_cols(tab[, c("lon", "lat", "nsi")], .)

      tab[is.nan(tab$m),]$m <- NA
      tab[tab$tot == 0,]$tot <- NA

      tab$m <- round(tab$m, digits = 4)
      tab$nsi <- round(tab$nsi, digits = 4)
      tab$tot <- round(tab$tot, digits = 4)

      m <- paste0(year, "_m")
      tot <- paste0(year, "_tot")
      nsi <- paste0(year, "_nsi")

      tab <- tab %>% dplyr::rename({{m}} := m, 
                                   {{tot}} := tot, 
                                   {{nsi}} := nsi)

    } else if (any(par %in% c("chl", "par", "k490"))) {

      tab <- tab %>% dplyr::select(-c("lon", "lat")) %>% 
                     dplyr::mutate(m = rowMeans(., na.rm = TRUE)) %>%
                     dplyr::select(m) %>%
                     dplyr::bind_cols(tab[, c("lon", "lat")], .)

      tab[is.nan(tab$m),]$m <- NA

      tab$m <- round(tab$m, digits = 4)

      m <- paste0(year, "_m")

      tab <- tab %>% dplyr::rename({{m}} := m)

    } else if (par == "sst") {

      tab <- tab %>% dplyr::select(-c("lon", "lat")) %>% 
                     dplyr::mutate(m = rowMeans(., na.rm = TRUE), 
                                   max = do.call(pmax, c(., na.rm = TRUE)),
                                   min = do.call(pmin, c(., na.rm = TRUE))) %>%
                     dplyr::mutate(r = max - min) %>%
                     dplyr::select(m, r) %>%
                     dplyr::bind_cols(tab[, c("lon", "lat")], .)

      tab[is.nan(tab$m),]$m <- NA

      tab$m <- round(tab$m, digits = 4)
      tab$r <- round(tab$r, digits = 4)

      m <- paste0(year, "_m")
      r <- paste0(year, "_r")

      tab <- tab %>% dplyr::rename({{m}} := m, 
                                   {{r}} := r)

    }

    if (!file.exists(paste0("variables/", par, "_final"))) {

      readr::write_delim(tab, 
                         file = paste0("variables/", par, "_final"),
                         delim = "\t",
                         na = "NA",
                         escape = "none",
                         append = FALSE)

      rm(tab)

    } else {

      final_tab <- readr::read_delim(file = paste0("variables/", par, "_final"), 
                                     delim = "\t",
                                     col_names = TRUE,
                                     show_col_types = FALSE)

      final_tab <- dplyr::bind_cols(final_tab, tab[, -c(1, 2)])

      readr::write_delim(final_tab, 
                         file = paste0("variables/", par, "_final"),
                         delim = "\t",
                         na = "NA",
                         escape = "none",
                         append = FALSE)

      rm(final_tab)

    }

  }

}