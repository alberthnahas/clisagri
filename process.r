### CLISAGRI ANALYSIS FOR KUPANG - FINAL SCRIPT WITH ALL PLOTS
# This script contains all necessary fixes and workflows, including the definitive hazard map fix.

### ===================================================================
### 1. SETUP AND LOAD DATA
### ===================================================================

# -- Load Required Packages --
library(zoo)
library(parallel)
library(ggplot2)
library(dplyr)
library(lubridate)
library(RColorBrewer)
library(scales)
library(metR)
library(tidyr)
library(gridExtra) # Added for combining plots

# -- Set Working Directory --
# setwd("/path/to/your/CLISAGRI/folder")

# -- Source All CLISAGRI R Functions --
tryCatch({
  sapply(list.files("source/", full.names = TRUE), source)
}, error = function(e) {
  stop("Could not find the 'source/' directory. Please make sure your working directory is set correctly.")
})


# -- Read and Clean Meteorological Data --
meteo <- read.csv("data/MeteoKupang.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
meteo$DAY <- as.Date(meteo$DAY, format = "%d/%m/%Y")
meteo <- meteo[!duplicated(meteo$DAY), ]
meteo$TEMPERATURE_MAX <- na.approx(meteo$TEMPERATURE_MAX, na.rm = FALSE)
meteo$TEMPERATURE_MIN <- na.approx(meteo$TEMPERATURE_MIN, na.rm = FALSE)
meteo$TEMPERATURE_AVG <- na.approx(meteo$TEMPERATURE_AVG, na.rm = FALSE)
meteo <- na.omit(meteo)

# -- Read Parameter and Sowing Files --
parameters <- read.csv("data/ParametersKupang.csv", header = TRUE, sep = ",")
sowing <- read.csv("data/SowingKupang.csv", header = TRUE, sep = ",")
sowing$DAY <- as.Date(sowing$DAY, format = "%d/%m/%Y")

# -- Adjust Parameters for Kupang's Tropical Climate --
parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE == "IDSL"] <- 0
if (!"TSUM1" %in% parameters$PARAMETER_CODE) {
  parameters <- rbind(parameters, data.frame(PARAMETER_CODE="TSUM1", PARAMETER_XVALUE=850, PARAMETER_DESCRIPTION="Temp sum emergence to anthesis"))
}
if (!"TSUM2" %in% parameters$PARAMETER_CODE) {
  parameters <- rbind(parameters, data.frame(PARAMETER_CODE="TSUM2", PARAMETER_XVALUE=1600, PARAMETER_DESCRIPTION="Temp sum anthesis to maturity"))
}


### ===================================================================
### 2. CORE MODEL AND FUNCTION PATCHES
### ===================================================================

# --- Patch 1: Override the flawed GDD calculation ---
dtsmtb <- function(temp) {
  t_base <- 10.0
  daily_sum <- ifelse(temp > t_base, temp - t_base, 0)
  daily_sum <- ifelse(daily_sum > 30, 30, daily_sum)
  return(daily_sum)
}

# --- Patch 2: Create a corrected phenology wrapper ---
phenology_fixed <- function(meteo, f.variety, sowing_dates, lat) {
  meteo$DVS <- NA
  for (sow_date in sowing_dates) {
    start_index <- which(meteo$DAY == sow_date)
    if (length(start_index) > 0) {
      end_index <- min((start_index + 365), nrow(meteo))
      meteo_season <- meteo[start_index:end_index,]
      
      self <- initialize(list(), sow_date, f.variety, start.type="sowing", end.type="maturity", vernalisation = FALSE)
      
      days <- 1
      drv <- list(LAT = lat)
      
      while (self$states$DVS <= 2 && days <= nrow(meteo_season) && !is.na(meteo_season$TEMPERATURE_AVG[days])) {
        current_day <- meteo_season$DAY[days]
        drv$TEMP <- meteo_season$TEMPERATURE_AVG[days]
        
        if (self$states$STAGE == "emerging") {
          dtsume <- drv$TEMP - self$params$TBASEM
          dtsume[dtsume < 0] <- 0
          self$rates$DTSUME <- dtsume
          self$rates$DVR <- 0
        } else {
          self$rates$DTSUME <- 0
          dtsum <- dtsmtb(drv$TEMP)
          if (self$states$DVS < 1) {
            self$rates$DVR <- dtsum / self$params$TSUM1
          } else {
            self$rates$DVR <- dtsum / self$params$TSUM2
          }
        }
        
        self <- integrate.dvs(self, current_day)
        meteo$DVS[meteo$DAY == current_day] <- self$states$DVS
        days <- days + 1
      }
    }
  }
  return(meteo)
}

# --- Patch 3: Override the hazard.map and dvs.plot functions ---
hazard.map = function(data, type, year=NULL)
{
  # Define colors and breaks based on type
  if(type >= 2 && type <= 6){
    colors = brewer.pal(n=7,name="RdBu")
    breaks = c(-2,-1.5,-1,1,1.5,2)
    labels=c(-2,"",-1,1,"",2)
    limits=c(-max(abs(data$value), 2, na.rm=TRUE), max(abs(data$value), 2, na.rm=TRUE))
  } else {
    colors = brewer.pal(n=9, name="YlOrRd")
    breaks = pretty_breaks()
    labels = waiver()
    limits = c(min(data$value, na.rm=TRUE), max(data$value, na.rm=TRUE))
  }
  
  if(is.null(year)) year = unique(data$year)[1]
  data_year = data[data$year == year,]
  
  if(nrow(data_year) < 10) { # Need enough points for loess/contours
    warning(paste("Not enough valid data points to generate a complete map for year", year))
    return(NULL)
  }
  
  # Use loess to get smoothed data for contour lines
  data_year$maturity.s = loess(maturity ~ tsum1 * tsum2, data = data_year, span=0.75, na.action=na.exclude)$fitted
  data_year$flow.s = loess(flower ~ tsum1 * tsum2, data = data_year, span=0.75, na.action=na.exclude)$fitted
  
  g1 = ggplot(data_year, aes(x=tsum1, y=tsum2)) + 
    # DEFINITIVE FIX: Use geom_contour_fill to interpolate from sparse data
    geom_contour_fill(aes(z=value), na.fill=TRUE) +
    # Overlay the contour lines for phenology
    geom_contour(aes(z=maturity.s), colour="black", na.rm=TRUE) +
    geom_contour(aes(z=flow.s), colour="grey40", na.rm=TRUE) +
    # Apply the color scale to the fill
    scale_fill_gradientn(breaks=breaks, colours=colors, labels=labels, limits=limits, na.value="transparent") +
    # Add labels to the contour lines
    geom_text_contour(aes(z=maturity.s), stroke=0.2, size=4, check_overlap=TRUE, na.rm=TRUE) +
    geom_text_contour(aes(z=flow.s), stroke=0.2, size=4, color="grey40", check_overlap=TRUE, na.rm=TRUE) +
    coord_fixed() +
    theme_bw() +
    labs(x="TSUM1 [GDD]", y="TSUM2 [GDD]", fill="Value", title=paste("Hazard Map, Type:", type, "Year:", year))
  
  print(g1)
  return(g1)
}

dvs.plot = function(data)
{
  tsum1 = unique(data$tsum1)
  tsum2 = unique(data$tsum2)
  
  z <- outer(seq(-1,1,length=length(tsum1)), seq(-1,1,length=length(tsum2)), FUN = fun_xy)
  
  Z = cbind(expand.grid(tsum1 = tsum1, tsum2 = tsum2), color = as.vector(z))
  gz = ggplot(Z, aes(x=tsum1, y=tsum2)) +
    geom_tile(aes(fill=color)) +
    scale_fill_identity() +
    labs(x="TSUM1", y="TSUM2") +
    theme_minimal(base_size=8) +
    theme(axis.text=element_text(size=6))
  
  data$color = NA
  for(i in 1:length(tsum1)) {
    for(j in 1:length(tsum2)) {
      data$color[which(data$tsum1==tsum1[i] & data$tsum2==tsum2[j])] = z[i,j]
    }
  }
  
  dvs.stages = data.frame(DVS=c(0.01, 0.1, 0.34, 0.72, 0.82, 0.92, 1.16, 2.0), 
                          NAME=c("emergence", "tillering", "stem elongation", "booting", "heading", "flowering", "grain filling", "maturity"))
  
  g1 = ggplot(data, aes(x=DAY, y=DVS, group=color)) +
    geom_path(aes(color=color), size=0.5, alpha=0.6) +
    scale_color_identity() + 
    theme_bw() +
    labs(x="Day", y="Development Stage (DVS)", title="DVS Trajectories for a Single Season") +
    scale_x_date(breaks = pretty_breaks(n=10), date_minor_breaks="1 month") + 
    scale_y_continuous(breaks=dvs.stages$DVS, labels=dvs.stages$NAME) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 2), xlim=c(min(data$DAY), max(data$DAY)), expand=FALSE)
  
  final_plot <- grid.arrange(g1, gz, ncol=2, widths = c(3, 1))
  
  return(final_plot)
}


### ===================================================================
### 3. WORKFLOW DEFINITIONS (HAZARD + DVS PLOT)
### ===================================================================

par.fun.hazard <- function(X, meteo, latitude, parameters, combin, sowing, type) {
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE == "TSUM1"] <- combin[X, 1]
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE == "TSUM2"] <- combin[X, 2]
  meteo_pheno <- phenology_fixed(meteo, parameters, sowing$DAY, latitude)
  meteo_pheno$DVS <- DVS2BBCH(meteo_pheno$DVS, dvs.end = 2)
  dvs_results <- clisagri(meteo = meteo_pheno, types = type, lat = latitude, sowing = sowing$DAY)
  if (is.null(dvs_results) || nrow(dvs_results) == 0) return(NULL)
  stage_dates <- sowing %>%
    mutate(year = ifelse(yday(DAY) > 180, year(DAY) + 1, year(DAY))) %>%
    rowwise() %>%
    mutate(
      season_data = list(meteo_pheno %>% filter(DAY >= DAY, DAY <= DAY + years(1))),
      flower = if(any(season_data$DVS >= 65, na.rm=T)) yday(season_data$DAY[which.min(abs(season_data$DVS - 65))]) else NA,
      maturity = if(any(season_data$DVS >= 89, na.rm=T)) yday(season_data$DAY[which.min(abs(season_data$DVS - 89))]) else NA
    ) %>%
    ungroup() %>%
    select(year, flower, maturity)
  dvs_results <- left_join(dvs_results, stage_dates, by = "year")
  return(dvs_results)
}

hazard.calc <- function(meteo, sowing, parameters, latitude, r.tsum1, r.tsum2, type, parallel = FALSE, ncores = 2) {
  combin <- expand.grid(tsum1 = r.tsum1, tsum2 = r.tsum2)
  par.fun.safe <- function(i) {
    try(par.fun.hazard(i, meteo, latitude, parameters, combin, sowing, type), silent = TRUE)
  }
  results <- if (parallel) {
    mclapply(1:nrow(combin), par.fun.safe, mc.cores = ncores)
  } else {
    lapply(1:nrow(combin), par.fun.safe)
  }
  bind_rows(lapply(1:length(results), function(i) {
    res <- results[[i]]
    if (is.data.frame(res) && nrow(res) > 0 && "year" %in% names(res)) {
      res$tsum1 <- combin$tsum1[i]
      res$tsum2 <- combin$tsum2[i]
      return(res)
    }
    return(NULL)
  }))
}

par.fun.dvs <- function(X, meteo, latitude, parameters, combin, sowing_df) {
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE == "TSUM1"] <- combin[X, 1]
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE == "TSUM2"] <- combin[X, 2]
  meteo_pheno <- phenology_fixed(meteo, parameters, sowing_df$DAY, latitude)
  result <- meteo_pheno %>% 
    filter(!is.na(DVS)) %>%
    select(DAY, DVS)
  if (nrow(result) > 0) { return(result) }
  return(NULL)
}

dvs.calc <- function(meteo, sowing_df, parameters, latitude, r.tsum1, r.tsum2, parallel = FALSE, ncores = 2) {
  combin <- expand.grid(tsum1 = r.tsum1, tsum2 = r.tsum2)
  par.fun.safe <- function(i) {
    try(par.fun.dvs(i, meteo, latitude, parameters, combin, sowing_df), silent = TRUE)
  }
  results <- if (parallel) {
    mclapply(1:nrow(combin), par.fun.safe, mc.cores = ncores)
  } else {
    lapply(1:nrow(combin), par.fun.safe)
  }
  bind_rows(lapply(1:length(results), function(i) {
    res <- results[[i]]
    if (is.data.frame(res) && nrow(res) > 0) {
      res$tsum1 <- combin$tsum1[i]
      res$tsum2 <- combin$tsum2[i]
      return(res)
    }
    return(NULL)
  }))
}


### ===================================================================
### 4. EXECUTE ANALYSIS AND PLOT RESULTS
### ===================================================================

# -- Define the parameter ranges to test --
tsum1_range <- seq(700, 1000, by = 50)
tsum2_range <- seq(1400, 1800, by = 50)

# -- Run the hazard calculation --
message("Running hazard calculation for type 6...")
type.6.kupang <- hazard.calc(
  meteo = meteo,
  sowing = sowing,
  parameters = parameters,
  latitude = -10.2,
  r.tsum1 = tsum1_range,
  r.tsum2 = tsum2_range,
  type = 6,
  parallel = TRUE
)

# -- Plot hazard results if successful --
if (!is.null(type.6.kupang) && nrow(type.6.kupang) > 0) {
  
  # DEFINITIVE FIX: Find the year with the MOST valid data points to ensure a good plot
  best_year_df <- type.6.kupang %>%
    filter(!is.na(value) & !is.na(flower) & !is.na(maturity)) %>%
    count(year, name = "valid_points") %>%
    arrange(desc(valid_points))
  
  if(nrow(best_year_df) > 0){
    target_year <- best_year_df$year[1]
    message(paste("Best year for plotting hazard map found:", target_year, "with", best_year_df$valid_points[1], "valid points."))
  } else {
    target_year <- NULL
    warning("No single year produced enough valid data for a complete map.")
  }
  
  if(!is.null(target_year)){
    try(hazard.map(type.6.kupang, type = 6, year = target_year), silent = TRUE)
  }
  
  message("Generating hazard time series plot...")
  try(clim.plot(type.6.kupang, type = 6), silent = TRUE)
  
} else {
  warning("Hazard calculation failed to produce any results.")
}


# --- Run DVS data generation for a single representative season ---
message("Generating data for DVS plot...")

target_sowing_date <- sowing$DAY[1] 
sowing_for_plot <- data.frame(DAY = target_sowing_date)

dvs_plot_data <- dvs.calc(
  meteo = meteo,
  sowing_df = sowing_for_plot,
  parameters = parameters,
  latitude = -10.2,
  r.tsum1 = tsum1_range,
  r.tsum2 = tsum2_range,
  parallel = TRUE
)

if (!is.null(dvs_plot_data) && nrow(dvs_plot_data) > 0) {
  message("Generating combined DVS plot and legend...")
  try(dvs.plot(dvs_plot_data), silent = TRUE)
} else {
  warning("DVS data generation failed to produce any results.")
}