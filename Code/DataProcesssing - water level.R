# Data is water level per wetland in the six regions. 
# There are four documents:
  # 1) survey_dates_1997_2025.csv
  # 2) daily stage in each wetland.csv
  # 3) wetlands and demes.csv
  # 4) snailkite counts by wetland 1997_2025.csv

# And we are aiming to generate three matrices to use in the second objective

  # 1) the counts (three surveys as replicates)
  # 2) the mean weekly water level (including the week before starting each survey)
  # 3) the %CV of the weekly water level (SD/mean * 100). 
      # This allow to compare across surveys with different length

library(tidyverse)

# Reading the data

# surveys dates to link water level data ####
survey_date <- read.csv("Code/survey_dates_1997_2025.csv") |>
  mutate(start = ymd(start),
         end = ymd(end),
         start.w = floor_date(start-6,"week"),
         end.w = floor_date(end, "week"))

# water level data per wetland - daily value in long format ####
daily_wetland_water <- read.csv("Code/daily stage in each wetland.csv") |>
  mutate(date = dmy(Daily.Date)) |>
  pivot_longer(cols = !c("date", "Daily.Date"),
               names_to = "wetland",
               values_to = "daily_water")

# Filter water level at the weekly scale for the dates of surveys
www_filt <- pmap(survey_date, function(start.w, end.w, ...) {
  daily_wetland_water |>
    filter(date >= start.w & date <= end.w)
})

www_filt_df <- bind_rows(www_filt, .id = "survey_id")
www_filt_df$X <- as.numeric(www_filt_df$survey_id)

# relationship of wetlands and regions ####
wetland_region <- read.csv("Code/wetlands and demes.csv") |>
  # adjust name "regions"
  rename(region = deme) |>
  # unify names (spaces and dash as points)
  mutate(wetland = gsub("[- ]", ".", wetland))

# count of number of wetlands per region
table(wetland_region$region)

# Counts per wetland in 6 regions ####
snki_N_wetland <- read.csv("Code/snailkite counts by wetland 1997_2025.csv")
names(snki_N_wetland)

# long format with water level per survey
snki_N_wetland_region <- snki_N_wetland |>
  # move wetland columns to a single, having the count reported
  pivot_longer(cols = !c(X, year, survey_num), 
               names_to = "wetland", 
               values_to = "count") |> 
  # merge with the daily water wetland level
  left_join(www_filt_df) |>
  # extract survey mean and CV
  group_by(wetland, year, survey_id, survey_num) |>
  summarise(
    count = mean(count, na.rm = TRUE),
    mu_water = mean(daily_water, na.rm = TRUE),
    sd_water = sd(daily_water, na.rm = TRUE),
    CV_water = (sd_water/mu_water) * 100
            ) |>
  # assign region from the wetland name
  left_join(wetland_region)

# Unify counts with survey mean water and variation ####
# 30 wetlands have water level information
# 46 wetlands have counts in 6 surveys

# which wetland in each region has more information?
snki_N_wetland_region |>
  mutate(region = ifelse(wetland == "STA.Lakeside.Ranch","OTHER",
                         region)) |>
  # to summarize, `group_by()`
  group_by(region, wetland) |>
  # remove NA surveys
  drop_na(count) |>
  # how many surveys per wetland?
  count() |> 
  # sort the information (for Table 1 in the report)
  arrange(region, -n) |> 
  filter(region != "OTHER") |>
  as.data.frame()

# select the wetlands with ≥ 75% of the maximum counts per region

snki_N_water <- snki_N_wetland_region |>
  ungroup() |>
  # here we can select the wetlands with more data to represent each region
  filter(wetland %in% c(
    "Grassy.Waters.Preserve", # EAST
    "WCA.3A", # EVER
    "WCA.2B", # EVER
    "Loxahatchee.NWR", # EVER
    "WCA.3B", # EVER
    "East.Lake.Tohopekaliga", # KRV
    "Lake.Kissimmee", # KRV
    "Lake.Istokpoga", # KRV
    "Lake.Okeechobee", # OKEE
    "Paynes.Prairie", # PP
    "SJM" # SJM
    )) |>
  # unify in a single name
  mutate(region.wetland = paste(region, wetland, sep = "_")) |>
  as.data.frame()

# with this, we can generate three matrices per wetland-region: 
  # 1) the counts
  # 2) the mean monthly water level
  # 3) the %CV of the monthly water level

snki_counts <- snki_N_water |> 
  dplyr::select(year, survey_num, region.wetland, count) |>
  pivot_wider(names_from = region.wetland, values_from = count) |>
  as.data.frame()

snki_water <- snki_N_water |> 
  dplyr::select(year, survey_num, region.wetland, mu_water) |>
  pivot_wider(names_from = region.wetland, values_from = mu_water) |>
  as.data.frame()

snki_water_CV <- snki_N_water |> 
  dplyr::select(year, survey_num, region.wetland, CV_water) |>
  pivot_wider(names_from = region.wetland, values_from = CV_water) |>
  as.data.frame()

names.snkicounts <- names(snki_counts)
wet.names <- c(
  "EAST_Grassy.Waters.Preserve", # EAST
  "EVER_WCA.3A", # EVER
  "EVER_WCA.2B", # EVER
  "EVER_Loxahatchee.NWR", # EVER
  "EVER_WCA.3B", # EVER
  "KRV_East.Lake.Tohopekaliga", # KRV
  "KRV_Lake.Kissimmee", # KRV
  "KRV_Lake.Istokpoga", # KRV
  "OKEE_Lake.Okeechobee", # OKEE
  "PP_Paynes.Prairie", # PP
  "SJM_SJM" # SJM
  )
cols.with.counts <- match(wet.names, names.snkicounts)
years <- unique(snki_counts$year)
nyears <- length(years)
year.col <- which(names.snkicounts=="year",arr.in=TRUE)
surv.num <- unique(snki_counts$survey_num)
tt <- years - years[1]

S4mat <- S5mat <- S6mat <- matrix(0, nrow=6, ncol=nyears)

Counts.list <- list(S4mat = S4mat, S5mat=S5mat, S6mat=S6mat)
Water.list <- list(W4mat =  S4mat, W5mat=S5mat, W6mat=S6mat)
WaterCV.list <- list(Wcv4mat =  S4mat, Wcv5mat=S5mat, Wcv6mat=S6mat)

for(i in 1:3){
  SN <- surv.num[(i+3)]
  ith.rows   <- which(snki_counts$survey_num==SN, arr.ind=TRUE)
  
  ith.datmat <- t(snki_counts[ith.rows, cols.with.counts])
  ith.watmat <- t(snki_water[ith.rows, cols.with.counts])
  ith.watCVmat <- t(snki_water_CV[ith.rows, cols.with.counts])
  colnames(ith.datmat) <- colnames(ith.watmat) <- colnames(ith.watCVmat) <- tt
  
  Counts.list[[i]] <- ith.datmat
  Water.list[[i]] <- ith.watmat
  WaterCV.list[[i]] <- ith.watCVmat
}
