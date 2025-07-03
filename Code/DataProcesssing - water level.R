# Data is water level per wetland in the six regions. 
# There are three documents:
  # 1) wetlands and demes.csv
  # 2) snailkite counts by wetland 1997_2025.csv
  # 3) daily stage in each wetland.csv

# And we are aiming to generate three matrices to use in the second objectiv

  # 1) the counts (three surveys as replicates)
  # 2) the mean monthly water level (should be corrected with dates of surveys?)
  # 3) the %CV of the monthly water level (also to be updated)

library(tidyverse)

# Reading the data

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

# long format
snki_N_wetland_region <- snki_N_wetland |>
  # move wetland columns to a single, having the count reported
  pivot_longer(cols = !c(X, year, survey_num), 
               names_to = "wetland", 
               values_to = "count") |>
  # add a potential date per survey (at the level of month)
  mutate(month = ifelse(survey_num == 1,"March",
                 ifelse(survey_num == 2,"April",
                 ifelse(survey_num == 3,"May",
                 ifelse(survey_num == 4,"June",
                 ifelse(survey_num == 5,"July",
                 ifelse(survey_num == 6,"August",
                               "none")))))),
         date = ymd(paste0(year,"-",month,"01"))) |>
  # assign region from the wetland name.
  left_join(wetland_region)

# water level data per wetland - daily to monthly to match survey ####
daily_wetland_water <- read.csv("Code/daily stage in each wetland.csv") |>
  mutate(date = dmy(Daily.Date),
         date = floor_date(date,"month"))
summary(daily_wetland_water)

# long format - summarized by month
monthly_wetland_water <- daily_wetland_water |>
  pivot_longer(cols = !c("date", "Daily.Date"),
               names_to = "wetland", 
               values_to = "daily_water") |>
  group_by(date, wetland) |>
  summarise(mu_month_water = mean(daily_water, na.rm = T),
            sd_month_water = sd(daily_water, na.rm = T),
            CV_month_water = (sd_month_water/mu_month_water)*100)

# Unify counts with Month mean water and variation ####
# 30 wetlands have water level information
# 46 wetlands have counts in 6 surveys

# which wetland in each region has more information?
snki_N_wetland_region |>
  # joint the monthly level of water per wetland
  left_join(monthly_wetland_water) |> 
  # to summarize, `group_by()`
  group_by(region,wetland) |>
  # remove NA surveys
  drop_na(count) |>
  # how many surveys per wetland?
  count() |> 
  # sort the information (for Table 1 in the report)
  arrange(region, -n) |> 
  as.data.frame()

snki_N_water <- snki_N_wetland_region |>
  left_join(monthly_wetland_water) |>
  ungroup() |>
  # here we can select the wetlands with more data, or representative of the demes
  filter(wetland %in% c(
    "East.Lake.Tohopekaliga", # KRV
    "SJM", # SJM
    "WCA.3A", # EVER
    "Lake.Okeechobee", # OKEE
    "Grassy.Waters.Preserve", # EAST
    "Paynes.Prairie" # PP
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
  dplyr::select(year, survey_num, region.wetland, mu_month_water) |>
  pivot_wider(names_from = region.wetland, values_from = mu_month_water) |>
  as.data.frame()

snki_water_CV <- snki_N_water |> 
  dplyr::select(year, survey_num, region.wetland, CV_month_water) |>
  pivot_wider(names_from = region.wetland, values_from = CV_month_water) |>
  as.data.frame()

names.snkicounts <- names(snki_counts)
wet.names <- c("EAST_Grassy.Waters.Preserve",
               "EVER_WCA.3A",
               "KRV_East.Lake.Tohopekaliga",
               "OKEE_Lake.Okeechobee",
               "PP_Paynes.Prairie",
               "SJM_SJM")
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

# table of region, wetland, latitude, longitude, and relationship of water level