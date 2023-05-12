## Process environmental data
##
## This script processes the environmental data, including:
##  - SST   monthly anomalies
##  - PDO   multi-month averages
##  - NPGO  multi-month averages


## Load raw SST data ---------------------------------------
sst.raw <- read.csv("./data/sst_raw.csv")

head(sst.raw)
tail(sst.raw)
sapply(sst.raw, class)
summary(sst.raw)

## Plot centers of SST grid cells
# dat1 <- sst.raw[sst.raw$year == 2015 & sst.raw$month == 3, ]
# plot(dat1$lon, dat1$lat, pch = 16)
# map('world2', add = TRUE, col = "grey50")



## Calc SST anomalies --------------------------------------
sst.anom <- sst_anomaly(sst.raw, ref.years = 1950:2010)

head(sst.anom)
tail(sst.anom)
summary(sst.anom)
sapply(sst.anom, class)

## Convert longitude to match salmon data
sst.anom$lon <- ifelse(sst.anom$lon > 180, sst.anom$lon - 360, sst.anom$lon)

## Exclude data west of date line
sst.anom <- sst.anom[sst.anom$lon < 0, ]

## Plot centers of SST grid cells with new lon
# dat1 <- sst.anom[sst.anom$year == 2015 & sst.anom$month == 3, ]
# plot(dat1$lon, dat1$lat, pch = 16)
# plot(countriesLow, add = TRUE, border = "grey50")



## PDO -----------------------------------------------------
pdo <- read.csv("./data/pdo.csv")

## Calculate PDO index
pdo.winter <- enviro_avg_months(pdo, 12, 3, "pdo")
pdo.spring <- enviro_avg_months(pdo, 4, 5, "pdo")
pdo.summer <- enviro_avg_months(pdo, 6, 8, "pdo")



## NPGO ----------------------------------------------------
npgo <- read.csv("./data/npgo.csv")

## Calculate NPGO index
npgo.winter <- enviro_avg_months(npgo, 12, 3, "npgo")
npgo.spring <- enviro_avg_months(npgo, 4, 5, "npgo")
npgo.summer <- enviro_avg_months(npgo, 6, 8, "npgo")



## Save outputs --------------------------------------------
save(sst.anom, file = "./output/sst_anom.RData")
save(pdo.winter, file = "./output/pdo_winter.RData")
save(pdo.spring, file = "./output/pdo_spring.RData")
save(pdo.summer, file = "./output/pdo_summer.RData")
save(npgo.winter, file = "./output/npgo_winter.RData")
save(npgo.spring, file = "./output/npgo_spring.RData")
save(npgo.summer, file = "./output/npgo_summer.RData")
