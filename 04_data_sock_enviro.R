## Combine sockeye + environmental data
##
## This script assigns annual environmental data to each stock specific brood
## year. Because sockeye populations can enter the ocean over multiple years,
## weighted means of the annual environmental data are calculated to index ocean
## conditions for a brood year, where the weights are equal to the proportion of
## sockeye enter the ocean in a given year.


# Prep data ------------------------------------------------
sock.covar  <- brood.table



## Calc stock specific SST averages ------------------------
sst.stock <- sst_averager(sock.info, sst.anom, distance = 400)
names(sst.stock) <- c("year","sst_raw", "sst_anom", "stock_id")



## Add enviro indices --------------------------------------

## SST
sock.covar$sst_anom <- wgt_mean_sst(sock.covar,
                                    sst.stock,
                                    sst.var = "sst_anom")
sock.covar$sst_raw  <- wgt_mean_sst(sock.covar,
                                    sst.stock,
                                    sst.var = "sst_raw")


## PDO
sock.covar$pdo_winter <- wgt_mean_pdo(sock.covar, pdo.winter)
sock.covar$pdo_spring <- wgt_mean_pdo(sock.covar, pdo.spring)
sock.covar$pdo_summer <- wgt_mean_pdo(sock.covar, pdo.summer)


## NPGO
sock.covar$npgo_winter <- wgt_mean_pdo(sock.covar, npgo.winter)
sock.covar$npgo_spring <- wgt_mean_pdo(sock.covar, npgo.spring)
sock.covar$npgo_summer <- wgt_mean_pdo(sock.covar, npgo.summer)



## Add derived columns -------------------------------------
sock.covar <- plyr::ddply(sock.covar, .(stock_no), transform,
                          rps = recruits / spawners,
                          rps_stnd = scale(recruits / spawners)[ , 1],
                          ln_rps = log(recruits / spawners),
                          rk_resid = residuals(lm(log(recruits / spawners) ~ spawners)),
                          sst_anom_stnd = scale(sst_anom)[ , 1],
                          sst_raw_stnd = scale(sst_raw)[ , 1],
                          pdo_winter_stnd = scale(pdo_winter)[ , 1],
                          pdo_spring_stnd = scale(pdo_spring)[ , 1],
                          pdo_summer_stnd = scale(pdo_summer)[ , 1],
                          npgo_winter_stnd = scale(npgo_winter)[ , 1],
                          npgo_spring_stnd = scale(npgo_spring)[ , 1],
                          npgo_summer_stnd = scale(npgo_summer)[ , 1])
sock.covar$stock <- factor(sock.covar$stock,
                           levels = unique(sock.covar$stock))



## Save output ---------------------------------------------
save(sst.stock, file = "./output/sst_stock.RData")
save(sock.covar, file = "./output/sock_covar.RData")
