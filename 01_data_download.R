## Download data for project
##
## This script downloads the data needed for the project and writes it to CSV
## files for processing in other scripts. All downloaded data are saved in the
## "./data" directory.


## Master brood table --------------------------------------
## doi:10.5063/F18S4N6X
bt.url <- "https://cn.dataone.org/cn/v2/resolve/urn:uuid:c747eee7-89fd-44ba-bce1-a14d0670792d"
download.file(bt.url, "./data/raw_brood_table_2019_03_31.csv")

bi.url <- "https://cn.dataone.org/cn/v2/resolve/urn:uuid:178c61cc-6119-42c2-9387-5b0b602324d4"
download.file(bi.url, "./data/raw_stock_info_2019_03_31.csv")



## SST raw -------------------------------------------------
ersst::sst_download(years = 1950:2016,
                    months = 1:12,
                    save.dir = "./data/sst_raw/",
                    version = 5)

sst.raw.full <- ersst::sst_load(years = 1950:2016,
                                months = 1:12,
                                read.dir = "./data/sst_raw/",
                                version = 5)

sst.raw.np <- ersst::sst_subset_space(sst.raw.full,
                                      lat.min = 36,
                                      lat.max = 80,
                                      lon.min = 170,
                                      lon.max = 250)

sst.raw.df <- ersst::sst_dataframe(sst.raw.np)

write.csv(sst.raw.df, "./data/sst_raw.csv", row.names = FALSE)



## PDO + NPGO ----------------------------------------------
years <- 1950:2016
pdo   <- get_pdo(years)
npgo  <- get_npgo(years)

write.csv(pdo, "./data/pdo.csv", row.names = FALSE)
write.csv(npgo, "./data/npgo.csv", row.names = FALSE)
