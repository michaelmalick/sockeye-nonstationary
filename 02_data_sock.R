## Process sockeye data
##
## This script cleans and processes the raw downloaded sockeye data. The output of
## this script is a master brood table for the analysis and a summary info table.


## Read in downloaded data ---------------------------------
s.brood.raw <- read.table("./data/raw_brood_table_2019_03_31.csv",
                          sep = ",", skip = 0, stringsAsFactors = FALSE,
                          header = TRUE)

head(s.brood.raw)
tail(s.brood.raw)
nrow(s.brood.raw)
sapply(s.brood.raw, class)
summary(s.brood.raw)



## Clean-up column names -----------------------------------
r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood.raw), value = TRUE)

s.brood.nor <- s.brood.raw[ , !names(s.brood.raw) %in% r.cols]
names(s.brood.nor) <- tolower(names(s.brood.nor))
names(s.brood.nor) <- gsub(".", "_", names(s.brood.nor), fixed = TRUE)
names(s.brood.nor)[names(s.brood.nor) == "broodyear"] <- "brood_yr"
names(s.brood.nor)[names(s.brood.nor) == "totalescapement"] <- "spawners"
names(s.brood.nor)[names(s.brood.nor) == "totalrecruits"] <- "recruits"

## add ordering column
s.brood.nor$read_order <- 1:nrow(s.brood.nor)

s.brood <- cbind(s.brood.nor, s.brood.raw[ , names(s.brood.raw) %in% r.cols])


## Update ocean region
# o.region <- rep("GOA", nrow(s.brood))
# o.region <- ifelse(s.brood$region == "Fraser River", "WC", o.region)
# o.region <- ifelse(s.brood$region %in% c("Bristol Bay", "AYK"), "BS", o.region)
# s.brood$ocean_region <- o.region



## Clean-up master brood table -----------------------------

# NA values need to be replaced with 0s for years with no recruits
# This replacement is only done for the "recruit" columns
r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood), value = TRUE)
s.brood[ , r.cols][is.na(s.brood[ , r.cols])] <- 0


## Don't use Nushagak prior to 1985
# s.brood[s.brood$stock == "Nushagak", ]
s.brood$useflag[s.brood$stock == "Nushagak" & s.brood$brood_yr < 1985] <- 0


## 2011 Weaver total returns are suspect
# s.brood[s.brood$stock == "Weaver", ]
s.brood$useflag[s.brood$stock == "Weaver" & s.brood$brood_yr == 2011] <- 0


## Don't use L. Washington
## Questionable methods w/ brood table creation
s.brood$useflag[s.brood$stock == "Washington"] <- 0


## Subset usable data points
s.brood.use <- s.brood[s.brood$useflag == 1, ]


## Add R2.5 age-class if it doesn't exist
if(!"R2.5" %in% names(s.brood)) {
    s.brood.use$R2.5 <- 0
}

head(s.brood.use)
tail(s.brood.use)
nrow(s.brood.use)
sapply(s.brood.use, class)
summary(s.brood.use)



## Calculate age proportions -------------------------------

## Create data frame with estimates of recruits, spawners and proportion of
## recruits that entered ocean at age 0, 1 and 2 and from each age class
bt <- plyr::ddply(s.brood.use, c("stock_id", "brood_yr"),function(x) {
                      R <- sum(x[ , grep("^R[[:digit:]]\\.[[:digit:]]", names(x))],
                               na.rm = TRUE) # total recruits
                      S <- x$spawners # spawners
                      ocean_0 <- sum(x[ ,grep("^R0\\.", names(x))], na.rm = TRUE) / R
                      ocean_1 <- sum(x[ ,grep("^R1\\.", names(x))], na.rm = TRUE) / R
                      ocean_2 <- sum(x[ ,grep("^R2\\.", names(x))], na.rm = TRUE) / R
                      ocean_3 <- sum(x[ ,grep("^R3\\.", names(x))], na.rm = TRUE) / R
                      ocean_4 <- sum(x[ ,grep("^R4\\.", names(x))], na.rm = TRUE) / R
                      R0.1 <- x[ , "R0.1"] / R
                      R0.2 <- x[ , "R0.2"] / R
                      R0.3 <- x[ , "R0.3"] / R
                      R0.4 <- x[ , "R0.4"] / R
                      R0.5 <- x[ , "R0.5"] / R
                      R1.1 <- x[ , "R1.1"] / R
                      R1.2 <- x[ , "R1.2"] / R
                      R1.3 <- x[ , "R1.3"] / R
                      R1.4 <- x[ , "R1.4"] / R
                      R1.5 <- x[ , "R1.5"] / R
                      R2.1 <- x[ , "R2.1"] / R
                      R2.2 <- x[ , "R2.2"] / R
                      R2.3 <- x[ , "R2.3"] / R
                      R2.4 <- x[ , "R2.4"] / R
                      R2.5 <- x[ , "R2.5"] / R
                      R3.1 <- x[ , "R3.1"] / R
                      R3.2 <- x[ , "R3.2"] / R
                      R3.3 <- x[ , "R3.3"] / R
                      R3.4 <- x[ , "R3.4"] / R
                      R4.1 <- x[ , "R4.1"] / R
                      R4.2 <- x[ , "R4.2"] / R
                      R4.3 <- x[ , "R4.3"] / R
                      data.frame(stock = x$stock,
                                 jurisdiction = x$jurisdiction,
                                 region = x$region,
                                 ocean_region = x$ocean_region,
                                 sub_region = x$sub_region,
                                 lat = x$lat,
                                 lon = x$lon,
                                 recruits = R,
                                 spawners = S,
                                 ocean_0,
                                 ocean_1,
                                 ocean_2,
                                 ocean_3,
                                 ocean_4,
                                 R0.1, R0.2, R0.3, R0.4, R0.5,
                                 R1.1, R1.2, R1.3, R1.4, R1.5,
                                 R2.1, R2.2, R2.3, R2.4, R2.5,
                                 R3.1, R3.2, R3.3, R3.4,
                                 R4.1, R4.2, R4.3,
                                 stringsAsFactors = FALSE)
                      })

## Check that all Rx.x columns were included
r1 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood.use), value = TRUE)
r2 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(bt), value = TRUE)
if(!all.equal(sort(r1), sort(r2)))
    stop("Missing ages in brood table cleaning")

# Drop years with missing data
bt.out.1 <- bt[complete.cases(bt), ]


## Subsetting years
##
## In the data processing that Jeanette conducted, small age classes were filled
## in using the mean of the 2 previous if that age class was less than 10% of
## the total returns in the two previous years.
##
## From Jeanette on infilling small age classes:
##   "If age-specific abundance is NA and abundance is less than 10% of the
##   total brood return in the previous 2 brood years, then the abundance for
##   that age group is estiamted using the mean of the age-specific values from
##   the previous 2 brood years. If the age-specific value that is to be filled
##   in is greater than 10% of the total brood return, then that value remains
##   NA."
##
## From Jeanette on the UseFlag column:
##   "Set the UseFlag based on whether a row has a complete age class
##   estimation. A year of data is considered complete when all of the major
##   age classes have real values (not NA). A major age class is defined for
##   these purposes as any age class where the long term mean is greater than
##   1% of the total recruits for that population."
##
# bt.out.2 <- bt.out.1[bt.out.1$BY >= 1950 & bt.out.1$BY <= 2007, ]


## Fill in missing years that fall w/in min and max BY for each stock
bt.out.3 <- fill_time_series(bt.out.1)


## Limit time series start/end years within an ocean region
##   This limits the influence of one or a few stocks that
##   have data at before/after the majority of stocks in
##   an ocean region.
# bt.out.3wc  <- bt.out.3[bt.out.3$ocean_region == "WC"  &
#                         bt.out.3$brood_yr %in% 1950:2011, ]
# bt.out.3goa <- bt.out.3[bt.out.3$ocean_region == "GOA" &
#                         bt.out.3$brood_yr %in% 1965:2010, ]
# bt.out.3bs  <- bt.out.3[bt.out.3$ocean_region == "BS"  &
#                         bt.out.3$brood_yr %in% 1963:2012, ]
# bt.out.3a  <- rbind(bt.out.3wc, bt.out.3goa, bt.out.3bs)


## Trim time series of NA values by selecting the longest
## string of consecutive brood years for stocks with NA values.
bt.out.4 <- plyr::ddply(bt.out.3, .(stock_id), function(i) {
                        BY <- i$brood_yr
                        R  <- i$recruits
                        BY.na <- ifelse(is.na(R), NA, BY)
                        ind <- consec_years(BY.na)
                        return(i[i$brood_yr %in% ind, ])
                       })
# unique(bt.out.4$Stock)
# unique(bt.out.3$Stock[which(is.na(bt.out.3$R))])
# bt.out.4[which(is.na(bt.out.4$R))]


## Remove stocks w/ less than 20 brood years
nyr <- plyr::ddply(bt.out.4, .(stock_id), summarize,
                   n = sum(!is.na(recruits)))
bt.out.5 <- bt.out.4[!bt.out.4$stock_id %in%
                     nyr$stock_id[which(nyr$n < 20)], ]



## Create consecutive stock number -------------------------
bt.split <- split(bt.out.5, bt.out.5$stock_id)
bt.out.6 <- lapply(1:length(bt.split), function(i) {
                       dat.i <- bt.split[[i]]
                       dat.i$stock_no <- i
                       return(dat.i)
                    })
bt.out.6 <- plyr::rbind.fill(bt.out.6)



## Scale spawners and recruits -----------------------------
bt.out.7 <- bt.out.6
names(bt.out.7)[names(bt.out.7) == "spawners"] <- "spawners_raw"
names(bt.out.7)[names(bt.out.7) == "recruits"] <- "recruits_raw"
bt.out.7$spawners <- bt.out.7$spawners_raw / 1e5
bt.out.7$recruits <- bt.out.7$recruits_raw / 1e5



## Reorder columns -----------------------------------------
r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(bt.out.7), value = TRUE)
o.cols <- grep("^ocean\\_[[:digit:]]", names(bt.out.7), value = TRUE)
m.cols <- c("stock_id",
            "stock_no",
            "stock",
            "jurisdiction",
            "region",
            "ocean_region",
            "sub_region",
            "lat",
            "lon",
            "brood_yr",
            "recruits",
            "spawners",
            "recruits_raw",
            "spawners_raw")
bt.out.8 <- cbind(bt.out.7[ , m.cols],
                  bt.out.7[ , o.cols],
                  bt.out.7[ , r.cols])

head(bt.out.8)
tail(bt.out.8)
sapply(bt.out.8, class)



## Final brood table ---------------------------------------
brood.table <- bt.out.8
head(brood.table)
tail(brood.table)
sapply(brood.table, class)
summary(brood.table)



## Create stock info table ---------------------------------
sock.info <- plyr::ddply(brood.table, .(stock_id), summarize,
                         stock_no = unique(stock_no),
                         stock = unique(stock),
                         jurisdiction = unique(jurisdiction),
                         region = unique(region),
                         ocean_region = unique(ocean_region),
                         sub_region = unique(sub_region),
                         lat = unique(lat),
                         lon = unique(lon),
                         na_count = sum(is.na(recruits)),
                         n_years = sum(!is.na(recruits)),
                         yr_start = min(brood_yr),
                         yr_end = max(brood_yr))
sock.info$stock <- factor(sock.info$stock, levels = unique(sock.info$stock))



## Save outputs --------------------------------------------
save(brood.table, file = "./output/brood_table.RData")
save(sock.info, file = "./output/sock_info.RData")
