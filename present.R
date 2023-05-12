## Presentation figures


if(!dir.exists("./figures/present/"))
    dir.create("./figures/present/")


col.stock <- RColorBrewer::brewer.pal(8, "Dark2")
col.sst <- chroma::dpal(500, hue = c(240, 0), chroma = 70, power = 1.0)

## Map data
cl <- rnaturalearth::ne_states(country = c("United States of America",
                                           "Mexico",
                                           "Canada"))
namerica.state.sp <- sp::spTransform(cl, sp::CRS("+init=epsg:4326"))
namerica.state.df  <- suppressMessages(ggplot2::fortify(namerica.state.sp))
dat.map <- ocean_region_lab(sock.info)


## Fig: R/S time series ------------------------------------

dat <- sock.covar[sock.covar$stock_id == 102 & sock.covar$brood_yr > 1951, ]
g <- ggplot(dat) +
    xlab("Brood year") +
    ylab("Productivity (R/S)") +
    geom_line(aes(x = brood_yr, y = rps), color = "steelblue") +
    geom_point(aes(x = brood_yr, y = rps), color = "steelblue") +
    theme_sleek(base_size = 14)
print(g)
ggsave("./figures/present/early_stuart_rps.pdf", width = 5, height = 3)
ggsave("./figures/present/early_stuart_rps.png", width = 5, height = 3)
ggsave("./figures/present/early_stuart_rps.jpg", width = 5, height = 3)

dat <- sock.covar[sock.covar$stock_id == 156 & sock.covar$brood_yr > 1964, ]
g <- ggplot(dat) +
    xlab("Brood year") +
    ylab("Productivity (R/S)") +
    geom_line(aes(x = brood_yr, y = rps), color = "steelblue") +
    geom_point(aes(x = brood_yr, y = rps), color = "steelblue") +
    theme_sleek(base_size = 14)
print(g)
ggsave("./figures/present/egegik_rps.pdf", width = 5, height = 3)
ggsave("./figures/present/egegik_rps.png", width = 5, height = 3)
ggsave("./figures/present/egegik_rps.jpg", width = 5, height = 3)



dat <- sock.covar[sock.covar$stock_id == 102, ]
    g <- ggplot(dat) +
        xlab("SST anomaly") +
        ylab("Productivity (R/S)") +
        geom_point(aes(x = sst_anom_stnd, y = rps), color = "tomato") +
        theme_sleek(base_size = 14)
    print(g)
ggsave("./figures/present/early_stuart_sst.pdf", width = 5, height = 4)
ggsave("./figures/present/early_stuart_sst.png", width = 5, height = 4)


## Fig: random walk ----------------------------------------
set.seed(129)
rw <- cumsum(rnorm(100))

pdf("./figures/present/random_walk.pdf", width = 5, height = 3)
    plot(rw, type = "l",
         col = "grey40",
         lwd = 6,
         axes = FALSE,
         xlab = "",
         ylab = "")
dev.off()



## Fig: Map ------------------------------------------------

range(sock.info$n_years)
range(sock.covar$brood_yr)
r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(sock.covar), value = TRUE)
sum(sapply(sock.covar[ , r.cols], sum) > 0)


## TODO handle multiple stocks w/ same ocean entry location
## [X] Fraser River
## [ ] Upper Station
## [ ] Karluk
## [ ] Chignik

n.fr <- sum(sock.info$region == "Fraser River")
fr <- data.frame(lon = dat.map$lon[1],
                 lat = seq(dat.map$lat[1], dat.map$lat[1] + 4, length.out = n.fr))
n.us <- sum(sock.info$region == "Fraser River")



g <- ggplot(dat.map) +
    geom_map(data = namerica.state.df,
             map = namerica.state.df,
             aes(map_id = id),
             fill = "grey90", color = "grey60", size = 0.1) +
    geom_point(data = fr, aes(x = lon, y = lat),
               color = "grey40", size = 0.5) +
    geom_point(aes(x = lon, y = lat,
                   fill = ocean_region_lab,
                   shape = ocean_region_lab),
               size = 3) +
    geom_point(aes(x = lon, y = lat,
                   shape = ocean_region_lab,
                   color = ocean_region_lab),
               size = 3, alpha = 0.8) +
    xlab(expression(paste("Longitude (", degree, "E)"))) +
    ylab(expression(paste("Latitude (", degree, "N)"))) +
    scale_color_manual(values = col.stock) +
    scale_fill_manual(values = scales::alpha(col.stock, 0.5)) +
    scale_shape_manual(values = c(21, 22, 24)) +
    labs(color =  "Ecosystem",
         fill = "Ecosystem",
         shape = "Ecosystem") +
    theme_sleek(base_size = 14) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.85, 0.6))
print(g)

ggsave("./figures/present/map_ocean_entry.pdf", width = 10, height = 5.5)
ggsave("./figures/present/map_ocean_entry.png", width = 10, height = 5.5)



## Fig: NPGO + PDO -----------------------------------------

d.npgo <- data.frame(year = npgo.winter$year,
                     value = npgo.winter$index,
                     var = "NPGO")
d.pdo <- data.frame(year = pdo.winter$year,
                    value = pdo.winter$index,
                    var = "PDO")
dat.ev <- rbind(d.npgo, d.pdo)
dat.ev <- dat.ev[dat.ev$year > 1950 &
                 dat.ev$year < 2012, ]


g <- ggplot(dat.ev) +
    geom_hline(yintercept = 0, col = "grey70", linetype = 2) +
    geom_vline(xintercept = 1976, col = "grey70", linetype = 2) +
    geom_vline(xintercept = 1989, col = "grey70", linetype = 2) +
    geom_line(aes(x = year, y = value), color = "steelblue", size = 1) +
    labs(x = "Year",
         y = "Value") +
    theme_sleek(base_size = 14) +
    facet_wrap( ~ var, ncol = 1)
print(g)
ggsave("./figures/present/npgo_pdo_ts.pdf", width = 7, height = 6)
ggsave("./figures/present/npgo_pdo_ts.png", width = 7, height = 6)



## Fig: sst map --------------------------------------------

sst.anom.sub <- sst.anom[sst.anom$lon <= -120 & sst.anom$lon >= -165 &
                         sst.anom$lat <= 62 & sst.anom$lat >= 48 &
                         sst.anom$month == 5 & sst.anom$year == 1984, ]


    g <- ggplot(dat.map) +
        geom_raster(data = sst.anom.sub, aes(x = lon, y = lat, fill = sst_anom)) +
        geom_map(data = namerica.state.df,
                 map = namerica.state.df,
                 aes(map_id = id),
                 fill = "grey90", color = "grey60", size = 0.1) +
        geom_point(data = fr, aes(x = lon, y = lat), color = "grey40", size = 0.5) +
        geom_point(aes(x = lon, y = lat,
                       color = ocean_region_lab,
                       shape = ocean_region_lab),
                   size = 3, alpha = 1) +
        xlab(expression(paste("Longitude (", degree, "E)"))) +
        ylab(expression(paste("Latitude (", degree, "N)"))) +
        scale_color_manual(values = col.stock, labels = c("Apr-Jul", "May-Aug", "Jun-Sep")) +
        scale_shape_manual(values = c(16, 15, 17), labels = c("Apr-Jul", "May-Aug", "Jun-Sep")) +
        scale_fill_gradientn(colours = col.sst, limits = c(-1.25, 1.25)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(color =  "Average",
             fill = "SST",
             shape = "Average") +
        theme_sleek(base_size = 14)
    print(g)

ggsave("./figures/present/map_sst.pdf", width = 10, height = 5.5)
ggsave("./figures/present/map_sst.png", width = 10, height = 5.5)



## Fig: sst stock time series ------------------------------

g <- ggplot(sock.covar[sock.covar$stock_id == 102, ]) +
    aes(x = brood_yr, y = sst_anom) +
    geom_hline(yintercept = 0, col = "grey70", linetype = 2) +
    geom_line() +
    geom_point() +
    labs(x = "Brood year",
         y = "SST anomaly") +
    theme_sleek(base_size = 12)
print(g)
ggsave("./figures/present/sst_ts_stock.pdf", width = 6, height = 3)
ggsave("./figures/present/sst_ts_stock.png", width = 6, height = 3)



## Fig: era by region --------------------------------------
ss <- plyr::ddply(hbm.gamma.era.stock, .(era, ocean_region, ocean_region_lab,
                                         var, stock), summarize,
                  median = median(value))
sig <- ddply(hbm.gamma.era, .(era, ocean_region_lab, var), summarize,
             median = median(value),
             upper95 = quantile(value, 0.975),
             lower95 = quantile(value, 0.025))
sig$sig <- ((sig$lower95 < 0 & sig$upper95 < 0) |
            (sig$lower95 > 0 & sig$upper95 > 0))
hbm.gamma.era2 <- plyr::join(hbm.gamma.era, sig,
                             by = c("era", "var", "ocean_region_lab"))

col_stock  <- rev(chroma::qpal(7, alpha = 0.4)[c(1, 4, 6)])
col_region <- rev(chroma::qpal(7, luminance = 40)[c(1, 4, 6)])
reg <- as.character(unique(hbm.gamma.era2$ocean_region))

for(i in seq_along(reg)) {

    reg_i  <- reg[i]
    file_i <- paste0("./figures/present/hbm_era_", tolower(reg_i), ".pdf")
    file2_i <- paste0("./figures/present/hbm_era_", tolower(reg_i), ".png")
    dk_i   <- col_region[i]
    lt_i   <- col_stock[i]

    dat_i <- hbm.gamma.era2[hbm.gamma.era2$ocean_region == reg_i, ]
    ss_i  <- ss[ss$ocean_region == reg_i, ]

    g <- ggplot() +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
        geom_violin(data = dat_i[!dat_i$sig, ],
                    aes(era, value, color = era, fill = NULL),
                    adjust = 1.5, width = 0.9) +
        geom_violin(data = dat_i[dat_i$sig, ],
                    aes(era, value, color = era, fill = era),
                    adjust = 1.5, width = 0.9) +
        ylim(-0.7, 0.85) +
        geom_point(data = ss_i, aes(era, median, color = era), size = 0.2) +
        labs(x = "Era", y = "Coefficient") +
        scale_color_manual(values = c(rep(dk_i, 3), "grey40"), guide = FALSE, drop = FALSE) +
        scale_fill_manual(values = c(rep(lt_i, 3), "grey80"), guide = FALSE, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        theme_sleek() +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt")) +
        facet_grid(var ~ ocean_region_lab)
    print(g)

    ggsave(file_i, width = 3.5, height = 5)
    ggsave(file2_i, width = 3.5, height = 5)

}


## Fig: kalman by region -----------------------------------
dat <- hbm.gamma.diff
dat <- ocean_region_lab(dat)
cor.mean <- plyr::ddply(dat, .(var, ocean_region, ocean_region_lab, brood_yr), summarize,
                        cor_mean = mean(gamma),
                        n_cor = length(var))
## should have at least stocks with data for a year to average
cor.mean$cor_mean <- ifelse(cor.mean$n_cor <= 3, NA, cor.mean$cor_mean)


col_stock  <- rev(chroma::qpal(7, alpha = 0.4)[c(1, 4, 6)])
col_region <- rev(chroma::qpal(7, luminance = 40)[c(1, 4, 6)])
reg <- as.character(unique(dat$ocean_region))

for(i in seq_along(reg)) {

    reg_i  <- reg[i]
    file_i <- paste0("./figures/present/hbm_diff_", tolower(reg_i), ".pdf")
    file2_i <- paste0("./figures/present/hbm_diff_", tolower(reg_i), ".png")
    dk_i   <- col_region[i]
    lt_i   <- col_stock[i]

    dat_i <- dat[dat$ocean_region == reg_i, ]
    cor_i <- cor.mean[cor.mean$ocean_region == reg_i, ]

    g <- ggplot(dat_i) +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
        geom_line(aes(x = brood_yr, y = gamma, group = stock),
                  color = lt_i, size = 0.4) +
        geom_line(data = cor_i, aes(x = brood_yr, y = cor_mean),
                  color = dk_i, size = 1.0) +
        labs(x = "Brood year",
             y = "Coefficient") +
        theme_sleek() +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt")) +
        facet_grid(var ~ ocean_region_lab)
    print(g)

    ggsave(file_i, width = 3.5, height = 5)
    ggsave(file2_i, width = 3.5, height = 5)
}


## Fig: kalman proportion ----------------------------------
dat <- hbm.gamma.diff
dat <- ocean_region_lab(dat)
dat.p <- plyr::ddply(dat, .(ocean_region_lab, var, brood_yr), summarize,
                     ocean_region = unique(ocean_region),
                     n_stocks = length(var),
                     n_pos = sum(gamma > 0),
                     n_neg = sum(gamma <= 0),
                     prop_pos = sum(gamma > 0) / length(stock),
                     prop_neg = sum(gamma <= 0) / length(stock))
dat.p$neg <- dat.p$prop_pos - 1

dat.ph <- dat.p[dat.p$n_stocks >= 2, ]
dat.pg <- dat.p[dat.p$n_stocks == 1, ]

reg <- as.character(unique(dat.p$ocean_region))

for(i in seq_along(reg)) {

    reg_i  <- reg[i]
    file_i <- paste0("./figures/present/hbm_prop_", tolower(reg_i), ".pdf")
    file2_i <- paste0("./figures/present/hbm_prop_", tolower(reg_i), ".png")
    dat_ph_i <- dat.ph[dat.ph$ocean_region == reg_i, ]
    dat_pg_i <- dat.pg[dat.pg$ocean_region == reg_i, ]

    g <- ggplot(dat_ph_i) +
        geom_segment(aes(x = brood_yr, xend = brood_yr, y = 0, yend = prop_pos),
                     color = "tomato", size = 1.2) +
        geom_segment(aes(x = brood_yr, xend = brood_yr, y = 0, yend = neg),
                     color = "steelblue", size = 1.2)

    if(nrow(dat_pg_i) > 0) {
        g <- g + geom_segment(data = dat_pg_i,
                              aes(x = brood_yr, xend = brood_yr, y = 0, yend = prop_pos),
                              color = "grey60", size = 1.2) +
            geom_segment(data = dat_pg_i,
                         aes(x = brood_yr, xend = brood_yr, y = 0, yend = neg),
                         color = "grey60", size = 1.2)
    }

    g <- g + geom_hline(yintercept = 0.0, col = "grey50", linetype = 1) +
        labs(x = "Brood year",
             y = "Proportion of stocks") +
        scale_y_continuous(limits = c(-1, 1),
                           breaks = seq(-1, 1, 0.5),
                           labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0")) +
        theme_sleek() +
        ## theme(panel.spacing.y = unit(-0.5, "pt"),
        ##       panel.spacing.x = unit(-0.5, "pt"),
        ##       axis.text.y = element_text(hjust = 4.5),
        ##       axis.title.y = element_text(vjust = 4.5)) +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt")) +
        facet_grid(var ~ ocean_region_lab)

    print(g)

    ggsave(file_i, width = 3.5, height = 5)
    ggsave(file2_i, width = 3.5, height = 5)
}
## Fig: NPGO era -------------------------------------------
ss <- plyr::ddply(hbm.gamma.era.stock, .(era, ocean_region, ocean_region_lab,
                                         var, stock), summarize,
                  median = median(value))
sig <- ddply(hbm.gamma.era, .(era, ocean_region_lab, var), summarize,
             median = median(value),
             upper95 = quantile(value, 0.975),
             lower95 = quantile(value, 0.025))
sig$sig <- ((sig$lower95 < 0 & sig$upper95 < 0) |
            (sig$lower95 > 0 & sig$upper95 > 0))
hbm.gamma.era2 <- plyr::join(hbm.gamma.era, sig,
                             by = c("era", "var", "ocean_region_lab"))

col_stock  <- rev(chroma::qpal(7, alpha = 0.4)[c(1, 4, 6)])
col_region <- rev(chroma::qpal(7, luminance = 40)[c(1, 4, 6)])
reg <- as.character(unique(hbm.gamma.era2$ocean_region))

for(i in seq_along(reg)) {

    reg_i  <- reg[i]
    file_i <- paste0("./figures/present/hbm_era_npgo_", tolower(reg_i), ".png")
    dk_i   <- col_region[i]
    lt_i   <- col_stock[i]

    dat_i <- hbm.gamma.era2[hbm.gamma.era2$ocean_region == reg_i &
                            hbm.gamma.era2$var == "NPGO", ]
    ss_i  <- ss[ss$ocean_region == reg_i & ss$var == "NPGO", ]

    g <- ggplot() +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
        geom_violin(data = dat_i[!dat_i$sig, ],
                    aes(era, value, color = era, fill = NULL),
                    adjust = 1.5, width = 0.9) +
        geom_violin(data = dat_i[dat_i$sig, ],
                    aes(era, value, color = era, fill = era),
                    adjust = 1.5, width = 0.9) +
        ylim(-0.7, 0.85) +
        geom_point(data = ss_i, aes(era, median, color = era), size = 0.2) +
        labs(x = "Era", y = "Coefficient") +
        scale_color_manual(values = c(rep(dk_i, 3), "grey40"), guide = FALSE, drop = FALSE) +
        scale_fill_manual(values = c(rep(lt_i, 3), "grey80"), guide = FALSE, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        theme_sleek()
    print(g)

    ggsave(file_i, width = 4, height = 2.3)

}


## Fig: NPGO kalman ----------------------------------------
dat <- hbm.gamma.diff
dat <- ocean_region_lab(dat)
cor.mean <- plyr::ddply(dat, .(var, ocean_region, ocean_region_lab, brood_yr), summarize,
                        cor_mean = mean(gamma),
                        n_cor = length(var))
## should have at least stocks with data for a year to average
cor.mean$cor_mean <- ifelse(cor.mean$n_cor <= 3, NA, cor.mean$cor_mean)


col_stock  <- rev(chroma::qpal(7, alpha = 0.4)[c(1, 4, 6)])
col_region <- rev(chroma::qpal(7, luminance = 40)[c(1, 4, 6)])
reg <- as.character(unique(dat$ocean_region))

for(i in seq_along(reg)) {

    reg_i  <- reg[i]
    file_i <- paste0("./figures/present/hbm_diff_npgo_", tolower(reg_i), ".png")
    dk_i   <- col_region[i]
    lt_i   <- col_stock[i]

    dat_i <- dat[dat$ocean_region == reg_i & dat$var == "NPGO", ]
    cor_i <- cor.mean[cor.mean$ocean_region == reg_i & cor.mean$var == "NPGO", ]

    g <- ggplot(dat_i) +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
        geom_line(aes(x = brood_yr, y = gamma, group = stock),
                  color = lt_i, size = 0.4) +
        geom_line(data = cor_i, aes(x = brood_yr, y = cor_mean),
                  color = dk_i, size = 1.0) +
        ylim(-0.8, 1) +
        labs(x = "Brood year",
             y = "Coefficient") +
        theme_sleek()
    print(g)

    ggsave(file_i, width = 4, height = 2.3)
}

