## Exploratory data graphics
##
## This script creates exploratory data graphics for the sockeye and
## environmental data sets.


if(!dir.exists("./figures/data_explore/"))
    dir.create("./figures/data_explore/")


## Plot monthly SST anomalies ------------------------------

## This code takes a while to run and is only marginally useful.
## It is commented out to decrease the runtime of reproducing the
## project during development.

# cols <- chroma::dpal(500, hue = c(240, 0), chroma = 70, power = 1.0)
# at <- seq(-4.5, 4.5, 0.5)
# if(!dir.exists("./figures/data_explore/sst-anomaly-maps/"))
#     dir.create("./figures/data_explore/sst-anomaly-maps/")
# sst_map(data = sst.anom,
#         plot.dir = "./figures/data_explore/sst-anomaly-maps/",
#         progress.bar = TRUE,
#         col.regions = cols,
#         par.settings = par_mjm(),
#         at = at)



## Map ocean entry -----------------------------------------
pdf("./figures/data_explore/sock_map_ocean_entry.pdf", width = 7, height = 6)

    plot(sock.info$lon, sock.info$lat,
         xlim = c(-170, -120),
         ylim = c(45, 70),
         main = "Unique ocean entry locations",
         type = "n",
         las = 1,
         xlab = "Longitude",
         ylab = "Latitude")
    plot(countriesLow, add = TRUE, border = "grey50")
    points(sock.info$lon, sock.info$lat, pch = 16, col = "red3")

dev.off()



## Productivity time series --------------------------------
pdf("./figures/data_explore/productivity_by_stock.pdf", width = 19, height = 9)

    g <- xyplot(rps ~ brood_yr | stock, data = sock.covar,
                type = "l",
                scales = list(y = list(relation = "free")),
                ylab = "Productivity (R/S)",
                xlab = "Brood year",
                main = "Productivity by stock",
                par.settings = par_mjm(),
                panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                })
    print(g)

dev.off()

pdf("./figures/data_explore/productivity_by_stock_stnd.pdf", width = 19, height = 9)

    g <- xyplot(rps_stnd ~ brood_yr | stock, data = sock.covar,
                type = "l",
                scales = list(y = list(relation = "free")),
                ylab = "Standardized productivity",
                xlab = "Brood year",
                main = "Standardized productivity by stock",
                par.settings = par_mjm(),
                panel = function(x, y, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.xyplot(x, y, ...)
                })
    print(g)

dev.off()



## Spawner/Recruit -----------------------------------------
pdf("./figures/data_explore/spawner_recruit.pdf", width = 19, height = 9)

    g <- xyplot(recruits ~ spawners | stock, data = sock.covar,
                type = "p",
                pch = 16, cex = 0.5,
                scales = list(relation = "free"),
                ylab = "Recruits",
                xlab = "Spawners",
                main = "Recruits vs. spawners w/ fitted Ricker curve",
                par.settings = par_mjm(),
                panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                    rk <- lm(log(y/x) ~ x)
                    a  <- exp(coef(rk)[1])
                    b  <- coef(rk)[2]
                    panel.curve(a*x*exp(b*x),
                                from = min(x, na.rm = TRUE),
                                to = max(x, na.rm = TRUE))
                })
    print(g)

dev.off()


pdf("./figures/data_explore/spawner_recruit_ln.pdf", width = 19, height = 9)

    g <- xyplot(ln_rps ~ spawners | stock, data = sock.covar,
                type = "p",
                pch = 16, cex = 0.5,
                scales = list(relation = "free"),
                ylab = "ln(R/S)",
                xlab = "Spawners",
                main = "ln(R/S) vs. spawners w/ fitted Ricker line",
                par.settings = par_mjm(),
                panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                    panel.abline(lm(y ~ x))
                })
    print(g)

dev.off()



## SST indices ---------------------------------------------
pdf("./figures/data_explore/sst_raw_stock.pdf", width = 19, height = 9)

    g <- xyplot(sst_raw ~ brood_yr | stock, data = sock.covar,
                type = "l",
                par.settings = par_mjm(),
                xlab = "Brood year",
                ylab = "SST raw index",
                main = "SST raw index by stock",
                panel = function(x, y, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.loess(x, y, col = "red2", lwd = 2)
                    panel.xyplot(x, y, ...)
                })
    print(g)

dev.off()


pdf("./figures/data_explore/sst_anom_stock.pdf", width = 19, height = 9)

    g <- xyplot(sst_anom ~ brood_yr | stock, data = sock.covar,
                type = "l",
                par.settings = par_mjm(),
                xlab = "Brood year",
                ylab = "SST anomaly index",
                main = "SST anomaly index by stock",
                panel = function(x, y, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.loess(x, y, col = "red2", lwd = 2)
                    panel.xyplot(x, y, ...)
                })
    print(g)

dev.off()



## PDO/NPGO indices ----------------------------------------
pdf("./figures/data_explore/pdo_npgo_ts.pdf", width = 10, height = 6)

    tmp <- rbind(pdo.winter, pdo.spring, pdo.summer,
                 npgo.winter, npgo.spring, npgo.summer)
    g <- xyplot(index ~ year | type, data = tmp,
                groups = months,
                auto.key = list(space = "right"),
                type = "l",
                layout = c(1, 2),
                par.settings = par_mjm(),
                xlab = "Year",
                ylab = "Index",
                main = "PDO and NPGO indices",
                panel = function(x, y, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.xyplot(x, y, ...)
                })
    print(g)

dev.off()



## PDO/NPGO splom  -----------------------------------------
pdf("./figures/data_explore/pdo_npgo_splom.pdf", width = 10, height = 10)

    tmp <- data.frame(pdo.winter = pdo.winter$index,
                      pdo.spring = pdo.spring$index[-1],
                      pdo.summer = pdo.summer$index[-1],
                      npgo.winter = npgo.winter$index,
                      npgo.spring = npgo.spring$index[-1],
                      npgo.summer = npgo.summer$index[-1])
    g <- splom(~ tmp, col = "grey50", pch = 16, cex = 0.8,
               par.settings = par_mjm(),
               panel = function(x, y, ...) {
                   panel.splom(x, y, ...)
                   panel.loess(x, y, lwd = 2)
                   panel.text(min(x), max(y),
                              paste("r =", round(cor(x, y), 2)),
                              cex = 0.8, adj = 0, col = "red2")
               })
    print(g)

dev.off()
