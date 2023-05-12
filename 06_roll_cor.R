## Rolling correlations
##
## This script runs the stationary and rolling correlation analysis between
## sockeye productivity and the environmental variables.


if(!dir.exists("./figures/roll_cor/"))
    dir.create("./figures/roll_cor/")


## Stationary correlations ---------------------------------
cor.sst <- cor_stock(sock.covar,
                     x = "rk_resid",
                     y = "sst_anom_stnd",
                     var = "SST",
                     level = 0.90,
                     ar1.correction = TRUE,
                     n.lags = 5)

cor.pdo <- cor_stock(sock.covar,
                     x = "rk_resid",
                     y = "pdo_winter_stnd",
                     var = "PDO",
                     level = 0.90,
                     ar1.correction = TRUE,
                     n.lags = 5)

cor.npgo <- cor_stock(sock.covar,
                      x = "rk_resid",
                      y = "npgo_winter_stnd",
                      var = "NPGO",
                      level = 0.90,
                      ar1.correction = TRUE,
                      n.lags = 5)


save(cor.sst, file = "./output/cor_sst.RData")
save(cor.pdo, file = "./output/cor_pdo.RData")
save(cor.npgo, file = "./output/cor_npgo.RData")


pdf("./figures/roll_cor/cor_stationary.pdf", width = 10, height = 8)

    cor.enviro <- rbind(cor.sst, cor.pdo, cor.npgo)
    g <- xyplot(stock ~ cor | var, groups = ocean_region,
                data = cor.enviro,
                lower = cor.enviro[ , "ci_lower"],
                upper = cor.enviro[ , "ci_upper"],
                layout = c(3, 1),
                xlim = c(-1, 1),
                par.settings = par_mjm(),
                ylab = "",
                xlab = "Correlation",
                auto.key = list(space = "right"),
                panel = panel.superpose,
                panel.groups = function(x, y, subscripts, col.line,
                                        group.number, upper, lower, ...) {
                    panel.abline(v = 0, col = "grey50", lty = 2)
                    upper <- upper[subscripts]
                    lower <- lower[subscripts]
                    panel.xyplot(x, y, ...)
                    panel.segments(x0 = lower, x1 = upper,
                                   y0 = y, y1 = y, col = col.line)
                })
    print(g)

dev.off()



## Rolling correlations ------------------------------------
cor.roll.sst <- cor_roll_stock(sock.covar,
                               x = "rk_resid",
                               y = "sst_anom_stnd",
                               window = 15,
                               var = "SST",
                               level = 0.90,
                               ar1.correction = TRUE,
                               n.lags = 5)

cor.roll.pdo <- cor_roll_stock(sock.covar,
                               x = "rk_resid",
                               y = "pdo_winter_stnd",
                               window = 15,
                               var = "PDO",
                               level = 0.90,
                               ar1.correction = TRUE,
                               n.lags = 5)

cor.roll.npgo <- cor_roll_stock(sock.covar,
                                x = "rk_resid",
                                y = "npgo_winter_stnd",
                                window = 15,
                                var = "NPGO",
                                level = 0.90,
                                ar1.correction = TRUE,
                                n.lags = 5)


save(cor.roll.sst, file = "./output/cor_roll_sst.RData")
save(cor.roll.pdo, file = "./output/cor_roll_pdo.RData")
save(cor.roll.npgo, file = "./output/cor_roll_npgo.RData")


## By stock
pdf("./figures/roll_cor/cor_roll_sst.pdf", width = 14, height = 8)
    plot_roll_ts(cor ~ mid | stock, data = cor.roll.sst,
                 groups = ocean_region,
                 data.mean = cor.sst,
                 ylab = "Correlation",
                 xlab = "Middle year",
                 ylim = c(-1, 1))
dev.off()

pdf("./figures/roll_cor/cor_roll_pdo.pdf", width = 14, height = 8)
    plot_roll_ts(cor ~ mid | stock, data = cor.roll.pdo,
                 groups = ocean_region,
                 data.mean = cor.pdo,
                 ylab = "Correlation",
                 xlab = "Middle year",
                 ylim = c(-1, 1))
dev.off()

pdf("./figures/roll_cor/cor_roll_npgo.pdf", width = 14, height = 8)
    plot_roll_ts(cor ~ mid | stock, data = cor.roll.npgo,
                 groups = ocean_region,
                 data.mean = cor.npgo,
                 ylab = "Correlation",
                 xlab = "Middle year",
                 ylim = c(-1, 1))
dev.off()



## By ocean region
pdf("./figures/roll_cor/cor_roll_sst_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(cor ~ mid | ocean_region, data = cor.roll.sst,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                ylim = c(-1, 1),
                layout = c(1, 3),
                ylab = "Correlation",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()

pdf("./figures/roll_cor/cor_roll_pdo_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(cor ~ mid | ocean_region, data = cor.roll.pdo,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                ylim = c(-1, 1),
                layout = c(1, 3),
                ylab = "Correlation",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()

pdf("./figures/roll_cor/cor_roll_npgo_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(cor ~ mid | ocean_region, data = cor.roll.npgo,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                ylim = c(-1, 1),
                layout = c(1, 3),
                ylab = "Correlation",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()
