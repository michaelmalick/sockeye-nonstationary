## Rolling Ricker models
##
## This script runs the stationary and rolling Ricker model analysis.


if(!dir.exists("./figures/roll_rk/"))
    dir.create("./figures/roll_rk/")


## Stationary Ricker models --------------------------------
rk.sst.full   <- rk_stock(sock.covar, ln_rps ~ spawners + sst_anom_stnd,
                          level = 0.90,
                          correlation = corAR1(),
                          var = "SST")
rk.pdo.full   <- rk_stock(sock.covar, ln_rps ~ spawners + pdo_winter_stnd,
                          level = 0.90,
                          correlation = corAR1(),
                          var = "PDO")
rk.npgo.full  <- rk_stock(sock.covar, ln_rps ~ spawners + npgo_winter_stnd,
                          level = 0.90,
                          correlation = corAR1(),
                          var = "NPGO")
rk.alpha.full <- rk_stock(sock.covar, ln_rps ~ spawners,
                          level = 0.90,
                          correlation = corAR1(),
                          var = NULL)

rk.sst   <- rk.sst.full[["sst_anom_stnd"]]
rk.pdo   <- rk.pdo.full[["pdo_winter_stnd"]]
rk.npgo  <- rk.npgo.full[["npgo_winter_stnd"]]
rk.npgo  <- rk.npgo.full[["npgo_winter_stnd"]]
rk.alpha <- rk.alpha.full[["(Intercept)"]]
rk.beta  <- rk.alpha.full[["spawners"]]


save(rk.sst.full, file = "./output/rk_sst_full.RData")
save(rk.pdo.full, file = "./output/rk_pdo_full.RData")
save(rk.npgo.full, file = "./output/rk_npgo_full.RData")
save(rk.alpha.full, file = "./output/rk_alpha_full.RData")

save(rk.sst, file = "./output/rk_sst.RData")
save(rk.pdo, file = "./output/rk_pdo.RData")
save(rk.npgo, file = "./output/rk_npgo.RData")
save(rk.alpha, file = "./output/rk_alpha.RData")
save(rk.beta, file = "./output/rk_beta.RData")


pdf("./figures/roll_rk/rk_stationary.pdf", width = 10, height = 8)

    rk.enviro <- rbind(rk.sst, rk.pdo, rk.npgo)
    g <- xyplot(stock ~ estimate | var, groups = ocean_region,
                data = rk.enviro,
                lower = rk.enviro[ , "ci_lower"],
                upper = rk.enviro[ , "ci_upper"],
                layout = c(3, 1),
                xlim = c(min(rk.enviro[ , "ci_lower"]),
                         max(rk.enviro[ , "ci_upper"])),
                par.settings = par_mjm(),
                ylab = "",
                xlab = "Coefficient",
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



## Rolling Ricker models -----------------------------------
rk.roll.sst.full <- rk_roll_stock(sock.covar, ln_rps ~ spawners + sst_anom_stnd,
                                  window = 15,
                                  correlation = corAR1(),
                                  level = 0.90)
rk.roll.sst <- rk_roll_getparam(rk.roll.sst.full,
                                param = "sst_anom_stnd",
                                var = "SST")


rk.roll.pdo.full <- rk_roll_stock(sock.covar, ln_rps ~ spawners + pdo_winter_stnd,
                                  window = 15,
                                  correlation = corAR1(),
                                  level = 0.90)
rk.roll.pdo <- rk_roll_getparam(rk.roll.pdo.full,
                                param = "pdo_winter_stnd",
                                var = "PDO")


rk.roll.npgo.full <- rk_roll_stock(sock.covar, ln_rps ~ spawners + npgo_winter_stnd,
                                   window = 15,
                                   correlation = corAR1(),
                                   level = 0.90)
rk.roll.npgo <- rk_roll_getparam(rk.roll.npgo.full,
                                 param = "npgo_winter_stnd",
                                 var = "NPGO")

rk.roll.alpha.full <- rk_roll_stock(sock.covar, ln_rps ~ spawners,
                                    window = 15,
                                    correlation = corAR1(),
                                    level = 0.90)
rk.roll.alpha <- rk_roll_getparam(rk.roll.alpha.full,
                                  param = "(Intercept)",
                                  var = "Alpha")
rk.roll.beta <- rk_roll_getparam(rk.roll.alpha.full,
                                 param = "spawners",
                                 var = "Beta")


save(rk.roll.sst.full, file = "./output/rk_roll_sst_full.RData")
save(rk.roll.pdo.full, file = "./output/rk_roll_pdo_full.RData")
save(rk.roll.npgo.full, file = "./output/rk_roll_npgo_full.RData")
save(rk.roll.alpha.full, file = "./output/rk_roll_alpha_full.RData")

save(rk.roll.sst, file = "./output/rk_roll_sst.RData")
save(rk.roll.pdo, file = "./output/rk_roll_pdo.RData")
save(rk.roll.npgo, file = "./output/rk_roll_npgo.RData")
save(rk.roll.alpha, file = "./output/rk_roll_alpha.RData")
save(rk.roll.beta, file = "./output/rk_roll_beta.RData")



## By stock
pdf("./figures/roll_rk/rk_roll_sst.pdf", width = 14, height = 8)
    plot_roll_ts(estimate ~ mid | stock, data = rk.roll.sst,
                 groups = ocean_region,
                 data.mean = rk.sst,
                 ylab = "Coefficient",
                 xlab = "Middle year",
                 ylim = range(c(rk.roll.sst[ , "ci_lower"],
                                rk.roll.sst[ , "ci_upper"])))
dev.off()

pdf("./figures/roll_rk/rk_roll_pdo.pdf", width = 14, height = 8)
    plot_roll_ts(estimate ~ mid | stock, data = rk.roll.pdo,
                 groups = ocean_region,
                 data.mean = rk.pdo,
                 ylab = "Coefficient",
                 xlab = "Middle year",
                 ylim = range(c(rk.roll.pdo[ , "ci_lower"],
                                rk.roll.pdo[ , "ci_upper"])))
dev.off()

pdf("./figures/roll_rk/rk_roll_npgo.pdf", width = 14, height = 8)
    plot_roll_ts(estimate ~ mid | stock, data = rk.roll.npgo,
                 groups = ocean_region,
                 data.mean = rk.npgo,
                 ylab = "Coefficient",
                 xlab = "Middle year",
                 ylim = range(c(rk.roll.npgo[ , "ci_lower"],
                                rk.roll.npgo[ , "ci_upper"])))
dev.off()

pdf("./figures/roll_rk/rk_roll_alpha.pdf", width = 14, height = 8)
    plot_roll_ts(estimate ~ mid | stock, data = rk.roll.alpha,
                 groups = ocean_region,
                 data.mean = rk.alpha,
                 ylab = "Coefficient",
                 xlab = "Middle year",
                 ylim = range(c(rk.roll.alpha[ , "ci_lower"],
                                rk.roll.alpha[ , "ci_upper"])))
dev.off()

pdf("./figures/roll_rk/rk_roll_beta.pdf", width = 14, height = 8)
    plot_roll_ts(estimate ~ mid | stock, data = rk.roll.beta,
                 groups = ocean_region,
                 data.mean = rk.beta,
                 ylab = "Coefficient",
                 xlab = "Middle year",
                 scales = list(y = list(relation = "free")),
                 ylim = range(c(rk.roll.beta[ , "ci_lower"],
                                rk.roll.beta[ , "ci_upper"])))
dev.off()


## By ocean region
pdf("./figures/roll_rk/rk_roll_sst_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(estimate ~ mid | ocean_region, data = rk.roll.sst,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                layout = c(1, 3),
                ylab = "Coefficient",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()

pdf("./figures/roll_rk/rk_roll_pdo_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(estimate ~ mid | ocean_region, data = rk.roll.pdo,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                layout = c(1, 3),
                ylab = "Coefficient",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()

pdf("./figures/roll_rk/rk_roll_npgo_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(estimate ~ mid | ocean_region, data = rk.roll.npgo,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                layout = c(1, 3),
                ylab = "Coefficient",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()

pdf("./figures/roll_rk/rk_roll_alpha_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(estimate ~ mid | ocean_region, data = rk.roll.alpha,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                layout = c(1, 3),
                ylab = "Coefficient",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()

pdf("./figures/roll_rk/rk_roll_beta_oceanregion.pdf", width = 6, height = 8)
    g <- xyplot(estimate ~ mid | ocean_region, data = rk.roll.beta,
                groups = stock,
                type = "l",
                lwd = 2,
                col = "steelblue",
                layout = c(1, 3),
                ylab = "Coefficient",
                xlab = "Middle year",
                abline = list(h = 0, col = "grey50", lty = 2),
                par.settings = par_mjm())
    print(g)
dev.off()
