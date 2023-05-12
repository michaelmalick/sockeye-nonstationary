## MARSS univariate models
##
## This script fits dynamic linear models (DLM) using
## the KFAS package. The models are formulated as Ricker
## models with time varying parameters.
##
## TODO could run optim over multiple inits

if(!dir.exists("./figures/dlm/"))
    dir.create("./figures/dlm/")


## Ricker: Alpha varying -----------------------------------
dlm.alpha.fits <- lapply(split(sock.covar, sock.covar$stock), function(x) {
        y  <- y_matrix(x, "ln_rps", id.var = "stock", kfas = TRUE)
        x1 <- y_matrix(x, "spawners", id.var = "stock", kfas = TRUE)
        mod1 <- kfas_model(y, x1 = x1, x2 = NULL,
                           Qa = "diagonal and unequal",
                           Qb = "zero",
                           P1 = c(0, 0), P1inf = c(1, 1))
        fit <- kfas_fit(mod1)
        fit
})
save(dlm.alpha.fits, file = "./output/dlm_alpha_fits.RData")
dlm.alpha.lst <- lapply(dlm.alpha.fits, kfas_tidy)
dlm.alpha <- plyr::rbind.fill(dlm.alpha.lst)
save(dlm.alpha, file = "./output/dlm_alpha.RData")

pdf("./figures/dlm/dlm_alpha_region.pdf", width = 9, height = 8)
    g <- xyplot(alpha ~ brood_yr | sub_region,
                data = dlm.alpha, groups = stock, col = "steelblue",
                type = "l", par.settings = par_mjm(),
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "Alpha",
                abline = list(h = 0, col = "grey50", lty = 2))
    print(g)
dev.off()

pdf("./figures/dlm/dlm_alpha_stock.pdf", width = 14, height = 8)
    g <- xyplot(alpha ~ brood_yr | stock,
                data = dlm.alpha, type = "l",
                groups = ocean_region,
                lower = dlm.alpha$alpha - qnorm(0.975) * dlm.alpha$se_alpha,
                upper = dlm.alpha$alpha + qnorm(0.975) * dlm.alpha$se_alpha,
                lwd = 2,
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "Alpha",
                par.settings = par_mjm(),
                panel = function(x, y, subscripts, upper, lower, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.polygon(x = c(x, rev(x)),
                                  y = c(lower[subscripts], rev(upper[subscripts])),
                                  border = NA, col = "grey85")
                    panel.xyplot(x, y, subscripts = subscripts, ...)
                })
    print(g)
dev.off()



## Ricker: SST varying -------------------------------------
dlm.sst.fits <- lapply(split(sock.covar, sock.covar$stock), function(x) {
        y  <- y_matrix(x, "ln_rps", id.var = "stock", kfas = TRUE)
        x1 <- y_matrix(x, "spawners", id.var = "stock", kfas = TRUE)
        x2 <- y_matrix(x, "sst_anom_stnd", id.var = "stock", kfas = TRUE)
        mod1 <- kfas_model(y, x1 = x1, x2 = x2,
                           Qa = "zero",
                           Qb = "zero",
                           Qg = "diagonal and unequal",
                           P1 = c(0, 0, 0), P1inf = c(1, 1, 1))
        fit <- kfas_fit(mod1, tol = .Machine$double.eps^0.5)
        fit
})
save(dlm.sst.fits, file = "./output/dlm_sst_fits.RData")
dlm.sst.lst <- lapply(dlm.sst.fits, kfas_tidy)
dlm.sst <- plyr::rbind.fill(dlm.sst.lst)
save(dlm.sst, file = "./output/dlm_sst.RData")

pdf("./figures/dlm/dlm_sst_region.pdf", width = 6, height = 8)
    g <- xyplot(gamma ~ brood_yr | ocean_region,
                data = dlm.sst, groups = stock, col = "steelblue",
                type = "l", par.settings = par_mjm(),
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "SST",
                layout = c(1, 3),
                abline = list(h = 0, col = "grey50", lty = 2))
    print(g)
dev.off()

pdf("./figures/dlm/dlm_sst_stock.pdf", width = 14, height = 8)
    g <- xyplot(gamma ~ brood_yr | stock,
                data = dlm.sst, type = "l",
                groups = ocean_region,
                lower = dlm.sst$gamma - qnorm(0.975) * dlm.sst$se_gamma,
                upper = dlm.sst$gamma + qnorm(0.975) * dlm.sst$se_gamma,
                lwd = 2,
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "SST",
                par.settings = par_mjm(),
                panel = function(x, y, subscripts, upper, lower, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.polygon(x = c(x, rev(x)),
                                  y = c(lower[subscripts], rev(upper[subscripts])),
                                  border = NA, col = "grey85")
                    panel.xyplot(x, y, subscripts = subscripts, ...)
                })
    print(g)
dev.off()



## Ricker: NPGO varying ------------------------------------
dlm.npgo.fits <- lapply(split(sock.covar, sock.covar$stock), function(x) {
        y  <- y_matrix(x, "ln_rps", id.var = "stock", kfas = TRUE)
        x1 <- y_matrix(x, "spawners", id.var = "stock", kfas = TRUE)
        x2 <- y_matrix(x, "npgo_winter_stnd", id.var = "stock", kfas = TRUE)
        mod1 <- kfas_model(y, x1 = x1, x2 = x2,
                           Qa = "zero",
                           Qb = "zero",
                           Qg = "diagonal and unequal",
                           P1 = c(0, 0, 0), P1inf = c(1, 1, 1))
        fit <- kfas_fit(mod1, tol = .Machine$double.eps^0.5)
        fit
})
save(dlm.npgo.fits, file = "./output/dlm_npgo_fits.RData")
dlm.npgo.lst <- lapply(dlm.npgo.fits, kfas_tidy)
dlm.npgo <- plyr::rbind.fill(dlm.npgo.lst)
save(dlm.npgo, file = "./output/dlm_npgo.RData")

pdf("./figures/dlm/dlm_npgo_region.pdf", width = 6, height = 8)
    g <- xyplot(gamma ~ brood_yr | ocean_region,
                data = dlm.npgo, groups = stock, col = "steelblue",
                type = "l", par.settings = par_mjm(),
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "NPGO",
                layout = c(1, 3),
                ylim = c(-2, 1.5),
                abline = list(h = 0, col = "grey50", lty = 2))
    print(g)
dev.off()

pdf("./figures/dlm/dlm_npgo_stock.pdf", width = 14, height = 8)
    g <- xyplot(gamma ~ brood_yr | stock,
                data = dlm.npgo, type = "l",
                groups = ocean_region,
                lower = dlm.npgo$gamma - qnorm(0.975) * dlm.npgo$se_gamma,
                upper = dlm.npgo$gamma + qnorm(0.975) * dlm.npgo$se_gamma,
                lwd = 2,
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "NPGO",
                ylim = c(-2, 1.5),
                par.settings = par_mjm(),
                panel = function(x, y, subscripts, upper, lower, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.polygon(x = c(x, rev(x)),
                                  y = c(lower[subscripts], rev(upper[subscripts])),
                                  border = NA, col = "grey85")
                    panel.xyplot(x, y, subscripts = subscripts, ...)
                })
    print(g)
dev.off()


## Ricker: PDO varying -------------------------------------
dlm.pdo.fits <- lapply(split(sock.covar, sock.covar$stock), function(x) {
        y  <- y_matrix(x, "ln_rps", id.var = "stock", kfas = TRUE)
        x1 <- y_matrix(x, "spawners", id.var = "stock", kfas = TRUE)
        x2 <- y_matrix(x, "pdo_winter_stnd", id.var = "stock", kfas = TRUE)
        mod1 <- kfas_model(y, x1 = x1, x2 = x2,
                           Qa = "zero",
                           Qb = "zero",
                           Qg = "diagonal and unequal",
                           P1 = c(0, 0, 0), P1inf = c(1, 1, 1))
        fit <- kfas_fit(mod1, tol = .Machine$double.eps^0.5)
        fit
})
save(dlm.pdo.fits, file = "./output/dlm_pdo_fits.RData")
dlm.pdo.lst <- lapply(dlm.pdo.fits, kfas_tidy)
dlm.pdo <- plyr::rbind.fill(dlm.pdo.lst)
save(dlm.pdo, file = "./output/dlm_pdo.RData")

pdf("./figures/dlm/dlm_pdo_region.pdf", width = 6, height = 8)
    g <- xyplot(gamma ~ brood_yr | ocean_region,
                data = dlm.pdo, groups = stock, col = "steelblue",
                type = "l", par.settings = par_mjm(),
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "PDO",
                layout = c(1, 3),
                abline = list(h = 0, col = "grey50", lty = 2))
    print(g)
dev.off()

pdf("./figures/dlm/dlm_pdo_stock.pdf", width = 14, height = 8)
    g <- xyplot(gamma ~ brood_yr | stock,
                data = dlm.pdo, type = "l",
                groups = ocean_region,
                lower = dlm.pdo$gamma - qnorm(0.975) * dlm.pdo$se_gamma,
                upper = dlm.pdo$gamma + qnorm(0.975) * dlm.pdo$se_gamma,
                lwd = 2,
                ylab = "Coefficient",
                xlab = "Brood year",
                main = "PDO",
                par.settings = par_mjm(),
                panel = function(x, y, subscripts, upper, lower, ...) {
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    panel.polygon(x = c(x, rev(x)),
                                  y = c(lower[subscripts], rev(upper[subscripts])),
                                  border = NA, col = "grey85")
                    panel.xyplot(x, y, subscripts = subscripts, ...)
                })
    print(g)
dev.off()
