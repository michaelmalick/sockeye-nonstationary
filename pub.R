## Publication figures + tables


if(!dir.exists("./pub/figures/"))
    dir.create("./pub/figures/")
if(!dir.exists("./pub/tables/"))
    dir.create("./pub/tables/")

min(sock.info$yr_start)
max(sock.info$yr_end)
range(sock.info$n_years)
median(sock.info$n_years)

labs <- expand.grid(ocean_region = c("WC", "GOA", "BS"), var = c("NPGO", "PDO", "SST"))
labs <- ocean_region_lab(labs)
labs$label <- paste0("(", letters[1:nrow(labs)], ")")
labs$era <- "early"


## Table: data info ----------------------------------------
dat.info <- ocean_region_lab(sock.info)
dat.info$stock <- as.character(dat.info$stock)
dat.info$brood_yrs <- paste0(dat.info$yr_start, "--", dat.info$yr_end)
dat.info.t <- data.frame("No" = 1:nrow(dat.info),
                         Ecosystem = dat.info$ocean_region_lab,
                         Region = dat.info$region,
                         Stock = dat.info$stock,
                         N = dat.info$n_years,
                         Brood.years = dat.info$brood_yrs)
names(dat.info.t) <- gsub(".", " ", names(dat.info.t), fixed = TRUE)
# dat.info.t$Ecosystem <- replace_dup(dat.info.t$Ecosystem)
# dat.info.t$Region <- replace_dup(dat.info.t$Region)
dat.info.t

write.table(dat.info.t, file = "./pub/tables/data_info.csv",
            sep = ",", row.names = FALSE)


## Table: era model coef -----------------------------------
era.tab <- plyr::ddply(hbm.gamma.era, .(var, ocean_region_lab, era), summarize,
                       lower = quantile(value, probs = 0.025),
                       median = median(value),
                       upper = quantile(value, probs = 0.975))
era.tab$median <- formatC(era.tab$median, digits = 2, format = "f")
era.tab$lower <- formatC(era.tab$lower, digits = 2, format = "f")
era.tab$upper <- formatC(era.tab$upper, digits = 2, format = "f")
names(era.tab) <- c("Variable", "Ecosystem", "Period",
                    "Lower 95% CI", "Median", "Upper 95% CI")
write.table(era.tab, file = "./pub/tables/era_coef.csv",
            sep = ",", row.names = FALSE)


## S2
erat.tab.s2 <- plyr::ddply(hbm.gamma.era.s2, .(var, ocean_region_lab, era), summarize,
                           lower = quantile(value, probs = 0.025),
                           median = median(value),
                           upper = quantile(value, probs = 0.975))
erat.tab.s2$median <- formatC(erat.tab.s2$median, digits = 2, format = "f")
erat.tab.s2$lower <- formatC(erat.tab.s2$lower, digits = 2, format = "f")
erat.tab.s2$upper <- formatC(erat.tab.s2$upper, digits = 2, format = "f")
names(erat.tab.s2) <- c("Variable", "Ecosystem", "Period",
                        "Lower 95% CI", "Median", "Upper 95% CI")
write.table(erat.tab.s2, file = "./pub/tables/era_coef_s2.csv",
            sep = ",", row.names = FALSE)



## Fig 1: map + covars -------------------------------------
## Map data
cl <- rnaturalearth::ne_states(country = c("United States of America",
                                           "Mexico",
                                           "Canada"))
namerica.state.sp <- sp::spTransform(cl, sp::CRS("+init=epsg:4326"))
namerica.state.df  <- suppressMessages(ggplot2::fortify(namerica.state.sp))
dat.map <- ocean_region_lab(sock.info)

col.stock <- RColorBrewer::brewer.pal(8, "Dark2")
col.light <- RColorBrewer::brewer.pal(8, "Set2")
dat.map <- ocean_region_lab(sock.info)
dm <- plyr::ddply(dat.map, .(ocean_region_lab, lon, lat), summarize,
                  n = length(stock_no),
                  no_min = min(stock_no),
                  no_max = max(stock_no))
dm$label <- ifelse(dm$no_min == dm$no_max, dm$no_min, paste0(dm$no_min, "-", dm$no_max))
dm$x <- dm$lon + 0.5
dm$y <- dm$lat + 0.5
dm$y[dm$no_min == 26] <- dm$lat[dm$no_min == 26] - 0.45
dm$y[dm$no_min == 27] <- dm$lat[dm$no_min == 27]
dm$x[dm$no_min == 27] <- dm$lon[dm$no_min == 27] + 0.8
dm$y[dm$no_min == 28] <- dm$lat[dm$no_min == 28]
dm$x[dm$no_min == 28] <- dm$lon[dm$no_min == 28] - 0.8
dm$y[dm$no_min == 29] <- dm$lat[dm$no_min == 29] - 0.5
dm$x[dm$no_min == 29] <- dm$lon[dm$no_min == 29]
dm$x[dm$no_min == 35] <- dm$lon[dm$no_min == 35] - 0.45
dm$x[dm$no_min == 37] <- dm$lon[dm$no_min == 37] - 0.6
dm$y[dm$no_min == 37] <- dm$lat[dm$no_min == 37] + 0.2
dm$y[dm$no_min == 38] <- dm$lat[dm$no_min == 38] + 0
dm$x[dm$no_min == 38] <- dm$lon[dm$no_min == 38] + 1
dm$label[dm$no_min == 39] <- ""
dm$label[dm$no_min == 40] <- "39-40"
dm$x[dm$no_min == 40] <- dm$lon[dm$no_min == 40] + 1.2
dm$label[dm$no_min == 41] <- ""
dm$label[dm$no_min == 42] <- ""
dm$label[dm$no_min == 43] <- "41-43"
dm$y[dm$no_min == 43] <- dm$lat[dm$no_min == 43] + 0.8

## Salmon map
g.map <- ggplot(dm) +
    geom_map(data = namerica.state.df,
                map = namerica.state.df,
                aes(map_id = id),
                fill = "grey90", color = "grey60", size = 0.1) +
    geom_text(aes(x = x, y = y, label = label), size = 2) +
    geom_point(aes(x = lon, y = lat,
                    fill = ocean_region_lab,
                    shape = ocean_region_lab),
                size = 3) +
    geom_point(aes(x = lon, y = lat,
                    shape = ocean_region_lab,
                    color = ocean_region_lab),
                size = 3) +
    xlab(expression(paste("Longitude (", degree, "E)"))) +
    ylab(expression(paste("Latitude (", degree, "N)"))) +
    scale_color_manual(values = col.stock) +
    scale_fill_manual(values = col.light) +
    scale_shape_manual(values = c(21, 22, 24)) +
    labs(color =  "",
            fill = "",
            shape = "") +
    theme_sleek() +
    theme(legend.justification = c(0, 0),
            legend.position = c(0.25, 0.1))
print(g.map)

## Inset map
g.in <- ggplot() +
    geom_map(data = namerica.state.df,
             map = namerica.state.df,
             aes(map_id = id),
             fill = "grey50", color = "grey50", size = 0.1) +
    geom_rect(data = mtcars[1,],
              aes(xmin = -165, xmax = -121, ymin = 48, ymax = 62),
              color = "grey20", fill = NA, size = 0.8,
              inherit.aes = FALSE, linetype = 1) +
    ylim(10, 75) +
    xlim(-170, -60) +
    theme_sleek() +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
print(g.in)



## Covariate time series
pdo <- read.csv("./data/pdo.csv")
npgo <- read.csv("./data/npgo.csv")
names(pdo)[names(pdo) == "pdo"] <- "value"
names(npgo)[names(npgo) == "npgo"] <- "value"
pdo$var <- "PDO"
npgo$var <- "NPGO"
sst.anom.sub <- sst.anom[sst.anom$id == 474, ]
sst <- data.frame(year = sst.anom.sub$year,
                  month = sst.anom.sub$month,
                  value = sst.anom.sub$sst_anom)
sst$var <- "SST"
sst$time <- 1:nrow(sst)
pdo$time <- 1:nrow(pdo)
npgo$time <- 1:nrow(npgo)
dat.ev <- rbind(npgo, pdo, sst)
dat.ev <- dat.ev[dat.ev$year > 1950 &
                 dat.ev$year < 2013, ]
dat.ev$date <- as.Date(paste0(dat.ev$year, "-", dat.ev$month, "-", "1"))
labs.ev <- data.frame(var = c("NPGO", "PDO", "SST"), year = 1951)
#
ma <- function(x, n = 5) {filter(x, rep(1 / n, n), sides = 2)}
dat.ev <- plyr::ddply(dat.ev, .(var), transform,
                      ma = as.vector(ma(value, 24)))
#
g.ev <- ggplot(dat.ev) +
    geom_vline(xintercept = as.Date("1977-06-01"), color = "grey50", linetype = 2) +
    geom_vline(xintercept = as.Date("1988-06-01"), color = "grey50", linetype = 2) +
    geom_line(aes(x = date, y = value, color = value), color = "grey50", size = 0.3) +
    geom_path(aes(x = date, y = ma, color = ma), size = 1, lineend = "round", na.rm = TRUE) +
    #
    ## geom_path(aes(x = date, y = value, color = value), lineend = "round") +
    ## geom_line(aes(x = date, y = ma), na.rm = TRUE, color = "grey50") +
    scale_colour_gradient2(low = "steelblue",
                           mid = "grey70",
                           high = "tomato",
                           midpoint = 0,
                           space = "Lab",
                           na.value = "green",
                           guide = FALSE) +
    geom_text(data = labs.ev, aes(x = as.Date("1952-01-01"), y = Inf, label = var),
              hjust = 0, vjust = 2.1, color = "grey30", size = 2.5) +
    scale_x_date(expand = c(0, 0)) +
    facet_wrap( ~ var, ncol = 1) +
    theme_sleek() +
    labs(x = "Year",
         y = "Index value") +
    theme(strip.text.x = element_blank(),
          panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt"))
print(g.ev)


## Add inset map to main map
g.map <- g.map +
    annotation_custom(ggplotGrob(g.in),
                      xmin = -132, xmax = -121.2, ymin = 55, ymax = 62)
print(g.map)



## Save full figure
pdf("./pub/figures/map_ocean_entry.pdf", width = 7, height = 7)
g <- gridExtra::grid.arrange(g.map, g.ev, layout_matrix = rbind(c(1,1,1),
                                                                c(1,1,1),
                                                                c(2,2,2),
                                                                c(2,2,2)))
dev.off()

jpeg("./pub/figures/map_ocean_entry.jpg", width = 7, height = 7, units = "in", res = 500)
g <- gridExtra::grid.arrange(g.map, g.ev, layout_matrix = rbind(c(1,1,1),
                                                                c(1,1,1),
                                                                c(2,2,2),
                                                                c(2,2,2)))
dev.off()

postscript("./pub/figures/1_map_ocean_entry.eps", width = 7, height = 7,
           horizontal = FALSE, onefile = FALSE, paper = "special")
g <- gridExtra::grid.arrange(g.map, g.ev, layout_matrix = rbind(c(1,1,1),
                                                                c(1,1,1),
                                                                c(2,2,2),
                                                                c(2,2,2)))
dev.off()



## Fig: gamma scenarios ------------------------------------
N <- 50
# seed <- sample(1:50000, 1)
seed <- 12744
set.seed(seed)


## Simulate gamma series
s.const <- rep(0.5, N)
s.era <- ts_step(N, k = 2,
                 up_mean = 0.8,
                 down_mean = -0.8,
                 up_sd = 0,
                 down_sd = 0)
s.rw <- ts_ar1(N, x_0 = 0, rho = 0.8, sd = 0.2, mean = 0)
s.line <- seq(0.8, -0.8, length.out = N) + s.rw
# plot(0, type = "n", ylim = c(-2, 2), xlim = c(1, N))
# lines(s.const, type = "o", col = 2)
# lines(s.era, type = "o", col = 3)
# lines(s.line, type = "o", col = 5)


## Simulate response
d1 <- sim_ss(n_years = N, sigma = 0.5, gamma = s.const)
d2 <- sim_ss(n_years = N, sigma = 0.5, gamma = s.era)
d3 <- sim_ss(n_years = N, sigma = 0.5, gamma = s.line)
# plot(d1$data$y, d1$data$x2)
# plot(d2$data$y, d2$data$x2)
# plot(d3$data$y, d3$data$x2)

## Combine datasets
col.bi <- c(rep("tomato", N/2), rep("steelblue", N/2))
pch.bi <- c(rep(16, N/2), rep(17, N/2))
lst <- list(data.frame(x = 1:N, y = s.const,
                       col = col.bi, pch = pch.bi,
                       stringsAsFactors = FALSE),
            data.frame(x = d1$data$x2, y = d1$data$y,
                       col = col.bi, pch = pch.bi,
                       stringsAsFactors = FALSE),
            data.frame(x = 1:N, y = s.era,
                       col = col.bi, pch = pch.bi,
                       stringsAsFactors = FALSE),
            data.frame(x = d2$data$x2, y = d2$data$y,
                       col = col.bi, pch = pch.bi,
                       stringsAsFactors = FALSE),
            data.frame(x = 1:N, y = s.line,
                       col = col.bi, pch = pch.bi,
                       stringsAsFactors = FALSE),
            data.frame(x = d3$data$x2, y = d3$data$y,
                       col = col.bi, pch = pch.bi,
                       stringsAsFactors = FALSE))

pdf("./pub/figures/gamma_scenario.pdf", width = 6, height = 5)
    par(mfrow = c(3, 2),
        mar = c(0, 0, 0, 0),
        oma = c(4.0, 4.2, 0.5, 4.2),
        ps = 12)
    for(i in 1:6) {
        x <- lst[[i]]$x
        y <- lst[[i]]$y
        if(i %in% c(1, 3, 5)) {
            xlim <- c(1, 50)
            ylim <- c(-0.9, 0.9)
        }
        if(i %in% c(2, 4, 6)) {
            xlim <- c(-1.3, 1.3)
            ylim <- range(y)
        }
        plot(0, type = "n",
             xlim = xlim,
             ylim = ylim,
             xlab = "",
             ylab = "",
             axes = FALSE)
        box(col = "grey50", lwd = 0.7)
        malick::add_label(paste0("(", letters[i], ")"),
                          xfrac = 0.9, col = "grey25")
        corners <- par("usr")
        if(i %in% 5:6)
            axis(side = 1, lwd = 0, lwd.tick = 1,
                 col = "grey50", col.axis = "grey25",
                 tcl = -0.2, padj = -1)
        if(i %in% c(1, 3, 5))
            axis(side = 2, lwd = 0, lwd.tick = 1, las = 1,
                 col = "grey50", col.axis = "grey25",
                 tcl = -0.2, hadj = 0.6)
        if(i %in% c(2, 4, 6))
            axis(side = 4, lwd = 0, lwd.tick = 1, las = 1,
                 col = "grey50", col.axis = "grey25",
                 tcl = -0.2, hadj = 0.4)
        if(i == 3)
            mtext("Covariate effect", side = 2, line = 2.5,
                  cex = 0.8, col = "grey25")
        if(i == 4)
            text(x = corners[2]+.45, y = mean(corners[3:4]), "Response",
                 srt = 270, xpd = NA, cex = 1.17, col = "grey25")
        if(i == 5)
            mtext("Year", side = 1, line = 2.0, cex = 0.8, col = "grey25")
        if(i == 6)
            mtext("Covariate", side = 1, line = 2.0, cex = 0.8, col = "grey25")
        if(i %in% c(1, 3, 5))
            lines(x, y, lwd = 0.7, col = "grey70")
        points(x, y, pch = lst[[i]]$pch, col = lst[[i]]$col)
        if(i %in% c(2, 4, 6)) {
            fit1 <- lm(y[1:(N/2)] ~ x[1:(N/2)])
            fit2 <- lm(y[(N/2+1):N] ~ x[(N/2+1):N])
            x.pred1 <- x[1:(N/2)]
            y.pred1 <- predict(fit1)
            x.pred2 <- x[(N/2+1):N]
            y.pred2 <- predict(fit2)
            lines(x.pred1, y.pred1, col = "tomato")
            lines(x.pred2, y.pred2, col = "steelblue")
        }
    }
dev.off()



## Fig: priors ---------------------------------------------

x <- seq(0, 6, length.out = 1000)
p1 <- dt.scaled(x, 3, 0, 0.5)
p2 <- dt.scaled(x, 1, 0, 1)
d1 <- "half-t(3, 0, 0.5)"
d2 <- "half-Cauchy(0, 1)"
df <- data.frame(dist = c(rep(d1, 1000), rep(d2, 1000)),
                 x = c(x, x),
                 y = c(p1, p2))

pdf("./pub/figures/prior_sigma_gamma.pdf", width = 5, height = 4)

    g <- ggplot(df) +
        geom_line(aes(x = x, y = y, linetype = dist, color = dist), size = 1) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_color_manual(values = c("tomato", "steelblue")) +
        scale_linetype_manual(values = c(1, 2)) +
        labs(color = "",
             linetype = "",
             y = "Probability density",
             x = expression(sigma[w])) +
        theme_sleek() +
        theme(legend.justification = c(0, 0),
              legend.position = c(0.6, 0.8),
              legend.key.width =  unit(0.5, "in"))
    print(g)

dev.off()


## Fig: random walk sim ------------------------------------


N <- 46   ## years
K <- 100  ## trials
set.seed(129)
lst <- vector("list", K)
for(i in 1:K) {
    # rt <- abs(rt.scaled(1, 1, 0, 1))  ## draw from prior
    rt <- abs(rt.scaled(1, 3, 0, 0.5))  ## draw from prior
    r1 <- c(rnorm(1, 0, 3), rnorm(N - 1, 0, rt))
    rw <- cumsum(r1)
    lst[[i]] <- data.frame(trial = rep(i, N), year = 1:N, rw = rw)
}
df <- plyr::rbind.fill(lst)
df.avg <- plyr::ddply(df, .(year), summarize,
                      mean = mean(rw),
                      lower = mean(rw) - sd(rw),
                      upper = mean(rw) + sd(rw))


pdf("./pub/figures/rw_sim.pdf", width = 6, height = 4)

    g <- ggplot(df) +
        # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -3, ymax = 3),
        #           fill = "black", alpha = 0.1) +
        geom_line(aes(x = year, y = rw, group = trial),
                  color = alpha("grey50", 0.5)) +
        geom_hline(yintercept = 5, linetype = 2, color = "grey10", size = 1) +
        geom_hline(yintercept = -5, linetype = 2, color = "grey10", size = 1) +
        # geom_line(data = df.avg, aes(x = year, y = mean),
        #           color = "red3", size = 0.75) +
        # geom_line(data = df.avg, aes(x = year, y = lower),
        #           color = "red3", linetype = 2, size = 0.75) +
        # geom_line(data = df.avg, aes(x = year, y = upper),
        #           color = "red3", linetype = 2, size = 0.75) +
        # ylim(-10, 10) +
        scale_x_continuous(expand = c(0, 0)) +
        labs(x = "Year",
             y = expression(gamma)) +
        theme_sleek()
    print(g)

dev.off()



## Fig: Prior predictive -----------------------------------
load_rdata("./output/hbm_prior_pred/")

yrep <- rstan::extract(pr.hbm5.sst, pars = "yrep")[[1]]
y <- rep(stan.dat.sst$y, 4)
ind <- c(1, 2, 4, 6)
tr <- yrep[ind, ]
for(i in 1:nrow(tr)) tr[i, ] <- i
df   <- data.frame(y = y, yrep = as.vector(yrep[ind, ]), trial = as.vector(tr))
labs.pr <- data.frame(trial = 1:4, label = c("(a)", "(b)", "(c)", "(d)"))

pdf("./pub/figures/prior_pred.pdf", width = 6, height = 5)

    g <- ggplot(df) +
        geom_point(aes(x = yrep, y = y), color = "steelblue", alpha = 0.3) +
        geom_abline(aes(intercept = 0, slope = 1), color = "grey30", linetype = 2) +
        facet_wrap( ~ trial) +
        labs(y = "Observed log(R/S)",
             x = "Prior simulated log(R/S)") +
        geom_text(data = labs.pr, aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.4, vjust = 1.6, color = "grey30") +
        theme_sleek() +
        theme(panel.spacing.x = unit(-0.5, "pt"),
              panel.spacing.y = unit(-0.5, "pt"),
              strip.background = element_blank(),
              strip.text.x = element_blank())
    print(g)

dev.off()


## Fig: rolling correlations -------------------------------
dat <- rbind(cor.roll.sst, cor.roll.npgo, cor.roll.pdo)
dat <- ocean_region_lab(dat)
cor.mean <- plyr::ddply(dat, .(var, ocean_region_lab, mid), summarize,
                        cor_mean = mean(cor),
                        n_cor = length(var))
## should have at least stocks with data for a year to average
cor.mean$cor_mean <- ifelse(cor.mean$n_cor <= 3, NA, cor.mean$cor_mean)
xlim  <- range(dat$mid)
xseq  <- seq(1962, 2002, 10)
xlabs <- paste0(xseq - 7, "\n", xseq + 7)  ## window size = 15


g <- ggplot(dat) +
    geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
    geom_line(aes(x = mid, y = cor, group = stock),
              color = "grey75",
              size = 0.4) +
    geom_line(data = cor.mean, aes(x = mid, y = cor_mean),
              color = "grey30",
              size = 1.0) +
    labs(x = "Brood year range",
         y = "Correlation") +
    scale_x_continuous(limits = c(xlim[1], xlim[2]),
                       expand = c(0.01, 0.01),
                       breaks = xseq,
                       labels = xlabs) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.4, vjust = 1.6, color = "grey30") +
    ylim(-1, 1) +
    theme_sleek() +
    theme(panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_grid(var ~ ocean_region_lab)
print(g)


pdf("./pub/figures/roll_cor.pdf", width = 8, height = 5)
print(g)
dev.off()

jpeg("./pub/figures/roll_cor.jpg", width = 8, height = 5, units = "in", res = 500)
print(g)
dev.off()


postscript("./pub/figures/roll_cor.eps", width = 8, height = 5,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(g)
dev.off()



## Fig 2: era ----------------------------------------------

## Dot
g <- ggplot(hbm.gamma.era) +
    aes(era, value, color = era) +
    geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
    stat_summary(fun.y = median,
                 fun.ymin = function(x) quantile(x, probs = 0.025),
                 fun.ymax = function(x) quantile(x, probs = 0.975),
                 geom = "pointrange",
                 size = 0.25) +
    labs(x = "Era",
         y = "Coefficient") +
    scale_color_manual(values = c(rep("steelblue", 3), "grey20"),
                       guide = FALSE) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.4, vjust = 1.6, color = "grey30") +
    theme_sleek() +
    theme(panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_grid(var ~ ocean_region_lab)
print(g)

pdf("./pub/figures/hbm_era_dot.pdf", width = 6, height = 4)
print(g)
dev.off()


## Violin
ss <- plyr::ddply(hbm.gamma.era.stock, .(era, ocean_region_lab, var, stock), summarize,
                  median = median(value))
sig <- ddply(hbm.gamma.era, .(era, ocean_region_lab, var), summarize,
             median = median(value),
             upper95 = quantile(value, 0.975),
             lower95 = quantile(value, 0.025))
sig$sig <- ((sig$lower95 < 0 & sig$upper95 < 0) |
            (sig$lower95 > 0 & sig$upper95 > 0))
hbm.gamma.era2 <- plyr::join(hbm.gamma.era, sig,
                             by = c("era", "var", "ocean_region_lab"))
g <- ggplot(hbm.gamma.era2) +
    geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
    geom_violin(data = hbm.gamma.era2[hbm.gamma.era2$sig, ],
                aes(era, value, color = era, fill = era),
                adjust = 1.5) +
    geom_violin(data = hbm.gamma.era2[!hbm.gamma.era2$sig, ],
                aes(era, value, color = era, fill = NULL),
                adjust = 1.5) +
    geom_point(data = ss,
               aes(era, median, color = era), size = 0.2) +
    labs(x = "Era",
         y = "Coefficient") +
    scale_color_manual(values = c(rep("steelblue", 3), "grey40"),
                       guide = FALSE) +
    scale_fill_manual(values = c(rep("slategray2", 3), "grey80"),
                      guide = FALSE) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.4, vjust = 1.6, color = "grey30") +
    theme_sleek() +
    theme(panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_grid(var ~ ocean_region_lab)
print(g)

pdf("./pub/figures/hbm_era_violin.pdf", width = 7, height = 4.7)
print(g)
dev.off()

jpeg("./pub/figures/hbm_era_violin.jpg", width = 7, height = 4.7, units = "in", res = 500)
print(g)
dev.off()

postscript("./pub/figures/2_hbm_era_violin.eps", width = 7, height = 4.7,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(g)
dev.off()


## Violin S2
pdf("./pub/figures/hbm_era_violin_s2.pdf", width = 6, height = 4)

    ss <- plyr::ddply(hbm.gamma.era.stock.s2, .(era, ocean_region_lab, var, stock), summarize,
                      median = median(value))
    sig <- ddply(hbm.gamma.era.s2, .(era, ocean_region_lab, var), summarize,
                 median = median(value),
                 upper95 = quantile(value, 0.975),
                 lower95 = quantile(value, 0.025))
    sig$sig <- ((sig$lower95 < 0 & sig$upper95 < 0) |
                (sig$lower95 > 0 & sig$upper95 > 0))
    hbm.gamma.era2 <- plyr::join(hbm.gamma.era.s2, sig,
                                 by = c("era", "var", "ocean_region_lab"))
    g <- ggplot(hbm.gamma.era2) +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
        geom_violin(data = hbm.gamma.era2[hbm.gamma.era2$sig, ],
                    aes(era, value, color = era, fill = era),
                    adjust = 1.5) +
        geom_violin(data = hbm.gamma.era2[!hbm.gamma.era2$sig, ],
                    aes(era, value, color = era, fill = NULL),
                    adjust = 1.5) +
        geom_point(data = ss,
                   aes(era, median, color = era), size = 0.2) +
        labs(x = "Era",
             y = "Coefficient") +
        scale_color_manual(values = c(rep("steelblue", 3), "grey40"),
                           guide = FALSE) +
        scale_fill_manual(values = c(rep("slategray2", 3), "grey80"),
                          guide = FALSE) +
        geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.4, vjust = 1.6, color = "grey30") +
        theme_sleek() +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt")) +
        facet_grid(var ~ ocean_region_lab)
    print(g)

dev.off()



## Fig: gamma different ------------------------------------

pdf("./pub/figures/hbm_tv_diff.pdf", width = 8, height = 5)

    g <- ggplot(hbm.gamma.diff) +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2, size = 0.25) +
        geom_line(aes(x = brood_yr, y = gamma, group = stock, color = gamma),
                  size = 0.5) +
        scale_colour_gradient2(low = "steelblue",
                               mid = "grey70",
                               high = "tomato",
                               midpoint = 0,
                               space = "Lab",
                               na.value = "green") +
        labs(x = "Brood year",
             y = expression("Coefficient " (gamma[t])),
             color = "") +
        geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.4, vjust = 1.6, color = "grey30") +
        theme_sleek() +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              legend.justification = c(0, 0),
              legend.position = c(0.81, -0.01),
              legend.key.size = unit(5, "mm"),
              legend.direction = "horizontal",
              legend.text = element_text(size = 7)) +
        facet_grid(var ~ ocean_region_lab)
    print(g)

dev.off()


pdf("./pub/figures/hbm_tv_diff_s1.pdf", width = 8, height = 5)

    g <- ggplot(hbm.gamma.diff.s1) +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2, size = 0.25) +
        geom_line(aes(x = brood_yr, y = gamma, group = stock, color = gamma),
                  size = 0.5) +
        scale_colour_gradient2(low = "steelblue",
                               mid = "grey70",
                               high = "tomato",
                               midpoint = 0,
                               space = "Lab",
                               na.value = "green") +
        labs(x = "Brood year",
             y = expression("Coefficient " (gamma[t])),
             color = "") +
        geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.4, vjust = 1.6, color = "grey30") +
        theme_sleek() +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              legend.justification = c(0, 0),
              legend.position = c(0.81, -0.01),
              legend.key.size = unit(5, "mm"),
              legend.direction = "horizontal",
              legend.text = element_text(size = 7)) +
        facet_grid(var ~ ocean_region_lab)
    print(g)

dev.off()



## Fig 3: gamma different mean -----------------------------
dat <- hbm.gamma.diff
dat <- ocean_region_lab(dat)
cor.mean <- plyr::ddply(dat, .(var, ocean_region_lab, brood_yr), summarize,
                        cor_mean = mean(gamma),
                        n_cor = length(var))
## should have at least stocks with data for a year to average
cor.mean$cor_mean <- ifelse(cor.mean$n_cor <= 3, NA, cor.mean$cor_mean)

pdf("./pub/figures/hbm_tv_diff_mean.pdf", width = 8, height = 5)

    g <- ggplot(dat) +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
        geom_line(aes(x = brood_yr, y = gamma, group = stock),
                  color = "grey75",
                  size = 0.4) +
        geom_line(data = cor.mean, aes(x = brood_yr, y = cor_mean),
                  color = "grey30",
                  size = 1.0) +
        labs(x = "Brood year",
             y = "Coefficient") +
        geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.4, vjust = 1.6, color = "grey30") +
        theme_sleek() +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt")) +
        facet_grid(var ~ ocean_region_lab)
    print(g)

dev.off()


g <- ggplot(dat) +
    geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
    geom_line(aes(x = brood_yr, y = gamma, group = stock),
              color = "grey75",
              size = 0.4) +
    geom_line(data = cor.mean, aes(x = brood_yr, y = cor_mean),
              color = "grey30",
              size = 1.0) +
    labs(x = "Brood year",
         y = "Coefficient") +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.4, vjust = 1.6, color = "grey30") +
    theme_sleek() +
    theme(panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_grid(var ~ ocean_region_lab, scales = "free")
print(g)

pdf("./pub/figures/hbm_tv_diff_mean_ax.pdf", width = 7, height = 4.7)
print(g)
dev.off()

jpeg("./pub/figures/hbm_tv_diff_mean_ax.jpg", width = 7, height = 4.7, units = "in", res = 500)
print(g)
dev.off()

postscript("./pub/figures/3_hbm_tv_diff_mean_ax.eps", width = 7, height = 4.7,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(g)
dev.off()



## Fig 4: gamma diff coef +/- proportions ------------------
df.hbm5.pp <- plyr::ddply(hbm.gamma.diff, .(ocean_region_lab, var, brood_yr), summarize,
						  n_stocks = length(var),
                          n_pos = sum(gamma > 0),
                          n_neg = sum(gamma <= 0),
                          prop_pos = sum(gamma > 0) / length(stock),
                          prop_neg = sum(gamma <= 0) / length(stock))
df.hbm5.pp$neg <- df.hbm5.pp$prop_pos - 1

df.hbm5.pp.h <- df.hbm5.pp[df.hbm5.pp$n_stocks >= 2, ]
df.hbm5.pp.g <- df.hbm5.pp[df.hbm5.pp$n_stocks == 1, ]

df.hbm5.pp[df.hbm5.pp$var == "NPGO", ]
x <- df.hbm5.pp[df.hbm5.pp$var == "SST" &
                df.hbm5.pp$ocean_region_lab == "West Coast", ]
min(x$prop_neg)
which.min(x$prop_neg)

x <- df.hbm5.pp[df.hbm5.pp$var == "SST" &
                df.hbm5.pp$ocean_region_lab == "Gulf of Alaska", ]
min(x$prop_neg)
which.min(x$prop_neg)


g <- ggplot(df.hbm5.pp.h) +
    geom_segment(aes(x = brood_yr, xend = brood_yr, y = 0, yend = prop_pos),
                 color = "tomato", size = 0.6) +
    geom_segment(aes(x = brood_yr, xend = brood_yr, y = 0, yend = neg),
                 color = "steelblue", size = 0.6) +
    geom_segment(data = df.hbm5.pp.g,
                 aes(x = brood_yr, xend = brood_yr, y = 0, yend = prop_pos),
                 color = "grey60", size = 0.6) +
    geom_segment(data = df.hbm5.pp.g,
                 aes(x = brood_yr, xend = brood_yr, y = 0, yend = neg),
                 color = "grey60", size = 0.6) +
    geom_hline(yintercept = 0.0, col = "grey50", linetype = 1) +
    labs(x = "Brood year",
         y = "Proportion of stocks") +
    scale_y_continuous(limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.5),
                       labels = c("1.0", "0.5", "0.0", "0.5", "1.0")) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.4, vjust = 1.6, color = "grey30") +
    theme_sleek() +
    theme(panel.spacing.y = unit(-0.5, "pt"),
          panel.spacing.x = unit(-0.5, "pt")) +
    facet_grid(var ~ ocean_region_lab)
print(g)



pdf("./pub/figures/hbm_tv_diff_prop.pdf", width = 7, height = 4.7)
print(g)
dev.off()

jpeg("./pub/figures/hbm_tv_diff_prop.jpg", width = 7, height = 4.7, units = "in", res = 500)
print(g)
dev.off()

postscript("./pub/figures/4_hbm_tv_diff_prop.eps", width = 7, height = 4.7,
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(g)
dev.off()



## Fig: gamma same -----------------------------------------

pdf("./pub/figures/hbm_tv_same.pdf", width = 8, height = 5)

    g <- ggplot(hbm.gamma.same) +
        geom_hline(yintercept = 0, col = "grey50", linetype = 2) +
        geom_ribbon(aes(x = brood_yr, ymin = lower, ymax = upper),
                    fill = "grey50", alpha = 0.3) +
        geom_line(aes(x = brood_yr, y = gamma),
                  color = "grey20",
                  size = 0.4) +
        labs(x = "Brood year",
             y = "Coefficient") +
        geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.4, vjust = 1.6, color = "grey30") +
        theme_sleek() +
        theme(panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt")) +
        facet_grid(var ~ ocean_region_lab)
    print(g)

dev.off()



## Fig: Stock-specific gamma diff --------------------------

drk <- RColorBrewer::brewer.pal(8, "Dark2")

## SST
pdf("./pub/figures/hbm_tv_diff_ind_sst.pdf", width = 8, height = 7.5)
    g <- ggplot(hbm.gamma.diff[hbm.gamma.diff$var == "SST", ]) +
        geom_hline(yintercept = 0, color = "grey50",
                   linetype = 2, size = 0.25) +
        geom_ribbon(aes(x = brood_yr, ymin = lower, ymax = upper,
                        fill = ocean_region_lab), alpha = 0.2) +
        geom_line(aes(x = brood_yr, y = gamma, color = ocean_region_lab),
                  size = 1) +
        geom_point(data = hbm.gamma.diff[hbm.gamma.diff$var == "SST" &
                                         hbm.gamma.diff$sig, ],
                  aes(x = brood_yr, y = gamma),
                  size = 0.5, color = "grey10") +
        geom_text(data = sock.info, aes(x = 1949, y = 0.9, label = stock),
                  hjust = 0, color = "grey40", size = 2.7) +
        labs(color = "",
             fill = "",
             y = expression("SST coefficient " (gamma[t])),
             x = "Brood year") +
        coord_cartesian(ylim = c(-1.05, 1.05)) +
        theme_sleek() +
        scale_color_manual(values = drk[1:3]) +
        scale_fill_manual(values = drk[1:3]) +
        theme(legend.position = "top",
              strip.background = element_blank(),
              axis.text = element_text(colour = "grey30", size = 8),
              panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              strip.text.x = element_blank()) +
        facet_wrap( ~ stock, as.table = FALSE, ncol = 5)
    print(g)
dev.off()

## NPGO
pdf("./pub/figures/hbm_tv_diff_ind_npgo.pdf", width = 8, height = 7.5)
    g <- ggplot(hbm.gamma.diff[hbm.gamma.diff$var == "NPGO", ]) +
        geom_hline(yintercept = 0, color = "grey50",
                   linetype = 2, size = 0.25) +
        geom_ribbon(aes(x = brood_yr, ymin = lower, ymax = upper,
                        fill = ocean_region_lab), alpha = 0.2) +
        geom_line(aes(x = brood_yr, y = gamma, color = ocean_region_lab),
                  size = 1) +
        geom_point(data = hbm.gamma.diff[hbm.gamma.diff$var == "NPGO" &
                                         hbm.gamma.diff$sig, ],
                  aes(x = brood_yr, y = gamma),
                  size = 0.5, color = "grey10") +
        geom_text(data = sock.info, aes(x = 1949, y = 0.9, label = stock),
                  hjust = 0, color = "grey40", size = 2.7) +
        labs(color = "",
             fill = "",
             y = expression("NPGO coefficient " (gamma[t])),
             x = "Brood year") +
        coord_cartesian(ylim = c(-1.05, 1.05)) +
        theme_sleek() +
        scale_color_manual(values = drk[1:3]) +
        scale_fill_manual(values = drk[1:3]) +
        theme(legend.position = "top",
              strip.background = element_blank(),
              axis.text = element_text(colour = "grey30", size = 8),
              panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              strip.text.x = element_blank()) +
        facet_wrap( ~ stock, as.table = FALSE, ncol = 5)
    print(g)
dev.off()


## PDO
pdf("./pub/figures/hbm_tv_diff_ind_pdo.pdf", width = 8, height = 7.5)
    g <- ggplot(hbm.gamma.diff[hbm.gamma.diff$var == "PDO", ]) +
        geom_hline(yintercept = 0, color = "grey50",
                   linetype = 2, size = 0.25) +
        geom_ribbon(aes(x = brood_yr, ymin = lower, ymax = upper,
                        fill = ocean_region_lab), alpha = 0.2) +
        geom_line(aes(x = brood_yr, y = gamma, color = ocean_region_lab),
                  size = 1) +
        geom_point(data = hbm.gamma.diff[hbm.gamma.diff$var == "PDO" &
                                         hbm.gamma.diff$sig, ],
                  aes(x = brood_yr, y = gamma),
                  size = 0.5, color = "grey10") +
        geom_text(data = sock.info, aes(x = 1949, y = 0.9, label = stock),
                  hjust = 0, color = "grey40", size = 2.7) +
        labs(color = "",
             fill = "",
             y = expression("PDO coefficient " (gamma[t])),
             x = "Brood year") +
        coord_cartesian(ylim = c(-1.05, 1.05)) +
        theme_sleek() +
        scale_color_manual(values = drk[1:3]) +
        scale_fill_manual(values = drk[1:3]) +
        theme(legend.position = "top",
              strip.background = element_blank(),
              axis.text = element_text(colour = "grey30", size = 8),
              panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              strip.text.x = element_blank()) +
        facet_wrap( ~ stock, as.table = FALSE, ncol = 5)
    print(g)
dev.off()



## Fig: Stock-specific rolling cors ------------------------
drk <- RColorBrewer::brewer.pal(8, "Dark2")
rkr.sst  <- ocean_region_lab(rk.roll.sst)
rkr.npgo <- ocean_region_lab(rk.roll.npgo)
rkr.pdo  <- ocean_region_lab(rk.roll.pdo)
xlim  <- range(rkr.sst$mid)
xseq  <- seq(1962, 2002, 10)
xlabs <- paste0(xseq - 7, "\n", xseq + 7)  ## window size = 15

pdf("./pub/figures/roll_cor_sst.pdf", width = 8, height = 7.5)
    g <- ggplot(rkr.sst) +
        geom_hline(yintercept = 0, color = "grey50",
                   linetype = 2, size = 0.25) +
        geom_ribbon(aes(x = mid, ymin = ci_lower, ymax = ci_upper,
                        fill = ocean_region_lab), alpha = 0.2) +
        geom_line(aes(x = mid, y = estimate, color = ocean_region_lab),
                  size = 1) +
        geom_point(data = rkr.sst[rkr.sst$sig, ],
                   aes(x = mid, y = estimate),
                  size = 0.5, color = "grey10") +
        geom_text(data = sock.info, aes(x = 1958, y = 0.9, label = stock),
                  hjust = 0, color = "grey40", size = 2.7) +
        labs(color = "",
             fill = "",
             y = "Correlation",
             x = "Brood year range") +
        coord_cartesian(ylim = c(-1.05, 1.05)) +
        theme_sleek() +
        scale_color_manual(values = drk[1:3]) +
        scale_fill_manual(values = drk[1:3]) +
        scale_x_continuous(limits = c(xlim[1], xlim[2]),
                           expand = c(0.01, 0.01),
                           breaks = xseq,
                           labels = xlabs) +
        theme(legend.position = "top",
              strip.background = element_blank(),
              axis.text = element_text(colour = "grey30", size = 8),
              panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              strip.text.x = element_blank()) +
        facet_wrap( ~ stock, as.table = FALSE, ncol = 5)
    print(g)
dev.off()


pdf("./pub/figures/roll_cor_npgo.pdf", width = 8, height = 7.5)
    g <- ggplot(rkr.npgo) +
        geom_hline(yintercept = 0, color = "grey50",
                   linetype = 2, size = 0.25) +
        geom_ribbon(aes(x = mid, ymin = ci_lower, ymax = ci_upper,
                        fill = ocean_region_lab), alpha = 0.2) +
        geom_line(aes(x = mid, y = estimate, color = ocean_region_lab),
                  size = 1) +
        geom_point(data = rkr.npgo[rkr.npgo$sig, ],
                   aes(x = mid, y = estimate),
                  size = 0.5, color = "grey10") +
        geom_text(data = sock.info, aes(x = 1958, y = 0.9, label = stock),
                  hjust = 0, color = "grey40", size = 2.7) +
        labs(color = "",
             fill = "",
             y = "Correlation",
             x = "Brood year range") +
        coord_cartesian(ylim = c(-1.05, 1.05)) +
        theme_sleek() +
        scale_color_manual(values = drk[1:3]) +
        scale_fill_manual(values = drk[1:3]) +
        scale_x_continuous(limits = c(xlim[1], xlim[2]),
                           expand = c(0.01, 0.01),
                           breaks = xseq,
                           labels = xlabs) +
        theme(legend.position = "top",
              strip.background = element_blank(),
              axis.text = element_text(colour = "grey30", size = 8),
              panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              strip.text.x = element_blank()) +
        facet_wrap( ~ stock, as.table = FALSE, ncol = 5)
    print(g)
dev.off()

pdf("./pub/figures/roll_cor_pdo.pdf", width = 8, height = 7.5)
    g <- ggplot(rkr.pdo) +
        geom_hline(yintercept = 0, color = "grey50",
                   linetype = 2, size = 0.25) +
        geom_ribbon(aes(x = mid, ymin = ci_lower, ymax = ci_upper,
                        fill = ocean_region_lab), alpha = 0.2) +
        geom_line(aes(x = mid, y = estimate, color = ocean_region_lab),
                  size = 1) +
        geom_point(data = rkr.pdo[rkr.pdo$sig, ],
                   aes(x = mid, y = estimate),
                  size = 0.5, color = "grey10") +
        geom_text(data = sock.info, aes(x = 1958, y = 0.9, label = stock),
                  hjust = 0, color = "grey40", size = 2.7) +
        labs(color = "",
             fill = "",
             y = "Correlation",
             x = "Brood year range") +
        coord_cartesian(ylim = c(-1.05, 1.05)) +
        theme_sleek() +
        scale_color_manual(values = drk[1:3]) +
        scale_fill_manual(values = drk[1:3]) +
        scale_x_continuous(limits = c(xlim[1], xlim[2]),
                           expand = c(0.01, 0.01),
                           breaks = xseq,
                           labels = xlabs) +
        theme(legend.position = "top",
              strip.background = element_blank(),
              axis.text = element_text(colour = "grey30", size = 8),
              panel.spacing.y = unit(-0.5, "pt"),
              panel.spacing.x = unit(-0.5, "pt"),
              strip.text.x = element_blank()) +
        facet_wrap( ~ stock, as.table = FALSE, ncol = 5)
    print(g)
dev.off()
