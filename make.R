## Reproduce project
##
## Sourcing this file will reproduce the entire project.


## Clean-up workspace
rm(list=ls())
graphics.off()
source("./functions.R")


## Specify project script order
m.scripts <- c("load.R",
               "02_data_sock.R",
               "03_data_enviro.R",
               "04_data_sock_enviro.R",
               "05_data_explore.R",
               "06_roll_cor.R",
               "07_roll_rk.R",
               "08_dlm.R",
               "09_hbm_fit.R",
               "10_hbm_inf.R",
               "present.R",
               "pub.R")


## Specify dynamic dirs
m.remove <- c("output",
              "figures")


## Run make
m.time <- make(scripts = m.scripts, remove = m.remove)
# m.time <- make(scripts = m.scripts, remove = NULL)


## Write log file
make_log(dir = "./logs",
         file = "make.log",
         session.info = TRUE,
         run.time = paste(round(m.time, 2), names(m.time), collapse = "   "))
