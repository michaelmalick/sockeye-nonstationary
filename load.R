## Load functions and libraries for analysis


## Load external packages
library(lattice)
library(nlme)
library(plyr)
library(KFAS)
library(MARSS)
library(RColorBrewer)
library(malick)    # see: https://github.com/michaelmalick/r-malick
library(chroma)    # see: https://github.com/michaelmalick/r-chroma
library(codatools) # see: https://github.com/michaelmalick/r-codatools
library(ersst)     # see: https://github.com/michaelmalick/r-ersst
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(bayesplot))
suppressPackageStartupMessages(library(maps))
suppressPackageStartupMessages(library(rworldmap))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(loo))


## Load datasets
data(countriesLow)


## Source scripts
source("./functions.R")
source("./bin/tsim.R")


## Setup directories
if(!dir.exists("./figures/"))
    dir.create("./figures/")
if(!dir.exists("./output/"))
    dir.create("./output/")
if(!dir.exists("./data/"))
    dir.create("./data/")
if(!dir.exists("./logs/"))
    dir.create("./logs/")


## Load output data
for(i in list.files(path = "./output/", pattern = "*.RData$")) {
    load(paste("./output/", i, sep = ""))
}


## Set bayesplot theme
bayesplot::bayesplot_theme_set(new = theme_sleek())
