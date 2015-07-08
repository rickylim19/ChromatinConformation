## ----setup, include=FALSE, cache=TRUE------------------------------------
# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',
               dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)

## ------------------------------------------------------------------------
sessionInfo()

## ----knitIt, include=FALSE, cache=TRUE, message=FALSE, warning=FALSE-----
library(knitr)
purl("chromatinConformation.Rnw" ) # compile to tex
purl("chromatinConformation.Rnw", documentation = 0) # extract R code only
knit2pdf("chromatinConformation.Rnw")

