
## ----setup, include=FALSE, cache=FALSE-----------------------------------
# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',
               dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)


## ----Rcode---------------------------------------------------------------
source('/home/ricky/Rlim/ChromatinConformation/MotifCalls/MotifCalls.R')
work_dir = '/home/ricky/Rlim/ChromatinConformation/MotifCalls/P53/'


## ----CentDistP53---------------------------------------------------------
NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_direct <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_direct/result.txt'))
NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_indirect <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_indirect/result.txt'))


## ------------------------------------------------------------------------
sessionInfo()


## ----knitIt,include=FALSE, cache=TRUE, results='hide', message=FALSE, warning=FALSE----
library(knitr)
purl("motifAnalysis.Rnw" ) # compile to tex
purl("motifAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("motifAnalysis.Rnw")


