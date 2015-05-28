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

## ----Rcodes--------------------------------------------------------------
work_dir = '/home/ricky/Rlim/ChromatinConformation/PeakCalls/CebpE/'
source('/home/ricky/Rlim/ChromatinConformation/PeakCalls/PeakCalls.R')

## ----dataPrep, cache=TRUE------------------------------------------------
KoefflerLab_BM_ChIPseq_CebpE_3000 <- read.table(paste0(work_dir, 
                           'Input/3000bin-KoefflerLab_BM_ChIPseq_CebpE_Input.bed'),
                                                header=T)
KoefflerLab_BM_ChIPseq_CebpE_300 <- read.table(paste0(work_dir, 
                           'Input/300bin-KoefflerLab_BM_ChIPseq_CebpE_Input.bed'),
                                               header=T)

head(KoefflerLab_BM_ChIPseq_CebpE_3000)
KoefflerLab_BM_ChIPseq_CebpE_300_hmm <- KoefflerLab_BM_ChIPseq_CebpE_300[, c(1,4,5)]
KoefflerLab_BM_ChIPseq_CebpE_3000_hmm <- KoefflerLab_BM_ChIPseq_CebpE_3000[, c(1,4,5)]

## ----fitjaHMM, cache=TRUE------------------------------------------------
set.seed(12345)
fit_300 <- jahmm(KoefflerLab_BM_ChIPseq_CebpE_300_hmm)
fit_3000 <- jahmm(KoefflerLab_BM_ChIPseq_CebpE_3000_hmm)

plotFit(KoefflerLab_BM_ChIPseq_CebpE_300_hmm, fit_300, 1:100000, 
        paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_300_Rjahmm.pdf'))

plotFit(KoefflerLab_BM_ChIPseq_CebpE_3000_hmm, fit_3000, 1:10000, 
        paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_3000_Rjahmm.pdf'))

no_fit_300 <- table(fit_300$path)
no_fit_300 <- as.data.frame(no_fit_300)
colnames(no_fit_300) <- c('state', 'freq')
no_fit_3000 <- table(fit_3000$path)
no_fit_3000 <- as.data.frame(no_fit_3000)
colnames(no_fit_3000) <- c('state', 'freq')


## ----peakTable, echo=FALSE, results='asis'-------------------------------
print(xtable(no_fit_300, caption='Frequency of States in 300bp Peak Resolution'))
print(xtable(no_fit_3000, caption='Frequency of States in 3000bp Peak Resolution'))

## ----pdfTable, echo=FALSE------------------------------------------------
createPDFTable(no_fit_300, 'Stat Freq in 300bp', paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_300_Rjahmm_Freq.pdf'))
createPDFTable(no_fit_3000, 'Stat Freq in 3000bp', paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_3000_Rjahmm_Freq.pdf'))

## ----peakCallsJahmm, cache=TRUE------------------------------------------
KoefflerLab_BM_ChIPseq_CebpE_300$path <- fit_300$path
KoefflerLab_BM_ChIPseq_CebpE_3000$path <- fit_3000$path
head(KoefflerLab_BM_ChIPseq_CebpE_300)

KoefflerLab_BM_ChIPseq_CebpE_300_targets <- with(KoefflerLab_BM_ChIPseq_CebpE_300, 
                                subset(KoefflerLab_BM_ChIPseq_CebpE_300, path == 2))
KoefflerLab_BM_ChIPseq_CebpE_3000_targets<- with(KoefflerLab_BM_ChIPseq_CebpE_3000, 
                                subset(KoefflerLab_BM_ChIPseq_CebpE_3000, path == 2))

targets_300 <- nrow(KoefflerLab_BM_ChIPseq_CebpE_300_targets)
targets_3000 <- nrow(KoefflerLab_BM_ChIPseq_CebpE_3000_targets)

write.table(KoefflerLab_BM_ChIPseq_CebpE_300_targets, 
            paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_300_targets.bed'),  
            row.names=F, quote=F, sep='\t')
write.table(KoefflerLab_BM_ChIPseq_CebpE_3000_targets, 
            paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_3000_targets.bed'),  
            row.names=F, quote=F, sep='\t')

## ----gRangesPrep, cache=TRUE---------------------------------------------

macs_targets <- read.table(paste0(work_dir,
                    'Output/300bin-KoefflerLab_BM_ChIPseq_CebpE_mm10_summits_filtered.bed'))
colnames(macs_targets) <- c('chr', 'start', 'end', 'peakScore')
# bam_scores for number of reads for target and input
bam_scores <- read.table(paste0(work_dir,
                    'Input/300bin-KoefflerLab_BM_ChIPseq_CebpE_Input.bed'), header=T)

# substitute peakScore with the bam score
macs_targets <- merge(macs_targets, bam_scores, by=c('chr', 'start', 'end') )
macs_targets <- macs_targets[, c(1,2,3,5,6)]


jahmm_targets <- KoefflerLab_BM_ChIPseq_CebpE_300_targets[,c(1,2,3,4,5)]
rownames(jahmm_targets) <- 1:nrow(jahmm_targets)
macs_targets <- data.table(macs_targets, key=c('chr', 'start', 'end'))
jahmm_targets <- data.table(jahmm_targets, key=c('chr', 'start', 'end'))

macs_targetsOnly <- macs_targets[!jahmm_targets]
jahmm_targetsOnly <- jahmm_targets[!macs_targets]
macs_jahmm <- merge(macs_targets, jahmm_targets,
                    by=c('chr', 'start', 'end'))

head(macs_jahmm)
macs_jahmm_targets <- macs_jahmm[, 1:5, with=F]
head(macs_jahmm_targets)
setnames(macs_jahmm_targets, c('chr', 'start', 'end', 'Input.x', 'CebpE.x'),
                             c('chr', 'start', 'end', 'Input', 'CebpE'))

## ----plotMacsVsJahmm, cache=TRUE-----------------------------------------
plotTargetInput(macs_targetsOnly, 1:500, ylim= c(0,100), 
                'Output/KoefflerLab_BM_ChIPseq_CebpE_macs_targets.pdf',
                title_f = 'macs Targets Only')     

plotTargetInput(jahmm_targetsOnly, 1:500, ylim = c(0,100),  
                'Output/KoefflerLab_BM_ChIPseq_CebpE_jahmm_targets.pdf',     
                title_f = 'jahmm Targets Only')     

plotTargetInput(macs_jahmm_targets, 1:500, ylim=c(0,100), 
                'Output/KoefflerLab_BM_ChIPseq_CebpE_macs_jahmm_targets.pdf',     
                title_f = 'macs and jahmm Targets')     

# plot venn diagram
pdf(paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_macs_jahmm_venn.pdf'),
    pointsize=21)
v <- draw.pairwise.venn(nrow(macs_targets), nrow(jahmm_targets), nrow(macs_jahmm),
                        c('macs', 'jahmm'), 
                        fill=c('red', 'blue'), lty='blank', 
                        cex = 1, lwd = 1)
dev.off()

## ----readsAlignment, cache=TRUE, echo = F--------------------------------
#library(ggbio)
#library(GenomicAlignments)
#set.seed(12345)

#KoefflerLab_BM_ChIPseq_CebpE_300_jahmm_granges <- 
#    makeGRangesFromDataFrame(KoefflerLab_BM_ChIPseq_CebpE_300_jahmm,
#                             keep.extra.columns = TRUE,
#                             ignore.strand = TRUE)
#
#KoefflerLab_BM_ChIPseq_CebpE_300_macs_granges <- 
#    makeGRangesFromDataFrame(KoefflerLab_BM_ChIPseq_CebpE_300_macs,
#                             keep.extra.columns = TRUE,
#                             ignore.strand = TRUE)
#
#KoefflerLab_BM_ChIPseq_CebpE_summits_macs_granges <- 
#    makeGRangesFromDataFrame(KoefflerLab_BM_ChIPseq_CebpE_summits_macs,
#                             keep.extra.columns = TRUE,
#                             ignore.strand = TRUE)

plotReadsOnTargets <- function(targets, n){
 
    total_targets <- length(targets)
    #select randomnly n targets to be plotted with the aligned reads
    ix <- sample(1:total_targets, n, replace = T)
    selected_targets <- targets[ix]
    for (i in 1:n){
        which <- selected_targets[i]
        param <- ScanBamParam(which = which)
        alignments <- readGAlignments(bam_file, index = bam_index, param = param)

        alignmentPlot <- autoplot(alignments)
        coveragePlot <- ggplot(as(alignments, 'GRanges')) + 
                            stat_coverage(color='gray40', fill = 'skyblue')

        pdf(paste0(work_dir, 'Output/KoefflerLab_BM_ChIPseq_CebpE_bam_coverage_', i ,'.pdf'))
        print(tracks(alignmentPlot, coveragePlot))
        dev.off()
    }
}

#plotReadsOnTargets(KoefflerLab_BM_ChIPseq_CebpE_300_jahmm_granges, 3)



## ------------------------------------------------------------------------
sessionInfo()

## ----knitIt, message=FALSE, warning=FALSE, cache=TRUE, include=FALSE-----
library(knitr)
purl("peakAnalysis.Rnw" ) # compile to tex
purl("peakAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("peakAnalysis.Rnw")

