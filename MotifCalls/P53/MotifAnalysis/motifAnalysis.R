
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


## ----boxplotProfile, cache=TRUE, echo=FALSE, results=FALSE---------------

NCMLS_SAOS2_P53_group1_compSorted4 <- read.csv(sep='\t',paste0(work_dir, 
    'Input/NCMLS_SAOS2_ChIPseq_P53_GMM_ModelAssignment_log_group1_compSorted4_ProfileMatrix.txt'),
    header=TRUE, row.names=1, check.names=F)
NCMLS_SAOS2_P53_group2_compSorted4 <- read.csv(sep='\t',paste0(work_dir, 
    'Input/NCMLS_SAOS2_ChIPseq_P53_GMM_ModelAssignment_log_group2_compSorted4_ProfileMatrix.txt'),
    header=TRUE, row.names=1, check.names=F)
NCMLS_SAOS2_P53_group3_compSorted4 <- read.csv(sep='\t',paste0(work_dir, 
    'Input/NCMLS_SAOS2_ChIPseq_P53_GMM_ModelAssignment_log_group3_compSorted4_ProfileMatrix.txt'),
    header=TRUE, row.names=1, check.names=F)
NCMLS_SAOS2_P53_group4_compSorted4 <- read.csv(sep='\t',paste0(work_dir, 
    'Input/NCMLS_SAOS2_ChIPseq_P53_GMM_ModelAssignment_log_group4_compSorted4_ProfileMatrix.txt'),
    header=TRUE, row.names=1, check.names=F)

pdf('figs/NCMLS_SAOS2_P53_BoxplotProfile.pdf', 
    useDingbats=FALSE)
boxplot(NCMLS_SAOS2_P53_group1_compSorted4, cex=0.1, col='red', bty='n')
boxplot(NCMLS_SAOS2_P53_group2_compSorted4, cex=0.1, col='green', add=T, bty='n')
boxplot(NCMLS_SAOS2_P53_group3_compSorted4, cex=0.1, col='yellow', add=T, bty='n')
boxplot(NCMLS_SAOS2_P53_group4_compSorted4, cex=0.1, col='blue', add=T, bty='n')
dev.off()


## ----CentDistP53Group, echo=FALSE, cache=TRUE----------------------------
NCMLS_SAOS2_ChIPseq_P53_group1_compSorted4_Centdist <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_ModelAssignment_GMM_group1_compSorted4_0.5kb/result.txt'))
NCMLS_SAOS2_ChIPseq_P53_group2_compSorted4_Centdist <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_ModelAssignment_GMM_group2_compSorted4_0.5kb/result.txt'))
NCMLS_SAOS2_ChIPseq_P53_group3_compSorted4_Centdist <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_ModelAssignment_GMM_group3_compSorted4_0.5kb/result.txt'))
NCMLS_SAOS2_ChIPseq_P53_group4_compSorted4_Centdist <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_ModelAssignment_GMM_group4_compSorted4_0.5kb/result.txt'))


## ----echo=FALSE, results=FALSE-------------------------------------------
CTCF_TP53_TP63_noDuplicate <- read.table(paste0(work_dir, 
                                    'Output/ChipMeme_NCMLS_SAOS2/',
                                    'Motifs/CTCF_TP53_TP63_Motif_Groups_noDuplicates.bed'))

colnames(CTCF_TP53_TP63_noDuplicate) <- c('Chr', 'Start', 'End', 'Motif', 'Score', 'Comp')
CTCF_noDuplicate <- subset(CTCF_TP53_TP63_noDuplicate, Motif=='CTCF')
TP53_noDuplicate <- subset(CTCF_TP53_TP63_noDuplicate, Motif=='TP53')
TP63_noDuplicate <- subset(CTCF_TP53_TP63_noDuplicate, Motif=='TP63')

# Total Component Frequency
Freq_groups <- read.table(paste0(work_dir,
                    'Output/ChipMeme_NCMLS_SAOS2/',
                    'Motifs/NCMLS_SAOS2_ChIPseq_P53_GMM_ModelAssignment_log_compSorted4_0.5kb.bed'))
total_freq <- as.data.frame(table(Freq_groups$V4))$Freq

CTCF_noDuplicate <- as.data.frame(table(CTCF_noDuplicate$Comp))
colnames(CTCF_noDuplicate) <- c('Comp', 'Freq')
CTCF_noDuplicate$Total <- total_freq
CTCF_noDuplicate$Percentage <- CTCF_noDuplicate$Freq/CTCF_noDuplicate$Total*100.0

TP53_noDuplicate <- as.data.frame(table(TP53_noDuplicate$Comp))
colnames(TP53_noDuplicate) <- c('Comp', 'Freq')
TP53_noDuplicate$Total <- total_freq
TP53_noDuplicate$Percentage <- TP53_noDuplicate$Freq/TP53_noDuplicate$Total*100.0

TP63_noDuplicate <- as.data.frame(table(TP63_noDuplicate$Comp))
colnames(TP63_noDuplicate) <- c('Comp', 'Freq')
TP63_noDuplicate$Total <- total_freq
TP63_noDuplicate$Percentage <- TP63_noDuplicate$Freq/TP63_noDuplicate$Total*100.0





## ----results='asis', echo=FALSE------------------------------------------
print(xtable(CTCF_noDuplicate, 
      caption='CTCF Discovered by Centrimo in Groups of 4-Component Log-Transformed ChIPseq'), include.rownames=FALSE)
print(xtable(TP53_noDuplicate, 
      caption='TP53 Discovered by Centrimo in Groups in 4-Component Log-Transformed ChIPseq'), include.rownames=FALSE)
print(xtable(TP63_noDuplicate, 
      caption='TP63 Discovered by Centrimo in Groups of 4-Component Log-Transformed ChIPseq:no duplicate'), include.rownames=FALSE)


## ----CentDistP53, echo=FALSE---------------------------------------------
NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_direct <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_direct/result.txt'))
NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_indirect <- 
    read.delim(paste0(work_dir,'Output/Centdist/',
    'NCMLS_SAOS2_ChIPseq_P53_Bicluster_compSorted4_dist50kb_indirect/result.txt'))


## ----TP53TP63CTCF, cache=TRUE, echo=FALSE--------------------------------
memeMotif_dir=paste0('/home/ricky/Rlim/ChromatinConformation/MotifCalls/P53/',
                     'Output/ChipMeme_NCMLS_SAOS2/Motifs/')

meme_msite <- read.table(paste0(memeMotif_dir,'TP53_TP63_CTCF.msite' ))

TP53_motif <- meme_msite[meme_msite$V2 == 'TP53', ]
TP63_motif <- meme_msite[meme_msite$V2 == 'TP63', ]
CTCF_motif <- meme_msite[meme_msite$V2 == 'CTCF', ]

peak_0.5kb <- read.table(paste0(work_dir, 
                                'Input/NCMLS_SAOS2_ChIPseq_P53_Peaks_0.5kb.bed'))
colnames(peak_0.5kb) <- c('chr', 'start', 'end', 'score')

TP53_score <- getPeakScore(TP53_motif, peak_0.5kb)
TP63_score <- getPeakScore(TP63_motif, peak_0.5kb)
CTCF_score <- getPeakScore(CTCF_motif, peak_0.5kb)
allmotif_score <- getPeakScore(meme_msite, peak_0.5kb)

colnames(TP53_score) <- c('chr', 'start', 'end', 'id', 'score')
colnames(TP63_score) <- c('chr', 'start', 'end', 'id', 'score')
colnames(CTCF_score) <- c('chr', 'start', 'end', 'id', 'score')
colnames(allmotif_score) <- c('chr', 'start', 'end', 'id', 'score')




## ----plotTP53TP63CTCF, cache=TRUE, results='hide', echo=FALSE------------

p <- ggplot(TP53_score, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of TP53 Motif')
ggsave(file='figs/NCMLS_SAOS2_ChIPseq_P53_Score_TP53Motif.pdf')

p <- ggplot(TP63_score, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of TP63 Motif')
ggsave(file='figs/NCMLS_SAOS2_ChIPseq_P63_Score_TP53Motif.pdf')

p <- ggplot(CTCF_score, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CTCF Motif')
ggsave(file='figs/NCMLS_SAOS2_ChIPseq_CTCF_Score_TP53Motif.pdf')

p <- ggplot(allmotif_score, aes(score, colour=id, group=id)) +
        geom_density()  
ggsave(file='figs/NCMLS_SAOS2_ChIPseq_Allmotif_Score_TP53Motif.pdf')



## ------------------------------------------------------------------------
sessionInfo()


## ----knitIt,cache=TRUE, include=FALSE, results='hide', message=FALSE, warning=FALSE----
library(knitr)
purl("motifAnalysis.Rnw" ) # compile to tex
purl("motifAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("motifAnalysis.Rnw")


