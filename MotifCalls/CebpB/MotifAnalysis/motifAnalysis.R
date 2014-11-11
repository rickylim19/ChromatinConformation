
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
work_dir = '/home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpB/'


## ----CentDistP53---------------------------------------------------------
work_dir_centdist = '/home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpB/'
HudsonAlpha_GM12878_CebpB_Bicluster_compSorted4_dist50kb_direct <- 
    read.delim(paste0(work_dir_centdist,'Output/CentDist/',
    'HudsonAlpha_GM12878_CebpB_Bicluster_compSorted4_dist50kb_direct/',
    'result.txt'))
HudsonAlpha_GM12878_CebpB_Bicluster_compSorted4_dist50kb_indirect <- 
    read.delim(paste0(work_dir_centdist,'Output/CentDist/',
    'HudsonAlpha_GM12878_CebpB_Bicluster_compSorted4_dist50kb_indirect/',
    'result.txt'))

#HCT116
HudsonAlpha_HCT116_CebpB_Bicluster_compSorted4_dist50kb_direct <- 
    read.delim(paste0(work_dir_centdist,'Output/CentDist/',
    'HudsonAlpha_HCT116_CebpB_Bicluster_compSorted4_dist50kb_direct/',
    'result.txt'))

HudsonAlpha_HCT116_CebpB_Bicluster_compSorted4_dist50kb_indirect <- 
    read.delim(paste0(work_dir_centdist,'Output/CentDist/',
    'HudsonAlpha_HCT116_CebpB_Bicluster_compSorted4_dist50kb_indirect/',
    'result.txt'))

#K562
HudsonAlpha_K562_CebpB_Bicluster_compSorted4_dist50kb_direct <- 
    read.delim(paste0(work_dir_centdist,'Output/CentDist/',
    'HudsonAlpha_K562_CebpB_Bicluster_compSorted4_dist50kb_direct/',
    'result.txt'))

HudsonAlpha_K562_CebpB_Bicluster_compSorted4_dist50kb_indirect <- 
    read.delim(paste0(work_dir_centdist,'Output/CentDist/',
    'HudsonAlpha_K562_CebpB_Bicluster_compSorted4_dist50kb_indirect/',
    'result.txt'))


## ----MemeGM12878---------------------------------------------------------
meme_chipScore_GM12878 <- getCebpAllMotif(peakScore_f=paste0(work_dir,
                                          'Input/HudsonAlpha_GM12878_CebpB_Peaks_0.5kb.bed'), 
                                          allmotif_f=paste0(work_dir,
                                          'Output/ChipMeme_HudsonAlpha_GM12878/Motif/allmotif.msite'),
                                          cebp_motif= 'motif_3')
head(meme_chipScore_GM12878$allmotif)
head(meme_chipScore_GM12878$cebpmotif)
p <- ggplot(meme_chipScore_GM12878$cebpmotif, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif: GM12878-MemeChip')
ggsave(file='figs/HudsonAlpha_GM12878_CebpB_Score_CebpMotif.pdf')

p <- ggplot(meme_chipScore_GM12878$allmotif, aes(score, colour=id, group=id)) +
        geom_density()  
ggsave(file='figs/HudsonAlpha_GM12878_CebpB_Score_AllMotif.pdf')



## ----MemeHCT116----------------------------------------------------------
meme_chipScore_HCT116 <- getCebpAllMotif(peakScore_f=paste0(work_dir,
                                          'Input/HudsonAlpha_HCT116_CebpB_Peaks_0.5kb.bed'), 
                                          allmotif_f=paste0(work_dir,
                                          'Output/ChipMeme_HudsonAlpha_HCT116/Motif/allmotif.msite'),
                                          cebp_motif= 'motif_1')
p <- ggplot(meme_chipScore_HCT116$cebpmotif, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif: HCT116-MemeChip')
ggsave(file='figs/HudsonAlpha_HCT116_CebpB_Score_CebpMotif.pdf')

p <- ggplot(meme_chipScore_HCT116$allmotif, aes(score, colour=id, group=id)) +
        geom_density()  
ggsave(file='figs/HudsonAlpha_HCT116_CebpB_Score_AllMotif.pdf')



## ----MemeK562------------------------------------------------------------
meme_chipScore_K562 <- getCebpAllMotif(peakScore_f=paste0(work_dir,
                                          'Input/HudsonAlpha_K562_CebpB_Peaks_0.5kb.bed'), 
                                          allmotif_f=paste0(work_dir,
                                          'Output/ChipMeme_HudsonAlpha_K562/Motif/allmotif.msite'),
                                          cebp_motif= 'motif_1')
p <- ggplot(meme_chipScore_K562$cebpmotif, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif: K562-MemeChip')
ggsave(file='figs/HudsonAlpha_K562_CebpB_Score_CebpMotif.pdf')

p <- ggplot(meme_chipScore_K562$allmotif, aes(score, colour=id, group=id)) +
        geom_density()  
ggsave(file='figs/HudsonAlpha_K562_CebpB_Score_AllMotif.pdf')



## ----DataPreProcessingMotif, cache=TRUE----------------------------------


rsat_GM12878 <- getCebpNon_rsat(input_f=paste0(work_dir, 
    'Output/Rsat/HudsonAlpha_GM12878_CebpB_Peaks_1kb/results/sites/peak_motifs_sites.tab'),
                               peak_f=paste0(work_dir, 
                                          'Input/HudsonAlpha_GM12878_CebpB_Peaks_1kb.bed'), 
                               jaspar_f=paste0(work_dir, 
    'Output/Rsat/HudsonAlpha_GM12878_CebpB_Peaks_1kb/results',
    '/discovered_vs_db/peak-motifs_db_jaspar.tab'))

rsat_HCT116 <- getCebpNon_rsat(input_f=paste0(work_dir, 
    'Output/Rsat/HudsonAlpha_HCT116_CebpB_Peaks_1kb/results/sites/peak_motifs_sites.tab'),
                               peak_f=paste0(work_dir, 
                                          'Input/HudsonAlpha_HCT116_CebpB_Peaks_1kb.bed'), 
                               jaspar_f=paste0(work_dir, 
    'Output/Rsat/HudsonAlpha_HCT116_CebpB_Peaks_1kb/results',
    '/discovered_vs_db/peak-motifs_db_jaspar.tab'))

rsat_K562 <- getCebpNon_rsat(input_f=paste0(work_dir, 
    'Output/Rsat/HudsonAlpha_K562_CebpB_Peaks_1kb/results/sites/peak_motifs_sites.tab'),
                               peak_f=paste0(work_dir, 
                                          'Input/HudsonAlpha_K562_CebpB_Peaks_1kb.bed'), 
                               jaspar_f=paste0(work_dir, 
    'Output/Rsat/HudsonAlpha_K562_CebpB_Peaks_1kb/results',
    '/discovered_vs_db/peak-motifs_db_jaspar.tab'))



## ----plotRsat------------------------------------------------------------

## GM12878
p <- ggplot(rsat_GM12878$nonCEBP, 
            aes(score, colour=motifName, group=motifName )) +
        geom_density()
ggsave(file='figs/HudsonAlpha_GM12878_CebpB_Score_nonCEBP_rsat.pdf')

# plot motifs in chunks, each for 10 motifs
plotNMotifScore(rsat_GM12878$nonCEBP, 
                output_f ='figs/HudsonAlpha_GM12878_CebpB_Score_nonCEBP_rsat_chunk', 
                N=5)

p <- ggplot(rsat_GM12878$CEBP, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif:GM12878')
ggsave(file='figs/HudsonAlpha_GM12878_CebpB_Score_CEBP_rsat')


## HCT116
p <- ggplot(rsat_HCT116$nonCEBP, 
            aes(score, colour=motifName, group=motifName )) +
        geom_density()
ggsave(file='figs/HudsonAlpha_HCT116_CebpB_Score_nonCEBP_rsat.pdf')
# plot motifs in chunks, each for 10 motifs
plotNMotifScore(rsat_HCT116$nonCEBP, 
                output_f ='figs/HudsonAlpha_HCT116_CebpB_Score_nonCEBP_rsat_chunk', 
                N=5)
p <- ggplot(rsat_HCT116$CEBP, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif:HCT116')
ggsave(file='figs/HudsonAlpha_HCT116_CebpB_Score_CEBP_rsat.pdf')

## K562
p <- ggplot(rsat_K562$nonCEBP, 
            aes(score, colour=motifName, group=motifName )) +
        geom_density()
ggsave(file='figs/HudsonAlpha_K562_CebpB_Score_nonCEBP_rsat.pdf')
# plot motifs in chunks, each for 10 motifs
plotNMotifScore(rsat_K562$nonCEBP, 
                output_f ='figs/HudsonAlpha_K562_CebpB_Score_nonCEBP_rsat_chunk', 
                N=5)
p <- ggplot(rsat_K562$CEBP, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif:K562')
ggsave(file='figs/HudsonAlpha_K562_CebpB_Score_CEBP_rsat.pdf')





## ------------------------------------------------------------------------
sessionInfo()


## ----knitIt,cache=TRUE, results='hide', message=FALSE, warning=FALSE-----
library(knitr)
purl("motifAnalysis.Rnw" ) # compile to tex
purl("motifAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("motifAnalysis.Rnw")


