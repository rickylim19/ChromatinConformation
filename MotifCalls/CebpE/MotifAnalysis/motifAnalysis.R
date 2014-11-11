
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


## ----Rcodes--------------------------------------------------------------
source('/home/ricky/Rlim/ChromatinConformation/MotifCalls/MotifCalls.R')
work_dir= '/home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/'


## ----DataPrep, echo=FALSE, cache=TRUE------------------------------------

Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist3kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance=3000, filterOutSingle=FALSE, 
    output_f = paste0(work_dir,
    'Output/Koeffler_BM_CebpE_NBM_BiclusterAssignment_compSorted4_dist3kb'))

Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist5kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance=5000,filterOutSingle=FALSE,
    paste0(work_dir,
    'Output/Koeffler_BM_CebpE_NBM_BiclusterAssignment_compSorted4_dist5kb'))

Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist8kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance=8000,filterOutSingle=FALSE,
    paste0(work_dir,
    'Output/Koeffler_BM_CebpE_NBM_BiclusterAssignment_compSorted4_dist8kb'))

Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist20kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance=20000,filterOutSingle=FALSE,
    paste0(work_dir,
    'Output/Koeffler_BM_CebpE_NBM_BiclusterAssignment_compSorted4_dist20kb'))

Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist50kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance=50000,filterOutSingle=FALSE,
    paste0(work_dir,
    'Output/Koeffler_BM_CebpE_NBM_BiclusterAssignment_compSorted4_dist50kb'))

# filter out the clusters with single peaks only
Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist3kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance = 3000, filterOutSingle = TRUE, 
    paste0(work_dir,
    'Output/',
    'Koeffler_BM_CebpE_NBM_BiclusterAssignment_SinglePeakFilteredOut_compSorted4_dist3kb'))
Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist5kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance = 5000, filterOutSingle = TRUE, 
    paste0(work_dir,
    'Output/',
    'Koeffler_BM_CebpE_NBM_BiclusterAssignment_SinglePeakFilteredOut_compSorted4_dist5kb'))

Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist8kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance = 8000, filterOutSingle = TRUE, 
    paste0(work_dir,
    'Output/',
    'Koeffler_BM_CebpE_NBM_BiclusterAssignment_SinglePeakFilteredOut_compSorted4_dist8kb'))

Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist20kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance = 20000, filterOutSingle = TRUE, 
    paste0(work_dir,
    'Output/',
    'Koeffler_BM_CebpE_NBM_BiclusterAssignment_SinglePeakFilteredOut_compSorted4_dist20kb'))

Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist50kb <-
 assignBiclustering(paste0(work_dir,
    'Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed'),
    distance = 50000, filterOutSingle = TRUE, 
    paste0(work_dir,
    'Output/',
    'Koeffler_BM_CebpE_NBM_BiclusterAssignment_SinglePeakFilteredOut_compSorted4_dist50kb'))


## ----DataExplorePlot, echo=FALSE, results='hide', cache=TRUE-------------

plotCountCluster(Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist3kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist3kb', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist3kb_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist5kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist5kb', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist5kb_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist8kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist8kb', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist8kb_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist20kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist20kb', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist20kb_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist50kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist50kb', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist50kb_countCluster.pdf'))

# Singlepeaks filteredout
plotCountCluster(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist3kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist3kb: SinglePeaks Filteredout', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist3kb_SinglePeakFilteredOut_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist5kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist50kb: SinglePeaks Filteredout', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist5kb_SinglePeakFilteredOut_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist8kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist50kb: SinglePeaks Filteredout', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist8kb_SinglePeakFilteredOut_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist20kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist50kb: SinglePeaks Filteredout', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist20kb_SinglePeakFilteredOut_countCluster.pdf'))
plotCountCluster(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist50kb$result$local, 
        title_p = 'Koeffler_BM_CebpE_NBM_comp4_dist50kb: SinglePeaks Filteredout', 
        output_f = paste0(work_dir, 
        'Output/Koeffler_BM_CebpE_NBM_comp4_dist50kb_SinglePeakFilteredOut_countCluster.pdf'))


## ----DataExploreTable----------------------------------------------------

## Table
Koeffler_BM_CebpE_Bicluster_NBM_table <- matrix(nrow=4, ncol=5) 
rownames(Koeffler_BM_CebpE_Bicluster_NBM_table) <- c('Direct', 'Indirect', 
                                                     'Noise', 'Total')
colnames(Koeffler_BM_CebpE_Bicluster_NBM_table) <- c('Dist3kb', 'Dist5kb', 'Dist8kb', 
                                                 'Dist20kb', 'Dist50kb')
Koeffler_BM_CebpE_Bicluster_NBM_table[, 1] <- 
    Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist3kb$numbindings 
Koeffler_BM_CebpE_Bicluster_NBM_table[, 2] <- 
    Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist5kb$numbindings 
Koeffler_BM_CebpE_Bicluster_NBM_table[, 3] <- 
    Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist8kb$numbindings 
Koeffler_BM_CebpE_Bicluster_NBM_table[, 4] <- 
    Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist20kb$numbindings 
Koeffler_BM_CebpE_Bicluster_NBM_table[, 5] <- 
    Koeffler_BM_CebpE_Bicluster_NBM_comp4_dist50kb$numbindings 

## Table: single peaks filtered out
Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table <- matrix(nrow=4, ncol=5) 
rownames(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table) <- 
    c('Direct', 'Indirect', 'Noise', 'Total')
colnames(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table) <- 
    c('Dist3kb', 'Dist5kb', 'Dist8kb', 'Dist20kb', 'Dist50kb')

Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table[, 1] <- 
    Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist3kb$numbindings 
Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table[, 2] <- 
    Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist5kb$numbindings 
Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table[, 3] <- 
    Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist8kb$numbindings 
Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table[, 4] <- 
    Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist20kb$numbindings 
Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table[, 5] <- 
    Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_comp4_dist50kb$numbindings 



## ----biclusteringNBMTable, echo=FALSE, results='asis'--------------------

print(xtable(Koeffler_BM_CebpE_Bicluster_NBM_table, 
             caption='Clustering Peaks in Different Cluster Distance'))
print(xtable(Koeffler_BM_CebpE_BiclusterSinglePeakFilteredOut_NBM_table, 
             caption='Clustering Peaks in Different Cluster Distance: Single Peaks Filteredout'))


## ----MotifDataPrep, cache=TRUE-------------------------------------------
# awk -F'\t' '{print $1"\t"$2-500"\t"$3+500"\t"$5}' 
# Input/Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed 
# > Output/Koeffler_BM_CebpE_Peaks_1kb.bed


## ----CentDistP53---------------------------------------------------------
work_dir_centdist = '/home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/'
Koeffler_BM_ChIPseq_CebpE_Bicluster_compSorted_dist50kb_direct <- 
    read.delim(paste0(work_dir,'Output/CentDist/',
    'Koeffler_BM_CebpE_NBM_BiclusterAssignment_SinglePeakFilteredOut_compSorted4_dist50kb_direct/',
    'result.txt'))
Koeffler_BM_ChIPseq_CebpE_Bicluster_compSorted_dist50kb_indirect <- 
    read.delim(paste0(work_dir,'Output/CentDist/',
    'Koeffler_BM_CebpE_NBM_BiclusterAssignment_SinglePeakFilteredOut_compSorted4_dist50kb_indirect/',
    'result.txt'))


## ----engine='bash'-------------------------------------------------------
# pwd: /home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/Output/
#      Rsat/peak-motifs.2014-11-03.083858_2014-11-03.083858_XPVNnt/results/sites

# skip the lines starting with '#' and ';'
# awk 'NF && $1!~/^#|^;/' peak-motifs_all_motifs_seqcoord.tab > peak_motifs_sites.tab

# for the jaspar db
# pwd: /home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/Output/Rsat/
#      peak-motifs.2014-11-03.083858_2014-11-03.083858_XPVNnt/results/discovered_vs_db

# skip the lines starting with '#' and ';'
# awk 'NF && $1!~/^#|^;/' peak-motifs_motifs_vs_db_jaspar_core_vertebrates.tab > 
#   peak-motifs_db_jaspar.tab


## ----DataPreProcessingMotif, cache=TRUE----------------------------------
#motif_dir = '/home/ricky/Rlim/ChromatinConformation/MotifCalls/'

Koeffler_BM_CebpE_Peaks_MotifSites <- read.table(paste0(work_dir, 
    'Output/Rsat/peak-motifs.2014-11-03.083858_2014-11-03.083858_XPVNnt/results/',
    'sites/peak_motifs_sites.tab'))

# remove the duplicated sites
Koeffler_BM_CebpE_Peaks_MotifSites_Unique <- Koeffler_BM_CebpE_Peaks_MotifSites[
            !duplicated(Koeffler_BM_CebpE_Peaks_MotifSites$V1),] 
Koeffler_BM_CebpE_Peaks_MotifSites_Unique <- 
            Koeffler_BM_CebpE_Peaks_MotifSites_Unique[, c('V1', 'V3')] 

Koeffler_BM_CebpE_Peaks_MotifSites_Unique <- 
            addChrCoordinates(Koeffler_BM_CebpE_Peaks_MotifSites_Unique)

head(Koeffler_BM_CebpE_Peaks_MotifSites_Unique)

# get the peak score
Koeffler_BM_CebpE_Peaks_1kb  <- read.table(paste0(work_dir, 
                                             'Input/Koeffler_BM_CebpE_Peaks_1kb.bed'))
colnames(Koeffler_BM_CebpE_Peaks_1kb) <- c('chr', 'start', 'end', 'score')
Koeffler_BM_CebpE_Motifs<- 
    merge(Koeffler_BM_CebpE_Peaks_MotifSites_Unique, Koeffler_BM_CebpE_Peaks_1kb, 
          by=c('chr', 'end'))
head(Koeffler_BM_CebpE_Motifs)
Koeffler_BM_CebpE_Motifs <- Koeffler_BM_CebpE_Motifs[, 
          c('chr', 'start', 'end', 'score', 'V3')]
colnames(Koeffler_BM_CebpE_Motifs) <- 
    c('chr', 'start', 'end', 'score', 'motifId')

# get the motif name from jaspar database
Koeffler_BM_CebpE_Jaspar <- read.table(paste0(work_dir, 
    'Output/Rsat/peak-motifs.2014-11-03.083858_2014-11-03.083858_XPVNnt/results/',
    '/discovered_vs_db/peak-motifs_db_jaspar.tab'))
Koeffler_BM_CebpE_Jaspar <- Koeffler_BM_CebpE_Jaspar[, c('V1', 'V4')]
colnames(Koeffler_BM_CebpE_Jaspar) <- c('motifId', 'motifName')
Koeffler_BM_CebpE_Motifs <- merge(Koeffler_BM_CebpE_Motifs, Koeffler_BM_CebpE_Jaspar, 
                        by=c('motifId','motifId'))
Koeffler_BM_CebpE_Motifs <- Koeffler_BM_CebpE_Motifs[, 
                        c('chr', 'start', 'end', 'score', 'motifId', 'motifName')]

# get only CEBP motif
CEBP_motif <- Koeffler_BM_CebpE_Motifs$motifName %in% c('CEBPA', 'CEBPB', 'Cebpa')
Koeffler_BM_CebpE_CEBP_motif <- Koeffler_BM_CebpE_Motifs[CEBP_motif, 
                                                          c('motifName', 'score')] 
Koeffler_BM_CebpE_CEBP_motif[, 'motifName'] <- rep('CEBP', 
                                                   nrow(Koeffler_BM_CebpE_CEBP_motif))
Koeffler_BM_CebpE_nonCEBP_motif <- Koeffler_BM_CebpE_Motifs[!CEBP_motif, 
                                                          c('motifName', 'score')] 




## ----plotNonCEBPMotif, cache=TRUE, results='hide', echo=FALSE------------
# the whole motifs in one plot
p <- ggplot(Koeffler_BM_CebpE_nonCEBP_motif, 
            aes(score, colour=motifName, group=motifName )) +
        geom_density()
ggsave(file='figs/Koeffler_BM_CebpE_Score_nonCEBP.pdf')

# plot motifs in chunks, each for 10 motifs
plotNMotifScore(Koeffler_BM_CebpE_nonCEBP_motif, 
                output_f ='figs/Koeffler_BM_CebpE_Score_nonCEBP_chunk', 
                N=5)



## ----plotCEBPMotifRsat, cache=TRUE, results='hide', echo=FALSE-----------
p <- ggplot(Koeffler_BM_CebpE_CEBP_motif, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif')
ggsave(file='figs/Koeffler_BM_CebpE_Score_CEBP.pdf')



## ----engine='bash'-------------------------------------------------------

## for a list of motifs (10 motifs output of meme)
# for i in $(ls motif_*_fasta.txt); do awk '/^>/' $i | sed 's/^>//g' | awk '{print $1}'> 
# ''`basename $i .txt`.site''; done

## add filename without _fasta.site to each line
# for i in $(ls *.site);do awk '{print $1"\t"FILENAME}' $i | sed 's/_fasta.site//g' > 
# "`basename $i .site`.msite"; done

## concatenate all the motif sites
# cat *.msite > allmotif.msite

## get the peak with width of 0.5kb
# awk -F``\t'' '{print $1``\t''$2-250``\t''$3+249``\t''$5}' 
#  Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed>Koeffler_BM_CebpE_Peaks_0.5kb.bed
# awk -F"\t" '{print $1"\t"$2-500"\t"$3+499"\t"$5}' 
#  Koeffler_BM_CebpE_NBM_ModelAssignment_compSorted4.bed>Koeffler_BM_CebpE_Peaks_1kb.bed



## ----memeDataOutput, cache=TRUE------------------------------------------
# note that ChIP meme uses ChIP regions of 0.5kb
peak_0.5kb <- read.table(paste0(work_dir, 
                                'Input/Koeffler_BM_CebpE_Peaks_0.5kb.bed'))
colnames(peak_0.5kb) <- c('chr', 'start', 'end', 'score')


all_motifs <- read.table(paste0(work_dir, 
                             'Output/ChipMeme/MotifMeme/allmotif.msite'))
head(all_motifs)
cebp_motif <- all_motifs[all_motifs$V2 == 'motif_2',]

allmotif_score <- getPeakScore(all_motifs, peak_0.5kb)
cebpmotif_score <- getPeakScore(cebp_motif, peak_0.5kb)
colnames(allmotif_score) <- c('chr', 'start', 'end', 'id', 'score')
colnames(cebpmotif_score) <- c('chr', 'start', 'end', 'id', 'score')

#motif_1_score <- getPeakScore(motif_1, peak_0.5kb)
#motif_2_score <- getPeakScore(motif_2, peak_0.5kb)
head(allmotif_score)
head(cebpmotif_score)



## ----plotCEBPMotifMeme, cache=TRUE, results='hide', echo=FALSE-----------
p <- ggplot(cebpmotif_score, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of CEBP Motif: Meme-Chip')
ggsave(file='figs/Koeffler_BM_CebpE_Score_CEBPmotif.pdf')

p <- ggplot(allmotif_score, aes(score, colour=id, group=id)) +
        geom_density()  
ggsave(file='figs/Koeffler_BM_CebpE_Score_Allmotif.pdf')


## ----junChIPSites, engine = 'bash'---------------------------------------
## convert jun sites in gtf to bed
# ./gtfMeme2Bed.py jun_sites.gtf > jun_sites.bed  

## intersect the bed files with the ChIPseq bed
# bedtools intersect 
# -b /home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/Output/
#    ChipMeme/jun_sites.bed 
# -a /home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/Input/
#    Koeffler_BM_CebpE_Peaks_0.5kb.bed 
# -u > /home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/Output/
#    ChipMeme/MotifMeme/jun_ChIPSites.bed

## ----plotJun-------------------------------------------------------------

jun_chipScore <- read.table(paste0(work_dir, 'Output/ChipMeme/MotifMeme/jun_ChIPSites.bed'))
head(jun_chipScore)
colnames(jun_chipScore) <- c('chr', 'start', 'end', 'score')
p <- ggplot(jun_chipScore, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of JUN Motif: Meme-Chip')
ggsave(file='figs/Koeffler_BM_CebpE_Score_JUN_MemeChip.pdf')


## ----engine='bash'-------------------------------------------------------
## grep only unique sites for homo and heterodimer motif sites
# grep -v -x -f cebpPerfectHeterodimer.site cebpPerfectHomodimer.site > 
# cebpPerfectHomodimerOnly.site
# grep -v -x -f cebpPerfectHomodimer.site cebpPerfectHeterodimer.site > 
# cebpPerfectHeterodimerOnly.site

## add the filename to the 2nd column
# for i in $(ls *.site);do awk '{print $1``\t''FILENAME}' $i | sed 's/.site//g' > ```basename $i .site`.msite''; done

## combine the homo and heterodimer motif sites
# cat cebpPerfectHeterodimerOnly.msite cebpPerfectHomodimerOnly.msite > 
# cebpPerfectHomoHeterodimerOnly.msite


## ----homoheterodimerDataPrep, cache=TRUE---------------------------------
homohetero_dir=paste0('/home/ricky/Rlim/ChromatinConformation/MotifCalls/CebpE/',
                      'Output/ChipMeme/MotifMeme/HomoHetero/')

homoheterodimer_msite <- read.table(paste0(homohetero_dir, 
                                    'cebpPerfectHomoHeterodimerOnly.msite' ))

homodimer_motif <- homoheterodimer_msite[
                        homoheterodimer_msite$V2 == 'cebpPerfectHomodimerOnly',]
nHomodimer <- nrow(homodimer_motif)
heterodimer_motif <- homoheterodimer_msite[
                        homoheterodimer_msite$V2 == 'cebpPerfectHeterodimerOnly',]
nHeterodimer <- nrow(heterodimer_motif)


homoheterodimer_score <- getPeakScore(homoheterodimer_msite, peak_0.5kb)
homodimer_score <- getPeakScore(homodimer_motif, peak_0.5kb)
heterodimer_score<- getPeakScore(heterodimer_motif, peak_0.5kb)
colnames(homoheterodimer_score) <- c('chr', 'start', 'end', 'id', 'score')
colnames(homodimer_score) <- c('chr', 'start', 'end', 'id', 'score')
colnames(heterodimer_score) <- c('chr', 'start', 'end', 'id', 'score')



## ----plotHomoHeteroCebp, cache=TRUE, results='hide', echo=FALSE----------
p <- ggplot(homodimer_score, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of Perfect Homodimer CEBP Motif')
ggsave(file='figs/Koeffler_BM_CebpE_Score_PerfectHomoCebp.pdf')

p <- ggplot(heterodimer_score, aes(x=score))+ 
        geom_histogram(aes(y = ..density..), binwidth=1, fill='#969696') + 
        geom_line(stat='density', alpha=0.6, size=1.5) + 
        ggtitle('ChIPseq Score Distribution of Perfect Heterodimer CEBP Motif')
ggsave(file='figs/Koeffler_BM_CebpE_Score_PerfectHeteroCebp.pdf')

p <- ggplot(homoheterodimer_score, aes(score, colour=id, group=id)) +
        geom_density()  
ggsave(file='figs/Koeffler_BM_CebpE_Score_HomoHeteroCebp.pdf')


## ------------------------------------------------------------------------
sessionInfo()


## ----knitIt, cache=TRUE, results='hide', message=FALSE, warning=FALSE----
library(knitr)
purl("motifAnalysis.Rnw" ) # compile to tex
purl("motifAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("motifAnalysis.Rnw")


