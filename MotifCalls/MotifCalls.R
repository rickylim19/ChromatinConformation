library(ggplot2)
library(scales)
library(xtable)
library(plyr)

assignCluster <- function(bed_f, distance=100, filterOutSingle=FALSE){

    ###################################################################
    #                                                                 #
    # Returns bed table with the cluster ids (V6)                     #
    #  Clusters based on their interdistance                          #
    #  Interdistance: cur.start - prev.end                            #
    #                                                                 #
    # > assignCluster('Test/inputCluster.bed') (dummy example         #
    #                                                                 #
    #    V1  V2   V3    v4  v5  V6                                    #
    # 1 chr1 100  150   1   30    1                                   #
    # 2 chr1 200  250   3   100   1                                   #
    # 3 chr1 501 1000   2   75    2                                   #
    # 4 chr2 100  150   3   200   3                                   #
    # 5 chr2 200  250   2   80    3                                   #
    # 6 chr2 501 1000   2   85    4                                   #
    #                                                                 #
    # # chr start end comp score local                                #
    # local: cluster num                                              #
    #                                                                 #
    #                                                                 #
    ###################################################################

    bed <- read.table(bed_f)
    cluster_list <- 1 #start with the cluster_id of 1
    for (i in 2: nrow(bed)){
        cur_start <- bed[,2][i]; prev_end<- bed[,3][i-1]
        cur_chr <- bed[,1][i]; prev_chr <- bed[,1][i-1]

        cluster_dist <- (cur_start - prev_end)
        if (cluster_dist <= distance && prev_chr == cur_chr){
            cluster_list[i] <- cluster_list[i-1]
        }
        else{
            cluster_list[i] <- cluster_list[i-1] + 1
        }
    }
    bed$V6 <- cluster_list
    if (filterOutSingle){
        count_tab <- as.data.frame(table(bed$V6))
        numRow <- nrow(bed)
        isSingle <- c()
        for (i in 1:numRow){
            if (count_tab[count_tab$Var1 == bed[i,6], 'Freq'] == 1){
                isSingle <- append(isSingle, TRUE)}
            else {isSingle<- append(isSingle, FALSE)}

        }
        bed$isSingle <- isSingle
        bed <- bed[bed$isSingle == FALSE, ]
        bed <- subset(bed, select=-isSingle)
        return (bed)
    }
    else {
        return (bed)
    }
}

assignBiclustering  <- function(bed, distance=100, filterOutSingle, output_f){

    #############################################################################
    #                                                                           #
    # Write 4 tables after assigning the direct(D), indirect(I), and            #
    #  noise clusters and returns a list containing :                           #
    #  - result: direct and indirect bindings                                   #
    #  - result_direct: direct bindings only                                    #
    #  - result_indirect: indirect bindings only                                #
    #  - result_noise: noise bindings only                                      #
    #  - result_numbindings: # for direct, indirect, noise, and total bindings  #
    #                                                                           #
    #  Header: #chr start end bicluster                                         #
    # 1). Table with both direct and indirect bindings                          #
    # 2). Table with only peaks assigned to direct bindings                     #
    # 3). Table with only peaks assigned to indirect bindings                   #
    # 4). Table with only peaks assigned to noise bindings                      #
    #                                                                           #
    # The biclustering is based on:                                             #
    #  - Local clustering: cluster peaks that are in close proximity            #
    #  - Assignment of the higher component within the locals with direct and   #
    #    the lower component with indirect bindings                             #
    #  - If the cluster contains only the lowest comp (1) then N is assigned    #
    #    *N: noise                                                              #
    #                                                                           #
    #  e.g:                                                                     #
    # #chr    start           end             comp    score   local   bicluster #
    # chr1    30832929        30832930        1       20      65      I         #
    # chr1    30872925        30872926        1       31      65      I         #
    # chr1    30910456        30910457        2       49      65      D         #
    # chr1    30920978        30920979        1       35      65      I         #
    # chr1    30921389        30921390        2       40      65      D         #
    # chr1    30930339        30930340        2       61      65      D         #
    # chr1    30948870        30948871        1       28      65      I         #
    # chr1    30950137        30950138        2       43      65      D         #
    #                                                                           #
    # > bicluster_3kb<-assignBiclustering(input.bed ,distance=3000, output.bed  #
    #                                                                           #
    #############################################################################

    # local clustering: cluster the bed according the distance in between
    #comp_local <- bedClusterPeaks(bed=bed, distance=distance)
    comp_local <- assignCluster(bed_f=bed, distance=distance, filterOutSingle)
    colnames(comp_local) <- c('chr', 'start', 'end', 'comp', 'score', 'local')
    comp_locallist<- split(comp_local, comp_local$local)
    f <- function(x){
            #get the highest component for each local cluster
            max_comp <- max(x$comp)
            apply(x, 1, function(x) ifelse((max_comp) == 1, 'N',
                                           ifelse(x[4]==max_comp, 'D', 'I')))
            }
    bicluster <- unlist(lapply(comp_locallist, FUN=f))
    comp_local$bicluster<- bicluster

    #result <- comp_local[, c('chr', 'start', 'end', 'bicluster')]
    result_list <- vector('list', length = 4)

    result <- comp_local
    result_direct<- result[result$bicluster=='D',]
    result_indirect <- result[result$bicluster=='I',]
    result_noise<- result[result$bicluster=='N',]

    write.table(result, paste0(output_f, '.bed'), quote=F, sep='\t',
                row.names=F, col.names=F)
    write.table(result_direct, paste0(output_f,'_direct', '.bed'), quote=F, sep='\t',
                row.names=F, col.names=F)
    write.table(result_indirect, paste0(output_f, '_indirect', '.bed'), quote=F, sep='\t',
                row.names=F, col.names=F)
    write.table(result_noise, paste0(output_f, '_noise', '.bed'), quote=F, sep='\t',
                row.names=F, col.names=F)

    # return the number of direct, indirect, and noise bindings
    list_numbindings <- c(nrow(result_direct),
                         nrow(result_indirect),
                         nrow(result_noise),
                         nrow(result))

    result_list$result<- result
    result_list$direct <- result_direct
    result_list$indirect <- result_indirect
    result_list$noise <- result_noise
    result_list$numbindings <- list_numbindings


    return(result_list)

}


plotCountCluster <- function(cluster, title_p, output_f){

    ###################################################################################
    # Returns a frequency plot (in percent) of each number of peak in a cluster       #
    #                                                                                 #
    # > plotCountCluster(Koeffler_BM_Cebpe_Bicluster_NBM_comp4_dist3kb$result$local,  #
    #    title_p = 'Koeffler_BM_Cebpe_NBM_comp4_dist3kb',                             #
    #    output_f = paste0(Input_dir,                                                 #
    #    'Output/Plots/Koeffler_BM_Cebpe_NBM_comp4_dist3kb_countCluster.pdf'))        #
    ###################################################################################

    countCluster <- as.data.frame(table(as.data.frame(table(cluster))$Freq))
    p <- ggplot(countCluster,
                aes(x = Var1, y = Freq/sum(Freq))) + geom_bar(stat='identity')
    p + scale_y_continuous(labels=percent) + xlab('Peak') + ylab('') + ggtitle(title_p)
    ggsave(output_f)

}

createMotifClouds <- function(motif, motif_freq, output_n){

    ###########################################################################
    #                                                                         #
    # Produces the [Motif]wordcloud.pdf                                       #
    #                                                                         #
    # > createMotifClouds(Koeffler_BM_Cebpe_CentDistFamily_direct300,         #
    #                     Koeffler_BM_Cebpe_CentDistScore_direct300,          #
    #                    'Koeffler_BM_Cebpe_CentDist_direct300_WordCloud')    #
    #                                                                         #
    ###########################################################################

    pdf(paste0('MotifCalls/Output/Plots/', output_n, '.pdf'), 
        useDingbats=FALSE)
    wordcloud(motif, motif_freq,
              random.order=FALSE, random.color=FALSE, 
              colors=brewer.pal(8, 'Dark2'), rot.per=0)
    dev.off()
}

addChrCoordinates <- function(x){

    ###################################################################
    # add chr end columns to dataframe x                              #
    # Note that the fasta format for the motif sites, e.g:            #
    # - mm10_chr10_81266663_81267163_+_site_1                         #
    #                                                                 #
    # > addChrCoordinates(Koeffler_BM_Cebpe_Peaks_MotifSites_Unique)  #
    #                                                                 #
    ###################################################################

    x <- as.matrix(x)
    chr_ <- apply(x,1, function(x) ldply(strsplit(x['V1'], '_'))$V2)
    end_ <- apply(x,1, function(x) ldply(strsplit(x['V1'], '_'))$V4)
    x <- as.data.frame(x)
    x$chr <- chr_
    x$end <- end_
    return(x)
}

plotNMotifScore <- function(motif_tab, output_f, N){

    #############################################################################
    # Returns incrementally each a plot of ChIPseq score for N motifs           #
    #                                                                           #
    # > plotNMotifScore(Koeffler_BM_Cebpe_nonCEBP_motif,                        #
    #        output_f = paste0('/home/ricky/Rlim/ChromatinConformation/',       #
    #        'MotifCalls/Output/Plots/Koeffler_BM_Cebpe_Score_nonCEBP_chunk'),  #
    #        N=5)                                                               #
    #                                                                           #
    #############################################################################

    motif_names <- unique(motif_tab$motifName)
    motif_chunks <- split(motif_names, ceiling(seq_along(motif_names)/N))
    for (i in 1:length(motif_chunks)){
        tab <- motif_tab[motif_tab$motifName %in% motif_chunks[[i]],]
        # add cebp motif
        tab <- rbind(tab, Koeffler_BM_CebpE_CEBP_motif)
        p <- ggplot(tab, aes(score, colour=motifName, group=motifName )) +
                geom_density()
        ggsave(file = paste0(output_f,i,'.pdf'))
    }
}

getPeakScore <- function(motif, peak){

    ############################################################
    # Returns the table with ChIPseq score                     #
    #                                                          #
    # Input:                                                   #
    # - motif                                                  #
    #    Table, each row looks like:                           #
    #    mm10_chr4_117886122_117886622_+_site_1                #
    # - peak (bed file)                                        #
    #    chr start end score                                   #
    #                                                          #
    # Output (dummpy example):                                 #
    #   chr start end sore                                     #
    #   chr10   100 120 21                                     #
    #   chr10   200 220 34                                     #
    #                                                          #
    # > motif_1_score <- getPeakScore(motif_1, peak_0.5kb)     #
    ############################################################

    result <- addChrCoordinates(motif)
    result <- merge(result, peak, by = c('chr', 'end'))
    result <- result[, c('chr', 'start', 'end', 'V2', 'score')]
    colnames(result) <- c('chr', 'start', 'end', 'V2', 'score')
    return (result)

}
getCebpAllMotif <- function(peakScore_f, allmotif_f, cebp_motif){

    #####################################################################################
    # Returns                                                                           #
    #                                                                                   #
    # peakScore_f : path to the bed file with the score                                 #
    # e.g of peak_score_f:                                                              #
    # chr1    1711747 1712248 80                                                        #
    # chr1    3573209 3573710 66                                                        #
    # chr1    3593444 3593945 63                                                        #
    #                                                                                   #
    # allmotif_f : path to the all motif file from chip meme                            #
    # hg19_chr15_60927876_60928376_+_site_1 motif_1                                     #
    # hg19_chr3_101617466_101617966_+_site_1  motif_1                                   #
    # hg19_chr9_114840750_114841250_+_site_1  motif_1                                   #
    #                                                                                   #
    # cebp_motif: the number of the motif reflects the cebp (TGA--TCA)                  #
    # e.g: motif_3 means motif number 3                                                 #
    #                                                                                   #
    #                                                                                   #
    #                                                                                   #
    # meme_chipScore_GM12878 <- getCebpAllMotif(peakScore_f=paste0(work_dir,            #
    #                    'Input/HudsonAlpha_GM12878_CebpB_Peaks_0.5kb.bed'),            #
    #                    allmotif_f=paste0(work_dir,                                    #
    #                    'Output/ChipMeme_HudsonAlpha_GM12878/Motif/allmotif.msite'),   #
    #                    cebp_motif= 'motif_3')                                         #
    #                                                                                   #
    #####################################################################################

    result_list <- vector('list', length=2)

    peak <- read.table(peakScore_f)
    colnames(peak) <- c('chr', 'start', 'end', 'score')
    all_motifs<- read.table(allmotif_f)
                                
    cebp_motif <- all_motifs[all_motifs$V2 == cebp_motif,]

    allmotif_score<- getPeakScore(all_motifs, peak)
    cebpmotif_score<- getPeakScore(cebp_motif, peak)
    colnames(allmotif_score) <- c('chr', 'start', 'end', 'id', 'score')
    colnames(cebpmotif_score) <- c('chr', 'start', 'end', 'id', 'score')

    result_list$allmotif <- allmotif_score 
    result_list$cebpmotif <- cebpmotif_score
    return (result_list)


}

getCebpNon_rsat <- function(input_f, peak_f, jaspar_f){

    ##############################################################################################
    # input_f: awk 'NF && $1!~/^#|^;/' peak-motifs_all_motifs_seqcoord.tab`                      #
    # peak_f: bed file with ChIPscore                                                            #
    # jaspar_f: awk 'NF && $1!~/^#|^;/' peak-motifs_motifs_vs_db_jaspar_core_vertebrates.tab     #
    #                                                                                            #
    # > getCebpNon_rsat(input_f=paste0(work_dir,                                                 #
    #    'Output/Rsat/HudsonAlpha_HCT116_CebpB_Peaks_1kb/results/sites/peak_motifs_sites.tab'),  #
    #                            peak_f=paste0(work_dir,                                         #
    #                                        'Input/HudsonAlpha_HCT116_CebpB_Peaks_1kb.bed'),    #
    #                            jaspar_f=paste0(work_dir,                                       #
    #    'Output/Rsat/HudsonAlpha_HCT116_CebpB_Peaks_1kb/results',                               #
    #    '/discovered_vs_db/peak-motifs_db_jaspar.tab'))                                         #
    #                                                                                            #
    ##############################################################################################

    result_list <- vector('list', length=2)

    MotifSites <- read.table(input_f)

    # remove the duplicated sites
    MotifSites_Unique <- MotifSites[!duplicated(MotifSites$V1),] 
    MotifSites_Unique <- MotifSites_Unique[, c('V1', 'V3')] 

    MotifSites_Unique <- addChrCoordinates(MotifSites_Unique)

    # get the peak score
    Peaks_1kb  <- read.table(peak_f)
    colnames(Peaks_1kb) <- c('chr', 'start', 'end', 'score')
    Motifs<- merge(MotifSites_Unique, Peaks_1kb, by=c('chr', 'end'))
    Motifs <- Motifs[, c('chr', 'start', 'end', 'score', 'V3')]
    colnames(Motifs) <- c('chr', 'start', 'end', 'score', 'motifId')

    # get the motif name from jaspar database
    Jaspar <- read.table(jaspar_f)
    Jaspar <- Jaspar[, c('V1', 'V4')]
    colnames(Jaspar) <- c('motifId', 'motifName')
    Motifs <- merge(Motifs, Jaspar, by=c('motifId','motifId'))
    Motifs <- Motifs[, c('chr', 'start', 'end', 'score', 'motifId', 'motifName')]


    # get only CEBP motif
    CEBP <- Motifs$motifName %in% c('CEBPA', 'CEBPB', 'Cebpa', 'Cebpb')
    CEBP_motif <- Motifs[CEBP, c('motifName', 'score')] 
    CEBP_motif[, 'motifName'] <- rep('CEBP', nrow(CEBP_motif))
    nonCEBP_motif <- Motifs[!CEBP, c('motifName', 'score')] 

    result_list$CEBP <- CEBP_motif
    result_list$nonCEBP <- nonCEBP_motif

    return(result_list)
}



