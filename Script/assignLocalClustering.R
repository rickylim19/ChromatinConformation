#!/usr/bin/env Rscript

options(options = 3)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mixtools"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("plyr"))

### main arguments ###
option_list <- list(
    make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                help = "Print extra output [default]"),
    make_option(c("-q", "--quietly"), action = "store_false",
                dest = "verbose", help = "Print little output"),
    make_option(c("-s", "--seed"), type = "integer", default = 12345,
                help = "set the seed number for reproducibility",
                metavar = "number"),
    make_option("--core", type = "integer", default = 4,
                help = "set the number of cores for parallelization",
                metavar = "number"),
    make_option("--distance", type = "integer", 
                help = "distance in bp within a cluster", default= '3000',
                metavar = "number"),
    make_option("--filter", type = "logical", 
                help = "filter out cluster that contains only the lowest score", default= T,
                metavar = "logical"),
    make_option("--input_dir", action = "store", type = "character", 
                help = "input directory",
                metavar = "character"),
    make_option("--output_dir", action = "store", type = "character", 
                help = "output directory",
                metavar = "character")
    )

assignCluster <- function(bed, distance=3000, filterOutSingle=TRUE){

    ###################################################################
    #                                                                 #
    # Returns bed table with the cluster ids (V6)                     #
    #  Clusters based on their interdistance                          #
    #  Interdistance: cur.start - prev.end                            #
    #                                                                 #
    # > assignCluster('Test/inputCluster.bed', distance=50)           #
    #                                                                 #
    #    chrom  start   end    score  comp  cluster                   #
    # 1 chr1 100  150   30  1  1                                      #
    # 2 chr1 200  250   100 3  1                                      #
    # 3 chr1 501 1000   75  2  2                                      #
    # 4 chr2 100  150   200 3  3                                      #
    # 5 chr2 200  250   80  2  3                                      #
    # 6 chr2 501 1000   85  2  4                                      #
    #                                                                 #
    # # chr start end comp score local                                #
    # local: cluster num                                              #
    #                                                                 #
    #                                                                 #
    ###################################################################

    #bed <- read.table(bed_f)
    cluster_list <- 1 #start with the cluster_id of 1
    for (i in 2: nrow(bed)){
        cur_start <- bed[,'start'][i]; prev_end<- bed[,'end'][i-1]
        cur_chr <- bed[,'chrom'][i]; prev_chr <- bed[,'chrom'][i-1]

        cluster_dist <- (cur_start - prev_end)
        if (cluster_dist <= distance && prev_chr == cur_chr){
            cluster_list[i] <- cluster_list[i-1]
        }
        else{
            cluster_list[i] <- cluster_list[i-1] + 1
        }
    }
    bed$cluster <- cluster_list
    if (filterOutSingle){
        count_tab <- as.data.frame(table(bed$cluster))
        numRow <- nrow(bed)
        isSingle <- c()
        for (i in 1:numRow){
            if (count_tab[count_tab$Var1 == bed[i,'cluster'], 'Freq'] == 1){
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

apply_localClustering <- function(x, distance=3000, filterOutSingle=T, input_dir, output_dir){

    #############################################################################
    #                                                                           #
    # Write 1 table after assigning the direct(D), indirect(I), and             #
    #  noise clusters and returns a list containing :                           #
    #  - result: direct and indirect bindings                                   #
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

    filename <- removeExt(x)
    #inputFile <- read.table(paste0(input_dir, x), header=T, sep='\t') 
    inputFile <- read.table(paste0(input_dir, x), sep='\t') 
    colnames(inputFile) <- c('chrom', 'start', 'end', 'score', 'comp')
    # local clustering: cluster the bed according the distance in between
    #comp_local <- bedClusterPeaks(bed=bed, distance=distance)
    comp_local <- assignCluster(bed=inputFile, distance=distance, filterOutSingle)
    colnames(comp_local) <- c('chrom', 'start', 'end', 'score', 'comp', 'local')
    comp_locallist<- split(comp_local, comp_local$local)
    f <- function(x){
            #get the highest component for each local cluster
            max_comp <- max(x$comp)
            apply(x, 1, function(x) ifelse((max_comp) == 1, 'N',
                                           ifelse(x[5]==max_comp, 'D', 'I')))
            }
    bicluster <- unlist(lapply(comp_locallist, FUN=f))
    comp_local$bicluster<- bicluster

    #result <- comp_local[, c('chr', 'start', 'end', 'bicluster')]
    result_list <- vector('list', length = 4)

    result <- comp_local
    result_direct<- result[result$bicluster=='D',]
    result_indirect <- result[result$bicluster=='I',]
    result_noise<- result[result$bicluster=='N',]

    write.table(result, paste0(output_dir, filename, '_', distance, 'bp_localClusters.bed'), 
                quote=F, sep='\t',row.names=F, col.names=F)
    write.table(result_direct, 
                paste0(output_dir, filename, '_', distance, 'bp_localClusters_direct.bed'), 
                quote=F, sep='\t',row.names=F, col.names=F)
    write.table(result_indirect, 
                paste0(output_dir, filename, '_', distance, 'bp_localClusters_indirect.bed'), 
                quote=F, sep='\t',row.names=F, col.names=F)
    write.table(result_noise, 
                paste0(output_dir, filename, '_', distance, 'bp_localClusters_noise.bed'), 
                quote=F, sep='\t',row.names=F, col.names=F)
}

### call the command with apply in parallel)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(file.path(opt$output_dir), showWarnings = F, recursive = T)
#filePattern = paste0(opt$resolution, 'bin-.*peaks.bed')
files <- list.files(path = opt$input_dir, pattern = '*_components.bed') 
set.seed(opt$seed)
invisible(mclapply(files, apply_localClustering, opt$distance, opt$filter, opt$input_dir, opt$output_dir, mc.cores = opt$core))
