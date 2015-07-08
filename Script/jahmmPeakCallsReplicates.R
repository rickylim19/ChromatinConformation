#!/usr/bin/env Rscript

options(options = 3)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Rjahmm"))
suppressPackageStartupMessages(library("dplyr"))

### main arguments ###
option_list <- list(
    make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                help = "Print extra output [default]"),
    make_option(c("-q", "--quietly"), action = "store_false",
                dest = "verbose", help = "Print little output"),
    make_option(c("-s", "--seed"), type = "integer", default = 12345,
                help = "set the seed number for reproducibility",
                metavar = "number"),
    make_option("--xrange", type = "integer", default = 100000,
                help = "set the number of rows plotted in the plotfit",
                metavar = "number"),
    make_option("--resolution", type = "character", default = "300",
                help = "bin resolution",
                metavar = "character"),
    make_option("--input_file", action = "store", type = "character", 
                help = "input file .jahmm",
                metavar = "character"),
    make_option("--output_dir", action = "store", type = "character", 
                help = "output directory",
                metavar = "character")
    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list))

dat <- read.table(opt$input_file, sep = '\t', header = T)
datDim <- dim(dat)
coordinates <- dat[, 1:3] #chr start end
readDat <- dat[, c(1, 4:datDim[2])] #chr protX ...
rm(dat)

colnames(readDat) <- toupper(colnames(readDat))
proteins <- colnames(readDat)

proteins <- unique(gsub("_REP[0-9]", "", proteins[grepl("_REP[0-9]", proteins)]))

plotFit <- function(chipProfile, fit, xrange, output_dir){

    ###########################################################################################
    #                                                                                         #
    # Returns the pdf plot of the chipPeaks coloured with identified targets                  #
    #                                                                                         #
    # Requires:                                                                               #
    # - chipProfile looks like:                                                               #
    # chr Input CebpE                                                                         #
    # chr1 0    0                                                                             #
    # chr1 0    0                                                                             #
    #                                                                                         #
    # - fit                                                                                   #
    # fit object returned by jahmm                                                            #
    # - xrange                                                                                #
    # range of the chromosome positions (range of the rows of the dataframe from chipProfile) #
    # - output_dir                                                                            #
    # output to return the pdf coupled with its name                                          #
    #                                                                                         #
    ###########################################################################################

    pdf(paste0(output_dir))
    par(mfrow=c(2,1))
    # plot the protein
    plot(chipProfile[xrange, 3], 
        #col = as.factor(fit$path), type = 'h', 
        col = as.factor(fit), type = 'h', 
        xlab='Position on Chr1 (kb)', ylab='read', main='Protein')
    legend('topright', c('non-target', 'intermediate', 'target'), 
        lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = c('black', 'red', 'green'), bty='n')

    # plot the input
    plot(chipProfile[xrange, 2], 
        #col = as.factor(fit$path), type = 'h', 
        col = as.factor(fit), type = 'h', 
        xlab='Position on Chr1 (kb)', ylab='read', main='Input')
    legend('topright', c('non-target', 'intermediate', 'target'), 
        lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = c('black', 'red', 'green'), bty='n')
    dev.off()

}

apply_jahmm <- function(x, readDat, coordinates, xrange, resolution,  output_dir){
    chipProfile <- readDat %>% 
                    select(CHR, INPUT, contains(x)) 
    fit <- jahmm(chipProfile)$path

    prots <- colnames(chipProfile)
    reps <- prots[grepl('_REP[0-9]', prots)]
    for (i in reps){
        plotFit(select(chipProfile, CHR, INPUT, contains(i)),
                fit, 1:xrange, 
                paste0(output_dir, i, '_', resolution, '_plotFit.pdf'))
    }
    peaks <- cbind(coordinates, fit) %>%
                filter(fit == 2)
    write.table(peaks, 
            paste0(output_dir, x, '_', resolution,'_peaks.bed'),  
            row.names=F, quote=F, sep='\t')

}

lapply(proteins, apply_jahmm, readDat, coordinates, opt$xrange, opt$resolution, opt$output_dir)
