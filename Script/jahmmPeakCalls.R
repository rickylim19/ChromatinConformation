#!/usr/bin/env Rscript

options(options = 3)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Rjahmm"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("parallel"))

legendFit <- function(path){
    ########################################################################
    #                                                                      #
    # >legendFit(path)                                                     #
    # input: path object returned by jahmm$path                            #
    #                                                                      #
    # Return the legend based on the number of states returned by jahmm    #
    #                                                                      #
    ########################################################################
    numStates = length(levels(as.factor(path)))
    if (numStates == 3) {
        return(legend('topright', c('non-target', 'intermediate', 'target'), 
            lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = 1:3, bty='n'))
    } else if ( (numStates == 2) & (all(levels(as.factor(path)) == c("0", "2"))) ){
        return(legend('topright', c('non-target', 'target'), 
            lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = 1:2, bty='n'))
    
    } else if (numStates == 2 & all(levels(as.factor(path)) == c("0", "1"))){
        return(legend('topright', c('non-target', 'intermediate'), 
            lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = 1:2, bty='n'))
    
    } else {
        return(legend('topright', 'non-target', 
            lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = 1, bty='n'))
    }
}

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
    # fit$path object returned by jahmm                                                            #
    # - xrange                                                                                #
    # range of the chromosome positions (range of the rows of the dataframe from chipProfile) #
    # - output_dir                                                                            #
    # output to return the pdf coupled with its name                                          #
    #                                                                                         #
    ###########################################################################################

    pdf(paste0(output_dir))
    par(mfrow=c(2,1))
    # plot the protein
    plot(chipProfile[1:xrange, 3], 
        col = as.factor(fit), type = 'h', 
        xlab='Position on Chr1 (kb)', ylab='read', main='Protein')
    legendFit(fit)

    # plot the input
    plot(chipProfile[1:xrange, 2], 
        col = as.factor(fit), type = 'h', 
        xlab='Position on Chr1 (kb)', ylab='read', main='Input')
    legendFit(fit)
    dev.off()

}

apply_jahmm <- function(x, input, coordinates, xrange, resolution,  input_dir, output_dir){
    filename <- file_path_sans_ext(x) %>%
                    basename()
    chip <- read.table(paste0(input_dir, x), sep="\t", header = F) 
    colnames(chip) <- c("chrom", "start", "end", filename)
    inputJahmm <- cbind(chip, input) %>%
                    select(chrom, input, contains(filename))
    fit <- jahmm(inputJahmm)$path
    #fit <- invisible(jahmm(inputJahmm))$path

    plotFit(inputJahmm, fit, xrange, 
            paste0(output_dir, filename, '_plotFit.pdf'))

    score <- chip[, filename]
    peaks <- cbind(coordinates, score, fit) %>%
                filter(fit == 2 & score > 0)
    write.table(peaks[, 1:4], 
            paste0(output_dir, filename, '_peaks.bed'),  
            row.names=F, quote=F, sep='\t')

}

### main arguments ###
option_list <- list(
    make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                help = "Print extra output [default]"),
    make_option(c("-q", "--quietly"), action = "store_false",
                dest = "verbose", help = "Print little output"),
    make_option(c("-s", "--seed"), type = "integer", default = 12345,
                help = "set the seed number for reproducibility",
                metavar = "number"),
    make_option(c("-c", "--core"), type = "integer", default = 4,
                help = "set the number of cores for parallelization",
                metavar = "number"),
    make_option("--xrange", type = "integer", default = 100000,
                help = "set the number of rows plotted in the plotfit",
                metavar = "number"),
    make_option("--resolution", type = "character", default = "300",
                help = "bin resolution",
                metavar = "character"),
    make_option("--input_dir", action = "store", type = "character", 
                help = "input directory",
                metavar = "character"),
    make_option("--output_dir", action = "store", type = "character", 
                help = "output directory",
                metavar = "character")
    )

### call the command with apply in parallel)
opt <- parse_args(OptionParser(option_list = option_list))
files <- list.files(path = opt$input_dir, pattern = opt$resolution) 
inputFile <- files[grepl("Input", ignore.case=T, files)]
samples <- files[files != inputFile]
input <- read.table(paste0(opt$input_dir, inputFile), sep="\t", header=F)
colnames(input) <- c('chrom', 'start', 'end', 'input')
coordinates <- input[, 1:3]
input <- input$input
xrange <- opt$xrange
resolution <- opt$resolution
input_dir <- opt$input_dir
output_dir <- opt$output_dir

invisible(mclapply(samples, apply_jahmm, input, coordinates, xrange, resolution, 
                    input_dir, output_dir, mc.cores = opt$core))
