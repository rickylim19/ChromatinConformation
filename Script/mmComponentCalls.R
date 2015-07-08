#!/usr/bin/env Rscript

############################################################################################
# Pipeline to call components using mixture models                                         #
# Currently only gaussian mixture model has been applied                                   #
#                                                                                          #
# $ Rscript Script/mmComponentCalls.R --model='GMM'\                                       #
# --ncomp=3 --oneComp=TRUE\                                                                #
# --input_dir='/home/ricky/Rlim/ChromatinConformation/Input/ComponentCalls/Jahmm/Test/'\   #
# --output_dir='/home/ricky/Rlim/ChromatinConformation/Output/ComponentCalls/Jahmm/Test/'\ #
#                                                                                          #
############################################################################################

options(options = 3)
options(warn = 2)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mixtools"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("plyr"))

assessMM <- function(fit_model, k, d=1){

    ######################################################################
    #                                                                    #
    # Return BIC to assess the simplicity of all models that fit         #
    #  The smaller the BIC, the better the fit                           #
    #                                                                    #
    # Input:                                                             #
    #  - fitted model from normalmixEM (mixtools)                        #
    #  - d: dimension of the normalmix model, default:1                  #
    # Output:                                                            #
    #  - bic                                                             #
    #                                                                    #
    # e.g:                                                               #
    #  > fit <- normalmixEM(Koeffler_BM_Cebpe$pileup, k=k, maxit=5000)   #
    #  > assessMM(fit, k)                                               #
    #                                                                    #
    ######################################################################
    
    # P: number of parameters
    # Parameters of a GMM: mu(mean), sigma(stdev), and lambda(mixing probability)
    # e.g for 2 component-GMM, P=5, i.e, mu1 + mu2 + sigma1 + sigma2 + lambda 
    # Formula for the number of free parameters, referenced to
    #  R.C. Veltkamp, et al. 2001 on page 182, section 8.4.1
    #   State-of-the-Art in Content-Based Image and Video Retrival, 

    P <- (k-1) + k*d + k*((d*(d+1))/2)
    L <- fit_model$loglik 
    n <- length(fit_model$x) # number of obs

    bic <- (-2*L) + (log(n)*P) 
    #aic <- (-2*L) + (2*P)
    #result <- list(bic=bic, aic=aic)
    return (bic)
}

sortComponent <- function(inputData){

    ############################################################
    # Sort the assigned components based on the median scores  #
    # note that the higher the component the higher the score  #
    #                                                          #
    # > testCase                                               #
    #   chrom start end score component                        #
    # 1  chr1     0 300     1         3                        #
    # 2  chr1   301 600    10         1                        #
    # 3  chr1   601 900     5         2                        #
    # > sortComponent(testCase)                                #
    #   chrom start end score component                        #
    # 1  chr1     0 300     1         1                        #
    # 2  chr1   301 600    10         3                        #
    # 3  chr1   601 900     5         2                        #
    #                                                          #
    ############################################################

    # get the sorted component assignment according to the median
    unsortedAssignment <- ddply(inputData, 'component', 
                                summarise, medscore=median(score)) 
    sortedAssignment <- unsortedAssignment[with(unsortedAssignment, 
                                                order(medscore)), ]
    sortedAssignment$sortedComp <- 1:nrow(sortedAssignment)
    
    # substitute the unsorted comps to the sorted ones
    sortedComp <- mapvalues(inputData$component, from=sortedAssignment$component, 
                                      to=sortedAssignment$sortedComp)
    result <- mutate(inputData, component = sortedComp)
    return(result)
}



apply_mm <- function(x, ncomp, input_dir, output_dir, model = 'GMM', oneComp=TRUE){
    filename <- removeExt(x)

    if (model == 'GMM'){
        # log-transformation was required to normalize (gaussian) the input 
        inputFile <- read.table(paste0(input_dir, x), header=T, sep='\t') %>%
                        mutate(score = log(score))
        # fit once with only n-comp 
        if (oneComp){
            fit <- normalmixEM(inputFile$score, k = ncomp, maxit = 5000)
            inputFile$component <- apply(fit$posterior, 1, which.max)

        # fit multiple times from 2 up to n:comp
        } else{
            fit_list <- lapply(2:ncomp, function(x){
                                        fit <- normalmixEM(inputFile$score, k = x, maxit = 5000)
                                        fit$assessment <- assessMM(fit, k = x)
                                        return (fit)})
            bic <- lapply(1: length(fit_list), function(x) {
                                                return(fit_list[[x]]$assessment)})
            # return the component when the decrease is significant
            # note 2 was addedd as the component starts with 2 
            best_component <- which.min(diff(unlist(bic))) + 1 
            inputFile$component <- apply(fit_list[[best_component]]$posterior, 1, which.max)
            ncomp <- best_component
        }
        #sort the component based on median score and reverse the log-score by exp
        inputFile <- sortComponent(inputFile) %>%
                        mutate(score = round(exp(score))) 

        # write for each n component assigned
        invisible(lapply(seq(ncomp), function(x) inputFile %>%
                            filter(component == x) %>%
                            write.table(paste0(output_dir, filename, '_component', x, '.bed'),
                                        quote=F, row.names=F, sep='\t', col.names=F)))

        # write all components assigned
        write.table(inputFile, paste0(output_dir, filename, '_components.bed'),
                    quote=F, row.names=F, sep='\t', col.names=F)
                        
    } else {return (NULL)}
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
    make_option("--core", type = "integer", default = 4,
                help = "set the number of cores for parallelization",
                metavar = "number"),
    make_option("--model", type = "character", 
                help = "mixture model", default= 'GMM',
                metavar = "character"),
    make_option("--ncomp", type = "integer", 
                help = "number of components tried", default= '4',
                metavar = "number"),
    make_option("--oneComp", type = "logical", 
                help = "Fit with only ncomp not trying to find best fit from 2:ncomp", default= TRUE,
                metavar = "logical"),
    make_option("--input_dir", action = "store", type = "character", 
                help = "input directory",
                metavar = "character"),
    make_option("--output_dir", action = "store", type = "character", 
                help = "output directory",
                metavar = "character")
    )



### call the command with apply in parallel)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(file.path(opt$output_dir), showWarnings = F, recursive = T)
#filePattern = paste0(opt$resolution, 'bin-.*peaks.bed')
files <- list.files(path = opt$input_dir, pattern = '*_peaks.bed') 
invisible(mclapply(files, apply_mm, model = opt$model, ncomp = opt$ncomp, 
                   oneComp = opt$oneComp, input_dir = opt$input_dir, 
                   output_dir = opt$output_dir, mc.cores = opt$core))
