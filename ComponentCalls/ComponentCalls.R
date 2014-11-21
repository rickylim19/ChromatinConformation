library(mixtools)
library(ggplot2)
library(reshape2)
library(xtable)
library(scales)
library(RColorBrewer)
library(ggbio)
library(rtracklayer)
library(BSgenome)
library(wordcloud)
library(broom) 
library(gridExtra)
library(plyr)

assessGMM <- function(fit_model, k, d=1){

    ######################################################################
    #                                                                    #
    # Return BIC and AIC to assess the simplicity of all models that fit #
    #  The smaller the BIC or AIC, the better the fit                    #
    #                                                                    #
    # Input:                                                             #
    #  - fitted model from normalmixEM (mixtools)                        #
    #  - d: dimension of the normalmix model, default:1                  #
    # Output:                                                            #
    #  - list of bic and aic                                             #
    #                                                                    #
    # e.g:                                                               #
    #  > fit <- normalmixEM(Koeffler_BM_Cebpe$pileup, k=k, maxit=5000)   #
    #  > assessGMM(fit, k)                                               #
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
    aic <- (-2*L) + (2*P)
    result <- list(bic=bic, aic=aic)
    return (result)
}

fitMMs <- function(X, max_k, model){

    ##########################################################################
    #                                                                        #
    # Returns the list of the fitted max_k-component-models                  #
    #  fit the data with mixture models starting from 2 components to        # 
    #  max_# k components                                                    #
    #  the list for each fit contains:                                       #
    #  - the normalmixEM object                                              #
    #  - the assessment criterion (BIC and AIC)                              #
    #  - the assignment component based on the most probable posterior.prob  #
    #                                                                        #
    # e.g:                                                                   #
    #  > GMMs_list <- fitMMs(Koeffler_BM_Cebpe$pileup, 7, model='GMM')       #
    #                                                                        #
    ##########################################################################

    n_component = length(seq(2,max_k))
    fit_list <- vector('list', length=n_component)

    for (k in seq(2,max_k)){
        if (model == 'GMM'){
            fit <- normalmixEM(X, k=k, maxit=5000)
            fit$assessment <- assessGMM(fit, k)
            fit$assignment <- apply(fit$posterior, 1, which.max)
        }
        if (model == 'NBM'){ 
            fit <- fit.nbmix(data=X, K=k)
        }
        fit_list[[k-1]] <- fit
    }
    return (fit_list)

} 

# Probability density function for a Gaussian mixture
  # Presumes the mixture object has the structure used by mixtools
dnormalmix <- function(x,mixture,log=FALSE) {
    lambda <- mixture$lambda
    k <- length(lambda)
    # Calculate share of likelihood for all data for one component
    like.component <- function(x,component) {
    lambda[component]*dnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
    }
    # Create array with likelihood shares from all components over all data
    likes <- sapply(1:k,like.component,x=x)
    # Add up contributions from components
    d <- rowSums(likes)
    if (log) {
    d <- log(d)
    }
    return(d)
}


visualizeFitMMs <- function(X, MMs_list, max_k, output_n, model='GMM', titleName){

    ##################################################################################
    #                                                                                #
    # Produces plots that fit the GMMs on the ChIP-seq Peaks                         #
    #  Two type of plots:                                                            #
    #  1. Fitting of each component in the mixture on the emperical distribution     #
    #  2. Fitting of all components in the mixture on the emperical distribution     #
    #     (Empirical vs Model)                                                       #
    #                                                                                #
    #                                                                                #
    # e.g:                                                                           #
    #  > visualizeFitMMs(X=Koeffler_BM_Cebpe$pileup, MMs_list=GMMs_list,             #
    #                     max_k=7, output_n='Koeffler_BM_Cebpe_GMM',                 #
    #                     model = 'GMM')                                             #
    #                                                                                #
    ##################################################################################

    #ifelse(logTransf, x <- X$pileup,
    #                  x <- seq(min(X$pileup)-1, max(X$pileup)+1, 1))
    x <- seq(min(X$pileup)-1, max(X$pileup)+1, 1)

    for (MM in MMs_list){
        #cbPalette <- brewer.pal(max_k, 'Spectral')
        #cbPalette <- brewer.pal(max_k, 'Set3')
        cbPalette <- colorRampPalette(brewer.pal(10,'Spectral'))(max_k)

        plot_table <- data.frame(x=x)

        if (model == 'GMM'){
        k <- length(MM$mu)
        for (i in seq(k)){
            plot_table[, paste0('fit',i)]<-dnorm(x, MM$mu[i], MM$sigma[i])*MM$lambda[i]
            }
        plot_table$fitSum <- apply(plot_table[, -1], 1, sum)
        }
        
        if (model == 'NBM'){
        k <- length(MM$prop)
        for (i in seq(k)){
            plot_table[, paste0('fit',i)]<-dnbinom(x, size = MM$size[i], 
                                                      prob = MM$prob[i])*MM$prop[i]
            }
        plot_table$fitSum <- apply(plot_table[, -1], 1, sum)
        }

        # fitting of each component in the mixture
        output_f <- paste0(output_n, '_comp', k, '.pdf')   
        p <- ggplot(X, aes(x=pileup))+
            #geom_histogram(aes(y=..density..), binwidth=1, fill='#969696') +
            geom_histogram(aes(y=..density..), fill='#969696') +
            geom_line(stat='density', alpha=0.6, size=1.5)
            for (i in seq(k)){
                p <- p + geom_area(data=plot_table, 
                                aes_string(x='x', y=paste0('fit', i)), 
                                #colour='#084594', size=1)
                                fill=cbPalette[i], size=1, alpha=.5)
      
            }

        p + ggtitle(paste0(k,'-Component ', model, '-', titleName)) 
        ggsave(file=output_f)

        # fitting of all components in the mixture
        # (comparing density estimates, empirical vs fitting)
        output_f2 <- paste0(output_n, '_comp', k, '_empVsModel', '.pdf')   
        p2 <- ggplot(X) + geom_line(stat='density', aes(x=pileup),  size=1)    
        p2 + geom_line(data=plot_table, aes(x = x, y = fitSum), 
                    alpha=0.8, linetype='dashed', colour='red', size=0.6) +
             ggtitle(paste0(k, '-Component ', model,'-', titleName))
        ggsave(file=output_f2)
    } 
}

plotFitMM <- function(X, MM_list, output_n, titleName, model='GMM', ymax){

    ###############################################################################
    #                                                                             #
    # Returns two types of plots in pdf format                                    #
    #  1. Fitting of each component in the mixture on the empirical distribution  #
    #  2. Fitting of all components in the mixture on the empirical distribution  #
    #     (Empirical vs. Mixture Model)                                           #
    #                                                                             #
    # >plotFitMM(Koeffler_BM_CebpE_log$pileup, GMMs_list_Koeffler_BM_CebpE_log,   #
    #            output_n='figs/Koeffler_BM_CebpE_GMM_ModelVisualization_log',    #
    #            titleName='Koeffler_BM_CebpE_log', model='GMM')                  #
    #                                                                             #
    ###############################################################################

    if (model=='GMM'){
        for (k in seq_along(MM_list)){

            # fit each component of the mixture model
            pdf(paste0(output_n, '_comp', k+1, '.pdf'), 
                useDingbats=FALSE)
            plot(MM_list[[k]], which=2, 
                main2=paste0(titleName, k+1, '.pdf'))
            lines(density(X, lty=2, lwd=2))
            dev.off()

            # empirical vs gaussian
            pdf(paste0(output_n, '_comp', k+1, '_empVsModel.pdf'), 
                useDingbats=FALSE)
            density_val<- density(X)
            plot(density_val,lty=1,ylim=c(0,max(density_val$y)), 
                main=paste("Comparison of density estimates\n",
                            "Empirical vs. Gaussian mixture\n",
                            titleName, '_comp', k+1, '.pdf'),bty='n')
            curve(dnormalmix(x,MM_list[[k]]),lty=2,add=TRUE)
            legend('topright', lty=c(1,2),c('= Empirical', '= Gaussian mixture'), bty='n')
            dev.off()
        }

    }
}

getModelAssessment <- function(MMs_list, max_k, output_n, model='GMM', titleName){

    ################################################################################
    #                                                                              #
    #  Produces table and plot of BIC and AIC for each component-MM                #
    #                                                                              #
    #  Input:                                                                      #
    #   - MMs_list: list of all the fitted Mixture model (MM) (output from fitMMs) #
    #   - output_n: the filename of the output                                     #
    #                                                                              #
    #  e,g:                                                                        #
    #  > getModelAssessment(GMMs_list, 7, 'Koeffler_BM_Cebpe_GMM_ModelAssessment', #
    #                        model='GMM')                                          #
    #   # the output of plots: ComponentCalls/Output/Plots/<output_n>              #
    #   # the output of tables: ComponentCallsOutput/Tables/<output_n>             #
    #                                                                              #
    ################################################################################

    n_component = length(seq(2,max_k))
    assessment_table <- data.frame(component=seq(2,max_k))
    bic <- c()
    aic <- c()
    #negLoglik <- c()


    for (MM in MMs_list){
        if (model == 'GMM'){
            bic <- append(bic, MM$assessment$bic)
            aic <- append(aic, MM$assessment$aic)
            #negLoglik <- append(negLoglik, -(MM$loglik))
        }

        if (model == 'NBM'){
            bic <- append(bic, MM$BIC)
            aic <- append(aic, MM$AIC)
        }
    }
    assessment_table$BIC <- bic
    assessment_table$AIC <- aic
    #assessment_table$negLoglik <- negLoglik

    # write table for BIC and AIC
    write.table(assessment_table, 
                file=paste0(output_n, '.txt'), 
                sep='\t', row.names=FALSE,col.names=colnames(assessment_table))

    # plot the AIC and BIC
    assessment_AICBIC <- assessment_table[, c('component', 'AIC', 'BIC')]
    assessment_AICBIC <- melt(assessment_AICBIC, id='component')

    pdf(paste0(output_n,'_AICBIC','.pdf'), 
        useDingbats=FALSE)
    p <- ggplot(data=assessment_AICBIC, 
                aes(x=component, y=value, colour=variable))+
            geom_line() + scale_colour_discrete('Criterion')+
            #coord_cartesian(ylim=c(488000, 505000)) + 
            ggtitle(paste0('Model Assessment:', model, 
                           '-', titleName))
    plot(p)
    dev.off()

}


assignComponentMMs <- function(X, MMs_list, output_n, model, sort.comp=FALSE){

    #########################################################################
    #                                                                       #
    # Produces tables (in BED format).                                      #
    #  sort.comp: sorting the components based on the range of pileup score #
    #   * sorting based on median values of each component score range      #
    #   * higher component corresponds to the higher score range            #
    #   e.g: comp.1 comp.2 comp.3                                           #
    #          20    50     30     (median score range)                     #
    #                                                                       #
    # The table looks like:                                                 #
    #  chr  start   end     pileup  comp # header                           #
    #   1   73832   73833   20  2                                           #
    #                                                                       #
    # > assignComponentMMs(Koeffler_BM_Cebpe, MMs_list,                     # 
    #                      'Koeffler_BM_Cebpe_ModelAssignment_NBM',         #
    #                      model = 'NBM', sort.comp=TRUE)                   #
    #                                                                       #
    #########################################################################             

    # col header: chr, start, end 

    options(scipen=999)
    bed <- data.frame(chr=X$chr, start=X$abs_summit, 
                      end=X$abs_summit+1, pileup=X$pileup) 
    #for (MM in MMs_list){
    for (i in 1: length(MMs_list)){
        if (model == 'GMM'){
            max_comp <- MMs_list[[i]]$assignment
            #k <- length(MM$mu)
            k <- length(MMs_list[[i]]$mu)
        }
        if (model == 'NBM'){
            max_comp <- class.nbmix(MMs_list[[i]], data = X$pileup) 
            k <- length(MMs_list[[i]]$prob)
        }

        bed$component <- max_comp

        if (sort.comp){

            # get the sorted component assignment according to the median
            unsortedAssignment <- ddply(bed, 'component', 
                                        summarise, medscore=median(pileup)) 
            sortedAssignment <- unsortedAssignment[with(unsortedAssignment, 
                                                   order(medscore)), ]
            sortedAssignment$sortedComp <- 1:nrow(sortedAssignment)
            
            # substitute the unsorted comps to the sorted ones
            bed$sortedComp <- mapvalues(bed$component, from=sortedAssignment$component, 
                                        to=sortedAssignment$sortedComp)
            result <- data.frame(bed$chr, bed$start, bed$end, bed$sortedComp, bed$pileup)
            compSorted <- '_compSorted'
            
        }
        else{
            result <- data.frame(bed$chr, bed$start, bed$end, bed$component, bed$pileup)
            compSorted <- '_comp'
        }

        colnames(result) <- c('chr', 'start', 'end', 'comp', 'pileup')

        # write per group file: e.g Koeffler_BM_Cebpe_ModelAssignment_group1_comp2.bed
        for (n in seq(1,k)){
            componentCluster <- result[result$comp==n, ]
            write.table(componentCluster, 
                        file=paste0( 
                        output_n,'_group', n, compSorted, k,".bed"), 
                        col.names=FALSE, row.names=FALSE, 
                        sep="\t", quote=FALSE)
        }

        write.table(result, 
                file=paste0(output_n, compSorted, k, '.bed'),
                col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
    }
}

getCompFreqTable <- function(X, k, n, model){

    ##########################################################
    #                                                        #
    # Returns a frequency table for each assigned component  #
    #  k -> k components                                     #
    #  n -> n times (repeat the EM fit algorithm n times)    #
    # > getCompFreqTable(Koeffler_BM_Cebpe$pileup, 15, 5,    #
    #                    model='GMM')                        #
    #                                                        #
    ##########################################################

    compFreq <- matrix(ncol=n, nrow=k)
    for (i in 1:n){
        if (model == 'GMM'){
            fit <- normalmixEM(X, k=k, maxit=5000)
            comp <- sort(table(apply(fit$posterior, 1, which.max)))
        }
        if (model == 'NBM'){
            fit <- fit.nbmix(data=X, K=k)
            comp <- sort(table(class.nbmix(fit, X)))
        }
        compFreq[,i] <- comp 

    }
    colnames(compFreq) <- sapply(seq(n), function(x) paste0('Rep.', x))
    rownames(compFreq) <- sapply(seq(k), function(x) paste0('Comp.', x))
    return(compFreq)
}

visualizeCompAssign <- function(bed_f, chrom, title_f, output_f){

    #################################################################################
    # Returns a pdf file that plots the peak score(pileup) on Y-axis                #
    #  and the genomic location on X-axis,                                          #
    #  each peak coloured according to the component assignment.                    #
    #  NOTE:                                                                        #
    #   - Only the first 50 peaks that are plotted here                             #
    #   - the range of the genomic location is according                            #
    #         to a certain chromosome, such as chr14                                #
    #                                                                               #
    # Input:                                                                        #
    #  - bed_f: bed file                                                            #
    #     The format of the bed file: chr start end comp score                      #
    #  - chrom: chromosome of the genomic location on X-axis                        #
    #  - title_f: title of the plot according to the filename                       #
    #  - output_f: output filename                                                  #
    #                                                                               #
    # > visualizeCompAssign('Koeffler_BM_Cebpe_ModelAssignment_comp3.bed',          #
    #                        chrom='chr14',                                         #
    #                        title_f='Koeffler_BM_Cebpe_ModelAssignment_comp3.bed', #
    #                        output_f='Koeffler_BM_Cebpe_ModelAssignment_comp3.pdf' #
    #                        )                                                      #
    #                                                                               #
    #################################################################################

    range_data <- import(bed_f) 
    chrData <- range_data[seqnames(range_data)==chrom]

    # get only the first 50 peaks
    chrData <- chrData[1:100]
    p <- autoplot(chrData, xlab=chrom, title=title_f, 
                  geom='bar', aes(col=name) , alpha=.5) + 
                  scale_colour_discrete(name='Component') + 
                  ggtitle(title_f)
    ggsave(file=output_f, width=15)
}

visualizeBedPeaks <- function(bed_f, chrom, title_f, output_f, max_peak){

    ####################################################################################################
    # Returns a pdf file that plots the score (Y-axis) and the genomic location (X-axis)               #
    #                                                                                                  #
    # Input:                                                                                           #
    #  - bed_f: #chr start end name score                                                              #
    #                                                                                                  #
    #     chr1    6214409 6216409 CEBPB   48                                                           #
    #     chr1    6466192 6468192 CTCF    86                                                           #
    #     chr1    6466192 6468192 CEBPB   86                                                           #
    #     chr1    7088912 7090912 CEBPB   31                                                           #
    #     chr1    7108615 7110615 CEBPB   29                                                           #
    # - chrom: coordinates chromosome you wish to visualize                                            #
    # - title_f: title of the plot                                                                     #
    # - output_f: the name of the output, (without the extension, as automatically in pdf)             #
    #   *Note that we appended _chunk in the output filename*                                          #
    # - max_peak: maximum number of peaks in each chunk                                                #
    #                                                                                                  #
    #                                                                                                  #
    # >visualizeBedPeaks(bed_f = paste0(work_dir,                                                      #
    #                   'Output/ChipMeme_Koeffler_BM_CebpE_GMM_ModelAssignment_log_compSorted4_2kb/',  #
    #                   'Motif/CEBPB_CTCF_Motif.bed'),                                                 #
    #                chrom='chr14', title_f = 'Koeffler_BM_CebpE_GMM_ModelAssignment_log_CEBPB_CTCF',  #
    #                output_f = 'figs/Koeffler_BM_CebpE_GMM_ModelAssignment_log_CEBPBCTCF',            #
    #                max_peak=50)                                                                      #
    #                                                                                                  #
    ####################################################################################################

    range_data <- import(bed_f)
    chrData <- range_data[seqnames(range_data)==chrom]
    x <- seq_along(chrData)
    chrData_chunks <- split(chrData, ceiling(x/max_peak))
    for (i in 1:length(chrData_chunks)){
        chrData_chunk <- chrData_chunks[[i]]
        p <- autoplot(chrData_chunk, xlab=chrom, title=title_f, 
                        geom='bar', aes(col=name) , alpha=.3) + 
                        ggtitle(title_f)
        ggsave(file=paste0(output_f, '_chunk', i,'_', chrom, '.pdf'), width=15)
    }
}


createPDFTable <- function(tabData, title_tab, output_f){
    
    ####################################################
    # Generate PDF from data.frame                     #
    #                                                  #
    # > createPDFTable(compFreq3, '3-Component GMM',   #
    #                  'Output/Tables/compFreq3.pdf')  # 
    #                                                  #
    ####################################################

    #pdf(output_f, useDingbats=FALSE, paper=paper, width=100)
    pdf(output_f, useDingbats=FALSE)
    #grid.newpage()
    tab <- tableGrob(tabData)
    h <- grobHeight(tab)
    w <- grobWidth(tab)
    #title <- textGrob(title_tab, 
    #                  y=unit(0.5,"npc") + 0.5*h, 
    #                  vjust=0, gp=gpar(fontsize=20))

    #gt <- gTree(children=gList(tab, title))
    gt <- gTree(children=gList(tab))
    grid.draw(gt)
    dev.off()
}

boxplot_compAssignment <- function(bed_f, title_f, output_f){

    ###################################################################################
    # Returns a boxplot of scores by components                                       #
    #                                                                                 #
    # Input Bed file with the format as follows:                                      #
    #  chr start end comp score e.g:                                                  #
    #  chr1 4391326 4391327    1    19                                                #
    #                                                                                 #
    #                                                                                 #
    # >boxplot_compAssignment('ComponentCalls/Output/Tables/                          #
    #                          Koeffler_BM_Cebpe_ModelAssignment_comp3.bed',          #
    #                          title_f='Distribution of Scores in 3-Component GMM',   #
    #                          'Output/Plots/Koeffler_BM_Cebpe_Boxplot_comp3.pdf')    #
    #                                                                                 #
    ###################################################################################

    compData <- read.table(bed_f)
    colnames(compData) <- c('chr', 'start', 'end', 'comp', 'score')
    ggplot(compData, aes(factor(comp), score, fill=factor(comp)))+
        geom_boxplot() + 
        scale_fill_discrete(name='Components') + 
        ggtitle(title_f)
    ggsave(file=output_f)
}


## Mixture NB ##

# function to fit mixture of negative binomial 
# (Poisson with Gamma overdispersion) distribution
# K = number of components
fit.nbmix <- function(data,K) {
   lambda <- rep(0,K-1)
   prob.init   <- seq(0.2,0.8,len=K)
   lsize.init <- log(rep(round(mean(data)),K))
   beta.init <- log(prob.init/(1-prob.init))
   par.init  <- c(lambda,beta.init,lsize.init)
    
   # fit model
   iter <- 1
   cat('Iteration ',iter,'\n',sep="")
   out <- optim(p=par.init,fn=nlogl.nbmix,data=data,K=K,method='Nelder-Mead')
   if(out$conv!=0) {
      out <- optim(p=out$p,fn=nlogl.nbmix,data=data,K=K,method='Nelder-Mead')
   }
   while(out$conv!=0) {
     iter <- iter+1
     cat('Iteration ',iter,'\n',sep="")
     out <- optim(p=out$p,fn=nlogl.nbmix,data=data,K=K,method='L-BFGS-B')
   }
   # collate results
   lambda = c(0,out$par[1:(K-1)])
   prop = exp(lambda)/sum(exp(lambda))

   prob = exp(out$par[K:(2*K-1)])/(1+exp(out$par[K:(2*K-1)]))
   size = exp(out$par[(2*K):(3*K-1)])

   AIC = 2*out$v + 2*length(out$p)
   BIC = 2*out$v + length(out$p)*log(length(data))

   result = list(prop=prop,prob=prob,size=size,AIC=AIC,BIC=BIC)

   result
}

#negative log-likelihood function for mixture of NB
nlogl.nbmix <- function(p,data,K) {
    
    lambda = c(0,p[1:(K-1)])
    prop = exp(lambda)/sum(exp(lambda))
    prob = exp(p[K:(2*K-1)])/(1+exp(p[K:(2*K-1)]))
    size = exp(p[(2*K):(3*K-1)])

    lik = t(apply(as.matrix(data),1,dnbinom,prob=prob,size=size)) %*% as.matrix(prop)
    #lik = ifelse(lik==0,1e-200,lik)
    logl = sum(log(lik))
    nlogl= - logl
nlogl
}

#class prediction 
class.nbmix <- function(fit,data) {
   K = length(fit$prop)
   prop = fit$prop
   prob = fit$prob
   size = fit$size

   mar.prob = t(apply(as.matrix(data),1,dnbinom,
                      prob=prob,size=size)) %*% as.matrix(diag(prop))
   mar.prob = mar.prob/apply(mar.prob,1,sum)

   max.index <- which(mar.prob==apply(mar.prob,1,max),arr.ind=TRUE)
   class <- rep(0,length(data))
   class[max.index[,1]] <- max.index[,2]

   class
}




## Local Clustering for direct and indirect bindings ##

bedClusterPeaks <- function(functionString='bedtools cluster', distance=100, bed){
    
    ###################################################################################
    #                                                                                 #
    # Calls the bedtools cluster to cluster peaks within the distance                 #
    #  Default distance is 100b                                                       #
    #                                                                                 #
    # > comp3_bedCluster <- read.table('Output/Tables/                                #
    #                                   Koeffler_BM_Cebpe_ModelAssignment_comp3.bed') #
    #                                                                                 #
    ###################################################################################

    options(scipen=99)
    bed.file = tempfile()
    out = tempfile()

    bed <- read.table(bed)
    write.table(bed, file=bed.file, quote=F, sep='\t', col.names=F, row.names=F)

    command=paste(functionString, '-d', distance, '-i', bed.file, '>', out, sep= ' ')
    cat(command, '\n')
    try(system(command))

    result <- read.table(out, header=F)
    unlink(bed.file); unlink(out)
    return (result)
}

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




