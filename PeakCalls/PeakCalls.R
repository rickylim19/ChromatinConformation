library(Rjahmm)
library(xtable)
library(gridExtra)
library(GenomicRanges)
library(data.table)
library(reshape2)
library(VennDiagram)

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
        col = as.factor(fit$path), type = 'h', 
        xlab='Position on Chr1 (kb)', ylab='read', main='Protein')
    legend('topright', c('non-target', 'intermediate', 'target'), 
        lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = c('black', 'red', 'green'), bty='n')

    # plot the input
    plot(chipProfile[xrange, 2], 
        col = as.factor(fit$path), type = 'h', 
        xlab='Position on Chr1 (kb)', ylab='read', main='Input')
    legend('topright', c('non-target', 'intermediate', 'target'), 
        lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col = c('black', 'red', 'green'), bty='n')
    dev.off()

}

plotTargetInput <- function(target_df, ranges, ylim, output_f, title_f){

    ##################################################################
    # Input:                                                         #
    # - target_df:                                                   #
    #   chr start end Input CebpE                                    #
    #   chr1 5    10   5    10                                       #
    # - ranges:e.                                                    #
    #   number of targets plotted, e.g: 1:1000 means 1000 targets    #
    # - output_f                                                     #
    #   output filename                                              #
    #                                                                #
    # > plotTargetInput(macs_jahmm, 1:500,                           #
    #           'KoefflerLab_BM_ChIPseq_CebpE_macs_targets.pdf')     #
    ################################################################## 
    

    target_df$targets <- 1: nrow(target_df) 
    target.long <- melt(target_df[ranges],
                        id.vars = c('chr', 'start', 'end', 'targets'))

    p <- ggplot(data = target.long,
        aes(x=targets, y=value, colour=variable)) +
        stat_smooth() + 
        coord_cartesian(ylim = ylim) + 
        ggtitle(title_f) + theme_bw()
    ggsave(paste0(work_dir, output_f), p)

}

