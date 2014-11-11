
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
library(plyr)
library(xtable)
library(ggplot2)
library(reshape)
work_dir='/home/ricky/Rlim/ChromatinConformation/DimerCebpMotif/'
#work_dir = '/Users/RickyLim/Desktop/CebpE/'

Barplot_Chunk <- function(dataTable, N, vars, output_n, ylabel, TERM=TRUE){

    #############################################################################
    # Returns N barplots that chunks the dataTable                              #
    #                                                                           #
    # dataTable looks like (tab-delimited)                                      #
    #                                                                           #
    # Term Function Percent.Homo Fold.Homo Bonferroni.Homo ..                   #
    #  1    binding  1.5          1.44     1.075                                #
    #                                                                           #
    #                                                                           #
    #                                                                           #
    # > Barplot_Chunk(dataTable=homohetero_GO, N=12,                            #    
    #                    vars=c('Term', 'Benjamini.Homo', 'Benjamini.Hetero'),  #    
    #                    ylabel='- Log(P Value)',                               #    
    #                    output_n= 'figs/homohetero_BenjaminiGO')               #    
    #                                                                           #
    #############################################################################

    chunks <- split(dataTable,                                         
                    (as.numeric(rownames(dataTable))-1) %/% N )        
    for (i in 1: length(chunks)){                                      
        chunk.m <- melt(chunks[[i]][, vars], id.vars=1)                
                                                                        
        if (TERM){
            p <- ggplot(chunk.m, aes(x=Term, y=value, fill=variable)) +    
                    geom_bar(stat='identity', position='dodge') +          
                    theme(axis.text.x = element_text(angle=90, hjust=1))+  
                    ylab(ylabel)+xlab('')                                  
            }
        else{
        #p <- ggplot(chunk.m, aes(x=Term, y=value, fill=variable)) +    
            p <- ggplot(chunk.m, aes(x=Function, y=value, fill=variable)) +    
                    geom_bar(stat='identity', position='dodge') +          
                    theme(axis.text.x = element_text(angle=90, hjust=1))+  
                    ylab(ylabel)+xlab('')                                  
        }
        p                                                              
        ggsave(paste0(output_n,'_', i, '.pdf'))                        
    }                                                                  
}                                                                      



## ----GODavid, cache=TRUE-------------------------------------------------
homodimer_GO <- read.table(paste0(work_dir,
                    'CebpE/Output/ChipMeme/Motif/',
                    'homo_GOChart_Percent_Fold_Bonferroni_Benjamini.txt'),
                    header=FALSE,sep='\t')

colnames(homodimer_GO) <- c('Term', 'Function',
                            'Percent', 'Fold',
                            'Bonferroni','Benjamini')
homodimer_GO[, 'Benjamini'] <- (-log(homodimer_GO$Benjamini))
homodimer_GO[, 'Bonferroni'] <- (-log(homodimer_GO$Bonferroni))
head(homodimer_GO)

heterodimer_GO <- read.table(paste0(work_dir,
                    'CebpE/Output/ChipMeme/Motif/',
                    'hetero_GOChart_Percent_Fold_Bonferroni_Benjamini.txt'),
                    header=FALSE,sep='\t')
colnames(heterodimer_GO) <- c('Term', 'Function',
                            'Percent', 'Fold',
                            'Bonferroni','Benjamini')
heterodimer_GO[, 'Benjamini'] <- (-log(heterodimer_GO$Benjamini))
heterodimer_GO[, 'Bonferroni'] <- (-log(heterodimer_GO$Bonferroni))

# combine the term
homohetero_GO <- merge(homodimer_GO, heterodimer_GO, by.x='Term', by.y='Term')
homohetero_GO <- subset(homohetero_GO, select=-Function.y)
colnames(homohetero_GO) <- c('Term', 'Function',
                             'Percent.Homo', 'Fold.Homo',
                             'Bonferroni.Homo','Benjamini.Homo',
                             'Percent.Hetero', 'Fold.Hetero',
                             'Bonferroni.Hetero','Benjamini.Hetero')

head(homohetero_GO)    

# Benjamini 
Barplot_Chunk(dataTable=homohetero_GO, N=12,
              vars=c('Function', 'Benjamini.Homo', 'Benjamini.Hetero'), 
              ylabel='- Log(P Value)',
              output_n= 'figs/homohetero_BenjaminiGO', TERM=FALSE)

# Percent
Barplot_Chunk(dataTable=homohetero_GO, N=12,
              vars=c('Function', 'Percent.Homo', 'Percent.Hetero'), 
              ylabel='Percent',
              output_n= 'figs/homohetero_PercentGO', TERM=FALSE)



homodimerOnly_GO <- subset(homodimer_GO, !homodimer_GO$Term %in% heterodimer_GO$Term)
heterodimerOnly_GO <- subset(heterodimer_GO, !heterodimer_GO$Term %in% homodimer_GO$Term)
rownames(heterodimerOnly_GO) <- 1:nrow(heterodimerOnly_GO)
rownames(homodimerOnly_GO) <- 1:nrow(homodimerOnly_GO)
homodimerOnly_GO$Function <- substr(homodimerOnly_GO$Function, 1, 50) 
head(homodimerOnly_GO)
dim(heterodimerOnly_GO)


# homodimerOnly Benjamini
Barplot_Chunk(dataTable=homodimerOnly_GO, N=20,
              vars=c('Function', 'Benjamini'), 
              ylabel='- Log(P Value)',
              output_n='figs/homodimerOnly_BenjaminiGO', TERM=FALSE)

# homodimerOnly Percent
Barplot_Chunk(dataTable=homodimerOnly_GO, N=20,
              vars=c('Function', 'Percent'), 
              ylabel='Percent',
              output_n='figs/homodimerOnly_PercentGO', TERM=FALSE)

# heterodimerOnly Benjamini 
Barplot_Chunk(dataTable=heterodimerOnly_GO, N=20,
              vars=c('Function', 'Benjamini'), 
              ylabel='- Log(P Value)',
              output_n='figs/heteroOnly_BenjaminiGO', TERM=FALSE)

# heterodimerOnly Percent
Barplot_Chunk(dataTable=heterodimerOnly_GO, N=20,
              vars=c('Function', 'Percent'), 
              ylabel='Percent',
              output_n='figs/heteroOnly_PercentGO', TERM=FALSE)


annot_GO_homohetero <- subset(homohetero_GO, select=c('Term', 'Function'))
annot_GO_homohetero_chunks <-  split(annot_GO_homohetero,
                    (as.numeric(rownames(annot_GO_homohetero))-1) %/% 35 )        
length(annot_GO_homohetero_chunks)

annot_GO_homo <- subset(homodimerOnly_GO, select=c('Term', 'Function'))
annot_GO_homo_chunks <-  split(annot_GO_homo,
                    (as.numeric(rownames(annot_GO_homo))-1) %/% 35 )        
length(annot_GO_homo_chunks)

annot_GO_hetero <- subset(heterodimerOnly_GO, select=c('Term', 'Function'))
annot_GO_hetero_chunks <-  split(annot_GO_hetero,
                    (as.numeric(rownames(annot_GO_hetero))-1) %/% 35 )        
length(annot_GO_hetero_chunks)



## ----GOannotHomoHetero1, echo=FALSE, results='asis'----------------------
print(xtable(annot_GO_homohetero_chunks[[1]], caption='GO Functional Annotation in Homo and Heterodimer (1)'))


## ----GOannotHomoHetero2, echo=FALSE, results='asis'----------------------
print(xtable(annot_GO_homohetero_chunks[[2]], caption='GO Functional Annotation in Homo and Heterodimer (2)'))


## ----GOannotHomo1, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[1]], caption='GO Functional Annotation in Homodimer (1)'))


## ----GOannotHomo2, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[2]], caption='GO Functional Annotation in Homodimer (2)'))


## ----GOannotHomo3, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[3]], caption='GO Functional Annotation in Homodimer (3)'))


## ----GOannotHomo4, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[4]], caption='GO Functional Annotation in Homodimer (4)'))


## ----GOannotHomo5, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[5]], caption='GO Functional Annotation in Homodimer (5)'))


## ----GOannotHomo6, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[6]], caption='GO Functional Annotation in Homodimer (6)'))


## ----GOannotHomo7, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[7]], caption='GO Functional Annotation in Homodimer (7)'))


## ----GOannotHomo8, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[8]], caption='GO Functional Annotation in Homodimer (8)'))


## ----GOannotHomo9, echo=FALSE, results='asis'----------------------------
print(xtable(annot_GO_homo_chunks[[9]], caption='GO Functional Annotation in Homodimer (9)'))


## ----GOannotHomo10, echo=FALSE, results='asis'---------------------------
print(xtable(annot_GO_homo_chunks[[10]], caption='GO Functional Annotation in Homodimer (10)'))


## ----GOannotHetero1, echo=FALSE, results='asis'--------------------------
print(xtable(annot_GO_hetero_chunks[[1]], caption='GO Functional Annotation in Heterodimer (1)'))


## ----GOannotHetero2, echo=FALSE, results='asis'--------------------------
print(xtable(annot_GO_hetero_chunks[[2]], caption='GO Functional Annotation in Heterodimer (2)'))


## ----GOannotHetero3, echo=FALSE, results='asis'--------------------------
print(xtable(annot_GO_hetero_chunks[[3]], caption='GO Functional Annotation in Heterodimer (3)'))


## ----GOPathway, cache=TRUE-----------------------------------------------
homodimer_path <- read.table(paste0(work_dir,
                    'CebpE/Output/ChipMeme/Motif/',
                    'homo_Pathway_Percent_Fold_Bonferroni_Benjamini.txt'),
                    header=FALSE,sep='\t')

colnames(homodimer_path) <- c('Term', 'Function',
                            'Percent', 'Fold',
                            'Bonferroni','Benjamini')
homodimer_path[, 'Benjamini'] <- (-log(homodimer_path$Benjamini))
homodimer_path[, 'Bonferroni'] <- (-log(homodimer_path$Bonferroni))
head(homodimer_path)

heterodimer_path <- read.table(paste0(work_dir,
                    'CebpE/Output/ChipMeme/Motif/',
                    'hetero_Pathway_Percent_Fold_Bonferroni_Benjamini.txt'),
                    header=FALSE,sep='\t')
colnames(heterodimer_path) <- c('Term', 'Function',
                            'Percent', 'Fold',
                            'Bonferroni','Benjamini')
heterodimer_path[, 'Benjamini'] <- (-log(heterodimer_path$Benjamini))
heterodimer_path[, 'Bonferroni'] <- (-log(heterodimer_path$Bonferroni))

# combine the term
homohetero_path <- merge(homodimer_path, heterodimer_path, by.x='Term', by.y='Term')
homohetero_path <- subset(homohetero_path, select=-Function.y)
head(homohetero_path)
colnames(homohetero_path) <- c('Term', 'Function',
                             'Percent.Homo', 'Fold.Homo',
                             'Bonferroni.Homo','Benjamini.Homo',
                             'Percent.Hetero', 'Fold.Hetero',
                             'Bonferroni.Hetero','Benjamini.Hetero')
dim(homohetero_path)
head(homohetero_path)
# Benjamini 
Barplot_Chunk(dataTable=homohetero_path, N=2,
              vars=c('Function', 'Benjamini.Homo', 'Benjamini.Hetero'), 
              ylabel='- Log(P Value)',
              output_n='figs/homohetero_BenjaminiPath', TERM=FALSE)
# Percent
Barplot_Chunk(dataTable=homohetero_path, N=2,
              vars=c('Function', 'Percent.Homo', 'Percent.Hetero'), 
              ylabel='Percent',
              output_n='figs/homohetero_PercentPath', TERM=FALSE)

# homodimerOnly or heterodimerOnly
homodimerOnly_path <- subset(homodimer_path, !homodimer_path$Term %in% heterodimer_path$Term)
heterodimerOnly_path <- subset(heterodimer_path, !heterodimer_path$Term %in% homodimer_path$Term)
rownames(heterodimerOnly_path) <- 1:nrow(heterodimerOnly_path)
rownames(homodimerOnly_path) <- 1:nrow(homodimerOnly_path)

# homodimerOnly Benjamini
Barplot_Chunk(dataTable=homodimerOnly_path, N=20,
              vars=c('Function', 'Benjamini'), 
              ylabel='- Log(P Value)',
              output_n='figs/homodimerOnly_BenjaminiPath', TERM=FALSE)

# homodimerOnly Percent
Barplot_Chunk(dataTable=homodimerOnly_path, N=20,
              vars=c('Function', 'Percent'), 
              ylabel='Percent',
              output_n='figs/homodimerOnly_PercentPath', TERM=FALSE)

# heterodimerOnly Benjamini 
Barplot_Chunk(dataTable=heterodimerOnly_path, N=7,
              vars=c('Function', 'Benjamini'), 
              ylabel='- Log(P Value)',
              output_n='figs/heteroOnly_BenjaminiPath', TERM=FALSE)

# heterodimerOnly Percent
Barplot_Chunk(dataTable=heterodimerOnly_path, N=7,
              vars=c('Function', 'Percent'), 
              ylabel='Percent',
              output_n='figs/heteroOnly_PercentPath', TERM=FALSE)

annot_Path <- subset(homohetero_path, select=c('Term', 'Function'))
annot_Path_homoOnly <- subset(homodimerOnly_path, select=c('Term', 'Function'))
annot_Path_heteroOnly <- subset(heterodimerOnly_path, select=c('Term', 'Function'))


## ----PathwayAnnotHomoHetero, echo=FALSE, results='asis'------------------
print(xtable(annot_Path, caption='GO Pathway Annotation (Homo-Heterodimer)'))


## ----PathwayAnnotHeteroOnly, echo=FALSE, results='asis'------------------
print(xtable(annot_Path_heteroOnly, caption='GO Pathway Annotation (Heterodimer Only)'))


## ----PathwayAnnotHomoOnly, echo=FALSE, results='asis'--------------------
print(xtable(annot_Path_homoOnly, caption='GO Pathway Annotation (Homodimer Only)'))


## ----InteractDavid, cache=TRUE-------------------------------------------
homodimer_Interact <- read.table(paste0(work_dir,
                    'CebpE/Output/ChipMeme/Motif/',
                    'homo_InteractChart_Percent_Fold_Bonferroni_Benjamini.txt'),
                    header=FALSE,sep='\t')

colnames(homodimer_Interact) <- c('Term', 'Percent', 'Fold', 'Bonferroni','Benjamini')
homodimer_Interact[, 'Bonferroni'] <- (-log(homodimer_Interact$Bonferroni))
homodimer_Interact[, 'Benjamini'] <- (-log(homodimer_Interact$Benjamini))
head(homodimer_Interact)

heterodimer_Interact <- read.table(paste0(work_dir,
                    'CebpE/Output/ChipMeme/Motif/',
                    'hetero_InteractChart_Percent_Fold_Bonferroni_Benjamini.txt'),
                    header=FALSE,sep='\t')
colnames(heterodimer_Interact) <- c('Term', 'Percent', 'Fold', 'Bonferroni','Benjamini')
heterodimer_Interact[, 'Bonferroni'] <- (-log(heterodimer_Interact$Bonferroni))
heterodimer_Interact[, 'Benjamini'] <- (-log(heterodimer_Interact$Benjamini))

# combine the term
homohetero_Interact <- merge(homodimer_Interact, heterodimer_Interact, by.x='Term', by.y='Term')
colnames(homohetero_Interact) <- c('Term', 
                             'Percent.Homo','Fold.Homo', 'Bonferroni.Homo','Benjamini.Homo',
                             'Percent.Hetero','Fold.Hetero', 'Bonferroni.Hetero','Benjamini.Hetero')

head(homohetero_Interact)    
dim(homohetero_Interact)

# Benjamini 
Barplot_Chunk(dataTable=homohetero_Interact, N=12,
              vars=c('Term', 'Benjamini.Homo', 'Benjamini.Hetero'), 
              ylabel='- Log(P Value)',
              output_n= 'figs/homohetero_BenjaminiInteract', TERM=TRUE)

# Percent
Barplot_Chunk(dataTable=homohetero_Interact, N=12,
              vars=c('Term', 'Percent.Homo', 'Percent.Hetero'), 
              ylabel='Percent',
              output_n= 'figs/homohetero_PercentInteract', TERM=TRUE)



homodimerOnly_Interact <- subset(homodimer_Interact, !homodimer_Interact$Term %in% heterodimer_Interact$Term)
heterodimerOnly_Interact <- subset(heterodimer_Interact, !heterodimer_Interact$Term %in% homodimer_Interact$Term)
rownames(heterodimerOnly_Interact) <- 1:nrow(heterodimerOnly_Interact)
rownames(homodimerOnly_Interact) <- 1:nrow(homodimerOnly_Interact)
head(homodimerOnly_Interact)


# homodimerOnly Benjamini
Barplot_Chunk(dataTable=homodimerOnly_Interact, N=21,
              vars=c('Term', 'Benjamini'), 
              ylabel='- Log(P Value)',
              output_n='figs/homodimerOnly_BenjaminiInteract', TERM=TRUE)

# homodimerOnly Percent
Barplot_Chunk(dataTable=homodimerOnly_Interact, N=21,
              vars=c('Term', 'Percent'), 
              ylabel='Percent',
              output_n='figs/homodimerOnly_PercentInteract', TERM=TRUE)

# heterodimerOnly Benjamini 
Barplot_Chunk(dataTable=heterodimerOnly_Interact, N=20,
              vars=c('Term', 'Benjamini'), 
              ylabel='- Log(P Value)',
              output_n='figs/heteroOnly_BenjaminiInteract', TERM=TRUE)

# heterodimerOnly Percent
Barplot_Chunk(dataTable=heterodimerOnly_Interact, N=20,
              vars=c('Term', 'Percent'), 
              ylabel='Percent',
              output_n='figs/heteroOnly_PercentInteract', TERM=TRUE)





## ------------------------------------------------------------------------
sessionInfo()


## ----knitIt, cache=TRUE, results='hide', message=FALSE, warning=FALSE----
library(knitr)
purl("dimerAnalysis.Rnw" ) # compile to tex
purl("dimerAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("dimerAnalysis.Rnw")


