# Chromatin Conformation Prediction

## Pipeline

- ReadAligns `*.SRA -> *.bam` 
 
    - Aligning the ChIPseq reads with the ref.genome using bowtie2
    
    getBamFromSRA 
    Usage -c [no.core] -a [annot.file] -l [list.file] -i [input_dir] -o [output_dir]
    [no.core]: for parallel processing
    [annot.file]: conversion of sra number to file descriptor
    [input_dir] : list of .rsa files
    [output_dir}: storage directory
    Prepare in advance the annot.txt and also double check for the list of files created (step no.3 in this script)
    e.g:
    time getBamFromSRA -c 4 -g /home/ricky/Rlim/Biotools/Genomes/hg19/bowtie2/hg19 -a CebpB/Input/annot.txt -i CebpB/Input -o CebpB/Output 2> CebpB/log_CebpB.txt &

- PeakCalls `*.bam -> *.xls`
    
    - Calling ChIPseq Peaks using macs2
     
    macsPeakCalling
    genome for human: hs and for mouse: mm
    macsPeakCalling -c 4 -g hs -i Input/bam -o Output/macs 2> log/macsPeakCall.txt
    e.g:
    time macsPeakCalling -c 4 -g hs -i CebpB/Input -o CebpB/Output 2> CebpB/log_CebpB.txt &
    
- ComponentCalls

    - Fitting the ChIPseq Peak Counts with Mixture Models (GMM AND/OR NBM) using R Code `ComponentCalls.R`
    - The R analysis is at ComponentAnalysis/
    
- MotifCalls

    - Using Centdist for motif database discovery
    - Using ChIP-meme for DeNovo motif discovery
    - Using rsat-peak-motifs for DeNovo motif discovery
