#! /bin/bash
#
# Copyright (C) 2015 ricky <ricky@tblab-Precision-T1650>
#
#

###############################################
# Peak calling pipeline using macs            #
#                                             #
# Input: *.bam                                #
#                                             #
# Output: *_peak.bed                          #
#                                             #
###############################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
    echo "./macsPeakCalling.sh -c 4 -g hs -i Input/ -o Output/ 2> log/jahmmPeakCall.txt"
exit 1
fi

while getopts ":c:g:s:i:o:" opt; do
    case $opt in
        c)
            echo "# core: $OPTARG" >&2
            core="$OPTARG"
            ;;
        g)
            echo "# genome: $OPTARG" >&2
            genome="$OPTARG"
            ;;
        i)
            echo "# inputDir: $OPTARG" >&2
            inputDir="$OPTARG"
            ;;
        o)
            echo "# outputDir: $OPTARG" >&2
            outputDir="$OPTARG"
            ;;
        \?)
            echo "# invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "# option -$OPTARG requires an argument" >&2
            exit 1
            ;;
    esac
done

# check if the directory exists
if [ ! -d $outputDir ]; then
    mkdir -p $outputDir
fi

# check if the bam files are sorted 
[[ $(ls -A $inputDir*.sorted.bam) ]] && echo "contains sorted bam files" >&2 || exit 1

# call macs peak depending if there is a control or not
control=( $(find $inputDir -maxdepth 1 -iname '*Input.sorted.bam') )
if [ -z "$control" ];
    then 
        echo "no Control for peak calling" >&2
        find $inputDir ! -name "*Input.sorted.bam" -type f | parallel -j $core "macs2 callpeak -g $genome -p 0.001 --call-summits -t {} -n {.}"   
    else
        find $inputDir ! -name "*Input.sorted.bam" -type f | parallel -j $core "macs2 callpeak -g $genome -p 0.001 --call-summits -t {} -c $control -n {.}"   
fi

# move the results into the output directory
find $inputDir ! -name "*.bam" -type f | parallel -j $core "mv {} $outputDir"

# add header and get only columns of chr, start, end, score
n=0
maxjobs=$core
for f in `ls $outputDir*_summits.bed`;
do
    fileName="${f%.*}"
    outputFile="$fileName"_peaks.bed""
    (awk 'BEGIN{print "chrom\tstart\tend\tscore"}1 {print $1"\t"$2"\t"$3"\t"$5}' $f > $outputFile) &
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
