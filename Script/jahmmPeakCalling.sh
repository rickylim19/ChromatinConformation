#! /bin/bash
#
# peakCallJahmm.sh
# Copyright (C) 2015 ricky <ricky@tblab-Precision-T1650>
#
#

###############################################
# Peak calling pipeline using jahmm           #
#                                             #
# input: *.bam                                #
#                                             #
# output: *_peak.bed                          #
#                                             #
# Scripts used:                               #
# -bedtools bamtobed                          #
# -./binitBed.py                              #
# -Rscript jahmmPeakCalls.R                   #
#                                             #
###############################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
    echo "./jahmmPeakCalling.sh -c 4 -g hg19 -r 300 -i Input/ -o Output/ 2> log/jahmmPeakCall.txt"
    echo "sampleFile.csv contains the filename but with *.bed extension instead of *.bam !!!" 
exit 1
fi

while getopts ":c:g:r:s:i:o:" opt; do
    case $opt in
        c)
            echo "# core: $OPTARG" >&2
            core="$OPTARG"
            ;;
        g)
            echo "# genome: $OPTARG" >&2
            genome="$OPTARG"
            ;;
        r)
            echo "# resolution: $OPTARG" >&2
            resolution="$OPTARG"
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

# check if the bam files are sorted and there is the sorted input file (for control during peak calling)
#[[ $(ls -A $inputDir*.sorted.bam) && $(ls -A $inputDir*Input.sorted.bam) ]] && echo "contains sorted bam files" >&2 || exit 1
[[ $(ls -A $inputDir*.sorted.bam) && $(ls -A $inputDir*Input*.sorted.bam) ]] && echo "contains sorted bam files" >&2 || exit 1

# convert bam to bed
ls $inputDir*.sorted.bam | parallel -j $core "bedtools bamtobed -i {} > {.}.bed" && echo "Bam to bed conversion: OK" >&2;

# bin the bed tools, output in 300bin-*.bed
ls $inputDir*.sorted.bed | parallel -j $core binitBed.py -b $resolution -l $genome -F 'bed' -n $core -od $inputDir {}  && echo "Binning Bed in $resolution bp: OK" >&2;

# peak calling
Rscript Script/jahmmPeakCalls.R --core=$core --input_dir=$inputDir --output_dir=$outputDir --resolution=$resolution && echo "Jahmm Peak calling in $resolution bp: OK" >&2 && rm $inputDir*.bed;
wait
# removal of unwanted files
#rm $inputDir*.bed
