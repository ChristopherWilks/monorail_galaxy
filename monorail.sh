set -o pipefail -o errexit 

path2=$PATH

#e.g. /data7/references/ath10/star_idx
ref=$1
#e.g. 4
threads=$2
#e.g. 10
min_uniq_qual=$3
#e.g. ./ath10/gtf/exons.bed
annotation=$4
#e.g. ./tmp2
tmp=$5
#e.g. tmp/SRR8505407_SRP182756_ath10_sra_1.fastq.gz
reads1=$6
#e.g. tmp/SRR8505407_SRP182756_ath10_sra_2.fastq.gz, just leave off if not paired
reads2=$7

#align
STAR --runMode alignReads --runThreadN $threads --genomeDir ${ref} --readFilesIn $reads1 $reads2 --readFilesCommand zcat --twopassMode None --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMreadID Number --readFilesType Fastx --outTmpDir ./star_tmp --outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH MD --chimOutType Junctions --chimOutJunctionFormat 1 --chimSegmentReadGapMax 3 --chimJunctionOverhangMin 12 --chimSegmentMin 12 > star.log 2>&1

sout="Aligned.out.bam"

#sort
sorted_bam="sorted.bam"
samtools sort -T ./samtools_temp -@ $threads -m 64M -o $sorted_bam $sout > stools.log 2>&1
#index
samtools index -@ $threads $sorted_bam >> stools.log 2>&1

#get coverage summaries
bamcount $sorted_bam --threads $threads --coverage --no-head --require-mdz --min-unique-qual $min_uniq_qual --frag-dist bc --bigwig bc --annotation $annotation bc --auc bc --alts bc > bc.log 2>&1 
