set -o pipefail -o errexit 
#get this script's path
#p=`perl -e '@f=split("/","'$0'"); pop(@f); print "".join("/",@f)."\n";'`

#/data7/miniconda2/envs/_galaxy_18.09/bin:/data7/miniconda2/condabin:/home/cwilks/.nimble/bin:/data/pypy3-v6.0.0-linux64/bin:/home/cwilks/zstd_dir:/data/bedtools2/bin:/usr/local/cuda-8.0/bin:/home/cwilks/Python3/bin:/data2/sqlite-autoconf-3110000.fts5:/data/java/bin:/data/ant/bin:/home/cwilks/anaconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/cwilks/samtools-1.2:/home/cwilks/samtools-1.2/htslib-1.2.1:/data/WiggleTools/bin:/data/kent_tools:/data2/liftover:/home/cwilks/bin

path2=$PATH
#starp=`/bin/which STAR`
#stp=`/bin/which samtools`

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


#racc=`perl -e '$r1="'${reads1}'"; @f=split(/\//,$r1); $r1=pop(@f); ($n1)=split(/[\._]/,$r1); print "$n1\n";'`

#cd ${tmp}
#align
STAR --runMode alignReads --runThreadN $threads --genomeDir ${ref} --readFilesIn $reads1 $reads2 --readFilesCommand zcat --twopassMode None --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMreadID Number --readFilesType Fastx --outTmpDir ./star_tmp --outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH MD --chimOutType Junctions --chimOutJunctionFormat 1 --chimSegmentReadGapMax 3 --chimJunctionOverhangMin 12 --chimSegmentMin 12 > star.log 2>&1

sout="Aligned.out.bam"

#sort
sorted_bam="sorted.bam"
samtools sort -T ./samtools_temp -@ $threads -m 64M -o $sorted_bam $sout > stools.log 2>&1
#index
samtools index -@ $threads $sorted_bam >> stools.log 2>&1

#get coverage summaries
#bamcount $sorted_bam --threads $threads --coverage --no-head --require-mdz --min-unique-qual $min_uniq_qual --frag-dist bc --bigwig bc --annotation $annotation bc --auc bc --lats bc > bc.log 2>&1 

#echo $PATH > ${sorted_bam}.bai
#echo $starp >> ${sorted_bam}.bai
#echo $stp >> ${sorted_bam}.bai
