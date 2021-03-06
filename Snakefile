"""
Parameters:
- fastq_dump_args: arguments to pass to fastq dumping tool
- fastq_dump_retries: number of retry attempts before dying
- fastq_dump_tool: name of tool to run, e.g. 'fastq-dump' or 'fasterq-dump'
- star: arguments to pass to STAR aligner
- salmon_args: arguments to pass to salmon quant
- unique_qual: minimum MAPQ needed to be counted in unique BW [default: 10]
- bw_bed: name of BED file to use with bwtool
- max_unalign: maximum number of unaligned reads to save per run accession
"""

STEPS = ['download', 'fastq_check', 'align', 'sort',
         'bamcount',
         'salmon',
         'align_unmapped',
         'extract_jx']

FILES = ['sjout.zst', 'fastq_check.tsv.zst',
         'unmapped.bam', 'unmapped.fastq.zst',
         'bamcount_nonref.csv.zst',
         'bamcount_auc.tsv',
         'bamcount_frag.tsv',
         'Chimeric.out.junction.zst',
         'all.exon_bw_count.zst', 'unique.exon_bw_count.zst',
         'all.bw.zst', 'unique.bw.zst',
         'salmon.tsv.zst',
         'jx_bed.zst',
         'manifest'] + list(map(lambda x: x + '.log', STEPS))

#get rid of the pesky warnings about current paths
#and make sure we're not assuming any wrong relative paths
import os
config['output']=os.path.abspath(config['output'])
config['temp']=os.path.abspath(config['temp'])
config['exon_bed']=os.path.abspath(config['exon_bed'])
config['ref']=os.path.abspath(config['ref'])

import subprocess
def run_command(cmd_args):
    cmd_args = ' '.join(cmd_args)
    try:
        subprocess.check_call(args=cmd_args, shell=True) 
    except subprocess.CalledProcessError as cpe:
        sys.stderr.write("error in run_command for command: %s\n" % cmd_args)
        raise cpe

import re
STEPS_FILES_FILTER=re.compile(r'(unmapped)|(download)|(salmon)|(extract_jx)|(jx_bed)|(manifest)')
def remove_steps_files():
    #modify STEP and FILES
    #so we don't run download or unmapped steps
    global FILES
    global STEPS
    for i in reversed(range(0,len(STEPS))):
        f = STEPS[i]
        if STEPS_FILES_FILTER.search(f):
            del(STEPS[i])
    for i in reversed(range(0,len(FILES))):
        f = FILES[i]
        if STEPS_FILES_FILTER.search(f):
            del(FILES[i])


FASTQ_PATT=re.compile(r'([^_\.]+)(_(\d+))?\.((fastq)|fq)(\.gz)?$')
def prep_for_galaxy_run():
    remove_steps_files()
    try:
        os.mkdir(config['temp'])
    except OSError as ose:
        pass
    fastqs = config['inputs'].split(',')
    m = FASTQ_PATT.search(fastqs[0])
    run_acc = 'sample'
    if m is not None:
        run_acc = m.group(1)
    study_acc = 'study'
    if 'study' in config:
        study_acc = config['study']
    genome = 'hg38'
    if 'genome' in config:
        genome = config['genome']
    method = 'sra'
    # SRR,SRP,genome
    # e.g. SRR1557855,SRP045778,ce10
    #create links which will be used in the align step
    i = 1
    for f in fastqs:
        newf = '%s/%s_%s_%s_%s_%d.fastq' % (config['temp'], run_acc, study_acc, genome, method, i)
        if 'compressed' in config:
            run_command(['zcat',f,'>',newf])
        else:
            os.symlink(os.path.abspath(f), newf)
        #create fastq 0
        if i == 1:
            os.symlink(os.path.abspath(newf), '%s/%s_%s_%s_%s_%d.fastq' % (config['temp'], run_acc, study_acc, genome, method, 0))
        i += 1
    #create fastq 2 if not paired
    if i == 2:
        open('%s/%s_%s_%s_%s_%d.fastq' % (config['temp'], run_acc, study_acc, genome, method, 2), "w").close() 
    return([run_acc, study_acc, genome, method])


def get_accessions(wildcards):
    """
    Grouping of SRRs with the same SRP could happen here
    """
    #if running under galaxy where the user will input the 
    #FASTQs, study, and genome directly
    if 'inputs' in config:
        (run_acc, study_acc, genome, method) = prep_for_galaxy_run()
        for ext in FILES:
            yield os.path.join(config['output'], '%s_%s_%s_%s.%s' % (run_acc, study_acc, genome, method, ext))
        #here to get the make_galaxy_links rule to fire
        yield os.path.join(config['output'], '%s_%s_%s_%s.done' % (run_acc, study_acc, genome, method))
    #not running under galaxy, takes a list of accessions
    else:
        for fn in config['input'].split():
            with open(fn, 'r') as fh:
                for ln in fh:
                    if ln.count(',') < 2:
                        continue
                    toks = ln.rstrip().split(',')
                    assert 3 <= len(toks) <= 4
                    method = 'sra'
                    if len(toks) == 4:
                        method = toks[3]
                    # SRR,SRP,genome
                    # e.g. SRR1557855,SRP045778,ce10
                    for ext in FILES:
                        yield os.path.join(config['output'], '%s_%s_%s_%s.%s' % (toks[0], toks[1], toks[2], method, ext))

rule all:
    input:
        get_accessions

rule make_manifest:
    input:
        config['output'] + '/{quad}.salmon.tsv.zst',
        config['output'] + '/{quad}.sjout.zst',
        config['output'] + '/{quad}.jx_bed.zst',
        config['output'] + '/{quad}.Chimeric.out.junction.zst',
        config['output'] + '/{quad}.unmapped.bam',
        config['output'] + '/{quad}.unmapped.fastq.zst',
        config['output'] + '/{quad}.bamcount_nonref.csv.zst',
        config['output'] + '/{quad}.bamcount_auc.tsv',
        config['output'] + '/{quad}.bamcount_frag.tsv',
        config['output'] + '/{quad}.fastq_check.tsv.zst',
        config['output'] + '/{quad}.all.exon_bw_count.zst',
        config['output'] + '/{quad}.unique.exon_bw_count.zst',
        config['output'] + '/{quad}.all.bw.zst',
        config['output'] + '/{quad}.unique.bw.zst',
        config['output'] + '/{quad}.align.log',
        config['output'] + '/{quad}.extract_jx.log',
        config['output'] + '/{quad}.bamcount.log',
        config['output'] + '/{quad}.align_unmapped.log',
        config['output'] + '/{quad}.download.log',
        config['output'] + '/{quad}.fastq_check.log',
        config['output'] + '/{quad}.sort.log',
        config['output'] + '/{quad}.salmon.log'
    output:
        config['output'] + '/{quad}.manifest'
    params:
        quad=lambda wildcards: wildcards.quad
    run:
        with open(output[0], 'wt') as fh:
            for fn in FILES:
                fh.write(params.quad + "." + fn + '\n')

def galaxy_link_files(op):
    remove_steps_files()
    a = [op + '/' + f for f in FILES]
    return a

rule make_galaxy_links:
    input: 
        config['output'] + '/{quad}.sjout.zst',
        config['output'] + '/{quad}.fastq_check.tsv.zst',
        config['output'] + '/{quad}.bamcount_nonref.csv.zst',
        config['output'] + '/{quad}.bamcount_auc.tsv',
        config['output'] + '/{quad}.bamcount_frag.tsv',
        config['output'] + '/{quad}.Chimeric.out.junction.zst',
        config['output'] + '/{quad}.all.exon_bw_count.zst',
        config['output'] + '/{quad}.unique.exon_bw_count.zst',
        config['output'] + '/{quad}.all.bw.zst',
        config['output'] + '/{quad}.unique.bw.zst',
        config['output'] + '/{quad}.fastq_check.log',
        config['output'] + '/{quad}.align.log',
        config['output'] + '/{quad}.sort.log',
        config['output'] + '/{quad}.bamcount.log',
    output:
        config['output'] + '/{quad}.done'
    params:
        quad=lambda wildcards: wildcards.quad,
        out=galaxy_link_files(config['output'])
    run:
        for (i,fn) in enumerate(input):
            os.symlink(os.path.abspath(fn), params.out[i])
        os.symlink(os.path.abspath(input[-1]), output[0])


rule bamcount:
    input:
        bam=config['temp'] + '/{quad}-sorted.bam',
        bamidx=config['temp'] + '/{quad}-sorted.bam.bai',
        bed=lambda wildcards: '%s/%s/gtf/%s' % (config['ref'], wildcards.quad.split('_')[2], config.get('bw_bed', 'exons.bed'))
    output:
        nonref=config['output'] + '/{quad}.bamcount_nonref.csv.zst',
        auc=config['output'] + '/{quad}.bamcount_auc.tsv',
        frag=config['output'] + '/{quad}.bamcount_frag.tsv',
        all_bw=config['output'] + '/{quad}.all.bw.zst',
        unique_bw=config['output'] + '/{quad}.unique.bw.zst',
        all_bw_count=config['output'] + '/{quad}.all.exon_bw_count.zst',
        unique_bw_count=config['output'] + '/{quad}.unique.exon_bw_count.zst'
    log:
        config['output'] + '/{quad}.bamcount.log'
    params:
        srr=lambda wildcards: wildcards.quad.split('_')[0],
        uniq_qual=config.get('unique_qual', 10)
    threads: 4
    shell:
        """
        TMP={config[temp]}/{params.srr}_bamcount
        bamcount {input.bam} \
            --threads {threads} \
            --coverage \
            --no-head \
            --require-mdz \
            --min-unique-qual {params.uniq_qual} \
            --frag-dist ${{TMP}} \
            --bigwig ${{TMP}} \
            --annotation {input.bed} ${{TMP}} \
            --auc ${{TMP}} \
            --alts ${{TMP}} \
            2>&1 | tee -a {log}

        #
        # --alts
        #

        (time zstd ${{TMP}}.alts.tsv -o {output.nonref}) 2>&1 | tee -a {log}
        size=$(wc -c < {output.nonref})
        echo "COUNT_NonrefSize ${{size}}"
        rm -f ${{TMP}}.alts.tsv

        #
        # --auc
        #
        mv ${{TMP}}.auc.tsv {output.auc}
        size=$(wc -c < {output.auc})
        echo "COUNT_AucSize ${{size}}"
        rm -f ${{TMP}}.auc.tsv

        #
        # --frag-dist
        #
        mv ${{TMP}}.frags.tsv {output.frag}
        size=$(wc -c < {output.frag})
        echo "COUNT_FragDistSize ${{size}}"
        rm -f ${{TMP}}.frags.tsv

        #
        # --bigwig
        #

        (time zstd ${{TMP}}.all.bw -o {output.all_bw}) 2>&1 | tee -a {log}
        size=$(wc -c < {output.all_bw})
        echo "COUNT_BwSize ${{size}}"
        rm -f ${{TMP}}.all.bw

        (time zstd ${{TMP}}.unique.bw -o {output.unique_bw}) 2>&1 | tee -a {log}
        size=$(wc -c < {output.unique_bw})
        echo "COUNT_BwSize ${{size}}"
        rm -f ${{TMP}}.unique.bw

        #
        # --annotation
        #

        (time zstd ${{TMP}}.all.tsv -o {output.all_bw_count}) 2>&1 | tee -a {log}
        size=$(wc -c < {output.all_bw_count})
        echo "COUNT_BwQuantSize ${{size}}"
        rm -f ${{TMP}}.all.tsv

        (time zstd ${{TMP}}.unique.tsv -o {output.unique_bw_count}) 2>&1 | tee -a {log}
        size=$(wc -c < {output.unique_bw_count})
        echo "COUNT_BwQuantSize ${{size}}"
        rm -f ${{TMP}}.unique.tsv

        # Check that all temporaries were properly purged
        set +o pipefail ; num_files=$(ls -d ${{TMP}}* 2>/dev/null | wc -l)
        if (( $num_files > 0 )) ; then
            echo "Failed to purge files (ignore . and ..): $(ls -ad ${{TMP}}*)"
            exit 1
        fi

        echo "COUNT_BamcountComplete 1"
        """

rule bw_zstd:
    input:
        config['temp'] + '/{prefix}.bw'
    output:
        config['output'] + '/{prefix}.bw.zst'
    shell:
        """
        zstd {input} -o {output}
        size=$(wc -c < {output})
        echo "COUNT_BwBytes ${{size}}"
        echo "COUNT_BwZstdComplete 1"
        """

rule sort:
    input:
        config['temp'] + '/{quad}.bam'
    wildcard_constraints:
        quad="[^-]+"
    output:
        bam=temp(config['temp'] + '/{quad}-sorted.bam'),
        bai=temp(config['temp'] + '/{quad}-sorted.bam.bai')
    log:
        config['output'] + '/{quad}.sort.log'
    params:
        srr=lambda wildcards: wildcards.quad.split('_')[0]
    threads: 8
    shell:
        """
        TMP="{config[temp]}/sort_temp.{params.srr}"
        mkdir -p ${{TMP}}
        time samtools sort \
            -T ${{TMP}}/samtools_temp \
            -@ {threads} \
            -m 64M \
            -o {output.bam} {input} 2>&1 | tee -a {log}
        rm -rf ${{TMP}}
        size=$(wc -c < {output.bam})
        echo "COUNT_SortedBAMBytes ${{size}}"

        time samtools index -@ {threads} {output.bam} 2>&1 | tee -a {log}
        echo "COUNT_SortComplete 1"
        """

rule salmon:
    input:
        reads0=config['temp'] + '/{quad}_0.fastq',
        reads1=config['temp'] + '/{quad}_1.fastq',
        reads2=config['temp'] + '/{quad}_2.fastq',
        index1=lambda wildcards: '%s/%s/salmon_index/hash.bin' % (config['ref'], wildcards.quad.split('_')[2]),
        index2=lambda wildcards: '%s/%s/salmon_index/sa.bin' % (config['ref'], wildcards.quad.split('_')[2])
    output:
        config['output'] + '/{quad}.salmon.tsv.zst'
    log:
        config['output'] + '/{quad}.salmon.log'
    params:
        index_base=lambda wildcards: '%s/%s/salmon_index' % (config['ref'], wildcards.quad.split('_')[2]),
        salmon_args=config.get('salmon_args', '')
    threads: 8
    shell:
        """
        READ_FILES="-r {input.reads0}"
        if [[ -s {input.reads2} ]] ; then
            READ_FILES="-1 {input.reads1} -2 {input.reads2}"
        fi
        if set -o pipefail && time salmon quant \
            --libType U \
            --quiet \
            --validateMappings \
            -i {params.index_base} \
            -p {threads} \
            {params.salmon_args} \
            ${{READ_FILES}} \
            --output salmon_quant \
            --minAssignedFrags 1 \
            2>&1 | tee -a {log}
        then
            echo "COUNT_SalmonSuccess 1"

            time zstd salmon_quant/quant.sf -o {output} 2>&1 | tee -a {log}
            size=$(wc -c < {output})
            echo "COUNT_SalmonQuantBytes ${{size}}"
        else
            touch {output}
            echo "COUNT_SalmonFailure 1"
        fi
                
        rm -rf salmon_quant
        echo "COUNT_SalmonComplete 1"
        """

rule align_unmapped:
    input:
        unmapped1=config['temp'] + '/{quad}_1.unmappedfastq',
        unmapped2=config['temp'] + '/{quad}_2.unmappedfastq',
        index=lambda wildcards: '%s/%s/unmapped_hisat2_idx/genome.1.ht2' % (config['ref'], wildcards.quad.split('_')[2])
    output:
        bam=config['output'] + '/{quad}.unmapped.bam',
        sample=config['output'] + '/{quad}.unmapped.fastq.zst'
    log:
        config['output'] + '/{quad}.align_unmapped.log'
    params:
        index_base=lambda wildcards: '%s/%s/unmapped_hisat2_idx/genome' % (config['ref'], wildcards.quad.split('_')[2]),
        srr=lambda wildcards: wildcards.quad.split('_')[0],
        hisat2_params=config.get('hisat2', ''),
        max_unalign=config.get('max_unalign', 100000)
    threads: 16
    shell:
        """
        TMP="{config[temp]}/align_unmapped_temp.{params.srr}"
        READ_FILES="-1 {input.unmapped1} -2 {input.unmapped2}"
        if [[ ! -s {input.unmapped2} ]] ; then
            READ_FILES="-U {input.unmapped1}"
        fi
        time hisat2 \
            $READ_FILES \
            -t --mm \
            -x {params.index_base} \
            --threads {threads} \
            {params.hisat2_params} \
            --un ${{TMP}}.un \
            --un-conc ${{TMP}}.un_conc \
            -S ${{TMP}}.sam \
            2>&1 | tee -a {log}

        #
        # Make BAM file out of aligned reads
        #
        time samtools view \
            -b -F 4 \
            -o {output.bam} \
            ${{TMP}}.sam \
            2>&1 | tee -a {log}
        rm -f ${{TMP}}.sam

        #
        # Save a subset of the doubly unaligned reads
        #
        if [[ ! -s {input.unmapped2} ]] ; then
            # unpaired
            fastq-sample \
                -n {params.max_unalign} \
                -o ${{TMP}}.samp \
                ${{TMP}}.un
            rm -f ${{TMP}}.un
        else
            # paired-end
            fastq-sample \
                -n {params.max_unalign} \
                -o ${{TMP}}.samp \
                ${{TMP}}.1.un_conc ${{TMP}}.2.un_conc
            rm -f ${{TMP}}.1.un_conc ${{TMP}}.2.un_conc
            
            # interleave the two output fastqs into single file
            paste ${{TMP}}.samp.1.fastq ${{TMP}}.samp.2.fastq \
                | paste - - - - \
                | awk -v OFS="\n" -v FS="\t" '{{print($1,$3,$5,$7,$2,$4,$6,$8)}}' \
                > ${{TMP}}.samp.fastq
            rm -f ${{TMP}}.samp.1.fastq ${{TMP}}.samp.2.fastq
        fi

        test -f ${{TMP}}.samp.fastq
        time zstd ${{TMP}}.samp.fastq -o {output.sample} 2>&1 | tee -a {log}
        size=$(wc -c < {output.sample})
        echo "COUNT_UnmappedSampleBytes ${{size}}"
        rm -f ${{TMP}}.samp.fastq
    
        size=$(wc -c < {output.bam})
        echo "COUNT_UnmappedBamBytes ${{size}}"

        echo "COUNT_AlignUnmappedComplete 1"
        """

rule extract_jx:
    input:
        bam=config['temp'] + '/{quad}-sorted.bam',
        bamidx=config['temp'] + '/{quad}-sorted.bam.bai',
        fa=lambda wildcards: '%s/%s/fasta/genome.fa' % (config['ref'], wildcards.quad.split('_')[2]),
        gtf=lambda wildcards: '%s/%s/gtf/genes.gtf' % (config['ref'], wildcards.quad.split('_')[2])
    output:
        config['output'] + '/{quad}.jx_bed.zst'
    params:
        srr=lambda wildcards: wildcards.quad.split('_')[0]
    log:
        config['output'] + '/{quad}.extract_jx.log'
    shell:
        """
        TMP="{config[temp]}/extract_jx.{params.srr}"
        time regtools junctions extract \
            -i 20 -a 1 \
            -o ${{TMP}}.jx_tmp \
            {input.bam} 2>&1 | tee -a {log}
        time zstd ${{TMP}}.jx_tmp -o {output} 2>&1 | tee -a {log}
        rm -f ${{TMP}}.jx_tmp

        size=$(wc -c < {output})
        echo "COUNT_ExtractJxBytes ${{size}}"

        echo "COUNT_ExtractJxComplete 1"
        """

rule align:
    input:
        reads0=config['temp'] + '/{quad}_0.fastq',
        reads1=config['temp'] + '/{quad}_1.fastq',
        reads2=config['temp'] + '/{quad}_2.fastq',
        index1=lambda wildcards: '%s/%s/star_idx/SAindex' % (config['ref'], wildcards.quad.split('_')[2]),
        index2=lambda wildcards: '%s/%s/star_idx/SA' % (config['ref'], wildcards.quad.split('_')[2])
    wildcard_constraints:
        quad="[^-]+"
    output:
        bam=temp(config['temp'] + '/{quad}.bam'),
        jxs=config['output'] + '/{quad}.sjout.zst',
        chimeric=config['output'] + '/{quad}.Chimeric.out.junction.zst',
        unmapped1=config['temp'] + '/{quad}_1.unmappedfastq',
        unmapped2=config['temp'] + '/{quad}_2.unmappedfastq'
    log:
        config['output'] + '/{quad}.align.log'
    params:
        index_base=lambda wildcards: '%s/%s/star_idx' % (config['ref'], wildcards.quad.split('_')[2]),
        srr=lambda wildcards: wildcards.quad.split('_')[0],
        star_params="%s %s" % (config.get('star', ''), '--genomeLoad LoadAndRemove' if 'inputs' not in config else '')
    threads: 16
    shell:
        """
        READ_FILES="{input.reads0}"
        if [[ -s {input.reads2} ]] ; then
            READ_FILES="{input.reads1} {input.reads2}"
        fi
        TMP="{config[temp]}/align_temp.{params.srr}"
        rm -rf ${{TMP}}
        time STAR \
            {params.star_params} \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.index_base} \
            --readFilesIn ${{READ_FILES}} \
            --twopassMode None \
            --outTmpDir ${{TMP}} \
            --outReadsUnmapped Fastx \
            --outMultimapperOrder Random \
            --outSAMreadID Number \
            --outSAMtype BAM Unsorted \
            --outSAMmode NoQS \
            --outSAMattributes NH MD \
            --chimOutType Junctions \
            --chimOutJunctionFormat 1 \
            --chimSegmentReadGapMax 3 \
            --chimJunctionOverhangMin 12 \
            --chimSegmentMin 12 2>&1 | tee -a {log}
   
        # Full set of output files:
        #
        # Aligned.out.bam
        # Chimeric.out.junction
        # Log.final.out
        # Log.out
        # Log.progress.out
        # SJ.out.tab
        # Unmapped.out.mate1
        # Unmapped.out.mate2 (if any reads were paired-end)

        #
        # Logs
        #
        rm -rf ${{TMP}}
        cat Log.out >> {log}
        cat Log.final.out >> {log}
        rm -f Log*.out

        #
        # Junctions
        #
        test -f SJ.out.tab
        time zstd SJ.out.tab -o {output.jxs} 2>&1 | tee -a {log}
        rm -f SJ.out.tab
        size=$(wc -c < {output.jxs})
        echo "COUNT_CompressedJxBytes ${{size}}"

        #
        # Chimerics
        #
        test -f Chimeric.out.junction
        test -s Chimeric.out.junction
        sort -k1,1 -n -k2,2 Chimeric.out.junction > Chimeric.out.junction.sorted
        time zstd Chimeric.out.junction.sorted -o {output.chimeric} 2>&1 | tee -a {log}
        rm -f Chimeric.out.junction Chimeric.out.junction.sorted
        size=$(wc -c < {output.chimeric})
        echo "COUNT_ChimericBytes ${{size}}"

        #
        # Unmapped
        #
        touch {output.unmapped2}
        test -f Unmapped.out.mate1
        mv Unmapped.out.mate1 {output.unmapped1}
        if [[ -f Unmapped.out.mate2 ]] ; then
            mv Unmapped.out.mate2 {output.unmapped2}
        fi

        #
        # Alignments
        #
        size=$(wc -c < Aligned.out.bam)
        echo "COUNT_BAMBytes ${{size}}"
        mv Aligned.out.bam {output.bam}

        echo "COUNT_AlignComplete 1"
        """

rule fastq_check:
    input:
        reads0=config['temp'] + '/{quad}_0.fastq',
        reads1=config['temp'] + '/{quad}_1.fastq',
        reads2=config['temp'] + '/{quad}_2.fastq'
    output:
        config['output'] + '/{quad}.fastq_check.tsv.zst'
    log:
        config['output'] + '/{quad}.fastq_check.log'
    params:
        srr=lambda wildcards: wildcards.quad.split('_')[0]
    shell:
        """
        TMP="{config[temp]}/fastq_check-{params.srr}.tsv"
        touch ${{TMP}}
        if [[ -s {input.reads0} ]] ; then
            time seqtk fqchk -q0 {input.reads0} >>${{TMP}} 2>>{log}
        fi
        if [[ -s {input.reads1} ]] ; then
            time seqtk fqchk -q0 {input.reads1} >>${{TMP}} 2>>{log}
        fi
        if [[ -s {input.reads2} ]] ; then
            time seqtk fqchk -q0 {input.reads2} >>${{TMP}} 2>>{log}
        fi
        time zstd ${{TMP}} -o {output} 2>&1 | tee -a {log}
        size=$(wc -c < {output})
        echo "COUNT_FastqCheckBytes ${{size}}"

        echo "COUNT_FastqCheckComplete 1"
        """

rule download:
    output:
        temp(config['temp'] + '/{quad}_0.fastq'),
        temp(config['temp'] + '/{quad}_1.fastq'),
        temp(config['temp'] + '/{quad}_2.fastq')
    log:
        config['output'] + '/{quad}.download.log'
    params:
        srr=lambda wildcards: wildcards.quad.split('_')[0],
        method=lambda wildcards: wildcards.quad.split('_')[3],
        fd_params=config.get('fastq_dump_args', ''),
        retries=config.get('fastq_dump_retries', '5')
    threads: 4
    shell:
        """
                set -x
        SUCCESS=0
        TIMEOUT=10
        PARAMS=""
        TMP="{config[temp]}/dl-{params.srr}"
        ! test -d ${{TMP}}
        
        if [[ {params.method} == "sra" ]] ; then
            USE_FASTERQ=1
            for i in {{ 1..{params.retries} }} ; do
                if time prefetch -t fasp -O ${{TMP}} -L info {params.srr} 2>&1 >> {log} ; then
                    SUCCESS=1
                    break
                else
                    echo "COUNT_SraRetries 1"
                    TIMEOUT=$((${{TIMEOUT}} * 2))
                    sleep ${{TIMEOUT}}
                fi
            done
            if (( $SUCCESS == 0 )) ; then
                echo "COUNT_SraFailures 1"
                exit 1
            fi
            test -f ${{TMP}}/*.sra
            size=$(cat ${{TMP}}/*.sra | wc -c)
            echo "COUNT_SraBytesDownloaded ${{size}}"
            if (( ${{USE_FASTERQ}} == 1 )) ; then
                time fasterq-dump ${{TMP}}/*.sra \
                    -e {threads} \
                    -t ${{TMP}} \
                    -L info \
                    --split-files \
                    --skip-technical \
                    -o {params.srr}.fastq \
                    2>&1 >> {log}
                test -f {params.srr}_2.fastq || mv {params.srr}.fastq {params.srr}_0.fastq
            else
                time fastq-dump ${{TMP}}/*.sra \
                    -L info \
                    --split-files \
                    --skip-technical \
                    -O . \
                    2>&1 >> {log}
                test -f {params.srr}_2.fastq || mv {params.srr}_1.fastq {params.srr}_0.fastq
            fi
            rm -rf ${{TMP}}
        elif [[ {params.method} == "gdc" ]] ; then
            TOKEN=~/gdc/latest.txt
            if [[ ! -f ${{TOKEN}} ]] ; then
                echo "ERROR: no GDC token file found at ${{TOKEN}}"
                exit 1
            fi
            for i in {{ 1..{params.retries} }} ; do
                if time gdc-client download \
                    -t ${{TOKEN}} \
                    --log-file {{TMP}}/log.txt \
                    -d {{TEMP}} \
                    {params.srr} 2>&1 >> {log}
                then
                    SUCCESS=1
                    break
                else
                    echo "COUNT_GdcRetries 1"
                    TIMEOUT=$((${{TIMEOUT}} * 2))
                    sleep ${{TIMEOUT}}
                fi
            done
            if (( $SUCCESS == 0 )) ; then
                echo "COUNT_GdcFailures 1"
                exit 1
            fi
            test -d {{TEMP}}/{params.srr}
            test -f {{TEMP}}/{params.srr}/*.tar.gz
            
            echo "=== gdc-client log.txt begin ===" >> {log}
            cat {{TMP}}/log.txt >> {log}
            echo "=== gdc-client log.txt end===" >> {log}
            
            size=$(cat {{TEMP}}/{params.srr}/*.tar.gz | wc -c)
            echo "COUNT_GdcBytesDownloaded ${{size}}"

            tar zxvf {{TEMP}}/{params.srr}/*.tar.gz
            rm -rf {{TEMP}}
            
            num_1s=$(ls -1 *_1.fastq | wc -l)
            num_2s=$(ls -1 *_2.fastq | wc -l)
            if (( ${{num_1s}} == 0 )) ; then
                echo "ERROR: No _1.fastq files output"
                exit 1
            fi
            if (( ${{num_1s}} > 1 )) ; then
                echo "ERROR: More than one _1.fastq file found"
                exit 1
            fi
            if (( ${{num_2s}} == 0 )) ; then
                # unpaired
                mv *_1.fastq {params.srr}_0.fastq
            else
                # paired-end
                mv *_1.fastq {params.srr}_1.fastq
                mv *_2.fastq {params.srr}_2.fastq
            fi
        fi
        
        # Next chunk expects the FASTQ files to exist in the current directory
        # named {{params.srr}}_{{0,1,2}}.fastq
        size=0
        for i in {{0..2}} ; do
            fn={params.srr}_${{i}}.fastq
            if [[ -f ${{fn}} ]] ; then
                echo "COUNT_FilesDownloaded 1"
            else
                touch ${{fn}}
            fi
            size=$((${{size}} + $(wc -c < ${{fn}})))
            mv ${{fn}} {config[temp]}/{wildcards.quad}_${{i}}.fastq
        done
        echo "COUNT_BytesDownloaded ${{size}}"
        echo "COUNT_DownloadComplete 1"
        """
