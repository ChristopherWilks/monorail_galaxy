"""
Parameters:
- star: arguments to pass to STAR aligner
- unique_qual: minimum MAPQ needed to be counted in unique BW [default: 10]
"""

STEPS = ['align', 'sort', 'bamcount']

FILES = ['sjout.zst',
         'bamcount_nonref.csv.zst',
         'bamcount_auc.tsv',
         'bamcount_frag.tsv',
         'Chimeric.out.junction.zst',
         'all.exon_bw_count.zst', 'unique.exon_bw_count.zst',
         'manifest']

import subprocess
def run_command(cmd_args):
    cmd_args = ' '.join(cmd_args)
    try:
        subprocess.check_call(args=cmd_args, shell=True) 
    except subprocess.CalledProcessError as cpe:
        sys.stderr.write("error in run_command for command: %s\n" % cmd_args)
        raise cpe


import re
FASTQ_PATT=re.compile(r'([^_\.]+)_(\d+)\.((fastq)|fq)(\.gz)?$')
import os
def prep_for_galaxy_run():
    try:
        os.mkdir(config['temp'])
    except OSError as ose:
        pass
    fastqs = config['inputs'].split(',')
    m = FASTQ_PATT.search(fastqs[0])
    run_acc = m.group(1)
    study_acc = run_acc
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
        run_command(['zcat',f,'>',newf])
        #create fastq 0
        if i == 1:
            try:
                os.symlink(os.path.abspath(newf), '%s/%s_%s_%s_%s_%d.fastq' % (config['temp'], run_acc, study_acc, genome, method, 0))
            except FileExistsError as fee:
                pass
        i += 1
    #create fastq 2 if not paired
    if i == 2:
        open('%s/%s_%s_%s_%s_%d.fastq' % (config['temp'], run_acc, study_acc, genome, method, 2), "w").close() 
    #create expected file structure for annotated exon bed file & reference index
    ref = config['ref']
    config['ref'] = '.'
    os.makedirs('%s/%s' % (config['ref'], genome))
    os.symlink(ref, '%s/%s/star_idx' % (config['ref'], genome))
    os.makedirs('%s/%s/gtf' % (config['ref'], genome))
    os.symlink(config['exon_bed'], 'exons.tmp')
    os.symlink('../../exons.tmp', '%s/%s/gtf/%s' % (config['ref'], genome, config.get('bw_bed', 'exons.bed')))
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
        config['output'] + '/{quad}.sjout.zst',
        config['output'] + '/{quad}.Chimeric.out.junction.zst',
        config['output'] + '/{quad}.bamcount_nonref.csv.zst',
        config['output'] + '/{quad}.bamcount_auc.tsv',
        config['output'] + '/{quad}.bamcount_frag.tsv',
        config['output'] + '/{quad}.all.exon_bw_count.zst',
        config['output'] + '/{quad}.unique.exon_bw_count.zst',
    output:
        config['output'] + '/{quad}.manifest'
    params:
        quad=lambda wildcards: wildcards.quad
    run:
        with open(output[0], 'wt') as fh:
            for fn in FILES:
                fh.write(params.quad + "." + fn + '\n')

rule bamcount:
    input:
        bam=config['temp'] + '/{quad}-sorted.bam',
        bamidx=config['temp'] + '/{quad}-sorted.bam.bai',
        #exe='/bamcount/bamcount',
        exe='bamcount',
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
        {input.exe} {input.bam} \
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
