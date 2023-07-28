import tempfile
import pandas as pd
from itertools import product

configfile: 'config/config.yaml'
MANIFEST = config.get('manifest', 'config/manifest.tab')
TMPDIR = tempfile.mkdtemp(prefix='flagger-snakemake.',suffix='.tmp')

manifest_df = pd.read_csv(MANIFEST, sep='\t')
manifest_df['tech'] = manifest_df['tech'].str.lower()
manifest_df = manifest_df.set_index(["sample", "tech"], drop=False)

flagger_container = 'docker://mobinasri/flagger:v0.2'

wildcard_constraints:
    sample='|'.join(manifest_df.index.get_level_values('sample')),
    tech='hifi|ont',
    cov_calculate_input='filtered_bam|corrected_bam',
    filter_target='no_mapq_filter|high_mapq_filter',
    flagger_types='error|duplicated|haploid|collapsed',

parts=16 # used for splitting an aligned bam
n_of_batches=3 # splits the read run into bins of 3
scattergather:
    split=parts,
    batches=n_of_batches

def get_both_asm(wildcards):
    return manifest_df.loc[(wildcards.sample, wildcards.tech), ['hap1','hap2']].to_dict()

def get_sample_run(what_type='fq'):
    def inner(wildcards):
        fofn_df = pd.read_csv(manifest_df.at[(wildcards.sample, wildcards.tech), 'fofn'], sep='\t', header=None, names=['file'])
        if what_type == 'fq':
            return fofn_df.at[int(wildcards.run_number), 'file']
        elif what_type == 'fai':
            return fofn_df.at[int(wildcards.run_number), 'file']+'.fai'
    return inner

def consolidate_parts(what_type='bam'):
    def inner(wildcards):
        fofn_df = pd.read_csv(manifest_df.at[(wildcards.sample, wildcards.tech), 'fofn'],sep='\t',header=None,names=['file'])
        run_number = fofn_df.index
        if what_type == 'log':
            pattern = 'results/{{sample}}/{{tech}}/secphase/phase_logs/phasing_out_run-{run_number}_batch-{scatteritem}.log'
        else:
            pattern = 'results/{{sample}}/{{tech}}/alignment/bam_batches/run-{run_number}/batch-{scatteritem}.bam'
        return gather.batches(pattern, run_number=run_number)
    return inner

def variant_calling_input(wildcards):
    # decision based on the tech of sample
    if wildcards.tech == 'hifi':
        return 'results/{sample}/{tech}/calls/deepvariant/{sample}.vcf.gz'
    else:
        return 'results/{sample}/{tech}/calls/pepper_margin_dv/{sample}.vcf.gz'

flagger_types = ["error", "duplicated", "haploid", "collapsed"]
def get_final_output(wildcards):
    cov_calculate_input = ["filtered_bam"]
    # cov_calculate_input = ["filtered_bam", "corrected_bam"]
    filter_target = ["high_mapq_filter"]
    # filter_target = ["no_mapq_filter","high_mapq_filter"]
    targets = [cov_calculate_input, filter_target]
    all = ["results/{{sample}}/{{tech}}/flagger/bed/{}/{}/{{sample}}.window_corrected_final.bed".format(*item) for item in product(*targets)]
    return expand(all, zip, sample=manifest_df.index.get_level_values('sample'), tech=manifest_df.index.get_level_values('tech'))

rule all:
    input:
        get_final_output

# -------- Step 0: Split the sequences -------- #
rule get_batch_ids:
    input:
        fai = get_sample_run(what_type='fai')
    output:
        batches = temp(scatter.batches('results/{{sample}}/{{tech}}/reads/run-{{run_number}}_batch-{scatteritem}.txt'))
    threads: 1
    resources:
        mem = 4,
        hrs = 8
    run:
        NIDS = len(output.batches)
        batch_dict = {}

        for i in range(NIDS):
            batch_dict[i] = []

        with open(input.fai, "r") as infile:
            fai_list = [line.split("\t")[0] for line in infile]

        for j in range(len(fai_list)):
            batch_dict[j % NIDS].append(fai_list[j])

        outs = [open(f,"w+") for f in output.batches]

        for i in range(NIDS):
            outs[i].write("\n".join(batch_dict[i]) + "\n")
            outs[i].close()

rule make_batch_fqs:
    input:
        fastq = get_sample_run(what_type='fq'),
        batch_file = 'results/{sample}/{tech}/reads/run-{run_number}_batch-{scatteritem}.txt'
    output:
        batch_fastq = temp('results/{sample}/{tech}/reads/run-{run_number}_batch-{scatteritem}.fastq'),
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.12",
    threads: 1
    resources:
        mem = 4,
        hrs = 8
    shell:
        '''
        samtools fqidx {input.fastq} -r {input.batch_file} > {output.batch_fastq}
        '''

# -------- Step 1: Align long reads to the diploid assembly -------- #
rule combine_asm:
    input:
        unpack(get_both_asm)
    output:
        combined = temp('results/{sample}/{tech}/both_haps.asm'),
        new_index =  temp('results/{sample}/{tech}/both_haps.asm.fai'),
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.12",
    threads: 1
    resources:
        mem = 4,
        hrs = 8,
        tmpdir=TMPDIR
    shell:
        '''
        starting_cmd=cat
        if [[ {input.hap1} == *.gz  ]]; then
            starting_cmd=zcat
        fi
        $starting_cmd {input.hap1} {input.hap2} > {output.combined}
        samtools faidx {output.combined}
        '''

rule make_kmers:
    input:
        asm = rules.combine_asm.output.combined,
    output:
        kmer_freq = 'results/{sample}/{tech}/repetitive_k15.txt',
        kmer_library = temp(directory('results/{sample}/{tech}/merylDB/'))
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "winnowmap/2.03",
    benchmark: 'benchmarks/{sample}_{tech}_make-kmers_benchmark.txt'
    threads: 1
    resources:
        mem = 4,
        hrs = 12
    shell:
        '''
        meryl count k=15 output {output.kmer_library} {input.asm} \
            && meryl print greater-than distinct=0.9998 {output.kmer_library} > {output.kmer_freq}
        '''

rule align_long_single_reads:
    # aligned reads will be sorted by query name as required for secphase input
    input:
        asm = rules.combine_asm.output.combined,
        kmer_freq = rules.make_kmers.output.kmer_freq,
        run_batch = 'results/{sample}/{tech}/reads/run-{run_number}_batch-{scatteritem}.fastq'
    output:
        qname_sorted_split_bam = 'results/{sample}/{tech}/alignment/bam_batches/run-{run_number}/batch-{scatteritem}.bam'
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "winnowmap/2.03",
        "samtools/1.12",
    benchmark: 'benchmarks/{sample}_{tech}-run_{run_number}-batch_{scatteritem}-align_long_single_reads_benchmark.txt'
    params:
        preset = lambda wildcards: config[wildcards.tech]['mapping_args']
    threads: 8
    resources:
        mem = 4,
        hrs = 24,
        tmpdir=TMPDIR
    shell:
        '''
        winnowmap -W {input.kmer_freq} -ax {params.preset} -Y -L --eqx --cs -I8g -t {threads} {input.asm} {input.run_batch} | samtools view -h - | samtools sort -T {resources.tmpdir} -n -@ {threads} -o {output.qname_sorted_split_bam} -
        '''

# -------- Step 2: Phase and relocalize the reads with secondary alignments using secphase -------- #
rule secphase_phase_reads:
    input:
        sorted_bam = rules.align_long_single_reads.output.qname_sorted_split_bam,
        fasta = rules.combine_asm.output.combined,
    output:
        flagged_regions = 'results/{sample}/{tech}/secphase/phase_logs/phasing_out_run-{run_number}_batch-{scatteritem}.log'
    container:
        'docker://mobinasri/secphase:v0.1'
    benchmark: 'benchmarks/{sample}_{tech}-{run_number}_batch-{scatteritem}_secphase_phase_reads_benchmark.txt'
    threads: 2
    resources:
        mem = 8,
        hrs = 12,
        tmpdir=TMPDIR
    params:
        secphase_args = lambda wildcards: config[wildcards.tech]['secphase_args']
    shell:
        '''
        secphase -q -c {params.secphase_args} \
            -i {input.sorted_bam} \
            -f {input.fasta} > {output.flagged_regions}
        '''

rule combine_aligned_bam:
    # combine the qname sorted splitted bams and coordinate sort
    input:
        aligned_bams = consolidate_parts(what_type='bam')
    output:
        qname_sorted_combined_bam = temp('results/{sample}/{tech}/alignment/qname_{sample}.bam'),
        coord_sorted_combined_bam = temp('results/{sample}/{tech}/alignment/coord_{sample}.bam')
    container: 'docker://quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2'
    benchmark: 'benchmarks/{sample}_{tech}-combine_aligned_bam_benchmark.txt'
    threads: 8
    resources:
        mem=2,
        hrs=12,
        tmpdir=TMPDIR
    shell:
        '''
        sambamba merge --nthreads={threads} {output.qname_sorted_combined_bam} {input.aligned_bams} && sambamba sort --show-progress --tmpdir={resources.tmpdir} --nthreads={threads} --out={output.coord_sorted_combined_bam} {output.qname_sorted_combined_bam}
        '''

rule secphase_post_correct:
    input:
        sorted_bam = rules.combine_aligned_bam.output.coord_sorted_combined_bam,
        phasing_log = consolidate_parts(what_type='log')
    output:
        corrected_bam = 'results/{sample}/{tech}/alignment/{sample}-coord_sorted-corrected.bam',
    container: 'docker://mobinasri/secphase:v0.1'
    benchmark: 'benchmarks/{sample}_{tech}-secphase_post_correct_benchmark.txt'
    threads: 8
    resources:
        mem = 2,
        hrs = 12,
        tmpdir=TMPDIR
    params:
        max_div=lambda wildcards: config[wildcards.tech]['max_divergence']
    shell:
        '''
        correct_bam \
            --primaryOnly \
            --maxDiv {params.max_div} \
            -n {threads} \
            -i {input.sorted_bam} \
            -P <(cat {input.phasing_log}) \
            -o {output.corrected_bam}

        samtools index -@ {threads} {output.corrected_bam}
        '''

# -------- Step 3: Call and filter variants  -------- #
rule get_split_beds:
    input:
        fai = rules.combine_asm.output.new_index,
    output:
        split_bed = temp(scatter.split('results/{{sample}}/{{tech}}/split_bed/{{sample}}_{scatteritem}.bed')),
    threads: 1
    resources:
        mem = 2,
        hrs=8
    container: flagger_container
    shell:
        '''
        out_dir=$(dirname {output.split_bed[0]}) \

        ## make a bed file that covers the whole assembly & split the bed file of the whole assembly into multiple bed files
        python3 /home/programs/src/split_bed_contig_wise.py \
            --bed <( cat {input.fai} | awk '{{print $1"\\t"0"\\t"$2}}' ) \
            --n {parts} \
            --dir $out_dir \
            --prefix {wildcards.sample}

        sleep 2

        ## change output names
        ls $out_dir | while read line; do new_name=${{line/./-of-{parts}.}}; mv "${{out_dir}}/${{line}}" "${{out_dir}}/${{new_name}}"; done
        '''

rule split_bam_into_contigs:
    input:
        bam = rules.secphase_post_correct.output.corrected_bam,
        bed = 'results/{sample}/{tech}/split_bed/{sample}_{scatteritem}.bed'
    output:
        split_bam = temp('results/{sample}/{tech}/split_bam/{sample}_{scatteritem}.bam'),
        split_bam_bai = temp('results/{sample}/{tech}/split_bam/{sample}_{scatteritem}.bam.bai')
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.12",
    benchmark: 'benchmarks/{sample}_{tech}_{scatteritem}-split_bam_into_contigs_benchmark.txt'
    threads: 4
    resources:
        mem=1,
        hrs=12
    shell:
        '''
        samtools view -@ {threads} -h -b -L {input.bed} {input.bam} > {output.split_bam}
        samtools index -@ {threads} {output.split_bam}
        '''

rule deepvariant:
    input:
        asm = rules.combine_asm.output.combined,
        bam = rules.split_bam_into_contigs.output.split_bam,
        bai= rules.split_bam_into_contigs.output.split_bam_bai
    output:
        vcf = temp('results/{sample}/{tech}/calls/deepvariant/{sample}_{scatteritem}.vcf'),
        vcf_gz = temp('results/{sample}/{tech}/calls/deepvariant/{sample}_{scatteritem}.vcf.gz'),
    threads: 8
    resources:
        mem = 2,
        hrs=12,
        tmpdir=TMPDIR
    container: 'docker://google/deepvariant:1.3.0'
    benchmark: 'benchmarks/{sample}_{tech}-deepvariant_benchmark-{scatteritem}.txt'
    shell:
        '''
        /opt/deepvariant/bin/run_deepvariant \
            --model_type="PACBIO" \
            --ref={input.asm} \
            --reads={input.bam} \
            --output_vcf={output.vcf} \
            --make_examples_extra_args="keep_supplementary_alignments=true,min_mapping_quality=0" \
            --call_variants_extra_args="use_openvino=true" \
            --num_shards={threads} \
            --dry_run=false

        gzip -c {output.vcf} > {output.vcf_gz}
        '''

rule pepper_margin_dv:
    input:
        asm = rules.combine_asm.output.combined,
        bam = rules.split_bam_into_contigs.output.split_bam,
        bai= rules.split_bam_into_contigs.output.split_bam_bai
    output:
        vcf_gz= temp('results/{sample}/{tech}/calls/pepper_margin_dv/{scatteritem}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz'),
    threads: 8
    resources:
        mem = 16,
        hrs=12,
        tmpdir=TMPDIR
    container: 'docker://kishwars/pepper_deepvariant:r0.8'
    benchmark: 'benchmarks/{sample}_{tech}_{scatteritem}-pepper_margin_dv_benchmark.txt'
    shell:
        '''
        run_pepper_margin_deepvariant call_variant \
          -b {input.bam} \
          -f {input.asm} \
          -o $( dirname {output.vcf_gz} ) \
          -t {threads} \
          --ont_r9_guppy5_sup \
          --pepper_include_supplementary \
          --dv_min_mapping_quality 0 \
          --pepper_min_mapping_quality 0 \
        '''

rule gather_vcfs:
    input:
        vcf_gz = lambda wildcards: gather.split('results/{{sample}}/{{tech}}/calls/{{caller}}/{{sample}}_{scatteritem}.vcf.gz') if wildcards.tech == 'hifi' else gather.split('results/{{sample}}/{{tech}}/calls/{{caller}}/{scatteritem}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz')
    output:
        merged_vcf = temp('results/{sample}/{tech}/calls/{caller}/{sample}-merged.vcf'),
        merged_vcf_gz = 'results/{sample}/{tech}/calls/{caller}/{sample}.vcf.gz'
    threads: 2
    resources:
        mem = 4,
        hrs=72
    container: 'docker://quay.io/biocontainers/bcftools:1.12--h45bccc9_1'
    benchmark: 'benchmarks/{sample}_{tech}_{caller}-gather_vcfs_benchmark.txt'
    shell:
        '''
        ## make a header for the merged vcf file
        zcat {input.vcf_gz[0]} | awk 'substr($0,1,1) == "#"' \
            | awk -v colname={wildcards.sample} '{{if ($1 == "#CHROM"){{ for(i =1; i < 10; i++){{printf("%s\\t",$i)}}; printf("%s\\n",colname)}} else {{print $0}}}}' > {output.merged_vcf} \
        && \
        zcat {input.vcf_gz} | awk 'substr($0,1,1) != "#"' >> {output.merged_vcf}

        ## sort the merged vcf file and produce the final gzipped vcf file
        bcftools sort -o {output.merged_vcf_gz} -Oz {output.merged_vcf}
        '''

# -------- Use biallelic SNVs to remove the alignments with alternative alleles -------- #
rule filter_biallelic_snps:
    input:
        vcf = variant_calling_input,
        bam = rules.secphase_post_correct.output.corrected_bam,
    output:
        filtered_vcf = 'results/{sample}/{tech}/calls/bi_snps-passed.vcf',
        filtered_bam = 'results/{sample}/{tech}/alignment/alt_filtered.bam',
        reads_removed_bam = 'results/{sample}/{tech}/alignment/alt.bam',
    params:
        filter_alt_reads = config['flagger']['filter_alt_reads'],
        gq = lambda wildcards: config[wildcards.tech]['genotype_quality']
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-filter_biallelic_snps_benchmark.txt'
    threads: 8
    resources:
        mem=8,
        hrs=72
    shell:
        '''
        ## Get the biallelic snps
        bcftools view -Ov -f PASS -m2 -M2 -v snps {input.vcf} -e 'FORMAT/VAF<0.5 | FORMAT/GQ<{params.gq}' > {output.filtered_vcf}
        
        ## Filter the reads that contain the alternative alleles of the snps with given VCF
        filter_alt_reads -i {input.bam} -o {output.filtered_bam} -f {output.reads_removed_bam} -v {output.filtered_vcf} -t {threads} {params.filter_alt_reads}
        
        ## Index
        samtools index -@ {threads} {output.filtered_bam}
        samtools index -@ {threads} {output.reads_removed_bam}
        '''

rule cov_calculate:
    input:
        bam = lambda wildcards: rules.filter_biallelic_snps.output.filtered_bam if wildcards.cov_calculate_input == 'filtered_bam' else rules.secphase_post_correct.output.corrected_bam,
        asm_index = rules.combine_asm.output.new_index,
    output:
        bam_cov = temp('results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/cov_calculate_{filter_target}.cov'),
        bam_cov_gz = 'results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/cov_calculate_{filter_target}.cov.tar.gz',
        bam_cov_counts = 'results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}.counts',
        bam_cov_mean = 'results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}-cov_mean.txt',
        bam_cov_sd = 'results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}-cov_sd.txt',
    params:
        min_mapq = lambda wildcards: 0 if wildcards.filter_target == 'no_mapq_filter' else 20
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-cov_calculate_{cov_calculate_input}_{filter_target}_benchmark.txt'
    threads: 1
    resources:
        mem=16,
        hrs=12
    shell:
        '''
        # Convert depth into compressed format
        depth2cov -d <( samtools depth -aa -Q {params.min_mapq} {input.bam} ) \
            -f {input.asm_index} -o {output.bam_cov} \
        && \
        # Convert cov to counts
        cov2counts -i {output.bam_cov} -o {output.bam_cov_counts}
        
        # Calculate mean and standard deviation
        python3 $CALC_MEAN_SD_PY \
            --countsInput {output.bam_cov_counts} \
            --meanOutput {output.bam_cov_mean} \
            --sdOutput {output.bam_cov_sd}
            
        tar -czf {output.bam_cov_gz} {output.bam_cov}
        '''

rule cov_calculate_by_window:
    input:
        bam_cov='results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/cov_calculate_{filter_target}.cov',
        asm_index=rules.combine_asm.output.new_index,
    output:
        covs=temp(directory('results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}_covs_windowed')),
        covs_gz='results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}_covs_windowed.tar.gz',
        counts=temp(directory('results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}_counts_windowed')),
        counts_gz='results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}_counts_windowed.tar.gz',
        windows='results/{sample}/{tech}/flagger/cov/{cov_calculate_input}/{filter_target}_windows.txt',
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-cov_calculate_{cov_calculate_input}_{filter_target}_benchmark.txt'
    threads: 1
    resources:
        mem=16,
        hrs=72
    shell:
        '''
        mkdir -p {output.covs}
        mkdir -p {output.counts}
        # Make a separate cov file for each contig
        split_cov_by_window -c {input.bam_cov} -f {input.asm_index} -p {output.covs}/{wildcards.sample} > {output.windows}
        
        # Count each window-specific cov file
        for f in $(readlink -f {output.covs}/*); do cov2counts -i $f -o ${{f/.cov/.counts}}; echo $f" finished";done
        mv {output.covs}/*.counts {output.counts}
        sleep 2
        
        tar -czf {output.covs_gz} --directory {output.covs} .
        tar -czf {output.counts_gz} --directory {output.counts} .
        '''

# -------- Flagger steps -------- #
rule fit_mixture_model:
    input:
        counts = rules.cov_calculate.output.bam_cov_counts,
        expected_cov = rules.cov_calculate.output.bam_cov_mean,
    output:
        table = 'results/{sample}/{tech}/flagger/bed/{cov_calculate_input}/{filter_target}_read_alignment.table'
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-{cov_calculate_input}_{filter_target}_fit_mixture_model_benchmark.txt'
    threads: 2
    resources:
        mem=8,
        hrs=72
    shell:
        '''
        python3 $FIT_GMM_PY \
            --counts {input.counts} \
            --cov $(cat {input.expected_cov}) \
            --output {output.table}
        '''

rule fit_mixture_model_windowed:
    input:
        counts = rules.cov_calculate_by_window.output.counts,
        windows_txt = rules.cov_calculate_by_window.output.windows,
        expected_cov = rules.cov_calculate.output.bam_cov_mean
    output:
        tables = directory('results/{sample}/{tech}/flagger/bed/{cov_calculate_input}/{filter_target}_tables')
    params:
        min_contig_size = 5*(10**6)
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-{cov_calculate_input}_{filter_target}_fit_mixture_model_windowed_benchmark.txt'
    threads: 16
    resources:
        mem=2,
        hrs=72
    shell:
        '''
        mkdir -p {output.tables}
        cat {input.windows_txt} | awk '{params.min_contig_size} <= $3' | xargs -n 3 -P {threads} sh -c 'python3 $FIT_GMM_PY --cov $(cat {input.expected_cov}) --counts {input.counts}/{wildcards.sample}.$0_$1_$2.counts --output {output.tables}/{wildcards.sample}.$0_$1_$2.table'
        '''

rule find_blocks_from_table:
    input:
        table = rules.fit_mixture_model.output.table,
        cov = rules.cov_calculate.output.bam_cov,
    output:
        beds = expand('results/{{sample}}/{{tech}}/flagger/bed/{{cov_calculate_input}}/{{filter_target}}_{{sample}}.{flagger_types}.bed', flagger_types=flagger_types)
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-{cov_calculate_input}_{filter_target}-find_blocks_from_table_benchmark.txt'
    threads: 2
    resources:
        mem=4,
        hrs=72
    shell:
        '''
         find_blocks_from_table \
            -c {input.cov} \
            -t {input.table} \
            -p $( dirname {output.beds[0]} )/{wildcards.filter_target}_{wildcards.sample}
        '''

rule find_blocks_by_windows:
    input:
        tables = rules.fit_mixture_model_windowed.output.tables,
        covs = rules.cov_calculate_by_window.output.covs,
        windows_txt = rules.cov_calculate_by_window.output.windows
    output:
        beds = expand('results/{{sample}}/{{tech}}/flagger/bed/{{cov_calculate_input}}/{{filter_target}}_beds/window_based.{flagger_types}.bed',flagger_types=flagger_types)
    params:
        min_win_size = 5 * (10 ** 6)
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-{cov_calculate_input}_{filter_target}-find_blocks_by_windows_benchmark.txt'
    threads: 8
    resources:
        mem=10,
        hrs=72
    shell:
        '''
        beds_dir=$( dirname {output.beds[0]} )
        mkdir -p $beds_dir
        cat {input.windows_txt} | awk '{params.min_win_size} <= $3' | xargs -n3 -P {threads} sh -c 'find_blocks_from_table -c {input.covs}/{wildcards.sample}.$0_$1_$2.cov -t {input.tables}/{wildcards.sample}.$0_$1_$2.table -p $( dirname {output.beds[0]} )/{wildcards.sample}.$0_$1_$2'
        
        cat $beds_dir/*.error.bed | bedtools sort -i - | bedtools merge -i - > $beds_dir/window_based.error.bed
        cat $beds_dir/*.duplicated.bed | bedtools sort -i - | bedtools merge -i - > $beds_dir/window_based.duplicated.bed
        cat $beds_dir/*.haploid.bed | bedtools sort -i - | bedtools merge -i - > $beds_dir/window_based.haploid.bed
        cat $beds_dir/*.collapsed.bed | bedtools sort -i - | bedtools merge -i - > $beds_dir/window_based.collapsed.bed
        
        find $beds_dir -name {wildcards.sample}.*.bed -delete 
        '''

rule combine_beds:
    input:
        blocks=rules.find_blocks_from_table.output.beds,
        windowed_blocks=rules.find_blocks_by_windows.output.beds,
    output:
        window_corrected = expand('results/{{sample}}/{{tech}}/flagger/bed/{{cov_calculate_input}}/{{filter_target}}_window_corrected/{flagger_types}.bed', flagger_types=flagger_types),
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-combine_beds_{cov_calculate_input}_{filter_target}_benchmark.txt'
    threads: 1
    resources:
        mem=16,
        hrs=72
    shell:
        '''
        out_dir=$(dirname {output.window_corrected[0]})
        first=$( dirname {input.blocks[0]} )
        second=$( dirname {input.windowed_blocks[0]} )
        fms=$out_dir/first_minus_second
        mkdir -p $fms
        cat $second/*.bed | sort -k1,1 -k2,2n | bedtools merge -i - > $second/second_all.bed 
        for c in error duplicated haploid collapsed
        do
            bedtools subtract -sorted -a <( sort -k1,1 -k2,2n $first/{wildcards.filter_target}_{wildcards.sample}.${{c}}.bed ) -b $second/second_all.bed > $fms/{wildcards.sample}.${{c}}.bed
            cat $fms/{wildcards.sample}.${{c}}.bed $second/*.${{c}}.bed | sort -k1,1 -k2,2n | bedtools merge -i - > $out_dir/${{c}}.bed
        done
        
        rm -rf $fms
        '''

rule get_final_bed:
    input:
        beds=rules.combine_beds.output.window_corrected
    output:
        final_bed = 'results/{sample}/{tech}/flagger/bed/{cov_calculate_input}/{filter_target}/{sample}.window_corrected_final.bed'
    container: flagger_container
    benchmark: 'benchmarks/{sample}_{tech}-get_final_bed_{cov_calculate_input}_{filter_target}_benchmark.txt'
    threads: 1
    resources:
        mem=8,
        hrs=12
    shell:
        '''
        tar_name="$( dirname {output.final_bed} )/window_corrected_beds.tar.gz"
        tar -czf $tar_name --directory $( dirname {input.beds[0]} ) .
        bash /home/scripts/combine_comp_beds.sh \
            -b $tar_name \
            -m /home/scripts/colors.txt \
            -t {wildcards.sample}_window_corrected \
            -o {output.final_bed}
        '''

onsuccess:
    shell("rm -rf results/*/*/alignment/bam_batches")