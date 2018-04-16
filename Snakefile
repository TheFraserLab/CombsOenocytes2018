configfile: "config.json"


mel_release = "r5.57_FB2014_03"
sim_release = "r2.01_FB2016_01"
sec_release = "r1.3_FB2015_01"

mel_version, mel_date = mel_release.split('_', 1)
sim_version, sim_date = sim_release.split('_', 1)
sec_version, sec_date = sec_release.split('_', 1)

num_mel_windows = 17
dates = {'mel': mel_date, 'sim': sim_date, 'sec': sec_date}
versions = {'mel': mel_version, 'sim': sim_version, 'sec': sec_version}

module = '''module () {
        eval `$LMOD_CMD bash "$@"`
        }'''

from os import path

def getreads(readnum):
    if readnum == 0:
        formatstr = 'sequence/{}.fastq.gz'
    else:
        formatstr = 'sequence/{{}}_{}.fastq.gz'.format(readnum)
    def retfun(wildcards):
        return [ancient(formatstr.format(srr))
                for srr in config['samples'][path.basename(wildcards.sample)]]
    return retfun

def getreads_nowc(readnum):
    if readnum == 0:
        formatstr = 'sequence/{}.fastq.gz'
    else:
        formatstr = 'sequence/{{}}_{}.fastq.gz'.format(readnum)
    def retfun(sample):
        return [formatstr.format(srr)
                for srr in config['samples'][path.basename(sample)]]
    return retfun

def getreadscomma(readnum):
    if readnum == 0:
        formatstr = 'sequence/{}.fastq.gz'
    else:
        formatstr = 'sequence/{{}}_{}.fastq.gz'.format(readnum)
    def retfun(wildcards):
        return ",".join([formatstr.format(srr)
                        for srr in config['samples'][path.basename(wildcards.sample)]])
    return retfun

def interleave_reads(wildcards):
    retval =  " ".join(" ".join([a, b])
            for a, b
            in zip(
                getreads_nowc(1)(path.basename(wildcards.sample)),
                getreads_nowc(2)(path.basename(wildcards.sample))
                )
            )
    return retval



rule all:
    input:
        expand('analysis/{target}/figure-specificity/oefemale.svg',
                target=['mel', 'sim', 'sec'])


rule prepare_emase:
    input:
        ancient("Reference/{species}/"),
        fasta="prereqs/d{species}.fasta",
        gtf="Reference/{species}_good_exonids.gtf"
    output:
        "Reference/{species}/emase.transcripts.fa",
        "Reference/{species}/emase.transcripts.info",
        "Reference/{species}/emase.gene2transcripts.tsv",
    shell: """
    source activate emase
    prepare-emase -G {input.fasta} -g {input.gtf} -o Reference/{wildcards.species}/ -m --no-bowtie-index


    """

rule bam_to_emase:
    input:
        bam="analysis/{genome}/{sample}/mapped.bam",
        transcriptome_info="Reference/{genome}/emase.transcripts.info",
    output:
        "analysis/{genome}/{sample}/aligned.transcriptome.h5"
    shell:"""
    source activate emase
    bam-to-emase -a {input.bam} \
            -i {input.transcriptome_info} \
            -s sim,sec \
            -o {output}
    """

rule run_emase:
    input:
        h5 = "analysis/{target}/{sample}/aligned.transcriptome.h5",
        transcripts = "Reference/{target}/emase.gene2transcripts.tsv",
        transcript_info = "Reference/{target}/emase.pooled.transcripts.info",
    output:
        "analysis/{target}/{sample}/emase.isoforms.effective_read_counts",
        "analysis/{target}/{sample}/emase.isoforms.tpm",
        "analysis/{target}/{sample}/emase.isoforms.alignment_counts",
        "analysis/{target}/{sample}/emase.genes.effective_read_counts",
        "analysis/{target}/{sample}/emase.genes.tpm",
        "analysis/{target}/{sample}/emase.genes.alignment_counts",
    shell: """ source activate emase
    run-emase -i {input.h5} \
            -g {input.transcripts} \
            -L {input.transcript_info} \
            -M 4 \
            -o analysis/{wildcards.target}/{wildcards.sample}/emase \
            -r 101 \
            -c
            """


rule block_gzip:
    input:
        vcf="{file}.gvcf"
    output:
        gz="{file}.gvcf.gz",
        gzi="{file}.gvcf.gz.gzi",
        tbi="{file}.gvcf.gz.tbi",
    shell:""" source activate emase
    bgzip --index --stdout {input.vcf} > {output.gz}
    tabix  -p vcf -f {output.gz}
    """

rule add_exonids:
    input: gtf="{file}.gtf", code="AddExonIDs.py"
    output: "{file}_exonids.gtf"
    shell: "python AddExonIDs.py {input.gtf} {output}"

rule vcf_to_genome:
    input:
        vcf_indels="Reference/{genome}/simsec_variants_on_{genome}_indels.gvcf.gz",
        vcf_snps="Reference/{genome}/simsec_variants_on_{genome}_snps.gvcf.gz",
        reffasta="Reference/d{genome}.fa",
        gtf="Reference/{genome}_good_exonids.gtf",
    output:
        chain="Reference/{genome}/ref-to-{target}.chain",
        patched="Reference/{genome}/{target}.patched.fa",
        transformed="Reference/{genome}/{target}.fa",
        gtf="Reference/{genome}/{target}.gtf",
        gtfdb="Reference/{genome}/{target}.gtf.db",
    shell:"""
    source activate emase
    g2gtools vcf2chain --pass -f {input.reffasta} -i {input.vcf_indels} -s {wildcards.target}_gdna -o {output.chain}
    g2gtools patch -i {input.reffasta} -s {wildcards.target}_gdna -v {input.vcf_snps} -o {output.patched}
    g2gtools transform -i {output.patched} -c {output.chain} -o {output.transformed}
    g2gtools convert -c {output.chain} -i {input.gtf} -f gtf -o {output.gtf}
    g2gtools gtf2db -i {output.gtf} -o {output.gtfdb}
    """

rule build_transcriptome:
    input:
        sim="Reference/{target}/sim.fa",
        sec="Reference/{target}/sec.fa",
    output:
           "Reference/{target}/emase.pooled.transcripts.fa",
           "Reference/{target}/emase.pooled.transcripts.info",
    shell:""" source activate emase
    prepare-emase \
            -G {input.sim},{input.sec} \
            -s sim,sec -o Reference/{wildcards.target} \
            --no-bowtie-index
    """

ruleorder: tophat_transcriptome > build_transcriptome > prepare_emase > vcf_to_genome

rule orig_seq_star_ref:
    input:
        fasta="Reference/{target}/simsec_corrected.fasta",
        gtf="Reference/{target}_good.gtf",
    output:
        outdir="Reference/{target}/orig",
        outfile="Reference/{target}/orig/Genome"
    priority: 50
    shell:""" {module}; module load STAR
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
	STAR --runMode genomeGenerate --genomeDir {output.outdir} \
        --runThreadN 12 \
        --limitGenomeGenerateRAM 48705389440 \
		--outTmpDir {output.outdir}/_tmp/ \
        --sjdbGTFfile {input.gtf} \
		--genomeFastaFiles {input.fasta}
    """

# I don't actually need to mask for the WASP pipeline, but it's safer to use
# the program that I know already works to generate the bed file.
rule mask_and_make_bed:
    input:
            var_tab='Reference/{target}/simsec_variants.tsv',
            ref_fasta="prereqs/d{target}.fasta",
    output:
            fasta="Reference/{target}/simsec_masked.fasta",
            corr_fasta="Reference/{target}/simsec_corrected.fasta",
            bed="Reference/{target}/simsec_variant.bed",
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python MaskReferenceFromGATKTable.py \
            --target-species sim \
            --emit-bed {output.bed} \
            --emit-corrected {output.corr_fasta} \
            --outfasta {output.fasta} \
                {input.ref_fasta} \
                {input.var_tab}
                """

rule diploid_star_ref:
    input:
        fasta="Reference/{target}/emase.pooled.transcripts.fa",
        gtf="Reference/{target}_good.gtf",
    output:
        outdir="Reference/{target}/emase.pooled",
        outfile="Reference/{target}/emase.pooled/Genome"
    priority: 50
    shell:""" {module}; module load STAR
    rm -rf {output.outdir}
    mkdir -p {output.outdir}
	STAR --runMode genomeGenerate --genomeDir {output.outdir} \
        --runThreadN 12 \
        --limitGenomeGenerateRAM 48705389440 \
		--outTmpDir {output.outdir}/_tmp/ \
		--genomeFastaFiles {input.fasta}
    """


rule star_map_orig:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        genome="Reference/{target}/orig/Genome",
        genomedir="Reference/{target}/orig/",
    params:
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
    output: "analysis/{target}/{sample}/orig_mapped.bam"
    threads: 16
    shell: """{module}; module load STAR
    STAR \
    --runThreadN 16 \
    --runMode alignReads \
    --readFilesCommand zcat \
    --genomeDir {input.genomedir} \
    --outFileNamePrefix analysis/{wildcards.target}/{wildcards.sample}/orig_ \
    --outSAMattributes All \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn {params.r1s} {params.r2s}
    if [ -s  analysis/{wildcards.target}/{wildcards.sample}/orig_Aligned.sortedByCoord.out.bam ]; then
          mv analysis/{wildcards.target}/{wildcards.sample}/orig_Aligned.sortedByCoord.out.bam {output};
    fi
    """

rule make_snpdir:
    input:
        vcf="Reference/{target}/simsec_variants_on_{target}.gvcf.gz"
    output:
        dir="Reference/{target}/snpdir",
        file="Reference/{target}/snpdir/all.txt.gz",
    shell:"""
    mkdir -p {output.dir}
    gzip -dc {input.vcf}             \
            | grep -v "^#"             \
            | awk 'BEGIN {{OFS="\t"}}; length($4) == 1 && length($5) == 1 {{print $1,$2,$4,$5}};' \
            | gzip -c  \
            > {output.file}
"""

rule dedup:
    input: "{sample}.bam"
    output: "{sample}.dedup.bam"
    log: "{sample}_dedup.log"
    shell: """{module}; module load picard
    picard MarkDuplicates \
        SORTING_COLLECTION_SIZE_RATIO=.05 \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        MAX_RECORDS_IN_RAM=2500000 \
		READ_NAME_REGEX=null \
		REMOVE_DUPLICATES=true \
		DUPLICATE_SCORING_STRATEGY=RANDOM \
		INPUT={input} OUTPUT={output} METRICS_FILE={log}
        """


rule wasp_find_snps:
    input:
        bam="analysis/{genome}/{sample}/{prefix}.dedup.bam",
        bai="analysis/{genome}/{sample}/{prefix}.dedup.bam.bai",
        snpdir="Reference/{genome}/snpdir",
        snpfile="Reference/{genome}/snpdir/all.txt.gz"
    output:
        temp("analysis/{genome}/{sample}/{prefix}.dedup.remap.fq1.gz"),
        temp("analysis/{genome}/{sample}/{prefix}.dedup.remap.fq2.gz"),
        temp("analysis/{genome}/{sample}/{prefix}.dedup.keep.bam"),
        temp("analysis/{genome}/{sample}/{prefix}.dedup.to.remap.bam"),

    shell:
        """python ~/FWASP/mapping/find_intersecting_snps.py \
            --progressbar \
            --phased --paired_end \
            {input.bam} {input.snpdir}
        """


rule wasp_remap:
    input:
        R1="analysis/{genome}/{sample}/{prefix}.remap.fq1.gz",
        R2="analysis/{genome}/{sample}/{prefix}.remap.fq2.gz",
        genome="Reference/{genome}/orig/Genome",
        genomedir="Reference/{genome}/orig/"
    output:
        temp("analysis/{genome}/{sample}/{prefix}.remap.bam")
    threads: 16
    shell: """{module}; module load STAR;
    rm -rf analysis/{wildcards.genome}/{wildcards.sample}/STARtmp
    STAR \
            --genomeDir {input.genomedir} \
            --outFileNamePrefix analysis/{wildcards.genome}/{wildcards.sample}/remap \
            --outSAMattributes MD NH --clip5pNbases 6 \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis/{wildcards.genome}/{wildcards.sample}/STARtmp \
            --limitBAMsortRAM 20000000000 \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn {input.R1} {input.R2}
    mv analysis/{wildcards.genome}/{wildcards.sample}/remapAligned.out.bam {output}
            """

rule wasp_keep:
    input:
        toremap="analysis/{file}.to.remap.bam",
        remapped="analysis/{file}.remap.bam",
    output:
        temp("analysis/{file}.remap.kept.bam"),
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python ~/FWASP/mapping/filter_remapped_reads.py \
            -p \
            {input.toremap} {input.remapped} \
            {output} """

rule wasp_merge:
    input:
        "analysis/{file}.remap.kept.bam",
        "analysis/{file}.keep.bam",
    output:
        temp("analysis/{file}.keep.merged.bam")
    shell:
        "{module}; module load samtools; samtools merge {output} {input}"


rule star_map_emase:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        genome="Reference/{target}/emase.pooled/Genome",
        genomedir="Reference/{target}/emase.pooled/",
    params:
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
    output: "analysis/{target}/{sample}/mapped.bam"
    threads: 6
    log: "{sample}/assigned_dmelR.log"
    shell: """{module}; module load STAR
    STAR \
    --runThreadN 16 \
    --runMode alignReads \
    --readFilesCommand zcat \
    --genomeDir {input.genomedir} \
    --outFileNamePrefix analysis/{wildcards.target}/{wildcards.sample}/ \
    --outSAMattributes All \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn {params.r1s} {params.r2s}
    if [ -s  analysis/{wildcards.target}/{wildcards.sample}/Aligned.sortedByCoord.out.bam ]; then
          mv analysis/{wildcards.target}/{wildcards.sample}/Aligned.sortedByCoord.out.bam {output};
    fi
    """


rule makedir:
    output: "{prefix}.log", "{prefix}/"
    shell: "touch {wildcards.prefix}.log; mkdir -p {wildcards.prefix}"

rule kallisto_index_unsplit:
    input:
        fasta="Reference/{target}/simsec_corrected_transcripts.fa"
    output:
        'Reference/{target}/kallisto_unsplit'
    shell:"""~/Downloads/kallisto/kallisto index\
        --index {output} \
        {input.fasta}
        """

rule bowtie_build_corrected:
    input:
        fasta="Reference/{target}/simsec_corrected.fasta"
    output:
        expand("Reference/{{target}}/simsec_corrected.{n}.bt2", n=[1,2,3,4])
    shell:"""
    {module}
    module load bowtie2
    bowtie2-build --threads 16 {input.fasta} Reference/{wildcards.target}/simsec_corrected
    """

rule tophat_transcriptome:
    input:
        "Reference/{target}/simsec_corrected.1.bt2",
        gtf="Reference/{target}_good.gtf"
    output:
        "Reference/{target}/simsec_corrected_transcripts.fa",
        "Reference/{target}/simsec_corrected_transcripts.fa.tlst",
        "Reference/{target}/simsec_corrected_transcripts.gff",
    shell:""" {module}
    module load tophat bowtie2
    tophat2 \
            --GTF {input.gtf} \
            --num-threads 16 \
            --transcriptome-index Reference/{wildcards.target}/simsec_corrected_transcripts \
            Reference/{wildcards.target}/simsec_corrected
    """


rule kallisto_index:
    input:
        fasta="Reference/{target}/emase.pooled.transcripts.fa"
    output:
        'Reference/{target}/kallisto'
    shell:"""~/Downloads/kallisto/kallisto index\
        --index {output} \
        {input.fasta}

        """


rule kallisto_quant_unsplit:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        index='Reference/{target}/kallisto_unsplit',
        dir=ancient('analysis/{target}/{sample}/unsplit/')
    priority: 50
    params:
        reads=interleave_reads
    output: "analysis/{target}/{sample}/unsplit/abundance.h5", "analysis/{target}/{sample}/unsplit/abundance.tsv"
    shell:"""
    mkdir -p {wildcards.sample}
    ~/Downloads/kallisto/kallisto quant \
        -i {input.index} \
        -o analysis/{wildcards.target}/{wildcards.sample}/unsplit \
        {params.reads}
    """


rule kallisto_quant:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        index='Reference/{target}/kallisto',
        dir=ancient('analysis/{target}/{sample}/')
    priority: 50
    params:
        reads=interleave_reads
    output: "analysis/{target}/{sample}/abundance.h5", "analysis/{target}/{sample}/abundance.tsv"
    shell:"""
    mkdir -p {wildcards.sample}
    ~/Downloads/kallisto/kallisto quant \
        -i {input.index} \
        -o analysis/{wildcards.target}/{wildcards.sample}/ \
        {params.reads}
    """


rule get_combined_variants:
    input:
        ref_fasta="prereqs/d{parent}.fasta",
        sim_gvcf="Reference/{parent}/sim_gdna_raw_variants_uncalibrated.p.g.vcf",
        sec_gvcf="Reference/{parent}/sec_gdna_raw_variants_uncalibrated.p.g.vcf",

    output:
        tsv="Reference/{parent}/simsec_variants.tsv",
        vcf="Reference/{parent}/simsec_variants_on_{parent}.gvcf"
    params:
        dir="Reference/{parent}/"
    shell: """
    {module}; module load java
	gatk -T GenotypeGVCFs \
		-R {input.ref_fasta} \
		-V {input.sim_gvcf} \
		-V {input.sec_gvcf} \
		-o {output.vcf}
	gatk -T VariantsToTable \
		-R {input.ref_fasta} \
		-V {output.vcf} \
		-F CHROM -F POS -F REF -F ALT -F QUAL \
		-F HET -F HOM-REF -F HOM-VAR -F NCALLED \
		-GF GT \
		-o {output.tsv}
    """

rule split_variants:
    input:
        "{variants}.gvcf"
    output:
        indels="{variants}_indels.gvcf",
        snps="{variants}_snps.gvcf"
    shell:"""
    {module}; module load bcftools

    bcftools view --types snps   -o {output.snps}   {input}
    bcftools view --types indels -o {output.indels} {input}

    """



rule map_gdna:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        ancient("Reference/{reference}/"),
        bt2_index="prereqs/d{reference}.1.bt2"
    output:
        "Reference/{reference}/{sample}_bowtie2.bam"
    params:
        index="prereqs/d{reference}",
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
        outdir= lambda wildcards, output: path.dirname(output[0])
    threads: 12
    shell: """{module}; module load samtools/1.3 bowtie2
    bowtie2 \
		--very-sensitive-local \
		-p 11 \
		--rg-id {wildcards.sample} \
		--rg "SM:{wildcards.sample}" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x {params.index} \
		-1 {params.r1s} \
		-2 {params.r2s} \
		| samtools view -b \
		| samtools sort -o {output} -T {params.outdir}/{wildcards.sample}_bowtie2_sorting
        """

rule sort_bam:
    input: "{sample}.bam"
    output: "{sample}.sort.bam"
    shell: "{module}; module load samtools; samtools sort -o {output} {input}"

rule index_bam:
    input: "{sample}.bam"
    output: "{sample}.bam.bai"
    log: "{sample}.bam.bai_log"
    shell: "{module}; module load samtools; samtools index {input}"


rule bowtie2_build:
    input: "{base}.fasta"
    output: "{base}.1.bt2"
    log:    "{base}.bt2.log"
    shell: "{module}; module load bowtie2; bowtie2-build --offrate 3 {input} {wildcards.base}"

ruleorder: bowtie_build_corrected > bowtie2_build


rule call_variants:
    input:
        ref_fasta="prereqs/d{reference}.fasta",
        ref_fai="prereqs/d{reference}.fasta.fai",
        ref_dict="prereqs/d{reference}.dict",
        bam="Reference/{reference}/{species}_gdna_bowtie2.dedup.bam",
        bai="Reference/{reference}/{species}_gdna_bowtie2.dedup.bam.bai",
    output:
        "Reference/{reference}/{species}_gdna_raw_variants_uncalibrated.p.g.vcf"
    threads: 4
    shell: """
    {module}; module load java
	gatk -T HaplotypeCaller \
		-R {input.ref_fasta} \
		-I {input.bam} \
		-nct 16 \
		--genotyping_mode DISCOVERY \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o {output}
    """

rule get_sra:
    output:
        "sequence/{srr}_1.fastq.gz",
        "sequence/{srr}_2.fastq.gz"
    #log: "sequence/{srr}.log"
    resources: max_downloads=1
    shell: """{module}; module load sra-tools

    fastq-dump --gzip --split-3 --outdir sequence {wildcards.srr}
    """

rule extra_fasta: # Some things DEMAND a .fa extension. Fuckers.
    output: "Reference/d{species}.fa"
    input: "prereqs/d{species}.fasta"
    shell: "cp {input} {output} """

rule ref_genome:
    output: "prereqs/d{species}.fasta"
    log: "prereqs/d{species}.log"
    params:
        date=lambda wildcards: dates[wildcards.species],
        version=lambda wildcards: versions[wildcards.species]
    shell: """{module}; #module load wget
    mkdir -p prereqs
	wget -O {output}.gz ftp://ftp.flybase.org/releases/{params.date}/d{wildcards.species}_{params.version}/fasta/d{wildcards.species}-all-chromosome-{params.version}.fasta.gz
	gunzip --force {output}.gz
    """

rule fadict:
    input: "{file}.fasta"
    output: "{file}.dict"
    shell: "{module}; module load picard; picard CreateSequenceDictionary R={input} O={output}"

rule index_fasta:
    input: "{file}.fasta"
    output: "{file}.fasta.fai"
    shell: "samtools faidx {input}"

rule ref_gff:
    output: "prereqs/d{species}.gff"
    log: "prereqs/d{species}.log"
    params:
        date=lambda wildcards: dates[wildcards.species],
        version=lambda wildcards: versions[wildcards.species]

    shell: """{module}; #module load wget
    mkdir -p prereqs
	wget -O {output}.gz ftp://ftp.flybase.org/releases/{params.date}/d{wildcards.species}_{params.version}/gff/d{wildcards.species}-all-{params.version}.gff.gz
	gunzip --force {output}.gz
    """

rule all_gtf:
    input: "prereqs/d{species}.gff"
    output: "Reference/{species}_all.gtf"
    log: "Reference/{species}_all.log"
    shell: """{module}; module load cufflinks
    gffread {input} -C -E -T -o- \
        > {output}
    """

rule good_gtf:
    input: "Reference/{species}_all.gtf"
    output: "Reference/{species}_good.gtf"
    log: "Reference/{species}_good.log"
    shell: """
	cat {input} \
		| grep -vP '(snoRNA|CR[0-9]{{4}}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		| grep 'gene_id' \
		> {output}
        """

rule bad_gtf:
    input: "Reference/{species}_all.gtf"
    output: "Reference/{species}_bad.gtf"
    log: "Reference/{species}_bad.log"
    shell: """
    echo {threads}
	cat {input} \
		| grep -vP '(snoRNA|CR[0-9]{{4}}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:|His.*:)' \
		> {output}
        """

rule sample_gene_ase_pval:
    input:
        bam="analysis/{target}/{sample}/orig_mapped.dedup.keep.merged.sort.bam",
        bai="analysis/{target}/{sample}/orig_mapped.dedup.keep.merged.sort.bam.bai",
        variants="Reference/{target}/simsec_variant.bed",
        gtf="Reference/{target}_good.gtf",
        sentinel=path.join("analysis", 'recalc_ase')
    threads: 1
    output:
        "analysis/{target}/{sample}/gene_ase_pval.tsv"
    shell: """ export PYTHONPATH=$HOME/ASEr/
    python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --id-name gene_name \
        --ase-function log10pval \
        --min-reads-per-allele 0 \
        {input.variants} \
        {input.gtf} \
        {input.bam}
    """

rule sample_gene_ase:
    input:
        bam="analysis/{target}/{sample}/orig_mapped.dedup.keep.merged.sort.bam",
        bai="analysis/{target}/{sample}/orig_mapped.dedup.keep.merged.sort.bam.bai",
        variants="Reference/{target}/simsec_variant.bed",
        gtf="Reference/{target}_good.gtf",
        sentinel=path.join("analysis", 'recalc_ase')
    threads: 1
    output:
        "analysis/{target}/{sample}/gene_ase_by_read.tsv"
    shell: """ export PYTHONPATH=$HOME/ASEr/
    python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --assign-all-reads \
        --id-name gene_name \
        --ase-function log2offset \
        --min-reads-per-allele 0 \
        {input.variants} \
        {input.gtf} \
        {input.bam}
    """

rule nowasp_gene_ase:
    input:
        bam="analysis/{target}/{sample}/orig_mapped.dedup.sort.bam",
        bai="analysis/{target}/{sample}/orig_mapped.dedup.sort.bam.bai",
        variants="Reference/{target}/simsec_variant.bed",
        gtf="Reference/{target}_good.gtf",
        sentinel=path.join("analysis", 'recalc_ase')
    threads: 1
    output:
        "analysis/{target}/{sample}/gene_ase_nowasp.tsv"
    shell: """ export PYTHONPATH=$HOME/ASEr/
    python ~/ASEr/bin/GetGeneASEbyReads.py \
        --outfile {output} \
        --assign-all-reads \
        --id-name gene_name \
        --ase-function pref_index \
        --min-reads-per-allele 0 \
        {input.variants} \
        {input.gtf} \
        {input.bam}
    """

rule ase_summary:
    input:
        lambda wc: expand('analysis/{target}/{sample}/gene_ase_by_read.tsv',
                target=[wc.target],
                sample=[c for c in config['samples'] if 'gdna' not in c])
    output:
        'analysis/{target}/summary_ase.tsv'
    log:
        'analysis/{target}/mst.log'
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python MakeSummaryTable.py \
	   --filename gene_ase_by_read.tsv \
	   --key gene \
       --out-basename summary_ase \
	   --column ase_value \
		analysis/{wildcards.target}/ \
		| tee {log}
        """

rule ase_refalt_summary:
    input:
        lambda wc: expand('analysis/{target}/{sample}/gene_ase_by_read.tsv',
                target=[wc.target],
                sample=[c for c in config['samples'] if 'gdna' not in c])
    output:
        'analysis/{target}/summary_ase_refalt.tsv'
    log:
        'analysis/{target}/mst_refalt.log'
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python MakeSummaryTable.py \
	   --filename gene_ase_by_read.tsv \
	   --key gene \
       --refalt \
       --out-basename summary_ase_refalt \
       --float-format %5.0f \
	   --column ref_counts \
		analysis/{wildcards.target}/ \
		| tee {log}
        """


rule nowasp_ase_summary:
    input:
        lambda wc: expand('analysis/{target}/{sample}/gene_ase_nowasp.tsv', target=[wc.target], sample=[c for c in config['samples'] if 'gdna' not in c])
    output:
        'analysis/{target}/summary_ase_nowasp.tsv'
    log:
        'analysis/{target}/mst.log'
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python MakeSummaryTable.py \
	   --filename gene_ase_nowasp.tsv \
	   --key gene \
       --out-basename summary_ase_nowasp \
	   --column ase_value \
		analysis/{wildcards.target}/ \
		| tee {log}
        """

rule pval_summary:
    input:
        lambda wc: expand('analysis/{target}/{sample}/gene_ase_pval.tsv', target=[wc.target], sample=[c for c in config['samples'] if 'gdna' not in c])
    output:
        'analysis/{target}/summary_ase_pvals.tsv'
    log:
        'analysis/{target}/mst.log'
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python MakeSummaryTable.py \
	   --filename gene_ase_pval.tsv \
	   --key gene \
       --out-basename summary_ase_pvals \
	   --column ase_value \
		analysis/{wildcards.target}/ \
		| tee {log}
        """

rule tissue_pvals:
    input:
        lambda wc: expand('analysis/{target}/{sample}/gene_ase_by_read.tsv', target=[wc.target], sample=[c for c in config['samples'] if 'gdna' not in c])
    output:
        'analysis/{target}/summary_ase_tissue_pvals.tsv'
    shell:"""
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python CombinedASEbySample.py -o {output} {input}
    """


rule kallisto_summary:
    input:
        lambda wc: expand('analysis/{target}/{sample}/unsplit/abundance.tsv', target=[wc.target], sample=[c for c in config['samples'] if 'gdna' not in c])
    output:
        'analysis/{target}/summary_kallisto_in_unsplit.tsv'
    log:
        'analysis/{target}/mst.log'
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python MakeSummaryTable.py \
	   --filename abundance.tsv \
       --in-subdirectory unsplit \
	   --key target_id \
       --out-basename summary_kallisto \
	   --column tpm \
		analysis/{wildcards.target}/ \
		| tee {log}
        """

rule kallisto_summary_renamed:
    input:
        expr = 'analysis/{target}/summary_kallisto_in_unsplit.tsv',
        info = 'Reference/{target}/simsec_corrected_transcripts.gff',
        tlst = 'Reference/{target}/simsec_corrected_transcripts.fa.tlst',
    output: 'analysis/{target}/summary_kallisto_genes.tsv'
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python CombineExprByGenes.py \
        --cufflinks-tlst {input.tlst} \
        --use-gtf --use-gene-name \
        --outfile {output} \
        {input.info} {input.expr}
        """

rule draw_specificity:
    input:
        code="DrawSpecificity.py",
        outdir=ancient('analysis/{target}/figure-specificity/'),
        expr='analysis/{target}/summary_kallisto_genes.tsv',
        ase='analysis/{target}/summary_ase.tsv',
        pvals='analysis/{target}/deseq_pvals.tsv',
    output:
        expand('analysis/{{target}}/figure-specificity/{tissue}{sex}.svg',
        tissue=['oe', 'fb'], sex=['male', 'female'],)
    shell:  """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python DrawSpecificity.py \
            --pvals {input.pvals} \
            --try-orthologs prereqs \
            {input.expr} {input.ase} {input.outdir}
        """


rule alignment:
    input:
        A="analysis/targets/{gene}/{A}.fasta",
        B="analysis/targets/{gene}/{B}.fasta",
    wildcard_constraints:
        A="[a-z][a-z][a-z]+",
        B="[a-z][a-z][a-z]+",
    output:
        "analysis/targets/{gene}/{A}_{B}.needleall"
    log:
        "analysis/targets/{gene}/{A}_{B}.needleall.log"
    shell: """{module}; module load  EMBOSS
    needleall -aseq {input.A} -bseq {input.B} \
            -aformat3 srspair -gapopen 10.0 -gapextend 0.5 \
            -outfile {output}"""

rule combined_fasta:
    input:
        A="analysis/targets/{gene}/{A}.fasta",
        B="analysis/targets/{gene}/{B}.fasta",
    wildcard_constraints:
        A="[a-z][a-z][a-z]+",
        B="[a-z][a-z][a-z]+",
    output:
        "analysis/targets/{gene}/{A}_{B}.fasta"
    log:
        "analysis/targets/{gene}/{A}_{B}.fasta.log"
    shell:
        "cat {input.A} {input.B} > {output}"

rule melsimsec_fastas:
    input:
        dir=ancient("analysis/targets/{region}/mss/"),
        mel="analysis/targets/{region}/mel.fasta",
        sim="analysis/targets/{region}/sim.fasta",
        sec="analysis/targets/{region}/sec.fasta",
    output:
        #dynamic("analysis/targets/{region}/mss/region_{rnum}.fasta")
    shell:"""
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    python PrepareForClustal.py \
            --outdir {input.dir} \
            {input.mel} {input.sim} {input.sec}
            """

rule clustalo_align:
    input:
        "analysis/targets/{region}/mss/region_{rnum}.fasta"
    output:
        "analysis/targets/{region}/mss/region_{rnum}.clu"
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
    clustalo \
            --in {input} --out {output} --outfmt=clustal

    """

rule mel_intergenic_bed:
    input:
        gtf="Reference/mel_all.gtf",
        dir=ancient("analysis/targets/{upstream}_{downstream}_intergenic/"),
    output:
        outbed="analysis/targets/{upstream}_{downstream}_intergenic/mel.bed"
    run:
        starts = []
        stops = []
        active = False
        with open(output.outbed, 'w') as outf:
            for line in open(input.gtf):
                data = line.split('\t')
                if len(data) < 8:
                    import sys
                    sys.stderr.write(line)
                    continue
                if data[2] != 'CDS':
                    continue
                start = int(data[3])
                stop = int(data[4])
                if wildcards.upstream in data[-1]:
                    active=True
                    if stops:
                        stops[0] = max(stop, stops[0])
                    else:
                        stops.append(stop)
                elif wildcards.downstream in data[-1] and active:
                    starts.append(start)
                    active = False
                    i_start = 0
                    i_stop = 0
                    region_n = 1
                    while i_stop < len(stops) and i_start < len(starts):
                        start = starts[i_start]
                        stop = stops[i_stop]
                        if start <= stop:
                            i_start += 1
                        else:
                            print(data[0],
                                  stop, start,
                                  '{}_{}_{}'.format(
                                      wildcards.upstream,
                                      wildcards.downstream,
                                      region_n,
                                      ),
                                  '.',  '+',
                                  stop, start,
                                  sep='\t',
                                  file=outf,
                                  )
                            region_n += 1
                            i_stop += 1
                elif active:
                    starts.append(start)
                    stops.append(stop)
                else:
                    pass


rule non_mel_bed:
    input:
        fa="analysis/targets/{gene}/mel.fasta",
        blastdb="Reference/d{species}.fa.nhr"
    output: "analysis/targets/{gene}/{species}.bed"
    #wildcard_constraints: species="^(?!mel$).*$"
    shell: """{module}; module load blast bioawk
    blastn -db Reference/d{wildcards.species}.fa \
            -outfmt "6 sseqid sstart send qseqid evalue sstrand length qlen slen qstart qend" \
            -gapextend 0\
            -query {input.fa} \
        | awk '!_[$4]++' \
        | bioawk -t '$10 > 1 && $6 ~ /minus/ {{$2 += $10 + 1}}; \
                    $10 > 1 && $6 ~ /plus/ {{$2 -= $10 + 1}}; \
                    $11 < $8 && $6 ~ /minus/ {{$3 -= ($8 - $11) + 1}}; \
                    $2 > $3 {{gsub("mel", "{wildcards.species}", $4); print $1,$3+1,$2,$4,$7/($8+1),"-", $3, $2 }}; \
                    $2 < $3 {{gsub("mel", "{wildcards.species}", $4); print $1,$2,$3+1,$4,$7/($8+1),"+", $2, $3 }}; '\
        > {output}
        """

rule mel_bed:
    input:
        gtf="Reference/mel_good.gtf",
        melsize="Reference/dmel.chr.sizes",
        oreganno="Reference/oreganno.prepend.bed",
        dnase="Reference/binding/dnase_peaks_prepend.bed",
    output:
        "analysis/targets/{gene}/mel.bed"
    shell: """{module}; module load bioawk bedtools
        mkdir -p `dirname output`
		grep '"{wildcards.gene}"' {input.gtf} \
		| bedtools sort \
		| bedtools merge  \
		| bedtools window -w 10000 -b - -a {input.dnase} \
        | bioawk -t '{{print $1, $2, $3, "mel_{wildcards.gene}_" NR, "0", "+", $2, $3}}' \
        | uniq --skip-fields 4 \
		> {output}

		grep '"{wildcards.gene}"' {input.gtf} \
		| bedtools sort \
		| bedtools merge  \
		| bedtools window -w 10000 -b - -a {input.dnase} \
        | bioawk -t '{{print $1, $2, $3, "mel_oreganno_{wildcards.gene}_" NR, "0", "+", $2, $3}}' \
        | uniq --skip-fields 4 \
        >> {output}
    """

rule nonmel_fasta_from_bed:
    input:
        bed='analysis/targets/{gene}/{species}.bed',
        full_fasta='Reference/d{species}.fa',
    output:
        'analysis/targets/{gene}/{species}.fasta'
    shell: """ {module}; module load bedtools
    bedtools getfasta -fi {input.full_fasta} -bed {input.bed} \
            -fo {output} -s -name
    """
rule make_blastdb:
    input: "Reference/{file}.fa"
    output: "Reference/{file}.fa.nhr"
    shell: """{module}; module load blast; makeblastdb -dbtype nucl -in {input}"""

ruleorder: mel_intergenic_bed > mel_bed > non_mel_bed

rule deseq2:
    input:
        ancient('analysis/{target}/combined/'),
        code='DESeq.R',
        data="analysis/{target}/summary_ase_refalt.tsv",
    output:
        expand('analysis/{{target}}/combined/{file}_deseq.tsv',
                file=['oefemale', 'oemale', 'fbfemale', 'fbmale'])
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate deseq
    Rscript {input.code} {input.data}
    """

rule combine_deseq_pvals:
    input:
        expand('analysis/{{target}}/combined/{file}_deseq.tsv',
                file=['oefemale', 'oemale', 'fbfemale', 'fbmale'])
    output:
        'analysis/{target}/deseq_pvals.tsv'
    run:
        import pandas as pd
        from os import path
        outdict = {}
        for fname in input:
            indf = pd.read_table(fname, index_col=0)
            colname = fname.split('/')[-1].split('_')[0]
            outdict[colname] = pd.np.log10(indf.padj) * pd.np.sign(indf.log2FoldChange)
        outdf = pd.DataFrame(outdict)
        outdf.to_csv(output[0], sep='\t', na_rep='0')




