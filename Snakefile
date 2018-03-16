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
        return [formatstr.format(srr)
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

ruleorder: build_transcriptome > vcf_to_genome

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


rule star_map:
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

rule kallisto_quant:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        index='Reference/dmel_5.57_kallisto',
        dir=ancient('{sample}/')
    priority: 50
    params:
        reads=interleave_reads
    output: "{sample}/abundance.h5", "{sample}/abundance.tsv"
    shell:"""
    mkdir -p {wildcards.sample}
    ~/Downloads/kallisto/kallisto quant \
        -i {input.index} \
        -o {wildcards.sample} \
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

rule index_bam:
    input: "{sample}.bam"
    output: "{sample}.bam.bai"
    log: "{sample}.bam.bai_log"
    shell: "{module}; module load samtools; samtools index {input}"

rule dedup:
    input: "{sample}.bam"
    output: ("{sample}_dedup.bam")
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

rule bowtie2_build:
    input: "{base}.fasta"
    output: "{base}.1.bt2"
    log:    "{base}.bt2.log"
    shell: "{module}; module load bowtie2; bowtie2-build --offrate 3 {input} {wildcards.base}"


rule call_variants:
    input:
        ref_fasta="prereqs/d{reference}.fasta",
        ref_fai="prereqs/d{reference}.fasta.fai",
        ref_dict="prereqs/d{reference}.dict",
        bam="Reference/{reference}/{species}_gdna_bowtie2_dedup.bam",
        bai="Reference/{reference}/{species}_gdna_bowtie2_dedup.bam.bai",
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

