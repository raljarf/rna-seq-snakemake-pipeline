SAMPLES = glob_wildcards(
    "../../learning/real-fastq-practice/{sample}_1.fastq.gz"
).sample

READS = ["1", "2"]


rule all:
    input:
        # FastQC reports for raw reads
        expand(
            "fastqc_results/{sample}_{read}_fastqc.html",
            sample=SAMPLES,
            read=READS
        ),
        expand(
            "fastqc_results/{sample}_{read}_fastqc.zip",
            sample=SAMPLES,
            read=READS
        ),

        # Trimmed FASTQ files
        expand(
            "trimmed/{sample}_{read}.trimmed.fastq.gz",
            sample=SAMPLES,
            read=READS
        ),

        # FastQC reports for trimmed reads
        expand(
            "fastqc_trimmed/{sample}_{read}.trimmed_fastqc.html",
            sample=SAMPLES,
            read=READS
        ),
        expand(
            "fastqc_trimmed/{sample}_{read}.trimmed_fastqc.zip",
            sample=SAMPLES,
            read=READS
        ),

        # STAR alignment: Generate coordinate-sorted BAM files for all samples
        expand(
            "aligned/{sample}_Aligned.sortedByCoord.out.bam",
            sample=SAMPLES
        ),

        # BAM indexing: Generate .bai index files for all samples
        expand(
            "aligned/{sample}_Aligned.sortedByCoord.out.bam.bai",
            sample=SAMPLES
        ),
 
       # Gene quantification: Generate gene count tables for all samples
        expand(
            "counts/{sample}_gene_counts.txt",
            sample=SAMPLES
        ),

        # Final combined gene count matrix
        "counts/gene_count_matrix.tsv",


rule fastqc:
    input:
        "../../learning/real-fastq-practice/{sample}_{read}.fastq.gz"
    output:
        html="fastqc_results/{sample}_{read}_fastqc.html",
        zip="fastqc_results/{sample}_{read}_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc {input} --outdir fastqc_results"


rule trim:
    input:
        r1="../../learning/real-fastq-practice/{sample}_1.fastq.gz",
        r2="../../learning/real-fastq-practice/{sample}_2.fastq.gz"
    output:
        r1="trimmed/{sample}_1.trimmed.fastq.gz",
        r2="trimmed/{sample}_2.trimmed.fastq.gz",
        html="trimmed/{sample}_fastp.html",
        json="trimmed/{sample}_fastp.json"
    conda:
        "envs/fastp.yaml"
    shell:
        """
        mkdir -p trimmed
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --html {output.html} \
            --json {output.json}
        """

rule fastqc_trimmed:
    input:
        "trimmed/{sample}_{read}.trimmed.fastq.gz"
    output:
        html="fastqc_trimmed/{sample}_{read}.trimmed_fastqc.html",
        zip="fastqc_trimmed/{sample}_{read}.trimmed_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p fastqc_trimmed
        fastqc {input} --outdir fastqc_trimmed
        """
rule genome_index:
    input:
        fasta="../../datasets/reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
        gtf="../../datasets/reference/Saccharomyces_cerevisiae.R64-1-1.115.gtf"
    output:
        "star_index/Genome"
    conda:
        "envs/star.yaml"
    shell:
        """
        mkdir -p star_index

        STAR --runMode genomeGenerate \
             --runThreadN 1 \
             --genomeDir star_index \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --genomeSAindexNbases 10
        """
rule align:
    input:
        r1="trimmed/{sample}_1.trimmed.fastq.gz",
        r2="trimmed/{sample}_2.trimmed.fastq.gz",
        index="star_index/Genome"
    output:
        bam="aligned/{sample}_Aligned.sortedByCoord.out.bam"
    conda:
        "envs/star.yaml"
    threads: 1
    shell:
        """
        mkdir -p aligned

        STAR \
            --runThreadN {threads} \
            --genomeDir star_index \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix aligned/{wildcards.sample}_
        """
# Index coordinate-sorted BAM files for efficient genomic access
rule index_bam:
    input:
        bam="aligned/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bai="aligned/{sample}_Aligned.sortedByCoord.out.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam} {output.bai}
        """
# Gene quantification: Count paired-end fragments assigned to genes
rule count_genes:
    input:
        bam="aligned/{sample}_Aligned.sortedByCoord.out.bam",
        bai="aligned/{sample}_Aligned.sortedByCoord.out.bam.bai",
        gtf="../../datasets/reference/Saccharomyces_cerevisiae.R64-1-1.115.gtf"
    output:
        counts="counts/{sample}_gene_counts.txt"
    conda:
        "envs/subread.yaml"
    shell:
        """
        mkdir -p counts

        featureCounts \
            -T 1 \
            -p \
            --countReadPairs \
            -t exon \
            -g gene_id \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam}
        """

rule combine_counts:
    input:
        expand(
            "counts/{sample}_gene_counts.txt",
            sample=SAMPLES
        )
    output:
        "counts/gene_count_matrix.tsv"
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/combine_counts.py"
