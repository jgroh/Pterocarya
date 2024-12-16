import pandas as pd

# local paths to conda environments
SAMTOOLS='/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml'
PIXY="/home/jgroh/heterodichogamy/conda_envs/pixy.yaml"
BEDTOOLS="/home/jgroh/heterodichogamy/conda_envs/bedtools.yaml"

SRA_PAIRED_META=pd.read_table("sra_paired_samples.txt", header = None, names = ['Run', 'Note'])
SRA_PAIRED_SAMPLES=SRA_PAIRED_META['Run'].tolist()

ALL_PTERO_META=pd.read_table("all_Pterocarya_WGS_phenotypes.txt", header = None, names = ['ID', 'type'])
ALL_PTERO_SAMPLES=ALL_PTERO_META['ID'].tolist()

PSTE_META=pd.read_table("Pste_WGS_phenotypes.txt", header = None, names = ['ID', 'type'])
PSTE_SAMPLES=PSTE_META['ID'].tolist()

PSTE_UCD1470_HAP1='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1.fasta'
PSTE_UCD1470_HAP2='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta'
PSTE_UCD1470_HAP1_LIFTOFF='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1_liftoff.gff'
PSTE_UCD1470_HAP2_LIFTOFF='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2_liftoff.gff'
PSTE_UCD1470_HAP1_BWT='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1.fasta.bwt'
PSTE_UCD1470_HAP2_BWT='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta.bwt'

PSTE_UCD1470_HAP1_FAI='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1.fasta.fai'
PSTE_UCD1470_HAP2_FAI='/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta.fai'

HAP1_FAI_TBL = pd.read_table(PSTE_UCD1470_HAP1_FAI, header = None, names = ['chr','w','x','y','z'])
HAP2_FAI_TBL = pd.read_table(PSTE_UCD1470_HAP2_FAI, header = None, names = ['chr','w','x','y','z'])
HAP1_CHRS = HAP1_FAI_TBL['chr'].tolist()
HAP2_CHRS = HAP2_FAI_TBL['chr'].tolist()


IMNU_FASTA='/home/jgroh/Pterocarya/genome_assemblies/IMNU/Pterocarya_stenoptera_genome.fasta'
IMNU_GFF='/home/jgroh/Pterocarya/genome_assemblies/IMNU/Pterocarya_stenoptera_gene.gff3'
IMNU_BWT='/home/jgroh/Pterocarya/genome_assemblies/IMNU/Pterocarya_stenoptera_genome.fasta.bwt'
IMNU_FAI='/home/jgroh/Pterocarya/genome_assemblies/IMNU/Pterocarya_stenoptera_genome.fasta.fai'
IMNU_FAI_TBL = pd.read_table(IMNU_FAI, header = None, names = ['chr','w','x','y','z'])
IMNU_CHRS = IMNU_FAI_TBL['chr'].tolist()

BNU_FASTA='/home/jgroh/Pterocarya/genome_assemblies/BNU/PST_2.0.assembly.fasta'
BNU_BWT='/home/jgroh/Pterocarya/genome_assemblies/BNU/PST_2.0.assembly.fasta.bwt'
BNU_FAI='/home/jgroh/Pterocarya/genome_assemblies/BNU/PST_2.0.assembly.fasta.fai'
BNU_FAI_TBL = pd.read_table(BNU_FAI, header = None, names = ['chr','w','x','y','z'])
BNU_CHRS = BNU_FAI_TBL['chr'].tolist()

PFRA_SAMPLES = ['PFRA_W14_2', 'PFRA_W18_2']
PMAC_SAMPLES = ['PMAC_SO4_S63']
PRHO_SAMPLES = ['PRHO_SO1_S61', 'PRHO_SO6_S64']

ISOSEQ_SAMPLES_STEVENS = ['SRR25617144', 'SRR25617145', 'SRR25617146', 'SRR25617147', 'SRR25617149', 'SRR25617150']
ISOSEQ_SAMPLES_ZHANG =  ['SRR26994970']

IMNU_RNASEQ_SAMPLES = ['SRR26994971', 'SRR26994972', 'SRR26994973', 'SRR26994974']

RNASEQ_MYSAMPLES = ['DV_136', 'DV_138_F','DV_138_M', 'DV_145', 'DV_146.5', 'DV_149_F','DV_149_M', 'DV_150', 'PSTE_UCD1', 'PSTE_WS_11.01', 'PSTE_WS_2.08','PSTE_WS_2.10']

ruleorder: align_to_IMNU > align_to_IMNU_SRA
ruleorder: align_to_HAP1 > align_to_HAP1_SRA
ruleorder: align_to_HAP2 > align_to_HAP2_SRA
ruleorder: index_bam_general > index_bam
ruleorder: index_bam_general > index_bam_IsoSeq

rule all:
  input:
    expand("results/coverage_hap1_RNAseq/mydata/{sample}_Gloc_RNAseq.txt.gz", sample = RNASEQ_MYSAMPLES),
    expand("results/coverage_hap2_RNAseq/mydata/{sample}_Gloc_RNAseq.txt.gz", sample = RNASEQ_MYSAMPLES),
    expand("results/coverage_hap1_RNAseq/published_data/{sample}_Gloc_RNAseq.txt.gz", sample = IMNU_RNASEQ_SAMPLES),
    expand("results/coverage_hap2_RNAseq/published_data/{sample}_Gloc_RNAseq.txt.gz", sample = IMNU_RNASEQ_SAMPLES),
#    "denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1_liftoff.gff",
#    "denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2_liftoff.gff",
    #"denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1_liftoff.gff", 
#    "denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2_liftoff.gff", 
    expand("output/PSTE_UCD1470_HAP{N}/Pste.assoc.txt", N = ['1', '2']), 
    expand("calls/PSTE_UCD1470_HAP{N}/Pste.bed", N=['1','2']),
#    expand("calls/PSTE_UCD1470_HAP{N}/biallelic_SNPs_filtered.vcf.gz", N=['1','2']),
#    "calls/PSTE_UCD1470_HAP2/allsites_unfiltered.vcf.gz",
#    "calls/PSTE_UCD1470_HAP1/allsites_unfiltered.vcf.gz",
    expand("results/coverage_UCD1470_HAP{N}/{sample}_Gloc.txt.gz", N=['1','2'], sample= PSTE_SAMPLES),
    expand("results/coverage_UCD1470_HAP1/{sample}_l{n}.txt.gz", sample= PSTE_SAMPLES, n = ['1','2']),
#    expand("alignment_files/PSTE_UCD1470_HAP{N}/{sample}.cram", N=['1','2'],sample = PSTE_SAMPLES),
    expand("alignment_files/{hap}_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam", hap = ['hap1', 'hap2'], sample = RNASEQ_MYSAMPLES), 
    expand("alignment_files/{hap}_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam", hap = ['hap1', 'hap2'], sample = IMNU_RNASEQ_SAMPLES),
#    "results/pixy/1kb_pi.txt",
#    "results/pixy/UCD2_SNPs_pi.txt",
#    "calls/IMNU/UCD2_heterozygous_SNPs.bed",
    expand("results/coverage_IMNU/{sample}_chr{N}.txt.gz", sample=PSTE_SAMPLES, N = ['12','14','16']),
#j    expand("alignment_files/IMNU/{sample}_focal_region.bam", sample = PFRA_SAMPLES+PMAC_SAMPLES+PRHO_SAMPLES),
    expand("sra_fasta_paired/{sample}_{N}.fasta.gz", sample=SRA_PAIRED_SAMPLES, N = ['1','2']),
    "calls/IMNU/allsites_filtered_Chr11.vcf.gz",
#    "HiFi_data_IMNU/SRR26994977.fasta",
    "calls/IMNU/allsites_filtered_Chr11.vcf.gz",
    BNU_FASTA+".mod.EDTA.TEanno.gff3",
    IMNU_FASTA+".mod.EDTA.TEanno.gff3",
    expand("results/coverage_IMNU/{sample}_Gloc.txt.gz", sample = ALL_PTERO_SAMPLES),
    expand("results/coverage_IMNU_RNAseq/{sample}_Gloc.txt.gz", sample = IMNU_RNASEQ_SAMPLES),
    expand("results/coverage_IsoSeq_hap1/{sample}_Gloc.txt.gz", sample = ISOSEQ_SAMPLES_STEVENS),
    expand("results/coverage_IsoSeq_hap1/{sample}_Gloc.txt.gz", sample = ISOSEQ_SAMPLES_ZHANG),
    expand("results/coverage_IsoSeq_hap2/{sample}_Gloc.txt.gz", sample = ISOSEQ_SAMPLES_STEVENS),
    expand("results/coverage_IsoSeq_hap2/{sample}_Gloc.txt.gz", sample = ISOSEQ_SAMPLES_ZHANG),
#    expand("calls/{genome}/Pste.bed", genome = ['IMNU', 'BNU']),   # final version complete, remain commented
#    expand("results/coverage_IMNU/{sample}_Gloc.txt.gz", sample = ALL_PTERO_SAMPLES)


rule index_assemblies:
  input:
    "{assembly}.fasta"
  output:
    "{assembly}.fasta.bwt"
  conda:
    SAMTOOLS
  shell:
    "bwa index {input}"


rule align_to_HAP1:
  input:
    fa=PSTE_UCD1470_HAP1,
    bwt=PSTE_UCD1470_HAP1_BWT,
    R1="fastq/{sample}_R1.fastq.gz",
    R2="fastq/{sample}_R2.fastq.gz",
  output:
    "alignment_files/PSTE_UCD1470_HAP1/{sample}.cram"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    runtime=24*60,
    mem_mb=30000
  threads:
    24
  shell:
    """
    mkdir -p tmp_sorted_reads/PSTE_UCD1470_HAP1/{wildcards.sample} &&
    bwa mem -t {threads} -M \
        -R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' \
        {input.fa} {input.R1} {input.R2} | \
        samtools sort -O bam -l 0 -T tmp_sorted_reads/IMNU/{wildcards.sample} - | \
        samtools view -T {input.fa} -C -o {output} -
    rm -r tmp_sorted_reads/PSTE_UCD1470_HAP1/{wildcards.sample}
    """

rule align_to_HAP1_SRA:
  input:
    fa=PSTE_UCD1470_HAP1,
    bwt=PSTE_UCD1470_HAP1_BWT,
    R1="sra_fasta_paired/{sample}_1.fasta.gz",
    R2="sra_fasta_paired/{sample}_1.fasta.gz",
  output:
    "alignment_files/PSTE_UCD1470_HAP1/{sample}.cram"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    runtime=24*60,
    mem_mb=30000
  threads:
    24
  wildcard_constraints: 
    sample="(?!SRR26994971$)(?!SRR26994972$)(?!SRR26994973$)(?!SRR26994974$).+" # Exclude the IMNU RNAseq samples from this list
  shell:
    """
    mkdir -p tmp_sorted_reads/HAP1/{wildcards.sample} &&
    bwa mem -t {threads} -M \
        -R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' \
        {input.fa} {input.R1} {input.R2} | \
        samtools sort -O bam -l 0 -T tmp_sorted_reads/HAP1/{wildcards.sample} - | \
        samtools view -T {input.fa} -C -o {output} -
    rm -r tmp_sorted_reads/HAP1/{wildcards.sample}
    """

rule align_to_HAP2:
  input:
    fa=PSTE_UCD1470_HAP2,
    bwt=PSTE_UCD1470_HAP2_BWT,
    R1="fastq/{sample}_R1.fastq.gz",
    R2="fastq/{sample}_R2.fastq.gz",
  output:
    "alignment_files/PSTE_UCD1470_HAP2/{sample}.cram"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    runtime=24*60,
    mem_mb=30000
  threads:
    24
  shell:
    """
    mkdir -p tmp_sorted_reads/PSTE_UCD1470_HAP2/{wildcards.sample} &&
    bwa mem -t {threads} -M \
        -R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' \
        {input.fa} {input.R1} {input.R2} | \
        samtools sort -O bam -l 0 -T tmp_sorted_reads/IMNU/{wildcards.sample} - | \
        samtools view -T {input.fa} -C -o {output} -
    rm -r tmp_sorted_reads/PSTE_UCD1470_HAP2/{wildcards.sample}
    """

rule align_to_HAP2_SRA:
  input:
    fa=PSTE_UCD1470_HAP2,
    bwt=PSTE_UCD1470_HAP2_BWT,
    R1="sra_fasta_paired/{sample}_1.fasta.gz",
    R2="sra_fasta_paired/{sample}_1.fasta.gz",
  output:
    "alignment_files/PSTE_UCD1470_HAP2/{sample}.cram"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    runtime=24*60,
    mem_mb=30000
  threads:
    24
  wildcard_constraints: 
    sample="(?!SRR26994971$)(?!SRR26994972$)(?!SRR26994973$)(?!SRR26994974$).+" # Exclude the IMNU RNAseq samples from this list
  shell:
    """
    mkdir -p tmp_sorted_reads/HAP2/{wildcards.sample} &&
    bwa mem -t {threads} -M \
        -R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' \
        {input.fa} {input.R1} {input.R2} | \
        samtools sort -O bam -l 0 -T tmp_sorted_reads/HAP2/{wildcards.sample} - | \
        samtools view -T {input.fa} -C -o {output} -
    rm -r tmp_sorted_reads/HAP2/{wildcards.sample}
    """






rule bcftools_call_by_chrom_HAP1:
  input:
    fa=PSTE_UCD1470_HAP1,
    fai=PSTE_UCD1470_HAP1+'.fai',
    crams=expand("alignment_files/PSTE_UCD1470_HAP1/{sample}.cram", sample = ALL_PTERO_SAMPLES),
    crai=expand("alignment_files/PSTE_UCD1470_HAP1/{sample}.cram.crai", sample = ALL_PTERO_SAMPLES),
  output:
    temp("calls/PSTE_UCD1470_HAP1/{chr}_allsites.vcf.gz")
  resources:
    runtime = 48*60,
    mem_mb=50000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    bcftools mpileup -Ou -f {input.fa} -r {wildcards.chr} \
        --annotate "AD,DP,INFO/AD" {input.crams} | bcftools call -m -f GQ,GP -o {output}
    """

rule bcftools_call_by_chrom_HAP2:
  input:
    fa=PSTE_UCD1470_HAP2,
    fai=PSTE_UCD1470_HAP2+'.fai',
    crams=expand("alignment_files/PSTE_UCD1470_HAP2/{sample}.cram", sample = ALL_PTERO_SAMPLES),
    crai=expand("alignment_files/PSTE_UCD1470_HAP2/{sample}.cram.crai", sample = ALL_PTERO_SAMPLES),
  output:
    temp("calls/PSTE_UCD1470_HAP2/{chr}_allsites.vcf.gz")
  resources:
    runtime = 48*60,
    mem_mb=50000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    bcftools mpileup -Ou -f {input.fa} -r {wildcards.chr} \
        --annotate "AD,DP,INFO/AD" {input.crams} | bcftools call -m -f GQ,GP -o {output}
    """

rule bcftools_concat_variant_HAP1:
  input:
    expand("calls/PSTE_UCD1470_HAP1/{chr}_allsites.vcf.gz", chr = HAP1_CHRS)
  output:
    "calls/PSTE_UCD1470_HAP1/allsites_unfiltered.vcf.gz"
  resources:
    runtime = 24*60,
    mem_mb = 30000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools concat {input} -o {output}"

rule bcftools_concat_variant_HAP2:
  input:
    expand("calls/PSTE_UCD1470_HAP2/{chr}_allsites.vcf.gz", chr = HAP2_CHRS)
  output:
    "calls/PSTE_UCD1470_HAP2/allsites_unfiltered.vcf.gz"
  resources:
    runtime = 24*60,
    mem_mb = 30000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools concat {input} -o {output}"


rule filter_variants_allsites_1:
  input:
    vcf="calls/PSTE_UCD1470_HAP{N}/allsites_unfiltered.vcf.gz",
    fa="denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap{N}/Pste_UCD1470_hap{N}.fasta",
  output:
    temp("calls/PSTE_UCD1470_HAP{N}/allsites.norm.bcf")
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools norm -f {input.fa} {input.vcf} -Ob -o {output}"

rule filter_variants_allsites_2:
  input:
    bcf="calls/PSTE_UCD1470_HAP{N}/allsites.norm.bcf",
  output:
    "calls/PSTE_UCD1470_HAP{N}/allsites.norm.flt-indels.bcf"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools filter --IndelGap 5 {input.bcf} -Ob -o {output}"


rule filter_variants_biallelic:
  input:
    "calls/PSTE_UCD1470_HAP{N}/allsites_unfiltered.vcf.gz"
  output:
    "calls/PSTE_UCD1470_HAP{N}/biallelic_SNPs_filtered.vcf.gz"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools view -i 'QUAL>=30 & FMT/GQ>30 & FMT/DP>10 & FMT/DP<150' -m2 -M2 -v snps -o {output} {input}"




rule coverage_PSTE_UCD1470_HAP1:
  input:
    cram="alignment_files/PSTE_UCD1470_HAP1/{sample}.cram",
    crai="alignment_files/PSTE_UCD1470_HAP1/{sample}.cram.crai",
  output:
    focal="results/coverage_UCD1470_HAP1/{sample}_Gloc.txt.gz",
    norm="results/coverage_UCD1470_HAP1/{sample}_norm.txt",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.cram} -r Chr11:3780925-4061815 | gzip > {output.focal}
    samtools depth -a {input.cram} -r Chr01 | awk '{{ sum += $3 }} END {{ mean = sum/ NR; print mean }}' > {output.norm}
    """

rule coverage_PSTE_UCD1470_HAP1_peaks:
  input:
    cram="alignment_files/PSTE_UCD1470_HAP1/{sample}.cram",
    crai="alignment_files/PSTE_UCD1470_HAP1/{sample}.cram.crai",
  output:
    l1="results/coverage_UCD1470_HAP1/{sample}_l1.txt.gz",
    l2="results/coverage_UCD1470_HAP1/{sample}_l2.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.cram} -r Chr04:12994795-12996810 | gzip > {output.l1}
    samtools depth -a {input.cram} -r Chr04:16459732-16461732 | gzip > {output.l2}
    """


rule coverage_PSTE_UCD1470_HAP2:
  input:
    cram="alignment_files/PSTE_UCD1470_HAP2/{sample}.cram",
    crai="alignment_files/PSTE_UCD1470_HAP2/{sample}.cram.crai",
  output:
    focal="results/coverage_UCD1470_HAP2/{sample}_Gloc.txt.gz",
    norm="results/coverage_UCD1470_HAP2/{sample}_norm.txt"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.cram} -r Chr11:4266881-4534617 | gzip > {output.focal}
    samtools depth -a {input.cram} -r Chr01 | awk '{{ sum += $3 }} END {{ mean = sum/ NR; print mean }}' > {output.norm}
    """

rule RNAseq_coverage_mydata_PSTE_UCD1470_HAP1:
  input:
    bam="alignment_files/hap1_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam",
    bai="alignment_files/hap1_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam.bai",
  output:
    focal="results/coverage_hap1_RNAseq/mydata/{sample}_Gloc_RNAseq.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.bam} -r Chr11:3780925-4061815 | gzip > {output.focal}
    """

rule RNAseq_coverage_mydata_PSTE_UCD1470_HAP2:
  input:
    bam="alignment_files/hap2_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam",
    bai="alignment_files/hap2_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam.bai",
  output:
    focal="results/coverage_hap2_RNAseq/mydata/{sample}_Gloc_RNAseq.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.bam} -r Chr11:4266881-4534617 | gzip > {output.focal}
    """

rule RNAseq_coverage_published_data_PSTE_UCD1470_HAP1:
  input:
    bam="alignment_files/hap1_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam",
    bai="alignment_files/hap1_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam.bai",
  output:
    focal="results/coverage_hap1_RNAseq/published_data/{sample}_Gloc_RNAseq.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.bam} -r Chr11:3780925-4061815 | gzip > {output.focal}
    """

rule RNAseq_coverage_published_data_PSTE_UCD1470_HAP2:
  input:
    bam="alignment_files/hap2_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam",
    bai="alignment_files/hap2_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam.bai",
  output:
    focal="results/coverage_hap2_RNAseq/published_data/{sample}_Gloc_RNAseq.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.bam} -r Chr11:4266881-4534617 | gzip > {output.focal}
    """

# ------ analyses with previous assemblies -----
rule download_HiFi_reads_IMNU:
  input:
  output:
    "HiFi_data_IMNU/SRR26994977.fasta"
  resources:
    mem_mb=5000,
    runtime=60
  threads:
    8
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/sra-tools.yaml"
  shell:
    "fasterq-dump --outdir HiFi_data_IMNU --fasta --threads {threads} --skip-technical --verbose SRR26994977"

rule download_sra_paired:
  input:
  output:
    r1=temp("sra_fasta_paired/{sample}_1.fasta"),
    r2=temp("sra_fasta_paired/{sample}_2.fasta"),
  resources:
    mem_mb=40000,
    runtime=1440
  threads:
    8
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/sra-tools.yaml"
  shell:
    "fasterq-dump --outdir sra_fasta_paired --split-files --fasta --threads {threads} --skip-technical --verbose {wildcards.sample}"


rule download_sra_single:
  input:
  output:
    r1="sra_fasta_single/{sample}.fasta",
  resources:
    mem_mb=40000,
    runtime=1440
  threads:
    8
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/sra-tools.yaml"
  shell:
    "fasterq-dump --outdir sra_fasta_single --fasta --threads {threads} --skip-technical --verbose {wildcards.sample}"

#rule gzip_sra_fasta_paired:
#  input:
#    "sra_fasta_paired/{sample}_{N}.fasta"
#  output:
#    "sra_fasta_paired/{sample}_{N}.fasta.gz"
#  shell:
#    "gzip {input}"

rule gzip_sra_fasta_single:
  input:
    "sra_fasta_single/{sample}.fasta"
  output:
    "sra_fasta_single/{sample}.fasta.gz"
  shell:
    "gzip {input}"


rule trim_sra_fasta_paired:
  input:
    r1="sra_fasta_paired/{sample}_1.fasta",
    r2="sra_fasta_paired/{sample}_2.fasta",
  output:
    r1="sra_fasta_paired/{sample}-trimmed-pair1.fastq.gz",
    r2="sra_fasta_paired/{sample}-trimmed-pair2.fastq.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/skewer.yaml"
  shell:
    """
    skewer -m any --threads {threads} --compress --output sra_fasta_paired/{wildcards.sample} {input.r1} {input.r2}
    """

rule align_IsoSeq_Stevens_etal_2018_data_to_hap1:
  input:
    fasta="/home/jgroh/heterodichogamy/gene_expression/IsoSeq/fasta/{sample}.fasta",
    ref=PSTE_UCD1470_HAP1
  output:
    "alignment_files/IsoSeq_hap1/{sample}.bam"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  threads:
    10
  resources:
    mem_mb=50000,
    runtime=1440,
  wildcard_constraints:
    sample="(?!SRR26994970$).+"
  shell:
    "minimap2 -ax splice -uf -C5 {input.ref} {input.fasta} | samtools sort -O bam -o {output} -"

rule align_IsoSeq_Stevens_etal_2018_data_to_hap2:
  input:
    fasta="/home/jgroh/heterodichogamy/gene_expression/IsoSeq/fasta/{sample}.fasta",
    ref=PSTE_UCD1470_HAP2
  output:
    "alignment_files/IsoSeq_hap2/{sample}.bam"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  threads:
    10
  resources:
    mem_mb=50000,
    runtime=1440,
  wildcard_constraints:
    sample="(?!SRR26994970$).+"
  shell:
    "minimap2 -ax splice -uf -C5 {input.ref} {input.fasta} | samtools sort -O bam -o {output} -"

rule align_IsoSeq_Zhang_etal_2024_data_hap1:
  input:
    fasta="sra_fasta_single/SRR26994970.fasta",
    ref=PSTE_UCD1470_HAP1
  output:
    "alignment_files/IsoSeq_hap1/SRR26994970.bam"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    mem_mb=50000,
    runtime=1440
  threads:
    10
  shell:
    "minimap2 -ax splice -uf -C5 {input.ref} {input.fasta} | samtools sort -O bam -o {output} -"

rule align_IsoSeq_Zhang_etal_2024_data_hap2:
  input:
    fasta="sra_fasta_single/SRR26994970.fasta",
    ref=PSTE_UCD1470_HAP2
  output:
    "alignment_files/IsoSeq_hap2/SRR26994970.bam"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    mem_mb=50000,
    runtime=1440
  threads:
    10
  shell:
    "minimap2 -ax splice -uf -C5 {input.ref} {input.fasta} | samtools sort -O bam -o {output} -"

rule align_to_IMNU:
  input:
    fa=IMNU_FASTA,
    bwt=IMNU_BWT,
    R1="fastq/{sample}_R1.fastq.gz",
    R2="fastq/{sample}_R2.fastq.gz",
  output:
    "alignment_files/IMNU/{sample}.cram"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    runtime=24*60,
    mem_mb=30000
  threads:
    24
  shell:
    """
    mkdir -p tmp_sorted_reads/IMNU/{wildcards.sample} &&
    bwa mem -t {threads} -M \
        -R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' \
        {input.fa} {input.R1} {input.R2} | \
        samtools sort -O bam -l 0 -T tmp_sorted_reads/IMNU/{wildcards.sample} - | \
        samtools view -T {input.fa} -C -o {output} -
    rm -r tmp_sorted_reads/IMNU/{wildcards.sample}
    """

rule align_to_IMNU_SRA:
  input:
    fa=IMNU_FASTA,
    bwt=IMNU_BWT,
    R1="sra_fasta_paired/{sample}_1.fasta.gz",
    R2="sra_fasta_paired/{sample}_1.fasta.gz",
  output:
    "alignment_files/IMNU/{sample}.cram"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    runtime=24*60,
    mem_mb=30000
  threads:
    24
  wildcard_constraints: 
    sample="(?!SRR26994971$)(?!SRR26994972$)(?!SRR26994973$)(?!SRR26994974$).+" # Exclude the IMNU RNAseq samples from this list
  shell:
    """
    mkdir -p tmp_sorted_reads/IMNU/{wildcards.sample} &&
    bwa mem -t {threads} -M \
        -R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' \
        {input.fa} {input.R1} {input.R2} | \
        samtools sort -O bam -l 0 -T tmp_sorted_reads/IMNU/{wildcards.sample} - | \
        samtools view -T {input.fa} -C -o {output} -
    rm -r tmp_sorted_reads/IMNU/{wildcards.sample}
    """

rule align_to_BNU:
  input:
    fa=BNU_FASTA,
    bwt=BNU_BWT,
    R1="fastq/{sample}_R1.fastq.gz",
    R2="fastq/{sample}_R2.fastq.gz",
  output:
    "alignment_files/BNU/{sample}.cram"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  resources:
    runtime=24*60,
    mem_mb=30000
  threads:
    24
  shell:
    """
    mkdir -p tmp_sorted_reads/BNU/{wildcards.sample} &&
    bwa mem -t {threads} -M \
        -R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' \
        {input.fa} {input.R1} {input.R2} | \
        samtools sort -O bam -l 0 -T tmp_sorted_reads/BNU/{wildcards.sample} - | \
        samtools view --require-flags 2 --min-MQ 40 -T {input.fa} -C -o {output} -
    rm -r tmp_sorted_reads/BNU/{wildcards.sample}
    """

rule index_bam_general:
  input:
    "{file}.bam"
  output:
    "{file}.bam.bai"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "samtools index {input}"

rule index_bam:
  input:
    "alignment_files/IMNU/{sample}.bam"
  output:
    "alignment_files/IMNU/{sample}.bam.bai"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "samtools index {input}"

rule index_bam_RNAseq:
  input:
    "alignment_files/IMNU_RNAseq/{sample}.bam"
  output:
    "alignment_files/IMNU_RNAseq/{sample}.bam.bai"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "samtools index {input}"

rule index_bam_IsoSeq:
  input:
    "alignment_files/IsoSeq_{hap}/{sample}.bam"
  output:
    "alignment_files/IsoSeq_{hap}/{sample}.bam.bai"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "samtools index {input}"

rule index_cram:
  input:
    "alignment_files/{genome}/{sample}.cram"
  output:
    "alignment_files/{genome}/{sample}.cram.crai"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "samtools index {input}"

rule bcftools_call_by_chrom_IMNU:
  input:
    fa=IMNU_FASTA,
    fai=IMNU_FAI,
    crams=expand("alignment_files/IMNU/{sample}.cram", sample = ALL_PTERO_SAMPLES),
    crai=expand("alignment_files/IMNU/{sample}.cram.crai", sample = ALL_PTERO_SAMPLES),
  output:
    temp("calls/IMNU/{chr}_allsites.vcf.gz")
  resources:
    runtime = 48*60,
    mem_mb=50000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    bcftools mpileup -Ou -f {input.fa} -r {wildcards.chr} \
        --annotate "AD,DP,INFO/AD" {input.crams} | bcftools call -m -f GQ,GP -o {output}
    """

rule bcftools_call_by_chrom_BNU:
  input:
    fa=BNU_FASTA,
    fai=BNU_FAI,
    crams=expand("alignment_files/BNU/{sample}.cram", sample = ALL_PTERO_SAMPLES),
    crai=expand("alignment_files/BNU/{sample}.cram.crai", sample = ALL_PTERO_SAMPLES),
  output:
    temp("calls/BNU/{chr}_allsites.vcf.gz")
  resources:
    runtime = 48*60,
    mem_mb=50000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    bcftools mpileup -Ou -f {input.fa} -r {wildcards.chr} \
        --annotate "AD,DP,INFO/AD" {input.crams} | bcftools call -m -f GQ,GP -o {output}
    """

rule bcftools_concat_variant_IMNU:
  input:
    expand("calls/IMNU/{chr}_allsites.vcf.gz", chr = IMNU_CHRS)
  output:
    "calls/IMNU/allsites_unfiltered.vcf.gz"
  resources:
    runtime = 24*60,
    mem_mb = 30000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools concat {input} -o {output}"


rule bcftools_concat_variant_BNU:
  input:
    expand("calls/BNU/{chr}_allsites.vcf.gz", chr = BNU_CHRS)
  output:
    "calls/BNU/allsites_unfiltered.vcf.gz"
  resources:
    runtime = 24*60,
    mem_mb = 30000
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools concat {input} -o {output}"

rule subset_variant:
  input:
    "calls/{genome}/allsites_unfiltered.vcf.gz"
  output:
    "calls/{genome}/variant_unfiltered.vcf.gz"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/vcftools.yaml"
  shell:
    "vcftools --gzvcf {input} --mac 1 --recode --stdout | bgzip -c > {output}"

rule filter_variant:
  input:
    "calls/{genome}/variant_unfiltered.vcf.gz"
  output:
    "calls/{genome}/variant_filtered.vcf.gz"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/vcftools.yaml"
  shell:
    """
    vcftools --gzvcf {input} \
        --remove-indels \
        --mac 2 \
        --max-missing 0.8 \
        --minQ 30 \
        --minGQ 20 \
        --minDP 10 \
        --maxDP 200 \
        --min-alleles 2 \
        --max-alleles 2 \
        --recode --stdout | bgzip -c > {output}
    """

rule subset_and_filter_invariant_Chr11:
  input:
    "calls/{genome}/allsites_unfiltered.vcf.gz"
  output:
    "calls/{genome}/invariant_filtered_Chr11.vcf.gz"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/vcftools.yaml"
  shell:
    """
    vcftools --gzvcf {input} \
        --chr Chr11 \
        --max-maf 0 \
        --remove-indels \
        --max-missing 0.8 \
        --minGQ 20 \
        --recode --stdout | bgzip -c > {output}
    """

rule index_filtered_vcfs:
  input:
    "calls/{genome}/{sites}_filtered.vcf.gz"
  output:
    "calls/{genome}/{sites}_filtered.vcf.gz.csi"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools index {input}"

rule index_filtered_vcfs_Chr11:
  input:
    "calls/{genome}/{sites}_filtered_Chr11.vcf.gz"
  output:
    "calls/{genome}/{sites}_filtered_Chr11.vcf.gz.csi"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools index {input}"

rule subset_and_filter_variants_on_Chr11:
  input:
    "calls/{genome}/allsites_unfiltered.vcf.gz"
  output:
    "calls/{genome}/variant_filtered_Chr11.vcf.gz"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    vcftools --gzvcf {input} \
        --chr Chr11 \
        --remove-indels \
        --mac 2 \
        --minGQ 20 \
        --max-alleles 2 \
        --min-alleles 2 \
        --max-missing 0.8 \
        --recode --stdout | bgzip -c > {output}
    """

rule concat_vcfs_Chr11:
  input:
    variant="calls/{genome}/variant_filtered_Chr11.vcf.gz",
    variant_csi="calls/{genome}/variant_filtered_Chr11.vcf.gz.csi",
    invariant="calls/{genome}/invariant_filtered_Chr11.vcf.gz",
    invariant_csi="calls/{genome}/invariant_filtered_Chr11.vcf.gz.csi",
  output:
    "calls/{genome}/allsites_filtered_Chr11.vcf.gz"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "bcftools concat --allow-overlaps {input.variant} {input.invariant} -O z -o {output}"

rule index_vcf:
  input:
    "calls/{genome}/{sites}_filtered.vcf.gz"
  output:
    "calls/{genome}/{sites}_filtered.vcf.gz.tbi"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "tabix {input}"


rule vcf_to_plink:
  input:
    "calls/{genome}/biallelic_SNPs_filtered.vcf.gz"
  output:
    bed="calls/{genome}/Pste.bed",
    bim="calls/{genome}/Pste.bim",
    fam="calls/{genome}/Pste.fam",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/plink_1.90.yaml"
  shell:
    "plink --vcf {input} --remove no_pheno_samples_plink.txt --const-fid 0 --allow-extra-chr --make-bed --out calls/{wildcards.genome}/Pste"

rule gemma_relatedness_matrix: # only run after manual edit of fam
  input:
    bed="calls/{genome}/Pste.bed",
    bim="calls/{genome}/Pste.bim",
    fam="calls/{genome}/Pste.fam", 
  output:
    "output/{genome}/Pste.cXX.txt"
  resources:
    mem_mb=20000
  shell:
    "gemma -bfile calls/{wildcards.genome}/Pste -gk 1 -o {wildcards.genome}/Pste"


rule gemma_lmm:
  input:
    bed="calls/{genome}/Pste.bed",
    bim="calls/{genome}/Pste.bim",
    fam="calls/{genome}/Pste.fam", # AFTER MANUAL EDIT
    mat="output/{genome}/Pste.cXX.txt"
  output:
    "output/{genome}/Pste.assoc.txt"
  resources:
    runtime=60,
    mem_mb=50000
  shell:
    "gemma -bfile calls/{wildcards.genome}/Pste -k {input.mat} -lmm 2 -miss 0.1 -maf 0.1 -o {wildcards.genome}/Pste"

rule coverage_IMNU:
  input:
    cram="alignment_files/IMNU/{sample}.cram",
    crai="alignment_files/IMNU/{sample}.cram.crai"
  output:
    focal="results/coverage_IMNU/{sample}_Gloc.txt.gz",
    norm="results/coverage_IMNU/{sample}_norm.txt"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.cram} -r Chr11:4300000-4750000 | gzip > {output.focal}
    samtools depth -a {input.cram} -r Chr10 | awk '{{ sum += $3 }} END {{ mean = sum/ NR; print mean }}' > {output.norm}
    """

rule coverage_IMNU_other_peaks:
  input:
    cram="alignment_files/IMNU/{sample}.cram",
    crai="alignment_files/IMNU/{sample}.cram.crai"
  output:
    chr12="results/coverage_IMNU/{sample}_chr12.txt.gz",
    chr14="results/coverage_IMNU/{sample}_chr14.txt.gz",
    chr16="results/coverage_IMNU/{sample}_chr16.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.cram} -r Chr12:1430000-1440000 | gzip > {output.chr12}
    samtools depth -a {input.cram} -r Chr14:13600000-13650000 | gzip > {output.chr14}
    samtools depth -a {input.cram} -r Chr16:17150000-17200000 | gzip > {output.chr16}
    """

rule coverage_BNU:
  input:
    cram="alignment_files/BNU/{sample}.cram",
    crai="alignment_files/BNU/{sample}.cram.crai"
  output:
    focal="results/coverage_BNU/{sample}_Gloc.txt.gz",
    norm="results/coverage_BNU/{sample}_norm.txt"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.cram} -r CTG_841:1500000-1967201 | gzip > {output.focal}
    samtools depth -a {input.cram} -r CTG_1 | awk '{{ sum += $3 }} END {{ mean = sum/ NR; print mean }}' > {output.norm}
    """

rule coverage_IsoSeq_hap1:
  input:
    bam="alignment_files/IsoSeq_hap1/{sample}.bam",
    bai="alignment_files/IsoSeq_hap1/{sample}.bam.bai"
  output:
    focal="results/coverage_IsoSeq_hap1/{sample}_Gloc.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.bam} -r Chr11:3780925-4061815 | gzip > {output.focal}
    """

rule coverage_IsoSeq_hap2:
  input:
    bam="alignment_files/IsoSeq_hap2/{sample}.bam",
    bai="alignment_files/IsoSeq_hap2/{sample}.bam.bai"
  output:
    focal="results/coverage_IsoSeq_hap2/{sample}_Gloc.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.bam} -r Chr11:4266881-4534617 | gzip > {output.focal}
    """

rule coverage_IMNU_RNAseq:
  input:
    bam="alignment_files/IMNU_RNAseq/{sample}_filtered.bam",
    bai="alignment_files/IMNU_RNAseq/{sample}_filtered.bam.bai"
  output:
    focal="results/coverage_IMNU_RNAseq/{sample}_Gloc.txt.gz",
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    """
    samtools depth -a {input.bam} -r Chr11:4300000-4750000 | gzip > {output.focal}
    """

# Annotate TEs in the IMNU assembly
rule EDTA_IMNU:
  input:
    IMNU_FASTA
  output:
    IMNU_FASTA+".mod.EDTA.TEanno.gff3"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/EDTA.yml"
  threads:
    32
  resources:
    runtime=72*60,
    mem_mb=50000
  shell:
    """
    perl /home/jgroh/heterodichogamy/EDTA/EDTA.pl \
        --genome /home/jgroh/Pterocarya/genome_assemblies/IMNU/Pterocarya_stenoptera_genome.fasta \
        --anno 1 \
        --sensitive 1 \
        --overwrite 1 \
        --threads {threads} \
        --u 1.5e-9
    """


# Annotate TEs in the BNU assembly
rule EDTA_BNU:
  input:
    BNU_FASTA
  output:
    BNU_FASTA+".mod.EDTA.TEanno.gff3"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/EDTA.yml"
  threads:
    32
  resources:
    runtime=72*60,
    mem_mb=50000
  shell:
    """
    perl /home/jgroh/heterodichogamy/EDTA/EDTA.pl \
        --genome /home/jgroh/Pterocarya/genome_assemblies/BNU/PST_2.0.assembly.fasta \
        --anno 1 \
        --sensitive 1 \
        --overwrite 1 \
        --threads {threads} \
        --u 1.5e-9
    """


rule liftoff_IMNU_to_hap1:
  input:
    ref=IMNU_FASTA,
    gff=IMNU_GFF,
    target=PSTE_UCD1470_HAP1
  output:
    "denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1_liftoff.gff"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/liftoff.yaml"
  resources:
    mem_mb=40000,
    runtime=12*60
  shell:
    """
    liftoff -g {input.gff} \
        -dir "/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/" \
        -o {output} \
        -u "/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/liftoff_unmapped.txt" \
        {input.target} {input.ref}
    """

rule liftoff_IMNU_to_hap2:
  input:
    ref=IMNU_FASTA,
    gff=IMNU_GFF,
    target=PSTE_UCD1470_HAP2
  output:
    "denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2_liftoff.gff"
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/liftoff.yaml"
  resources:
    mem_mb=40000,
    runtime=12*60
  shell:
    """
    liftoff -g {input.gff} \
        -dir "/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/" \
        -o {output} \
        -u "/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/liftoff_unmapped.txt" \
        {input.target} {input.ref}
    """




# This step is a prerequisite to run STAR. 
rule STAR_index_hap1:
  input:
    fa=PSTE_UCD1470_HAP1,
    gff="denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/P.stenoptera_UCD1470_hap1_liftoff.gff"
  output:
    "/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/SA"
  threads:
    20
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  resources:
    mem_mb=50000,
    runtime=60
  shell:
    """
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir /home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1 \
        --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 \
         --outTmpDir hap1_STAR_index \
        --genomeSAindexNbases 13
    """

# This step is a prerequisite to run STAR. 
rule STAR_index_hap2:
  input:
    fa=PSTE_UCD1470_HAP2,
    gff="denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2_liftoff.gff"
  output:
    "/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/SA"
  threads:
    20
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  resources:
    mem_mb=50000,
    runtime=60
  shell:
    """
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir /home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2 \
        --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 \
        --outTmpDir hap2_STAR_index \
        --genomeSAindexNbases 13
    """


# This step is a prerequisite to run STAR. 
rule STAR_index_IMNU:
  input:
    fa=IMNU_FASTA,
    gff=IMNU_GFF
  output:
    "/home/jgroh/Pterocarya/genome_assemblies/IMNU/SA"
  threads:
    20
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  resources:
    mem_mb=50000,
    runtime=60
  shell:
    """
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir /home/jgroh/Pterocarya/genome_assemblies/IMNU \
        --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 \
        --genomeSAindexNbases 13
    """

# Align RNAseq reads from published datasets
rule STAR_align_IMNU_published_data:
  input:
    fa=IMNU_FASTA,
    index="/home/jgroh/Pterocarya/genome_assemblies/IMNU/SA",
    r1="sra_fasta_paired/{sample}_1.fastq.gz",
    r2="sra_fasta_paired/{sample}_2.fastq.gz",
  output:
    temp("alignment_files/IMNU_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam")
  threads:
    30
  resources:
    mem_mb=50000,
    runtime=12*60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  shell:
    """
    STAR --runThreadN {threads} \
        --genomeDir /home/jgroh/Pterocarya/genome_assemblies/IMNU \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix alignment_files/IMNU_RNAseq/published_data/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
    """


rule STAR_align_hap1_published_data:
  input:
    fa=PSTE_UCD1470_HAP1,
    index="/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/SA",
    r1="sra_fasta_paired/{sample}-trimmed-pair1.fastq.gz",
    r2="sra_fasta_paired/{sample}-trimmed-pair2.fastq.gz",
  output:
    "alignment_files/hap1_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam"
  threads:
    30
  resources:
    mem_mb=50000,
    runtime=12*60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  shell:
    """
    STAR --runThreadN {threads} \
        --genomeDir /home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1 \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix alignment_files/hap1_RNAseq/published_data/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMstrandField intronMotif 
    """ # last option for compatibility with BRAKER3

rule STAR_align_hap2_published_data:
  input:
    fa=PSTE_UCD1470_HAP2,
    index="/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/SA",
    r1="sra_fasta_paired/{sample}-trimmed-pair1.fastq.gz",
    r2="sra_fasta_paired/{sample}-trimmed-pair2.fastq.gz",
  output:
    "alignment_files/hap2_RNAseq/published_data/{sample}.Aligned.sortedByCoord.out.bam"
  threads:
    30
  resources:
    mem_mb=50000,
    runtime=12*60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  shell:
    """
    STAR --runThreadN {threads} \
        --genomeDir /home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2 \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix alignment_files/hap2_RNAseq/published_data/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMstrandField intronMotif 
    """

rule STAR_align_hap1_mydata:
  input:
    fa=PSTE_UCD1470_HAP1,
    index="/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1/SA",
    r1="RNAseq/trimmed_reads/{sample}-trimmed-pair1.fastq.gz",
    r2="RNAseq/trimmed_reads/{sample}-trimmed-pair2.fastq.gz",
  output:
    "alignment_files/hap1_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam"
  threads:
    30
  resources:
    mem_mb=50000,
    runtime=12*60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  shell:
    """
    STAR --runThreadN {threads} \
        --genomeDir /home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap1 \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix alignment_files/hap1_RNAseq/mydata/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMstrandField intronMotif 
    """

rule STAR_align_hap2_mydata:
  input:
    fa=PSTE_UCD1470_HAP2,
    index="/home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/SA",
    r1="RNAseq/trimmed_reads/{sample}-trimmed-pair1.fastq.gz",
    r2="RNAseq/trimmed_reads/{sample}-trimmed-pair2.fastq.gz",
  output:
    "alignment_files/hap2_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam"
  threads:
    30
  resources:
    mem_mb=50000,
    runtime=12*60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  shell:
    """
    STAR --runThreadN {threads} \
        --genomeDir /home/jgroh/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2 \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix alignment_files/hap2_RNAseq/mydata/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMstrandField intronMotif 
    """

rule STAR_align_IMNU_mydata:
  input:
    fa=IMNU_FASTA,
    index="/home/jgroh/Pterocarya/genome_assemblies/IMNU/SA",
    r1="RNAseq/FASTQ/{sample}_R1.fastq.gz",
    r2="RNAseq/FASTQ/{sample}_R2.fastq.gz",
  output:
    temp("alignment_files/IMNU_RNAseq/mydata/{sample}.Aligned.sortedByCoord.out.bam")
  threads:
    30
  resources:
    mem_mb=50000,
    runtime=12*60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/star.yaml"
  shell:
    """
    STAR --runThreadN {threads} \
        --genomeDir /home/jgroh/Pterocarya/genome_assemblies/IMNU \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix alignment_files/IMNU_RNAseq/mydata/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
    """

rule filter_STAR_alignment_published_data:
  input:
    "alignment_files/{genome}/published_data/{sample}.Aligned.sortedByCoord.out.bam"
  output:
    "alignment_files/{genome}/published_data/{sample}_filtered.bam"
  resources:
    mem_mb=10000,
    runtime=60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "samtools view --exclude-flags 2816 --min-MQ 255 {input} -o {output}"


rule filter_STAR_alignment_mydata:
  input:
    "alignment_files/{genome}/mydata/{sample}.Aligned.sortedByCoord.out.bam"
  output:
    "alignment_files/{genome}/mydata/{sample}_filtered.bam"
  resources:
    mem_mb=10000,
    runtime=60
  conda:
    "/home/jgroh/heterodichogamy/conda_envs/samtools_1.17.yaml"
  shell:
    "samtools view --exclude-flags 2816 --min-MQ 255 {input} -o {output}"
    

rule extract_focal_region_for_IGV_IMNU:
  input:
    "alignment_files/IMNU/{sample}_filtered.cram"
  output:
    "alignment_files/IMNU/{sample}_focal_region.bam"
  conda:
    SAMTOOLS
  shell:
    "samtools view {input} Chr11:4450000-4600000 -o {output}"


rule extract_UCD2_SNPs:
  input:
    "calls/IMNU/variant_filtered_Chr11.vcf.gz"
  output:
    "calls/IMNU/UCD2_heterozygous_SNPs.bed"
  conda:
    SAMTOOLS
  shell:
    """
    bcftools view -m2 -M2 -v snps -s "PSTE_UCD_2" {input} | \
        bcftools filter -i 'GT="het"' | \
        bcftools query -f '%CHROM\t%POS-1\\n' | \
        awk '{{print $1 "\t" $2-1 "\t" $2}}' > {output}
    """

rule exclude_SNPs_in_CNV_regions:
  input:
    "calls/IMNU/UCD2_heterozygous_SNPs.bed"
  output:
    "calls/IMNU/UCD2_heterozygous_SNPs_CNV_excluded.txt" # this is for input to pixy so should now be 1-based
  conda:
    BEDTOOLS
  shell:
    """
    bedtools intersect -v -a calls/IMNU/UCD2_heterozygous_SNPs.bed -b CNV_regions_to_exclude_IMNU.bed | \
        awk '{{print $1 "\t" $3}}' > {output}
    """



rule tabix_Chr11:
  input:
    "calls/IMNU/variant_filtered_Chr11.vcf.gz",
  output:
    "calls/IMNU/variant_filtered_Chr11.vcf.gz.tbi",
  conda:
    PIXY
  shell:
    "tabix {input}"


rule pi_at_UCD2_SNPs:
  input:
    vcf="calls/IMNU/variant_filtered_Chr11.vcf.gz",
    tbi="calls/IMNU/variant_filtered_Chr11.vcf.gz.tbi",
    pops="all_Ptero_pixy_popfile.txt",
    sites="calls/IMNU/UCD2_heterozygous_SNPs_CNV_excluded.txt"
  output:
    "results/pixy/UCD2_SNPs_pi.txt"
  conda:
    PIXY
  threads:
    50
  shell:
    """
    pixy --stats pi --vcf {input.vcf} \
        --populations {input.pops} \
        --sites_file {input.sites} \
        --window_size 1 \
        --output_folder "results/pixy" --output_prefix "UCD2_SNPs" \
        --n_cores {threads} \
        --bypass_invariant_check 'yes'
    """

rule pi_windows:
  input:
    vcf="calls/IMNU/allsites_filtered_Chr11.vcf.gz",
    tbi="calls/IMNU/allsites_filtered_Chr11.vcf.gz.tbi",
    pops="all_Ptero_pixy_popfile.txt",
  output:
    "results/pixy/1kb_pi.txt"
  conda:
    PIXY
  threads:
    50
  shell:
    """
    pixy --stats pi --vcf {input.vcf} \
        --populations {input.pops} \
        --chromosomes 'Chr11' \
        --interval_start 4400000 --interval_end 4700000 \
        --window_size 1000 \
        --output_folder "results/pixy" --output_prefix 1kb \
        --n_cores {threads} \
    """

