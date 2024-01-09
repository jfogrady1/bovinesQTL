
import pandas as pd
import os
import subprocess

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
rule all:
     input:
        expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{individual}_Log.final.out', individual = config["samples"]),
        expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/junction_files/{individual}.regtools_junc.txt.gz', individual = config["samples"]),
        '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.bed.gz',
        "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.gtf",
        "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.exons.txt.gz",
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/data/RNA_seq/Bos_taurus_exon.txt.gz',
        "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac_perind_numers.counts.gz",
        "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/leafcutter_ds_effect_sizes.txt",
        expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/{group}.leafcutter.bed.gz', group = config["groups"]),
        expand("/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}_splice.cis_qtl.txt.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}.cis_qtl_pairs.{chrn}.txt.gz",  group = config["groups"], chrn=autosomes),
        expand("/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}_splice.cis_independent_qtl.txt.gz", group = config["groups"])

rule Alignment:
    input:
        genome = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/",
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/{config["samples"][wildcards.individual]}_{{N}}.fq.gz', N=(1,2)),
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.FILTERED.dose.vcf.gz"
    params:
        prefix = lambda wildcards: f'{config["samples"][wildcards.individual]}'
    output:
        ind_vcf = "/home/workspace/jogrady/bovinesQTL/results/SNP_data/imputation/imputed/{individual}_CHR_Renamed.IMPUTED.FILTERED.dose.vcf.gz",
        finallog = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{individual}_Log.final.out',
        interlog = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{individual}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{individual}_Log.out',
        bams = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{individual}_Aligned.sortedByCoord.out.bam'
    threads: 40
    shell:
        '''

        vcftools --gzvcf {input.vcf} --indv {params.prefix} --recode --recode-INFO-all --stdout | bgzip -c > /home/workspace/jogrady/bovinesQTL/results/SNP_data/imputation/imputed/{params.prefix}_CHR_Renamed.IMPUTED.FILTERED.dose.vcf.gz

        STAR-2.7.1a  --genomeDir {input.genome} --runThreadN {threads} \
        --waspOutputMode SAMtag --varVCFfile <(zcat {output.ind_vcf} ) --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random \
        --twopassMode Basic --outSAMstrandField intronMotif --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand zcat \
        --outFileNamePrefix /home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{params.prefix}_ 
        
        '''


# Step 1. Converting bams to juncs
# We now advise that users create junction files using the regtools software as it much faster than our own implementation.

# regtools junctions extract uses the CIGAR strings in each bam file to quantify the usage of each intron.

rule generate_juncs:
    input:
        reads = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{individual}_Aligned.sortedByCoord.out.bam',
        regtools = '/home/workspace/jogrady/bovinesQTL/bin/regtools/build/regtools'
    output:
        filtered_bam = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/{individual}_Aligned_filtered.bam', 
        regtools_junctions = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/junction_files/{individual}.regtools_junc.txt.gz'
    params:
        id = '{individual}'


    singularity:
        "docker://francois4/leafcutter:latest"
    shell:
        '''
        samtools view -h -q 255 {input.reads} | grep -v "vW:i:[2-7]" | samtools view -b > {output.filtered_bam}
        samtools index {output.filtered_bam}
        regtools junctions extract -a 8 -m 50 -M 500000 -s 0 {output.filtered_bam} -o /home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/junction_files/{params.id}.regtools_junc.txt 
        
        gzip /home/workspace/jogrady/bovinesQTL/results/RNA-seq/Splicing_Alignment/junction_files/{params.id}.regtools_junc.txt 
        
        '''



# 3) Generating intron excision ratios with LeafCutter

# The cluster_prepare_fastqtl.py script wraps LeafCutter's leafcutter_cluster_regtools.py script, applies filters to remove introns with low counts or low complexity, and generates input files formatted for QTL mappers.

rule collapse_annotation:
    input:
        gtf = "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.gtf",
        script = "/home/workspace/jogrady/bovinesQTL/bin/collapse_annotation.py"
    output:
        gtf_collapsed = "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.gtf"
    singularity:
        "docker://francois4/leafcutter:latest"
    shell:
        '''
        python3 {input.script} {input.gtf} {output.gtf_collapsed} --collapse_only
        '''
rule exons_list:
    input:
        gtf_collapsed = "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.gtf",
        script = "/home/workspace/jogrady/bovinesQTL/bin/exon_list.py"
    output:
        gtf_exons = "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.exons.txt.gz"
    singularity:
        "docker://francois4/leafcutter:latest"
    shell:
        '''
        python3 {input.script} {input.gtf_collapsed} {output.gtf_exons} 
        '''


rule intron_excision:
    input:       
        script = "/home/workspace/jogrady/bovinesQTL/bin/cluster_prepare_fastqtl.py",
        junc_files_list = "/home/workspace/jogrady/bovinesQTL/data/junction_list.txt",
        collapsed_annotation = "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.gtf",
        exons_list = "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.exons.txt.gz",
        sample_participant_map = "/home/workspace/jogrady/bovinesQTL/data/participants_id.txt"
    params:
        prefix = 'Con_V_Reac'
    output:
        expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.bed.gz', individual = config["samples"]) ,
        expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.bed.gz.tbi', individual = config["samples"]),
        expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.txt', individual = config["samples"])

    shell:
        '''
        python3 {input.script} \
        {input.junc_files_list} \
        {input.exons_list} \
        {input.collapsed_annotation} \
        {params.prefix} \
        {input.sample_participant_map} \
        --leafcutter_dir /home/workspace/jogrady/bovinesQTL/bin/leafcutter \
        -o /home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/

        #grep -P -v '^\t0' /home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.txt > /home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.tmp.txt
        
        #mv /home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.tmp.txt /home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.txt

        '''


rule exon_gtf:
    input: 
        gtf =  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        script = "/home/workspace/jogrady/bovinesQTL/bin/splicing/gtf_to_exon.R"
        
    output:
        gtf_exon = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/data/RNA_seq/Bos_taurus_exon.txt.gz"

    shell:
        '''
        Rscript {input.script} {input.gtf} {output.gtf_exon}
        '''


rule junctionQC:
    input:
        "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac_perind_numers.counts.gz",
    output:
        "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac_filtered_perind_numers.counts.gz"

    params:
        script = "/home/workspace/jogrady/bovinesQTL/bin/leafcutter/scripts/Cluster_QC.R"
    shell:
        "Rscript {params.script} "
        "--outFolder /home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/ "
        "--dataCode Con_V_Reac "
        "--missingness 0.5 "
        "--minratio 0.05 "
        "--maxsize 10 "


rule leafcutter_ds:
    input:
        counts = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac_filtered_perind_numers.counts.gz",
        exon_file = "/home/workspace/jogrady/bovinesQTL/data/Bos_taurus.ARS-UCD1.2.110.genes.exons.txt.gz",
        script = "/home/workspace/jogrady/bovinesQTL/bin/leafcutter/scripts/leafcutter_ds.R",
        groups = "/home/workspace/jogrady/bovinesQTL/data/covariate_RNA_seq.txt"
    output:
        significance = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/leafcutter_ds_cluster_significance.txt",
        effect = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/leafcutter_ds_effect_sizes.txt"
    singularity:
        "docker://lindsayliang/leafcutter:latest"
    shell:
        """
        Rscript {input.script} -o /home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/leafcutter_ds \
        -e {input.exon_file} -p 8 {input.counts} {input.groups}
        """

rule prepare_bed_files_for_qtl:
    input:
        bed = '/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.bed.gz'

    output:
        beds =  expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/{group}.leafcutter.bed.gz', group = config["groups"]),
        index =  expand('/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/{group}.leafcutter.bed.gz.tbi', group = config["groups"])

    shell:
        '''
        zcat {input.bed} | bgzip -c > {output.beds[0]}
        zcat {input.bed} | cut -f 1-67 | bgzip -c > {output.beds[1]}
        zcat {input.bed} | cut -f 1-4,68-127 | bgzip -c > {output.beds[2]}

        # index
        tabix -p bed {output.beds[0]}
        tabix -p bed {output.beds[1]}
        tabix -p bed {output.beds[2]}
        '''

rule sqtl_mapping_permute:
    input:
        bed_4_tensor = multiext("/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/{group}.leafcutter", ".bed.gz", ".bed.gz.tbi"),# Phenotypes
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt", # Covariates
        phenotype_groups = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.txt"
    output:
        outputdir="/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}_splice.cis_qtl.txt.gz"
    resources:
        mem_mb = 6000,
        threads = 20
    params:
       plink_prefix_path="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor", # Genotypes
       group = "{group}"
    shell:
        '''
        mkdir -p /home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/

        prefix="/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{params.group}"

        python3 /home/workspace/jogrady/bovinesQTL/bin/run_tensorqtl.py \
        {params.plink_prefix_path} \
        {input.bed_4_tensor[0]} \
        ${{prefix}} \
        --covariates {input.covariates_file} \
        --groups {input.phenotype_groups} \
        --mode cis
        '''

rule tensorqtl_nominal_splice:
    """Get summary statistics for all tested cis-window SNPs per gene."""
    input:
        bed_4_tensor = expand("/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/{group}.leafcutter.bed.gz", group = config["groups"]),# Phenotypes
        covariates_file=expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt", group = config["groups"]), # Covariates
        phenotype_groups = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.txt",
        script = "/home/workspace/jogrady/bovinesQTL/bin/tensorqtl_nominal.py"
    output:
        parquet_files = expand("/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}.cis_qtl_pairs.{chrn}.parquet", group = config["groups"], chrn = autosomes)
    params:
        plink_prefix_path=expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor", group = config["groups"]), # Genotypes
        outdir = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/",
    shell:
        """
        python3 {input.script} {params.plink_prefix_path[0]} ALL {input.covariates_file[0]} {input.bed_4_tensor[0]} {input.phenotype_groups} {params.outdir}
        python3 {input.script} {params.plink_prefix_path[1]} CONTROL {input.covariates_file[1]} {input.bed_4_tensor[1]} {input.phenotype_groups} {params.outdir}
        python3 {input.script} {params.plink_prefix_path[2]} INFECTED {input.covariates_file[2]} {input.bed_4_tensor[2]} {input.phenotype_groups} {params.outdir}
        """
    
rule parquet_2_txt:
    input:
        parquet_files= "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}.cis_qtl_pairs.{chrn}.parquet",
        script = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/parquet2txt.py'
    output:
        txt_files = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}.cis_qtl_pairs.{chrn}.txt.gz"
    shell:
        '''
        python3 {input.script} {input.parquet_files} {output.txt_files}
        '''

rule sqtl_mapping_conditional:
    input:
        bed_4_tensor = multiext("/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/{group}.leafcutter", ".bed.gz", ".bed.gz.tbi"),# Phenotypes
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt", # Covariates
        phenotype_groups = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/leafcutter_output/Con_V_Reac.leafcutter.phenotype_groups.txt",
         cis = "/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}_splice.cis_qtl.txt.gz"
    output:
        outputdir="/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{group}_splice.cis_independent_qtl.txt.gz"
    resources:
        mem_mb = 6000,
        threads = 20
    params:
       plink_prefix_path="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor", # Genotypes
       group = "{group}"
    shell:
        '''
        mkdir -p /home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/

        prefix="/home/workspace/jogrady/bovinesQTL/results/RNA-seq/sQTL/{params.group}_splice.cis_independent_qtl.txt.gz"

        python3 /home/workspace/jogrady/bovinesQTL/bin/run_tensorqtl.py \
        {params.plink_prefix_path} \
        {input.bed_4_tensor[0]} \
        ${{prefix}} \
        --covariates {input.covariates_file} \
        --groups {input.phenotype_groups} \
        --cis_output {input.cis} \
        --mode cis_independent
        '''