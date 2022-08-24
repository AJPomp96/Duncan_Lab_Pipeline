#!/usr/bin/env nextflow

/*
Melinda Duncan Lab
RNA-Seq Pipeline. Started August 2022.
Anthony Pompetti
Adam Faranda
Suhotro Gorai
Patrick Dopler
Methodology adapted from https://github.com/SciLifeLab/NGI-RNAseq/blob/3ffd8fe92d4ba39dfc96e36f67156dc7679808e8/main.nf#L161-L168
*/

/*
Enable Nextflow DSL2
*/
nextflow.enable.dsl=2

/*
Define local params 
*/
params.outdir = "./results"
params.pubdir = "rn45s_sort"
params.singleEnd = false

process rn45s_sort {
    executor = 'slurm'
    memory '96 GB'
    cpus 8

    publishDir "${params.outdir}/${params.pubdir}"

    input:
    tuple val(file_id), file(sam)

    output:
    tuple val(file_id), path("*.fq")

    script:
    if( params.singleEnd )
    """
    samtools view \
	-@ ${task.cpus} \
	-bS ${sam} > ${file_id}_ribo_reads.bam

    samtools sort \
	-n \
	-@ ${task.cpus} \
	-m 12G \
	-o ${file_id}_byname_ribo.bam ${file_id}_ribo_reads.bam

    bedtools bamtofastq -i <(samtools view -h -f 4 ${file_id}_byname_ribo.bam) \
    -fq ./${file_id}_ribotrim.fq
    """

    else
    """
    samtools view \
	-@ ${task.cpus} \
	-bS ${sam} > ${file_id}_ribo_reads.bam

    samtools sort \
	-n \
	-@ ${task.cpus} \
	-m 12G \
	-o ${file_id}_byname_ribo.bam ${file_id}_ribo_reads.bam

    bedtools bamtofastq -i <(samtools view -h -f 13 ${file_id}_byname_ribo.bam) \
    -fq ./${file_id}_ribotrim_R1.fq \
    -fq2 ./${file_id}_ribotrim_R2.fq
    """
}