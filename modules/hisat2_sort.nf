#!/usr/bin/env nextflow

/*
Melinda Duncan Lab
RNA-Seq Pipeline. Started August 2022.

Contributors:
Anthony Pompetti <ajpompet@udel.edu>
Adam Faranda    <abf@udel.edu>
Suhotro Gorai <sugo@udel.edu>
Patrick Dopler <ptdopler@udel.edu>

Methodology adapted from: 
https://github.com/SciLifeLab/NGI-RNAseq/blob/3ffd8fe92d4ba39dfc96e36f67156dc7679808e8/main.nf#L161-L168
*/

/*
Enable Nextflow DSL2
*/
nextflow.enable.dsl=2

/*
Define local params 
*/
params.outdir = "./results"
params.pubdir = "sorted_bam"

process hisat2_sort {
    executor = 'slurm'
    memory '96 GB'
    cpus 8

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), file(sam)

    output:
    tuple val(file_id), path("*.sorted.bam")

    script:
    """
    samtools view \
	-@ ${task.cpus} \
	-bS ${sam} > ${file_id}.bam

    samtools sort \
	-@ ${task.cpus} \
	-m 12G \
	-o ${file_id}.sorted.bam \
    ${file_id}.bam
    
    samtools index ${file_id}.sorted.bam
    """
}