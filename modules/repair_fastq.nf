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
params.pubdir = "repair_fastq"

process repair_fastq {
    executor = 'slurm'
    memory '32 GB'
    cpus 2

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), file(reads)

    output:
    tuple val(file_id), path("*_repaired_*.fq")

    script:
    """
    repair.sh \
    in1=${reads[0]} \
    in2=${reads[1]} \
    out1=${file_id}_repaired_R1.fq \
    out2=${file_id}_repaired_R2.fq \
    outsingle=${file_id}_singletons.fq
    """
}