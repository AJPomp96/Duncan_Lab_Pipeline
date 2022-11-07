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
params.pubdir = "trim_galore"
params.singleEnd = false

/*
Run trim_galore on each read stored within the reads_ch channel
*/
process trim_galore {
    maxForks 3
    executor = 'slurm'
    memory '8 GB'
    cpus 8

    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    tuple val(file_id), file(reads)

    output:
    tuple val(file_id), path("*.fq.gz")

    script:
    if ( params.singleEnd )
    """
    trim_galore \
    --length 35 \
    --quality 28 \
    --phred33 \
    --cores ${task.cpus} \
    $reads
    """

    else
    """
    trim_galore \
    --length 35 \
    --quality 28 \
    --paired \
    --phred33 \
    --clip_R2 3 \
    --cores ${task.cpus} \
    $reads
    """
}