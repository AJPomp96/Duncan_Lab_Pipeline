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
params.pubdir = "fastqc"

/*
Run fastqc on fastq files
*/
process fastqc {
    executor = 'slurm'
    memory '8 GB'
    cpus 2
    
    publishDir "${params.outdir}/${params.pubdir}"

    input:
    tuple val(file_id), file(reads)

    output:
    file("*.{html,zip}")

    script:
    """
    fastqc -t 2 $reads
    """
}