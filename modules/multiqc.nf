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
params.pubdir = "multiqc"
params.multiqc_config = ""

process multiqc {
    executor = 'slurm'
    memory '32 GB'
    cpus 8

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'
    
    input:
    file multiqc_config
    file('multiqc_input/*')

    output:
    file('*multiqc_report.html')

    script:
    """
    multiqc multiqc_input/ --config ${multiqc_config}
    """
}