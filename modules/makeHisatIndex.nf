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
Configurable variables for process
*/
params.outdir = "./db"

process makeHisatIndex {
    executor = 'slurm'
    memory '256 GB'
    cpus 8

    publishDir "${params.outdir}", mode: 'link'
    
    input:
    tuple path(gtf), path(fasta), path(ss), path(exon)

    shell:
    '''
    hisat2-build -f \
    -p 8 \
    --ss !{ss} \
    --exon !{exon} \
    !{fasta} \
    !{params.genome}_!{params.ens_rls}
    '''
}
