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
params.singleEnd = false
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = false

process hisat2_align {
    executor = 'slurm'
    memory '256 GB'
    cpus 8

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), path(reads), val(ht2_base)

    output:
    tuple val(file_id), path("*.sam")
    path("*_align_report.txt"), emit: align_report

    script:
    /*
    Check for params that define rnastrandness, default is reverse_stranded
    */
    rnastrandness = ''
    if (params.forward_stranded && !params.unstranded){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    }
    else if (params.reverse_stranded && !params.unstranded){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }

    if( params.singleEnd )
    """
    hisat2 -p ${task.cpus} \
    --verbose \
    --phred33 \
    --dta \
    --summary-file ${file_id}_align_report.txt \
    -x ${ht2_base} \
    -U ${reads} \
	${rnastrandness} \
    -S ${file_id}.sam
    """

    else
    """
    hisat2 -p ${task.cpus} \
    --verbose \
    --phred33 \
    --dta \
    --summary-file ${file_id}_align_report.txt \
    -x ${ht2_base} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
	${rnastrandness} \
    -S ${file_id}.sam
    """
}