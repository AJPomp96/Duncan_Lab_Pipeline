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
params.pubdir = "rn45s_align"
params.singleEnd = false

process rn45s_align {
    executor = 'slurm'
    memory '128 GB'
    cpus 8

	publishDir "${params.outdir}/${params.pubdir}"

    input:
    tuple val(file_id), file(reads)

    output:
    tuple val(file_id), path("*_ribo_reads.sam")

    script:
    if( params.singleEnd )
    """
    hisat2 -p ${task.cpus} \
	--verbose \
	--phred33 \
	--dta \
	--fr \
	--summary-file ${file_id}_RiboAlignStat.txt \
	-x /work/ajpompet/EnsMm_grc39_104/RNA45SN5_Index \
	-U $reads \
	-S ${file_id}_ribo_reads.sam
    """

    else
    """
    hisat2 -p ${task.cpus} \
	--verbose \
	--phred33 \
	--dta \
	--fr \
	--summary-file ${file_id}_RiboAlignStat.txt \
	-x /work/ajpompet/EnsMm_grc39_104/RNA45SN5_Index \
	-1 ${reads[0]} \
	-2 ${reads[1]} \
	-S ${file_id}_ribo_reads.sam
    """
}