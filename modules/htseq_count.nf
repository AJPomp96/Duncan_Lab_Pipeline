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
params.pubdir = "htseq_count"

process htseq_count {
    executor = 'slurm'
    memory '32 GB'
    cpus 2

    publishDir "${params.outdir}/${params.pubdir}", mode: 'link'

    input:
    tuple val(file_id), file(bam)

    output:
    path("*_rf_GeneCount.txt")

    script:
    """
    htseq-count \
	-i gene_id \
	-r pos \
	-f bam \
	-s reverse \
	-m union \
	--type exon \
	${bam} \
	${params.db}*.gtf > ${file_id}_rf_GeneCount.txt
    """
}
