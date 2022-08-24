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
params.pubdir = "sortmerna"
params.singleEnd = false
params.rRNAdb = "/work/ajpompet/01Aug2022_Nextflow_Test/db/rRNA-l19.mus_musculus.fa"

/*
Run sortmerna on each read stored within the trim_galore channel
*/
process sortMeRNA {
    executor = 'slurm'
    memory '32 GB'
    cpus 8

    publishDir "${params.outdir}/${params.pubdir}"
    
    input:
    tuple val(file_id), file(reads)

    output:
    tuple val(file_id), path("*rmrRNA*.fq.gz")

    script:
    if( params.singleEnd )
    """
    echo $params.singleEnd
    
    sortmerna --ref ${params.rRNAdb} \
	--reads $reads \
	--fastx \
    --workdir . \
	--aligned ${file_id}_rRNA --other ${file_id}_rmrRNA \
	--threads 20
    """

    else 
    """
    echo $params.singleEnd

    sortmerna --ref ${params.rRNAdb} \
    --reads ${reads[0]} --reads ${reads[1]} \
    --fastx --pair_out --out2 \
    --workdir . \
    --aligned ${file_id}_rRNA --other ${file_id}_rmrRNA \
    --threads 20
    """
}