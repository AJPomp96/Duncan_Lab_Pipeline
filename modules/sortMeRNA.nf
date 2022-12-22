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
params.pubdir = "sortmerna"
params.singleEnd = false
params.rRNAdb = "${params.db}*.g19.fa"

/*
Run sortmerna on each read stored within the trim_galore channel
*/
process sortMeRNA {
    maxForks 3
    executor = 'slurm'
    memory '32 GB'
    cpus 4

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'
    
    input:
    tuple val(file_id), file(reads), file(rRNAdb)

    output:
    tuple val(file_id), path("*rmrRNA*.fq.gz")
    path("*.log"), emit: smr_log

    script:
    if( params.singleEnd )
    """
    echo $params.singleEnd
    
    sortmerna --ref ${rRNAdb} \
    --reads $reads \
    --fastx \
    --workdir . \
    --aligned ${file_id}_rRNA --other ${file_id}_rmrRNA \
    --threads 10 \
    --log
    """

    else 
    """
    echo $params.singleEnd

    sortmerna --ref ${rRNAdb} \
    --reads ${reads[0]} --reads ${reads[1]} \
    --fastx --paired_out --out2 \
    --workdir . \
    --aligned ${file_id}_rRNA --other ${file_id}_rmrRNA \
    --threads 10 \
    --log
    """
}
