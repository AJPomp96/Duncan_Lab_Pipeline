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
params.gtf = ""

process downloadGTF {
    executor = 'slurm'
    memory '64 GB'
    cpus 2

    publishDir "${params.outdir}", mode: 'copy'
    
    output:
    path("*.gtf")

    shell:
    '''
    wget "!{params.gtf}"
    if [ -f *.gz ]; then
        gunzip *.gz
    fi

    GTF_FILE=$(ls *.gtf)

    #fix Lim2/Gm52993 error
    #Lim2 is obscured by a nonsense mediated decay transcript
    #prune Gm52993 from gtf file
    tmpfile=$(mktemp)
    INTERFERING_ID=ENSMUSG00000093639

    wc -l $GTF_FILE

    cat $GTF_FILE\
        | grep -v $INTERFERING_ID > ${tmpfile}

    cat ${tmpfile} > $GTF_FILE

    rm -f ${tmpfile}

    wc -l $GTF_FILE
    '''
}