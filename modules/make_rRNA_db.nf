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
https://github.com/igordot/genomics/blob/master/workflows/rrna-ref.md
*/

/*
Enable Nextflow DSL2
*/
nextflow.enable.dsl=2

/*
Configurable variables for process
*/
params.outdir = "./db"
params.species = "mus_musculus"

process make_rRNA_db{
    executor = 'slurm'
    memory '64 GB'
    cpus 2

    publishDir "${params.outdir}", mode: 'link'

    output:
    path("*.g19.fa")

    shell:
    '''
    wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/12.0/sequences/rnacentral_species_specific_ids.fasta.gz
    
    #Convert multi-line sequences to single-line (uses fasta_formatter from FASTX-Toolkit)
    gzip -cd rnacentral_species_specific_ids.fasta.gz \
        | fasta_formatter -w 0 \
        | gzip \
        > rnacentral.nowrap.fasta.gz
    
    #Remove empty lines, replace spaces with underscores, and keep only rRNA sequences
    gzip -cd rnacentral.nowrap.fasta.gz \
        | sed '/^$/d' \
        | sed 's/\s/_/g' \
        | grep -E -A 1 "ribosomal_RNA|rRNA" \
        | grep -v "^--$" \
        | gzip \
        > rnacentral.ribosomal.nowrap.fasta.gz

    #Extract species-specific rRNA sequences
    zcat rnacentral.ribosomal.nowrap.fasta.gz \
        | grep -A 1 -F -i "!{params.species}" \
        | grep -v "^--$" \
        | fasta_formatter -w 80 \
        > rRNA.!{params.species}.fa
    
    #Remove sequences less than 19bp (uses reformat.sh from BBMap suite) (necessary for sortmerna process)
    reformat.sh \
        in=rRNA.!{params.species}.fa \
        out=rRNA.!{params.species}.g19.fa \
        minlength=19
    '''
}