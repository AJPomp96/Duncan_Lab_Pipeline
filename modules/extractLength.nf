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

process extractLength {
    executor = 'slurm'
    memory '32 GB'
    cpus 2

    publishDir "${params.outdir}", mode: 'link'
    
    input:
    path(gtf)

    output:
    path("*.tsv")

    shell:
    '''
    #!/usr/bin/env Rscript
    
    library(GenomicFeatures)

    #create txdb from gtf file
    txdb <- makeTxDbFromGFF("!{gtf}", format="gtf")
    
    #collect exons per gene id
    exons.list.per.gene <- exonsBy(txdb,by="gene")

    #for each gene reduce all exons to a set of non overlapping exons, calculate their lengths (widths) and sum them
    exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

    #to see them in table format, unlist them
    unlist_geneLength <- unlist(exonic.gene.sizes)
    write.table(unlist_geneLength, "EnsMm_!{params.genome}_!{params.ens_rls}_Length.tsv", quote=FALSE, col.names=FALSE)
    '''
}