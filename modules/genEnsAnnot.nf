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

process genEnsAnnot {
    executor = 'slurm'
    memory '32 GB'
    cpus 2

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(tsv)

    output:
    path("*.csv")

    shell:
    '''
    #!/usr/bin/env Rscript
    
    library(tidyverse)
    library(biomaRt)
    library(httr)
    
    new_config <- httr::config(ssl_verifypeer = FALSE)
    httr::set_config(new_config, override = FALSE)

    url <- filter(listEnsemblArchives(), version == !{params.ens_rls})$url
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",
                host = url)

    martdf <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version",
                            "external_gene_name", "chromosome_name", 
                            "description", "gene_biotype"), mart = mart)

    length_table <- read.table("!{tsv}", header = F, quote = "") 

    names(length_table) <- c("gene_id", "eu_length")

    print(paste("Rows in length table:", nrow(length_table)))

    names(martdf) <- c("gene_id", "gene_id_version", "SYMBOL", "SEQNAME", "DESCRIPTION", "GENEBIOTYPE")

    martdf <- martdf %>% rowwise() %>% mutate(gene_id_version=gsub(paste0(gene_id,"."),"",gene_id_version))

    print(paste("Rows in meta table:", nrow(martdf)))
    
    gene_annotations <- inner_join(length_table, martdf, by="gene_id")

    print(paste("Rows in annotation table:", nrow(gene_annotations)))

    print(head(gene_annotations))

    write.csv(gene_annotations, file="Mouse_Gene_Annotations.csv")
    '''
}