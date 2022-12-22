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

//Configurable variables for pipeline
params.species = "mus_musculus"
params.genome = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.ens_rls = 107
params.reads = "$PWD/fastq/*{_R,_}{1,2}*.fastq.gz"
params.singleEnd = false
params.outdir = "./results"
params.name = false
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = false
params.rmrRNA = false
params.aligner = 'hisat2'
params.db = "${workflow.projectDir.getParent()}/${params.genome}_${params.ens_rls}/"
params.multiqc_config = "${workflow.projectDir}/multiqc_config.yaml"

//Include modules to main pipeline
include { make_rRNA_db } from './modules/make_rRNA_db.nf' addParams(species: params.species, outdir: params.db)
include { downloadGTF } from './modules/downloadGTF.nf' addParams(outdir: params.db, gtf: params.gtf)
include { downloadFASTA } from './modules/downloadFASTA.nf' addParams(outdir: params.db, fasta: params.fasta)
include { extractExons } from './modules/extractExons.nf' addParams(outdir: params.db)
include { extractSpliceSites } from './modules/extractSpliceSites.nf' addParams(outdir: params.db)
include { extractLength } from './modules/extractLength.nf' addParams(outdir: params.db)
include { genEnsAnnot } from './modules/genEnsAnnot.nf' addParams(outdir: params.db)
include { dexSeqPrep } from './modules/dexSeqPrep.nf' addParams(outdir: params.db)

include { fastqc as pretrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'pretrim_fastqc')
include { fastqc as posttrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'posttrim_fastqc')
include { fastqc as sortmerna_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'sortmerna_fastqc')
include { fastqc as ribotrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'ribotrim_fastqc')
include { trim_galore } from './modules/trim_galore.nf'

include { dexSeqCount } from './modules/dexSeqCount.nf'

//Only include rmRNA processes if user wants ribosomal RNA filtering (default = true)
if( params.rmrRNA ){
    include { sortMeRNA } from './modules/sortMeRNA.nf'
    include { rn45s_align } from './modules/rn45s_align.nf'
    include { rn45s_sort } from './modules/rn45s_sort.nf'
    //Include repair step if reads are paired-end
    if( !params.singleEnd ){
        include { repair_fastq } from './modules/repair_fastq.nf'
    }
}

//Only include hisat2 processes if params.aligner is 'hisat2' (default = hisat2)
if( params.aligner == 'hisat2') {
    include { makeHisatIndex } from './modules/makeHisatIndex.nf' addParams(outdir: params.db)
    include { hisat2_align } from './modules/hisat2_align.nf'
    include { hisat2_sort } from './modules/hisat2_sort.nf'
}
include { htseq_count } from './modules/htseq_count.nf'
include { multiqc } from './modules/multiqc.nf' addParams(multiqc_config: params.multiqc_config)

workflow DLP {

    main:
    /*
    PREPROCESSING WORKFLOW:
    Download Fasta, 
    Download GTF, 
    Build HISAT2/STAR indexes, 
    Build BED file,
    */
    //Generate folder where gtf file, fasta file, rRNA fasta, and aligner files will be contained
    if(file(params.db).isEmpty()){
        file(params.db).mkdir()
    }

    //Generate rRNA db for specified species (default = mus musculus)
    if(file("${params.db}*.g19.fa").isEmpty()){
        make_rRNA_db()
        rRNAdb = make_rRNA_db.out
    }
    else{
        Channel
        .fromPath( "${params.db}*.g19.fa" )
        .set{ rRNAdb }
    }
    
    //Detect if gtf file exists, if not, download gtf
    if(file("${params.db}*.gtf").isEmpty()){
        downloadGTF()
        gtf = downloadGTF.out
    }
    else{
        Channel
        .fromPath( "${params.db}*.gtf" )
        .set{ gtf }
    }

    //Detect if fasta file exists, if not, download fasta
    if(file("${params.db}*primary_assembly.fa").isEmpty()){
        downloadFASTA()
        fasta = downloadFASTA.out
    }
    else{
        Channel
        .fromPath( "${params.db}*${params.genome}*.fa" )
        .set{ fasta }
    }

    //Detect if genome.ss exists, if not, extract splice sites from gtf
    if(file("${params.db}*.ss").isEmpty()){
        extractSpliceSites(gtf)
        ss = extractSpliceSites.out
    }
    else{
        Channel
        .fromPath( "${params.db}*.ss" )
        .set{ ss }
    }

    //Detect if genome.exon exists, if not, extract exons from gtf
    if(file("${params.db}*.exon").isEmpty()){
        extractExons(gtf)
        exon = extractExons.out
    }
    else{
        Channel
        .fromPath( "${params.db}*.exon" )
        .set{ exon }
    }

    //Detect if hisat indexes exist, if not, make hisat indexes
    if(file("${params.db}*.ht2").isEmpty()){
        makeHisatIndex(gtf.merge(fasta, ss, exon))
        ht2_base = makeHisatIndex.out.ht2_base
    }
    else{
        ht2_base = Channel.from("${params.db}${params.genome}_${params.ens_rls}")
    }

    //Detect if length file exists, if not, extract lengths
    if(file("${params.db}*Length.tsv").isEmpty()){
        extractLength(gtf)
        lengths = extractLength.out
    }
    else{
        Channel
        .fromPath( "${params.db}*Length.tsv" )
        .set{ lengths }
    }

    //Detect if annotation file exists, if not, generate annotation file
    if(file("${params.db}Mouse_Gene_Annotations.csv").isEmpty()){
        genEnsAnnot(lengths)
    }

    //Detect if gff file exists, if not run dexseq prep to create gff file
    if(file("${params.db}*.gff").isEmpty()){
        dexSeqPrep(gtf)
        gff = dexSeqPrep.out
    }
    else{
        Channel
        .fromPath( "${params.db}*.gff" )
        .set{ gff }
    }    

    /*
    TRIM_FILTER WORKFLOW:
    Trim reads with Trim Galore, 
    Filter reads for rRNA,
    Perform FastQC on reads before and after trim/filter
    */
    //Create channel to auto detect paired end data. Specifiy --singleEnd if your fastq is in single-end format
    Channel
    .fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set{ reads_ch }

    //Perform fastqc on pretrimmed reads, trim reads, and perform fastqc on trimmed reads
    pretrim_fastqc(reads_ch)
    trim_galore(reads_ch)
    posttrim_fastqc(trim_galore.out[0])

    //Steps if user specifies for samples to have ribosomal RNA filtering (default = true)
    //rRNA filter trimmed reads
    if( params.rmrRNA ) {
        sortMeRNA(trim_galore.out[0].combine(rRNAdb))
        //skip repair_fastq process if fastq is single-end
        if( params.singleEnd ) {
            sortmerna_fastqc(sortMeRNA.out[0])
            filtered_fastq = sortMeRNA.out[0]
        }
        else {
            repair_fastq(sortMeRNA.out[0])
            sortmerna_fastqc(repair_fastq.out)
            filtered_fastq = repair_fastq.out
        }
    }

    else {
        filtered_fastq = trim_galore.out[0]
    }
        

    /*
    ALIGN_COUNT WORKFLOW:
    Trim reads with Trim Galore, 
    Filter reads for rRNA,
    Perform FastQC on reads before and after trim/filter
    */
    //Conditional processes depending on aligner chosen (default = hisat2)
    if( params.aligner == 'hisat2') {
        hisat2_align(filtered_fastq.combine(ht2_base))
        sam = hisat2_align.out[0]
        hisat2_sort(hisat2_align.out[0])
        output = hisat2_sort.out
    }
    if( params.aligner == 'star') {
        //TODO: add processes for STAR aligner
    }
    dexSeqCount(sam.combine(gff))
    htseq_count(hisat2_sort.out)

    /*
    POST_PROCESS WORKFLOW:
    MultiQC on PreTrim_FastQC, Trimgalore_Report, 
    PostTrim_FastQC, SortMeRNA_Report, PostSortMeRNA_FastQC,
    HiSAT2_Report, HTSeq_Report
    */

    //Combine all output files necessary for multiqc report into one channel
    //Wait until htseq_count is finished to concat all the reports into a single channel
    multiqc_input = htseq_count.out.collect()
    multiqc_input.concat(trim_galore.out.trimming_report.collect(), posttrim_fastqc.out.collect(), sortMeRNA.out.smr_log.collect(), sortmerna_fastqc.out.collect(), hisat2_align.out.align_report.collect(), pretrim_fastqc.out.collect())
    
    multiqc_input.view()

    //Run MultiQC to create compiled QC report
    multiqc(params.multiqc_config, multiqc_input)
}

/*
MAIN: 
Workflow Execution
*/
workflow {
    DLP()
}
