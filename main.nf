#!/usr/bin/env nextflow
nextflow.enable.dsl=2

start_var = Channel.from("""
********* Start running nf-wgs_amr pipeline *********
nf-wgs_amr is a workflow for genomics qc, assembly, decontamination by
taxonomic classification, analysis of Antimicrobial Resistance Genes
**************************************
""")
start_var.view()

if (params.help) { exit 0, helpMSG() }

// Help Message
def helpMSG() {
    log.info """
    ********* Assembly and taxonomic classification workflow for (viral) metagenomics *********

        Usage example:
    nextflow run main.nf --illumina illumina/ -profile itrop --species Ecoli

        Input:
    --illumina                  path to the directory containing the illumina read file (fastq) (default: $params.illumina)
        Optional input:
    --k2nt_db                   path to the Kraken2 nucleotide database (e.g. nt) [default: $params.k2nt_db]
    --species                   bacterial species to assemble (e.g. Ecoli, Kpneumoniae, Salmonella, Ecloacae, Saureus) [default: $params.species]
      Output:
    --output                    path to the output directory (default: $params.output)

        Outputed directories:
    qc                          The reads file after qc, qc logs and host mapping logs
    assembly                    The spades assembly output directory
    taxonomic_classif           The taxonomic classifications at contigs level

        Basic Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB [default: $params.memory]

        Workflow Options:
    --mash_dataset              path to mash dataset prepared for the species [default: $params.mash_dataset]
    --busco_lineage             to specify is species is not among: Ecoli, Kpneumoniae, Salmonella, Ecloacae, Saureus [default: $params.busco_lineage]
        Nextflow options:
    -profile                    change the profile of nextflow both the engine and executor more details on github README
    -resume                     resume the workflow where it stopped
    """
}

workflow {

    // error handling
    if (
        workflow.profile.contains('itrop') ||
        workflow.profile.contains('ifb')
    ) { "executer selected" }
    else { exit 1, "No executer selected: -profile itrop/ifb"}


    //*************************************************
    // STEP 0 - Include needed modules
    //*************************************************

    include {fastp} from './modules/fastp.nf' params(output: params.output)
    // including assembler module
    include {spades} from './modules/spades.nf' params(output: params.output)
    // include quast, busco, checkm
    include {quast} from './modules/quast.nf' params(output: params.output)
    include {busco} from './modules/busco.nf' params(output: params.output)

    // including Kraken2 - nucleotide level
    include {kraken2nt_contigs} from './modules/kraken2.nf' params(output: params.output)


    //*************************************************
    // STEP 1 QC with fastp
    //*************************************************

    // DATA INPUT ILLUMINA
    illumina_input_ch = Channel
        .fromFilePairs( "${params.illumina}/*_R{1,2}*.fastq{,.gz}", checkIfExists: true)
        .view()

    // run fastp module
    fastp(illumina_input_ch)
    illumina_clean_ch = fastp.out[0]

    //*************************************************
    // STEP 2 - Assembly
    //*************************************************
    spades(illumina_clean_ch)
    contigs_ch = spades.out[0]

    //*************************************************
    // STEP 2 - Assembly QC Quast, Busco
    //*************************************************
    // QUAST Assembly QC
    quast(contigs_ch,illumina_clean_ch)

    // BUSCO completeness - only available for IFB or try Singularity
/*
    if (params.species == "Ecoli") || (params.species == "Kpneumoniae") ||
      (params.species == "Salmonella") || (params.species == "Ecloacae")
    {
      lineage = "enterobacterales_odb10"
    }
    else if (params.species == "Saureus"){
      lineage = "bacillales_odb10"
      //Bacteria; Terrabacteria group; Firmicutes; Bacilli; Bacillales; Staphylococcaceae; Staphylococcus
    }
    else {
      if (params.busco_lineage){
        lineage = params.busco_lineage
      }
      else {
        exit 1, "unknow lineage for this bacterial species! \
        Please provide lineage with --busco_lineage (e.g. enterobacterales_odb10)"
      }

    }
    busco(contigs_ch, lineage)
*/
    //*************************************************
    // STEP 3 - decontamination with Kraken2
    //*************************************************
    //kraken2nt contigs
    kraken2nt_contigs(contigs_ch, params.k2nt_db)
    // KrakenTools
//    deconta_contigs_ch = kraken_tools(contigs_ch, params.species)

    //*************************************************
    // STEP 4 - Find closest relative with Mash
    //*************************************************
    // using the mash dataset provided
    // grep the sort output for "complete genome" and select top hit
    //retrieve the fasta file from Antimicrobial

    //*************************************************
    // STEP 5 - Map to the closest with Minimap2
    //*************************************************
    // minimap2 with PAF Output - use all contigs: contigs_ch
    // compare the number of contigs mapped with the number of deconta_contigs_ch

    //*************************************************
    // STEP 6 - ARGs search: CARD RGI and AMRFinderPlus
    //*************************************************

    //*************************************************
    // STEP 7 - Bakta annotation
    //*************************************************
    // bakta annotation of deconta contigs (and mapped contigs)

    //*************************************************
    // STEP 8 - PlasmidFinder et al. Platon ? MGEFinder..
    //*************************************************

  //  PlasForest
  //  MOB-recon
}
