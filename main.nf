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
    --card_db                   path to the CARD json Database for Antimicrobial Resistance Genes prediction [default: $params.card_db]

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
    --species                   bacterial species to assemble (e.g. Ecoli, Kpneumoniae, Salmonella, Ecloacae, Saureus) [default: $params.species]
    --species_taxid             NCBI TaxID of the bacterial species to assemble [default: $params.species_taxid]
    --mash_dataset              path to mash dataset prepared for the species [default: $params.mash_dataset]
    --busco_lineage             to specify according to the bacterial species. e.g. enterobacterales_odb10, bacillales_odb10... check BUSCO [default: $params.busco_lineage]
    --amrfinder_organism        To specify for PointMutation detection if not Ecoli, Salmonella, Kpneumoniae or Saureus.
                                Can be among these: Acinetobacter_baumannii, Campylobacter,
                                Clostridioides_difficile, Enterococcus_faecalis, Enterococcus_faecium,
                                Escherichia, Klebsiella, Neisseria, Pseudomonas_aeruginosa,
                                Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius,
                                Streptococcus_agalactiae, Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae.
                                The amrfinderplus will be run if not specified, but no point mutations are detected.
                                [default: $params.amrfinder_organism]
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
    include {extract_kraken} from './modules/kraken2.nf' params(output: params.output)

    // include mash_screen
    if (params.mash_dataset){
      include {mash_screen} from './modules/mash.nf' params(output: params.output)
    }
    // AMR analysis modules
    include {amrfinderplus} from './modules/amrfinderplus.nf' params(output: params.output)
    include {card_rgi} from './modules/card.nf' params(output: params.output)


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
    quast(contigs_ch, illumina_clean_ch)

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
    krak_res = kraken2nt_contigs.out[0]
    krak_report = kraken2nt_contigs.out[1]

    // KrakenTools
// Maybe I should treat all this with parameters
    if(params.species_taxid){
      sp_taxid = params.species_taxid
    }
    else if (params.species == "Ecoli") {
      sp_taxid = '562'
    }
    else if (params.species == "Kpneumoniae"){
      sp_taxid = '573'
    }
    else if (params.species == "Salmonella"){
      sp_taxid = '590'
    }
    else if (params.species == "Ecloacae"){
      sp_taxid = '550'
    }
    else if (params.species == "Saureus"){
      sp_taxid = '1280'
    }
    else {
      exit 1, "No species or species taxid specified for retrieving the species of interest"
    }

    extract_kraken(contigs_ch,krak_res,krak_report,sp_taxid,params.krakentools_extract)
    deconta_contigs_ch = extract_kraken.out[0]

    //*************************************************
    // STEP 4 - Find closest relative with Mash
    //*************************************************
    // using the mash dataset provided
    if(params.mash_dataset){
      mash_screen(contigs_ch,params.species,params.mash_sketch)
      // later, mash_screen will output fasta file of closest relative for mapping
    }

    //*************************************************
    // STEP 5 - Map to the closest with Minimap2
    //*************************************************
    // minimap2 with PAF Output - use all contigs: contigs_ch
    // compare the number of contigs mapped with the number of deconta_contigs_ch

    // mapping could be done on same ref for all... need to select from mash results
    // might need a second pipeline for mappings and core/pan genomes analysis,
    // after the assemblies and decontamination are done...

    //*************************************************
    // STEP 6 - ARGs search: CARD RGI and AMRFinderPlus
    //*************************************************
    // on the deconta_contigs_ch
    // AMRFinderPlus NCBI

    if(params.species == "Ecoli"){
      organism = 'Escherichia'
      amrfinderplus(deconta_contigs_ch,organism)
    }
    else if (params.species == "Kpneumoniae"){
      organism = 'Klebsiella'
      amrfinderplus(deconta_contigs_ch,organism)
    }
    else if (params.species == "Salmonella"){
      organism = 'Salmonella'
      amrfinderplus(deconta_contigs_ch,organism)
    }
    else if (params.species == "Saureus"){
      organism = 'Staphylococcus_aureus'
      amrfinderplus(deconta_contigs_ch,organism)
    }
    else if (params.amrfinder_organism){
      amrfinderplus(deconta_contigs_ch,params.amrfinder_organism)
    }
    else {
      amrfinderplus_no_species(deconta_contigs_ch)
    }

    // CARD Resistance Genes Identifier

    card_rgi(deconta_contigs_ch,params.card_db)

    //*************************************************
    // STEP 7 - Bakta annotation
    //*************************************************
    // bakta annotation of deconta contigs (and mapped contigs)

    //*************************************************
    // STEP 8 - PlasmidFinder et al. Platon ? MGEFinder..
    //*************************************************
  //Platon
  //  PlasForest
  //  MOB-recon
  // on deconta_contigs_ch
}
