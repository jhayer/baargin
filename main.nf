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
    --plasmidfinder_db          path to the CGE PlasmidFinder database [default: $params.plasmidfinder_db]
    --bakta_db                  path to the bakta annotation database [default: $params.bakta_db]

      Output:
    --output                    path to the output directory (default: $params.output)
    --tmpdir                    path to the tmp directory (default: $params.tmpdir)

        Outputed directories:
    qc                          The reads file after qc, qc logs and host mapping logs
    assembly                    The spades assembly output directory
    taxonomic_classif           The taxonomic classifications at contigs level

        Basic Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB [default: $params.memory]

        Workflow Options:
    --genus                     Bacterial genus (Escherichia, Salmonella, Enterobacter, Klebsiella, Staphylococcus)  [default: $params.genus]
    --species                   bacterial species to assemble (e.g. coli, pneumoniae, cloacae, aureus) [default: $params.species]
    --species_taxid             NCBI TaxID of the bacterial species to assemble [default: $params.species_taxid]

    --mash_dataset              path to mash dataset prepared for the species [default: $params.mash_dataset]
    --busco_lineage             to specify according to the bacterial species. e.g. enterobacterales_odb10, bacillales_odb10... check BUSCO [default: $params.busco_lineage]
    --busco_db_offline          path to BUSCO datasets if user wants to run BUSCO offline [default: params.$busco_db_offline]
    --amrfinder_organism        To specify for PointMutation detection if not Ecoli, Salmonella, Kpneumoniae or Saureus.
                                Can be among these: Acinetobacter_baumannii, Campylobacter,
                                Clostridioides_difficile, Enterococcus_faecalis, Enterococcus_faecium,
                                Escherichia, Klebsiella, Neisseria, Pseudomonas_aeruginosa,
                                Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius,
                                Streptococcus_agalactiae, Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae.
                                The amrfinderplus will be run if not specified, but no point mutations are detected.
                                [default: $params.amrfinder_organism]
    --conda_card_rgi            path to existing conda env for CARD RGI [default: $params.conda_card_rgi]
    --conda_amrfinder           path to existing conda env for AMRFinderPlus [default: $params.conda_amrfinder]
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
    include {quast; quast as quast2} from './modules/quast.nf' params(output: params.output)
    include {busco; busco as busco2} from './modules/busco.nf' params(output: params.output)
    include {busco_auto_prok; busco_auto_prok as busco_auto_prok2} from './modules/busco.nf' params(output: params.output)

    // including Kraken2 - nucleotide level
    include {kraken2nt_contigs} from './modules/kraken2.nf' params(output: params.output)
    include {extract_kraken} from './modules/kraken2.nf' params(output: params.output)

    // include mash_screen
    if (params.mash_dataset){
      include {mash_screen} from './modules/mash.nf' params(output: params.output)
    }
    // AMR analysis modules
    include {amrfinderplus; amrfinderplus as amrfinderplus2} from './modules/amrfinderplus.nf' params(output: params.output)
    include {amrfinderplus_no_species; amrfinderplus_no_species as amrfinderplus_no_species2} from './modules/amrfinderplus.nf' params(output: params.output)
    include {card_rgi; card_rgi as card_rgi2} from './modules/card.nf' params(output: params.output)
    //plasmids
    include {plasmidfinder; plasmidfinder as plasmidfinder2} from './modules/plasmidfinder.nf' params(output: params.output)
    // MLST
    include {mlst; mlst as mlst2} from './modules/mlst.nf' params(output: params.output)
    // Annotation
    include {prokka} from './modules/prokka.nf' params(output: params.output)
    include {bakta} from './modules/bakta.nf' params(output: params.output)


    //*************************************************
    // STEP 1 QC with fastp
    //*************************************************

    // DATA INPUT ILLUMINA
    illumina_input_ch = Channel
        .fromFilePairs( "${params.illumina}/*_{1,2}*.fq{,.gz}", checkIfExists: true)
        .view()

    // run fastp module
    fastp(illumina_input_ch)
    illumina_clean_ch = fastp.out[0]

    //*************************************************
    // STEP 2 - Assembly
    //*************************************************
    spades(illumina_clean_ch)
    contigs_ch = spades.out.assembly
    contigs_w_reads = spades.out.quast

    //*************************************************
    // STEP 2 - Assembly QC Quast, Busco on raw assembly
    //*************************************************
    // QUAST Assembly QC
    //no taxonomic decontamination of the contigs yet
    //deconta= "raw"
    quast(contigs_w_reads, "raw")

    // BUSCO completeness - Singularity container
    if(params.busco_lineage){
      busco(contigs_ch, params.busco_lineage, "raw", params.busco_db_offline)
    }
    else {
      busco_auto_prok(contigs_ch, "raw")
    }

    //*************************************************
    // STEP 3 - plasmids prediction on all raw contigs
    // and MLST on raw contigs
    //*************************************************
    if(params.plasmidfinder_db){
      plasmidfinder(contigs_ch, params.plasmidfinder_db, "raw")
    }

    mlst(contigs_ch, "raw")

    //*************************************************
    // STEP 4 - ARGs search: CARD RGI and AMRFinderPlus
    //*************************************************
    // on the all contigs (raw)
    // AMRFinderPlus NCBI
    // For all AMRFinderPlus processes
    if(params.genus == "Escherichia"){
      organism = 'Escherichia'
    }
    else if (params.genus == "Klebsiella"){
      organism = 'Klebsiella'
    }
    else if (params.genus == "Salmonella"){
      organism = 'Salmonella'
    }
    if (params.genus="Staphylococcus" && params.species == "aureus"){
      organism = 'Staphylococcus_aureus'
    }

    // if amrfinder_organism is given in the params directly
    if (params.amrfinder_organism){
      amrfinderplus(contigs_ch,params.amrfinder_organism, "raw")
    }
    else{
      if(params.genus){
        amrfinderplus(contigs_ch,organism, "raw")
      }
      else{
        amrfinderplus_no_species(contigs_ch, "raw")
      }
    }

    // CARD Resistance Genes Identifier
    if (params.conda_card_rgi){
      card_rgi(contigs_ch,params.card_db, "raw")
    }

    //*************************************************
    // STEP 5 - decontamination with Kraken2
    //*************************************************
    //kraken2nt contigs
    kraken2nt_contigs(contigs_w_reads, params.k2nt_db)
    krak_res = kraken2nt_contigs.out.kn_results
    krak_report = kraken2nt_contigs.out.kn_report
    contigs_kn2 = kraken2nt_contigs.out.kn_reads_contigs

    // KrakenTools
// Maybe I should treat all this with parameters
    if(params.species_taxid){
      sp_taxid = params.species_taxid
    }
    else if (params.genus=="Escherichia" && params.species == "coli") {
      sp_taxid = '562'
    }
    else if (params.genus=="Klebsiella" && params.species == "pneumoniae"){
      sp_taxid = '573'
    }
    else if (params.genus == "Salmonella"){
      sp_taxid = '590'
    }
    else if (params.genus=="Enterobacter" && params.species == "cloacae"){
      sp_taxid = '550'
    }
    else if (params.genus=="Staphylococcus" && params.species == "aureus"){
      sp_taxid = '1280'
    }
    else {
      exit 1, "No species or species taxid specified for retrieving the species of interest"
    }

    extract_kraken(contigs_kn2,krak_res,krak_report,sp_taxid,params.krakentools_extract)
//    deconta_contigs_ch = extract_kraken.out.kn_contigs_deconta
//    deconta_for_quast = extract_kraken.out.kn_reads_contigs_deconta

    deconta_contigs_ch = extract_kraken.out[0]
    deconta_for_quast = extract_kraken.out[1]

    if (deconta_contigs_ch) {
      //*************************************************
      // STEP 6 -Assembly QC (bis) of decontaminated contigs
      //*************************************************
      // QUAST Assembly QC

      quast2(deconta_for_quast,"deconta")

      // BUSCO completeness - Singularity container
      if(params.busco_lineage){
        busco2(deconta_contigs_ch, params.busco_lineage, "deconta",params.busco_db_offline)
      }
      else {
        busco_auto_prok2(deconta_contigs_ch, "deconta")
      }

      //*************************************************
      // STEP 7 - ARGs search: CARD RGI and AMRFinderPlus
      //*************************************************
      // on the deconta_contigs_ch
      // AMRFinderPlus NCBI

      if (params.amrfinder_organism){
        amrfinderplus2(deconta_contigs_ch,params.amrfinder_organism, "deconta")
      }
      else{
        if(params.genus){
          amrfinderplus2(deconta_contigs_ch,organism, "deconta")
        }
        else{
          amrfinderplus_no_species2(deconta_contigs_ch, "deconta")
        }
      }

      // CARD Resistance Genes Identifier
      if (params.conda_card_rgi){
        card_rgi2(deconta_contigs_ch,params.card_db, "deconta")
      }

      //*************************************************
      // STEP 8 -  annotation
      //*************************************************
      // bakta annotation of deconta contigs (and mapped contigs)
  //    if(params.bakta_db){
  //      bakta(deconta_contigs_ch, params.bakta_db, params.genus, params.species)
  //    }

      // prokka
      prokka(deconta_contigs_ch, params.genus, params.species)

      //*************************************************
      // STEP 9 - PlasmidFinder et al. Platon ? MGEFinder..
      //*************************************************
      if(params.plasmidfinder_db){
        plasmidfinder2(deconta_contigs_ch, params.plasmidfinder_db, "deconta")
      }

      //Platon
      //  PlasForest
      //  MOB-recon
      // on deconta_contigs_ch

      //*************************************************
      // STEP 10 - MLST - Sequence typing
      //*************************************************
      mlst2(deconta_contigs_ch, "deconta")

      //*************************************************
      // STEP 11 - Find closest relative with Mash
      //*************************************************
      // using the mash dataset provided
  /*    if(params.mash_dataset){
        mash_screen(deconta_contigs_ch,params.species,params.mash_sketch)
        // later, mash_screen will output fasta file of closest relative for mapping
      }
*/
      //*************************************************
      // STEP 12 - Map to the closest with Minimap2
      //*************************************************
      // minimap2 with PAF Output - use all contigs: contigs_ch
      // compare the number of contigs mapped with the number of deconta_contigs_ch

      // mapping could be done on same ref for all... need to select from mash results
      // might need a second pipeline for mappings and core/pan genomes analysis,
      // after the assemblies and decontamination are done...


    }


}
