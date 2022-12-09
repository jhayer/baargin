#!/usr/bin/env nextflow
nextflow.enable.dsl=2

start_var = Channel.from("""
********* Start running nf-wgs_amr pipeline *********
nf-wgs_amr is a workflow for bacterial genomics qc, assembly, decontamination by
taxonomic classification, analysis of Antimicrobial Resistance Genes, plasmids detection
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
    --contigs                   path to the directory containing the already assembled contigs files (fasta) (default: $params.contigs)
    --hybrid_index              path to the csv file containing the mapping between sampleID, illuminaR1.fastq.gz, illuminaR2.fastq.gz, ont_read.fastq

        Optional input:
    --phred_type                phred score type. Specify if 33 (default and current) or 64 (ex. BGI, older...) [default: $params.phred_type]
    --k2nt_db                   path to the Kraken2 nucleotide database (e.g. MiniKraken, nt) [default: $params.k2nt_db]
    --card_db                   path to the CARD json Database for Antimicrobial Resistance Genes prediction [default: $params.card_db]
    --plasmidfinder_db          path to the CGE PlasmidFinder database [default: $params.plasmidfinder_db]

      Output:
    --output                    path to the output directory (default: $params.output)
    --tmpdir                    path to the tmp directory (default: $params.tmpdir)

        Outputed directories:
    qc                          The reads file after qc, qc logs and host mapping logs
    assembly                    The spades assembly output directory and all metrics related to assembly (busco, quast)
    taxonomic_classif           The taxonomic classifications at contigs level
    AMR                         The output directory for resistane genes analysis: ARMFinderPlus and CARD

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
    --amrfinder_organism        To specify for PointMutation detection.
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
        workflow.profile.contains('singu')
    ) { "executer selected" }
    else { exit 1, "No executer selected: -profile itrop/singu"}


    //*************************************************
    // STEP 0 - Include needed modules
    //*************************************************

    include {fastp} from './modules/fastp.nf' params(output: params.output)
    include {fastp_hybrid} from './modules/fastp.nf' params(output: params.output)
    // including assembler module
    include {spades} from './modules/spades.nf' params(output: params.output)
    include {unicycler} from './modules/unicycler.nf' params(output: params.output)
    // include quast, busco, checkm
    include {quast} from './modules/quast.nf' params(output: params.output)
    include {quast_contigs_only; quast_contigs_only as quast_contigs_only2} from './modules/quast.nf' params(output: params.output)
    include {quast_hybrid; quast_hybrid as quast_hybrid2} from './modules/quast.nf' params(output: params.output)
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
    //include {bakta} from './modules/bakta.nf' params(output: params.output)


    //*************************************************
    // STEP 1 QC with fastp
    //*************************************************

    // DATA INPUT ILLUMINA
    if(params.illumina){
      illumina_input_ch = Channel
          .fromFilePairs( "${params.illumina}/*{1,2}*.fastq{,.gz}", checkIfExists: true)
          .view()

      // run fastp module
      fastp(illumina_input_ch, params.phred_type)
      illumina_clean_ch = fastp.out[0]

      //*************************************************
      // STEP 2 - Assembly
      //*************************************************
      spades(illumina_clean_ch)
      contigs_ch = spades.out.assembly
      contigs_w_reads = spades.out.quast

      //*************************************************
      // STEP 2bis - Assembly QC Quast, Busco on raw assembly
      //*************************************************
      // QUAST Assembly QC
      //no taxonomic decontamination of the contigs yet
      //deconta= "raw"
      quast(contigs_w_reads, "raw")

    }
    else if(params.contigs){
      // DATA INPUT is Contigs from assembly
      contigs_files_ch = Channel
        .fromPath("${params.contigs}/*.{fasta,fa}", checkIfExists: true)
        .view()

        contigs_files_ch.map { file ->
          def id = ( file.baseName.toString() =~ /^[^._]*(?=\_)/ )[0] // first element in the name until underscore
          return tuple(id, file)
        }
          .set{ contigs_ch}

      quast_contigs_only(contigs_ch, "raw")
    }
    else if(params.hybrid_index){
      // the input is CSV file mapping sampleID to illumina and ont reads files
      //sampleId	read1	read2 ont
      //FC816RLABXX reads/FC816RLABXX_L1_R1.fastq.gz reads/FC816RLABXX_L1_R2.fastq.gz ont/FC816RLABXX.fastq
      hybrid_ch = Channel.fromPath(params.hybrid_index, checkIfExists: true) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleID, file(row.read1), file(row.read2), file(row.ont)) }

      // run fastp module on short reads
      fastp_hybrid(hybrid_ch, params.phred_type)
      trimmed_hybrid_ch = fastp_hybrid.out.trimmed_hybrid

      unicycler(trimmed_hybrid_ch)
      contigs_ch = unicycler.out.assembly
      quast_hybrid(contigs_ch, hybrid_ch, "raw")
    }
    else{
      //no input provided
    }


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

    // if amrfinder_organism is given in the params directly
    if (params.amrfinder_organism){
      amrfinderplus(contigs_ch,params.amrfinder_organism, "raw")
    }
    else{
      amrfinderplus_no_species(contigs_ch, "raw")
    }

    // CARD Resistance Genes Identifier
    if (params.card_db){
      card_rgi(contigs_ch,params.card_db, "raw")
    }

    //*************************************************
    // STEP 5 - decontamination with Kraken2
    //*************************************************
    //kraken2nt contigs
    kraken2nt_contigs(contigs_ch, params.k2nt_db)
    krak_res = kraken2nt_contigs.out.kn_results
    krak_report = kraken2nt_contigs.out.kn_report
    contigs_kn2 = kraken2nt_contigs.out.kn_contigs

    // KrakenTools
    if(params.species_taxid){
      sp_taxid = params.species_taxid
    }
    else {
      exit 1, "No TaxID specified for retrieving the species of interest - please provide a NCBI TaxID (from species or genus) with --species_taxid"
    }

    extract_kraken(contigs_kn2,krak_res,krak_report,sp_taxid,params.krakentools_extract)
//    deconta_contigs_ch = extract_kraken.out.kn_contigs_deconta
//    deconta_for_quast = extract_kraken.out.kn_reads_contigs_deconta

    deconta_contigs_ch = extract_kraken.out[0]

    if (deconta_contigs_ch) {
      //*************************************************
      // STEP 6 -Assembly QC (bis) of decontaminated contigs
      //*************************************************
      // QUAST Assembly QC

      quast_contigs_only2(deconta_contigs_ch,"deconta")

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
        amrfinderplus_no_species2(deconta_contigs_ch, "deconta")
      }

      // CARD Resistance Genes Identifier
      if (params.card_db){
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
