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
    --bakta_db                  path to the Bakta local database if the use prefer annotating the genomes with Bakta instead of Prokka [default: $params.bakta_db]
    --card_db                   path to the CARD json Database for Antimicrobial Resistance Genes prediction [default: $params.card_db]
    --amrfinder_db              path to a local AMRFinder Database for Antimicrobial Resistance Genes prediction [default: $params.amrfinder_db]
    --plasmidfinder_db          path to the CGE PlasmidFinder database [default: $params.plasmidfinder_db]

      Output:
    --output                    path to the output directory (default: $params.output)
    --tmpdir                    path to the tmp directory (default: $params.tmpdir)

        Outputed directories:
    qc                          The reads file after qc, qc logs and host mapping logs
    assembly                    The spades assembly output directory and all metrics related to assembly (busco, quast)
    taxonomic_classif           The taxonomic classifications at contigs level
    AMR                         The output directory for resistane genes analysis: ARMFinderPlus and CARD

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
include {quast_hybrid} from './modules/quast.nf' params(output: params.output)
include {compile_quast; compile_quast as compile_quast2} from './modules/quast.nf' params(output: params.output)
include {busco; busco as busco2} from './modules/busco.nf' params(output: params.output)
include {busco_auto_prok; busco_auto_prok as busco_auto_prok2} from './modules/busco.nf' params(output: params.output)
include {busco_proteins; busco_proteins_auto_prok} from './modules/busco.nf' params(output: params.output)
include {compile_busco; compile_busco as compile_busco2} from './modules/busco.nf' params(output: params.output)
include {compile_busco_prot} from './modules/busco.nf' params(output: params.output)

// including Kraken2 - nucleotide level
include {kraken2nt_contigs} from './modules/kraken2.nf' params(output: params.output)
include {extract_kraken} from './modules/kraken2.nf' params(output: params.output)

// AMR analysis modules
include {amrfinderplus; amrfinderplus as amrfinderplus2} from './modules/amrfinderplus.nf' params(output: params.output)
include {amrfinderplus_no_species; amrfinderplus_no_species as amrfinderplus_no_species2} from './modules/amrfinderplus.nf' params(output: params.output)
include {card_rgi; card_rgi as card_rgi2} from './modules/card.nf' params(output: params.output)
//plasmids
include {plasmidfinder; plasmidfinder as plasmidfinder2} from './modules/plasmidfinder.nf' params(output: params.output)
include {platon; platon as platon2} from './modules/platon.nf' params(output: params.output)
include {mefinder; mefinder as mefinder2} from './modules/mefinder.nf' params(output: params.output)

// MLST
include {mlst; mlst as mlst2} from './modules/mlst.nf' params(output: params.output)
include {compile_mlst; compile_mlst as compile_mlst2 } from './modules/mlst.nf' params(output: params.output)

// Annotation
include {prokka} from './modules/prokka.nf' params(output: params.output)
include {bakta} from './modules/bakta.nf' params(output: params.output)
// pangenome
include {roary} from './modules/roary.nf' params(output: params.output)
// compilation
include {compile_amrfinder; compile_amrfinder as compile_amrfinder2} from './modules/compile_amrfinder.nf'
include {compile_amrfinder_no_species; compile_amrfinder_no_species as compile_amrfinder_no_species2} from './modules/compile_amrfinder.nf'
include {compile_plasmidfinder; compile_plasmidfinder as compile_plasmidfinder2} from './modules/compile_plasmidfinder.nf'
include {compile_card; compile_card as compile_card2} from './modules/card.nf'



workflow {

    // error handling
    if (
        workflow.profile.contains('itrop') ||
        workflow.profile.contains('singu')
    ) { "executer selected" }
    else { exit 1, "No executer selected: -profile itrop/singu"}


    //*************************************************
    // STEP 1 QC with fastp
    //*************************************************

    // DATA INPUT ILLUMINA
    if(params.illumina){
      illumina_input_ch = Channel
          .fromFilePairs( "${params.illumina}/*R{1,2}*.fastq{,.gz}", checkIfExists: true)
          .view()
          .ifEmpty { exit 1, "Cannot find any reads in the directory: ${params.illumina}" }

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
      // STEP 2bis - Assembly QC Quast on raw assembly
      //*************************************************
      // QUAST Assembly QC
      //no taxonomic decontamination of the contigs yet
      //deconta= "raw"
      quast(contigs_w_reads, "raw")
      quast_collect = quast.out.quast_transpo.collect()
    }
    else if(params.contigs){
      // DATA INPUT is Contigs from assembly
      contigs_files_ch = Channel
        .fromPath("${params.contigs}/*.{fasta,fa}", checkIfExists: true)
        .view()
        .ifEmpty { exit 1, "Cannot find any contigs in the directory: ${params.contigs}" }

        contigs_files_ch.map { file ->
        //  def id = ( file.baseName.toString() =~ /^[^._]*(?=\_)/ )[0] // first element in the name until underscore
          def id = ( file.baseName.toString())
          return tuple(id, file)
        }
          .set{ contigs_ch}
      //*************************************************
      // STEP 2bis - Assembly QC Quast on raw assembly
      //*************************************************
      quast_contigs_only(contigs_ch, "raw")
      quast_collect = quast_contigs_only.out.quast_transpo.collect()
    }
    else if(params.hybrid_index){
      // the input is CSV file mapping sampleID to illumina and ont reads files
      //sampleId	read1	read2 ont
      //FC816RLABXX reads/FC816RLABXX_L1_R1.fastq.gz reads/FC816RLABXX_L1_R2.fastq.gz ont/FC816RLABXX.fastq
      hybrid_ch = Channel.fromPath(params.hybrid_index, checkIfExists: true) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleID, file(row.read1), file(row.read2), file(row.ont)) }

      //*************************************************
      // STEP 1 QC with fastp
      //*************************************************
      // run fastp module on short reads
      fastp_hybrid(hybrid_ch, params.phred_type)
      trimmed_hybrid_ch = fastp_hybrid.out.trimmed_hybrid
      //*************************************************
      // STEP 2 - Hybrid Assembly
      //*************************************************
      unicycler(trimmed_hybrid_ch)
      contigs_ch = unicycler.out.assembly
      //*************************************************
      // STEP 2bis - Assembly QC Quast on raw assembly
      //*************************************************
      quast_hybrid(contigs_ch, hybrid_ch, "raw")
      quast_collect = quast_hybrid.out.quast_transpo.collect()
    }
    else{
      //no input provided
      exit 1, "No input specified. Please provide input short reads with --illumina, \
        or input contigs with --contigs, \
        or provide a sample sheet for hybrid short and long reads with --hybrid_index"
    }

    //compile Quast results
    compile_quast(quast_collect, "raw")

    //*************************************************
    // STEP 2bis - Assembly QC Busco on raw assembly
    //*************************************************
    // BUSCO completeness - Singularity container
    if(params.busco_lineage){
      busco(contigs_ch, params.busco_lineage, "raw", params.busco_db_offline)
      busco_collect = busco.out.busco_sum.collect()
    }
    else {
      busco_auto_prok(contigs_ch, "raw", params.busco_db_offline)
      busco_collect = busco_auto_prok.out.busco_sum.collect()
    }
    // compile and format Busco results
    compile_busco(busco_collect, "raw")

    //*************************************************
    // STEP 3 - plasmids prediction on all raw contigs
    //*************************************************
    if(params.plasmidfinder_db){
      plasmidfinder(contigs_ch, params.plasmidfinder_db, "raw")
      compile_plasmidfinder(plasmidfinder.out.plasmidfinder_tab.collect(), "raw")
    }
    if(params.platon_db){
      platon(contigs_ch, params.platon_db, "raw")
    }

  //  mefinder(contigs_ch, "raw")
    //*************************************************
    // STEP 3bis - MLST on raw contigs
    //*************************************************
    mlst(contigs_ch, "raw")
    compile_mlst(mlst.out.mlst_tab.collect(), "raw")

    //*************************************************
    // STEP 4 - ARGs search: CARD RGI and AMRFinderPlus
    //*************************************************
    // on the all contigs (raw)
    // AMRFinderPlus NCBI

    // if amrfinder_organism is given in the params directly
    if (params.amrfinder_organism){
      amrfinderplus(contigs_ch,params.amrfinder_organism, params.amrfinder_db, "raw")
      compile_amrfinder(amrfinderplus.out.amrfile.collect(), amrfinderplus.out.amrfile_allmut.collect(), "raw")
    }
    else{
      amrfinderplus_no_species(contigs_ch, params.amrfinder_db, "raw")
      compile_amrfinder_no_species(compile_amrfinder_no_species.out.amrfile.collect(), "raw")
    }

    // CARD Resistance Genes Identifier
    if (params.card_db){
      card_rgi(contigs_ch,params.card_db, "raw")
      compile_card(card_rgi.out.card_json.collect(), "raw")
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
      extract_kraken(contigs_kn2,krak_res,krak_report,params.species_taxid,params.krakentools_extract)
      deconta_contigs_ch = extract_kraken.out[0]
    }
    else {
      exit 1, "No TaxID specified for retrieving the species of interest - please provide a NCBI TaxID (from species or genus) with --species_taxid"
    }

    if (deconta_contigs_ch) {
      //*************************************************
      // STEP 6 -Assembly QC (bis) of decontaminated contigs
      //*************************************************
      // QUAST Assembly QC
      quast_contigs_only2(deconta_contigs_ch,"deconta")
      //compile quast deconta results
      compile_quast2(quast_contigs_only2.out.quast_transpo.collect(), "deconta")

      // BUSCO completeness - Singularity container
      if(params.busco_lineage){
        busco2(deconta_contigs_ch, params.busco_lineage, "deconta", params.busco_db_offline)
        busco_collect2 = busco2.out.busco_sum.collect()
      }
      else {
        busco_auto_prok2(deconta_contigs_ch, "deconta", params.busco_db_offline)
        busco_collect2 = busco_auto_prok2.out.busco_sum.collect()
      }
      // compile and format Busco results
      compile_busco2(busco_collect2, "deconta")

      //*************************************************
      // STEP 7 - PlasmidFinder et al. Platon ? MGEFinder..
      //*************************************************
      // Platon
      if(params.plasmidfinder_db){
        plasmidfinder2(deconta_contigs_ch, params.plasmidfinder_db, "deconta")
        compile_plasmidfinder2(plasmidfinder2.out.plasmidfinder_tab.collect(), "deconta")
      }
      //Platon
      if(params.platon_db){
        platon2(deconta_contigs_ch, params.platon_db, "deconta")
      }

  //    mefinder2(contigs_ch, "deconta")

      //*************************************************
      // STEP 8 - MLST - Sequence typing
      //*************************************************
      mlst2(deconta_contigs_ch, "deconta")
      compile_mlst2(mlst2.out.mlst_tab.collect(), "deconta")

      //*************************************************
      // STEP 9 - ARGs search: CARD RGI and AMRFinderPlus
      //*************************************************
      // on the deconta_contigs_ch
      // AMRFinderPlus NCBI

      if (params.amrfinder_organism){
        amrfinderplus2(deconta_contigs_ch,params.amrfinder_organism, params.amrfinder_db, "deconta")
        compile_amrfinder2(amrfinderplus2.out.amrfile.collect(), amrfinderplus2.out.amrfile_allmut.collect(), "deconta")
      }
      else{
        amrfinderplus_no_species2(deconta_contigs_ch, params.amrfinder_db, "deconta")
        compile_amrfinder_no_species2(amrfinderplus_no_species2.out.amrfile.collect(), "deconta")
      }

      // CARD Resistance Genes Identifier
      if (params.card_db){
        card_rgi2(deconta_contigs_ch,params.card_db, "deconta")
        compile_card2(card_rgi2.out.card_json.collect(), "deconta")
      }

      //*************************************************
      // STEP 10 -  annotation and pangenome
      //*************************************************
      if(params.bakta_db){
        // Annotation with Batka if DB provided
        bakta(deconta_contigs_ch, params.bakta_db, params.genus, params.species)
        // proteins for Busco
        faa_annot = bakta.out.annot_faa
      
        // pangenome analysis with Roary using all gff outputs from bakta
        gff_annot = bakta.out.annot_gff
      }
      else { // if no Bakta DB provided, run Prokka as default
        // prokka
        prokka(deconta_contigs_ch, params.genus, params.species)
        // proteins for Busco
        faa_annot = prokka.out.annot_faa

        // GFF for pangenome analysis with Roary
        gff_annot = prokka.out.annot_gff
        
      }

      // Busco on annotation
      if(params.busco_lineage){
         busco_proteins(faa_annot, params.busco_lineage,params.busco_db_offline)
         busco_prot_collect = busco_proteins.out.busco_prot.collect()
      }
      else {
        busco_proteins_auto_prok(faa_annot, params.busco_db_offline)
        busco_prot_collect = busco_proteins_auto_prok.out.busco_prot.collect()
      }
      compile_busco_prot(busco_prot_collect)

      //*************************************************
      // STEP 11 -  pangenome with Roary
      //*************************************************
//      roary(prokka.out.prokka_gff.collect())
      // pangenome analysis with Roary using all gff outputs from Prokka or Bakta
      roary(gff_annot.collect())

    }

}
