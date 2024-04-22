#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

start_var = Channel.from("""
********* Start running baargin pipeline *********
baargin is a workflow for bacterial genomics qc, assembly, decontamination by
taxonomic classification, analysis of Antimicrobial Resistance Genes, plasmids detection
**************************************
""")
start_var.view()

if (params.help) { exit 0, helpMSG() }

// Help Message
def helpMSG() {
    log.info """
    ********* Bacterial Assembly and ARGs detection In Nextflow *********

        Usage example:
    nextflow run main.nf --reads_folder data --illumina_pattern "*R{1,2}_001_subs10000.fastq.gz" --genus Escherichia --species coli --species_taxid 562 -profile docker -resume
    --help                      prints the help section

        Input sequences:
    --illumina_pattern          pattern of the R1 and R2 illumina files paired. Ex: "*_{R1,R2}_001.fastq.gz" or "*_{1,2}.fastq.gz". Required with --reads_folder and must be quoted (default: $params.illumina_pattern)
    --reads_folder              path to the directory containing the illumina reads files (fastq.gz) (default: $params.reads_folder)
  OR
    --contigs                   path to the directory containing the already assembled contigs files (fasta) (default: $params.contigs)
  OR
    --hybrid_index              For users having both short and long reads:
                                path to the CSV file containing the mapping between sampleID, illuminaR1.fastq.gz, illuminaR2.fastq.gz, ont_read.fastq
                                Must have the header as follow:
                                sampleID,read1,read2,ont

                                Example of CSV index file:
                                sampleID,read1,read2,ont
                                124,test_illu_hybrid/124_1.fq,test_illu_hybrid/124_2.fq,test_ont/barcode05_concat.fastq
                                365,test_illu_hybrid/365_1.fq,test_illu_hybrid/365_2.fq,test_ont/barcode01_concat.fastq

        Output:
    --output                    path to the output directory (default: $params.output)
    --tmpdir                    path to the tmp directory (default: $params.tmpdir)

        Species mandatory options:
    --genus                     Bacterial genus (Escherichia, Salmonella, Enterobacter, Klebsiella, Staphylococcus)  [default: $params.genus]
    --species                   bacterial species to assemble - used for annotation (e.g. coli, pneumoniae, cloacae, aureus) [default: $params.species]
    --species_taxid             NCBI TaxID of the bacterial species to assemble - used for Kraken decontaminiation [default: $params.species_taxid]

        Databases path required (script provided for downloading them):
    --card_db                   path to the CARD json Database for Antimicrobial Resistance Genes prediction [default: $params.card_db]
    --kraken2_db                path to the local Kraken2 nucleotide database (e.g. MiniKraken, nt, standard) [default: $params.kraken2_db]
    --plasmidfinder_db          path to the CGE PlasmidFinder database [default: $params.plasmidfinder_db]

        Optional databases paths: if provided, the tool is run:
    --amrfinder_db              path to a local AMRFinder Database for Antimicrobial Resistance Genes prediction [default: $params.amrfinder_db] - a database is provided within the container

    --bakta_db                  path to the Bakta local database if the user prefers annotating the genomes with Bakta instead of Prokka [default: $params.bakta_db]
    --busco_db_offline          path to local BUSCO datasets if user wants to run BUSCO offline [default: $params.busco_db_offline]
    --platon_db                 path to the Platon local database

       Optional input:
    --phred_type                phred score type. Specify if 33 (default and current) or 64 (ex. BGI, older...) [default: $params.phred_type]
    --busco_lineage             to specify according to the bacterial species. e.g. enterobacterales_odb10, bacillales_odb10... check BUSCO [default: $params.busco_lineage]
                                If not provided, Busco will use prokaryotes database
    --amr_id_min                Minimum identity percentage for positive match with CARD and/or AMRFinder. (ex: 0.90, 0.85...) [default: $params.amr_id_min]
    --amr_cov_min               Minimum coverage of the reference protein for a positive match with CARD and/or AMRFinder. (ex: 0.60, 0.80...) [default: $params.amr_cov_min]
    --amrfinder_organism        To specify for PointMutation detection
                                Can be among these: Acinetobacter_baumannii, Campylobacter,
                                Clostridioides_difficile, Enterococcus_faecalis, Enterococcus_faecium,
                                Escherichia, Klebsiella, Neisseria, Pseudomonas_aeruginosa,
                                Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius,
                                Streptococcus_agalactiae, Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae.
                                The amrfinderplus will be run if not specified, but no point mutations are detected.
                                [default: $params.amrfinder_organism]
                                If not provided, resistance genes will be detected but not species-specific point mutations involved in AMR
        Nextflow options:
    -profile                    change the profile of nextflow both the engine and executor more details on github README
    -resume                     resume the workflow where it stopped

            Outputed directories:
    sample_ID
      AMR                       The output directory for resistance genes analysis: ARMFinderPlus and CARD
      annotation                The annotation directory containing Prokka or Bakta output if run
      assembly                  The spades assembly output directory with scaffolds.fasta files and all metrics related to assembly (busco, quast) and the "decontaminated" scaffolds
        |
         --taxonomic_classif    The taxonomic classification from Kraken2 at contigs/scaffolds level - the extracted scaffolds using provided TaxID are in the (upper) assembly directory
      plasmids                  The output directory for plasmids identification with PlasmidFinder and Platon
      qc                        The reads file after qc, qc logs and host mapping logs

    compile_results             The ouput directory for the summary files of all samples together, from all tools used. Presence/Absence (1/0) tabular (tsv) files
    pangenome                   The pangenome analysis output directory from Roary
    """
}

//*************************************************
// Check options
//*************************************************

// check plasmidfinder_db
if(params.plasmidfinder_db){
  File pfdb = new File("${params.plasmidfinder_db}");
  if(! pfdb.exists()){
    exit 1, "Mandatory plasmidfinder_db path does not exists! (${params.plasmidfinder_db})"   
  }
}
else {
  exit 1, "Mandatory plasmidfinder_db path missing!" 
}

// check card_db
if (params.card_db){
  File carddb = new File("${params.card_db}");
  if(carddb.exists()){
    if (carddb.isDirectory()){
      exit 1, "${params.card_db} is a folder, the card.json file is required here!" 
    }
  }
  else {
    exit 1, "Mandatory card_db path does not exists! (${params.card_db})"    
  }
}
else {
  exit 1, "Mandatory card_db path missing!" 
}

// check kraken2_db
if (params.kraken2_db){
  File kndb = new File("${params.kraken2_db}");
  if(! kndb.exists()){
    exit 1, "Mandatory kraken2_db path does not exists! (${params.kraken2_db})"
  }
}
else {
  exit 1, "Mandatory kraken2_db path missing!" 
}
//*************************************************
// STEP 0 - Include needed modules
//*************************************************

include {fastp} from './modules/fastp.nf'
include {fastp_hybrid} from './modules/fastp.nf'
// including assembler module
include {spades} from './modules/spades.nf'
include {unicycler} from './modules/unicycler.nf'
// include quast, busco, checkm
include {quast} from './modules/quast.nf'
include {quast_contigs_only as quast_contigs_only2; quast_contigs_only} from './modules/quast.nf'
include {quast_hybrid} from './modules/quast.nf'
include {compile_quast as compile_quast2; compile_quast} from './modules/quast.nf'
include {busco as busco2; busco} from './modules/busco.nf'
include {busco_proteins} from './modules/busco.nf'
include {compile_busco as compile_busco2; compile_busco} from './modules/busco.nf'
include {compile_busco_prot} from './modules/busco.nf'

// including Kraken2 - nucleotide level
include {kraken2nt_contigs} from './modules/kraken2.nf'
include {extract_kraken} from './modules/kraken2.nf'

// AMR analysis modules
include {amrfinderplus as amrfinderplus2; amrfinderplus} from './modules/amrfinderplus.nf'
include {amrfinderplus_no_db as amrfinderplus_no_db2; amrfinderplus_no_db} from './modules/amrfinderplus.nf'
include {amrfinderplus_no_species as amrfinderplus_no_species2; amrfinderplus_no_species} from './modules/amrfinderplus.nf'
include {amrfinderplus_no_species_no_db as amrfinderplus_no_species_no_db2; amrfinderplus_no_species_no_db} from './modules/amrfinderplus.nf'

include {card_rgi as card_rgi2; card_rgi} from './modules/card.nf'
//plasmids
include {plasmidfinder as plasmidfinder2; plasmidfinder} from './modules/plasmidfinder.nf'
include {platon as platon2; platon} from './modules/platon.nf'
include {platon_json2tsv as platon_json2tsv2; platon_json2tsv} from './modules/platon.nf'
include {compile_platon as compile_platon2; compile_platon} from './modules/platon.nf'
include {mefinder as mefinder2; mefinder} from './modules/mefinder.nf'

//split AMR into plasmids and chromosomes
include {split_amr_plas_chrom_amrfinder_sp as split_amr_plas_chrom_amrfinder_sp2; split_amr_plas_chrom_amrfinder_sp} from './modules/split_amr_plas_chrom.nf'
include {split_amr_plas_chrom_amrfinder_no_sp as split_amr_plas_chrom_amrfinder_no_sp2; split_amr_plas_chrom_amrfinder_no_sp} from './modules/split_amr_plas_chrom.nf'
include {split_card_plas_chrom as split_card_plas_chrom2; split_card_plas_chrom} from './modules/split_amr_plas_chrom.nf'

// MLST
include {mlst as mlst2; mlst} from './modules/mlst.nf'
include {compile_mlst as compile_mlst2; compile_mlst} from './modules/mlst.nf'

// Annotation
include {prokka; prokka_genus} from './modules/prokka.nf'
include {bakta; bakta_genus} from './modules/bakta.nf' 
// pangenome
include {roary} from './modules/roary.nf' params(output: params.output)
// compilation
include {compile_amrfinder as compile_amrfinder2; compile_amrfinder} from './modules/compile_amrfinder.nf'
include {compile_amrfinder_no_species as compile_amrfinder_no_species2; compile_amrfinder_no_species} from './modules/compile_amrfinder.nf'
include {compile_amrfinder_plasmid_split as compile_amrfinder_plasmid_split2; compile_amrfinder_plasmid_split} from './modules/compile_amrfinder.nf'
include {compile_amrfinder_no_sp_plasmid_split as compile_amrfinder_no_sp_plasmid_split2; compile_amrfinder_no_sp_plasmid_split} from './modules/compile_amrfinder.nf'

include {compile_plasmidfinder as compile_plasmidfinder2; compile_plasmidfinder} from './modules/compile_plasmidfinder.nf'
include {compile_card as compile_card2; compile_card} from './modules/card.nf'
include {compile_card_split_plasmid as compile_card_split_plasmid2; compile_card_split_plasmid} from './modules/compile_card_tsv.nf'



workflow {

  //  start_var.view()

    // error handling
    if (
        workflow.profile.contains('itrop') ||
        workflow.profile.contains('singularity') ||
        workflow.profile.contains('docker')
    ) { "executer selected" }
    else { exit 1, "No executer selected: -profile docker/singularity/itrop"}


    //*************************************************
    // STEP 1 QC with fastp
    //*************************************************
    samples_number = 0

    // DATA INPUT ILLUMINA
    if(params.reads_folder){

      // checking the number of samples in put
      def list_files = []
      File illumina_input = new File(params.reads_folder)
      if(illumina_input.exists()){
        if ( illumina_input.isDirectory()) {
          illumina_input.eachFileRecurse(FILES){
            if (it.name =~ ~/.*.f*q.*$/){
              list_files.add(it)
            }
          }
          nb_files = list_files.size()
          samples_number = nb_files/2
          log.info "${samples_number} samples in ${params.reads_folder}"
        }
        else {
          exit 1,  "The input ${params.reads_folder} is a file! A folder containing fastq(.gz) files pairs is expected!\n"
        }
      } else {
        exit 1, "The input folder ${params.reads_folder} does not exists!\n"
      }

      // creating channels from file pairs
      illumina_input_ch = Channel
          .fromFilePairs( "${params.reads_folder}/${params.illumina_pattern}", checkIfExists: true)
          .view()
          .ifEmpty { exit 1, "Cannot find any reads in the directory: ${params.reads_folder} matching the pattern ${params.illumina_pattern}\nNB: Path needs to be enclosed in quotes!" }


 /*     illumina_input_ch = Channel
        .fromFilePairs("${params.reads_folder}/${params.illumina_pattern}", size: 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.illumina}\nNB: Path needs to be enclosed in quotes!" }
        .map { row ->
          def meta = [:]
          meta.id = row[0]
          meta.group  = 0
          return [ meta, row[1] ] }
*/

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
      
      // checking the number of samples in put
      def list_files = []
      File fasta_input = new File(params.contigs)
      if(fasta_input.exists()){
        if ( fasta_input.isDirectory()) {
          fasta_input.eachFileRecurse(FILES){
            if (it.name =~ ~/.*.fasta/){
              list_files.add(it)
            }
          }
          samples_number = list_files.size()
          log.info "${samples_number} files in ${params.contigs}"
        }
        else {
          exit 1,  "The input ${params.contigs} is a file! A folder containing contigs in fasta is expected!\n"
        }
      } else {
        exit 1, "The input folder ${params.contigs} does not exists!\n"
      }

      // Creating contigs input channels from fasta files
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

      // checking number of samples
      // Count lines
      File hybrid_csv = new File(params.hybrid_index)
      if(hybrid_csv.exists()){
        if ( hybrid_csv.isDirectory()) {
          exit 1,  "The input ${params.hybrid_index} is a directory! A CSV file is expected!\n"
        }
        else {
          File file = new File(params.hybrid_index);
          LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(file));
          lineNumberReader.skip(Long.MAX_VALUE);
          int lines = lineNumberReader.getLineNumber();
          lineNumberReader.close();
          samples_number = lines - 1
          log.info "${samples_number} samples in ${params.hybrid_index}"
        }
      } else {
        exit 1, "The input CSV file ${params.hybrid_index} does not exists!\n"
      }
      // Creating the input hybrid channel containing illumina paired reads and ont reads
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
  //    quast_hybrid(contigs_ch, hybrid_ch, "raw")
      quast_hybrid(contigs_ch, trimmed_hybrid_ch, "raw")
      quast_collect = quast_hybrid.out.quast_transpo.collect()
    }
    else{
      //no input provided
      exit 1, "No input specified. Please provide input short reads with --illumina, \
        or input contigs with --contigs, \
        or provide a CSV sample sheet for hybrid: short and long reads with --hybrid_index"
    }

    //compile Quast results
    compile_quast(quast_collect, "raw")

    //*************************************************
    // STEP 2bis - Assembly QC Busco on raw assembly
    //*************************************************
    // BUSCO completeness - Singularity container
    add_cmd = ""
    if(params.busco_db_offline){
      File bu_offline = new File("${params.busco_db_offline}")
      if(bu_offline.exists()){
        add_cmd = add_cmd + " --offline --download_path ${params.busco_db_offline}"
      }
      else{
        exit 1, "${params.busco_db_offline} busco_db_offline path does not exists!"
      }
    }
    if(params.busco_lineage){
      add_cmd = add_cmd + " --lineage_dataset ${params.busco_lineage}"
    }
    else {
      add_cmd = add_cmd + " --auto-lineage-prok"
    }
    log.info "Here is the busco command line addition=${add_cmd}"
    busco(contigs_ch, "raw", add_cmd)
    // compile and format Busco results
    compile_busco(busco.out.busco_sum.collect(), "raw")

    //*************************************************
    // STEP 3 - plasmids prediction on all raw contigs
    //*************************************************
    
    plasmidfinder(contigs_ch, params.plasmidfinder_db, "raw")
    compile_plasmidfinder(plasmidfinder.out.plasmidfinder_tab.collect(), "raw")

    if(params.platon_db){
      File platondb = new File("${params.platon_db}");
      if(platondb.exists()){
        platon(contigs_ch, params.platon_db, "raw")
        if(platon.out.platon_json){
          platon_json2tsv(platon.out.platon_json, "raw")
          compile_platon(platon_json2tsv.out.platon_inc.collect(), platon_json2tsv.out.platon_plasmid.collect(), platon_json2tsv.out.platon_amr.collect(), "raw" )
        }
        else {
          log.info "Platon did not identify any plamid sequence"
        }

      }
      else {
        exit 1, "${params.platon_db} platon_db path does not exists!"
      }
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

    
      File amrdb = new File("${params.amrfinder_db}");
      if (amrdb.exists()){
        // if amrfinder_organism is given in the params directly
        if (params.amrfinder_organism){
          amrfinderplus(contigs_ch,params.amrfinder_organism, params.amrfinder_db, "raw", params.amr_id_min, params.amr_cov_min)

          if(platon.out.tp_platon_id_tsv){
            // split sp
            amrfinderplus.out.tp_id_amrf.join(platon.out.tp_platon_id_tsv).set{tp_id_amrf_platon}
            split_amr_plas_chrom_amrfinder_sp(tp_id_amrf_platon, "raw")
            // compile plas and chrom
            compile_amrfinder_plasmid_split(split_amr_plas_chrom_amrfinder_sp.out.amrf_plasmid.collect(),
                split_amr_plas_chrom_amrfinder_sp.out.amrf_allmut_plasmid.collect(),
                split_amr_plas_chrom_amrfinder_sp.out.amrf_chrom.collect(),
                split_amr_plas_chrom_amrfinder_sp.out.amrf_allmut_chrom.collect(), "raw")
          }
          else{
            compile_amrfinder(amrfinderplus.out.amrfile.collect(), amrfinderplus.out.amrfile_allmut.collect(), "raw")
          }

      //    compile_amrfinder(amrfinderplus.out.amrfile.collect(), amrfinderplus.out.amrfile_allmut.collect(), "raw")
        }
        else{
          amrfinderplus_no_species(contigs_ch, params.amrfinder_db, "raw", params.amr_id_min, params.amr_cov_min)
          // if platon results, split the amr file into plasmids and chrom
          if(platon.out.tp_platon_id_tsv){
            amrfinderplus_no_species.out.tp_id_amrf.join(platon.out.tp_platon_id_tsv).set{tp_id_amrf_platon}
            // split AMRFinder files
            split_amr_plas_chrom_amrfinder_no_sp(tp_id_amrf_platon, "raw")
            // compile plas and chrom
            compile_amrfinder_no_sp_plasmid_split(split_amr_plas_chrom_amrfinder_no_sp.out.amrf_plasmid.collect(),
                split_amr_plas_chrom_amrfinder_no_sp.out.amrf_chrom.collect(), "raw")
          }
          else{
            compile_amrfinder_no_species(amrfinderplus_no_species.out.amrfile.collect(), "raw")
          }
        }
      }
      else {
        // if amrfinder_organism is given in the params directly
        if (params.amrfinder_organism){
          amrfinderplus_no_db(contigs_ch,params.amrfinder_organism, "raw", params.amr_id_min, params.amr_cov_min)
          // if platon results, split the amr files into plasmids and chrom
          if(platon.out.tp_platon_id_tsv){
            // make joint channel platon file + amrfinder file
            amrfinderplus_no_db.out.tp_id_amrf.join(platon.out.tp_platon_id_tsv).set{tp_id_amrf_platon}
            // split AMRFinder files
            split_amr_plas_chrom_amrfinder_sp(tp_id_amrf_platon, "raw")
            // compile plas and chrom separately
            compile_amrfinder_plasmid_split(split_amr_plas_chrom_amrfinder_sp.out.amrf_plasmid.collect(),
                split_amr_plas_chrom_amrfinder_sp.out.amrf_allmut_plasmid.collect(),
                split_amr_plas_chrom_amrfinder_sp.out.amrf_chrom.collect(),
                split_amr_plas_chrom_amrfinder_sp.out.amrf_allmut_chrom.collect(), "raw")
          }
          else{
            compile_amrfinder(amrfinderplus_no_db.out.amrfile.collect(), amrfinderplus_no_db.out.amrfile_allmut.collect(), "raw")
          }
        }
        else{
          amrfinderplus_no_species_no_db(contigs_ch, "raw", params.amr_id_min, params.amr_cov_min)
          // if platon results, split the amr file into plasmids and chrom
          if(platon.out.tp_platon_id_tsv){
            // make joint channel platon file + amrfinder file
            amrfinderplus_no_species_no_db.out.tp_id_amrf.join(platon.out.tp_platon_id_tsv).set{tp_id_amrf_platon}
            // split AMRFinder files
            split_amr_plas_chrom_amrfinder_no_sp(tp_id_amrf_platon, "raw")
            // compile plas and chrom
            compile_amrfinder_no_sp_plasmid_split(split_amr_plas_chrom_amrfinder_no_sp.out.amrf_plasmid.collect(),
                split_amr_plas_chrom_amrfinder_no_sp.out.amrf_chrom.collect(), "raw")

          }
          else {
            compile_amrfinder_no_species(amrfinderplus_no_species_no_db.out.amrfile.collect(), "raw")
          }
          
        }

      }
    

    // CARD Resistance Genes Identifier
    
    card_rgi(contigs_ch,params.card_db, "raw", params.amr_id_min, params.amr_cov_min)
    if(platon.out.tp_platon_id_tsv){
      // split RGI_main.txt into plasmid and chromosome
      card_rgi.out.tp_id_card.join(platon.out.tp_platon_id_tsv).set{tp_id_card_platon}
      // split AMRFinder files
      split_card_plas_chrom(tp_id_card_platon, "raw")
      // compile with compile_card.py and not with rgi heatmap which is based on the json and not the tsv
      compile_card_split_plasmid(split_card_plas_chrom.out.rgi_plasmid.collect(), split_card_plas_chrom.out.rgi_chrom.collect(), "raw")
    }
    else{
      //compile all in one
      if (samples_number > 1) {
        compile_card(card_rgi.out.card_json.collect(), "raw")
      }
    }
 //   if (samples_number > 1) {
 //     compile_card(card_rgi.out.card_json.collect(), "raw")
 //   }

    //*************************************************
    // STEP 5 - decontamination with Kraken2
    //*************************************************
    //kraken2nt contigs
    kraken2nt_contigs(contigs_ch, params.kraken2_db)
    krak_res = kraken2nt_contigs.out.kn_results
    krak_report = kraken2nt_contigs.out.kn_report
    contigs_kn2 = kraken2nt_contigs.out.kn_contigs


    // KrakenTools
    if(params.species_taxid){
      extract_kraken(contigs_kn2,krak_res,krak_report,params.species_taxid)
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

      // BUSCO completeness
      add_cmd2 = ""
      if(params.busco_db_offline){
        File local_bu_db = new File("${params.busco_db_offline}")
        if(local_bu_db.exists()){
          add_cmd2 = add_cmd2 + " --offline --download_path ${params.busco_db_offline}"
        }
        else{
          exit 1, "${params.busco_db_offline} busco_db_offline path does not exists!"
        }
      }
      if(params.busco_lineage){
        add_cmd2 = add_cmd2 + " --lineage_dataset ${params.busco_lineage}"
      }
      else {
        add_cmd2 = add_cmd2 + " --auto-lineage-prok"
      }
      busco2(deconta_contigs_ch, "deconta", add_cmd2)
      // compile and format Busco results
      compile_busco2(busco2.out.busco_sum.collect(), "deconta")

      //*************************************************
      // STEP 7 - PlasmidFinder and Platon
      //*************************************************
      // PlasmidFinder
      if(params.plasmidfinder_db){
        plasmidfinder2(deconta_contigs_ch, params.plasmidfinder_db, "deconta")
        compile_plasmidfinder2(plasmidfinder2.out.plasmidfinder_tab.collect(), "deconta")
      }
      //Platon
      if(params.platon_db){
          platon2(deconta_contigs_ch, params.platon_db, "deconta")
          if(platon2.out.platon_json){
            platon_json2tsv2(platon2.out.platon_json, "deconta")
            compile_platon2(platon_json2tsv2.out.platon_inc.collect(), platon_json2tsv2.out.platon_plasmid.collect(), platon_json2tsv2.out.platon_amr.collect(), "deconta" )
          }
          else {
            log.info "Platon did not identify any plamid sequence"
          }
      }

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

      
        File amrdb2 = new File("${params.amrfinder_db}");
        if (amrdb2.exists()) {
          // if amrfinder_organism is given in the params directly, run the mutation search
          if (params.amrfinder_organism){
            amrfinderplus2(deconta_contigs_ch,params.amrfinder_organism, params.amrfinder_db, "deconta", params.amr_id_min, params.amr_cov_min)
            compile_amrfinder2(amrfinderplus2.out.amrfile.collect(), amrfinderplus2.out.amrfile_allmut.collect(), "deconta")
          }
          else {
            // no mutation search
            amrfinderplus_no_species2(deconta_contigs_ch, params.amrfinder_db, "deconta", params.amr_id_min, params.amr_cov_min)
            compile_amrfinder_no_species2(amrfinderplus_no_species2.out.amrfile.collect(), "deconta")
          }
        }
        else {
          //run AMRFinder with no localDB - will use the container DB
          // if amrfinder_organism is given in the params directly, run the mutation search
          if (params.amrfinder_organism){
            amrfinderplus_no_db2(deconta_contigs_ch,params.amrfinder_organism, "deconta", params.amr_id_min, params.amr_cov_min)
            compile_amrfinder2(amrfinderplus_no_db2.out.amrfile.collect(), amrfinderplus_no_db2.out.amrfile_allmut.collect(), "deconta")
          }
          else {
            // no mutation search
            amrfinderplus_no_species_no_db2(deconta_contigs_ch, "deconta", params.amr_id_min, params.amr_cov_min)
            compile_amrfinder_no_species2(amrfinderplus_no_species_no_db2.out.amrfile.collect(), "deconta")
          }
        }
      

      // CARD Resistance Genes Identifier
      if (params.card_db){
        card_rgi2(deconta_contigs_ch,params.card_db, "deconta", params.amr_id_min, params.amr_cov_min)
        if (samples_number > 1){
          compile_card2(card_rgi2.out.card_json.collect(), "deconta")
        }
      }

      //*************************************************
      // STEP 10 -  annotation and pangenome
      //*************************************************
      if(params.genus){
        log.info "Genus is: ${params.genus}"
      }
      else{
        exit 1, "Genus not provided, genus is needed for annotation!"
      }
      if(params.bakta_db){
        File baktadb = new File("${params.bakta_db}");
        if (baktadb.exists()) {
          // Annotation with Batka if DB provided
          if(params.species){
            bakta(deconta_contigs_ch, params.bakta_db, params.genus, params.species)
            // proteins for Busco
            faa_annot = bakta.out.annot_faa
            // pangenome analysis with Roary using all gff outputs from bakta
            gff_annot = bakta.out.annot_gff
          }
          else{
            bakta_genus(deconta_contigs_ch, params.bakta_db, params.genus)
            // proteins for Busco
            faa_annot = bakta_genus.out.annot_faa
            // pangenome analysis with Roary using all gff outputs from bakta
            gff_annot = bakta_genus.out.annot_gff
          }
        }
        else {
          exit 1, "${params.bakta_db} bakta_db path does not exists!"
        }
      }
      else { // if no Bakta DB provided, run Prokka as default
        // prokka
        if(params.species){
          prokka(deconta_contigs_ch, params.genus, params.species)
          // proteins for Busco
          faa_annot = prokka.out.annot_faa
          // GFF for pangenome analysis with Roary
          gff_annot = prokka.out.annot_gff
        }
        else{
          prokka_genus(deconta_contigs_ch, params.genus)
          // proteins for Busco
          faa_annot = prokka_genus.out.annot_faa
          // GFF for pangenome analysis with Roary
          gff_annot = prokka_genus.out.annot_gff
        }


      }

      // Busco on annotation
      add_p_cmd = ""
      if(params.busco_db_offline){
        File local_bup_db = new File("${params.busco_db_offline}")
        if(local_bup_db.exists()){
          add_p_cmd = add_p_cmd + " --offline --download_path ${params.busco_db_offline}"
        }
        else{
          exit 1, "${params.busco_db_offline} busco_db_offline path does not exists!"
        }
      }
      if(params.busco_lineage){
        add_p_cmd = add_p_cmd + " --lineage_dataset ${params.busco_lineage}"
      }
      else {
        add_p_cmd = add_p_cmd + " --auto-lineage-prok"
      }
      busco_proteins(faa_annot, add_p_cmd)
      compile_busco_prot(busco_proteins.out.busco_prot.collect())

      //*************************************************
      // STEP 11 -  pangenome with Roary
      //*************************************************

      // pangenome analysis with Roary using all gff outputs from Prokka or Bakta
      if (samples_number > 1){
        roary(gff_annot.collect())
      }

    }

}
