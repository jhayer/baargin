# baargin: Bacterial Assembly and Antimicrobial Resistance Genes In NextFlow
Workflow for analysis of Whole Genome Sequencing (WGS) data with AntiMicrobial Resistance (AMR) focus

<img src="doc/img/IRD.png" width="300" height="100" /> <img src="doc/img/MIVEGEC.png" width="150" height="100" />


## Table of Contents

   * [Prerequisites](#prerequisites)
   * [Installation](#installation)
      * [Using conda](#using-conda)
      * [Old school - Manually](#old-school---manually)
   * [Download databases](#download-databases)
     * [Mandatory databases](#mandatory-databases)
     * [Optional databases](#optional-databases)
   * [Test the workflow](#Test)
   * [Usage](#usage)
   * [Parameters](#parameters)
   * [Update](#update)
   * [Uninstall](#uninstall)
   * [Complete help and options](#complete)
   * [Citation](#citation)
   * [Author](#author-and-contributors)


## Prerequisites

You need to have installed Docker or Singularity as the workflow uses containers to run the different tools.


## Installation

### Using conda

Prerequisite: conda 

   <details>
      <summary>See here for conda installation</summary>
  
      Conda is installed by downloading and executing an installer from the Conda website, but which version you need depends on your operating system:

      ```
      # Install Miniconda3 for 64-bit Mac
      curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-MacOSX-x86_64.sh -O
      bash Miniconda3-4.7.12.1-MacOSX-x86_64.sh
      rm Miniconda3-4.7.12.1-MacOSX-x86_64.sh
      ```
      ```
      # Install Miniconda3 for 64-bit Linux
      curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O
      bash Miniconda3-4.7.12.1-Linux-x86_64.sh
      rm Miniconda3-4.7.12.1-Linux-x86_64.sh
      ```
      The installer will ask you questions during the installation:

      Do you accept the license terms? (Yes)
      Do you accept the installation path or do you want to choose a different one? (Probably yes)
      Do you want to run conda init to setup Conda on your system? (Yes)
      Restart your shell so that the settings in ~/.bashrc/~/.bash_profile can take effect. You can verify that the installation worked by running:

      ```
      conda --version
      ```
      Lastly, we will setup the default channels (from where packages will be searched for and downloaded if no channel is specified).

      ```
      conda config --add channels defaults
      conda config --add channels bioconda
      conda config --add channels conda-forge
      ```
   </details>

Clone the baargin repository and move into it
```
git clone https://github.com/jhayer/baargin.git
cd 
```

Create an environment conda with needed dependencies:
```
conda env create -f conda_environment_baargin.yml
```

Activate the environment to be ready to use baargin:
```
conda activate baargin
```

### Old school - Manually

Python modules for the setup and installation of minimal databases. 
Prerequisite: Python>=3.8.0  
```
pip install pyyaml gitpython requests biopython>=1.78 numpy>=1.22
```

Then, you need to install NextFlow version 22.04.0+ [https://www.nextflow.io](https://www.nextflow.io)

```
# verify Java version (at least version 8+)
java -version

# Setup nextflow (it will create a nextflow executable file in the current directory)
curl -s https://get.nextflow.io | bash

# Then clone the workflow repository
git clone https://github.com/jhayer/baargin.git
```
Then the repository folder appears in your local directory


## Download databases

### Mandatory databases

Go to the repository folder and run:
```
# download the mandatory databases to run the workflow
download_db.py
```

### Optional databases

Some databases are not mandatory can take a significant disk space. We do not provide a download script for those but they can be installed separately by the user who will then provide the path accordingly when running the workflow.

Please note that if path for these databases are not provided, the corresponding tool is not run by the workflow.

1. Platon database - for plasmids detection

Can take approx. 2.8G.
If you wish to download it, please visit:https://github.com/oschwengers/platon#database

The database can be downloaded without installing bakta. Just `curl` or `wget` the URL they provide.

2. Bakta database - for genome structural and functionnal annotation

Can take around 62 Gb.
If you wish to download it, please visit: [https://github.com/oschwengers/bakta](https://github.com/oschwengers/bakta#database)

If you do not want to install bakta for downloading its database, they provide a download link from Zenodo that you can directly `wget` or `curl` and then decompress. Example:

```
wget https://zenodo.org/record/7025248/files/db.tar.gz
tar -xzf db.tar.gz
```

If Bakta database is provided, the annotation will be performed by Bakta, otherwise Prokka will be used with its default database.


3. AMRFinderPlus => add if no auto update



## Usage

You can first check the available options and parameters by running:
`nextflow run /path/to/nf-wgs_amr/main.nf -profile singu --help`

You always need to select a profile to run the workflow with.

We provide `singu`, a profile using Singularity containers only with Slurm executor.
We also provide `itrop` an example of config with a module environment, where some tools are run from modules installed on a HPC environment.

The `local` provide is a config using Docker and without Slurm.

Feel free to add your own favourite config, in the `conf` folder.

## Test the workflow

We provide a test directory containing 3 illumina tests datasets, of *E. coli*,  that have been downsampled to be lighter. You can run this, from the directory of your choice, as long as you give the path to the nf-wgs_amr directory:

```
nextflow run /path/to/nf-wgs_amr/main.nf -profile singu \
  --illumina '/path/to/nf-wgs_amr/data' --genus 'Escherichia' --species 'coli' \
  --busco_lineage 'enterobacterales_odb10' --amrfinder_organism 'Escherichia' \
  --species_taxid '562' --output './results_test'
```


## Parameters

For running the workflow needs 3 main parameters:
1. the input datasets: 3 possible inputs:
  - directory containing paired-end short reads (Illumina type)
  - directory containing already assembled contigs/scaffolds
  - an index CSV file indicating path to short reads and long reads; for hybrid input requiring Unicycler hybrid assembly.
  The CSV index file should look as below and must include the columns headers:

```
sampleID,read1,read2,ont
124,test_illu_hybrid/124_1.fq,test_illu_hybrid/124_2.fq,test_ont/barcode05_concat.fastq
365,test_illu_hybrid/365_1.fq,test_illu_hybrid/365_2.fq,test_ont/barcode01_concat.fastq
```

2. Three mandatory databases: in the directory .... already in the nextflow.config. To overwrite in the command line if different

3. A TaxID (NCBI Taxonomy ID) to which extract from to get "decontaminated" scaffolds/contigs belonging to the expected bacterial taxon. It can be a TaxID corresponding to an *order*, a *genus* or a *species*, and all the contigs classified by Kraken2 under this specified taxon and lower in the taxonomy (children taxa) will be retrieved as decontaminated.


You can avoid writing all the parameters by providing a config file containing the parameters (e.g. paths to databases, busco lineage...)
here is an example config:

```
// Workflow parameters
params.output = "./results"
params.tmpdir = "./tmpdir"

//db
params.card_db = "/data2/projects/ARCAHE/databases/CARD_2022/card.json"
params.kraken2_db = "/data/projects/banks/kraken2/nt/21-09/nt/"
params.plasmidfinder_db = "database/plasmidfinder/plasmidfinder_db"

params.amrfinder_db = "/databases/amrfinder/latest"
params.bakta_db = "/databases/bakta_db_2208/db/"
params.busco_db_offline = "/databases/busco_downloads"
params.platon_db = "/databases/platon/db"

// Species options
params.amrfinder_organism = "Escherichia"
params.busco_lineage = "enterobacterales_odb10"
params.genus = "Escherichia"
params.species = "coli"
params.species_taxid = "562"

// Nextflow configuration options
workDir = './work'

profiles {
  itrop {
    process.executor = 'slurm'
    includeConfig "conf/itrop.config"
  }
  singu {
    process.executor = 'slurm'
    includeConfig "conf/singu.config"
  }
}

process {
    clusterOptions = '-p highmemplus --nodelist=node5'
    // You can also override existing process cpu or time settings here too.
}
```

If you have such a file, you can run the workflow that way:

```
nextflow run nf-wgs_amr/main.nf -profile singu \
  -c '/path_to_my_params/params_node5_slurm.config' \
  --illumina 'path/to/your/illumina/reads_folder' \
  --output 'results_Ecoli'
```


## Complete help and options

````
********* Workflow for bacterial genome assembly and detection of antimicrobial resistances and plasmids *********

    Usage example:
nextflow run main.nf --illumina short_reads_Ecoli --genus Escherichia --species coli --species_taxid 562 -profile singu -resume
--help                      prints the help section

    Input sequences:
--illumina                  path to the directory containing the illumina read file (fastq) (default: null)
--contigs                   path to the directory containing the already assembled contigs files (fasta) (default: null)
--hybrid_index              For users having both short and long reads:
                            path to the CSV file containing the mapping between sampleID, illuminaR1.fastq.gz, illuminaR2.fastq.gz, ont_read.fastq
                            Must have the header as follow:
                            sampleID,read1,read2,ont

                            Example of CSV index file:
                            sampleID,read1,read2,ont
                            124,test_illu_hybrid/124_1.fq,test_illu_hybrid/124_2.fq,test_ont/barcode05_concat.fastq
                            365,test_illu_hybrid/365_1.fq,test_illu_hybrid/365_2.fq,test_ont/barcode01_concat.fastq

    Output:
--output                    path to the output directory (default: ./results)
--tmpdir                    path to the tmp directory (default: ./tmpdir)

    Species mandatory options:
--genus                     Bacterial genus (Escherichia, Salmonella, Enterobacter, Klebsiella, Staphylococcus)  [default: null]
--species                   bacterial species to assemble (e.g. coli, pneumoniae, cloacae, aureus) [default: null]
--species_taxid             NCBI TaxID of the bacterial species to assemble [default: null]

    Databases path required (script provided for downloading them):
--card_db                   path to the CARD json Database for Antimicrobial Resistance Genes prediction [default: ]
--kraken2_db                path to the local Kraken2 nucleotide database (e.g. MiniKraken, nt, standard) [default: ]
--plasmidfinder_db          path to the CGE PlasmidFinder database [default: ]

    Optional databases paths: if provided, the tool is run:
--amrfinder_db              path to a local AMRFinder Database for Antimicrobial Resistance Genes prediction [default: ] - a database if provided within the container

--bakta_db                  path to the Bakta local database if the user prefers annotating the genomes with Bakta instead of Prokka [default: ]
--busco_db_offline          path to local BUSCO datasets if user wants to run BUSCO offline [default: null]
--platon_db                 path to the Platon local database

   Optional input:
--phred_type                phred score type. Specify if 33 (default and current) or 64 (ex. BGI, older...) [default: 33]
--busco_lineage             to specify according to the bacterial species. e.g. enterobacterales_odb10, bacillales_odb10... check BUSCO [default: null]
                            If not provided, Busco will use prokaryotes database
--amrfinder_organism        To specify for PointMutation detection
                            Can be among these: Acinetobacter_baumannii, Campylobacter,
                            Clostridioides_difficile, Enterococcus_faecalis, Enterococcus_faecium,
                            Escherichia, Klebsiella, Neisseria, Pseudomonas_aeruginosa,
                            Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius,
                            Streptococcus_agalactiae, Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae.
                            The amrfinderplus will be run if not specified, but no point mutations are detected.
                            [default: null]
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
````

## Update

## Uninstall

## Citation

## Author and contributors

Juliette Hayer  (@jhayer)
Jacques Dainat  (@Juke34)
