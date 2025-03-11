 ```
                                       _____
                                      /  . .\
                                      \    ---<
                                       \    /
_______________________________________/  /
            D E M U X   I D               /
_______________________________________ /
```
# DemuxID
### Demultiplexing Sorted Cluster Identification
##### Designed to supplement clweinrich/demux

DemuxID was created to perform a custom blast process for quickly identifying reads, or clusters of similar reads, from Oxford Nanopore Sequencing. Specifically, it was designed to work in tandem with clweinrich/demux, which outputs a set of folders documenting the processes of dereplication, filtering, clustering, chimera removal, and sorting clusters. **This tools is meant to identify the sorted clusters from the final output of demux and produce a summary of the relevant blast information for each sample.** 

It does so by creating an blast database (blastdb) from NCBI's GenBank. This precise sequences can be specified inside the primary script, demuxid.py, starting at line 184. The reason for including this as part of the main script is to minimize the additional work that needs to be done by the end user to create a stable sequence identification tool for an entire project. This is meant to be highly compressed such that only one script is needed to run the entire process. The hope is that this reduces the possiblity that individual users who may be less familiar with the techicalities of blast break the current system. 

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/guerbrown/demuxid.git
   cd demuxid
   ```

2. Install either Mamba (recommended) or Conda

3. Set up the Mamba or Conda environment and install dependencies
```
# This process has been automated for you as part of the demux.py script
./demux.py setup-env
```
Notes on installation: This was designed in a Windows Subsystem for Linux running Ubuntu. If you are running this on an HPC or a native linux interface, you might need to install additional packages

## Running
### 1. Configuration
The following script will create config.ini, which is required for any and all downstream processes.
```
./demuxid.py init --config config.ini
```

After creating config.ini you will need to update the parameters inside. Although this should be self-explantory, here is a detailed explanation just in case you are stuck:
- The first section is [Paths], which means the "address" of the files on your computer. The first is "script_dir", which is your demuxid.py script. If you are running your commands within the same working directory as demuxid.py, this can be left blank. However, if you are running the ./demuxid.py from any other directory, you will need to have the directory for demuxid.py listed here.
- Next are the input and output directories (input_dir and output_dir, respectively). The input directory should be the path to the directory with you sorted files from the output of your original demux run. If you are not running demux, this should be the folder with the .fasta sequences you would like to automate a blast process for. The output is where the annotated sequences and the taxonomic summary will be left
- db_path is the location you would like to store you blastdb, or the files you will search against to identify your sequences.
- It is recommended that you store only the files for certain processes in their respective folders with nothing else.

Customizing Output
The config.ini file allows you to customize which information is included in the FASTA headers:
- include_taxonomy: Include taxonomic identification (Y/N)
- include_accession: Include NCBI accession numbers (Y/N)
- include_identity: Include percent identity (Y/N)
- include_evalue: Include E-value information (Y/N)

### 2. Running the complete pipeline
Note: Before you do this, please check to line demuxid.py starting at line 184 that the taxa and genes listed are the ones *you* want to use. This was designed for in-house operation guerbrown. Usage outside of its original purpose requires updating this list.

```
./demuxid.py full --config config.ini
```

### 3. Running individual steps
Note: Whether you are rerunning this *after* making your initial blastdb or debugging, each of the steps can be run indepent of the others

```
# Create the database only
./demuxid.py create-db --config config.ini

# Process sequences only (requires a database)
./demux.py blast --config config.ini
```

## Prerequisites
- Python 3.9 or higher
- Git (for cloning the repository)
- Internet connection (for database creation)
- At least 10GB of free disk space (for the BLAST database)

## Output
DemuxID generates two main outputs:
1. **Annotated FASTA files** - Original sequences with headers updated to include taxonomic information
2. **Taxonomic summary file** - A text file (`taxonomic_summary.txt`) containing detailed information about each cluster, including:
   - Taxonomic identification
   - Sequence count
   - Percent identity
   - E-value
   - Accession numbers

## Troubleshooting
- **XML parsing errors** - These may appear during BLAST processing but typically don't affect results
- **Database creation failures** - Ensure you have a stable internet connection and sufficient disk space
- **Environment setup issues** - Try manually installing dependencies: `mamba install -n demuxid -c bioconda blast`

