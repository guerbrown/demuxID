 ```
                                       _____
                                      /  . .\
                                      \    ---<
                                       \   /
_______________________________________/  /
            D E M U X   I D              /
_______________________________________ /
```
# DemuxID
### Demultiplexing Sorted Cluster Identification
##### Designed to supplement clweinrich/demux

DemuxID was created to perform a custom blast process for quickly identifying reads, or clusters of similar reads, from Oxford Nanopore Sequencing. Specifically, it was designed to work in tandem with clweinrich/demux, which outputs a set of folders documenting the processes of dereplication, filtering, clustering, chimera removal, and sorting clusters. This tools is meant to identify the sorted clusters from the final output of demux. 

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


