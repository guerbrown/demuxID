#!/usr/bin/env python3
"""
DemuxID - Oak Gall Wasp Sequence Identification Pipeline

This script combines database creation and BLAST automation into a single tool.

Commands:
  setup-env    Set up a conda environment with required dependencies
  create-db    Create a BLAST database for oak gall wasps and parasites
  blast        Run BLAST automation on FASTA files
  full         Run the complete pipeline (create DB then BLAST)
  init         Create a default configuration file

Usage:
  demuxid.py setup-env
  demuxid.py create-db --email EMAIL --output_dir DIR [options]
  demuxid.py blast --input_dir DIR --output_dir DIR --db_path PATH [options]
  demuxid.py full --input_dir DIR --output_dir DIR [options]
  demuxid.py (-h | --help)

Options:
  -h --help                Show this help message
  --config FILE            Path to configuration file
  --email EMAIL            Email for NCBI Entrez/BLAST
  --input_dir DIR          Directory with input FASTA files
  --output_dir DIR         Directory for output files
  --db_path PATH           Path to BLAST database
  --db_dir DIR             Directory for database files
  --max_seq NUM            Maximum sequences for DB [default: 50000]
  --batch_size NUM         Batch size for downloads [default: 500]
  --delay NUM              Delay between BLAST requests [default: 10]
  --identity_threshold NUM Minimum percent identity [default: 80.0]
  --evalue_threshold NUM   Maximum E-value [default: 1e-10]
  --taxa TEXT              Custom taxa (format: "Family1:Genus1,Genus2;Family2")
  --gene TEXT              Custom genes (format: "COX1,COI")
  --input_file FILE        Process a single FASTA file instead of directory
"""

import os
import sys
import argparse
import time
import tempfile
import configparser
import re
import subprocess
from importlib.util import find_spec
import shutil

def check_environment():
    """
    Check if the required packages are installed and BLAST+ is available.
    If not, provide instructions for setting up the environment.
    
    Returns:
        bool: True if all requirements are met, False otherwise
    """
    missing_packages = []
    required_packages = ['Bio']
    
    # Check for required Python packages
    for package in required_packages:
        if find_spec(package) is None:
            missing_packages.append(package)
    
    # Check for BLAST+ binaries
    blast_binaries = ['makeblastdb', 'blastn']
    missing_binaries = []
    for binary in blast_binaries:
        if shutil.which(binary) is None:
            missing_binaries.append(binary)
    
    if not missing_packages and not missing_binaries:
        return True
    
    # Print environment setup instructions
    print("DemuxID Environment Setup Required")
    print("==================================")
    
    if missing_packages:
        print("\nMissing Python packages:")
        for package in missing_packages:
            print(f"  - {package}")
        
        print("\nTo install missing packages, you can use:")
        print(f"pip install {' '.join(missing_packages)}")
        print("\nOR set up a dedicated environment:")
        print("\nOption 1: Using conda/mamba (recommended):")
        print("mamba create -n demuxid python=3.9 biopython blast")
        print("conda activate demuxid")
        print("\nOption 2: Using Python venv:")
        print("python3 -m venv demuxid_env")
        print("source demuxid_env/bin/activate")
        print("pip install biopython")
        print("\nOption 3: Automated setup (this will create a conda environment):")
        print("./demuxid.py setup-env")
    
    if missing_binaries:
        print("\nMissing BLAST+ binaries:")
        for binary in missing_binaries:
            print(f"  - {binary}")
        
        print("\nTo install BLAST+:")
        print("Option 1: Using conda/mamba:")
        print("mamba install -c bioconda blast")
        print("\nOption 2: Using system package manager:")
        print("sudo apt-get install ncbi-blast+  # Debian/Ubuntu")
        print("sudo yum install blast  # CentOS/RHEL")
    
    return False

def setup_environment():
    """
    Set up a conda environment with all required dependencies.
    
    Returns:
        bool: True if setup was successful, False otherwise
    """
    print("Setting up DemuxID environment...")
    
    # Prefer mamba over conda
    conda_cmd = None
    if shutil.which('mamba'):
        conda_cmd = 'mamba'
    elif shutil.which('conda'):
        conda_cmd = 'conda'
    
    if conda_cmd is None:
        print("Error: Neither mamba nor conda found. Please install mamba first.")
        print("Visit: https://github.com/mamba-org/mamba")
        return False
    
    # Create a new environment with Python and Biopython
    env_name = 'demuxid'
    print(f"Creating '{env_name}' environment using {conda_cmd}...")
    
    try:
        # First create environment with Python and Biopython
        subprocess.run(
            [conda_cmd, 'create', '-y', '-n', env_name, 'python=3.9', 'biopython'],
            check=True
        )
        
        # Then install BLAST from bioconda
        print(f"Installing BLAST from bioconda...")
        subprocess.run(
            [conda_cmd, 'install', '-y', '-n', env_name, '-c', 'bioconda', 'blast'],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error creating environment: {e}")
        return False
    
    # Get the environment activation script path
    if os.name == 'nt':  # Windows
        activate_script = os.path.join(os.path.dirname(shutil.which(conda_cmd)), 'activate.bat')
        activate_cmd = f'call "{activate_script}" {env_name}'
    else:  # Unix/Linux
        activate_script = os.path.join(os.path.dirname(os.path.dirname(shutil.which(conda_cmd))), 'etc', 'profile.d', 'conda.sh')
        activate_cmd = f'source "{activate_script}" && conda activate {env_name}'
    
    print("\nEnvironment setup complete!")
    print(f"\nTo activate the environment, run:")
    print(f"{activate_cmd}")
    print("\nAfter activation, you can run DemuxID commands:")
    print("./demuxid.py init --config config.ini")
    
    return True

#########################################
# Database Creation Functions
#########################################

def get_oak_gall_wasp_config():
    """
    Return a predefined configuration for oak gall wasps and parasites.
    
    Returns:
        Dictionary with predefined taxa and genes
    """
    # Define all the taxa based on the comprehensive list provided
    taxa_dict = {
        # Main families
        "Cynipidae": [
            "Cynipini", 
            "Synergini", 
            "Ceroptresini", 
            "Aulacideini", 
            "Aylacini", 
            "Diastrophini"
        ],
        "Diplolepididae": [],
        "Figitidae": ["Euceroptres"],
        
        # Chalcidoidea superfamily and its families
        "Chalcidoidea": [],
        "Eurytomidae": [],
        "Torymidae": [],
        "Eulophidae": [],
        "Eupelmidae": [],
        "Megastigmidae": [],
        "Ormyridae": [],
        "Bethyliidae": [],
        "Chalcidae": [],
        "Crabonidae": [],
        "Pteromalidae": []
    }
    
    # Define specific genera that need special focus
    special_genera = [
        "Andricus",
        "Neuroterus", 
        "Cynips",
        "Biorhiza",
        "Dryocosmus",
        "Synergus",
        "Saphonecrus",
        "Euceroptres",
        "Ceroptres",
        "Torymus",
        "Ormyrus",
        "Eurytoma",
        "Sycophila",
        "Eupelmus"
    ]
    
    # Add the special genera to the taxa dictionary
    for genus in special_genera:
        taxa_dict[genus] = []
    
    # Define genes of interest - include all possible naming variations
    gene_list = [
        "COX1", "COI", "COXI", "CO1", 
        "cytochrome oxidase subunit 1", 
        "cytochrome oxidase subunit I",
        "cytochrome c oxidase subunit I",
        "cytochrome c oxidase subunit 1"
    ]
    
    return {"taxa": taxa_dict, "genes": gene_list}

def parse_taxa_input(taxa_input):
    """
    Parse the taxa input string into a structured dictionary.
    
    The input format is expected to be a comma-separated list of family names,
    with optional subfamilies/genera for each family specified as:
    Family1:Genus1,Genus2;Family2:Genus3,Genus4
    
    Args:
        taxa_input: Taxa input string
        
    Returns:
        Dictionary mapping families to lists of subfamilies/genera
    """
    taxa_dict = {}
    
    # Split by semicolon to get family groups
    family_groups = taxa_input.split(";")
    
    for group in family_groups:
        if ":" in group:
            # This family has specific subfamilies/genera
            family, subtaxa = group.split(":", 1)
            family = family.strip()
            taxa_dict[family] = [taxon.strip() for taxon in subtaxa.split(",")]
        else:
            # Just a family name
            family = group.strip()
            if family:
                taxa_dict[family] = []
    
    return taxa_dict

def build_search_term(taxa_dict, gene_list):
    """
    Build an NCBI search term for the specified taxa and genes.
    
    Args:
        taxa_dict: Dictionary of taxonomic groups organized by family
        gene_list: List of genes/regions
        
    Returns:
        NCBI search term string
    """
    # Build the taxa search term
    taxa_terms = []
    
    # Add all family names
    for family in taxa_dict.keys():
        taxa_terms.append(f"{family}[Organism]")
    
    # Add all subfamilies/genera
    for family, subtaxa in taxa_dict.items():
        for taxon in subtaxa:
            if taxon:  # Skip empty entries
                taxa_terms.append(f"{taxon}[Organism]")
    
    # Combine all taxa with OR
    taxa_term = " OR ".join(taxa_terms)
    
    # Create the gene search term
    gene_term = " OR ".join([f"{gene}[Gene]" for gene in gene_list])
    
    # Also search for gene names in the title as sometimes they're not properly annotated
    title_term = " OR ".join([f"{gene}[Title]" for gene in gene_list])
    
    return f"({taxa_term}) AND ({gene_term} OR {title_term})"

def download_sequences(search_term, output_file, email, max_seq=50000, batch_size=500):
    """
    Download sequences from NCBI using the specified search term.
    
    Args:
        search_term: NCBI search term
        output_file: Path to save the sequences
        email: Email for NCBI Entrez
        max_seq: Maximum number of sequences to download
        batch_size: Number of sequences to download in each batch
    """
    # Now that we know we need Bio, import it
    from Bio import Entrez
    
    Entrez.email = email
    
    # First, search for the IDs
    print(f"Searching NCBI with term: {search_term}")
    try:
        search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_seq)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        total_ids = len(id_list)
        
        if total_ids == 0:
            print("No sequences found with the specified search term.")
            return False
        
        print(f"Found {total_ids} sequences. Downloading...")
    except Exception as e:
        print(f"Error during Entrez search: {e}")
        return False
    
    # Download in batches to avoid large requests
    try:
        with open(output_file, 'w') as out_handle:
            for start in range(0, total_ids, batch_size):
                end = min(start + batch_size, total_ids)
                batch_ids = id_list[start:end]
                
                print(f"Downloading batch {start//batch_size + 1} ({end}/{total_ids} sequences)")
                
                # Fetch the sequences
                try:
                    fetch_handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="fasta", retmode="text")
                    out_handle.write(fetch_handle.read())
                    fetch_handle.close()
                    
                    # Add a small delay to avoid overloading NCBI servers
                    time.sleep(1)
                except Exception as e:
                    print(f"Error downloading batch: {e}")
                    continue
        
        print(f"Downloaded sequences saved to: {output_file}")
        return True
    except Exception as e:
        print(f"Error writing sequences to file: {e}")
        return False

def create_database(args):
    """Create a BLAST database for oak gall wasps and parasites."""
    # Use predefined configuration if no custom taxa specified
    if args.taxa is None:
        config = get_oak_gall_wasp_config()
        taxa_dict = config["taxa"]
        gene_list = config["genes"]
        print("Using predefined oak gall wasp and parasite configuration.")
    else:
        # Parse taxa input
        taxa_dict = parse_taxa_input(args.taxa)
        
        # Parse gene list
        gene_list = [gene.strip() for gene in args.gene.split(",")] if args.gene else ["COX1", "COI"]
    
    # Build search term
    search_term = build_search_term(taxa_dict, gene_list)
    
    # Print the search term for debugging
    print(f"\nUsing NCBI search term: {search_term}\n")
    
    # Determine where to store the database
    if hasattr(args, 'db_dir') and args.db_dir:
        # Use explicitly provided db_dir argument
        db_dir = args.db_dir
    elif hasattr(args, 'db_path') and args.db_path:
        # Use db_path from config directly
        db_dir = args.db_path
    else:
        # Fall back to output_dir
        db_dir = args.output_dir
    
    # Create directory if it doesn't exist
    print(f"Using database directory: {db_dir}")
    os.makedirs(db_dir, exist_ok=True)
    
    # Get database name
    db_name = getattr(args, 'db_name', "Oak_Gall_Wasps_Parasites_COX1_db")
    
    # Create paths for files
    fasta_file = os.path.join(db_dir, f"{db_name}.fasta")
    
    # Download sequences
    if download_sequences(
        search_term=search_term,
        output_file=fasta_file,
        email=args.email,
        max_seq=args.max_seq,
        batch_size=args.batch_size
    ):
        # Create BLAST database
        success = create_blast_database(
            input_file=fasta_file,
            output_dir=db_dir,
            db_name=db_name
        )
        
        if success:
            full_db_path = os.path.join(db_dir, db_name)
            print(f"\nDatabase created successfully at: {full_db_path}")
            print(f"Use this path with the 'blast' command: --db_path {full_db_path}")
            return 0
    
    print("\nError: Database creation failed.")
    return 1

#########################################
# BLAST Automation Functions
#########################################

def get_header_info(config, blast_result, cluster_num):
    """
    Get the header information based on config settings.
    
    Args:
        config: Configuration dictionary with header preferences
        blast_result: Dictionary with BLAST results
        cluster_num: Cluster number
        
    Returns:
        String with formatted header information
    """
    header_parts = [f"Cluster-{cluster_num}"]
    
    # Add taxonomic ID if requested
    if config.get('include_taxonomy', 'Y').upper() == 'Y':
        taxonomy = clean_taxonomy_string(blast_result['taxonomy'])
        header_parts.append(taxonomy)
    
    # Add accession number if requested
    if config.get('include_accession', 'N').upper() == 'Y' and 'accession' in blast_result:
        header_parts.append(blast_result['accession'])
    
    # Add percent identity if requested
    if config.get('include_identity', 'N').upper() == 'Y':
        header_parts.append(f"id={blast_result['percent_identity']:.2f}%")
    
    # Add E-value if requested
    if config.get('include_evalue', 'N').upper() == 'Y':
        header_parts.append(f"e={blast_result['e_value']:.2e}")
    
    return "_".join(header_parts)

def blast_sequence(sequence, email, local_blast=False, db_path=None, 
                   identity_threshold=80.0, evalue_threshold=1e-10):
    """
    Run BLAST search for a sequence.
    
    Args:
        sequence: The sequence to BLAST
        email: Email for NCBI BLAST
        local_blast: Whether to use local BLAST
        db_path: Path to local BLAST database
        identity_threshold: Minimum percent identity to accept a hit
        evalue_threshold: Maximum E-value to accept a hit
        
    Returns:
        The best hit taxonomic ID
    """
    # Import Biopython modules
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio.Blast.Applications import NcbiblastnCommandline
    
    if local_blast and db_path:
        # Use local BLAST
        try:
            # Create a temporary file for the sequence
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp:
                temp_filename = temp.name
                SeqIO.write(sequence, temp, "fasta")
            
            # Create a temporary file for the BLAST output
            output_file = tempfile.NamedTemporaryFile(suffix='.xml', delete=False).name
            
            # Run BLAST command
            blastn_cline = NcbiblastnCommandline(
                query=temp_filename,
                db=db_path,
                evalue=evalue_threshold,
                outfmt=5,  # XML output
                out=output_file,
                max_target_seqs=5
            )
            stdout, stderr = blastn_cline()
            
            # Parse the results
            with open(output_file) as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                blast_record = next(blast_records)
                
                # Get the best hit
                if blast_record.alignments:
                    best_alignment = blast_record.alignments[0]
                    best_hsp = best_alignment.hsps[0]
                    
                    # Extract taxonomic ID and accession from the alignment title
                    title_parts = best_alignment.title.split(' ', 1)
                    accession = title_parts[0] if len(title_parts) > 0 else "Unknown"
                    taxonomy = title_parts[1] if len(title_parts) > 1 else "Unknown"
                    
                    percent_identity = (best_hsp.identities / best_hsp.align_length) * 100
                    
                    # Check if it meets the threshold
                    if percent_identity >= identity_threshold and best_hsp.expect <= evalue_threshold:
                        return {
                            'taxonomy': taxonomy,
                            'accession': accession,
                            'percent_identity': percent_identity,
                            'e_value': best_hsp.expect
                        }
                    else:
                        return {
                            'taxonomy': 'Low-identity',
                            'accession': accession, 
                            'percent_identity': percent_identity, 
                            'e_value': best_hsp.expect
                        }
                else:
                    return {
                        'taxonomy': 'No-hit', 
                        'accession': 'None',
                        'percent_identity': 0, 
                        'e_value': 0
                    }
        except Exception as e:
            print(f"Local BLAST error: {e}")
            return {
                'taxonomy': f'BLAST-error: {str(e)}', 
                'accession': 'Error',
                'percent_identity': 0, 
                'e_value': 0
            }
        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_filename)
                os.unlink(output_file)
            except:
                pass
    else:
        # Use NCBI web BLAST
        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence.seq, entrez_query="", 
                                          expect=evalue_threshold, hitlist_size=5, email=email)
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
            # Get the best hit
            if blast_record.alignments:
                best_alignment = blast_record.alignments[0]
                best_hsp = best_alignment.hsps[0]
                
                # Extract taxonomic ID and accession from the alignment title
                title_parts = best_alignment.title.split('|')
                accession = title_parts[1] if len(title_parts) > 1 else "Unknown"
                taxonomy = title_parts[-1].strip() if len(title_parts) > 0 else "Unknown"
                
                percent_identity = (best_hsp.identities / best_hsp.align_length) * 100
                
                # Check if it meets the threshold
                if percent_identity >= identity_threshold and best_hsp.expect <= evalue_threshold:
                    return {
                        'taxonomy': taxonomy,
                        'accession': accession,
                        'percent_identity': percent_identity,
                        'e_value': best_hsp.expect
                    }
                else:
                    return {
                        'taxonomy': 'Low-identity',
                        'accession': accession, 
                        'percent_identity': percent_identity, 
                        'e_value': best_hsp.expect
                    }
            else:
                return {
                    'taxonomy': 'No-hit', 
                    'accession': 'None',
                    'percent_identity': 0, 
                    'e_value': 0
                }
        except Exception as e:
            print(f"NCBI BLAST error: {e}")
            return {
                'taxonomy': f'BLAST-error: {str(e)}', 
                'accession': 'Error',
                'percent_identity': 0, 
                'e_value': 0
            }

def clean_taxonomy_string(taxonomy):
    """
    Clean up the taxonomy string for use in a FASTA header.
    
    Args:
        taxonomy: Raw taxonomy string from BLAST
        
    Returns:
        Cleaned taxonomy string
    """
    # Replace spaces with hyphens
    cleaned = taxonomy.replace(' ', '-')
    
    # Remove any special characters that might cause issues in file names or headers
    cleaned = ''.join(c for c in cleaned if c.isalnum() or c in '-._()')
    
    return cleaned

def parse_fasta_header(header):
    """
    Parse a FASTA header to extract cluster size and other information.
    
    Args:
        header: FASTA header string
        
    Returns:
        Dictionary with parsed information
    """
    result = {'size': 1}  # Default size
    
    # Look for size information in the header
    size_match = re.search(r'size=(\d+)', header)
    if size_match:
        result['size'] = int(size_match.group(1))
    
    # Look for seqs information in the header
    seqs_match = re.search(r'seqs=(\d+)', header)
    if seqs_match:
        result['seqs'] = int(seqs_match.group(1))
    
    return result

def process_fasta_file(file_path, output_dir, email, local_blast=False, db_path=None, 
                       delay=10, identity_threshold=80.0, evalue_threshold=1e-10, header_config=None):
    """
    Process a single FASTA file.
    
    Args:
        file_path: Path to the FASTA file
        output_dir: Directory for output files
        email: Email for NCBI BLAST
        local_blast: Whether to use local BLAST
        db_path: Path to local BLAST database
        delay: Delay between BLAST requests in seconds
        identity_threshold: Minimum percent identity to accept a hit
        evalue_threshold: Maximum E-value to accept a hit
        header_config: Configuration for header information
    
    Returns:
        Dictionary with summary information
    """
    # Import SeqIO
    from Bio import SeqIO
    
    sample_id = os.path.splitext(os.path.basename(file_path))[0]
    output_file = os.path.join(output_dir, f"{sample_id}_blasted.fasta")
    
    try:
        sequences = list(SeqIO.parse(file_path, "fasta"))
        if not sequences:
            print(f"Warning: No sequences found in {file_path}")
            return {'sample_id': sample_id, 'clusters': []}
    except Exception as e:
        print(f"Error parsing file {file_path}: {e}")
        return {'sample_id': sample_id, 'clusters': []}
    
    updated_sequences = []
    clusters_info = []
    
    for i, seq in enumerate(sequences, 1):
        print(f"Processing {sample_id}, sequence {i} of {len(sequences)}")
        
        # Parse the original header for additional information
        header_info = parse_fasta_header(seq.description)
        
        # Run BLAST search
        blast_result = blast_sequence(
            seq, email, local_blast, db_path, 
            identity_threshold, evalue_threshold
        )
        
        # Use header configuration if provided
        if header_config:
            header_suffix = get_header_info(header_config, blast_result, i)
            new_id = f"{sample_id}_{header_suffix}"
        else:
            # Clean up taxonomic ID
            taxonomy = clean_taxonomy_string(blast_result['taxonomy'])
            new_id = f"{sample_id}_Cluster-{i}_{taxonomy}"
        
        # Update the sequence ID
        seq.id = new_id
        seq.description = ""  # Clear the description
        
        updated_sequences.append(seq)
        
        # Store cluster information for summary
        cluster_info = {
            'cluster_num': i,
            'size': header_info.get('size', 1),
            'seqs': header_info.get('seqs', 1),
            'taxonomy': blast_result['taxonomy'],
            'percent_identity': blast_result['percent_identity'],
            'e_value': blast_result['e_value'],
            'accession': blast_result.get('accession', 'Unknown')
        }
        clusters_info.append(cluster_info)
        
        # Add delay to avoid overloading NCBI server if using web BLAST
        if not local_blast and i < len(sequences):
            print(f"Waiting {delay} seconds before next BLAST request...")
            time.sleep(delay)
    
    # Write updated sequences to output file
    try:
        SeqIO.write(updated_sequences, output_file, "fasta")
        print(f"Processed file saved to: {output_file}")
    except Exception as e:
        print(f"Error writing to {output_file}: {e}")
    
    # Return summary information
    return {
        'sample_id': sample_id,
        'clusters': clusters_info
    }

def create_summary_file(output_dir, all_samples_summary):
    """
    Create a summary file with taxonomic information from all samples.
    
    Args:
        output_dir: Directory for output files
        all_samples_summary: List of dictionaries with sample information
    """
    summary_file = os.path.join(output_dir, "taxonomic_summary.txt")
    
    with open(summary_file, 'w') as f:
        for sample in all_samples_summary:
            sample_id = sample['sample_id']
            clusters = sample['clusters']
            
            f.write(f"{sample_id}\n")
            f.write(f"{len(clusters)} Clusters\n")
            
            for cluster in clusters:
                f.write(f"Cluster {cluster['cluster_num']}:\n")
                f.write(f"* {cluster.get('seqs', cluster.get('size', 1))} Seqs\n")
                f.write(f"* {cluster['taxonomy']}\n")
                f.write(f"* Percent Identity: {cluster['percent_identity']:.2f}%\n")
                if cluster['e_value'] != 0:
                    f.write(f"* E-value: {cluster['e_value']:.2e}\n")
                else:
                    f.write(f"* E-value: 0\n")
                f.write(f"* Accession: {cluster['accession']}\n")
                f.write("\n")
            
            f.write("\n")
    
    print(f"Summary file created at: {summary_file}")

def get_default_config():
    """
    Get default configuration for header information.
    
    Returns:
        Dictionary with default configuration
    """

def get_default_config():
    """
    Get default configuration for header information.
    
    Returns:
        Dictionary with default configuration
    """
    return {
        'include_taxonomy': 'Y',
        'include_accession': 'N',
        'include_identity': 'N',
        'include_evalue': 'N'
    }

def load_config_file(config_file):
    """
    Load configuration from a file.
    
    Args:
        config_file: Path to configuration file
        
    Returns:
        Dictionary with configuration
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    
    result = get_default_config()
    
    if 'Header' in config:
        for key in result:
            if key in config['Header']:
                result[key] = config['Header'][key]
    
    return result

def run_blast(args):
    """Run BLAST automation on FASTA files."""
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load header configuration
    header_config = get_default_config()
    if args.config:
        try:
            header_config = load_config_file(args.config)
            print("Loaded header configuration from config file.")
        except Exception as e:
            print(f"Error loading header configuration: {e}")
    
    # Process either a single file or a directory
    all_samples_summary = []
    
    if args.input_file:
        # Process a single file
        print(f"\nProcessing file: {args.input_file}")
        
        sample_summary = process_fasta_file(
            file_path=args.input_file,
            output_dir=args.output_dir,
            email=args.email,
            local_blast=True if args.db_path else False,
            db_path=args.db_path,
            delay=args.delay,
            identity_threshold=args.identity_threshold,
            evalue_threshold=args.evalue_threshold,
            header_config=header_config
        )
        
        all_samples_summary.append(sample_summary)
    else:
        # Process each FASTA file in the input directory
        fasta_files = [f for f in os.listdir(args.input_dir) if f.endswith(('.fasta', '.fa'))]
        total_files = len(fasta_files)
        
        if total_files == 0:
            print(f"No FASTA files found in {args.input_dir}")
            return 1
        
        print(f"Found {total_files} FASTA files to process")
        
        for i, file_name in enumerate(fasta_files, 1):
            file_path = os.path.join(args.input_dir, file_name)
            print(f"\nProcessing file {i} of {total_files}: {file_name}")
            
            sample_summary = process_fasta_file(
                file_path=file_path,
                output_dir=args.output_dir,
                email=args.email,
                local_blast=True if args.db_path else False,
                db_path=args.db_path,
                delay=args.delay,
                identity_threshold=args.identity_threshold,
                evalue_threshold=args.evalue_threshold,
                header_config=header_config
            )
            
            all_samples_summary.append(sample_summary)
    
    # Create summary file
    create_summary_file(args.output_dir, all_samples_summary)
    
    print("\nAll files processed successfully!")
    return 0

#########################################
# Command Line Interface
#########################################

def create_default_config_file(path):
    """
    Create a default configuration file at the specified path.
    
    Args:
        path: Path where to create the config file
    """
    config = configparser.ConfigParser()
    
    config['Paths'] = {
        'script_dir': '',
        'input_dir': '/path/to/fasta_files',
        'output_dir': '/path/to/results',
        'db_path': '/path/to/database',
        'db_name': 'Oak_Gall_Wasps_Parasites_COX1_db'  # Add this line
    }
    
    config['Blast'] = {
        'email': 'your.email@example.com',
        'delay': '10',
        'identity_threshold': '80.0',
        'evalue_threshold': '1e-10'
    }
    
    config['Header'] = {
        'include_taxonomy': 'Y',
        'include_accession': 'Y',
        'include_identity': 'Y',
        'include_evalue': 'N'
    }
    
    with open(path, 'w') as f:
        config.write(f)
    
    print(f"Created default configuration file at: {path}")

def parse_config_file(config_file):
    """
    Parse a configuration file and return a dictionary with settings.
    
    Args:
        config_file: Path to the configuration file
        
    Returns:
        Dictionary with settings
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    
    settings = {}
    
    # Paths section
    if 'Paths' in config:
        for key in ['script_dir', 'input_dir', 'output_dir', 'db_path']:
            if key in config['Paths']:
                settings[key] = config['Paths'][key]
    
    # Blast section
    if 'Blast' in config:
        if 'email' in config['Blast']:
            settings['email'] = config['Blast']['email']
        if 'delay' in config['Blast']:
            settings['delay'] = config.getint('Blast', 'delay')
        if 'identity_threshold' in config['Blast']:
            settings['identity_threshold'] = config.getfloat('Blast', 'identity_threshold')
        if 'evalue_threshold' in config['Blast']:
            settings['evalue_threshold'] = config.getfloat('Blast', 'evalue_threshold')
    
    # Header section
    if 'Header' in config:
        settings['header_config'] = {}
        for key in ['include_taxonomy', 'include_accession', 'include_identity', 'include_evalue']:
            if key in config['Header']:
                settings['header_config'][key] = config['Header'][key]
    
    return settings

def parse_cli_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='DemuxID - Oak Gall Wasp Sequence Identification Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Set up the required environment
  demuxid.py setup-env
  
  # Create a default config file
  demuxid.py init --config config.ini
  
  # Create a BLAST database
  demuxid.py create-db --email your.email@example.com --output_dir ./blastdb
  
  # Run BLAST automation
  demuxid.py blast --input_dir ./fasta_files --output_dir ./results --db_path ./blastdb/Oak_Gall_Wasps_Parasites_COX1_db
  
  # Run the complete pipeline
  demuxid.py full --config config.ini
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Setup environment command
    setup_parser = subparsers.add_parser('setup-env', help='Set up a conda environment with required dependencies')
    
    # Initialize with default config
    init_parser = subparsers.add_parser('init', help='Create a default configuration file')
    init_parser.add_argument('--config', required=True, help='Path to create configuration file')
    
    # Create database command
    db_parser = subparsers.add_parser('create-db', help='Create a BLAST database')
    db_parser.add_argument('--config', help='Path to configuration file')
    db_parser.add_argument('--email', help='Email for NCBI Entrez')
    db_parser.add_argument('--output_dir', help='Directory for output files')
    db_parser.add_argument('--db_dir', help='Directory for database files (overrides output_dir)')
    db_parser.add_argument('--max_seq', type=int, default=50000, help='Maximum number of sequences to download')
    db_parser.add_argument('--batch_size', type=int, default=500, help='Number of sequences to download in each batch')
    db_parser.add_argument('--taxa', help='Custom taxa (format: "Family1:Genus1,Genus2;Family2")')
    db_parser.add_argument('--gene', help='Custom genes (format: "COX1,COI")')
    
    # BLAST command
    blast_parser = subparsers.add_parser('blast', help='Run BLAST automation')
    blast_parser.add_argument('--config', help='Path to configuration file')
    blast_parser.add_argument('--input_dir', help='Directory containing FASTA files')
    blast_parser.add_argument('--output_dir', help='Directory for output files')
    blast_parser.add_argument('--input_file', help='Single FASTA file to process')
    blast_parser.add_argument('--email', help='Email for NCBI BLAST')
    blast_parser.add_argument('--db_path', help='Path to local BLAST database')
    blast_parser.add_argument('--delay', type=int, default=10, help='Delay between BLAST requests in seconds')
    blast_parser.add_argument('--identity_threshold', type=float, default=80.0, help='Minimum percent identity')
    blast_parser.add_argument('--evalue_threshold', type=float, default=1e-10, help='Maximum E-value')
    
    # Full pipeline command
    full_parser = subparsers.add_parser('full', help='Run the complete pipeline')
    full_parser.add_argument('--config', help='Path to configuration file')
    full_parser.add_argument('--input_dir', help='Directory containing FASTA files')
    full_parser.add_argument('--output_dir', help='Directory for output files')
    full_parser.add_argument('--email', help='Email for NCBI Entrez/BLAST')
    full_parser.add_argument('--db_dir', help='Directory for database files')
    
    return parser.parse_args()

def merge_config_args(args):
    """Merge command line arguments with config file settings."""
    if not args.config:
        return args
    
    # Parse config file
    try:
        config = parse_config_file(args.config)
    except Exception as e:
        print(f"Error parsing config file: {e}")
        return args
    
    # Create a new args namespace
    merged_args = argparse.Namespace()
    
    # Copy all attributes from args
    for attr in vars(args):
        setattr(merged_args, attr, getattr(args, attr))
    
    # Override with config file values if not provided on command line
    for key, value in config.items():
        if key == 'header_config':
            # Special case for header config
            merged_args.header_config = value
        elif not hasattr(merged_args, key) or getattr(merged_args, key) is None:
            setattr(merged_args, key, value)
    
    return merged_args

def main():
    """Main entry point for the script."""
    # Check for setup-env command before import checks
    if len(sys.argv) > 1 and sys.argv[1] == 'setup-env':
        if setup_environment():
            return 0
        else:
            return 1
    
    # Check environment for other commands
    if len(sys.argv) > 1 and sys.argv[1] != '--help' and sys.argv[1] != '-h':
        if not check_environment():
            return 1
    
    # Parse command line arguments
    args = parse_cli_args()
    
    if not args.command:
        print("Error: No command specified. Use --help to see available commands.")
        return 1
    
    # Initialize with default config
    if args.command == 'init':
        create_default_config_file(args.config)
        return 0
    
    # Merge command line arguments with config file
    if args.config:
        args = merge_config_args(args)
    
    # Run the requested command
    if args.command == 'create-db':
        # Validate required arguments
        if not args.email:
            print("Error: --email is required for create-db command.")
            return 1
        if not args.output_dir and not args.db_dir:
            print("Error: --output_dir or --db_dir is required for create-db command.")
            return 1
        
        return create_database(args)
    
    elif args.command == 'blast':
        # Validate required arguments
        if not args.email:
            print("Error: --email is required for blast command.")
            return 1
        if not args.output_dir:
            print("Error: --output_dir is required for blast command.")
            return 1
        if not args.input_dir and not args.input_file:
            print("Error: --input_dir or --input_file is required for blast command.")
            return 1
        
        return run_blast(args)
    
    elif args.command == 'full':
        # Validate required arguments
        if not args.email:
            print("Error: --email is required for full command.")
            return 1
        if not args.output_dir:
            print("Error: --output_dir is required for full command.")
            return 1
        if not args.input_dir:
            print("Error: --input_dir is required for full command.")
            return 1
        
        # Set up database directory
        db_dir = args.db_dir if args.db_dir else os.path.join(args.output_dir, 'blastdb')
        os.makedirs(db_dir, exist_ok=True)
        
        # Create database
        db_args = argparse.Namespace(
            command='create-db',
            email=args.email,
            output_dir=None,
            db_dir=db_dir,
            max_seq=50000,
            batch_size=500,
            taxa=None,
            gene=None,
            config=args.config
        )
        
        print("\n=== Creating BLAST Database ===\n")
        result = create_database(db_args)
        if result != 0:
            print("Error during database creation. Aborting pipeline.")
            return result
        
        # Set the database path for BLAST
        db_path = os.path.join(db_dir, "Oak_Gall_Wasps_Parasites_COX1_db")
        
        # Run BLAST automation
        blast_args = argparse.Namespace(
            command='blast',
            email=args.email,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            input_file=None,
            db_path=db_path,
            delay=10,
            identity_threshold=80.0,
            evalue_threshold=1e-10,
            config=args.config,
            header_config=getattr(args, 'header_config', None)
        )
        
        print("\n=== Running BLAST Automation ===\n")
        return run_blast(blast_args)
    
    else:
        print(f"Error: Unknown command '{args.command}'")
        return 1

if __name__ == "__main__":
    sys.exit(main())