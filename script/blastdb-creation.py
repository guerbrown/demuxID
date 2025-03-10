#!/usr/bin/env python3
"""
Oak Gall Wasp and Parasite BLAST Database Creator

This script automates the process of:
1. Downloading sequences from NCBI for oak gall wasps and their parasites
2. Creating a local BLAST database from these sequences

Requirements:
- Biopython
- NCBI Entrez Direct utilities (for downloading)
- NCBI BLAST+ command line tools (for creating the database)

Usage:
python blastdb-creation.py --use_config --email your.email@example.com --output_dir /path/to/database
"""

import os
import argparse
import time
import subprocess
from Bio import Entrez

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Download sequences and create a custom BLAST database.')
    
    # Create a mutually exclusive group for the two different modes
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--use_config', action='store_true', 
                          help='Use predefined oak gall wasp and parasite configuration')
    mode_group.add_argument('--taxa', 
                          help='Semicolon-separated list of taxonomic groups with optional subfamilies/genera '
                               '(e.g., "Cynipidae:Cynipini,Synergini;Figitidae")')
    
    # Make gene only required if taxa is specified
    parser.add_argument('--gene', help='Comma-separated list of genes/regions (e.g., "COX1,COI")')
    
    # Always required arguments
    parser.add_argument('--email', required=True, help='Email for NCBI Entrez')
    parser.add_argument('--output_dir', required=True, help='Directory for output files and database')
    
    # Optional arguments
    parser.add_argument('--max_seq', type=int, default=50000, help='Maximum number of sequences to download')
    parser.add_argument('--batch_size', type=int, default=500, help='Number of sequences to download in each batch')
    
    args = parser.parse_args()
    
    # Validate that gene is provided if taxa is specified
    if args.taxa and not args.gene:
        parser.error("--gene is required when --taxa is specified")
    
    return args

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

def create_blast_database(input_file, output_dir, db_name):
    """
    Create a BLAST database from the downloaded sequences.
    
    Args:
        input_file: Path to the FASTA file with sequences
        output_dir: Directory for the database files
        db_name: Name of the database
    """
    db_path = os.path.join(output_dir, db_name)
    
    print(f"Creating BLAST database: {db_name}")
    
    try:
        # Run makeblastdb command
        cmd = [
            "makeblastdb",
            "-in", input_file,
            "-dbtype", "nucl",
            "-out", db_path,
            "-title", db_name,
            "-parse_seqids"
        ]
        
        subprocess.run(cmd, check=True)
        print(f"BLAST database created at: {db_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST database: {e}")
        return False
    except FileNotFoundError:
        print("Error: makeblastdb command not found. Make sure BLAST+ is installed.")
        return False

def main():
    """Main function to download sequences and create a BLAST database."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Use predefined configuration if requested
    if args.use_config:
        config = get_oak_gall_wasp_config()
        taxa_dict = config["taxa"]
        gene_list = config["genes"]
        print("Using predefined oak gall wasp and parasite configuration.")
    else:
        # Parse taxa input
        taxa_dict = parse_taxa_input(args.taxa)
        
        # Parse gene list
        gene_list = [gene.strip() for gene in args.gene.split(",")]
    
    # Build search term
    search_term = build_search_term(taxa_dict, gene_list)
    
    # Print the search term for debugging
    print(f"\nUsing NCBI search term: {search_term}\n")
    
    # Create a descriptive name for files
    db_name = "Oak_Gall_Wasps_Parasites_COX1_db"
    fasta_file = os.path.join(args.output_dir, f"{db_name}.fasta")
    
    # Download sequences
    if download_sequences(
        search_term=search_term,
        output_file=fasta_file,
        email=args.email,
        max_seq=args.max_seq,
        batch_size=args.batch_size
    ):
        # Create BLAST database
        create_blast_database(
            input_file=fasta_file,
            output_dir=args.output_dir,
            db_name=db_name
        )

if __name__ == "__main__":
    main()
