#!/usr/bin/env python3
"""
Oak Gall Wasp Sequence BLAST Automation

This script automates the process of:
1. Reading sequences from FASTA files in a directory
2. Running BLAST searches against a local BLAST database
3. Updating the FASTA headers with taxonomic information
4. Saving the updated sequences to new files
5. Creating a summary report of all taxonomic information

Requirements:
- Biopython
- NCBI BLAST+ command line tools

Usage:
python blast-automation.py --input_dir /path/to/fasta_files --output_dir /path/to/output 
                          --email your.email@example.com --local_blast 
                          --db_path /path/to/database/Oak_Gall_Wasps_Parasites_COX1_db
                          --config /path/to/config.ini
"""

import os
import argparse
import time
import tempfile
import configparser
import re
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Automate BLAST searches and update FASTA headers.')
    parser.add_argument('--input_dir', help='Directory containing FASTA files')
    parser.add_argument('--output_dir', help='Directory for output files')
    parser.add_argument('--email', help='Email for NCBI BLAST (required by NCBI)')
    parser.add_argument('--local_blast', action='store_true', help='Use local BLAST instead of NCBI web service')
    parser.add_argument('--db_path', help='Path to local BLAST database (required if using local_blast)')
    parser.add_argument('--delay', type=int, default=10, help='Delay between BLAST requests in seconds (default: 10)')
    parser.add_argument('--identity_threshold', type=float, default=80.0, 
                        help='Minimum percent identity to accept a hit (default: 80.0)')
    parser.add_argument('--evalue_threshold', type=float, default=1e-10, 
                        help='Maximum E-value to accept a hit (default: 1e-10)')
    parser.add_argument('--config', help='Path to configuration file')
    parser.add_argument('--input_file', help='Single FASTA file to process (instead of directory)')
    
    args = parser.parse_args()
    
    # Load config file if specified
    if args.config:
        config = configparser.ConfigParser()
        try:
            config.read(args.config)
            # Override command line arguments with config file values if not provided
            if 'Blast' in config:
                if not args.email and 'email' in config['Blast']:
                    args.email = config['Blast']['email']
                if 'delay' in config['Blast']:
                    args.delay = config.getint('Blast', 'delay')
                if 'identity_threshold' in config['Blast']:
                    args.identity_threshold = config.getfloat('Blast', 'identity_threshold')
                if 'evalue_threshold' in config['Blast']:
                    args.evalue_threshold = config.getfloat('Blast', 'evalue_threshold')
            if 'Paths' in config:
                if not args.input_dir and 'input_dir' in config['Paths']:
                    args.input_dir = config['Paths']['input_dir']
                if not args.output_dir and 'output_dir' in config['Paths']:
                    args.output_dir = config['Paths']['output_dir']
                if not args.db_path and 'db_path' in config['Paths']:
                    args.db_path = config['Paths']['db_path']
                    args.local_blast = True
        except Exception as e:
            print(f"Error reading config file: {e}")
    
    # Validate that we have either input_dir or input_file
    if not args.input_dir and not args.input_file:
        parser.error("Either --input_dir or --input_file must be specified")
    
    # Validate that we have output_dir
    if not args.output_dir:
        parser.error("--output_dir must be specified")
        
    # Validate that we have an email
    if not args.email:
        parser.error("--email must be specified (either via command line or config file)")
    
    # Validate local_blast arguments
    if args.local_blast and not args.db_path:
        parser.error("--db_path is required when using --local_blast")
    
    return args

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

def main():
    """Main function to process all FASTA files in a directory."""
    args = parse_arguments()
    
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
            local_blast=args.local_blast,
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
            return
        
        print(f"Found {total_files} FASTA files to process")
        
        for i, file_name in enumerate(fasta_files, 1):
            file_path = os.path.join(args.input_dir, file_name)
            print(f"\nProcessing file {i} of {total_files}: {file_name}")
            
            sample_summary = process_fasta_file(
                file_path=file_path,
                output_dir=args.output_dir,
                email=args.email,
                local_blast=args.local_blast,
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

if __name__ == "__main__":
    main()
