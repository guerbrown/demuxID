#!/bin/bash
# Master script to run the oak gall wasp and parasite sequence identification pipeline

# Default parameters
CONFIG_FILE="config.ini"
CREATE_DB=false
BLAST_ONLY=false
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# Help function
function show_help {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help             Show this help message"
    echo "  -c, --config FILE      Specify config file (default: config.ini)"
    echo "  --create-db            Create a BLAST database only"
    echo "  --blast-only           Run BLAST automation only (skip database creation)"
    echo "  --script-dir DIR       Specify script directory (default: ./script)"
    echo ""
    echo "This script runs the oak gall wasp and parasite sequence identification pipeline."
    echo "It can create a custom BLAST database and process FASTA files to identify species."
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            show_help
            exit 0
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --create-db)
            CREATE_DB=true
            shift
            ;;
        --blast-only)
            BLAST_ONLY=true
            shift
            ;;
        --script-dir)
            SCRIPT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found!"
    exit 1
fi

# Extract values from config file
EMAIL=$(grep -A10 "^\[Blast\]" "$CONFIG_FILE" | grep "email" | cut -d'=' -f2 | tr -d ' ')
INPUT_DIR=$(grep -A10 "^\[Paths\]" "$CONFIG_FILE" | grep "input_dir" | cut -d'=' -f2 | tr -d ' ')
OUTPUT_DIR=$(grep -A10 "^\[Paths\]" "$CONFIG_FILE" | grep "output_dir" | cut -d'=' -f2 | tr -d ' ')
DB_PATH=$(grep -A10 "^\[Paths\]" "$CONFIG_FILE" | grep "db_path" | cut -d'=' -f2 | tr -d ' ')
CONFIG_SCRIPT_DIR=$(grep -A10 "^\[Paths\]" "$CONFIG_FILE" | grep "script_dir" | cut -d'=' -f2 | tr -d ' ')

# Use script_dir from config if provided
if [ -n "$CONFIG_SCRIPT_DIR" ]; then
    SCRIPT_DIR="$CONFIG_SCRIPT_DIR"
fi

# Check if script directory exists
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Error: Script directory '$SCRIPT_DIR' not found!"
    exit 1
fi

# Check if scripts exist
BLASTDB_SCRIPT="$SCRIPT_DIR/blastdb-creation.py"
BLAST_SCRIPT="$SCRIPT_DIR/blast-automation.py"

if [ ! -f "$BLASTDB_SCRIPT" ]; then
    echo "Error: BLAST database creation script not found at '$BLASTDB_SCRIPT'!"
    exit 1
fi

if [ ! -f "$BLAST_SCRIPT" ]; then
    echo "Error: BLAST automation script not found at '$BLAST_SCRIPT'!"
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Create database if requested
if [ "$CREATE_DB" = true ] && [ "$BLAST_ONLY" = false ]; then
    echo "Creating BLAST database..."
    DB_DIR=$(dirname "$DB_PATH")
    
    if [ ! -d "$DB_DIR" ]; then
        echo "Creating database directory: $DB_DIR"
        mkdir -p "$DB_DIR"
    fi
    
    python3 "$BLASTDB_SCRIPT" --use_config --email "$EMAIL" --output_dir "$DB_DIR"
    
    if [ $? -ne 0 ]; then
        echo "Error creating BLAST database!"
        exit 1
    fi
    
    echo "BLAST database created successfully."
fi

# Run BLAST automation
if [ "$CREATE_DB" = false ] || [ "$BLAST_ONLY" = true ]; then
    echo "Running BLAST automation..."
    
    python3 "$BLAST_SCRIPT" --input_dir "$INPUT_DIR" --output_dir "$OUTPUT_DIR" --email "$EMAIL" --local_blast --db_path "$DB_PATH" --config "$CONFIG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "Error running BLAST automation!"
        exit 1
    fi
    
    echo "BLAST automation completed successfully."
fi

echo "All tasks completed."
