#!/usr/bin/env python3
"""
TeloTracker: A comprehensive tool for telomere analysis
Main entry point with modular architecture
"""

import argparse
import sys
import os
from pathlib import Path
import logging
from typing import Optional, Dict, Any

from .version import __version__
from .log import setup_logging
from .utils import validate_input_files, create_output_directory
from .reference import ReferenceModule
from .track import TrackModule


def create_parser() -> argparse.ArgumentParser:
    """Create the main argument parser with subcommands."""
    
    parser = argparse.ArgumentParser(
        prog='telotracker',
        description='TeloTracker: Comprehensive telomere sequence analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create reference and annotations
  telotracker reference --genome genome.fasta --output reference_dir/
  
  # Analyze sequencing data
  telotracker track --reference reference_dir/ --reads reads.fastq --output results/
  
  # Full pipeline
  telotracker reference --genome genome.fasta --output ref/ && \\
  telotracker track --reference ref/ --reads reads.fastq --output results/
        """
    )
    
    parser.add_argument(
        '--version', 
        action='version', 
        version=f'TeloTracker {__version__}'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='count',
        default=0,
        help='Increase verbosity (use -v, -vv, or -vvv)'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress all output except errors'
    )
    
    parser.add_argument(
        '--log-file',
        type=str,
        help='Write logs to specified file'
    )
    
    # Create subparsers for the two main modules
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        metavar='COMMAND'
    )
    
    # Reference module subparser
    reference_parser = subparsers.add_parser(
        'reference',
        help='Create reference genome and annotations for telomere analysis',
        description='Build reference resources including telomere annotations, X-elements, and Y-prime sequences'
    )
    
    reference_parser.add_argument(
        '--genome', '-g',
        required=True,
        type=str,
        help='Input genome FASTA file'
    )
    
    reference_parser.add_argument(
        '--output', '-o',
        required=True,
        type=str,
        help='Output directory for reference files'
    )
    
    reference_parser.add_argument(
        '--telomere-seq',
        type=str,
        default='TTAGGG',
        help='Telomere repeat sequence (default: TTAGGG)'
    )
    
    reference_parser.add_argument(
        '--min-telomere-length',
        type=int,
        default=100,
        help='Minimum telomere tract length to annotate (default: 100)'
    )
    
    reference_parser.add_argument(
        '--threads', '-t',
        type=int,
        default=1,
        help='Number of threads to use (default: 1)'
    )
    
    reference_parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing output directory'
    )
    
    # Track module subparser
    track_parser = subparsers.add_parser(
        'track',
        help='Analyze sequencing data for telomere features',
        description='Analyze reads for telomere length, recombination events, and other features'
    )
    
    track_parser.add_argument(
        '--reference', '-r',
        required=True,
        type=str,
        help='Reference directory created by "telotracker reference"'
    )
    
    track_parser.add_argument(
        '--reads',
        required=True,
        type=str,
        nargs='+',
        help='Input sequencing reads (FASTQ format, single or paired-end)'
    )
    
    track_parser.add_argument(
        '--output', '-o',
        required=True,
        type=str,
        help='Output directory for analysis results'
    )
    
    track_parser.add_argument(
        '--sample-name', '-s',
        type=str,
        help='Sample name for output files (default: inferred from input)'
    )
    
    track_parser.add_argument(
        '--threads', '-t',
        type=int,
        default=1,
        help='Number of threads to use (default: 1)'
    )
    
    track_parser.add_argument(
        '--memory',
        type=str,
        default='8G',
        help='Maximum memory usage (default: 8G)'
    )
    
    track_parser.add_argument(
        '--min-mapq',
        type=int,
        default=20,
        help='Minimum mapping quality for analysis (default: 20)'
    )
    
    track_parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing output directory'
    )
    
    # Analysis-specific options
    track_parser.add_argument(
        '--skip-mapping',
        action='store_true',
        help='Skip read mapping (use existing BAM files)'
    )
    
    track_parser.add_argument(
        '--analyze-circles',
        action='store_true',
        help='Perform extrachromosomal circular DNA analysis'
    )
    
    track_parser.add_argument(
        '--find-recombination',
        action='store_true',
        help='Search for subtelomeric recombination events'
    )
    
    return parser


def validate_arguments(args: argparse.Namespace) -> bool:
    """Validate command-line arguments."""
    
    if args.command == 'reference':
        # Validate reference command arguments
        if not os.path.exists(args.genome):
            logging.error(f"Genome file not found: {args.genome}")
            return False
            
        # Check if output directory exists and handle --force
        if os.path.exists(args.output) and not args.force:
            logging.error(f"Output directory exists: {args.output}. Use --force to overwrite.")
            return False
            
    elif args.command == 'track':
        # Validate track command arguments
        if not os.path.exists(args.reference):
            logging.error(f"Reference directory not found: {args.reference}")
            return False
            
        # Check for required reference files
        required_ref_files = [
            'genome.fasta',
            'telomere_annotations.bed',
            'reference_info.json'
        ]
        
        for req_file in required_ref_files:
            ref_file = os.path.join(args.reference, req_file)
            if not os.path.exists(ref_file):
                logging.error(f"Required reference file missing: {ref_file}")
                return False
        
        # Validate read files
        for read_file in args.reads:
            if not os.path.exists(read_file):
                logging.error(f"Read file not found: {read_file}")
                return False
        
        # Check output directory
        if os.path.exists(args.output) and not args.force:
            logging.error(f"Output directory exists: {args.output}. Use --force to overwrite.")
            return False
    
    return True


def main() -> int:
    """Main entry point for TeloTracker."""
    
    parser = create_parser()
    
    # If no arguments provided, show help
    if len(sys.argv) == 1:
        parser.print_help()
        return 1
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = setup_logging(
        verbose=args.verbose,
        quiet=args.quiet,
        log_file=getattr(args, 'log_file', None)
    )
    
    logging.info(f"Starting TeloTracker {__version__}")
    logging.info(f"Command: {' '.join(sys.argv)}")
    
    # Validate arguments
    if not validate_arguments(args):
        return 1
    
    try:
        # Execute the appropriate module
        if args.command == 'reference':
            logging.info("Running reference module")
            reference_module = ReferenceModule()
            success = reference_module.run(args)
            
        elif args.command == 'track':
            logging.info("Running track module")
            track_module = TrackModule()
            success = track_module.run(args)
            
        else:
            logging.error(f"Unknown command: {args.command}")
            parser.print_help()
            return 1
        
        if success:
            logging.info("TeloTracker completed successfully")
            return 0
        else:
            logging.error("TeloTracker failed")
            return 1
            
    except KeyboardInterrupt:
        logging.warning("Analysis interrupted by user")
        return 130
    except Exception as e:
        logging.error(f"Unexpected error: {str(e)}")
        if args.verbose >= 2:
            logging.exception("Full traceback:")
        return 1


if __name__ == '__main__':
    sys.exit(main())
