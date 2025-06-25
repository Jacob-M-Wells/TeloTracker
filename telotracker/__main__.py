#!/usr/bin/env python3
"""
This is the main module for TeloTracker.

Copyright 2024-2025 Jacob M. Wells jacobwells1203@gmail.com
https://github.com/Jacob-M-Wells/TeloTracker

This file is part of TeloTracker. TeloTracker is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. TeloTracker is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with TeloTracker.
If not, see <http://www.gnu.org/licenses/>.
"""
#!/usr/bin/env python3

import argparse
import pathlib
import sys

from .basecall import basecall
from .reference import reference_creation_pipeline
from .track import track_telomeres_pipeline

from .log import bold, get_ascii_art
from .version import __version__

def main():
    """
    Entry point for the command line interface.
    
    Subcommands:
    - basecall: Run basecalling for TeloTracker (without trimming)
    - reference: Generate a telomere reference from input reads
    - track: Track telomeres using an existing reference
    """
    args = parse_main_args(sys.argv[1:])
    validate_args(args)
    # Run the appropriate function based on the subcommand
    args.func(args)

def parse_main_args(args):
    telotracker_description = f"""R|{get_ascii_art()}\n
    TeloTracker: a tool to track telomere/chromosome end structural changes from ONT sequencing data."""
    #parser = MyParser(description=tt_descr, formatter_class=MyHelpFormatter, add_help=False)

    # Create the main parser
    parser = argparse.ArgumentParser(
        description=telotracker_description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-v', '--version', action='version', version=f'TeloTracker v{__version__}')

    # Add subparsers for different commands
    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name', required=True)

    # Subcommands for TeloTracker - each subcommand has its own parser
    basecall_subparser(subparsers)
    reference_subparser(subparsers)
    track_subparser(subparsers)

    # Add help and version options
    #help_args = parser.add_argument_group('Help')
    #help_args.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    #help_args.add_argument('-v', '--version', action='version',
    #                       version=f'TeloTracker v{__version__}',
    #                      help='Show program version and exit')
    
    return parser.parse_args(args)

def validate_args(args):
    """
    Validate the parsed arguments.
    
    Raises:
        ValueError: If required arguments are missing or invalid.
    """
    if args.subparser_name == 'basecall':
        if not args.pod5 or not args.outdir:
            raise ValueError("Missing required arguments for basecalling. Use -p/--pod5 and -o/--outdir.")
    
    elif args.subparser_name == 'reference':
        if not args.reads or not args.outdir:
            raise ValueError("Missing required arguments for reference creation. Use -r/--reads and -o/--outdir.")
    
    elif args.subparser_name == 'track':
        if not args.input_file or not args.reference or not args.outdir:
            raise ValueError("Missing required arguments for tracking. Use -i/--input_file, -ref/--reference, and -o/--outdir.")

def basecall_subparser(subparsers):
    subcommand = subparsers.add_parser('basecall', description='Run basecalling for TeloTracker (without trimming)',
                                            help='Run basecalling for TeloTracker (without trimming)')
    
    required_args = subcommand.add_argument_group('Required arguments')
    required_args.add_argument('-p', '--pod5', type=str, required=True,
                               help='Raw input file for basecalling')
    required_args.add_argument('-o', '--outdir', type=pathlib.Path, required=True,
                               help='Output directory for basecalled reads')
    
    settings_args = subcommand.add_argument_group('Settings')
    settings_args.add_argument('-m', '--model', type=str, default='sup',
                               help="""
                                Basecalling model to use\n
                                (e.g. 'sup', 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0', etc.)
                                """)
    settings_args.add_argument('-k', '--kit-name', default="",
                            help="""
                            Optional: If sequencing was done with barcoding, provided the nanopore kit name. Choose from:\n
                            SQK-NBD114-24 SQK-NBD114-96 SQK-RBK114-24 SQK-RBK114-96 SQK-MLK114-96-XL\n
                            * Note: Other kits are not compatible with TeloTracker analysis.
                            """)
    settings_args.add_argument('--threads', type=int, default=4, help='Number of threads')
    settings_args.add_argument(
        '--resume-from', choices=['basecall', 'demultiplex'],
        help="""
        Resume from a specific step.\n
        Options: basecall, trim_reads\n
        1. basecall: Resume from the basecalling step.
        2. demultiplex: Resume from the demultiplexing step.
        """
    )   
        
    subcommand.set_defaults(func=basecall)

def reference_subparser(subparsers):
    subcommand = subparsers.add_parser('reference', description='Generates a telomere reference from input reads.', help='Generate a telomere reference')
    
    required_args = subcommand.add_argument_group('Required arguments')
    required_args.add_argument('-r', '--reads', type=str, required=True,
                               help='Input FASTQ reads')
    required_args.add_argument('-o', '--outdir', type=pathlib.Path, required=True,
                               help='Output directory')
    
    settings_args = subcommand.add_argument_group('Settings')
    settings_args.add_argument('-t', '--threads', type=int, default=4,
                                help='Number of threads')
    settings_args.add_argument('-e', '--existing_reference', type=pathlib.Path, default="",
                                help='Optional existing reference FASTA file')
    settings_args.add_argument(
        '--resume-from', choices=[
            'prep_reads', 'create_reference', 'map_reads', 'polish_reference', 'fix_telomeres', 'create_reference_annotations'
            ],
        help="""
        Resume from a specific step.\n
        Options: prep_reads, create_reference, map_reads, polish_reference, fix_telomeres, create_reference_annotations
        1. prep_reads: Prepare reads for reference creation.
        2. create_reference: Create the initial telomere reference.
        3. map_reads: Map reads to the reference.
        4. polish_reference: Polish the reference using mapped reads.
        5. fix_telomeres: Fix telomere sequences in the reference.
        6. create_reference_annotations: Create annotations for the reference.
        """
    )
    
    subcommand.set_defaults(func=reference_creation_pipeline)
    

def track_subparser(subparsers):
    subcommand = subparsers.add_parser('track', help='track data using an existing reference')
    
    required_args = subcommand.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--input_file', type=str, required=True,
                               help='Input FASTQ of reads to analyze/track telomeres')
    required_args.add_argument('-ref', '--reference', type=str, required=True,
                               help='Reference FASTA')
    required_args.add_argument('-o', '--outdir', type=pathlib.Path, required=True,
                               help='Output directory')
    
    settings_args = subcommand.add_argument_group('Settings')
    settings_args.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    settings_args.add_argument('--subtelomere_reference', required=False, help='Optional subtelomere reference file')
    settings_args.add_argument('--resume', action='store_true', help='Resume skipped steps')
    settings_args.add_argument(
        '--resume-from', choices=[
            'anchor_reads', 'detect_adapter', 'determine_telomere_lengths', 'count_y_primes', 'annotate_subtelomeres', 'call_recombinations'
            ],
        help="""
        Resume from a specific step.\n
        Options: anchor_reads, detect_adapter, determine_telomere_lengths, count_y_primes, annotate_subtelomeres, call_recombinations
        1. anchor_reads: Anchor reads to the reference.
        2. detect_adapter: Detect adapter sequences in the reads.
        3. determine_telomere_lengths: Determine telomere lengths from the reads.
        4. count_y_primes: Count Y-prime elements in the reads.
        5. annotate_subtelomeres: Annotate subtelomere regions in the reference.
        6. call_recombinations: Call recombination events in the reads.
        """
    )
    subcommand.set_defaults(func=track_telomeres_pipeline)


if __name__ == '__main__':
    main()
