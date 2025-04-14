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

import argparse
import sys
import pathlib

from .version import __version__

def main():
    """
    Entry point for the command line interface.
    """
    # Create the argument parser.
    parser = argparse.ArgumentParser(description='telotracker: a description of what your tool does',
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Add subparsers to handle different subcommands.
    subparsers = parser.add_subparsers(dest='subparser_name', required=True)
    
    # Example subcommand 1: 'function1'
    subparser_function1 = subparsers.add_parser('function1', help='Description of function1',
                                               description='Function1 description - more details',
                                               formatter_class=argparse.RawTextHelpFormatter)
    subparser_function1.add_argument('--input', type=str, required=True,
                                     help='Input file')
    subparser_function1.add_argument('--output', type=str, required=True,
                                     help='Output file')

    # You can define additional subcommands here
    # Example subcommand 2: 'function2'
    # subparser_function2 = subparsers.add_parser('function2', help='Description of function2')
    # subparser_function2.add_argument('param', type=str, help='Parameter for function2')

    # Help and version options for main parser
    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    help_args.add_argument('-v', '--version', action='version',
                           version='telotracker v' + __version__,
                           help='Show program version and exit')

    # Parse arguments
    args = parser.parse_args()

    # Handle subcommands and invoke corresponding function
    if args.subparser_name == 'function1':
        from . import function1
        function1.main(args)
    
    # Add more subcommand logic as needed
    # elif args.subparser_name == 'function2':
    #     from . import function2
    #     function2.main(args)

    else:
        # Handle unsupported subcommands
        print(f"Error: Unsupported subcommand '{args.subparser_name}'")
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
