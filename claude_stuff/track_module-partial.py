"""
TeloTracker Track Module
Handles sequencing data analysis for telomere tracking
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from .utils import run_command, create_output_directory, validate_fastq
from .anchor_reads import ReadAnchor
from .circle_search import CircleSearcher
from .find_subtelomeric_recombination import RecombinationFinder


class TrackModule:
    """
    Main class for sequencing data analysis and telomere tracking.
    
    This module handles:
    1. Read mapping to reference
    2. Telomere length analysis
    3. Subtelomeric recombination detection
    4. Circular DNA analysis
    5. Feature quantification and reporting
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.output_dir = None
        self.reference_dir = None
        self.read_files = []
        self.config = {}
        self.reference_info = {}
        
    def run(self, args) -> bool:
        """
        Main execution method for track module.
        
        Args:
            args: Parsed command-line arguments from argparse
            
        Returns:
            bool: True if successful, False otherwise
        """
        
        try:
            # Store arguments and setup
            self.output_dir = Path(args.output)
            self.reference_dir = Path(args.reference)
            self.read_files = args.reads
            
            # Create configuration dictionary
            self.config = {
                'reference_directory': str(self.reference_dir),
                'read_files': self.read_files,
                'sample_name': args.sample_name or self._infer_sample_name(),
                'threads': args.threads,
                'memory': args.memory,
                'min_mapq': args.min_mapq,
                'output_directory': str(self.output_dir),
                'force_overwrite': args.force,
                'skip_mapping': args.skip_mapping,
                'analyze_circles': args.analyze_circles,
                'find_recombination': args.find_recombination
            }
            
            self.logger.info(f"Analyzing sample: {self.config['sample_name']}")
            self.logger.info(f"Output directory: {self.output_dir}")
            
            # Step 1: Setup and validation
            if not self._setup_analysis():
                return False
            
            # Step 2: Read mapping (unless skipped)
            if not self.config['skip_mapping']:
                if not self._map_reads():
                    return False
            
            # Step 3: Core telomere analysis
            if not self._analyze_telomeres():
                return False
            
            # Step 4: Optional analyses
            if not self._run_optional_analyses():
                return False
            
            # Step 5: Generate reports
            if not self._generate_reports():
                return False
            
            self.logger.info("Track module completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Track module failed: {str(e)}")
            return False
    
    def _setup_analysis(self) -> bool:
        """Setup analysis environment and validate inputs."""
        
        try:
            # Create output directory structure
            if not self._setup_output_directory():
                return False
            
            # Load reference information
            if not self._load_reference_info():
                return False
            
            # Validate input files
            if not self._validate_inputs():
                return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"Analysis setup failed: {str(e)}")
            return False
    
    def _setup_output_directory(self) -> bool:
        """Setup output directory structure."""
        
        try:
            # Handle existing directory
            if self.output_dir.exists() and not self.config['force_overwrite']:
                self.logger.error(f"Output directory exists: {self.output_dir}")
                return False
            
            if self.config['force_overwrite'] and self.output_dir.exists():
                import shutil
                shutil.rmtree(self.output_dir)
            
            # Create directory structure
            self.output_dir.mkdir(parents=True, exist_ok=False)
            
            subdirs = [
                'mapping',           # Read mapping results
                'telomeres',        # Telomere analysis results
                'features',         # Feature quantification
                'circles',          # Circular DNA analysis
                'recombination',    # Recombination analysis
                'reports',          # Final reports and summaries
                'logs',             # Analysis logs
                'temp'              # Temporary files
            ]
            
            for subdir in subdirs:
                (self.output_dir / subdir).mkdir(exist_ok=True)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to setup output directory: {str(e)}")
            return False
    
    def _load_reference_info(self) -> bool:
        """Load reference information and validate reference directory."""
        
        try:
            # Check for reference info file
            ref_info_file = self.reference_dir / 'reference_info.json'
            if not ref_info_file.exists():
                self.logger.error(f"Reference info file not found: {ref_info_file}")
                return False
            
            # Load reference metadata
            with open(ref_info_file, 'r') as f:
                self.reference_info = json.load(f)
            
            # Validate required reference files
            required_files = [
                'genome.fasta',
                'annotations/telomeres.bed',
                'indices/genome.mmi'  # or appropriate index files
