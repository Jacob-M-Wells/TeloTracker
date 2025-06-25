"""
TeloTracker Reference Module
Handles reference genome processing and annotation creation
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from .utils import run_command, create_output_directory, validate_fasta
from .get_reference_features import FeatureAnnotator
from .create_reference import ReferenceBuilder


class ReferenceModule:
    """
    Main class for reference genome processing and annotation.
    
    This module handles:
    1. Genome indexing and validation
    2. Telomere sequence annotation
    3. X-element identification
    4. Y-prime sequence annotation
    5. Creation of analysis-ready reference files
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.output_dir = None
        self.genome_file = None
        self.config = {}
        
    def run(self, args) -> bool:
        """
        Main execution method for reference module.
        
        Args:
            args: Parsed command-line arguments from argparse
            
        Returns:
            bool: True if successful, False otherwise
        """
        
        try:
            # Store arguments and setup
            self.output_dir = Path(args.output)
            self.genome_file = args.genome
            
            # Create configuration dictionary
            self.config = {
                'genome_file': args.genome,
                'telomere_sequence': args.telomere_seq,
                'min_telomere_length': args.min_telomere_length,
                'threads': args.threads,
                'output_directory': str(self.output_dir),
                'force_overwrite': args.force
            }
            
            self.logger.info(f"Setting up reference in: {self.output_dir}")
            
            # Step 1: Setup output directory and validate input
            if not self._setup_output_directory():
                return False
                
            if not self._validate_input():
                return False
            
            # Step 2: Process genome file
            if not self._process_genome():
                return False
            
            # Step 3: Create annotations
            if not self._create_annotations():
                return False
            
            # Step 4: Generate summary and metadata
            if not self._generate_summary():
                return False
            
            self.logger.info("Reference module completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Reference module failed: {str(e)}")
            return False
    
    def _setup_output_directory(self) -> bool:
        """Setup and validate output directory."""
        
        try:
            # Create output directory
            if self.config['force_overwrite'] and self.output_dir.exists():
                self.logger.warning(f"Removing existing directory: {self.output_dir}")
                import shutil
                shutil.rmtree(self.output_dir)
            
            self.output_dir.mkdir(parents=True, exist_ok=False)
            
            # Create subdirectories
            subdirs = ['annotations', 'indices', 'logs', 'temp']
            for subdir in subdirs:
                (self.output_dir / subdir).mkdir(exist_ok=True)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to setup output directory: {str(e)}")
            return False
    
    def _validate_input(self) -> bool:
        """Validate input genome file."""
        
        try:
            if not os.path.exists(self.genome_file):
                self.logger.error(f"Genome file not found: {self.genome_file}")
                return False
            
            # Validate FASTA format
            if not validate_fasta(self.genome_file):
                self.logger.error(f"Invalid FASTA format: {self.genome_file}")
                return False
            
            self.logger.info(f"Input genome validated: {self.genome_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Input validation failed: {str(e)}")
            return False
    
    def _process_genome(self) -> bool:
        """Process genome file - copy, index, and prepare for analysis."""
        
        try:
            self.logger.info("Processing genome file...")
            
            # Copy genome to output directory
            genome_output = self.output_dir / 'genome.fasta'
            
            # Use symlink if possible, copy otherwise
            try:
                genome_output.symlink_to(os.path.abspath(self.genome_file))
                self.logger.info(f"Created symlink to genome: {genome_output}")
            except OSError:
                import shutil
                shutil.copy2(self.genome_file, genome_output)
                self.logger.info(f"Copied genome to: {genome_output}")
            
            # Create genome index for mapping tools (BWA, minimap2, etc.)
            self._create_genome_indices(genome_output)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Genome processing failed: {str(e)}")
            return False
    
    def _create_genome_indices(self, genome_file: Path) -> bool:
        """Create genome indices for various mapping tools."""
        
        indices_dir = self.output_dir / 'indices'
        
        # List of indexing commands to run
        index_commands = []
        
        # BWA index
        bwa_cmd = f"bwa index -p {indices_dir}/genome {genome_file}"
        index_commands.append(("BWA", bwa_cmd))
        
        # Minimap2 index
        minimap2_cmd = f"minimap2 -d {indices_dir}/genome.mmi {genome_file}"
        index_commands.append(("Minimap2", minimap2_cmd))
        
        # FASTA index
        faidx_cmd = f"samtools faidx {genome_file}"
        index_commands.append(("FASTA", faidx_cmd))
        
        # Run indexing commands in parallel
        with ThreadPoolExecutor(max_workers=min(3, self.config['threads'])) as executor:
            futures = {}
            
            for name, cmd in index_commands:
                future = executor.submit(self._run_indexing_command, name, cmd)
                futures[future] = name
            
            # Wait for completion
            for future in as_completed(futures):
                name = futures[future]
                try:
                    success = future.result()
                    if success:
                        self.logger.info(f"{name} index created successfully")
                    else:
                        self.logger.warning(f"{name} index creation failed")
                except Exception as e:
                    self.logger.error(f"{name} index creation error: {str(e)}")
        
        return True
    
    def _run_indexing_command(self, name: str, command: str) -> bool:
        """Run a single indexing command."""
        try:
            result = run_command(command, capture_output=True)
            return result.returncode == 0
        except Exception as e:
            self.logger.error(f"{name} indexing failed: {str(e)}")
            return False
    
    def _create_annotations(self) -> bool:
        """Create all telomere-related annotations."""
        
        try:
            self.logger.info("Creating annotations...")
            
            genome_file = self.output_dir / 'genome.fasta'
            annotations_dir = self.output_dir / 'annotations'
            
            # Initialize feature annotator
            annotator = FeatureAnnotator(
                genome_file=str(genome_file),
                output_dir=str(annotations_dir),
                config=self.config
            )
            
            # Create different types of annotations
            annotation_tasks = [
                ("telomeres", self._annotate_telomeres),
                ("x_elements", self._annotate_x_elements),
                ("y_primes", self._annotate_y_primes),
                ("subtelomeric_regions", self._annotate_subtelomeric_regions)
            ]
            
            # Run annotation tasks
            for task_name, task_func in annotation_tasks:
                self.logger.info(f"Creating {task_name} annotations...")
                if not task_func(annotator):
                    self.logger.error(f"Failed to create {task_name} annotations")
                    return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"Annotation creation failed: {str(e)}")
            return False
    
    def _annotate_telomeres(self, annotator) -> bool:
        """Annotate telomere sequences."""
        try:
            telomere_bed = annotator.find_telomeres(
                telomere_seq=self.config['telomere_sequence'],
                min_length=self.config['min_telomere_length']
            )
            
            output_file = self.output_dir / 'annotations' / 'telomeres.bed'
            annotator.write_bed_file(telomere_bed, output_file)
            
            self.logger.info(f"Found {len(telomere_bed)} telomere regions")
            return True
            
        except Exception as e:
            self.logger.error(f"Telomere annotation failed: {str(e)}")
            return False
    
    def _annotate_x_elements(self, annotator) -> bool:
        """Annotate X-element sequences."""
        try:
            # Implementation depends on your specific X-element detection logic
            x_elements = annotator.find_x_elements()
            
            output_file = self.output_dir / 'annotations' / 'x_elements.bed'
            annotator.write_bed_file(x_elements, output_file)
            
            self.logger.info(f"Found {len(x_elements)} X-element regions")
            return True
            
        except Exception as e:
            self.logger.error(f"X-element annotation failed: {str(e)}")
            return False
    
    def _annotate_y_primes(self, annotator) -> bool:
        """Annotate Y-prime sequences."""
        try:
            y_primes = annotator.find_y_primes()
            
            output_file = self.output_dir / 'annotations' / 'y_primes.bed'
            annotator.write_bed_file(y_primes, output_file)
            
            self.logger.info(f"Found {len(y_primes)} Y-prime regions")
            return True
            
        except Exception as e:
            self.logger.error(f"Y-prime annotation failed: {str(e)}")
            return False
    
    def _annotate_subtelomeric_regions(self, annotator) -> bool:
        """Annotate subtelomeric regions."""
        try:
            subtel_regions = annotator.find_subtelomeric_regions()
            
            output_file = self.output_dir / 'annotations' / 'subtelomeric_regions.bed'
            annotator.write_bed_file(subtel_regions, output_file)
            
            self.logger.info(f"Found {len(subtel_regions)} subtelomeric regions")
            return True
            
        except Exception as e:
            self.logger.error(f"Subtelomeric annotation failed: {str(e)}")
            return False
    
    def _generate_summary(self) -> bool:
        """Generate summary statistics and metadata."""
        
        try:
            self.logger.info("Generating reference summary...")
            
            # Collect statistics
            stats = self._collect_reference_stats()
            
            # Create metadata file
            metadata = {
                'telotracker_version': self._get_version(),
                'creation_date': self._get_timestamp(),
                'configuration': self.config,
                'statistics': stats,
                'files': self._list_output_files()
            }
            
            # Write metadata
            metadata_file = self.output_dir / 'reference_info.json'
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            # Write summary report
            self._write_summary_report(stats)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Summary generation failed: {str(e)}")
            return False
    
    def _collect_reference_stats(self) -> Dict:
        """Collect statistics about the created reference."""
        
        stats = {}
        annotations_dir = self.output_dir / 'annotations'
        
        # Count features in each annotation file
        for bed_file in annotations_dir.glob('*.bed'):
            feature_type = bed_file.stem
            try:
                with open(bed_file, 'r') as f:
                    count = sum(1 for line in f if not line.startswith('#'))
                stats[f'{feature_type}_count'] = count
            except Exception:
                stats[f'{feature_type}_count'] = 0
        
        return stats
    
    def _write_summary_report(self, stats: Dict) -> None:
        """Write human-readable summary report."""
        
        report_file = self.output_dir / 'REFERENCE_SUMMARY.txt'
        
        with open(report_file, 'w') as f:
            f.write("TeloTracker Reference Summary\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Genome file: {self.config['genome_file']}\n")
            f.write(f"Output directory: {self.config['output_directory']}\n")
            f.write(f"Creation date: {self._get_timestamp()}\n\n")
            
            f.write("Annotation Statistics:\n")
            f.write("-" * 20 + "\n")
            for key, value in stats.items():
                f.write(f"{key.replace('_', ' ').title()}: {value}\n")
            
            f.write(f"\nFiles created:\n")
            for file_path in self._list_output_files():
                f.write(f"  - {file_path}\n")
    
    def _get_version(self) -> str:
        """Get TeloTracker version."""
        try:
            from .version import __version__
            return __version__
        except ImportError:
            return "unknown"
    
    def _get_timestamp(self) -> str:
        """Get current timestamp."""
        from datetime import datetime
        return datetime.now().isoformat()
    
    def _list_output_files(self) -> List[str]:
        """List all created output files."""
        output_files = []
        for file_path in self.output_dir.rglob('*'):
            if file_path.is_file():
                output_files.append(str(file_path.relative_to(self.output_dir)))
        return sorted(output_files)

