"""
TeloTracker Utilities
Common utility functions for both reference and track modules
"""

import os
import sys
import subprocess
import logging
from pathlib import Path
from typing import Optional, Union, List, Dict, Any
import gzip
import mimetypes


def run_command(
    command: Union[str, List[str]], 
    shell: bool = False,
    capture_output: bool = False,
    cwd: Optional[str] = None,
    timeout: Optional[int] = None
) -> subprocess.CompletedProcess:
    """
    Run a system command with proper error handling and logging.
    
    Args:
        command: Command to run (string or list)
        shell: Whether to run through shell
        capture_output: Whether to capture stdout/stderr
        cwd: Working directory
        timeout: Timeout in seconds
        
    Returns:
        CompletedProcess object
    """
    
    logger = logging.getLogger(__name__)
    
    try:
        if isinstance(command, list):
            cmd_str = ' '.join(command)
        else:
            cmd_str = command
        
        logger.debug(f"Running command: {cmd_str}")
        
        result = subprocess.run(
            command,
            shell=shell,
            capture_output=capture_output,
            text=True,
            cwd=cwd,
            timeout=timeout,
            check=False
        )
        
        if result.returncode != 0:
            logger.warning(f"Command failed with return code {result.returncode}: {cmd_str}")
            if capture_output and result.stderr:
                logger.warning(f"Error output: {result.stderr}")
        
        return result
        
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout} seconds: {cmd_str}")
        raise
    except Exception as e:
        logger.error(f"Command execution failed: {cmd_str}, Error: {str(e)}")
        raise


def create_output_directory(path: Union[str, Path], force: bool = False) -> bool:
    """
    Create output directory with proper error handling.
    
    Args:
        path: Directory path to create
        force: Whether to overwrite existing directory
        
    Returns:
        bool: True if successful
    """
    
    logger = logging.getLogger(__name__)
    path = Path(path)
    
    try:
        if path.exists():
            if not force:
                logger.error(f"Directory already exists: {path}")
                return False
            else:
                logger.warning(f"Removing existing directory: {path}")
                import shutil
                shutil.rmtree(path)
        
        path.mkdir(parents=True, exist_ok=False)
        logger.info(f"Created directory: {path}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create directory {path}: {str(e)}")
        return False


def validate_fasta(file_path: Union[str, Path]) -> bool:
    """
    Validate FASTA file format.
    
    Args:
        file_path: Path to FASTA file
        
    Returns:
        bool: True if valid FASTA format
    """
    
    logger = logging.getLogger(__name__)
    file_path = Path(file_path)
    
    try:
        if not file_path.exists():
            logger.error(f"FASTA file not found: {file_path}")
            return False
        
        # Check file extension
        valid_extensions = {'.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn'}
        if file_path.suffix.lower() not in valid_extensions:
            logger.warning(f"Unexpected FASTA file extension: {file_path.suffix}")
        
        # Check if file is gzipped
        opener = gzip.open if file_path.suffix == '.gz' else open
        mode = 'rt' if file_path.suffix == '.gz' else 'r'
        
        # Read first few lines to validate format
        with opener(file_path, mode) as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                logger.error(f"FASTA file doesn't start with '>': {file_path}")
                return False
            
            # Check that we have sequence data
            sequence_found = False
            for i, line in enumerate(f):
                if i > 10:  # Check first 10 lines
                    break
                line = line.strip()
                if line and not line.startswith('>'):
                    sequence_found = True
                    # Basic sequence validation
                    valid_chars = set('ATCGNRYKMSWBDHV-')  # Standard IUPAC codes
                    if not all(c.upper() in valid_chars for c in line):
                        logger.warning(f"Non-standard characters in FASTA: {file_path}")
            
            if not sequence_found:
                logger.error(f"No sequence data found in FASTA: {file_path}")
                return False
        
        logger.debug(f"FASTA validation passed: {file_path}")
        return True
        
    except Exception as e:
        logger.error(f"FASTA validation failed for {file_path}: {str(e)}")
        return False


def validate_fastq(file_path: Union[str, Path]) -> bool:
    """
    Validate FASTQ file format.
    
    Args:
        file_path: Path to FASTQ file
        
    Returns:
        bool: True if valid FASTQ format
    """
    
    logger = logging.getLogger(__name__)
    file_path = Path(file_path)
    
    try:
        if not file_path.exists():
            logger.error(f"FASTQ file not found: {file_path}")
            return False
        
        # Check file extension
        valid_extensions = {'.fastq', '.fq', '.fastq.gz', '.fq.gz'}
        if not any(str(file_path).endswith(ext) for ext in valid_extensions):
            logger.warning(f"Unexpected FASTQ file extension: {file_path}")
        
        # Check if file is gzipped
        opener = gzip.open if str(file_path).endswith('.gz') else open
        mode = 'rt' if str(file_path).endswith('.gz') else 'r'
        
        # Validate FASTQ format (4 lines per record)
        with opener(file_path, mode) as f:
            lines_checked = 0
            while lines_checked < 20:  # Check first 5 records (20 lines)
                # Header line
                header = f.readline()
                if not header:
                    break
                if not header.startswith('@'):
                    logger.error(f"Invalid FASTQ header: {file_path}")
                    return False
                
                # Sequence line
                sequence = f.readline()
                if not sequence:
                    logger.error(f"Incomplete FASTQ record: {file_path}")
                    return False
                
                # Plus line
                plus = f.readline()
                if not plus or not plus.startswith('+'):
                    logger.error(f"Invalid FASTQ plus line: {file_path}")
                    return False
                
                # Quality line
                quality = f.readline()
                if not quality:
                    logger.error(f"Missing quality line: {file_path}")
                    return False
                
                # Check that sequence and quality have same length
                if len(sequence.strip()) != len(quality.strip()):
                    logger.error(f"Sequence/quality length mismatch: {file_path}")
                    return False
                
                lines_checked += 4
        
        logger.debug(f"FASTQ validation passed: {file_path}")
        return True
        
    except Exception as e:
        logger.error(f"FASTQ validation failed for {file_path}: {str(e)}")
        return False


def check_dependencies(required_tools: List[str]) -> Dict[str, bool]:
    """
    Check if required command-line tools are available.
    
    Args:
        required_tools: List of tool names to check
        
    Returns:
        Dict mapping tool names to availability status
    """
    
    logger = logging.getLogger(__name__)
    results = {}
    
    for tool in required_tools:
        try:
            result = subprocess.run(
                ['which', tool],
                capture_output=True,
                text=True
            )
            available = result.returncode == 0
            results[tool] = available
            
            if available:
                logger.debug(f"Found tool: {tool}")
            else:
                logger.warning(f"Tool not found: {tool}")
                
        except Exception as e:
            logger.error(f"Error checking tool {tool}: {str(e)}")
            results[tool] = False
    
    return results


def get_file_info(file_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Get information about a file.
    
    Args:
        file_path: Path to file
        
    Returns:
        Dictionary with file information
    """
    
    file_path = Path(file_path)
    
    info = {
        'path': str(file_path),
        'exists': file_path.exists(),
        'size': None,
        'is_gzipped': False,
        'mime_type': None
    }
    
    if file_path.exists():
        stat = file_path.stat()
        info['size'] = stat.st_size
        info['modified'] = stat.st_mtime
        info['is_gzipped'] = str(file_path).endswith('.gz')
        
        # Try to determine MIME type
        mime_type, _ = mimetypes.guess_type(str(file_path))
        info['mime_type'] = mime_type
    
    return info


def format_file_size(size_bytes: int) -> str:
    """
    Format file size in human-readable format.
    
    Args:
        size_bytes: Size in bytes
        
    Returns:
        Formatted size string
    """
    
    if size_bytes == 0:
        return "0 B"
    
    size_names = ["B", "KB", "MB", "GB", "TB"]
    i = 0
    size = float(size_bytes)
    
    while size >= 1024.0 and i < len(size_names) - 1:
        size /= 1024.0
        i += 1
    
    return f"{size:.1f} {size_names[i]}"


def parse_memory_string(memory_str: str) -> int:
    """
    Parse memory string (e.g., '8G', '512M') to bytes.
    
    Args:
        memory_str: Memory string
        
    Returns:
        Memory in bytes
    """
    
    memory_str = memory_str.upper().strip()
    
    if memory_str.endswith('G'):
        return int(float(memory_str[:-1])
