<p align="center"><img src="images/telotracker_banner_v1.png" alt="TeloTracker" width="70%"></p>



# TeloTracker
Yeast Telomere Identification in Nanopore Reads

`TeloTracker` is a Python package designed to track and visualize Oxford Nanopore Technology (ONT) Sequencing of yeast telomeres.
It was specifically developed to study Alternative Lengthening of Telomeres (ALT) in yeast, but it can be used for general telomere analysis in *Saccharomyces cerevisiae*.

## 🔧 Installation
To install `TeloTracker`, you can clone this repository by running the command below.

```bash
git clone https://github.com/Jacob-M-Wells/TeloTracker.git
cd TeloTracker
```



### Installing with conda (recommended)
```bash
conda env create -n telotracker -f environment.yml
conda activate telotracker
```

### Installing with pip
```bash
pip install -r requirements.txt
```
**Note:** Several required tools are not available on PyPI and must be installed separately if using pip.
Using conda is highly recommended to avoid missing dependencies.

To avoid missing dependencies, **installing with Conda is strongly recommended (see above).**

The following external tools are required:

- [`BLAST+`](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (version 2.16.0)
- [`RepeatMasker`](http://www.repeatmasker.org/) (version 4.1.7-p1)
- [`Porechop_ABI`](https://github.com/bonsai-team/Porechop_ABI) (≥0.5.0)
- [`Clustal Omega`](http://www.clustal.org/omega/) (version 1.2.4)
- [`Flye`](https://github.com/fenderglass/Flye) (version 2.9.5)
- [`Minimap2`](https://github.com/lh3/minimap2) (version 2.28)
- [`QUAST`](http://quast.sourceforge.net/) (version 5.3.0)
- [`Medaka`](https://github.com/nanoporetech/medaka) (version 2.0.1)
- [`MUMmer`](https://mummer4.github.io/) (version 3.23)
- [`MUSCLE`](https://www.drive5.com/muscle/) (version 5.3)
- [`SAMtools & BCFtools`](http://www.htslib.org/) (version 1.21)
- [`BEDtools`](https://bedtools.readthedocs.io/) (version 2.31.1)
- [`RMBlast`](http://www.repeatmasker.org/RMBlast.html) (version 2.14.1)

## Installing Dorado (Required for Basecalling)

[Dorado](https://github.com/nanoporetech/dorado) is a high-performance, open source basecaller for Oxford Nanopore reads. Dorado can be run with a docker container or installed locally.

Dorado must be used to basecall (rather than using MinKnow), as basecalling without trimming adapters/barcode sequences is required.

📝 Detailed instructions and usage examples are available in the [Dorado documentation](https://dorado-docs.readthedocs.io/en/latest/).

#### Option 1. Run with Docker Container (using singularity)
📦 Download the latest release for your platform from the [Dorado Docker Releases page](https://hub.docker.com/r/nanoporetech/dorado/tags).
```bash
singularity pull dorado.sif docker://nanoporetech/dorado:latest
singularity exec dorado.sif dorado --help
```

#### Option 2. Run with local installation
📦 Download the latest release for your platform from the [Dorado GitHub Releases page](https://github.com/nanoporetech/dorado/releases).
The current version as of 03/13/2025: Dorado 0.9.1 @ https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz
```bash
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz
gunzip https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz
tar -xvf https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz
mv dorado-0.9.1-linux-x64/ ~/bin/
ln -s ~/bin/dorado-0.9.1-linux-x64/bin/dorado ~/bin/dorado
source ~/.bashrc
dorado --help
```

## Usage

TeloTracker processes ONT waveform data (.pod5/.fast5 files) to basecall and analyze reads for their telomere and sub-telomere features. TeloTracker has two distinct pipelines: `reference` and `track`.

- reference: creates a telomere reference genome and annotation
- track: analyzes reads against a reference for telomere dynamics

## Quick Start

1. **Generate a reference**
   
    telotracker reference \
        -r input_reads.fastq \
        -o teloref_outdir \
        -t 16

2. **Track telomeres in new reads**

    telotracker track \
        -i new_reads.fastq \
        -ref teloref_outdir/reference.fasta \
        -o tracking_outdir \
        -t 16

3. **Optional: Basecall and demultiplex from .pod5**

    telotracker basecall \
        -p raw_input.pod5 \
        -o basecall_outdir \
        -m sup \
        -t 8

## Citing

If you use TeloTracker in your research, please cite this GitHub repository:

    Wells, Jacob M. TeloTracker. https://github.com/Jacob-M-Wells/TeloTracker

