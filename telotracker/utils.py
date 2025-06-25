from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary with headers as keys and sequences as values.
    Uses Biopython for parsing.
    """
    sequences = {}
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = str(record.seq)
    return sequences
    

def format_fasta(sequence, line_length=60):
    formatted_sequence = []
    for i in range(0, len(sequence), line_length):
        formatted_sequence.append(sequence[i:i+line_length])
    return '\n'.join(formatted_sequence)


def extract_fasta_from_bed(fasta_file, bed_file, output_fasta):
    sequences = {}
    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            sequences[record.id] = record.seq

    extracted_sequences = {}
    with open(bed_file, "r") as bed_handle:
        for line in bed_handle:
            print(line.strip().split("\t"))
            chromosome, start, end, identifier, strand, length = line.strip().split("\t")
            start, end = int(start), int(end)
            sequence = sequences[chromosome][start-1:end]
            if strand == "-":
                sequence = sequence.reverse_complement()
            extracted_sequences[identifier] = sequence

    # Write extracted sequences to a new output FASTA file
    with open(output_fasta, "w") as output_handle:
        for identifier, sequence in extracted_sequences.items():
            record = SeqRecord(sequence, id=identifier, description="")
            SeqIO.write(record, output_handle, "fasta")
            
def check_file_type(file_path):
    """
    Checks the file type based on its extension.
    """
    file_types = {
        "bam": (".bam",),
        "fastq": (".fastq", ".fastq.gz"),
        "fasta": (".fasta", ".fa"),
        "bed": (".bed",),
    }

    for file_type, extensions in file_types.items():
        if file_path.endswith(extensions):
            return file_type

    raise ValueError(
        "Unsupported file type. Please provide a .fasta/.fa, .fastq/.fastq.gz, .bam, or .bed file."
    )