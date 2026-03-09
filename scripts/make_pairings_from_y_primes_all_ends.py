from Bio import SeqIO
import pandas as pd
import os
import sys

# Usage: make_pairings_from_y_primes_all_ends.py <base_name> <anchor_set> <output_dir> <strain_number>
base_name  = sys.argv[1]
anchor_set = sys.argv[2]
output_dir = sys.argv[3]
strain_id  = sys.argv[4]

print("Starting make_pairings_from_y_primes_all_ends.py")
print(f'Opening {base_name}...')

blast_dir             = os.path.join(output_dir, 'blast')
chr_anchor_reads_dir  = os.path.join(blast_dir, 'chr_anchor_reads')
pairing_output_dir    = os.path.join(output_dir, 'paired_by_y_prime_reads')
good_end_rm_tsv       = os.path.join(output_dir, f'{base_name}_good_end_y_repeatmasker.tsv')

os.makedirs(pairing_output_dir, exist_ok=True)

if strain_id == '6991':
    chr_end_first_y_prime_id_dict = {
        'chr1L': None,  'chr1R': None,
        'chr2L': 'ID4', 'chr2R': None,
        'chr3L': None,  'chr3R': None,
        'chr4L': None,  'chr4R': 'ID2',
        'chr5L': 'ID6', 'chr5R': 'ID1',
        'chr6L': 'ID4', 'chr6R': None,
        'chr7L': None,  'chr7R': 'ID5',
        'chr8L': 'ID1', 'chr8R': 'ID1',
        'chr9L': 'ID6', 'chr9R': None,
        'chr10L': 'ID6', 'chr10R': None,
        'chr11L': None, 'chr11R': None,
        'chr12L': 'ID1', 'chr12R': 'ID1',
        'chr13L': 'ID1', 'chr13R': None,
        'chr14L': 'ID5', 'chr14R': 'ID6',
        'chr15L': None, 'chr15R': 'ID2',
        'chr16L': 'ID5', 'chr16R': 'ID1'
    }
elif strain_id == '7172':
    chr_end_first_y_prime_id_dict = {
        'chr1L': None,  'chr1R': None,
        'chr2L': 'ID4', 'chr2R': None,
        'chr3L': None,  'chr3R': None,
        'chr4L': None,  'chr4R': 'ID2',
        'chr5L': 'ID6', 'chr5R': 'ID1',
        'chr6L': 'ID4', 'chr6R': None,
        'chr7L': None,  'chr7R': 'ID5',
        'chr8L': 'ID1', 'chr8R': 'ID1',
        'chr9L': 'ID6', 'chr9R': None,
        'chr10L': 'ID6', 'chr10R': None,
        'chr11L': None, 'chr11R': None,
        'chr12L': 'ID1', 'chr12R': 'ID1',
        'chr13L': 'ID1', 'chr13R': None,
        'chr14L': 'ID5', 'chr14R': 'ID6',
        'chr15L': None, 'chr15R': 'ID2',
        'chr16L': 'ID5', 'chr16R': 'ID1'
    }
elif strain_id == '7302':
    chr_end_first_y_prime_id_dict = {
        'chr1L': None,  'chr1R': None,
        'chr2L': 'ID4', 'chr2R': None,
        'chr3L': None,  'chr3R': None,
        'chr4L': None,  'chr4R': 'ID2',
        'chr5L': 'ID6', 'chr5R': 'ID1',
        'chr6L': 'ID4', 'chr6R': None,
        'chr7L': None,  'chr7R': 'ID5',
        'chr8L': 'ID1', 'chr8R': 'ID1',
        'chr9L': 'ID6', 'chr9R': None,
        'chr10L': 'ID6', 'chr10R': None,
        'chr11L': None, 'chr11R': None,
        'chr12L': 'ID1', 'chr12R': 'ID1',
        'chr13L': 'ID7', 'chr13R': None,
        'chr14L': 'ID5', 'chr14R': 'ID6',
        'chr15L': None, 'chr15R': 'ID2',
        'chr16L': 'ID5', 'chr16R': 'ID1'
    }
else:
    raise ValueError(f'Unknown strain_id: "{strain_id}" (from base_name: "{base_name}")')


def calculate_y_prime_change(df_read):
    chr_end            = df_read.at[0, 'chr_end']
    ref_first_y_prime  = chr_end_first_y_prime_id_dict[chr_end]
    telomere_side      = df_read.at[0, 'telomere_side']

    if telomere_side == 'beginning':
        read_first_y_prime_color = df_read.at[(len(df_read) - 1), 'y_prime_id_and_color']
    else:
        read_first_y_prime_color = df_read.at[0, 'y_prime_id_and_color']

    read_first_y_prime = read_first_y_prime_color.split('_')[0]

    if read_first_y_prime == ref_first_y_prime:
        return None
    else:
        return f'{chr_end}_and_{read_first_y_prime}'


df = pd.read_csv(good_end_rm_tsv, sep='\t')
df['y_prime_id_and_color'] = df['y_prime_group'].apply(lambda x: x.split('/')[2])

all_reads        = df['read_id'].unique()
y_prime_change_dict = {}
for read_id in all_reads:
    df_read = df[df['read_id'] == read_id].reset_index(drop=True)
    y_prime_change_dict[read_id] = calculate_y_prime_change(df_read)

df['y_prime_change'] = df['read_id'].apply(lambda x: y_prime_change_dict[x])
df = df.dropna(subset=['y_prime_change'])

if df["y_prime_change"].unique().size == 0:
    print(
        "*"*40,  
        "\nNOTICE:\n"
        "No reads with y_prime changes found in the input TSV file.\n"
        "This likely means that all reads have the same y_prime at the relevant chromosome end, or that there are no reads with y_primes at all.\n"
        "Please check the input TSV file and the reference y_prime assignments for each chromosome end.\n"
        "Exiting!\n",
        "*"*40)
    sys.exit(0)

for pairing in df['y_prime_change'].unique():
    df_pair          = df[df['y_prime_change'] == pairing]
    pairing_read_set = set(df_pair['read_id'].unique())

    chr_end = pairing.split('_')[0]
    input_fasta  = os.path.join(chr_anchor_reads_dir,
                                f'{base_name}_blasted_{anchor_set}_{chr_end}_anchor_reads.fasta')
    output_fasta = os.path.join(pairing_output_dir, f'{pairing}.fasta')

    print(f'Creating file: {output_fasta}')

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            if record.description.split(' ')[0] in pairing_read_set:
                SeqIO.write(record, outfile, 'fasta')
