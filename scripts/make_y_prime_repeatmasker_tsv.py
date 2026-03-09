import pandas as pd
import os
import sys

# Usage: make_y_prime_repeatmasker_tsv.py <base_name> <output_dir>
base_name  = sys.argv[1]
output_dir = sys.argv[2]

print("Starting make_y_prime_repeatmasker_tsv.py")
print(f'Opening {base_name}...')

repeatmasker_dir      = os.path.join(output_dir, 'repeatmasker', 'y_primes')
post_y_prime_tsv      = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')

output_combined       = os.path.join(output_dir, f'{base_name}_combined_repeatmasker.tsv')
output_good_end       = os.path.join(output_dir, f'{base_name}_good_end_y_repeatmasker.tsv')
output_gained_y       = os.path.join(output_dir, f'{base_name}_gained_y_repeatmasker.tsv')

all_results_files = [
    os.path.join(repeatmasker_dir, f)
    for f in os.listdir(repeatmasker_dir)
    if f.endswith('_y_prime_repeatmasker.ssv')
]

df_all = pd.DataFrame()

print(f'Processing {len(all_results_files)} repeatmasker result files...')

for rm_file in all_results_files:
    corrected = rm_file.replace('.ssv', '_corrected.ssv')

    with open(rm_file, 'r') as orig:
        corrected_data = []
        for line_num, line in enumerate(orig.readlines()):
            if line_num == 0:
                corrected_data.append(line)
            else:
                line = line.strip()
                if line and line[-1] != '*':
                    line = f'{line} -'
                corrected_data.append(line)

    with open(corrected, 'w') as cf:
        for fixed_line in corrected_data:
            cf.write(f'{fixed_line}\n')

    try:
        df_single = pd.read_csv(corrected, sep=r'\s+')
        # Extract chr_end from filename: <base_name>_<chr_end>_y_prime_repeatmasker_corrected.ssv
        chr_end = os.path.basename(corrected)
        chr_end = chr_end.replace(f'{base_name}_', '').replace('_y_prime_repeatmasker_corrected.ssv', '')
        df_single['chr_end'] = chr_end
        df_all = pd.concat([df_all, df_single])
    except Exception:
        print(f'Error in {corrected}: {sys.exc_info()}')

df_all.loc[df_all['sub_match'] == '*', 'sub_match'] = True
df_all.loc[df_all['sub_match'] == '-', 'sub_match'] = False

print(df_all)
df_all.to_csv(output_combined, sep='\t')

# ---- Filter using post_y_prime_probe results ----
df_filter = pd.read_csv(post_y_prime_tsv, sep='\t')
df_filter = df_filter.dropna(subset=['repeat_length'])
df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

df_strand = df_filter[['read_id', 'Repeat_Type']]
df_all    = df_all.merge(df_strand, how='left', on='read_id')

def determine_telomere_side(chr_end, strand):
    return 'beginning' if strand == 'AC' else 'end'

df_all['telomere_side'] = df_all.apply(
    lambda r: determine_telomere_side(r['chr_end'], r['Repeat_Type']), axis=1)

good_reads = set(df_filter['read_id'])
df_all['good_ends'] = df_all['read_id'].apply(lambda x: x in good_reads)
print(df_all['good_ends'].value_counts())

df_good_end = df_all[df_all['good_ends'] == True]
df_good_end.to_csv(output_good_end, sep='\t')

# ---- Y prime gain subset ----
reads_gain = set(df_filter[df_filter['delta_y_prime_sign'] == '+']['read_id'])
df_all['gained_y'] = df_all['read_id'].apply(lambda x: x in reads_gain)
print(df_all['gained_y'].value_counts())

df_gained_y = df_all[df_all['gained_y'] == True]
df_gained_y.to_csv(output_gained_y, sep='\t')
