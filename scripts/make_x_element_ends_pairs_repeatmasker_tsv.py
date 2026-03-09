import pandas as pd
import os
import sys

# Usage: make_x_element_ends_pairs_repeatmasker_tsv.py <base_name> <output_dir> <strain_number>
base_name  = sys.argv[1]
output_dir = sys.argv[2]
strain_id  = sys.argv[3]

print("Starting make_x_element_ends_pairs_repeatmasker_tsv.py")
print(f'Opening {base_name}...')

repeatmasker_dir = os.path.join(output_dir, 'repeatmasker', 'x_element_ends')

output_all    = os.path.join(output_dir, f'{base_name}_paired_x_element_ends_repeatmasker.tsv')
output_good   = os.path.join(output_dir, f'{base_name}_good_x_element_ends_paired_repeatmasker.tsv')
output_gained = os.path.join(output_dir, f'{base_name}_good_gained_y_x_element_ends_paired_repeatmasker.tsv')

post_y_prime_tsv = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')

all_results_files = [
    os.path.join(repeatmasker_dir, f)
    for f in os.listdir(repeatmasker_dir)
    if f.endswith('_x_element_ends_repeatmasker.ssv')
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
        # Parse chr_end from filename
        chr_end = os.path.basename(corrected)
        chr_end = chr_end.replace('_x_element_ends_repeatmasker_results_corrected.ssv', '')
        chr_end = chr_end.replace(f'{base_name}_', '')
        chr_end = chr_end.split('_and_ID')[0]
        df_single['original_chr_end_anchor'] = chr_end
        df_single = df_single.dropna(axis=1, how='all')
        df_all = pd.concat([df_all, df_single])
    except Exception:
        print(f'Error in {corrected}: {sys.exc_info()}')

df_all.loc[df_all['sub_match'] == '*', 'sub_match'] = True
df_all.loc[df_all['sub_match'] == '-', 'sub_match'] = False

df_all.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
print(df_all)
df_all.to_csv(output_all, sep='\t')

# ---- Filter for good matches ----
df_good = df_all[df_all['sub_match'] == False]
df_good = df_good[df_good['SW_score'] >= 500]
df_good = df_good[df_good['divergence_percent'] <= 2]

df_filter = pd.read_csv(post_y_prime_tsv, sep='\t')
df_filter = df_filter.dropna(subset=['repeat_length'])
df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

good_reads = set(df_filter['read_id'])
df_good['good_read'] = df_good['read_id'].apply(lambda x: x in good_reads)
print(df_good['good_read'].value_counts())
df_good = df_good[df_good['good_read'] == True]

df_good.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
df_good.to_csv(output_good, sep='\t')

# ---- Y prime gain subset ----
reads_gain = set(df_filter[df_filter['delta_y_prime_sign'] == '+']['read_id'])
df_good['gained_y'] = df_good['read_id'].apply(lambda x: x in reads_gain)
print(df_good['gained_y'].value_counts())

df_gained = df_good[df_good['gained_y'] == True].sort_values(
    by=['original_chr_end_anchor', 'read_id', 'match_start_on_read']
)
df_gained.to_csv(output_gained, sep='\t')
