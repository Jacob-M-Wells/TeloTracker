import pandas as pd
import os
import sys
import math

# Usage: make_spacer_pairs_repeatmasker_tsv.py <base_name> <output_dir> <strain_number>
base_name  = sys.argv[1]
output_dir = sys.argv[2]
strain_id  = sys.argv[3]

print("Starting make_spacer_pairs_repeatmasker_tsv.py")
print(f'Opening {base_name}...')

# Strain reference files live under a known references directory
# The shell script passes STRAIN_REF_DIR; here we derive from strain_id via output_dir sibling
strain_features_dir  = os.path.join(os.path.dirname(output_dir), 'references', f'{strain_id}_features')
bed_file_of_features = os.path.join(strain_features_dir, f'{strain_id}_final_features.bed')

repeatmasker_dir = os.path.join(output_dir, 'repeatmasker', 'spacers')

output_all    = os.path.join(output_dir, f'{base_name}_paired_spacer_repeatmasker.tsv')
output_good   = os.path.join(output_dir, f'{base_name}_paired_good_spacer_repeatmasker.tsv')
output_gained = os.path.join(output_dir, f'{base_name}_paired_good_gained_spacer_repeatmasker.tsv')

post_y_prime_tsv = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')

all_results_files = [
    os.path.join(repeatmasker_dir, f)
    for f in os.listdir(repeatmasker_dir)
    if f.endswith('_spacer_repeatmasker.ssv')
]

df_all = pd.DataFrame()
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
        chr_end = os.path.basename(corrected)
        chr_end = chr_end.replace('_spacer_repeatmasker_results_corrected.ssv', '')
        chr_end = chr_end.replace(f'{base_name}_', '')
        chr_end = chr_end.split('_and_ID')[0]
        df_single['original_chr_end_anchor'] = chr_end
        df_single = df_single.dropna(axis=1, how='all')
        df_all = pd.concat([df_all, df_single])
    except Exception:
        print(f'Error in {corrected}: {sys.exc_info()}')

df_all.loc[df_all['sub_match'] == '*', 'sub_match'] = True
df_all.loc[df_all['sub_match'] == '-', 'sub_match'] = False

print(df_all)

# ---- Label anchor sections using features BED ----
bed_columns = ['Chr', 'Start', 'End', 'Name', 'Strand', 'Length']
df_bed = pd.read_csv(bed_file_of_features, sep='\t', names=bed_columns)

distances_to_anchor_dict = {}
for chr_end in df_all['original_chr_end_anchor'].unique():
    chopped_segment_length = 250
    space_rows   = df_bed[df_bed['Name'] == f'{chr_end}_space_between_anchor']
    anchor_rows  = df_bed[df_bed['Name'] == f'{chr_end}_anchor']
    if space_rows.empty or anchor_rows.empty:
        continue
    dist_start = space_rows['Length'].iloc[0]
    dist_end   = dist_start + anchor_rows['Length'].iloc[0]
    distances_to_anchor_dict[chr_end] = (dist_start / chopped_segment_length,
                                          dist_end   / chopped_segment_length)

def label_anchor_on_chopped_segments(chr_end, section_number, row):
    if chr_end not in distances_to_anchor_dict:
        return 'Unknown'
    section_of_anchor_start, section_of_anchor_end = distances_to_anchor_dict[chr_end]

    if section_number < math.floor(section_of_anchor_start):
        return 'Before Anchor'
    elif section_number == math.floor(section_of_anchor_start):
        return 'Partial Before Anchor' if section_number == section_of_anchor_start else 'Is Anchor'
    elif section_number <= section_of_anchor_end:
        return 'Is Anchor'
    elif section_number < math.ceil(section_of_anchor_end):
        return 'Partial After Anchor'
    else:
        return 'After Anchor'

df_all['anchor_label'] = df_all.apply(
    lambda r: label_anchor_on_chopped_segments(r['original_chr_end_anchor'], r['section_number'], r), axis=1)

df_all.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
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
