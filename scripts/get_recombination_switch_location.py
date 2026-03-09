import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# Usage: get_recombination_switch_location.py <base_name> <output_dir>
base_name  = sys.argv[1]
output_dir = sys.argv[2]

print("Starting get_recombination_switch_location.py")

spacer_output_dir = os.path.join(output_dir, 'graphs', 'spacer_outputs')
os.makedirs(spacer_output_dir, exist_ok=True)

spacer_tsv  = os.path.join(output_dir, f'{base_name}_paired_good_spacer_repeatmasker.tsv')
x_ends_tsv  = os.path.join(output_dir, f'{base_name}_good_x_element_ends_paired_repeatmasker.tsv')
y_prime_tsv         = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')
combined_masker_tsv = os.path.join(output_dir, f'{base_name}_combined_repeatmasker.tsv')


def calculate_index_of_switch(chr_end_matching_order_list, anchored_chr_end):
    conservate_switch_index = 0
    aggressive_switch_index = 0
    score         = 0
    max_score     = 0
    first_switch  = True
    switch_to_chr_end = anchored_chr_end

    for i, matching_chr_end in enumerate(chr_end_matching_order_list):
        if matching_chr_end != anchored_chr_end:
            if i == 0:
                switch_to_chr_end = matching_chr_end
            score += 1
            if score >= max_score:
                aggressive_switch_index = i
        else:
            if i == 0:
                conservate_switch_index = 0
                aggressive_switch_index = 0
                switch_to_chr_end = anchored_chr_end
                break
            else:
                score -= 1
                if first_switch:
                    first_switch = False
                    conservate_switch_index = i

    return switch_to_chr_end, conservate_switch_index, aggressive_switch_index


print(f'Opening {spacer_tsv}...')
df_paired_spacer_data = pd.read_csv(spacer_tsv,  sep='\t')
df_paired_x_ends_data = pd.read_csv(x_ends_tsv,  sep='\t')
df_y_prime_data       = pd.read_csv(y_prime_tsv,         sep='\t')  # reference_y_primes
df_combined_masker    = pd.read_csv(combined_masker_tsv,  sep='\t')  # y_prime_id

# Remove known problematic ends
for df in [df_paired_spacer_data, df_paired_x_ends_data]:
    df.drop(df[df['original_chr_end_anchor'].isin(['chr1L', 'chr1R'])].index, inplace=True)
    df.reset_index(drop=True, inplace=True)

df_y_prime_data = df_y_prime_data[~df_y_prime_data['chr_end'].isin(['1L', '1R'])].reset_index(drop=True)

reads_of_y_prime_switches        = df_paired_spacer_data['read_id'].unique()
reads_of_y_prime_switches_x_ends = df_paired_x_ends_data['read_id'].unique()

print(f'Reads with spacer switches: {len(reads_of_y_prime_switches)}')

conservate_switch_distance_list    = []
aggressive_switch_distance_list    = []
spacer_switching_chr_end_pair_list = []
has_x_end_list                     = []
has_its_in_donor_list              = []

ends_with_its = ['chr4R', 'chr5L', 'chr5R', 'chr6L', 'chr8L', 'chr12R', 'chr13L', 'chr14L', 'chr14R', 'chr16L']

for read_id in reads_of_y_prime_switches:
    try:
        df_read = df_paired_spacer_data[df_paired_spacer_data['read_id'] == read_id].reset_index(drop=True)

        anchored_chr_end = df_read['original_chr_end_anchor'].iloc[0]
        strand           = df_read[df_read['anchor_label'] == 'Is Anchor']['strand'].mode().iloc[0]
        df_read['matched_chr_end'] = df_read['chr_end_tract'].apply(lambda x: x.split('_')[0])
        chr_end_order = df_read['matched_chr_end'].to_list()

        def get_switch_distance_fwd(df_r, switch_idx):
            if switch_idx == 0:
                return 0
            mn = df_r['match_end_on_read'][switch_idx]
            mx = df_r['match_start_on_read'][switch_idx + 1]
            return (mx + mn) / 2 - df_r['match_start_on_read'][0]

        def get_switch_distance_rev(df_r, switch_idx):
            if switch_idx == 0:
                return 0
            L  = len(chr_end_order) - 1
            mn = df_r['match_start_on_read'][L - switch_idx]
            mx = df_r['match_end_on_read'][L - (switch_idx + 1)]
            return df_r['match_start_on_read'][L] - (mx + mn) / 2

        is_left = 'L' in anchored_chr_end
        is_plus = strand == '+'

        forward = (is_left and is_plus) or (not is_left and not is_plus)

        if not forward:
            chr_end_order = chr_end_order[::-1]

        switch_to, c_idx, a_idx = calculate_index_of_switch(chr_end_order, anchored_chr_end)

        if forward:
            conservate_switch_distance = get_switch_distance_fwd(df_read, c_idx)
            aggressive_switch_distance = get_switch_distance_fwd(df_read, a_idx)
        else:
            conservate_switch_distance = get_switch_distance_rev(df_read, c_idx)
            aggressive_switch_distance = get_switch_distance_rev(df_read, a_idx)

        spacer_switching_chr_end_pair = ('No Switch' if anchored_chr_end == switch_to
                                         else f'Anchored {anchored_chr_end} switch to {switch_to}')

        if read_id in reads_of_y_prime_switches_x_ends:
            df_x  = df_paired_x_ends_data[df_paired_x_ends_data['read_id'] == read_id].reset_index(drop=True)
            x_end = df_x['x_element_ends'].iloc[0].split('_')[0]
            x_element_read_outcome = "No Switch in X element" if x_end == anchored_chr_end else "Switch in X element"
        else:
            x_element_read_outcome = 'No Results'

        try:
            donor_chr_ends = df_combined_masker[df_combined_masker['read_id'] == read_id]['y_prime_id'].iloc[0]
            its_in_donor   = any(end in donor_chr_ends for end in ends_with_its)
        except IndexError:
            its_in_donor = 'No Result'

        conservate_switch_distance_list.append(conservate_switch_distance)
        aggressive_switch_distance_list.append(aggressive_switch_distance)
        spacer_switching_chr_end_pair_list.append(spacer_switching_chr_end_pair)
        has_x_end_list.append(x_element_read_outcome)
        has_its_in_donor_list.append(its_in_donor)

    except Exception as e:
        print(f'Error for read {read_id}: {e}')
        conservate_switch_distance_list.append(None)
        aggressive_switch_distance_list.append(None)
        spacer_switching_chr_end_pair_list.append(None)
        has_x_end_list.append(None)
        has_its_in_donor_list.append(None)

df_switches = pd.DataFrame({
    'base_name':                                base_name,
    'conservate_switch_distance':               conservate_switch_distance_list,
    'aggressive_switch_distance':               aggressive_switch_distance_list,
    'anchor_and_spacer_switch_to_chr_end_pair': spacer_switching_chr_end_pair_list,
    'has_x_end':                                has_x_end_list,
    'has_its_in_donor':                         has_its_in_donor_list,
})

output_tsv = os.path.join(output_dir, f'{base_name}_recombination_switch_location.tsv')
df_switches.to_csv(output_tsv, sep='\t', index=False)
print(f'Saved switch location results to {output_tsv}')

# ---- Summary stats ----
print(f'\n--- Summary for {base_name} ---')
print(df_switches['anchor_and_spacer_switch_to_chr_end_pair'].value_counts())
print()
print(df_switches['has_x_end'].value_counts())
print()

total = len(df_switches)
if total > 0:
    try:
        no_switch_spacer = df_switches['anchor_and_spacer_switch_to_chr_end_pair'].value_counts()['No Switch']
    except KeyError:
        no_switch_spacer = 0
    switch_spacer = total - no_switch_spacer
    print(f"Total 1st Y' switch events: {total}")
    print(f'Switch in Spacer:    {switch_spacer}  ({switch_spacer / total * 100:.1f}%)')
    print(f'No Switch in Spacer: {no_switch_spacer}  ({no_switch_spacer / total * 100:.1f}%)')

print(f'Finished {base_name}.')