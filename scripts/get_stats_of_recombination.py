import pandas as pd
import sys
import os

# Usage: get_stats_of_recombination.py <base_name> <output_dir> <strain_number>
base_name     = sys.argv[1]
output_dir    = sys.argv[2]
strain_id     = sys.argv[3]

print("Starting get_stats_of_recombination.py")
print(f'Opening {base_name}...')

post_y_prime_tsv  = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')
good_end_rm_tsv   = os.path.join(output_dir, f'{base_name}_good_end_y_repeatmasker.tsv')
output_recomb_tsv = os.path.join(output_dir, f'{base_name}_y_prime_recombination.tsv')

df_y_repeatmasker = pd.read_csv(good_end_rm_tsv, sep='\t')
df_input          = pd.read_csv(post_y_prime_tsv, sep='\t')

df_input = df_input.dropna(subset=['repeat_length'])
df_input = df_input[df_input['Adapter_After_Telomere'] == True]

# ---- Y' reference order per strain ----
if strain_id == '6991':
    y_prime_order_dict = {
        'chr1L': None,  'chr1R': None,
        'chr2L': ('ID4',), 'chr2R': None,
        'chr3L': None,  'chr3R': None,
        'chr4L': None,  'chr4R': ('ID2','ID2','ID2','ID2','ID2','ID2','ID2'),
        'chr5L': ('ID6',), 'chr5R': ('ID1',),
        'chr6L': ('ID4',), 'chr6R': None,
        'chr7L': None,  'chr7R': ('ID5',),
        'chr8L': ('ID1',), 'chr8R': ('ID1',),
        'chr9L': ('ID6',), 'chr9R': None,
        'chr10L': ('ID6',), 'chr10R': None,
        'chr11L': None, 'chr11R': None,
        'chr12L': ('ID1',), 'chr12R': ('ID1','ID2','ID2','ID2','ID2','ID2','ID2'),
        'chr13L': ('ID1',), 'chr13R': None,
        'chr14L': ('ID5','ID2','ID3','ID3','ID3'), 'chr14R': ('ID6',),
        'chr15L': None, 'chr15R': ('ID2',),
        'chr16L': ('ID5',), 'chr16R': ('ID1',)
    }
elif strain_id == '7172':
    y_prime_order_dict = {
        'chr1L': None,  'chr1R': None,
        'chr2L': ('ID4',), 'chr2R': None,
        'chr3L': None,  'chr3R': None,
        'chr4L': None,  'chr4R': ('ID2','ID2','ID2','ID2','ID2','ID2','ID2'),
        'chr5L': ('ID6',), 'chr5R': ('ID1',),
        'chr6L': ('ID4',), 'chr6R': None,
        'chr7L': None,  'chr7R': ('ID5',),
        'chr8L': ('ID1',), 'chr8R': ('ID1',),
        'chr9L': ('ID6',), 'chr9R': None,
        'chr10L': ('ID6',), 'chr10R': None,
        'chr11L': None, 'chr11R': None,
        'chr12L': ('ID1',), 'chr12R': ('ID1','ID2','ID2','ID2','ID2','ID2','ID2','ID2'),
        'chr13L': ('ID1',), 'chr13R': None,
        'chr14L': ('ID5','ID2','ID3'), 'chr14R': ('ID6',),
        'chr15L': None, 'chr15R': ('ID2',),
        'chr16L': ('ID5','ID2','ID2'), 'chr16R': ('ID1',)
    }
elif strain_id == '7302':
    y_prime_order_dict = {
        'chr1L': None,  'chr1R': None,
        'chr2L': ('ID4',), 'chr2R': None,
        'chr3L': None,  'chr3R': None,
        'chr4L': None,  'chr4R': ('ID2','ID2','ID2','ID2','ID2','ID2','ID2'),
        'chr5L': ('ID6',), 'chr5R': ('ID1',),
        'chr6L': ('ID4',), 'chr6R': None,
        'chr7L': None,  'chr7R': ('ID5',),
        'chr8L': ('ID1',), 'chr8R': ('ID1',),
        'chr9L': ('ID6',), 'chr9R': None,
        'chr10L': ('ID6',), 'chr10R': None,
        'chr11L': None, 'chr11R': None,
        'chr12L': ('ID1',), 'chr12R': ('ID1','ID2','ID2','ID2','ID2','ID2'),
        'chr13L': ('ID7','ID2','ID7','ID2'), 'chr13R': None,
        'chr14L': ('ID5','ID2','ID3','ID3','ID3'), 'chr14R': ('ID6',),
        'chr15L': None, 'chr15R': ('ID2',),
        'chr16L': ('ID5',), 'chr16R': ('ID1',)
    }
else:
    raise ValueError(f'Unknown strain_id: "{strain_id}" (derived from base_name: "{base_name}"). '
                     f'Add this strain to the y_prime_order_dict section.')

df_y_repeatmasker['y_prime_id_and_color'] = df_y_repeatmasker['y_prime_group'].apply(lambda x: x.split('/')[2])
df_y_repeatmasker['y_prime_id']           = df_y_repeatmasker['y_prime_id_and_color'].apply(lambda x: x.split('_')[0])
df_y_repeatmasker['y_prime_id_and_match'] = (df_y_repeatmasker['y_prime_id'].astype(str)
                                               + '-' + df_y_repeatmasker['match_id'].astype(str))


def calculate_y_prime_change(read_id, chr_end):
    y_prime_ref_order   = y_prime_order_dict[chr_end]
    recombination_events = {}

    if read_id not in df_y_repeatmasker['read_id'].values:
        y_prime_match_order             = []
        better_repeatmasker_y_prime_count = 0
    else:
        df_read = df_y_repeatmasker[df_y_repeatmasker['read_id'] == read_id].reset_index(drop=True)
        telomere_side = df_read.at[0, 'telomere_side']
        print(df_read)

        y_prime_match_order             = df_read['y_prime_id_and_match'].unique()
        better_repeatmasker_y_prime_count = len(y_prime_match_order)
        print(f'df_read {df_read}')

        if telomere_side == 'beginning':
            y_prime_match_order = y_prime_match_order[::-1]
        print(f'y_prime_match_order {y_prime_match_order}')

    for i, y_prime_and_match_id in enumerate(y_prime_match_order):
        current_y_in_read = y_prime_and_match_id.split('-')[0]

        if y_prime_ref_order is None or i >= len(y_prime_ref_order):
            expected_y_in_ref = None
        else:
            expected_y_in_ref = y_prime_ref_order[i]

        print(y_prime_ref_order)
        print(f'from {expected_y_in_ref} to {current_y_in_read}')

        if current_y_in_read != expected_y_in_ref:
            recombination_events[f'read_y_num_{i+1}'] = f'from {expected_y_in_ref} to {current_y_in_read}'

    if y_prime_ref_order is not None:
        if len(y_prime_match_order) < len(y_prime_ref_order):
            for i in range(len(y_prime_match_order), len(y_prime_ref_order)):
                recombination_events[f'ref_y_num_{i+1}'] = f'from {y_prime_ref_order[i]} to {None}'

    if len(recombination_events) == 0:
        recombination_status = 'No Change'
    elif 'read_y_num_1' in recombination_events or 'ref_y_num_1' in recombination_events:
        recombination_status = "1st Y' Change"
    else:
        recombination_status = "Y' Recombination"

    return recombination_events, recombination_status, better_repeatmasker_y_prime_count


all_reads = df_input['read_id'].unique()

y_prime_recomb_events_dict           = {}
y_prime_recomb_status_dict           = {}
better_repeatmasker_y_prime_count_dict = {}

for read_id in all_reads:
    chr_end = f"chr{df_input[df_input['read_id'] == read_id].iloc[0]['chr_end']}"
    recomb_events, recomb_status, rm_count = calculate_y_prime_change(read_id, chr_end)
    y_prime_recomb_events_dict[read_id]            = recomb_events
    y_prime_recomb_status_dict[read_id]            = recomb_status
    better_repeatmasker_y_prime_count_dict[read_id] = rm_count

df_input['y_prime_recombination_events'] = df_input['read_id'].apply(lambda x: y_prime_recomb_events_dict[x])
df_input['y_prime_recombination_status'] = df_input['read_id'].apply(lambda x: y_prime_recomb_status_dict[x])
df_input['better_repeatmasker_y_prime_count'] = df_input['read_id'].apply(lambda x: better_repeatmasker_y_prime_count_dict[x])
df_input['better_repeatmasker_y_primes_relative_to_ref'] = (
    df_input['better_repeatmasker_y_prime_count'] - df_input['reference_y_primes'])

def calc_delta_group(delta):
    if delta == 0:   return 'Same Number'
    elif delta < 0:  return 'Loss'
    elif delta == 1: return 'Gain 1'
    else:            return 'Gain Multiple'

df_input['y_prime_delta_group'] = df_input['better_repeatmasker_y_primes_relative_to_ref'].apply(calc_delta_group)

for chr_end, group in df_input.groupby('chr_end'):
    print(f'Chromosome End: {chr_end}, with {len(group)} reads')
    print(group['y_prime_recombination_status'].value_counts())
    pct = (group['y_prime_recombination_status'].value_counts(normalize=True) * 100).apply(lambda x: f'{x:.2f}%')
    print(pct)
    print('\n')

print(f'Total Reads: {len(df_input)}\n')
df_input.to_csv(output_recomb_tsv, sep='\t', index=False)
print(f'Finished {base_name}.')
