#!/bin/bash
set -euo pipefail

#######################################################################################
# TeloTracker - Analysis Pipeline
#
# Usage:
#   bash telomere_analysis.sh <input> <base_name> <output_dir> <strain_number> [strain_ref_dir]
#
# Arguments:
#   input            Path to input .bam, .fastq, or .fastq.gz file
#   base_name        Sample name used for output file naming (e.g. "6991_day0")
#   output_dir       Directory where all outputs will be written
#   strain_number    Strain identifier (e.g. 6991, 7172, 7302)
#   strain_ref_dir   (Optional) Path to strain-specific reference directory for Steps 7-9
#                    e.g. "references/6991_features" — omit or set to "none" for core analysis only
#
# Examples:
#   bash telomere_analysis.sh sample.bam 6991_day0 ./results 6991
#   bash telomere_analysis.sh sample.fastq.gz 6991_day0 ./results 6991
#   bash telomere_analysis.sh sample.fastq.gz 6991_day0 ./results 6991 references/6991_features
#
# Repository structure expected:
#   telotracker/
#   ├── telomere_analysis.sh
#   ├── scripts/        # all Python analysis scripts
#   └── references/
#       ├─── universal/     # general adapter and Y' sequence files
#       │   ├── nanopore_sqk-slk114_adapter_sequence_truncated.txt
#       │   └── y_prime_probe.fasta
#       ├── anchors/        # anchor sequence files
#       │   ├── telomerase_shutoff_anchors.fasta
#       │   └── telomerase_deletion_anchors.fasta
#       ├── 6991_features/      # strain-specific reference feature files
#       │   ├── 6991_final_features.bed
#       │   ├── repeatmasker_6991_all_y_primes.fasta
#       │   ├── 6991_x_element_ends_pairings/
#       │   └── 6991_spacer_pairings/
#       ├── 7172_features/      # strain-specific reference feature files
#       │   ├── 7172_final_features.bed
#       │   ├── repeatmasker_7172_all_y_primes.fasta
#       │   ├── 7172_x_element_ends_pairings/
#       │   └── 7172_spacer_pairings/
#       └── 7302_features/      # strain-specific reference feature files
#           ├── 7302_final_features.bed
#           ├── repeatmasker_7302_all_y_primes.fasta
#           ├── 7302_x_element_ends_pairings/
#           └── 7302_spacer_pairings/
#######################################################################################


# --- Arguments ---
INPUT=$1
BASE_NAME=$2
OUTPUT_DIR=$3
STRAIN_NUMBER=$4
STRAIN_REF_DIR=${5:-"none"}

# --- Configuration ---
THREADS=${THREADS:-80}

# Anchor set — choose one:
#   telomerase_shutoff_anchors  : anchors designed for telomerase shutoff samples (6991, and subsequent transformants)
#   all_telomerase_deletion_anchorsanchors  : anchors for telomerase deletion samples (6212)
ANCHOR_SET="telomerase_shutoff_anchors"

# --- Reference paths ---
SCRIPT_DIR="$(dirname "$0")/scripts"

UNIVERSAL_REF_DIR="$(dirname "$0")/references/universal"
ADAPTER_SEQ_FILE=$UNIVERSAL_REF_DIR/nanopore_sqk-slk114_adapter_sequence_truncated.txt
Y_PRIME_PROBE_DB=$UNIVERSAL_REF_DIR/y_prime_probe.fasta

ANCHOR_REF_DIR="$(dirname "$0")/references/anchors"
ANCHOR_DB=$ANCHOR_REF_DIR/${ANCHOR_SET}.fasta

# --- Output subdirectories ---
dir_blast=$OUTPUT_DIR/blast
dir_chr_anchors=$OUTPUT_DIR/blast/chr_anchor_reads
dir_y_prime_blast=$OUTPUT_DIR/blast/y_prime_probe
dir_porechop=$OUTPUT_DIR/porechop
dir_graphs=$OUTPUT_DIR/graphs
dir_repeatmasker_y=$OUTPUT_DIR/repeatmasker/y_primes
dir_repeatmasker_x=$OUTPUT_DIR/repeatmasker/x_element_ends
dir_repeatmasker_spacer=$OUTPUT_DIR/repeatmasker/spacers

porechop_log=$dir_porechop/${BASE_NAME}_porechop.log

#######################################################################################

echo "===== TeloTracker: Starting $BASE_NAME ====="
echo "  Input:         $INPUT"
echo "  Base name:     $BASE_NAME"
echo "  Output dir:    $OUTPUT_DIR"
echo "  Strain number: $STRAIN_NUMBER"
echo "  Anchor set:    $ANCHOR_SET"
echo "  Strain ref:    $STRAIN_REF_DIR"
echo ""

mkdir -p $dir_blast $dir_chr_anchors $dir_y_prime_blast $dir_porechop $dir_graphs

# -------------------------------------------------------------------------------------
echo "----- Step 0: Prep input files -----"
# -------------------------------------------------------------------------------------

raw_fastq=$OUTPUT_DIR/${BASE_NAME}_raw.fastq
filtered_fastq=$OUTPUT_DIR/${BASE_NAME}.fastq
input_fasta=$OUTPUT_DIR/${BASE_NAME}.fasta

# Step 0a: Convert input to raw fastq
if [[ "$INPUT" == *.bam ]]; then
    echo "  Converting BAM to FASTQ..."
    samtools fastq -T '*' $INPUT > $raw_fastq
elif [[ "$INPUT" == *.fastq.gz ]]; then
    echo "  Decompressing FASTQ..."
    zcat $INPUT > $raw_fastq
elif [[ "$INPUT" == *.fastq ]]; then
    echo "  Copying FASTQ..."
    cp $INPUT $raw_fastq
else
    echo "ERROR: Unrecognized input format: $INPUT"
    echo "       Expected .bam, .fastq, or .fastq.gz"
    exit 1
fi

# Step 0b: Filter reads by length (>= 2000bp) and qscore (qs >= 10)
echo "  Filtering reads (min length: 2000bp, min qscore: 10)..."
python $SCRIPT_DIR/filter_reads.py $raw_fastq $filtered_fastq

# Step 0c: Convert filtered fastq to fasta
# sed explanation:
#   1~4s/^@/>/  — on every 4th line starting at line 1 (the header), replace @ with >
#   2~4p        — print every 4th line starting at line 2 (the sequence)
sed -n '1~4s/^@/>/p;2~4p' $filtered_fastq > $input_fasta

# -------------------------------------------------------------------------------------
echo "----- Step 1: BLAST reads for chromosome anchors -----"
# -------------------------------------------------------------------------------------

blast_out=$dir_blast/${BASE_NAME}_blasted_${ANCHOR_SET}.tsv

echo -e "read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\tanchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue" \
    > $blast_out

blastn \
    -query $input_fasta \
    -db $ANCHOR_DB \
    -task dc-megablast \
    -perc_identity 85 \
    -min_raw_gapped_score 5000 \
    -num_threads $THREADS \
    -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" \
    >> $blast_out

# -------------------------------------------------------------------------------------
echo "----- Step 2: Filter BLAST results for anchored reads -----"
# -------------------------------------------------------------------------------------

python $SCRIPT_DIR/filter_for_reads_with_anchors.py $BASE_NAME $ANCHOR_SET $OUTPUT_DIR

# -------------------------------------------------------------------------------------
echo "----- Step 3: Split and label chromosome-anchored reads -----"
# -------------------------------------------------------------------------------------

python $SCRIPT_DIR/split_and_label_all_reads_include_anchor.py $BASE_NAME $ANCHOR_SET $OUTPUT_DIR

# -------------------------------------------------------------------------------------
echo "----- Step 4: Trim adapters with Porechop -----"
# -------------------------------------------------------------------------------------

cat $dir_chr_anchors/${BASE_NAME}_blasted_${ANCHOR_SET}_chr*_anchor_reads.fasta \
    > $OUTPUT_DIR/${BASE_NAME}_all_chr_anchored_reads.fasta

porechop_abi \
    -t $THREADS \
    --format fastq -ddb -v 3 --no_split \
    -cap $ADAPTER_SEQ_FILE \
    -i $OUTPUT_DIR/${BASE_NAME}_all_chr_anchored_reads.fasta \
    -o $dir_porechop/${BASE_NAME}_trimmed.fastq \
    > $porechop_log

sed -n '1~4s/^@/>/p;2~4p' $dir_porechop/${BASE_NAME}_trimmed.fastq \
    > $dir_porechop/${BASE_NAME}_trimmed.fasta

python $SCRIPT_DIR/check_for_adapters.py $BASE_NAME $porechop_log $OUTPUT_DIR
python $SCRIPT_DIR/compare_adapter_callers_dorado.py $BASE_NAME $ANCHOR_SET $OUTPUT_DIR $input_fasta
python $SCRIPT_DIR/fine_telomere_trimming.py $BASE_NAME $ANCHOR_SET $OUTPUT_DIR

# -------------------------------------------------------------------------------------
echo "----- Step 5: BLAST anchored reads for Y' probe sequences -----"
# -------------------------------------------------------------------------------------

for chr_num in {1..16}; do
    for side in L R; do
        echo "  Blasting chr${chr_num}${side}"

        query=$dir_chr_anchors/${BASE_NAME}_blasted_${ANCHOR_SET}_chr${chr_num}${side}_anchor_reads
        out=$dir_y_prime_blast/${BASE_NAME}_chr${chr_num}${side}_blasted_probe.tsv

        echo -e "read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\tanchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue" \
            > $out

        blastn \
            -query $query.fasta \
            -db $Y_PRIME_PROBE_DB \
            -perc_identity 90 \
            -num_threads $THREADS \
            -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" \
            >> $out
    done
done

# -------------------------------------------------------------------------------------
echo "----- Step 6: Perfrom Y' analysis -----"
# -------------------------------------------------------------------------------------

python $SCRIPT_DIR/y_prime_analysis.py $BASE_NAME $OUTPUT_DIR

# -------------------------------------------------------------------------------------
echo "----- Step 6: Create telomere graphs -----"
# -------------------------------------------------------------------------------------

python $SCRIPT_DIR/single_sample_plots.py $BASE_NAME $OUTPUT_DIR

# -------------------------------------------------------------------------------------
# Steps 8-10: Reference-based RepeatMasker analysis (Y', X element ends, spacers)
# Requires a strain-specific reference directory — pass as 4th argument
# -------------------------------------------------------------------------------------

if [[ "$STRAIN_REF_DIR" == "none" ]]; then
    echo ""
    echo "===== No strain reference directory provided — skipping Steps 7-9. ====="
    echo "===== TeloTracker: $BASE_NAME complete (core analysis only). ====="
    exit 0
fi

if [[ ! -d "$STRAIN_REF_DIR" ]]; then
    echo ""
    echo "WARNING: Strain reference directory '$STRAIN_REF_DIR' not found — skipping Steps 7-9."
    echo "===== TeloTracker: $BASE_NAME complete (core analysis only). ====="
    exit 0
fi

echo "  Using reference strain feature for: $STRAIN_NUMBER ($STRAIN_REF_DIR)"

# -------------------------------------------------------------------------------------
echo "----- Step 8: RepeatMasker for Y' sequences -----"
# -------------------------------------------------------------------------------------

mkdir -p "$dir_repeatmasker_y"

if [ -n "$(ls -A "$dir_repeatmasker_y")" ]; then
    rm -r "$dir_repeatmasker_y"/*
fi

for chr_num in {1..16}; do
    for side in L R; do
        echo "  Searching chr${chr_num}${side}"

        query=$dir_chr_anchors/${BASE_NAME}_blasted_${ANCHOR_SET}_chr${chr_num}${side}_anchor_reads
        out=$dir_repeatmasker_y/${BASE_NAME}_chr${chr_num}${side}_y_prime_repeatmasker.ssv

        RepeatMasker \
            $query.fasta \
            -lib $STRAIN_REF_DIR/repeatmasker_${STRAIN_NUMBER}_all_y_primes.fasta \
            -s -pa 14 --cutoff 1000 -no_is -norna -gff \
            -dir $dir_repeatmasker_y

        echo -e "SW_score\tdivergence_percent\tdeletion_percent\tinsertion_percent\tread_id\tmatch_start_on_read\tmatch_end_on_read\tleftover_on_read\tstrand\ty_prime_id\ty_prime_group\tmatch_start_on_y_prime\tmatch_end_on_y_prime\tleftover_on_y_prime\tmatch_id\tsub_match" \
            > $out

        grep "Y_Prime" $dir_repeatmasker_y/${BASE_NAME}_blasted_${ANCHOR_SET}_chr${chr_num}${side}_anchor_reads.fasta.out \
            >> $out || true    # prevents grep from exiting with 'done'
    done
done

python $SCRIPT_DIR/make_y_prime_repeatmasker_tsv.py $BASE_NAME $OUTPUT_DIR
python $SCRIPT_DIR/get_summary_stats_for_y_prime_repeatmasker.py $BASE_NAME $ANCHOR_SET $OUTPUT_DIR
python $SCRIPT_DIR/get_stats_of_recombination.py $BASE_NAME $OUTPUT_DIR $STRAIN_NUMBER
python $SCRIPT_DIR/make_pairings_from_y_primes_all_ends.py $BASE_NAME $ANCHOR_SET $OUTPUT_DIR $STRAIN_NUMBER

# -------------------------------------------------------------------------------------
echo "----- Step 9: RepeatMasker for X element ends -----"
# -------------------------------------------------------------------------------------

mkdir -p "$dir_repeatmasker_x"

if [ -n "$(ls -A "$dir_repeatmasker_x")" ]; then
    rm -r "$dir_repeatmasker_x"/*
fi

for pairing_file in $OUTPUT_DIR/paired_by_y_prime_reads/*.fasta; do
    pairing=$(basename "$pairing_file" .fasta)
    echo "  Processing $pairing"

    RepeatMasker \
        $pairing_file \
        -lib $STRAIN_REF_DIR/${STRAIN_NUMBER}_x_element_ends_pairings/${STRAIN_NUMBER}_paired_${pairing}.fasta \
        -s -pa 14 --cutoff 500 -no_is -norna -gff \
        -dir $dir_repeatmasker_x

    filtered=$dir_repeatmasker_x/${pairing}_x_element_ends_filtered.out
    grep "x_ends" $dir_repeatmasker_x/${pairing}.fasta.out > $filtered || true  # prevents grep from exiting with 'done'

    out=$dir_repeatmasker_x/${BASE_NAME}_${pairing}_x_element_ends_repeatmasker.ssv
    echo -e "SW_score\tdivergence_percent\tdeletion_percent\tinsertion_percent\tread_id\tmatch_start_on_read\tmatch_end_on_read\tleftover_on_read\tstrand\tx_element_ends\tsection_number\tmatch_start_on_chr_end_section\tmatch_end_on_chr_end_section\tleftover_on_chr_end_section\tmatch_id\tsub_match" \
        > $out
    cat $filtered >> $out
done

python $SCRIPT_DIR/make_x_element_ends_pairs_repeatmasker_tsv.py $BASE_NAME $OUTPUT_DIR $STRAIN_NUMBER

# -------------------------------------------------------------------------------------
echo "----- Step 10: RepeatMasker for spacer sequences -----"
# -------------------------------------------------------------------------------------

mkdir -p "$dir_repeatmasker_spacer"

if [ -n "$(ls -A "$dir_repeatmasker_spacer")" ]; then
    rm -r "$dir_repeatmasker_spacer"/*
fi

for pairing_file in $OUTPUT_DIR/paired_by_y_prime_reads/*.fasta; do
    pairing=$(basename "$pairing_file" .fasta)
    echo "  Processing $pairing"

    RepeatMasker \
        $pairing_file \
        -lib $STRAIN_REF_DIR/${STRAIN_NUMBER}_spacer_pairings/${STRAIN_NUMBER}_paired_${pairing}.fasta \
        -s -pa 14 --cutoff 500 -no_is -norna -gff \
        -dir $dir_repeatmasker_spacer

    filtered=$dir_repeatmasker_spacer/${pairing}_spacer_filtered.out
    grep "_from_repeat_to_plus_50kb" $dir_repeatmasker_spacer/${pairing}.fasta.out > $filtered || true  # prevents grep from exiting with 'done'
    out=$dir_repeatmasker_spacer/${BASE_NAME}_${pairing}_spacer_repeatmasker.ssv
    echo -e "SW_score\tdivergence_percent\tdeletion_percent\tinsertion_percent\tread_id\tmatch_start_on_read\tmatch_end_on_read\tleftover_on_read\tstrand\tchr_end_tract\tsection_number\tmatch_start_on_chr_end_section\tmatch_end_on_chr_end_section\tleftover_on_chr_end_section\tmatch_id\tsub_match" \
        > $out
    cat $filtered >> $out
done

python $SCRIPT_DIR/make_spacer_pairs_repeatmasker_tsv.py $BASE_NAME $OUTPUT_DIR $STRAIN_NUMBER

# -------------------------------------------------------------------------------------
echo "----- Step 10: Get the recombination switch location -----"
# -------------------------------------------------------------------------------------

python $SCRIPT_DIR/get_recombination_switch_location.py $BASE_NAME $OUTPUT_DIR $STRAIN_NUMBER

#######################################################################################
echo "===== TeloTracker: $BASE_NAME complete ====="