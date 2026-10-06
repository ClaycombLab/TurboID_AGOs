#!/bin/bash

# Deffine inputs tmp files
QUERY="Y66D12A_19.fasta"
DATABASE="all_species.fasta"
PREFIX="Y66D12A.19"
SRCPATH=/home/isaac/Data/TURBO_DATA/PHYLO/src
THREADS=8
QCOV=30

# BLAST
BLAST_EVALUE=1e-3
MAX_TARGET_SEQS=2000

MIN_SUBJ_COV=40
MIN_LEN_RATIO=0.70
MAX_LEN_RATIO=1.30


# HMM filter
FULL_SEQ_EVALUE=1e-10
BEST_DOM_IEVALUE=1e-5
MIN_DOM_SCORE=30
MIN_HMM_COV=0.35
MAX_BIAS_TO_SCORE=0.30


# outFiles
HITS_TSV=${PREFIX}.blastp.qcov${QCOV}.tsv
HIT_IDS=${PREFIX}.hit_ids.txt
OUT_FASTA=${PREFIX}.hits.qcov${QCOV}.fa
OUT_ALN=${PREFIX}.hits.qcov${QCOV}.aln
OUT_HMM=${PREFIX}.hits.qcov${QCOV}.hmm
HMM_BUILD_LOG=${OUT_HMM}.build.log
HMM_STATS=${PREFIX}.hits.qcov${QCOV}.hmmstat.txt
HMMSEARCH_TBL=${PREFIX}.hits.qcov${QCOV}.hmmsearch.tbl
HMMSEARCH_DOMTBL=${PREFIX}.hits.qcov${QCOV}.hmmsearch.domtbl
HMMSEARCH_TXT=${PREFIX}.hits.qcov${QCOV}.hmmsearch.txt
HMMSEARCH_HIT_IDS=${PREFIX}.hmmsearch.hit_ids.txt
POS_IDS=${PREFIX}.positive_ids.txt
RECOVERY_TSV=${PREFIX}.hits.qcov${QCOV}.recovery.tsv

# Filtered HMM
BEST_HMM_SUMMARY=${PREFIX}.best_hmm_summary.tsv
BEST_SEQ_IDS=${PREFIX}.best_hmm_targets.ids.txt
SPLIT_DOMAIN_IDS=${PREFIX}.split_domain_targets.ids.txt

# Combined set for phylogeny
COMBINED_IDS=${PREFIX}.combined_seed_plus_hmm.ids.txt
COMBINED_FASTA=${PREFIX}.combined_seed_plus_hmm.fa
COMBINED_ALN=${PREFIX}.combined_seed_plus_hmm.aln

# timal
TRIMMED_ALN=${PREFIX}.combined_seed_plus_hmm.trimmed.aln

# modeltest-ng
MODELTEST_LOG=${PREFIX}.modeltest-ng.log
BEST_MODEL_TXT=${PREFIX}.best_model.txt
RAXML_PREFIX=${PREFIX}.raxmlng


mkdir -p $PREFIX

chems() {
    echo "
  ⠀⠀⠀⠀⡀⣯⡭⠀⢟⢿⣷⣄⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⢠⣼⠏⠴⠶⠈⠘⠻⣘⡆⠀⠀⠀⠀⠀⠀
  ⠀⠀⣠⣾⡟⠁⡀⠀⠀⠀⡼⠡⠛⡄⠀⠀⠀⠀⠀
  ⠀⠀⠙⠻⢴⠞⠁⠀⠊⠀⠀⠀⠈⠉⢄⠀⠀⠀⠀
  ⠀⠀⠀⠀⢀⠀⠀⠀⢃⠄⠂⠀⠀⢀⠞⢣⡀⠀⠀
  ⠀⠀⠀⠀⡌⠁⠀⠀⠀⢀⠀⠐⠈⠀⠀⡺⠙⡄⠀
  ⠀⠀⠀⠀⡿⡀⠀⠀⠀⠁⠀⠴⠁⠀⠚⠀⡸⣷⠀
  ⠀⠀⠀⠀⢹⠈⠀⠀⠀⠀⠔⠁⠀⢀⠄⠀⠁⢻⣧
  ⠀⠀⠀⠀⣸⠀⢠⣇⠀⢘⣬⠀⠀⣬⣠⣦⣼⣿⠏
  ⡠⠐⢂⡡⠾⢀⡾⠋⠉⠉⡇⠀⢸⣿⣿⣿⡿⠃⠀
  ⠉⢉⡠⢰⠃⠸⠓⠒⠂⠤⡇⠀⡿⠟⠛⠁⠀⠀⠀
  ⠘⢳⡞⣅⡰⠁⠀⠀⠀⢀⠇⠀⡇⠀⠀⠀⠀⠀⠀
  ⠀⠀⠉⠉⠀⠀⠀⠀⢀⣌⢀⢀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠘⠊⠀
  mwaiting"
}

msg() {
    echo
    echo "=== $1 ==="
}

die() {
    echo "ERROR: $1" >&2
    exit 1
}
# check tools
check_tools() {
    for tool in blastp seqkit mafft hmmbuild hmmstat hmmsearch trimal modeltest-ng raxml-ng; do
        command -v "$tool" >/dev/null 2>&1 || die "required tool not found in PATH: $tool"
    done
}

check_inputs() {
    [[ -s "$QUERY" ]] || die "q file not found or empty: $QUERY"
    [[ -s "$DATABASE" ]] || die "Db files not found or empty: $DATABASE"
}


blast_search() {
    msg "Step 1: BLASTP search and qcov/length-aware filter"

    blastp \
      -query $QUERY \
      -db $DATABASE \
      -evalue $BLAST_EVALUE \
      -matrix BLOSUM45 \
      -comp_based_stats 2 \
      -seg yes \
      -max_target_seqs $MAX_TARGET_SEQS \
      -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen stitle' \
      -num_threads $THREADS \
    | awk -v qcov=$QCOV \
          -v min_scov=$MIN_SUBJ_COV \
          -v min_ratio=$MIN_LEN_RATIO \
          -v max_ratio=$MAX_LEN_RATIO '
        BEGIN { OFS="\t" }
        {
          aln_len   = $4
          query_cov = $13
          qlen      = $14
          slen      = $15

          subj_cov  = (slen > 0 ? 100 * aln_len / slen : 0)
          len_ratio = (qlen > 0 ? slen / qlen : 0)

          if (query_cov >= qcov &&
              subj_cov  >= min_scov &&
              len_ratio >= min_ratio &&
              len_ratio <= max_ratio) {
              print
          }
        }
      ' > $HITS_TSV

    cut -f2 $HITS_TSV | sort -u > $HIT_IDS

    if [[ ! -s "$HIT_IDS" ]]; then
        echo "No BLAST hits passed qcov >= ${QCOV}, subject coverage >= ${MIN_SUBJ_COV}, length ratio ${MIN_LEN_RATIO}-${MAX_LEN_RATIO}"
        exit 0
    fi
}

extract_blast_hits() {
    msg "Step 2: Extract BLAST hit sequences"

    seqkit grep -f $HIT_IDS $DATABASE > $OUT_FASTA

    local nseq
    nseq=$(grep -c '^>' "$OUT_FASTA" || true)
    echo "Sequences extracted from BLAST hits: $nseq"

    if [[ "$nseq" -lt 5 ]]; then
        echo "Warning: very few sequences for HMM building ($nseq)"
    fi
}


build_seed_alignment() {
    msg "Step 3: Build MAFFT alignment from BLAST hits"
    mafft --thread $THREADS $OUT_FASTA > $OUT_ALN
}


build_hmm() {
    msg "Step 4: Build HMM"
    hmmbuild $OUT_HMM $OUT_ALN > $HMM_BUILD_LOG
}


summarize_hmm() {
    msg "Step 5: Summarize HMM"
    hmmstat $OUT_HMM > $HMM_STATS
}

search_hmm() {
    msg "Step 6: Search HMM back against database"

    hmmsearch \
      --cpu $THREADS \
      --tblout $HMMSEARCH_TBL \
      --domtblout $HMMSEARCH_DOMTBL \
      $OUT_HMM $DATABASE > $HMMSEARCH_TXT
}

recovery_assessment() {
    msg "Step 7: Simple recovery assessment"

    cp $HIT_IDS $POS_IDS
    awk '!/^#/{print $1}' $HMMSEARCH_TBL | sort -u > $HMMSEARCH_HIT_IDS

    {
      echo -e "metric\tvalue"
      echo -e "training_sequences\t$(wc -l < "$POS_IDS")"
      echo -e "hmmsearch_hits_total\t$(wc -l < "$HMMSEARCH_HIT_IDS")"
      echo -e "training_sequences_recovered\t$(comm -12 <(sort "$POS_IDS") <(sort "$HMMSEARCH_HIT_IDS") | wc -l)"
      echo -e "training_sequences_missed\t$(comm -23 <(sort "$POS_IDS") <(sort "$HMMSEARCH_HIT_IDS") | wc -l)"
    } > "$RECOVERY_TSV"
}

filter_hmm_hits() {
    msg "Step 8: HMM filtering from domtblout"
    Rscript $SRCPATH/filter_hmm.R \
        $HMMSEARCH_DOMTBL \
        $BEST_HMM_SUMMARY \
        $BEST_SEQ_IDS \
        $FULL_SEQ_EVALUE \
        $BEST_DOM_IEVALUE \
        $MIN_DOM_SCORE \
        $MIN_HMM_COV \
        $MAX_BIAS_TO_SCORE

    if [[ ! -s "$BEST_SEQ_IDS" ]]; then
        echo "No HMM hits passed the balanced domain-aware filter"
        exit 0
    fi

    echo "Protein-level HMM hits retained after balanced filtering: $(wc -l < "$BEST_SEQ_IDS")"
}

flag_multidomain() {
    msg "Step 9: Flag proteins with multiple detected domains"

    awk '
    !/^#/ {count[$1]++}
    END {
      for (id in count) {
        if (count[id] > 1) print id
      }
    }
    ' "$HMMSEARCH_DOMTBL" | sort -u > "$SPLIT_DOMAIN_IDS"

    echo "Targets with >1 detected domain: $(wc -l < "$SPLIT_DOMAIN_IDS")"
}

combine_ids() {
    msg "Step 10: Combine original BLAST seeds with filtered HMM hits"

    cat $HIT_IDS $BEST_SEQ_IDS | sort -u > $COMBINED_IDS
    local combined_count
    combined_count=$(wc -l < "$COMBINED_IDS")
    echo "Combined nonredundant proteins: $combined_count"
}

extract_combined_fasta() {
    msg "Step 11: Extract combined protein sequences"

    seqkit grep -f $COMBINED_IDS $DATABASE > $COMBINED_FASTA

    local combined_nseq
    combined_nseq=$(grep -c '^>' $COMBINED_FASTA || true)
    echo "Combined sequences extracted: $combined_nseq"

    if [[ "$combined_nseq" -lt 3 ]]; then
        echo "Too few combined sequences for phylogenetic alignment"
        exit 0
    fi
}

build_combined_alignment() {
    msg "Step 12: Build combined alignment for phylogeny"
    mafft --thread $THREADS $COMBINED_FASTA > $COMBINED_ALN
}


trim_alignment() {
    msg "Step 13: Trim alignment with trimAl automated1"

    trimal -in $COMBINED_ALN -out $TRIMMED_ALN -automated1

    local trimmed_nseq
    trimmed_nseq=$(grep -c '^>' "$TRIMMED_ALN" || true)
    echo "Trimmed alignment sequences: $trimmed_nseq"
}

run_modeltest() {
    msg "Step 14: ModelTest-NG model selection"

    modeltest-ng \
      -i $TRIMMED_ALN \
      -d aa \
      -t ml \
      -p $THREADS > $MODELTEST_LOG 2>&1

}

final_summary() {
    echo 
    echo "BLAST hits table:              $HITS_TSV"
    echo "BLAST hit IDs:                 $HIT_IDS"
    echo "BLAST hit FASTA:               $OUT_FASTA"
    echo "Initial alignment:             $OUT_ALN"
    echo "HMM:                           $OUT_HMM"
    echo "HMM build log:                 $HMM_BUILD_LOG"
    echo "HMM summary stats:             $HMM_STATS"
    echo "hmmsearch table:               $HMMSEARCH_TBL"
    echo "hmmsearch domain table:        $HMMSEARCH_DOMTBL"
    echo "hmmsearch text report:         $HMMSEARCH_TXT"
    echo "Recovery summary:              $RECOVERY_TSV"
    echo "Balanced HMM summary:          $BEST_HMM_SUMMARY"
    echo "Filtered HMM IDs:              $BEST_SEQ_IDS"
    echo "Split-domain targets:          $SPLIT_DOMAIN_IDS"
    echo "Combined IDs:                  $COMBINED_IDS"
    echo "Combined FASTA:                $COMBINED_FASTA"
    echo "Combined alignment:            $COMBINED_ALN"
    echo "Trimmed alignment:             $TRIMMED_ALN"
    echo "ModelTest-NG log:              $MODELTEST_LOG"
}

check_tools
check_inputs

# Step 1
chems
echo "Running BLASTP seed search"
blast_search
sleep 2

# Step 2
chems
echo "Extracting BLAST hits"
extract_blast_hits
sleep 2

# Step 3
chems
echo "Building seed alignment with MAFFT"
build_seed_alignment
sleep 2

# Step 4
chems
echo "Building profile HMM"
build_hmm
sleep 2

# Step 5
chems
echo "Summarizing HMM"
summarize_hmm
sleep 2

# Step 6
chems
echo "Searching HMM against database"
search_hmm
sleep 2

# Step 7
chems
echo "Assessing recovery of seed proteins"
recovery_assessment
sleep 2

# Step 8
chems
echo "Applying balanced HMM filters"
filter_hmm_hits
sleep 2

# Step 9
chems
echo "Flagging multi-domain targets"
flag_multidomain
sleep 2

# Step 10
chems
echo "Combining BLAST seeds and filtered HMM hits"
combine_ids
sleep 2

# Step 11
chems
echo "Extracting combined FASTA"
extract_combined_fasta
sleep 2

# Step 12
chems
echo "Building final alignment"
build_combined_alignment
sleep 2

# Step 13
chems
echo "Trimming alignment"
trim_alignment
sleep 2

# Step 14
# chems
# echo "Selecting substitution model with ModelTest-NG"
# run_modeltest
# sleep 2

# Step 15
chems
final_summary
sleep 2

mv ${PREFIX}.* $PREFIX
