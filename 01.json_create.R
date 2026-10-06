# This script aims to build JSON files for Alphafold3. The idea is to use UP Ids
# for bot. query and subject hits, which will be chain A and chain B.

# Update: since I'm preparing this script for AGO predictions I'll also incorporate
# an option to make it able to accept RNA as chain R

# Update: I will also incorporate a flag to be able to use a single FASTA file 
# as query, so now the input does not depended only on the UP IDs 
pacman::p_load(optparse, dplyr, Biostrings, jsonlite)

opt_lst = list(
  make_option(c("-u", "--uniprot_path"), type = "character",
              default = "/home/isaac/backup/ALPHAFOLD3/PROTEINS/",
              help = "Path to directory containing UniProt FASTA [default = %default]"),
  make_option(c("-f", "--uniprot_file"), type = "character",
              default = "uniprotkb_database_WormBase_2025_04_03.fasta.gz",
              help = "UniProt FASTA filename (can be .gz) [default = %default]"),
  make_option(c("-t", "--target_path"), type = "character",
              default = "/home/isaac/backup/ALPHAFOLD3/TABLES/",
              help = "Path to directory containing targets file [default = %default]"),
  make_option(c("-l", "--targets_file"), type = "character",
              default = "D2005.1_targets.txt",
              help = "Targets filename (one ID per line or first column) [default = %default]"),
  make_option(c("-j", "--template_json_path"), type = "character",
              default = "/home/isaac/backup/ALPHAFOLD3/SRC/",
              help = "Path to directory containing template JSON [default = %default]"),
  make_option(c("-k", "--template_json_file"), type = "character",
              default = "template.json",
              help = "Template JSON filename [default = %default]"),
  make_option(c("-o", "--out_path"), type = "character",
              default = "/home/isaac/backup/ALPHAFOLD3/INPUTS/",
              help = "Outpuut directory for JSON files [default = %default]"),
  make_option(c("-q", "--query_protein"), type = "character", default = NULL,
              help = "Single query protein ID (used as protein A) [default = %default]"),
  make_option(c("-p", "--query_path"), type = "character", default = "/home/isaac/backup/ALPHAFOLD3/TABLES/",
              help = "Path to directory containing query file [default = %default]"),
  make_option(c("-s", "--query_file"), type = "character", default = NULL,
              help = "query filename (one ID per line or first column) [default = %default]"),
  make_option(c("-n", "--name_prefix"), type = "character", default = "CELE",
              help = "prefix for output JSON names [default = %default]"),
  make_option(c("--with_rna"), action = "store_true", default = FALSE,
              help = "If set, include an RNA sequence in the JSON [default = %default]"),
  make_option(c("--rna_sequence"), type = "character", default = NULL,
              help = "RNA sequence to include when --with_rna is set [default = %default]"),
  make_option(c("--rna_id"), type = "character", default = "R",
              help = "ID for RNA chain when --with_rna is set [default = %default]"),
  make_option(c("--overwrite"), action = "store_true", default = FALSE,
              help = "Overwrite existing JSON files [default = %default]"),
  make_option(c("--query_fasta"), type = "character", default = NULL,
              help = "FASTA file containing one or more query protein sequences [default = %default]"),
  make_option(c("--model_seeds"), type = "character", default = "5,10",
              help = "Commaseparated list of model seeds, e.g. 1,2,3,4,5 [default = %default]")
)

opt = parse_args(OptionParser(option_list = opt_lst))

seed_vec = as.integer(trimws(strsplit(opt$model_seeds, ",")[[1]]))

if (length(seed_vec) == 0 || any(is.na(seed_vec))) {
  stop("--model_seeds must be a commaseparated list of integers, 1,2,3,4")
}


flg_sngl_q  = !is.null(opt$query_protein) && nzchar(opt$query_protein)
flg_q_file  = !is.null(opt$query_path) && nzchar(opt$query_path) &&
  !is.null(opt$query_file) && nzchar(opt$query_file)
flg_q_fasta = !is.null(opt$query_fasta) && nzchar(opt$query_fasta)

n_qmodes = sum(c(flg_sngl_q, flg_q_file, flg_q_fasta))

if (n_qmodes == 0) {
  stop("I must provide exactly one of: --query_protein, --query_file, or --query_fasta")
}

if (n_qmodes > 1) {
  stop("use only one of: --query_protein, --query_file, or --query_fasta")
}

up_dir      = opt$uniprot_path
up_fname    = opt$uniprot_file
tgt_dir     = opt$target_path
tgt_fname   = opt$targets_file
tmpl_dir    = opt$template_json_path
tmpl_fname  = opt$template_json_file
jf_pfx      = opt$name_prefix
jf_outdir   = opt$out_path
flg_owrt    = opt$overwrite
flg_rna     = opt$with_rna
rna_seq     = opt$rna_sequence
rna_cid     = opt$rna_id

up_fp       = file.path(up_dir,   up_fname)
tgt_fp      = file.path(tgt_dir,  tgt_fname)
tmpl_fp     = file.path(tmpl_dir, tmpl_fname)

if (flg_q_file) {
  qf_fp = file.path(opt$query_path, opt$query_file)
} else {
  qf_fp = NULL
}

if (flg_q_fasta) {
  qfa_fp = opt$query_fasta
} else {
  qfa_fp = NULL
}

if (!file.exists(up_fp))   stop("UniProt FASTA not found: ",   up_fp)
if (!file.exists(tgt_fp))  stop("Targets file not found: ",    tgt_fp)
if (!file.exists(tmpl_fp)) stop("tempolate not found: ",   tmpl_fp)

if (!is.null(qf_fp) && !file.exists(qf_fp)) {
  stop("Q file not found: ", qf_fp)
}

if (!dir.exists(jf_outdir)) {
  dir.create(jf_outdir, recursive = TRUE, showWarnings = FALSE)
}

if (flg_rna) {
  if (is.null(rna_seq) || !nzchar(rna_seq)) {
    stop("--with_rna is set, you need to pass an RNA!")
  }
  
  rna_seq = toupper(gsub("\\s+", "", rna_seq))
  rna_seq = chartr("T", "U", rna_seq)
  
  if (grepl("[^ACGU]", rna_seq)) {
    stop("RNA sq contains invalid characters!")
  }
}

if (!is.null(qfa_fp) && !file.exists(qfa_fp)) {
  stop("Q FASTA not found: ", qfa_fp)
}

up_aa = readAAStringSet(up_fp)

names(up_aa) = names(up_aa) %>%
  strsplit("\\|") %>%
  sapply(function(x) {
    if (length(x) >= 2) x[2] else x[1]
  })

tgt_raw = read.delim(tgt_fp, sep = "\t", header = FALSE)
tgt_ids = trimws(as.character(tgt_raw[[1]]))
tgt_ids = tgt_ids[nzchar(tgt_ids)]

if (length(tgt_ids) == 0) stop("No valid T IDs found in targets file")

if (flg_sngl_q) {
  qry_ids = trimws(opt$query_protein)
  qry_ids = qry_ids[nzchar(qry_ids)]
  
  if (length(qry_ids) == 0) stop("No valid Q IDs found")
  
  miss_q = setdiff(qry_ids, names(up_aa))
  if (length(miss_q) > 0) {
    warning("These query IDs were not found in UniProt and will be ommited: ",
            paste(miss_q, collapse = ", "))
  }
  qry_ids = intersect(qry_ids, names(up_aa))
  
  if (length(qry_ids) == 0) stop("None of the Q IDs were found in the UniProt fasta")
  
  qry_aa = up_aa[qry_ids]
  
} else if (flg_q_file) {
  qf_raw  = read.delim(qf_fp, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  qry_ids = trimws(as.character(qf_raw[[1]]))
  qry_ids = qry_ids[nzchar(qry_ids)]
  
  if (length(qry_ids) == 0) stop("No valid Q IDs found")
  
  miss_q = setdiff(qry_ids, names(up_aa))
  
  if (length(miss_q) > 0) {
    warning("These query IDs were not found in UniProt and will be omited: ",
            paste(miss_q, collapse = ", "))
  }
  
  qry_ids = intersect(qry_ids, names(up_aa))
  qry_aa = up_aa[qry_ids]
  
} else if (flg_q_fasta) {
  qry_aa = readAAStringSet(qfa_fp)
  names(qry_aa) = names(qry_aa) %>%
    sapply(function(x) {
      x = trimws(x)
      if (is.na(x) || !nzchar(x)) return("query")
      strsplit(x, "\\s+")[[1]][1]
    })
  
  if (anyDuplicated(names(qry_aa))) {
    dup_nms = unique(names(qry_aa)[duplicated(names(qry_aa))])
    stop("Duplicate sequence names in query FASTA: ",
         paste(dup_nms, collapse = ", "),
         ". Make FASTA headers unique")
  }
}

miss_t = setdiff(tgt_ids, names(up_aa))
if (length(miss_t) > 0) {
  warning("These target IDs were not found in UniProt and will be skipped: ",
          paste(miss_t, collapse = ", "))
}
tgt_ids = intersect(tgt_ids, names(up_aa))
tmpl = fromJSON(tmpl_fp, simplifyVector = FALSE)

jf_written = c()

for (qa_id in names(qry_aa)) {
  seq_A = as.character(qry_aa[qa_id])
  for (tb_id in tgt_ids) {
    seq_B = as.character(up_aa[tb_id])
    jf      = tmpl
    jf$modelSeeds                    = as.list(seed_vec)
    jf$sequences[[1]]$protein$id     = "A"
    jf$sequences[[1]]$protein$sequence = seq_A
    jf$sequences[[2]]$protein$id     = "B"
    jf$sequences[[2]]$protein$sequence = seq_B
    
    if (flg_rna) {
      jf$sequences[[3]] = list(rna = list(id = rna_cid, sequence = rna_seq))
    } else {
      if (length(jf$sequences) > 2) jf$sequences = jf$sequences[1:2]
    }
    
    jf$name = paste0(jf_pfx, "_", qa_id, "_", tb_id)
    jf_out  = file.path(jf_outdir, paste0(jf$name, ".json"))
    
    if (file.exists(jf_out) && !flg_owrt) {
      warning("Skipping existing file: ", jf_out)
      next
    }
    
    write_json(jf, path = jf_out, auto_unbox = TRUE, pretty = TRUE, digits = NA)
    jf_written = c(jf_written, jf_out)
  }
}