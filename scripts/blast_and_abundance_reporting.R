library(data.table)
library(openxlsx)
library(fuzzyjoin)
##library(Biostrings)
##library(tidyverse)
## install.packages("R.utils") # For fread to process .gz, otherwise gunzip .gz to .tsv first.

#TODO: Turn this into package.

#TODO: Cleanup style, jumping between 'DT' as name and other descriptive names.

#TODO: Cleanup long parameter choice styles?

## data.table idioms used herein:
## https://cran.r-project.org/web/packages/data.table/vignettes/datatable-faq.html#why-do-i-have-to-type-dt-sometimes-twice-after-using-to-print-the-result-to-console
## https://github.com/Rdatatable/data.table/issues/869
##
## When using ':=', use '[]' after the last ':=' in functions.
##
## NO : DT[, V1 := NULL]
## YES: DT[, V1 := NULL][]

##' Get column names from file.
##'
##' @param fname A filename for file containing space separated column
##'   names.
##' @return vector of column names.
##' @export
get_col_names <- function(fname){
  col_names <- scan(fname,
                    sep= " ",
                    what= "character",
                    quiet = TRUE)
  make.names(col_names,
             unique = TRUE)
}


##' Get and name blast table from data and column name files.
##'
##' This read method is appropriate for blast results reported in
##' outfmt6 tabular format.
##'
##' Reading blast results is complicated by the occasional use of
##' different quote styles in some fields. And the possible presence
##' of common comment characters, such as '#'.
##'
##' We only tested with tsv seperated.
##'
##' @param blast_fname A filename for file with blast output in tsv
##'   format, such as obtained from blast outfmt 6 and derivatives.
##' @param col_names_fname A filename for file with single line of
##'   space separated column names.
##' @return data.table of blast results.k
##' @export
blast_files_to_dt <- function(fname, fname_columns){
  col_names <- get_col_names(fname_columns)

  ## data.table's fread magically handles the quotes and comment chars.
  ##
  ## BL <- fread(fname, col.names = col_names)
  ## t <- read.table(fname,
  ##                 col.names = col_names,
  ##                 comment.char = "",
  ##                 quote = "",
  ##                 sep = "\t",
  ##                 as.is = TRUE)
  ## compare(t, BL)
  ## compareEqual(t, BL)
  ## compareEqual(as.data.table(t), BL) # TRUE
  fread(fname, col.names = col_names)
}


##' Reduces blast tables to single best query rows.
##'
##' If the ordering of the blast table has not been changed from the
##' ordering that blast itself outputs, the first encountered
##' query-target is the best query-target.
##'
##' All fields except query are ignored, but included in
##' output. Therefore, even if there are multiple HSPs only the HSP
##' that is first encountered is kept.
##'
##' About query column names:
##' - qseqid is the query column for diamond blastx on NR.
##' - qaccver is the query column for blastn (-tasks blastn
##'   and dc-megablast) on NT.
##' - contig is also acceptable.
##'
##' @param DT blast data.table containing a "query" column which is of
##'   "qaccver" or "qseqid" or "contig"
##' @return deduplicated blast data.table containing only the
##'   "best hits".
##' @export
blast_to_best_hits <- function(DT,
                               dedup_on = c("qaccver",
                                            "qseqid",
                                            "contig")){
  dedup_on <- match.arg(dedup_on)
  DT[!duplicated(eval(as.symbol(dedup_on))),]
}


##' Combine blast hits and abundance.
##'
##' Merges on "target_id" in abundance to ("qaccver", "qseqid", or
##' "contig"). Only one of the latter can be present. The merge column
##' is renamed to "contig".
##'
##' Provides sensible names for various length columns, if present.
##'
##' @param abund data.table with abundance.
##' @param blast data.table with blast hits.
##' @return The merged data.table
##' @export
merge_abund_blast <- function(abund, blast){
  blast_merge_col_vals <- c("qseqid", "qaccver", "contig")

  found <- names(blast) %in% blast_merge_col_vals

  if ((!any(found)) || (!sum(found) == 1)) {
    stop(paste("Error: Expect a single value to be present in blast column names for merging which can be one and only one of :", paste(blast_merge_col_vals, collapse = " ")))
  }

  y <- names(blast)[found]

  DT <- merge(abund,
              blast,
              by.x = "target_id",
              by.y = y,
              all.x = TRUE) # 1st cols becomes "target_id", "length",
                            # and there is a match "length.y" which
                            # might be larger than qlen because of
                            # gaps in alignments.

  names(DT)[names(DT) ==
              "target_id"] <- "contig" # In blast the "target" is what
                                       # is hit. Rename to contig to
                                       # clarify this is not the blast
                                       # target.
  names(DT)[names(DT) ==
              "length.x"] <- "contig_length" # So we can disambiguate
                                             # the various lengths.
  if ("length.y" %in% names(DT)) {
    names(DT)[names(DT) ==
                "length.y"] <- "match_length" # Note, it might be
                                              # longer than "qlen"
                                              # because of gaps and
                                              # indels.
  }

  if ("qlen" %in% (names(DT))){ # Our extended blastn (outfmt6plus)
    DT[, qlen := NULL][] # Drop redundant column (contig_length == qlen)
  }

  DT # Note, no more "qseqid" or "qaccver"
}


##' Extracts sample names to new column.
##'
##' Assumes the sample name can be extracted by simple substitution of
##' "" for the pattern.
##'
##' e.g. applying the pattern "_TRINITY_.*" and replacement "" to
##' "1_GB3_B_TRINITY_DN4_c0_g1_i8" adds "1_GB3_B" to a "sample"
##' column.
##'
##' @param DT data.table containing column which has sample names
##'   embedded.
##' @param column The column name containing the sample name. Probably contig.
##' @param pattern The pattern to delete from the column in order to
##'   produce a sample name.
##' @param replacement The replacement for the pattern.
##' @return data.table with a new "sample" column.
##' @export
extract_add_sample_name <- function(DT,
                                    column,
                                    pattern,
                                    replacement){
  DT[, "sample" := sub(pattern, replacement, get(column, DT))][]
}

##' Sorts blast hits by sample, then tpm.
##'
##' @param DT data.table of abundance and hits.
##' @param abund_sort_by Column name of abundance to sort on. Default
##'   "tpm".
##' @return data.table sorted on sample name and abundance.
##' @export
sort_sample_abund_blast <- function(DT,
                                    abund_sort_by = c("tpm", "est_counts")){
  abund_sort_by <- match.arg(abund_sort_by)
  data.table::setorderv(DT, c("sample", abund_sort_by), c(1, -1))
}

##' Merge and sort abundance and blast data.tables.
##'
##' Sorting is by sample name, then abundance (tpm or est_counts).
##'
##' @param abund data.table with abundance.
##' @param blast data.table with blast hits
##' @param pattern Pattern to be replaced.
##' @param replacement Replacement for pattern.
##' @param col_w_sample_info The column containing the sample info.
##' @param abund_sort_by The column to sort on (decreasing sort).
##' @return data.table with merged abundance and blast sorted by
##' @export
merge_sort_abund_blast <- function(abund,
                                   blast,
                                   pattern = "_TRINITY_.*",
                                   replacement = "",
                                   col_w_sample_info = "contig",
                                   abund_sort_by = c("tpm", "est_counts")){

  abund_sort_by <- match.arg(abund_sort_by)

  DT <- merge_abund_blast(abund, blast)
  DT <- extract_add_sample_name(DT,
                                column = col_w_sample_info,
                                pattern = pattern,
                                replacement = replacement)
  sort_sample_abund_blast(DT, abund_sort_by)
}

##' Splits blast results by column.
##'
##' Useful to report by kingdom or by species reports.
##'
##' @param DT data.table
##' @param split_column
##' @return list of data.frames or data.tables.
##' @export
split_by_column <- function(DT, split_column){
  split(DT, get(split_column, DT))
}


##' Write a blast data.table to file.
##'
##' If ftype is "tsv", the file will be identical to how we get
##' blast tsv files directly from blast itself.
##'
##' If ftype is "csv", or "xlsx" the file will be easily opened by
##' excel.
##'
##' For excel introducing spaces into column names is useful so that
##' in excel, the names can wrap more nicely if 'wrap' is chosen
##' within excel.
##'
##' Does not add extension to filename. User is expected to provide
##' valid extension for the filetype chosen.
##'
##' @param dt Data.table containing data to be written.
##' @param file Filename to write to, including the extension.
##' @param ftype File type.
##' @param xlsx_sheet_name Name to use for excel sheet.
##' @param xlsx_replace_char_w_space Character to replace with a space
##'   for nicer column name wrapping.
##' @export
write_blast_dt <- function(dt,
                           file,
                           ftype = c("tsv", "csv", "xlsx"),
                           xlsx_sheet_name = "Sheet1",
                           xlsx_replace_char_w_space = "."){
  if (ftype == "tsv") { # identical formatting as our input.
    fwrite(dt,
           file = file,
           col.names = FALSE,
           sep = "\t",
           quote = FALSE)
  } else if (ftype == "csv") { # identical formating as write.csv
    fwrite(dt,
           file = file,
           row.names = TRUE,
           quote = TRUE)
  } else if (ftype == "xlsx") {
    old_col_names <- colnames(dt)
    colnames(dt) <- gsub(xlsx_replace_char_w_space,
                         " ",
                         colnames(dt),
                         fixed = TRUE)
    header_style <- createStyle(textDecoration = "bold")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb,
                           sheetName = xlsx_sheet_name)
    openxlsx::writeDataTable(wb,
                             xlsx_sheet_name,
                             x = dt,
                             rowNames = TRUE,
                             headerStyle = header_style)
    openxlsx::freezePane(wb,
                         xlsx_sheet_name,
                         firstRow = TRUE)
    openxlsx::setColWidths(wb,
                           xlsx_sheet_name,
                           cols = 1:ncol(dt),
                           widths = 12)
    openxlsx::saveWorkbook(wb,
                           file = file,
                           overwrite = TRUE)
    colnames(dt) <- old_col_names
  } else {
    stop("Supported file types are tsv, csv, xlsx")
  }
}


##' Tries to insert qstrand and sstrand info immediately after bitscore.
##'
##' If no qstrand is found, tries to calculate qstrand from qstart and
##' qend.
##'
##' If no sstrand is found, tries to calculate sstrand from sstart and
##' send.
##'
##' If no bitscore column is found, qstrand and sstrand, if they exist, will remain in place. Otherwise, they are added to the end.
##'
##' Data.table will be edited by reference.
##'
##' @param DT data.table possibly containing strand or position info.
##' @export
get_and_place_strand_cols <- function(DT){
  if (!"qstrand" %in% names(DT) & all(c("qend", "qstart") %in% names(DT))){
    DT[,qstrand := ifelse(qend - qstart > 0, "plus", "minus")][]
  }
  if (!"sstrand" %in% names(DT) & all(c("send", "sstart") %in% names(DT))){
    DT[,sstrand := ifelse(send - sstart > 0, "plus", "minus")][]
  }

  bitscore_pos <- match("bitscore", names(DT))
  qstrand_pos <- match("qstrand", names(DT))
  sstrand_pos <- match("sstrand", names(DT))
  if (!any(is.na(c(bitscore_pos, qstrand_pos, sstrand_pos)))) {
    before <- names(DT)[1:bitscore_pos]
    after <- names(DT)[(bitscore_pos + 1):ncol(DT)]
    after <- after[!after %in% c("qstrand", "sstrand")]
    setcolorder(DT,
                c(before,
                  "qstrand",
                  "sstrand",
                  after))
  }
}

##' Subset columns to a few that we want to report in short format.
##'
##' Optionally apply a prefix or suffix to help disambiguate various
##' blast reports if they will later be joined.
##'
##' Standard blast columns reported are:
##' c("match_length", "evalue",  "staxids", "sseqid")
##'
##' Additional blast columns reported (if present):
##'
##' - diamond : c("sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")
##' - blast :  c("ssciname", "sskingdom", "stitle")
##'
##' @param DT data.table of merged abundance and blast.
##' @param prefix Prefix to apply to blast columns.
##' @param suffix Sufffix to apply to blast columns.
##' @param diamond Is blast from diamond?
##' @param with_rank Is within group rank column desired?
##' @param cols_to_report vector of column names to report.
##' @return data.table subsetted to just the report columns.
##' @export
report_key_cols <- function(DT,
                            prefix = "",
                            suffix = "",
                            diamond = FALSE,
                            with_rank = FALSE){

  blast_cols <- c("sscinames",
                  "evalue",
                  "qstrand",
                  "sstrand",
                  "stitle",
                  "match_length",
                  "staxids")

  if (diamond) {
    blast_cols <- c(blast_cols,
                    "sseqid",
                    "sskingdoms",
                    "skingdoms",
                    "sphylums")
  } else {
    blast_cols <- c(blast_cols,
                    "saccver",
                    "sskingdom",
                    "sskingdoms")
  }
  if (with_rank) {
    blast_cols <- c(blast_cols, "rank_of_hit")
  }

  abund_cols <-  c("sample",
                   "contig",
                   "contig_length",
                   "est_counts",
                   "tpm")

  cols_to_report <- c(abund_cols, blast_cols)
  DT <- DT[, ..cols_to_report]

  new_cols <- paste0(prefix, blast_cols, suffix)
  setnames(DT, old = blast_cols, new = new_cols)
  DT
}


##' Pulls out the single best hit per query. (Possibly as prioritized
##' by selected superkingdom).
##'
##' This is a high level summary function for our dataset combining
##' blast and abundance data.
##'
##' Abundance is in a form such as the "abundance.tsv" trinityrnaseqs
##' util script "align_and_estimate_abundance.pl" produces.
##'
##' Matches as a pattern on selected_sskingdoms, so compound
##' sskingdoms are extracted together.
##'
##'   e.g. "Viruses" might retrieve the following compound and simple
##'   super kingdoms:
##'
##'   - 0;Viruses
##'   - Bacteria;Viruses
##'   - Eukaryota;Viruses
##'   - Viruses
##'
##' The rank of the first encountered hit is reported. The best hit is
##' 1. If a contig has a best hit that is not in the selected
##' sskingdoms, but the next best hit is in the selected sskingdoms,
##' the ran_of_hit is 2.
##'
##'   e.g. A contig might have 10 hits because we used
##'   (-max_target_seqs 10) as a blast parameter. The best 7 might be
##'   to the superkingdom "Bacteria". The next 2 to "Viruses", and the
##'   last to "Archaea". If "Viruses" is selected as the sskingdoms to
##'   report, only the hit in position 8 will be reported, along with
##'   its rank '8'. If "Archaea" is selected as the sskingdoms, only
##'   the hit at position '9' is reported, along with its rank
##'   '9'. Likewise, selecting "Bacteria" would return the hit at
##'   position 1.
##'
##' Selecting a superkingdom can be overridden using an appropriate
##' pattern.
##'
##'   e.g. The selected_sskingdoms patterns "." or "" will select all
##'   sskingdoms.
##'
##' Multiple superkingdoms can be selected using an appropriate
##' pattern.
##'
##'   e.g. The selected_sskingdoms pattern "Euk|Bac" will select the
##'   best hit among Eukaryotes and Bacteria.
##'
##' Diamond reports the query as qseqid. NCBI blast reports reports
##' the query as qaccver. Some column names we report are specific to
##' NCBI blast. Some to diamond. Therefore we switch between them.
##'
##' Prefix and suffix can be useful if one will later merge results of different
##' types of blast. Helping to distinguish which blast was the source.
##'
##' As a summary for our data, it may be too specific to our data,
##' such that it fails yours. If so, use it as a template to design
##' your own reporting functions. For example, we have sample names
##' encoded in our query names and we parse those sample names out by
##' trimming removing the "_Trinity_.*" suffix (the default with the
##' "merge_sort_abund_blast" function). If your sample names in the
##' queries are differently encoded, modify that call per the
##' documentation for "merge_sort_abund_blast".
##'
##' @param abund data.table with abundance. Such as from running
##'   kallisto trinityrnaseq.
##' @param blast data.table with blast hits.
##' @param diamond TRUE or FALSE. Sets query id column name and
##'   reporting columns that might be specific to diamond blast.
##' @param prefix Prefix to apply to blast column names.
##' @param suffix Suffix to apply to blast column names.
##' @param selected_sskingdoms Pattern to use to select
##'   sskingdoms. Defaults to all sskingdoms.
##' @param best TRUE or FALSE. Return best hit (for selected kingdom).
##' @return A data.table with abundance and certain blast column
##'   serving as a high level overview of the best blast hits per
##'   query for a given superkingdom.
##' @export
blast_abund_sskingdoms_key_cols <- function(abund,
                                            blast,
                                            diamond = FALSE,
                                            prefix = "",
                                            suffix = "",
                                            selected_sskingdoms = "",
                                            best = TRUE){
  b <- copy(blast)
  if (diamond) {
    contig_col = "qseqid"
  } else {
    contig_col = "qaccver"
  }
  b[, rank_of_hit := seq_len(.N), by = contig_col][]
  b <- b[sskingdoms %like% selected_sskingdoms]
  if (best) {
    b <- blast_to_best_hits(b, dedup_on = contig_col)
  }
  ab <-  merge_sort_abund_blast(abund = abund, blast = b)
  ab <- ab[!is.na(evalue)]
  report_key_cols(ab, with_rank = TRUE, prefix = prefix, diamond = diamond)
}

##' Collects 3 blast reports from an abundance file and a blast file.
##'
##' See return for details on the reports.
##'
##' Note, in producing blast results, a strand column(s) might be
##' calculated if not present, or repositioned if present.
##'
##' @param abund data.table of abundance, such as from kallisto.
##' @param blast_fname filename of the blast results tsv.
##' @param blast_col_name filename of the blast column names, space
##'   separated.
##' @param prefix string to apply as a prefix to certain report
##'   columns. Useful to distinguish source of results.
##' @param ssk ssuperkingdoms to subset the report on. Used as a
##'   pattern. So "Viruses" will recover both "Viruses" and
##'   "Eukaryota;Viruses". Converted to lower case in results.
##' @param diamond Boolean. TRUE specifies the blast results are from
##'   diamond blastx. Switches to certain reporting columns that are
##'   specific to diamond or where the column names are different from
##'   NCBI blasts.
##' @param best_only Boolean. TRUE specifies that we only want the
##'   first (best) encountered hit. Or, when the ssk is not "", the
##'   1st encountered ssuperkingdom matching ssk. For example, ssk =
##'   "Viruses" and best_only = TRUE returns the first encountered
##'   Virus hit for the summary_ssk results.
##' @param paste_char Character for building column names. default is "."
##' @return list of blast reports.
##' \itemize{
##'   \item blast - Full blast results. There is no abundance data attached.
##'   \item summary - As in blast, but with abundance data attached, fewer blast specific columns, and possibly reduced to a single best row for each query.
##'   \item summary_ssk - As in blast, but with abundance data attached, fewer blast specific columns, and possibly reduced to only a single best ssuperkingdom row for each query. For example, if ssk = "Viruses", the best Virus hit is reported, whether or not it was the overall best hit.
##' }
process_reports <- function(abund,
                            blast_fname,
                            blast_col_name,
                            prefix = "",
                            ssk = "Viruses",
                            diamond = FALSE,
                            best_only = TRUE,
                            paste_char = "."){
  b <- blast_files_to_dt(blast_fname, blast_col_name)
  get_and_place_strand_cols(b)
  p_all <- paste0(prefix, "all", paste_char)
  summary <- blast_abund_sskingdoms_key_cols(a,
                                             b,
                                             prefix = p_all,
                                             selected_sskingdoms = "",
                                             diamond = diamond,
                                             best = best_only)
  p_sel <- paste0(prefix, tolower(ssk), paste_char)
  summary_ssk <- blast_abund_sskingdoms_key_cols(a,
                                                 b,
                                                 prefix = p_sel,
                                                 selected_sskingdoms = ssk,
                                                 diamond = diamond,
                                                 best = best_only)
  list(blast = b, summary = summary, summary_ssk = summary_ssk)
}


##' For each contig report all hits with equivalent or better rank the
##' best ssuperkingdom hit.
##'
##' Use this to limit summary reports for review. Usually we want to
##' review hits better than the hit for the ssuperkingdoms and do not
##' care about hits that are worse.
##'
##' For example: Run this report. Remove the contigs for the
##' sskingdoms that were best hits in order to focus on more
##' questionably reported hits for that sskingdoms. Write the
##' remainder to an excel file for manual review. This should cut down
##' substantially the amount of lines for each contig needing to be
##' reviewed, making these manual reviews more efficient.
##'
##' copy(dc_all$
##' @param best data.table A summary report of the single best hit for
##'   the ssuperkingdoms. This should probably be a copy of a report
##'   data.table. See details. ex. copy(dc$summary_vir)
##' @param all data.table A summary report of all hits all
##'   irregardless of ssuperkingdoms. This should probably be a copy of a report data.table. See details.
##' @param ssk The ssuperkingdom used in the 'best' data.table
##' @return data.table of all hits equal with rank to or better than
##'   the best hit of the ssk.
##' @export
##' @examples
##' \dontrun{
##' a <- fread("kallisto_abund_file.tsv")
##' a[, eff_length := NULL] # Drop redundant col (qlen = eff_length).
##' # summaries contain best hits.
##' dc <- process_reports(a,
##'                       "dc_megablast_file.tsv",
##'                       "dc_megablast_blast_cols.txt")
##' # summaries contain all hits.
##' dc_all <- process_reports(a,
##'                           "dc_megablast_file.tsv",
##'                           "dc_megablast_blast_cols.txt",
##'                           best_only = FALSE)
##' best <- copy(dc$summary_ssk)
##' names(best) <- sub("Viruses", "", names(best))
##' all <- copy(dc_all$summary)
##' names(all) <- sub("all ", "", names(all))
##'
##' for_review <- get_ranks_gte_ssk(best, all)
##' for_review <- for_review[!(rank_of_hit == 1 & sskingdoms %like% "Viruses")]
##' write_blast_dt(for_review, file = "review_these.xlsx", ftype = "xlsx", sheet_name = "for_review")
##'}
get_ranks_gte_ssk <- function(best, all, ssk = "Viruses"){
  ## Subset all hits to the contigs hit by the ssk. Find the best hit
  ## for the ssk. Then non-equi join to get ranks (of anything) <= the
  ## best virus rank. Note, we copy the original "rank_of_hit" because
  ## the join loses it to the rollup value otherwise.
  check <- all[contig %in% best$contig]
  check_ssk <-  check[`sskingdoms` %like% ssk,
                      list(best_ssk = min(`rank_of_hit`)),
                      by = contig]
  check[, roll_rank := `rank_of_hit`]
  check <- check[check_ssk, on=c("contig", "roll_rank<=best_ssk")]
  check[,best_ssk := `rank_of_hit` == roll_rank][]
}

##' Extracts fasta set to just those with matching names (ids).
##'
##' Subset a fasta file based on (part of) their names (ids). In our case,
##' we match upto first literal space character (' ') of the fasta header.
##' After matches are found, we join back the rest of the fasta header,
##' using that to subset the fastas.
##'
##' @param fastas_fname Filename of fastas to be subsetted.
##' @param ids Data.table of ids to subset on and info to add to fasta
##'   description.
##' @export
##' @return Subsetted fastas, possibly with new info added to fasta
##'   descriptions.
##' @examples
##' \dontrun{
##' subsetted <- extract_fastas_by_ids(test.fasta, ids_info)
##' Biostrings::writeXStringSet(subsetted_fastas,
##'   file = "test.subsetted.fasta",
##'   compress = FALSE)
##' }
extract_fastas_by_ids <- function(fastas_fname,
                                  ids) {
  require("Biostrings")
  require("tidyverse")
  fa <- Biostrings::readDNAStringSet(fastas_fname)
  ## ids <- suppressMessages(read_csv(ids_fname,
  ##                                  comment = "",
  ##                                  col_names = c("ids")))
  ids <-
  fas <- tibble::tibble(names(fa))

  names(fas) <- c("name")
  fas <- fas %>%
    separate(col = name,
             into = c("ids", "info"),
             sep = " ",
             extra = "merge",
             fill = "right")
  fas <- left_join(ids, fas, by = "ids") # ordering is by ids$ids
  fas <- fas %>%
    unite(name, sep = " ")

  fa[fas$name]
}
