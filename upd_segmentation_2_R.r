#######################################################
# segmentation_2.R : segmentation script for UPD mode
# Created by Nahyun Kong (nahyun@wustl.edu)
# Last edited 07/2025.

## Annotated pipeline for segmenting child's VAF data per groupA and groupB snps, combine segments per groups, and classify UPd .
#######################################################

# --- 1. Load & Install Required Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("PSCBS", "DNAcopy", "optparse", "dplyr", "tidyverse", "GenomicRanges", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("PSCBS", "DNAcopy", "GenomicRanges")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}
suppressMessages(library(PSCBS))
suppressMessages(library(DNAcopy))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tibble))
suppressMessages(library(purrr))
options(warn=-1)

# --- 2. Parse Command-Line Arguments ---
option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, help="Input counts file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="./", help="Output directory (default: ./)", metavar="character"),
  make_option(c("-s", "--segment_length"), type="double", default=1000000, help="Minimum segment length in MB", metavar="integer")
  )
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

counts_path = opt$input_file
output_dir = opt$output_dir
segment_length = opt$segment_length

counts_df = read_delim(counts_path, delim='\t')
#counts_df = read_delim("/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20240716_final_submission/rerun/hg38/triomix_snp/1-00180.child_chr14.counts", delim='\t')

# --- 3. Helper Functions for Chromosome Conversion ---
convert_chr_to_integer <- function(x){
  if(str_detect(x, '[0-9]+')){
    as.integer(gsub(pattern = 'chr', replacement = '', x=x))
  }else{
    23 # chrX -> 23
  }
}
convert_chr_to_integer_vec = Vectorize(convert_chr_to_integer)
convert_integer_to_chrom <- function(x){
  if(x==23){
    return('chrX')
  }else{
    return(paste0('chr', x))
  }
}
convert_integer_to_chrom_vec = Vectorize(convert_integer_to_chrom)

# --- 4. Load and Prepare Input Data ---
counts_df = counts_df %>% mutate(chrom=convert_chr_to_integer_vec(chrom)) %>% dplyr::rename(chromosome=chrom, x=pos, y=vaf, depth=depth)

# --- 5. Extract Snp Groups by Inheritance Patterns ---
df_homoalt_father = counts_df %>% filter(father_vaf==1 & mother_vaf==0) %>% dplyr::select(chromosome, x, y, depth)
df_homoalt_mother = counts_df %>% filter(mother_vaf==1 & father_vaf==0) %>% dplyr::select(chromosome, x, y, depth)
df_het_father = counts_df %>% filter(father_vaf > 0.4 & father_vaf < 0.6 & mother_vaf ==0) %>% dplyr::select(chromosome, x, y, depth)
df_het_mother = counts_df %>% filter(mother_vaf > 0.4 & mother_vaf < 0.6 & father_vaf ==0) %>% dplyr::select(chromosome, x, y, depth)

# --- 6. Segmentation of VAF Data by Snp Groups ---
segments = data.frame()
for(chrom in unique(counts_df$chromosome)){
  df_homoalt = list(father=df_homoalt_father %>% filter(chromosome==chrom),
                    mother=df_homoalt_mother %>% filter(chromosome==chrom))
  df_het = list(father=df_het_father %>% filter(chromosome==chrom),
                mother=df_het_mother %>% filter(chromosome==chrom))
  df_het_no0 = list(father=df_het_father %>% filter(chromosome==chrom) %>% filter(y!=0),
                    mother=df_het_mother %>% filter(chromosome==chrom) %>% filter(y!=0))
  for(parent in c('father', 'mother')){
    df_homoalt_parent = df_homoalt[[parent]]
    df_het_parent = df_het[[parent]]
    df_het_parent_no0 = df_het_no0[[parent]]
    if(dim(df_homoalt_parent)[1]>1){
      df_homoalt_parent = dropSegmentationOutliers(df_homoalt_parent)
      gaps <- findLargeGaps(df_homoalt_parent, minLength=segment_length * 1e+06)
      knownSegments <- gapsToSegments(gaps)
      fit <- segmentByCBS(df_homoalt_parent, knownSegments = knownSegments, seed = 48879, verbose = -10, joinSegments= FALSE)
      segmented = getSegments(fit, simplify = TRUE)
      segmented = segmented %>% dplyr::select(-sampleName)
      segmented$parent = paste0(parent,"hom_alt")
      segmented$mean_depth = mapply(function(start, end) {
        mean(df_homoalt_parent$depth[df_homoalt_parent$x >= start & df_homoalt_parent$x <= end], na.rm=TRUE)
      }, segmented$start, segmented$end)
      segments = bind_rows(segments, segmented)
    }
    if(dim(df_het_parent)[1]>1){
      df_het_parent = dropSegmentationOutliers(df_het_parent)
      gaps <- findLargeGaps(df_het_parent, minLength=segment_length * 1e+06)
      knownSegments <- gapsToSegments(gaps)
      fit <- segmentByCBS(df_het_parent, knownSegments = knownSegments, seed = 48879, verbose = -10, joinSegments= FALSE)
      segmented = getSegments(fit, simplify = TRUE)
      segmented = segmented %>% dplyr::select(-sampleName)
      segmented$parent = paste0(parent,"het")
      segmented$mean_depth = mapply(function(start, end) {
        mean(df_het_parent$depth[df_het_parent$x >= start & df_het_parent$x <= end], na.rm=TRUE)
      }, segmented$start, segmented$end)
      

      segments = bind_rows(segments, segmented)
    }
    if(dim(df_het_parent_no0)[1]>1){
      df_het_parent_no0 = dropSegmentationOutliers(df_het_parent_no0)
      gaps <- findLargeGaps(df_het_parent_no0, minLength=segment_length * 1e+06)
      knownSegments <- gapsToSegments(gaps)
      fit <- segmentByCBS(df_het_parent_no0, knownSegments = knownSegments, seed = 48879, verbose = -10, joinSegments= FALSE)
      segmented = getSegments(fit, simplify = TRUE)
      segmented = segmented %>% dplyr::select(-sampleName)
      segmented$parent = paste0(parent,"het_no0")
      segmented$mean_depth = mapply(function(start, end) {
        mean(df_het_parent_no0$depth[df_het_parent_no0$x >= start & df_het_parent_no0$x <= end], na.rm=TRUE)
      }, segmented$start, segmented$end)
      segments = bind_rows(segments, segmented)
    }
  }
}
segments = segments %>% mutate(chromosome=convert_integer_to_chrom_vec(chromosome))

# --- 7. Build Atomic Segments ---
segments$id <- seq_len(nrow(segments))
gr <- with(segments, GRanges(
  seqnames = chromosome,
  ranges = IRanges(start, end),
  parent = parent,
  mean = mean,
  id = id
))
atomic_gr <- disjoin(gr, with.revmap=TRUE)
atomic_df <- as.data.frame(atomic_gr)
atomic_df$revmap <- as.character(atomic_df$revmap)

revmap_list <- as.list(atomic_gr$revmap)

# --- 8. Assign Means per Parental Pattern ---
assign_means <- function(revmap_str) {
  ids <- as.integer(unlist(strsplit(gsub("[c() ]", "", revmap_str), ",")))
  segments[ids, c("parent", "mean")] %>% tibble::deframe()
}
mean_lists <- lapply(atomic_df$revmap, assign_means)
parent_names <- unique(segments$parent)
get_mean_for_parent <- function(mlist, parent) {
  if (parent %in% names(mlist)) as.numeric(mlist[parent]) else NA
}
for (parent in parent_names) {
  atomic_df[[parent]] <- sapply(mean_lists, get_mean_for_parent, parent=parent)
}
names(atomic_df)[names(atomic_df) == "seqnames"] <- "chromosome"
atomic_df <- atomic_df %>% select(-revmap)

# --- 9. UPD Classification ---
atomic_df <- atomic_df %>%
  mutate(
    UPD_type = case_when(
      (chromosome != "chrX" & 
         (fatherhom_alt < 0.1 & motherhom_alt > 0.9)) ~ "maternal_UPD",
      
      (chromosome != "chrX" & 
         (motherhom_alt < 0.1 & fatherhom_alt > 0.9 )) ~ "paternal_UPD",
      
      TRUE ~ "not_UPD"
    ),
    maternal_UPD_subtype = case_when(
      UPD_type == "maternal_UPD" & !is.na(motherhet_no0) &
        motherhet_no0 >= 0.4 & motherhet_no0 <= 0.6 ~ "maternal_heterodisomy",
      UPD_type == "maternal_UPD" & motherhet_no0 > 0.8 ~ "maternal_isodisomy",
      TRUE ~ NA_character_
    ),
    paternal_UPD_subtype = case_when(
      UPD_type == "paternal_UPD" & !is.na(fatherhet_no0) &
        fatherhet_no0 >= 0.4 & fatherhet_no0 <= 0.6 ~ "paternal_heterodisomy",
      UPD_type == "paternal_UPD" & fatherhet_no0 > 0.9 ~ "paternal_isodisomy",
      TRUE ~ NA_character_
    ),
    UPD_final = coalesce(maternal_UPD_subtype, paternal_UPD_subtype, UPD_type, "not_UPD")
  )

# --- 10. Write Atomic Segments Table ---
output_path = file.path(output_dir, paste0(basename(counts_path), ".upd.segments.tsv"))
write.table(atomic_df, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote atomic segment table to: ", output_path, "\n")

# --- 11. Merge and Summarize UPD Regions ---
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")

UPD_result <- atomic_df %>% 
  filter(UPD_final != "not_UPD") %>%
  mutate(
    chromosome = factor(chromosome, levels = chrom_order)
  ) %>%
  arrange(chromosome, start) %>%
  mutate(
    group = cumsum(
      lag(chromosome, default = as.character(chromosome[1])) != as.character(chromosome) |
        lag(UPD_final, default = UPD_final[1]) != UPD_final |
        lag(end, default = start[1]) + 1 != start
    )
  )

UPD_merged <- UPD_result %>%
  group_by(chromosome, UPD_final, group) %>%
  summarise(
    start = min(start),
    end = max(end),
    width = end - start + 1,
    .groups = "drop"
  ) %>%
  select(chromosome, start, end, UPD_final, width)%>%
  arrange(chromosome, start)

# --- 12. Remove Deletions in UPD Cadidates by Coverage Check ---
get_median_depth <- function(chr, start, end) {
  # Extract just the number from chr ("chr8" -> 8, "chrX" -> 23, "chrY" -> 24)
  chr_num <- if (grepl("^chrX$", chr)) {
    23
  } else {
    as.numeric(gsub("^chr", "", chr))
  }
  snps_in_region <- counts_df %>%
    filter(chromosome == chr_num, x >= start, x <= end)
  median(snps_in_region$depth, na.rm = TRUE)
}
get_median_depth_vec <- Vectorize(get_median_depth)

UPD_merged <- UPD_merged %>%
  mutate(median_depth = get_median_depth_vec(chromosome, start, end)) %>%
  select(chromosome, start, end, UPD_final, width, median_depth)

# Compute global median depth, excluding chrX
global_median_depth <- counts_df %>%
  filter(chrom != "chrX") %>%
  pull(depth) %>%
  median(na.rm = TRUE)

# Filter out any regions with median_depth less than half the global median
UPD_merged <- UPD_merged %>%
  filter(median_depth > 0.5 * global_median_depth)

# --- 13. Write Final UPD bed file ---
output_path = file.path(output_dir, paste0(basename(counts_path), ".upd.classification.bed"))
write.table(UPD_merged , file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote atomic segment table to: ", output_path, "\n")
