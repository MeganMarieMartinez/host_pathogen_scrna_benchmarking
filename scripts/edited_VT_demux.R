args <- commandArgs(trailingOnly = TRUE)

suppressMessages(library(Rsubread))

# Loading parameters
Parameter_file_path <- args[1]
Parameters <- read.table(Parameter_file_path, header = FALSE, sep = "\t")
Parameters <- as.character(Parameters$V1)
Parameters <- strsplit(Parameters, split = "=", fixed = TRUE)

Parameters_names  <- unlist(lapply(Parameters, function(x) x[1]))
Parameters_values <- unlist(lapply(Parameters, function(x) x[2]))
names(Parameters_values) <- Parameters_names

Output_directory      <- Parameters_values["Output_directory"]
Name_run              <- make.names(Parameters_values["Name_run"])
Viral_annotation_file <- as.character(Parameters_values["Viral_annotation_file"])

# Loading list of files to process
Parameter_target_file <- args[2]
File_to_process <- read.table(Parameter_target_file, header = FALSE, sep = "\t")
File_to_process <- unique(as.character(File_to_process$V1))

List_names <- c()
for (k in File_to_process) {
  if (file.exists(k)) {
    name_target <- basename(k)
    name_target <- sub("\\.gz$", "", name_target)
    name_target <- sub("\\.fastq$", "", name_target)
    List_names  <- c(List_names, name_target)
  }
}
List_target_path <- file.path(Output_directory, List_names)
names(List_target_path) <- List_names

# Loading all viral fragments detected in each sample
Identified_viral_fragments <- c()

for (k in List_target_path) {
  report <- read.table(file.path(k, "QC_filtered.txt"), header = TRUE)
  if (nrow(report) > 0) {
    Identified_viral_fragments <- c(Identified_viral_fragments, rownames(report))
  }
}

Identified_viral_fragments <- unique(Identified_viral_fragments)
Identified_viral_fragments <- Identified_viral_fragments[!is.na(Identified_viral_fragments)]

### Find the GTF; if none found, create a "pseudo GTF"
Path_GTF <- file.path(Output_directory, Name_run, "Merged_GTF.txt")

if (!file.exists(Path_GTF)) {
  cat("No GTF file found. Creating a basic GTF file\n")

  # Get the length of each reference from the BAM headers
  get_length <- function(segment) {
    for (p in List_target_path) {
      bam <- file.path(p, "Viral_BAM_files", paste0(segment, ".bam"))
      if (file.exists(bam)) {
        hdr <- system(paste("samtools view -H", shQuote(bam)), intern = TRUE)
        sq  <- hdr[grepl("^@SQ", hdr)]
        sn  <- sub(".*SN:([^\t]+).*", "\\1", sq)
        ln  <- as.numeric(sub(".*LN:([0-9]+).*", "\\1", sq))
        if (segment %in% sn) return(ln[match(segment, sn)])
      }
    }
    NA_real_
  }

  Length_segments <- vapply(Identified_viral_fragments, get_length, numeric(1))

  if (any(is.na(Length_segments))) {
    stop("Could not find a length for: ",
         paste(Identified_viral_fragments[is.na(Length_segments)], collapse = ", "))
  }

  gtf_lines <- data.frame(
    seqname   = Identified_viral_fragments,
    source    = "Empty_Annotation",
    feature   = "transcript",
    start     = 1L,
    end       = as.integer(Length_segments),
    score     = ".",
    strand    = "+",
    frame     = ".",
    attribute = paste0('gene_id "', Identified_viral_fragments, '_1";'),
    stringsAsFactors = FALSE
  )

  dir.create(dirname(Path_GTF), recursive = TRUE, showWarnings = FALSE)
  write.table(format(gtf_lines, scientific = FALSE, trim = TRUE),
              file = Path_GTF, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

### Counting
for (sample_name in names(List_target_path)) {
  path_temp <- List_target_path[[sample_name]]
  cat("Demultiplexing reads from sample", sample_name, "\n")

  # Aggregate reads from detected viruses
  List_bam_files <- file.path(path_temp, "Viral_BAM_files",
                              paste0(Identified_viral_fragments, ".bam"))
  List_bam_files <- List_bam_files[file.exists(List_bam_files)]
  if (length(List_bam_files) == 0) {
    cat("No BAM files found for", sample_name, "- skipping\n")
    next
  }
  List_bam_files <- paste0("'", List_bam_files, "'")

  merged_bam <- file.path(path_temp, "Reads_to_demultiplex.bam")
  command_merge <- paste("samtools merge", paste0("'", merged_bam, "'"), "-f",
                         paste(List_bam_files, collapse = " "))
  system(command_merge)

  # Assign reads to transcripts with featureCounts
  x <- suppressMessages(featureCounts(
    files = merged_bam,
    annot.ext = Path_GTF, isGTFAnnotationFile = TRUE,
    GTF.featureType = "transcript", GTF.attrType = "gene_id",
    reportReads = "BAM", reportReadsPath = path_temp,
    verbose = FALSE, primaryOnly = TRUE, allowMultiOverlap = TRUE))

  # Sort and index
  assigned_bam <- file.path(path_temp, "Reads_to_demultiplex.bam.featureCounts.bam")
  sorted_bam   <- file.path(path_temp, "Assigned_sorted.bam")
  system(paste("samtools sort", paste0("'", assigned_bam, "'"),
               "-o", paste0("'", sorted_bam, "'")))
  system(paste("samtools index", paste0("'", sorted_bam, "'")))

  # UMI-tools
  command_umi_tools <- paste(
    "umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell",
    "-I", paste0("'", sorted_bam, "'"),
    "-S", paste0("'", file.path(path_temp, "Expression_table.tsv"), "'"),
    "--wide-format-cell-counts")
  suppressMessages(system(command_umi_tools))
}
