# tumour_mrna_deg_module_lib.R

get_btx_palette <- function() {
	BTX_COL <-
		list(
			blue_d_in = c(189, 197, 255),
			blue_d_out = c(0, 62, 193),
			blue_l_in = c(179, 223, 255),
			blue_l_out = c(0, 151, 229),
			green_in = c(195, 249, 195),
			green_out = c(0, 169, 0),
			pink_in = c(254, 195, 252),
			pink_out = c(184, 56, 180),
			purple_in = c(212, 178, 255),
			purple_out = c(151, 73, 254),
			gray_in = c(214, 214, 214),
			gray_out = c(0, 0, 0)	
		)
	tmp <- rep(NA, length(BTX_COL))	
	for (i in 1:length(BTX_COL)) {
	  tmp[i] <- 
		rgb(
			red = (BTX_COL[[i]][1]), 
			green = (BTX_COL[[i]][2]), 
			blue = (BTX_COL[[i]][3]),
			maxColorValue = 255
		)
	}
	names(tmp) <- names(BTX_COL)
	BTX_COL <- as.list(tmp)	
	return(BTX_COL)

}


write_heatmap <- function(x, file_out) {
	svglite::svglite(file_out, width = 32.47, height = 13.1)
	ComplexHeatmap::draw(
		x, 
		heatmap_legend_side = "top", 
		annotation_legend_side = "top"
	)
	invisible(dev.off())
	return(file_out)
}

write_svg <- function(x, file_out, w = 32.47, h = 13.1) {
	svglite::svglite(file_out, width = w, height = h)
	plot(x)
	invisible(dev.off())
	return(file_out)
}


get_ensembl_ids <- function(x) {

	# Convert from gene.symbol to ensembl.gene
	logger::log_warn("Old version of ENSEMBL being used by get_ensembl_ids(); update to using HGNC API")
	library(EnsDb.Hsapiens.v86)
	output <- 
		ensembldb::select(
			EnsDb.Hsapiens.v86, 
			keys = x, 
			keytype = "SYMBOL", 
			columns = c("SYMBOL", "GENEID")
		) |>
		dplyr::distinct() |>
		dplyr::filter(grepl("ENSG", GENEID))

	# Test outupt
	if(all(x %in% output$SYMBOL)) {
		return(output)
	}else{
		logger::log_error("Missing ensembl ID.")
		stop()
	}
}	




#' Download TCGA/GTEx data from UCSC Xena Toil Hub for specified genes, phenotype, and survival.
#'
#' This function queries the UCSC Xena platform, specifically the Toil Hub
#' (\url{https://toil.xenahubs.net:443}), to retrieve:
#' \enumerate{
#'   \item Gene expression RNAseq data (expected counts) for a combined TCGA and GTEx dataset.
#'   \item Phenotype data for the TCGA/GTEx cohort.
#'   \item TCGA overall survival data.
#' }
#' It fetches expression values for the given gene symbols using `UCSCXenaTools::fetch_dense_values`
#' and transposes the resulting matrix so that samples are rows and genes are columns.
#'
#' This function aims to support workflows related to the combined TCGA and GTEx
#' expression datasets, such as the one described in:
#' "Network analysis of TCGA and GTEx gene expression datasets for
#' identification of trait-associated biomarkers in human cancer"
#' (\url{https://pmc.ncbi.nlm.nih.gov/articles/PMC8841814/}).
#'
	#' @param ensembl_id A string containing a ENSMEBL gene identifier.
#'
#' @return A tibble
#'
#' @importFrom checkmate assert_data_frame assert_names assert_character
#' @import UCSCXenaTools
#' @export
ucsc_toil_gene <- function(ensembl_id) {

	HOST <- "https://toil.xenahubs.net"
	## Gene expression data
	dataset <- "TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2" # RSEM expected_count (DESeq2 standardized)
	library(UCSCXenaTools)
	all_gene_identifiers <- .p_dataset_field(HOST, dataset) 
	ensembl_id_toil <-
		grep(paste(unique(ensembl_id$GENEID), collapse = "|"), all_gene_identifiers, value = TRUE)
	rsem_deseq_norm_count <-
		UCSCXenaTools::fetch_dense_values(
		  host = HOST,
		  dataset = dataset,
		  identifiers = ensembl_id_toil, 
		  use_probeMap = FALSE # change if using HUGO
		)
	# Transpose each fetched matrix so that samples are rows and genes are columns
	# This makes the data more commonly structured for downstream analysis.	
	rsem_deseq_norm_count <- t(rsem_deseq_norm_count) 
	# Simplify column names
	colnames(rsem_deseq_norm_count) <- 
		tools::file_path_sans_ext(colnames(rsem_deseq_norm_count))
	# Clean
	rsem_deseq_norm_count <- 
		rsem_deseq_norm_count |>
			as.data.frame() |>
			tibble::rownames_to_column(var = "sample")
	# Return fetched gene expression data
	return(rsem_deseq_norm_count = rsem_deseq_norm_count)
}


#' Merge Clinical Phenotype Data from UCSC Xena Toil Hub
#'
#' @description
#' This function downloads the clinical phenotype data for the TCGA TARGET GTEx
#' cohort from the UCSC Xena Toil hub. It then merges this clinical data with a
#' user-provided data frame based on a shared "sample" column.
#'
#' @details
#' The function automates the process of fetching and preparing the Toil hub
#' phenotype data using the `UCSCXenaTools` package. It performs the following
#' steps:
#' 1. Filters the `XenaData` to select the "phenotype" dataset from the "toilHub".
#' 2. Downloads the data using `XenaGenerate`, `XenaQuery`, and `XenaDownload`.
#' 3. Prepares the downloaded data object with `XenaPrepare`.
#' 4. Executes a left join to merge the clinical data into the input data frame `x`.
#'
#' @param x A data frame or tibble. This data frame must contain a column named
#'   `sample` whose values are sample identifiers matching those in the
#'   `TcgaTargetGTEX_phenotype` dataset (e.g., "TCGA-OR-A5J2-01A").
#'
#' @return A new data frame that is the result of a left join between the input
#'   data frame `x` and the downloaded clinical phenotype data. Columns from the
#'   phenotype data are appended to the columns of `x`.
#'
#' @importFrom dplyr filter left_join 
#' @importFrom UCSCXenaTools XenaData XenaGenerate XenaQuery XenaDownload XenaPrepare
#'
#' @examples
#' \dontrun{
#' # Ensure you have the required packages installed:
#' # install.packages(c("dplyr", "UCSCXenaTools"))
#'
#' # Create a sample data frame with sample identifiers
#' my_sample_data <- data.frame(
#'   sample = c("TCGA-OR-A5J2-01A", "TCGA-OR-A5J3-01A"),
#'   gene_expression = c(12.3, 8.5),
#'   mutation_status = c("Positive", "Negative")
#' )
#'
#' # Merge the clinical phenotype data
#' clinical_merged_data <- ucsc_toil_pheno(my_sample_data)
#'
#' # View the first few columns of the merged data
#' head(clinical_merged_data)
#' }
#'
ucsc_toil_pheno <- function(x, y) {

	clinical_data <-
		UCSCXenaTools::XenaData |>
		dplyr::filter(XenaHostNames == "toilHub") |>
		dplyr::filter(DataSubtype == "phenotype") |>  # select clinical dataset
		dplyr::filter(XenaDatasets == "TcgaTargetGTEX_phenotype.txt") |>  # select dataset
		UCSCXenaTools::XenaGenerate() |>  # generate a XenaHub object
		UCSCXenaTools::XenaQuery() |>
		UCSCXenaTools::XenaDownload()
	clinical_data <- 
		UCSCXenaTools::XenaPrepare(clinical_data) |>
		janitor::clean_names() |>
		dplyr::mutate(primary_site = iconv(primary_site, from = "ISO-8859-1", to = "UTF-8")) |>
		dplyr::mutate(sample_type = iconv(sample_type, from = "ISO-8859-1", to = "UTF-8")) 
	# Annotate Core Indications of Interest
	key_indication <- c("Lung", "Esophagus", "Stomach", "Colon", "Head and Neck region", "Pancreas", "Breast", "Bladder", "Ovary", "Prostate", "Liver", "Uterus", "Soft tissue,Bone")	
	clinical_data <-
		clinical_data |>
		dplyr::mutate(priority = ifelse(primary_site %in% key_indication, "Core", "Ancillary")) |>
		dplyr::mutate(priority = factor(priority, levels = c("Core", "Ancillary")))
	# Rename adjacent normals in TCGA (as a group negative control)
	clinical_data <-
		clinical_data |>
		dplyr::mutate(primary_disease_or_tissue = ifelse(sample_type == "Solid Tissue Normal", "Negative control", primary_disease_or_tissue)) 	
	# Augment with BRCA annotations 
	tmp <-
		y |>
		dplyr::mutate(CLID = stringr::str_sub(CLID, end = -2))
	clinical_data <-
		clinical_data |>
		dplyr::left_join(tmp, by = c("sample" = "CLID")) 
	# Rename BRCA ER/HR/TNBC subtypes
	clinical_data <-
		clinical_data |>
		dplyr::mutate(brca = ifelse(`Triple Negative Status` == "Yes", "TNBC", NA)) |>
		dplyr::mutate(brca = ifelse(`er_status_by_ihc` == "Positive" & `HER2.newly.derived` == "Negative", "BRCA_ER+HER2-", brca)) |>
		dplyr::mutate(brca = ifelse(`er_status_by_ihc` == "Negative" & `HER2.newly.derived` == "Positive", "BRCA_ER-HER2+", brca)) |>
		dplyr::mutate(brca = ifelse(`er_status_by_ihc` == "Positive" & `HER2.newly.derived` == "Positive", "BRCA_ER+HER2+", brca)) 
	# Add purity (where possible)
	# TODO: future work, complete purity estimates for all cohorts
	tmp_purity <-
		UCSCXenaShiny::load_data("tcga_purity") |>
			dplyr::select(sample, CPE) |>
			dplyr::group_by(sample) |>
			dplyr::slice_max(n = 1, order_by = CPE)
	clinical_data <-
		clinical_data |>
		dplyr::left_join(tmp_purity, by = "sample")
	# Clean outputs
	v <- grep("ENS", colnames(x), value = TRUE)
	out <-
		dplyr::left_join(x, clinical_data, by = "sample") |>
		dplyr::filter(!grepl("Leukemia", primary_disease_or_tissue, ignore.case = TRUE)) |>
		dplyr::filter(!grepl("Lymphoma", primary_disease_or_tissue, ignore.case = TRUE)) |>
		dplyr::filter(primary_disease_or_tissue != "Negative control") |> # Pending PCA
		dplyr::filter(sample_type != "Cell Line") |>
		dplyr::mutate(primary_disease_or_tissue = ifelse(is.na(brca), primary_disease_or_tissue, brca)) |>
		dplyr::select(all_of(c("study", "sample", "priority", "primary_disease_or_tissue", "primary_site", "sample_type", "CPE", v))) |>
		tibble::set_tidy_names(syntactic = TRUE) |>
		tibble::as_tibble() |>
		dplyr::distinct()

		#dplyr::filter(sample_type != "Solid Tissue Normal") # remove adjacent normal tissues

	return(out)
}


get_n <- function(x) {
	out <- 
		data.frame(
			y = 0,
			label = sprintf("(n=%s)", scales::number(length(x)))
		)
   return(out)
}

plot_expression_ridge <- function(x, keys, target, plt_tumour = TRUE) {

	## --- Scaffolding ---
	# tar_load(exp_gene_pheno_tb2); x = exp_gene_pheno_tb2; keys = all_ensgids; target = "ERBB2"; plt_tumour = TRUE


	## --- Retrieve transcript identifier ---
	bait <-
		keys |>
		dplyr::filter(SYMBOL == target) |>
		dplyr::select(GENEID) |>
		dplyr::distinct() |>
		unlist() |>
		paste(collapse = "|")
	v <- grep(bait, colnames(x), value = TRUE)

	# Remove missing data
	x <- tidyr::drop_na(x, all_of(v))
		
	summary_stats_study <- 
		x |>
		dplyr::group_by(study) |>
		dplyr::summarise(med = median(.data[[v]], na.rm = TRUE), n = dplyr::n())
	annotation_n <-
		x |>
		dplyr::filter(study != "GTEX") |>	
		dplyr::group_by(primary_disease_or_tissue) |>
		dplyr::count()
	med_line <- 
		summary_stats_study |>
		dplyr::filter(study == "GTEX") |>
		dplyr::select("med") |>
		as.numeric()

	summary_stats_psite <- 
		x |>
		dplyr::group_by(study, primary_disease_or_tissue) |>
		dplyr::summarise(med = median(.data[[v]], na.rm = TRUE), n = dplyr::n()) |>
		dplyr::arrange(desc(med))
	tissue_order <- summary_stats_psite$primary_disease_or_tissue
	x_lims <- round(max(x[[v]], na.rm = TRUE) * 1.05)

	if (plt_tumour) {
		to_plot <- 
			x |>
			dplyr::filter(study %in% c("TCGA", "TARGET")) |>
			dplyr::filter(sample_type != "Normal Tissue")  # filter to primary/recur/met
	}else{
		to_plot <- 
			x |>
			dplyr::filter(study %in% c("GTEX")) |>
			dplyr::filter(sample_type == "Normal Tissue") # filter to healthy	
	}

	# Formatting
	btx_cols <- 
		list(
			turq = c(0, 115, 172),
			p1 = c(216, 136, 213),
			p2 = c(229, 176, 227),
			p3 = c(242, 215, 241)		
		)
	tmp <- rep(NA, length(btx_cols))	
	for (i in 1:length(btx_cols)) {
	  tmp[i] <- 
		rgb(
			red = (btx_cols[[i]][1]), 
			green = (btx_cols[[i]][2]), 
			blue = (btx_cols[[i]][3]),
			maxColorValue = 255
		)
	}
	btx_cols <- as.list(rev(tmp))
	font_family <- "lato"
	
	# Plotting
	ridge_plt <- 
		to_plot |>
		dplyr::arrange(desc(median(.data[[v]], na.rm = TRUE))) |>
		dplyr::mutate(
			primary_disease_or_tissue = factor(primary_disease_or_tissue, levels = tissue_order)
		) |>
		ggplot(aes(x = primary_disease_or_tissue, y = .data[[v]])) +
		ggdist::stat_halfeye(fill_type = "segments", alpha = 0.3) +
		ggdist::stat_interval() +
		stat_summary(geom = "point", fun = median) +
		stat_summary(
			fun.data = get_n, 
			geom = "text", 
			family = font_family, 
			size = 3, color = "grey75") +
		geom_hline(yintercept = med_line, col = "grey30", lty = "dashed") +
		annotate("text", x = 0.5, y = med_line, label = "Median healthy expression",
			   family = font_family, size = 3, hjust = 0) +
		scale_x_discrete(labels = toupper) +
		scale_y_continuous(limits = c(0, x_lims), breaks = seq(0, 50, 1)) +
		#paletteer::scale_color_paletteer_d("MetBrewer::Hokusai1", 1) +
		#paletteer::scale_color_paletteer_d("MetBrewer::VanGogh3", 1) +
		scale_color_manual(values = btx_cols) +
		coord_flip(clip = "off") +
		guides(col = "none") +
		facet_wrap(vars(priority), scales = "free_y") +
		  theme_minimal(base_family = font_family) +
		  theme(
			plot.background = element_rect(color = NA, fill = "grey97"),
			panel.grid = element_blank(),
			panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
			plot.title.position = "plot",
			plot.subtitle = ggtext::element_textbox_simple(
			  margin = margin(t = 4, b = 16), size = 10),
			plot.caption = ggtext::element_textbox_simple(
			  margin = margin(t = 12), size = 7
			),
			plot.caption.position = "plot",
			axis.text.y = element_text(hjust = 0, margin = margin(r = -10)),
			plot.margin = margin(4, 4, 4, 4)
	  ) 
	#stat_n_text(family = font_family, size = 3, color = "grey75") # stopped working
	
	# Sanity Check: all tissues present?
	test <- toupper(unique(to_plot$primary_disease_or_tissue))
	y_labs <- 
		c(
			ggplot_build(ridge_plt)$layout$panel_params[[1]]$y$get_labels(),
			ggplot_build(ridge_plt)$layout$panel_params[[2]]$y$get_labels()
		)
	if(!all(test %in% y_labs)){
		msg <- "Missing tissues in plot: plot_expression_ridge()"
		logger::log_error(msg)
		stop(msg)
	}
	
	return(ridge_plt)
}




get_comparison_plan <- function(x) {
	## --- Scaffolding ---
	# tar_load(exp_gene_pheno_tb); x = exp_gene_pheno_tb; tar_load(ebrt_tb); ebrt = ebrt_tb

	## --- Dependencies ---
	library(stringdist)
	library(readr)	

	## --- Parse inputs ---
	indications <- 
		x |>
		dplyr::filter(study != "GTEX") |>
		dplyr::distinct(`primary_disease_or_tissue`) |>
		unlist() |>
		sort()
	
	normal_tissue <- 
		x |>
		dplyr::filter(study == "GTEX") |> 
		dplyr::distinct(`primary_disease_or_tissue`) |>
		dplyr::pull(`primary_disease_or_tissue`)
	test_tissue <-
		grep(
			c("adrenal|heart|kidney|liver|lung|salivary|pancreas|pituitary|intestine"),
			normal_tissue,
			ignore.case = TRUE,
			value = TRUE
		)

	## For unknown encoding issues, this is not producing matching strings: Gemini errors!
	# ebrt_tissues <-
		# ebrt |>
		# dplyr::mutate(GTEX_xena = iconv(GTEX_xena, from = "latin1", to = "UTF-8", sub = "")) |>
		# tidyr::separate_rows(GTEX_xena, sep = "\\|") |>
		# dplyr::filter(EBRT_bicycle <= 40) |>
		# dplyr::filter(GTEX_xena != "N/A") |>
		# dplyr::filter(GTEX_xena != "Testis") |>
		# dplyr::filter(GTEX_xena != "Ovary") |>
		# dplyr::filter(GTEX_xena != "Fallopian tube") |>
		# dplyr::pull(GTEX_xena) |>
		# unlist() |>
		# sort()		
	# match(normal_tissue, ebrt_tissues)

	targets <- grep("ENS", colnames(x), value = TRUE)

	## --- Use tidyr::crossing to get all combinations ---
	out <- tidyr::crossing(indications, targets, test_tissue)

	return(out)
}


get_iqr <- function(x) {
	#' Calculate Median and IQR for Gene Expression by Group
	#'
	#' @description
	#' This function takes a data frame of expression data and calculates the median
	#' and Interquartile Range (IQR) for all gene columns (identified by containing
	#' "ENS" in their names). The calculations are performed for each group defined
	#' by the combination of `primary_disease_or_tissue` and `study`.
	#'
	#' @details
	#' The function follows a "group-summarise-tidy" workflow:
	#' 1.  First, it cleans column names and removes rows with any `NA` values.
	#' 2.  It then groups the data frame by the `primary_disease_or_tissue` and `study` columns.
	#' 3.  Using `dplyr::summarise` with `dplyr::across`, it computes the median and IQR
	#'     for every column whose name contains "ENS". `NA` values are ignored in
	#'     these calculations.
	#' 4.  Finally, it reshapes the resulting wide-format summary table into a tidy,
	#'     long-format tibble. The output has one row for each combination of
	#'     tissue, study, and gene, with dedicated columns for the `median` and `iqr`.
	#'
	#' @param x A data frame or tibble. It must contain the columns
	#'   `primary_disease_or_tissue`, `study`, and at least one column with a name
	#'   containing the string "ENS".
	#'
	#' @return A tidy tibble with the following columns:
	#'   \itemize{
	#'     \item `primary_disease_or_tissue`: The grouping tissue/disease type.
	#'     \item `study`: The grouping study (e.g., "TCGA", "GTEX").
	#'     \item `gene`: The name of the gene column (e.g., "ENSG001").
	#'     \item `median`: The calculated median expression for that gene within that group.
	#'     \item `iqr`: The calculated Interquartile Range for that gene within that group.
	#'   }
	#'
	#' @importFrom tibble set_tidy_names
	#' @importFrom tidyr drop_na pivot_longer pivot_wider
	#' @importFrom dplyr group_by summarise across contains
	#'
	#' @examples
	#' # Create a sample data frame for demonstration
	#' mock_expression_data <- data.frame(
	#'   primary_disease_or_tissue = c("Lung Cancer", "Lung Cancer", "Normal Lung", "Normal Lung"),
	#'   study = c("TCGA", "TCGA", "GTEX", "GTEX"),
	#'   ENSG001 = c(10, 12, 105, 108),
	#'   ENSG002 = c(50, 55, 20, 22),
	#'   some_other_col = 1:4 # This column will be ignored
	#' )
	#'
	#' # Calculate the grouped statistics
	#' summary_stats <- get_iqr(mock_expression_data)
	#'
	#' # View the tidy output
	#' print(summary_stats)
	#' #> # A tibble: 4 × 5
	#' #> # Groups:   primary_disease_or_tissue, study [2]
	#' #>   primary_disease_or_tissue study gene    median   iqr
	#' #>   <chr>                     <chr>   <chr>    <dbl> <dbl>
	#' #> 1 Lung Cancer               TCGA    ENSG001   11       1
	#' #> 2 Lung Cancer               TCGA    ENSG002   52.5     2.5
	#' #> 3 Normal Lung               GTEX    ENSG001  106.      1.5
	#' #> 4 Normal Lung               GTEX    ENSG002   21       1

	# Reference: all normal tissues ------------------------------------------
	## Split data into tumour/normal
	to_compare <- x 

	out_tb <- 
		to_compare |>
			#dplyr::filter(`primary_disease_or_tissue` == indication & study != "GTEX") |>
		dplyr::group_by(`primary_disease_or_tissue`, study) |>	
		# 1. Summarise across all columns that contain "ENS"
		dplyr::summarise(
			dplyr::across(.cols = contains("ENS"),
			# 2. Apply a named list of functions to each column
			#  We use purrr-style anonymous functions `~` to add `na.rm = TRUE`
			.fns = list(
				q25 = ~quantile(.x, probs = 0.25, na.rm = TRUE),
				q50 = ~quantile(.x, probs = 0.5, na.rm = TRUE),
				q75 = ~quantile(.x, probs = 0.75, na.rm = TRUE)
			),
		# 3. Define the naming pattern for the new columns
		.names = "{.col}_{.fn}"
			),
		.groups = "drop"
		) |>
	# 4. Reshape the data from wide to long format
	tidyr::pivot_longer(
		cols = contains("ENS"),
		names_to = c("gene", "statistic"),
		names_pattern = "(.*)_(q25|q50|q75)", # Regex to separate gene name and statistic
		values_to = "value"
	) |>
	# 5. Reshape again to give each statistic its own column (the final tidy format)
	tidyr::pivot_wider(
	  names_from = statistic,
	  values_from = value
	)
	return(out_tb)
}

get_reference_dist <- function(x) {
	# Given an annotated input expression matrix, find the normal tissue with highest gene expression distribution. 
	out <-
		x |>
			dplyr::filter(study == "GTEX") |>
			dplyr::group_by(gene) |>
			dplyr::slice_max(n = 1, q50, with_ties = FALSE) |>
			dplyr::select(gene, primary_disease_or_tissue)
	return(out)
}

#' Calculate Earth Mover's Distance (EMD) for a Gene Between Two Tissues
#'
#' @description
#' This function compares the expression distribution of a single gene between a
#' "tumour" tissue and a "normal" reference tissue. The primary comparison is
#' done using the Earth Mover's Distance (EMD) metric from the `EMDomics`
#' package.
#'
#' Additionally, it calculates several supplementary statistics:
#' - A Wilcoxon rank-sum test p-value (two sided, for simplicity).
#' - The difference in median expression between the two groups.
#' - The skewness of the expression values in the tumour cohort.
#'
#' The function is designed to work with a specific data structure where one
#' cohort (e.g., TCGA) is compared against a reference cohort (e.g., GTEx).
#'
#' @param x A tibble or data frame containing gene expression data. It must
#'   include columns for `study`, `sample`, `primary_disease_or_tissue`, and
#'   at least one column for gene expression (e.g., "ENSG...").
#' @param pair_indices A character vector of length 3 containing the identifiers
#'   for the comparison:
#'   1.  `character`: The name of the tumour tissue (e.g., "Liver Hepatocellular Carcinoma").
#'   2.  `character`: The ENSEMBL ID of the gene to analyze (e.g., "ENSG00000123838").
#'   3.  `character`: The name of the normal reference tissue from the GTEx study
#'       (e.g., "Heart - Atrial Appendage").
#' @param perms An integer specifying the number of permutations to use for the
#'   EMD p-value calculation. Defaults to `1`.
#'
#' @return
#' A tibble (a 1-row data frame) containing the aggregated results, including:
#' - `GENEID`: The gene ID.
#' - `emd`: The EMD score.
#' - `p.value`: The p-value from the EMD calculation.
#' - `d`: The difference in median expression (tumour - normal).
#' - `skew_tumour`: The skewness of the tumour samples' expression.
#' - `pval_wilcox`: The p-value from the Wilcoxon test.
#' - `primary_disease_or_tissue`: The name of the tumour tissue.
#' - `ref_tissue`: The name of the normal reference tissue.
#'
#' The function returns `NULL` if the gene has a median expression below 0.5
#' in the combined dataset or if the EMD calculation fails.
#'
#' @importFrom dplyr filter select bind_rows summarise mutate group_by across left_join distinct
#' @importFrom tidyr drop_na pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom EMDomics calculate_emd
#' @importFrom e1071 skewness
#' @importFrom broom tidy
get_emd <- function(x, pair_indices, perms = 1) {

	## --- Scaffolding ---
	# tar_load(exp_gene_pheno_tb); x = exp_gene_pheno_tb; perms = 1; use_max_ref = TRUE; pair_indices = c("Liver Hepatocellular Carcinoma", "ENSG00000123838", "Heart - Atrial Appendage")

	## --- Dependencies ---
	library(dplyr)

	## --- Reference: all normal tissues ---
	## Split data into tumour/normal
	tumour_tb <- 
		x |> 
		dplyr::filter(primary_disease_or_tissue == pair_indices[1]) |>
		dplyr::select(study, sample, primary_disease_or_tissue, contains(pair_indices[2]))
	
	normal_tb <- 
		x |>
		dplyr::filter(study == "GTEX") |> 
		dplyr::filter(primary_disease_or_tissue == pair_indices[3]) |>		
		dplyr::select(study, sample, primary_disease_or_tissue, contains(pair_indices[2]))			


	## --- Sanity check ---
	if(FALSE){
		test_plt <- 
			dplyr::bind_rows(normal_tb, tumour_tb) |>
			ggplot(aes(x = study, y = !!sym(pair_indices[2]))) +
			geom_boxplot() + 
			ggpubr::stat_compare_means()
	}
	
	## --- Calculate wilcoxon test 
	wilcox_test <- 
		dplyr::bind_rows(normal_tb, tumour_tb) |> 
		dplyr::summarise(
			broom::tidy(wilcox.test(!!sym(pair_indices[2]) ~ study))
		)	
	
	
	## --- Calculate EMD ---			
	test_me <- 
		rbind(tumour_tb, normal_tb) |>
		tidyr::drop_na()
	exp_data <- 
		test_me |>
			dplyr::select(contains("ENS")) |>	
			as.matrix()
	rownames(exp_data) <-
		test_me |>
			dplyr::select(sample) |>
			unlist() |>
			as.vector()	
	## --- Skip unexpressed genes (in general) ---
	filter_in <- 
		apply(X = exp_data, FUN = function(a) median(a) >= 0.5, MARGIN = 2)
	if(!filter_in){return(NULL)}
	# Label by GTEx vs TCGA/TARGET
	lbls <- 
		test_me |>
			dplyr::mutate(study = as.character(study != "GTEX")) |>
			dplyr::select(study) |>
			unlist() |>		
			as.vector()
	names(lbls) <- rownames(exp_data)
	## --- Calculate Earth Mover Distance ---
	# TODO: throwing non finite errors when using use_max_ref?
	tmp <- t(exp_data)
	results <- 
		tryCatch({
			EMDomics::calculate_emd(
				rbind(tmp, tmp), 
				lbls, 
				#binSize = 0.2,
				nperm = perms, # amend in production
				seq = FALSE, # data is already in log2 space
				verbose = FALSE,
				parallel = FALSE # amend in production?
			)
		}, error = function(e) {
			logger::log_error("Failed EMD:", e$message)
			return(NULL) # Return NULL if reading fails
		})
	if(is.null(results)){return(NULL)}
	
	## --- Calculate relative direction of change ---
	median_diff <-
		test_me |>
			dplyr::group_by(study) |>
			dplyr::select(study, contains("ENS")) |>	
			dplyr::summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop") |>
			tidyr::pivot_longer(!study, names_to = "GENEID", values_to = "median") |>
			dplyr::group_by(GENEID) |>
			dplyr::summarise(d = diff(median), .groups = "drop")
	## --- Calculate skewness ---
	skew_tumour <-
		tumour_tb |>
			dplyr::select(contains("ENS")) |>	
			dplyr::summarise(across(where(is.numeric), ~ e1071::skewness(.x, na.rm = TRUE))) |> 
			unlist()
	## --- Aggregate outputs ---
	out_tb <- 
		results$emd |>
		as.data.frame() |>
		dplyr::slice(1) |>
		tibble::rownames_to_column(var = "GENEID") |>
		dplyr::left_join(median_diff, by = "GENEID") |>
		dplyr::mutate(skew_tumour = skew_tumour) |>
		dplyr::mutate(primary_disease_or_tissue = pair_indices[1]) |>
		dplyr::mutate(ref_tissue = pair_indices[3]) |>
		dplyr::mutate(pval_wilcox = wilcox_test$p.value)

	# Re-add annotations
	tmp <-
		x |>
		dplyr::select(study, primary_disease_or_tissue,  priority) |>
		dplyr::distinct()	 
	out_tb <-
		out_tb |>
		dplyr::left_join(tmp, by ="primary_disease_or_tissue")

	return(out_tb)
}


plot_static_heatmap <- function(x, y) {

	## --- Scafoolding ---
	# tar_load(gene_emd_max); x = gene_emd_max; tar_load(ensembl_ids); y = ensembl_ids; 

	## --- Dependencies ---
	library(ComplexHeatmap)
	library(dplyr)

	## --- Prepare data ---
	# Lowest EMD value for earch tissue-tumour comparison
	to_plot <-
		x |>
		dplyr::left_join(y, by = c("GENEID")) |>
		dplyr::filter(d >= 0) |>
		dplyr::filter(!grepl("Leukemia", primary_disease_or_tissue)) |>
		dplyr::filter(!grepl("Lymphoma", primary_disease_or_tissue)) |>
		dplyr::select(priority, primary_disease_or_tissue, SYMBOL, emd, pval_wilcox) |>
		dplyr::ungroup() |>
		dplyr::group_by(primary_disease_or_tissue, SYMBOL) |>
		dplyr::slice_min(n = 1, order_by = emd) 
	
	## --- Meta-analysis of Wilcoxon p-values ---
	# Fisher Method: Best for asymmetric scenarios where you suspect at least one of your tests has a strong effect, even if others do not. 
	# to_test <- split(x = to_plot$pval_wilcox, f = to_plot$SYMBOL)
	# tmp <- lapply(X = to_test, FUN = metap::sumlog, log.p = FALSE)
	# tmp1 <- unlist(sapply(X = tmp, FUN = "[", i = "p")) 
	# meta_tb <- 
		# tibble::tibble(
			# SYMBOL = gsub("\\.p$", "", names(tmp1)),
			# pval = -log10(p.adjust(tmp1, method = "bonferroni"))
		# ) 
	# idx <- is.infinite(meta_tb$pval)
	# top_p <- max(meta_tb$pval[!idx])	
	# meta_tb$pval[idx] <- top_p	# winsorise Inf values

	## --- Extract matrix ---
	tmp <- 
		to_plot |>
		dplyr::select(-priority, -pval_wilcox) |>
		tidyr::pivot_wider(names_from = SYMBOL, values_from = emd) |>
		dplyr::ungroup() 
	m <- 
		tmp |>
		dplyr::select(-primary_disease_or_tissue) |>
		as.matrix()
	rownames(m) <- tmp$primary_disease_or_tissue
	m[is.na(m)] <- 0

	## --- Colours ---
	BTX_COL <- get_btx_palette()
	tealrose_heatmap <- 
		as.character(paletteer::paletteer_c(palette = "grDevices::TealRose", direction = -1, n = 100))
	# breaks_1 <- seq(0, round(max(m)), length.out = length(tealrose_heatmap)) 
	breaks_1 <- quantile(m, c(1:100/100), na.rm = TRUE) # histogram-weighted, like boxplot	
	tealrose_heatmap <- 
		circlize::colorRamp2(breaks = breaks_1, colors = tealrose_heatmap)

	# pval_heatmap <- 
		# as.character(paletteer::paletteer_c(palette = "grDevices::Mako", direction = -1, n = 100))	
	# breaks_1 <- quantile(meta_tb$pval, c(1:100/100), na.rm = TRUE)	
	# pval_heatmap <- 
		# circlize::colorRamp2(breaks = breaks_1, colors = pval_heatmap)		
	
	## --- Bottom target annotation ---
	# pval_col <- meta_tb$pval[match(colnames(m), meta_tb$SYMBOL)]
	# bottom_ha <- 
		# ComplexHeatmap::HeatmapAnnotation(
			# fisher_p = pval_col,
			# #col = list(fisher_p = pval_heatmap),
			# annotation_name_side = "left",
			# annotation_legend_param = 
				# list(
					# fisher_p = list(
						# title = "Meta p-value (-log10)", 
						# direction = "horizontal", 
						# nrow = 1)
				# )
		# )

	## Indication annotation
	grps_row <- to_plot$priority[match(rownames(m), to_plot$primary_disease_or_tissue)]
	row_ha <- 
		ComplexHeatmap::rowAnnotation(
			priority = grps_row, 
			col = list(priority = c("Core" = BTX_COL$blue_d_out, "Ancillary" = BTX_COL$gray_in))
		)
	
	a_plt <-
		ComplexHeatmap::Heatmap(
			m, 
			col = tealrose_heatmap,
			#row_km = 6, 
			#column_km = 6, 
			name = "EMD (positive only)",
			heatmap_legend_param = list(
				title = "EMD (positive only)",  
				direction = "horizontal"
			),			
			na_col = "black",
			right_annotation = row_ha,
			# bottom_annotation = bottom_ha,
			row_title = NULL,
			column_title = NULL,
			row_names_gp = gpar(fontsize = 18),
			row_names_max_width = unit(20, "cm"), 			
			column_names_rot = 45,
			column_names_gp = gpar(fontsize = 18, hjust = 1),
			rect_gp = gpar(col = "white", lwd = 1),
			show_row_names = TRUE,
			show_column_names = TRUE
		)
			
	return(a_plt)
}




ucsc_met500_download <- function(annotations) {

	# Integrative Clinical Genomics of Metastatic Cancer - PMC
	## https://pmc.ncbi.nlm.nih.gov/articles/PMC5995337/

	## Catalog
	# UCSCXenaTools::XenaData |> View()

	## --- Scaffolding ---
	# tar_load(all_ensgids); annotations = all_ensgids;

	## --- Checkmate ---
	checkmate::assert_data_frame(
		annotations, 
		min.rows = 1, 
		.var.name = "annotations"
	)
	checkmate::assert_names(
		names(annotations), 
		must.include = "SYMBOL", 
		.var.name = "annotations"
	)
	checkmate::assert_character(
		annotations$SYMBOL,
		min.len = 1,          # Ensure gene symbols are not empty strings
		any.missing = FALSE,  # Ensure no missing (NA) gene symbols
		unique = FALSE,        
		.var.name = "annotations$SYMBOL"
	)


	## --- Dependencies ---
	library(EnsDb.Hsapiens.v86)
	HOST <- "https://ucscpublic.xenahubs.net"

	## --- Get Probe map values from Xena ---
	xe <- UCSCXenaTools::XenaGenerate(subset = XenaCohorts == "MET500 (expression centric)")
	xe_query <- UCSCXenaTools::XenaQueryProbeMap(xe)
	xe_download <- UCSCXenaTools::XenaDownload(xe_query)
	probe_map <- UCSCXenaTools::XenaPrepare(xe_download)	
	
	## --- Get query ENSGIDs ---
	query_ensgid <-
		annotations |>
		dplyr::filter(grepl("ENSG", GENEID)) |>
		dplyr::select(GENEID) |>
		dplyr::distinct() |>
		unlist() 
	
	## --- Lookup ENSGIDs ---
	ensgids <-
		query_ensgid |>
		stringr::str_c(collapse = "|")
	idx <-
		grepl(
			ensgids,
			probe_map[["id"]]
		)
	lookup_ids <-
		probe_map[idx, ] |>
		dplyr::select(id) |>
		unlist()

	## --- Gene expression data ---
	# https://xenabrowser.net/datapages/?dataset=MET500%2FgeneExpression%2FM.mx.log2.txt&host=https%3A%2F%2Fucscpublic.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
	ge <-
		UCSCXenaTools::fetch_dense_values(
		  host = HOST,
		  dataset = "MET500/geneExpression/M.mx.log2.txt",
		  identifiers = lookup_ids, 
		  use_probeMap = FALSE # Use full ensembl ids
		)
	# Transpose each fetched matrix so that samples are rows and genes are columns
	# This makes the data more commonly structured for downstream analysis.	
	ge <-
		t(ge) 	
	# Simplify column names
	colnames(ge) <- tools::file_path_sans_ext(colnames(ge))
	
	# Retrieve clinical data	
	clinical_data <-
		UCSCXenaTools::XenaData |>
		dplyr::filter(XenaHostNames == "publicHub") |>
		dplyr::filter(DataSubtype == "phenotype") |>  # select clinical dataset
		dplyr::filter(grepl("MET500", XenaCohorts)) |>
		UCSCXenaTools::XenaGenerate() |>  # generate a XenaHub object
		UCSCXenaTools::XenaQuery() |>
		UCSCXenaTools::XenaDownload()
	clinical_data = UCSCXenaTools::XenaPrepare(clinical_data)
	
	# Return a list containing the metadata of the cohorts and the fetched gene expression data
	return(list(clinical_data = clinical_data, ge = ge))
}


plot_met500 <- function(x, y, id, annotation) {

	## --- Scaffolding ---
	# tar_load(met500_tb); x = met500_tb; tar_load(all_ensgids); y = all_ensgids;  id =  "ST14"

	## --- Dependencies ---
	library(ggplot2)
	library(dplyr)

	met500_ge <- 
		x$ge |>
		tibble::as_tibble(rownames = "sample") # genes
	met500_phenotype <-
		x$clinical_data |>
		janitor::clean_names() |>
		dplyr::mutate(tissue = iconv(tissue, from = "ISO-8859-1", to = "UTF-8")) |>
		dplyr::filter(!grepl("blood", tissue, ignore.case = TRUE)) |>
		dplyr::filter(!grepl("cell line", tissue, ignore.case = TRUE)) |>
		tidyr::drop_na(tissue) |>
		dplyr::rename("sample" = "sample_id")
		
	key_indication <- tolower(c("Lung", "Esophagus", "Stomach", "Colon", "Oral", "Pancreas", "Breast", "Bladder", "Ovary", "Prostate", "Liver", "Uterus", "Soft tissue,Bone"))	

	met500_merged <- 	
		met500_phenotype |>
		dplyr::left_join(met500_ge, by = "sample") |>
		dplyr::mutate(priority = ifelse(tissue %in% key_indication, "Core", "Ancillary")) |>
		dplyr::mutate(priority = factor(priority, levels = c("Core", "Ancillary"))) 
	## --- Median and MAD and sample size ---
	# Ideal is high median, low MAD
	to_plot <-
		met500_merged |>
		tidyr::pivot_longer(contains("ENSG"), names_to = "ENSGID", values_to = "tpm") |>
		dplyr::group_by(cohort, ENSGID) |>
		dplyr::summarise(
			n = dplyr::n(),
			tpm_median = median(tpm),
			tpm_mad = mad(tpm),
			prop_expressed = sum(tpm > 0) / length(tpm),
			.groups = "drop_last"
		) |>
		#dplyr::mutate(scr = tpm_median * (1 - tpm_mad))|>
		dplyr::ungroup() |>
		dplyr::arrange(desc(tpm_median), tpm_mad) 


	## --- Update SYMBOL ---
	library(EnsDb.Hsapiens.v86)
	annot <-
		ensembldb::select(
			EnsDb.Hsapiens.v86, 
			keys = unique(to_plot$ENSGID), 
			keytype = "GENEID", 
			columns = c("GENEID", "SYMBOL")
		) |>
		dplyr::rename("ENSGID" = "GENEID")
	to_plot <-
		to_plot |>
		dplyr::left_join(annot, by = "ENSGID") |>
		dplyr::mutate(scr = row_number()) |>
		dplyr::mutate(lbl = ifelse(scr <= 20, SYMBOL, "")) 
	
	## --- Colours ------------------------------------------------------------
	BTX_COL <- get_btx_palette()
	library(sysfonts)
	sysfonts::font_add_google("Lato", "lato")
	sysfonts::font_add_google("Open Sans", "opensans")
	showtext::showtext_auto()
	library(ggplot2)
	old_theme <- 
	  ggplot2::theme_set(ggplot2::theme_minimal()) +
	  ggplot2::theme_update(
		# legend.title = element_blank(),
		# legend.justification = c(0, 1), 
		# legend.position = c(.1, 1.075),
		legend.background = element_blank(),
		axis.title= element_text(family = "opensans", size = 10),
		plot.title = element_text(family = "lato", size = 20, margin = margin(b = 10)),
		plot.subtitle = element_text(family = "opensans", size = 10, color = "darkslategrey", margin = margin(b = 25)),
		plot.caption = element_text(family = "opensans", size = 8, margin = margin(t = 10), color = "grey70", hjust = 0)
	  )

	# palettes_d_names |>
	# dplyr::filter(length >= length(unique(to_plot$cohort)))

	## --- Unexpressed genes ---
	# TODO: Use a working conservative definition of (log_2(FPKM + 0.001) ~ 0) 

	## --- Scatter plot ---
	a_plt <- 
		to_plot |>
		ggplot(aes(x = tpm_median, y = tpm_mad, size = n, color = cohort, fill = cohort, label = lbl)) +
		geom_vline(xintercept = 0, lty = 2) +
		geom_point(alpha = 0.75, shape = 21) +
		ggrepel::geom_text_repel(size = 3.5, max.overlaps = Inf) +
		paletteer::scale_color_paletteer_d("colorBlindness::SteppedSequential5Steps") +
		paletteer::scale_fill_paletteer_d("colorBlindness::SteppedSequential5Steps") +
		labs(
			x = "Median expression (log2(fpkm+0.001)",
			y = "Expression variance (Median Absolulte Deviation)",
			size = "Sample size",
			color = "Cohort",
			fill = "Cohort",
			caption = "Data: MET500, Nature. 2017 Aug 17;548(7667):297-303. doi: 10.1038/nature23306."
		)
	a_plt <- 
		stylehaven::add_logo(
			a_plt, 
			image = "./IMG/2k B.png", 
			position = "right", 
			height = 0.025
		)
	
	## --- Barplot of proportion expressed ---
	b_plt <- 
		to_plot |>
		dplyr::filter(SYMBOL == id) |>
		dplyr::mutate(cohort = forcats::fct_reorder(cohort, prop_expressed, .fun = sum)) |>
		ggplot(aes(x = cohort, y = prop_expressed, color = cohort, fill = cohort)) +
	    geom_bar(position = "stack", stat = "identity") +
		#coord_flip() + # Flip the x and y coordinates
		#facet_wrap(vars(passed), nrow = 2, scales = "free_y") +
		geom_hline(yintercept = 1) +
		paletteer::scale_color_paletteer_d("colorBlindness::SteppedSequential5Steps") +
		paletteer::scale_fill_paletteer_d("colorBlindness::SteppedSequential5Steps") +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		labs(
			title = "",
			x = "",
			y = "Proportion of samples target is expressed in (0-1)"
		)	
	b_plt <- 
		stylehaven::add_logo(
			b_plt, 
			image = "./IMG/2k B.png", 
			position = "right", 
			height = 0.025
		)
	
	## --- Beeswarm of panel ---

	all_cols <- 
		paletteer::palettes_d_names |>
		dplyr::filter(length >= length(unique(to_plot$cohort))) |>
		dplyr::filter(type == "qualitative")
			
	library(ggbeeswarm)
	c_plt <-
		to_plot |>
		dplyr::filter(SYMBOL %in% annotation) |>
		dplyr::mutate(lbl = ifelse(prop_expressed < 1, cohort, "")) |>
		ggplot(aes(y = SYMBOL, x = prop_expressed, label = lbl)) +
		ggbeeswarm::geom_quasirandom(aes(fill = cohort, color = cohort), orientation = 'y', alpha = 0.75, shape = 21, size = 2) +
		# geom_point(alpha = 0.75) + 
		stat_summary(
			#orientation = "y",
			fun = "median",      # Specify the function to compute
			geom = "point",        # Use a point for the geometry
			shape = 124,            # Shape 95 is a horizontal bar
			size = 4,             # Make the bar thick and visible
			color = "black"          # Use a distinct color
		) +
		ggrepel::geom_text_repel(size = 3.5, max.overlaps = Inf) +
		paletteer::scale_color_paletteer_d("khroma::discreterainbow") +
		paletteer::scale_fill_paletteer_d("khroma::discreterainbow") +
		labs(
			title = glue::glue("Loss of expression in metastatic samples"),
			caption = "Data: MET500, Nature. 2017 Aug 17;548(7667):297-303. doi: 10.1038/nature23306.",
			x = "Proportion of samples target is expressed in (0-1)",
			y = "",
			color = "Cancer type",
			fill = "Cancer type" 
		) +
		scale_y_discrete(
			labels = generate_bold_labels(id)
		) + 
		theme_minimal() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1),
			legend.position = "top"
		)
	c_plt <- 
		stylehaven::add_logo(
			c_plt, 
			image = "./IMG/2k B.png", 
			position = "left", 
			height = 0.025
		)
	return(list(a_plt = a_plt, b_plt = b_plt, c_plt = c_plt))	

}


get_binary_exp <- function(x, threshold = 2) {

	# Note: For gene expression RNA-seq data represented as RSEM expected counts that have been standardized/normalized using DESeq2 and transformed via log2(x+1) , a reasonable cutoff to classify a gene as "not expressed" generally falls in the range of 0.1 to 1.0, depending on your desired stringency.

	## --- Scaffold ---
	# x = exp_gene_pheno_tb2
	
	y <- log2(threshold + 1)
	out <-
		x |>
		tidyr::pivot_longer(contains("ENSG"), names_to = "ENSGID", values_to = "tpm") |>
		dplyr::group_by(study, primary_disease_or_tissue, ENSGID) |>
		dplyr::summarise(
			n = dplyr::n(),
			tpm_median = median(tpm),
			tpm_mad = mad(tpm),
			prop_expressed = sum(tpm > y) / length(tpm),
			.groups = "drop_last"
		) |>
		#dplyr::mutate(scr = tpm_median * (1 - tpm_mad))|>
		dplyr::ungroup() |>
		dplyr::arrange(desc(tpm_median), tpm_mad) 

	return(out)
	
}



write_binary_exp <- function(x, y, filters, file_out) {
 
	# tar_load(exp_binary); x = exp_binary; tar_load(ensembl_ids); y = ensembl_ids; filters = consort_plt$data; file_out = "../../output/tumor_mrna_hpa.xlsx"
 
	library(dplyr)
	
	x |>
	dplyr::filter(study != "GTEX") |>
	dplyr::left_join(y[, c("ENSGID", "SYMBOL")], by = "ENSGID") |>
	dplyr::filter(SYMBOL %in% filters$SYMBOL) |>
	dplyr::group_by(ENSGID) |>
	dplyr::arrange(desc(tpm_median)) |>	
	dplyr::mutate(indication_rank = dplyr::row_number(-tpm_median)) |>	
	writexl::write_xlsx(file_out)
 
	return(file_out)

}


plot_binary_exp <- function(x, y, filters) {
 
	# tar_load(exp_binary); x = exp_binary; tar_load(ensembl_ids); y = ensembl_ids; filters = consort_plt$data; file_out = "../../output/tumor_mrna_hpa.svg"
 
	library(ggplot2)
	library(dplyr)
	
	to_plot <- 
		x |>
		dplyr::filter(study != "TARGET") |>
		dplyr::left_join(y[, c("ENSGID", "SYMBOL")], by = "ENSGID") |>
		dplyr::filter(SYMBOL %in% filters$SYMBOL) |>
		dplyr::group_by(ENSGID) |>
		dplyr::arrange(desc(tpm_median)) 
		#dplyr::mutate(SYMBOL = forcats::fct_reorder(SYMBOL, tpm_median))

	## --- Create the bubble plot ---
	a_plt <-
		to_plot |>
		ggplot(aes(y = primary_disease_or_tissue, x = SYMBOL, size = prop_expressed, fill = tpm_median, color = tpm_median)) +
		geom_point(alpha = 0.7) + 
		facet_wrap(vars(study), nrow = 2, scales = "free_y") +
		# scale_size(
			# name = "Expression (mRNA)", 
             # range = c(2, 12)
		# ) +      # Control the min and max bubble size
		paletteer::scale_fill_paletteer_c("grDevices::Mako", direction = -1) + 
		paletteer::scale_color_paletteer_c("grDevices::Mako", direction = -1) +
		labs(
			title = "",
			caption = "Data: TCGA/GTEx in UCSC Xena",
			x = "",
			y = "",
			size = "Proportion of samples \nwith non-zero expression", 
			color = "mRNA Expression (Median expected count, log2(x+1))",
			fill = "mRNA Expression (Median expected count, log2(x+1))" 
		) +
		theme_minimal() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1)
		)
	a_plt <- 
		stylehaven::add_logo(
			a_plt, 
			image = "./IMG/2k B.png", 
			position = "left", 
			height = 0.025
		)
	return(a_plt)

}


get_parafac <- function(x, verbose = FALSE, k = 6) {

	## Scaffolding
	# tar_load(gene_emd_max); x = gene_emd_max

	## --- Load necessary packages ---
	library(tibble)
	library(dplyr)
	library(tidyr)
	library(multiway)

	## --- Convert the tibble into a 3D array (tensor) ---
	dim_names <- 
		list(
			tissue_t = unique(x$primary_disease_or_tissue),
			tissue_n = unique(x$ref_tissue),
			GENEID = unique(x$GENEID)
		)
	tensor_arr <- 
		x |>
		dplyr::arrange(primary_disease_or_tissue, ref_tissue, GENEID) |>
		dplyr::pull(emd) |>
		array(
			dim = c(
				length(dim_names$tissue_t), 
				length(dim_names$tissue_n), 
				length(dim_names$GENEID)
			),
			dimnames = dim_names
		)

	## --- Run the PARAFAC analysis ---
	# 'nfac' is the number of factors to extract.
	# 'const' applies non-negativity constraints to the modes.
	#  aim for the highest number of components that retains a high core consistency.
	
	parafac_ls <- vector("list", k)
	
	for(i in 2:k) {
		parafac_mdl <- 
			multiway::parafac(
				tensor_arr, 
				nfac = i, 
				# parallel = TRUE,
				const = c("nonneg", "nonneg", "nonneg")
			)	
		cc <- 
			multiway::corcondia(
				tensor_arr, 
				parafac_mdl
			)
		parafac_ls[[i]] <- 
			list(
				mdl = parafac_mdl,
				cc = cc
			)
	}
	


	## --- Inspect the results ---
	if(verbose) {
		summary(parafac_ls[[i]]$mdl)

		cat("\n--- Loadings for Mode A (Indications) ---\n")
		print(parafac_ls[[i]]$mdl$A)

		cat("\n--- Loadings for Mode B (Healthy tissues) ---\n")
		print(parafac_ls[[i]]$mdl$B)

		cat("\n--- Loadings for Mode C (Gene) ---\n")
		print(parafac_ls[[i]]$mdl$C)	
	}

	return(parafac_ls)

}


#' Convert an renv.lock file into a Dockerfile
#'
#' @param renv_file Path to the renv.json or renv.lock file.
#' @param output_dockerfile Path where the Dockerfile will be saved. Default: "Dockerfile".
#' @param method Either "renv" (copies lockfile and runs renv::restore; recommended for reproducibility)
#'               or "explicit" (parses packages into explicit R install statements using pak/remotes).
#' @param base_image Custom base Docker image. If NULL, auto-selects 'rocker/r-ver:<R-version>'.
#' @param use_ppm_binaries Logical; if TRUE, configures Posit Package Manager (PPM) binary repository
#'                         for Ubuntu to enable fast, pre-compiled package installations instead of slow source builds.
#' @param extra_sysreqs Character vector of additional Ubuntu apt packages to install.
#' @param workdir Working directory inside the Docker container. Default: "/project".
#' @param entry_cmd Command executed when container runs, e.g. 'CMD ["R"]' or 'CMD ["Rscript", "00_run_targets.R"]'.
#' @return Invisibly returns the Dockerfile content as a character vector.
renv_to_dockerfile <- function(
  renv_file = "renv.json",
  output_dockerfile = "Dockerfile",
  use_ppm_binaries = TRUE,
  entry_cmd = 'CMD ["R"]'
) {

	## Scaffolding
	# renv_file = "renv.json"; use_ppm_binaries = TRUE; entry_cmd = 'CMD ["R"]'
	
	## Dependencies
	library(jsonlite)
	
	## Input data
	lock <- jsonlite::fromJSON(renv_file, simplifyVector = FALSE)

	## --- Detect R version ---
	r_version <- lock$R$Version
	if (is.null(r_version) || !nzchar(r_version)) {
		r_version <- "4.3.2" # Fallback if unspecified
		logger::log_warn("R version not found in lockfile. Defaulting to 4.3.2.")
	}
    base_image <- glue::glue("rocker/r-ver:{r_version}")

	## --- Standard system libraries --- 
	default_sysreqs <- 
		c(
			"build-essential",
			"libcurl4-openssl-dev",
			"libssl-dev",
			"libxml2-dev",
			"libgit2-dev",
			"zlib1g-dev",
			"libfontconfig1-dev",
			"libharfbuzz-dev",
			"libfribidi-dev",
			"libfreetype6-dev",
			"libpng-dev",
			"libtiff5-dev",
			"libjpeg-dev",
			"git",
			"wget"
		)
	sysreqs_str <- paste(default_sysreqs, collapse = " \\\n    ")

	## --- Build Dockerfile lines ---
	dfile <- 
		c(
	sprintf("# Generated automatically from %s", basename(renv_file)),
	sprintf("FROM %s", base_image),
	"",
	"# Set non-interactive environment for apt & renv",
	"#ENV DEBIAN_FRONTEND=noninteractive",
	"#ENV RENV_PATHS_CACHE=/root/.cache/R/renv"
	)

  # Configure binary package repository (RSPM/PPM) if enabled for fast installs
  if (use_ppm_binaries) {
    dfile <- c(
      dfile,
      "",
      "# Use Posit Package Manager for pre-compiled Linux binaries (speeds up builds dramatically)",
      "RUN echo 'options(repos = c(CRAN = \"https://packagemanager.posit.co/cran/__linux__/jammy/latest\"))' >> /usr/local/lib/R/etc/Rprofile.site"
    )
  }

  # Install Linux system dependencies
  dfile <- c(
    dfile,
    "",
    "# Install core system dependencies",
    "RUN apt-get update && apt-get install -y --no-install-recommends \\",
    sprintf("    %s \\", sysreqs_str),
    "    && rm -rf /var/lib/apt/lists/*",
    ""
  )

    # ----------------------------------------------------
    # Explicit Package Installation
    # ----------------------------------------------------
	packages <- lock$Packages
    cran_pkgs <- character()
    bioc_pkgs <- character()
    github_pkgs <- character()

	for (pkg_name in names(packages)) {
		pkg_meta <- packages[[pkg_name]]
		source_type <- pkg_meta$Source %||% "Repository"

		if (identical(source_type, "Bioconductor")) {
			bioc_pkgs <- c(bioc_pkgs, pkg_name)
		}else if (identical(source_type, "GitHub")) {
			user <- pkg_meta$RemoteUsername
			repo <- pkg_meta$RemoteRepo
			ref <- pkg_meta$RemoteRef %||% pkg_meta$RemoteSha %||% "HEAD"
			github_pkgs <- c(github_pkgs, sprintf("%s/%s@%s", user, repo, ref))
		}else{
		cran_pkgs <- c(cran_pkgs, pkg_name)
		}
	}

    dfile <- c(
      dfile,
      "",
      "# Install pak for robust, multi-threaded package installations with sysreq resolution",
      'RUN R -e "install.packages(\'pak\', repos = \'https://r-lib.github.io/p-pkg\')"'
    )

    if (length(cran_pkgs) > 0) {
      chunk_size <- 40
      chunks <- split(cran_pkgs, ceiling(seq_along(cran_pkgs) / chunk_size))
      for (chunk in chunks) {
        pkg_list_str <- paste(sprintf('\"%s\"', chunk), collapse = ", ")
        dfile <- c(
          dfile,
          sprintf('RUN R -e \'pak::pkg_install(c(%s))\'', pkg_list_str)
        )
      }
    }

    if (length(bioc_pkgs) > 0) {
      bioc_list_str <- paste(sprintf('\"%s\"', bioc_pkgs), collapse = ", ")
      dfile <- c(
        dfile,
        sprintf('RUN R -e \'pak::pkg_install(c(%s))\'', bioc_list_str)
      )
    }

    if (length(github_pkgs) > 0) {
      for (gh_repo in github_pkgs) {
        dfile <- c(
          dfile,
          sprintf('RUN R -e \'pak::pkg_install(\"%s\")\'', gh_repo)
        )
      }
    }
  
	writeLines(dfile, con = output_dockerfile)
	logger::log_info("Generated Dockerfile at: {output_dockerfile}")
	invisible(dfile)
}

