data_loader <- 
	list(
		# Input gene list
		tar_target(params_precedented_f, "target_list.json", format = "file"),
		tar_target(params, jsonlite::read_json(params_precedented_f, simplifyVector = TRUE)),
		# Aux. data
		tar_target(ebrt_f, file.path(params$data_commons, "tissue_specificity/wahl_et_al.csv"), format = "file"),
		tar_target(ebrt_tb, readr::read_csv(ebrt_f)), 		
		tar_target(brca_file, file.path(params$data_commons, "brca/1-s2.0-S2666979X21000835-mmc2.xlsx"), format = "file"),
		tar_target(brca_hist, readxl::read_excel(brca_file))
	)

data_collation <- 
	list(
		# Get biomart annotations (use get_biomart when backup)
		tar_target(ensembl_ids, get_ensembl_ids(params$targets))
	)

###############################################################################
# Efficacy: Target enrichment: mRNA (external)
tumour_mrna_deg_module <- 
	list(
		## --- Load data ---
		tar_target(exp_gene_tb, ucsc_toil_gene(ensembl_id = ensembl_ids)),
		tar_target(exp_gene_pheno_tb, ucsc_toil_pheno(x = exp_gene_tb, y = brca_hist)),
		## --- Calculate Earth Mover Distance ---
		tar_target(comparison_plan, get_comparison_plan(x = exp_gene_pheno_tb)),
		tar_target(
			name = gene_emd_max_i,
			command = 
				get_emd(
					pair_indices = c(comparison_plan$indications, comparison_plan$targets, comparison_plan$test_tissue),
					x = exp_gene_pheno_tb,
					perms = 1
				),
			pattern = map(comparison_plan),
			#resources = tar_resources(crew = tar_resources_crew(controller = "my_local_controller")),
			error = "null", # continue on error
			iteration = "list" # ensures the results are structured cleanly
		),
		tar_target(
			name = gene_emd_max,
			command = dplyr::bind_rows(gene_emd_max_i)
		),
		tar_target(
			name = parafac_mdl,
			command = get_parafac(x = gene_emd_max)),
		tar_target(gene_emd_max_heat_static_plt, plot_static_heatmap(x = gene_emd_max, y = ensembl_ids)),
		tar_target(name = tumour_mrna_emd, command = write_heatmap(x = gene_emd_max_heat_static_plt, file_out = "../../output/tumour_mrna_emd.svg"), format = "file"),
		tar_target(
			name = gene_emd_max_rds,
			command = export_tsv(x = gene_emd_max, file_out = "../../output/gene_emd_max.rds"),
			format = "file"
		)
	)

