# tar_make(script = "_targets.R")

## To execute this in an R session:
# tar_make(callr_function = NULL)



# Core libraries
library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)

# Source functions from the 'function_lib' directory
tar_source(
	files = "function_lib",
	envir = targets::tar_option_get("envir"),
	change_directory = FALSE
)

## --- Parelellisation ---
# https://books.ropensci.org/targets/crew.html
# 13.5 Heterogeneous workers: these will be deployed at the target level for efficiency on local compute
controller_local <- 
	crew::crew_controller_local(
		name = "my_local_controller",
		workers = 8, 
		host = "127.0.0.1",
		garbage_collection = TRUE,
		options_local = crew_options_local(log_directory = "./crew_logs", log_join = FALSE)
	)


## --- Get global package dependencies for targets ---

tar_option_set(
	packages = c("quarto", "logger", "readr", "readxl", "janitor", "dplyr", "tidyr", "checkmate", "ensembldb", "EnsDb.Hsapiens.v86", "org.Hs.eg.db", "AnnotationDbi", "ggplot2", "ggrepel", "paletteer", "plotly", "pander", "gtExtras", "ComplexHeatmap"),
	memory = "transient", # extra options for downstream workers
	storage = "worker", 
	retrieval = "worker",	
	garbage_collection = TRUE,
	workspace_on_error = TRUE
)

## --- Graphical parameters ---

library(sysfonts)
sysfonts::font_add_google("Lato", "lato")
sysfonts::font_add_google("Open Sans", "opensans")
showtext::showtext_auto()
library(ggplot2)
old_theme <- 
  ggplot2::theme_set(ggplot2::theme_minimal()) +
  ggplot2::theme_update(
    legend.title = element_blank(),
    # legend.justification = c(0, 1), 
    # legend.position = c(.1, 1.075),
    legend.background = element_blank(),
    axis.title= element_text(family = "opensans", size = 10),
    plot.title = element_text(family = "lato", size = 20, margin = margin(b = 10)),
    plot.subtitle = element_text(family = "opensans", size = 10, color = "darkslategrey", margin = margin(b = 25)),
    plot.caption = element_text(family = "opensans", size = 8, margin = margin(t = 10), color = "grey70", hjust = 0)
  )

## --- Workflow modules ---
list(
	data_loader,
	data_collation,
	tumour_mrna_deg_module
)

