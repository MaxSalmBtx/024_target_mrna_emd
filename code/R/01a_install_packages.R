# author      = "Max Salm"
# copyright   = "MIT"
# license     = "NA"
# version     = "1.0.0"
# status      = "Development"
# Aim         = Simple script to install all necessary R packages

### Check web access first
if ( Sys.info()["sysname"] == "Linux" ) {
	cmd <- "wget -p https://www.google.com/"
	test_connection <- 
		system(cmd, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
	if (length(test_connection) == 0) {
		cat("Connection to internet verified.\n")
		unlink("www.google.com", recursive = TRUE)
	}else{
		stop("No connection to internet found.\n")
	}
}

#########################
### Base dependencies ###
#########################
if (Sys.getenv("BIOENV_IMAGE") == "") {
    .libPaths(c("/home/m.salm/R/personal_lib", .libPaths()))
}
if (Sys.getenv("USER") == "m.salm") {
  personal_repo <- "../../code/R/r_lib"
} else {
  personal_repo <- .libPaths()[1]
}

### Helper functions
check_package <- function(
  package_names, 
  dest_dir, 
  fresh_install = FALSE, 
  repository = "CRAN"
){

	#' Check and install R packages from various repositories.
	#'
	#' This function checks if a set of R packages are installed in a specified
	#' local repository. If not, it installs them from CRAN, or BIOC.
	#' It also allows for a fresh install, removing previous versions if needed.
	#'
	#' @param package_names A character vector of package names to check and install.
	#' @param dest_dir The destination directory for the local repository.
	#' @param fresh_install A logical value indicating whether to remove previous
	#'   installations of the packages. Defaults to FALSE.
	#' @param repository A character string specifying the repository to use.
	#'   Must be one of "CRAN" or "BIOC". Defaults to "CRAN".
	#'
	#' @return None (invisibly). The function installs packages and loads them.
	#'
	#' @examples
	#' \dontrun{
	#' # Check and install packages from CRAN
	#' checkPackage(
	#'   package_names = c("ggplot2", "dplyr"),
	#'   dest_dir = "~/R/my_packages"
	#' )
	#'
	#' # Perform a fresh install from CRAN
	#' checkPackage(
	#'   package_names = c("MASS", "nlme"),
	#'   dest_dir = "~/R/my_packages",
	#'   fresh_install = TRUE,
	#'   repository = "CRAN"
	#' )
	#' }
	#'
	#' @import checkmate
	#' @export

	# Argument checks using checkmate
	checkmate::assertCharacter(package_names, min.len = 1, any.missing = FALSE)
	checkmate::assertDirectory(dest_dir, access = "w")
	checkmate::assertLogical(fresh_install, len = 1)
	checkmate::assertChoice(repository, choices = c("CRAN", "BIOC"))

	# Check if the local repository directory exists.
	# If it doesn't, create it and update the library paths.
    if (!dir.exists(dest_dir)){
        dir.create(dest_dir, recursive = TRUE, mode = "0777")
        .libPaths(dest_dir)  # Update library paths to include the new directory
    } 
        
	# Get a vector of installed packages in the local repository.
    installed_packages  <- installed.packages(lib.loc = dest_dir)[, "Package"]

	# If fresh_install is TRUE, remove previously installed packages.
    if (fresh_install) {
        cat("Removing previous installations.\n")
        packages_to_remove  <- 
			package_names[package_names %in% installed_packages]
        remove.packages(packages_to_remove, lib = dest_dir)
    }

	# Determine which packages need to be newly installed.
    new_packages <- package_names[!(package_names %in% installed_packages)]
    if (length(new_packages) > 0) {
        message("Installing:")
        message(paste0(new_packages, "\n"))
		# Install packages based on the specified repository.		
        if(repository == "CRAN"){
            install.packages(new_packages, 
                            dependencies = TRUE, 
                            repos = "https://cran.rstudio.com/", 
                            lib = dest_dir)      
        }else if(repository == "BIOC"){
		if (!require("BiocManager", quietly = TRUE))
			install.packages("BiocManager")
		BiocManager::install(version = "3.23") # compatible with R 4.6.0				
            BiocManager::install(
              pkgs = new_packages, 
              lib = dest_dir, 
              suppressUpdates=TRUE, 
              ask=FALSE
            )     
        }else{
            stop("checkPackage(): Select from CRAN/BIOC")
        }                
    }else{
        logger::log_info("All packages found.\n")
    }

	# Load all requested packages.   
	package_load_status <- 
		lapply(
			package_names, 
			function(pkg) {
				tryCatch(
					{
						library(pkg, character.only = TRUE, lib.loc = dest_dir)
						TRUE # Return TRUE if loading succeeds
					}, 
					error = function(e) {
						cat(paste("Failed to load ", pkg, ":", e$message, "\n"))
						FALSE # Return FALSE if loading fails
					}
				)
			}
		)

	if (any(!unlist(package_load_status))) {
		logger::log_warn("Some packages failed to load. Check output for details.")		
	}
}


### --- List all libraries in source code ---

#' Find all lines containing a specific text pattern in files within a folder.
#'
#' @param dir_path The path to the directory to search. Defaults to the current directory (".").
#' @param pattern The string pattern to search for (e.g., "library").
#' @param is_case_sensitive A boolean to control whether the search is case-sensitive. Defaults to FALSE.
#' @param search_recursively A boolean to control whether to search in subdirectories. Defaults to TRUE.
#'
#' @return A named list where each name is a file path and the value is a character
#'         vector of the lines in that file that matched the pattern.
#'
find_library_files <- function(dir_path = "./function_lib", pattern, is_case_sensitive = FALSE, search_recursively = TRUE) {

	# dir_path = "./function_lib"; pattern = "library"; is_case_sensitive = FALSE; search_recursively = TRUE

	# 1. Get a list of all files in the directory.
	#    full.names = TRUE gives the complete file path.
	all_files <- 
		list.files(
			path = dir_path,
			recursive = search_recursively,
			full.names = TRUE,
			all.files = FALSE, # Include hidden files
			no.. = TRUE      # Exclude '.' and '..' pseudo-directories
		)

	# Filter out directories, keeping only actual files
	file_paths <- all_files[!sapply(all_files, function(x) file.info(x)$isdir)]

	# 2. Use lapply to iterate over each file path.
	#    This will return a list of results.
	results_list <- 
		lapply(file_paths, function(file_path) {
			# Use a try-catch block to gracefully handle errors, such as trying
			# to read a binary file that can't be parsed as text.
			inputs <- tryCatch({
				# Read all lines from the current file
				readLines(file_path, warn = FALSE)
			}, error = function(e) {
				# If an error occurs, print a message and return NULL for this file
				message(sprintf("Could not read file: %s", file_path))
				return(NULL)
			})

			# If the file was unreadable, lines will be NULL, so we skip it
			if (is.null(inputs)) {
			  return(NULL)
			}

			## --- Use grep to find the lines that contain the pattern.
			matching_lines <- 
				grep(pattern, inputs, value = TRUE, ignore.case = !is_case_sensitive)
			
			## --- extract library names ---
			matches <- unique(stringr::str_match(matching_lines, "\\((.*?)\\)")[, 2])
			return(matches)
	})

	out <- unique(unlist(results_list[sapply(results_list, function(x) length(x) > 0)]))

	return(out)
}


# Search the current folder and all subfolders for the string "library" (case-insensitive)
found_libs <- find_library_files(pattern = "library")

check_package(
	package_names = found_libs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)




base_pkgs <- c("quarto", "RColorBrewer", "pander", "checkmate", "docstring", "logger", "testthat", "readxl", "janitor", "stringdist", "targets", "tarchetypes", "assertions", "overlapping", "tidylog")
check_package(
	package_names = base_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)






####################
### Bioconductor ###
####################
bioc_pkgs <- c("EnsDb.Hsapiens.v86", "ComplexHeatmap", "EMDomics", "org.Hs.eg.db", "AnnotationDbi")
check_package(
	package_names = bioc_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "BIOC"
)


############
### CRAN ###
############
mran_pkgs <- c("rmarkdown", "readr", "Hmisc", "pander", "R.utils", "ggrepel", "UCSCXenaShiny", "sparklyr", "sparklyr.nested", "httr", "httr2", "jsonlite", "DescTools", "corrr")
check_package(
	package_names = mran_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)

table_pkgs <- c("kableExtra", "pander", "gt", "gtExtras", "gtsummary")
check_package(
	package_names = table_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)
    
plot_pkgs <- 
	c("ggplot2", "viridis", "ggpubr",  "magick", "showtext", "gplots", "ggExtra", "ggdist", "MetBrewer", "EnvStats", "ggtext", "colorRamp2", "iheatmapr", "legendry", "tidyclust", "GGally", "gghighlight")
check_package(
	package_names = plot_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)


# remotes::install_github("ropensci/UCSCXenaTools") # Needed to get the most up-to-date version, bug fix
# remotes::install_github("jespermaag/gganatogram") # organ maps
# remotes::install_github("CT-Data-Haven/stylehaven") # Logo annotations
# remotes::install_github("alex-bio/UniProtExtractR") # Tidy access to UniProtKB
# remotes::install_github("hrbrmstr/ggchicklet") # stacked barplots
# remotes::install_github("steveneschrich/surfaceome") # surfaceome
# remotes::install_github("hughjonesd/ggmagnify")
# remotes::install_github("anhtr/HPAanalyze")
