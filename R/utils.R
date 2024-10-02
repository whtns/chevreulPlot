

#' Record Experiment Metadata
#'
#' Records miscellaneous data
#' @param object A object
#' @param experiment_name name of the experiment
#' @param organism human or mouse
#'
#' @return a SingleCellExperiment object
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' record_experiment_data(small_example_dataset)
#'
record_experiment_data <- function(object, 
                                   experiment_name = "default_experiment", 
                                   organism = "human") {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' needed for this function to work. 
             Please install it.",
            call. = FALSE
        )
    }

    organism <- metadata(object)[["experiment"]][["organism"]] %||% organism

    experiment_name <- 
        metadata(object)[["experiment"]][["experiment_name"]] %||% experiment_name

    message(glue("[{format(Sys.time(), '%H:%M:%S')} 
                 Logging Technical Details..."))
    experiment <- list(
        experiment_name = experiment_name,
        organism = organism
    )
    experiment$date_of_export <- Sys.Date()
    experiment$date_of_analysis <- Sys.Date()

    experiment$parameters <- list(
        gene_nomenclature = "gene_symbol",
        discard_genes_expressed_in_fewer_cells_than = 10,
        keep_mitochondrial_genes = TRUE,
        variables_to_regress_out = "nCount_RNA",
        number_PCs = 30,
        tSNE_perplexity = 30,
        cluster_resolution = seq(0.2, 2.0, by = 0.2)
    )
    experiment$filtering <- list(
        UMI_min = 50,
        genes_min = 10
    )
    experiment$sessionInfo <- list(
        capture.output(sessionInfo())
    )

    if (!is.null(objectVersion(object))) {
        experiment$SingleCellExperiment_version <- objectVersion(object)
    }

    experiment$chevreul_version <- packageVersion("chevreul")

    metadata(object)[["experiment"]] <- experiment

    return(object)
}


#' Update a database of chevreul projects
#'
#' Add new/update existing projects to the database by recursing fully
#'
#' @param projects_dir The project directory to be updated
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#' @return a sqlite database with SingleCellExperiment objects
update_project_db <- function(
        projects_dir = NULL,
        cache_location = "~/.cache/chevreul",
        sqlite_db = "single-cell-projects.db",
        verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    projects_tbl <-
        dir_ls(projects_dir, glob = "*.here", recurse = TRUE, 
               fail = FALSE, all = TRUE) %>%
        path_dir(.) %>%
        set_names(path_file(.)) %>%
        enframe("project_name", "project_path") %>%
        mutate(project_slug = str_remove(project_name, "_proj$")) %>%
        mutate(project_type = path_file(path_dir(project_path))) %>%
        identity()

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl") %>%
        filter(file.exists(project_path)) %>%
        filter(!project_path %in% projects_tbl$project_path) %>%
        bind_rows(projects_tbl) %>%
        distinct(project_path, .keep_all = TRUE)

    dbWriteTable(con, "projects_tbl", projects_tbl, overwrite = TRUE)

    dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Append projects to database
#'
#' @param new_project_path new project path
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#' @return a sqlite database with SingleCellExperiment objects

append_to_project_db <- function(
        new_project_path,
        cache_location = "~/.cache/chevreul",
        sqlite_db = "single-cell-projects.db",
        verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    projects_tbl <-
        new_project_path %>%
        set_names(path_file(.)) %>%
        enframe("project_name", "project_path") %>%
        mutate(project_slug = str_remove(project_name, "_proj$")) %>%
        mutate(project_type = path_file(path_dir(project_path))) %>%
        identity()

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl") %>%
        filter(file.exists(project_path)) %>%
        filter(!project_path %in% projects_tbl$project_path) %>%
        bind_rows(projects_tbl) %>%
        distinct(project_path, .keep_all = TRUE)

    dbWriteTable(con, "projects_tbl", current_projects_tbl, overwrite = TRUE)

    dbDisconnect(con)
}

#' Read a database of chevreul projects
#'
#' Reads database of chevreul projects to a data frame
#'
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#' @return a tibble with SingleCellExperiment objects
#

read_project_db <- function(
        cache_location = "~/.cache/chevreul",
        sqlite_db = "single-cell-projects.db",
        verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl")

    dbDisconnect(con)

    return(current_projects_tbl)
}

#' Make Bigwig Database
#'
#'
#' @param new_project Project directory
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db containing bw files
#'
#' @return a sqlite database of bigwig files for cells 
#' in a SingleCellExperiment object
make_bigwig_db <- function(new_project = NULL, 
                           cache_location = "~/.cache/chevreul/", 
                           sqlite_db = "bw-files.db") {
    new_bigwigfiles <- dir_ls(new_project, glob = "*.bw", recurse = TRUE) %>%
        set_names(path_file(.)) %>%
        enframe("name", "bigWig") %>%
        mutate(sample_id = 
                   str_remove(name, "_Aligned.sortedByCoord.out.*bw$")) %>%
        filter(!str_detect(name, "integrated")) %>%
        distinct(sample_id, .keep_all = TRUE) %>%
        identity()

    con <- dbConnect(SQLite(), dbname = path(cache_location, sqlite_db))

    all_bigwigfiles <-
        dbReadTable(con, "bigwigfiles") %>%
        bind_rows(new_bigwigfiles)

    dbWriteTable(con, "bigwigfiles", all_bigwigfiles, overwrite = TRUE)

    return(all_bigwigfiles)
}

#' Retrieve Metadata from Batch
#'
#' @param batch batch
#' @param projects_dir path to project dir
#' @param db_path path to .db file
#'
#' @return a tibble with cell level metadata from a SingleCellExperiment object
metadata_from_batch <- function(
        batch, projects_dir = "/dataVolume/storage/single_cell_projects",
        db_path = "single-cell-projects.db") {
    mydb <- dbConnect(SQLite(), path(projects_dir, db_path))

    projects_tbl <- dbReadTable(mydb, "projects_tbl") %>%
        filter(!project_type %in% c("integrated_projects", "resources"))

    dbDisconnect(mydb)

    metadata <-
        projects_tbl %>%
        filter(project_slug == batch) %>%
        pull(project_path) %>%
        path("data") %>%
        dir_ls(glob = "*.csv") %>%
        identity()
}

#' Clean Vector of Chevreul Names
#'
#' Cleans names of objects provided in a vector form
#'
#' @param myvec A vector of object names
#'
#' @return a clean vector of object names
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' make_chevreul_clean_names(colnames(
#' get_cell_metadata(small_example_dataset)))
make_chevreul_clean_names <- function(myvec) {
    myvec %>%
        set_names(str_to_title(str_replace_all(., 
                                               "[^[:alnum:][:space:]\\.]",
                                               " ")))
}

#' Get metadata from object
#'
#' Get metadata from the given object
#'
#' @param object a SingleCellExperiment object
#'
#' @return a tibble with metadata from a SingleCellExperiment object
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' metadata_from_object(small_example_dataset)
#'
metadata_from_object <- function(object) {
    colnames(colData(object))
}

