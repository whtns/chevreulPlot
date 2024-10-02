#' Run Differential Expression
#'
#'
#'
#' @param object SingleCellExperiment object
#' @param cluster1 cluster 1
#' @param cluster2 cluster 2
#' @param resolution resolution
#' @param diffex_scheme scheme for differential expression
#' @param featureType gene or transcript
#' @param tests t, wilcox, or bimod
#'
#'
#' @return a dataframe with differential expression information
run_object_de <- function(object, cluster1, cluster2, resolution = 0.2,
                          diffex_scheme = "louvain", featureType = "gene",
                          tests = c("t", "wilcox", "bimod")) {

    data_env <- new.env(parent = emptyenv())
    data("grch38", envir = data_env, package = "chevreul")
    data("grch38_tx2gene", envir = data_env, package = "chevreul")
    grch38 <- data_env[["grch38"]]
    grch38_tx2gene <- data_env[["grch38_tx2gene"]]

    match.arg(tests)

    if (featureType == "transcript") object <- altExp(object, "transcript")

    if (diffex_scheme == "louvain") {
        if (query_experiment(object, "integrated")) {
            active_experiment <- "integrated"
        } else {
            active_experiment <- "gene"
        }
        colLabels(object) <- colData(object)[[paste0(active_experiment,
                                                     "_snn_res.", resolution)]]
        object <- object[, colLabels(object) %in% c(cluster1, cluster2)]
        colLabels(object) <- factor(colLabels(object))
    } else if (diffex_scheme == "custom") {
        object <- object[, c(cluster1, cluster2)]
        keep_cells <- c(cluster1, cluster2)
        new_idents <- c(rep(1, length(cluster1)), rep(2, length(cluster2)))
        names(new_idents) <- keep_cells
        new_idents <- new_idents[colnames(object)]
        colLabels(object) <- new_idents
        cluster1 <- 1
        cluster2 <- 2
    }
    test_list <- vector("list", length(tests))
    for (test in tests) {
        message(test)
        de <- findMarkers(object, test.type = test)
        if (featureType == "transcript") {
            de_cols <- c("enstxp", "ensgene", "symbol", "p_val" = "p.value",
                         "avg_log2FC", "pct.1", "pct.2", "p_val_adj" = "FDR")
            de <- de[[1]] |>
                as.data.frame() |>
                rownames_to_column("enstxp") |>
                left_join(grch38_tx2gene, by = "enstxp") |>
                left_join(grch38, by = "ensgene")
            if ("summary.logFC" %in% colnames(de)) {
                de <- mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
            }
            de <- select(de, any_of(de_cols))
        } else if (featureType == "gene") {
            de_cols <- c("ensgene", "symbol", 
                         "p_val" = "p.value", "avg_log2FC",
                         "pct.1", "pct.2", "p_val_adj" = "FDR")
            de <- de[[1]] |>
                as.data.frame() |>
                rownames_to_column("symbol") |>
                left_join(grch38, by = "symbol")
            if ("summary.logFC" %in% colnames(de)) {
                de <- mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
            }
            de <- select(de, any_of(de_cols))
        }
        test_list[[match(test, tests)]] <- de
    }
    names(test_list) <- tests
    return(test_list)
}




#' Prep Slider Values
#'
#' @param default_val Provide a default value
#'
#' @noRd
prep_slider_values <- function(default_val) {
    min <- round(default_val * 0.25, digits = 1)
    max <- round(default_val * 2.0, digits = 1)
    step <- 10^((ceiling(log10(default_val))) - 1)
    value <- default_val
    return(list(min = min, max = max, value = value, step = step))
}

