test_that("plotting works for sce", {
    expect_error(
        plot_marker_features(small_example_dataset, group_by = "gene_snn_res.1"),
        NA
    )
})
