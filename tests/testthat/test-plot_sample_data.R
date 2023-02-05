test_that("warn if feature does not exist", {
  expect_warning(plot_feature_normalized_data(assay = "METAB",
                                              tissue = "ADRNL",
                                              feature_ID = "glucose"))
  expect_warning(plot_feature_logfc(assay = "METAB",
                                    tissue = "ADRNL",
                                    feature_ID = "glucose"))
})

test_that("feature is correctly passed for logfc", {
  expect_message(plot_feature_normalized_data(feature = "METAB;SKM-GN;meta-reg:GMP",
                                              add_adj_p = TRUE))
})

test_that("RefMet IDs are handled", {
  expect_message(plot_feature_normalized_data(assay = "METAB",
                                              tissue = "LIVER",
                                              feature_ID = "Lactic acid",
                                              add_adj_p = TRUE))
  expect_message(plot_feature_logfc(assay = "METAB",
                                    tissue = "LIVER",
                                    feature_ID = "Lactic acid",
                                    add_adj_p = TRUE))
})
