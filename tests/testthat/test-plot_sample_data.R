test_that("warn if feature does not exist", {
  expect_warning(plot_feature_normalized_data(assay = "METAB",
                                              tissue = "ADRNL",
                                              feature_ID = "glucose"))
  expect_warning(plot_feature_logfc(assay = "METAB",
                                    tissue = "ADRNL",
                                    feature_ID = "glucose"))
})
