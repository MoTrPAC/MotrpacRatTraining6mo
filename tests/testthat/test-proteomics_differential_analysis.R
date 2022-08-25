# # test for limma_res_extract_se()
# test_that("limma_res_extract_se() returns same output", {
#   expect_snapshot(cat(limma_res_extract_se()))
# })

# test for proteomics_timewise_da()
test_that("PROT timewise DA returns same p-values", {
  expect_snapshot(cat(proteomics_timewise_da("PROT", "HEART")$p_value))
})

# test for proteomics_training_da()
test_that("PROT training DA returns same p-values", {
  expect_snapshot(cat(proteomics_training_da("PROT", "HEART")$p_value))
})
