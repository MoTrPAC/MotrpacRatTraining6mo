input = data.frame(data.table::data.table(MotrpacRatTraining6moData::GRAPH_PW_ENRICH)[
  tissue %in% c("SKM-GN","HEART") &
  ome %in% c("PROT","TRNSCRPT") & 
  grepl(":8w_F1_M1", cluster)])

test_that("pw_enrich_res is missing columns", {
  expect_error(enrichment_network_vis(pw_enrich_res = input[,!colnames(input) == "adj_p_value"]))
  expect_error(enrichment_network_vis(pw_enrich_res = input[,!colnames(input) == "ome"]))
  expect_error(enrichment_network_vis(pw_enrich_res = input[,!colnames(input) == "tissue"]))  
  expect_error(enrichment_network_vis(pw_enrich_res = input[,!colnames(input) == "intersection"]))
  expect_error(enrichment_network_vis(pw_enrich_res = input[,!colnames(input) == "computed_p_value"]))
  expect_error(enrichment_network_vis(pw_enrich_res = input[,!colnames(input) == "term_name"]))
  expect_error(enrichment_network_vis(pw_enrich_res = input[,!colnames(input) == "term_id"]))
})

test_that("pw_enrich_res is provided with other args", {
  expect_error(enrichment_network_vis(input, tissues="SKM-GN"))
  expect_error(enrichment_network_vis(input, cluster="8w_F1_M1"))
})

test_that("pw_enrich_res has more than one cluster", {
  input = data.frame(data.table::data.table(MotrpacRatTraining6moData::GRAPH_PW_ENRICH)[
    tissue %in% c("SKM-GN","HEART") &
      ome %in% c("PROT","TRNSCRPT")])
  expect_error(enrichment_network_vis(input))
})

test_that("more than one cluster is provided", {
  expect_error(enrichment_network_vis(cluster=c("8w_F1_M1", "8w_F0_M1")))
})

test_that("not all args are specified", {
  expect_error(enrichment_network_vis(cluster="8w_F1_M1"))
  expect_error(enrichment_network_vis(tissue="SKM-GN"))
})

test_that("invalid tissue", {
  expect_error(enrichment_network_vis(tissue="SKMGN", cluster="8w_F1_M1"))
})

test_that("invalid assay", {
  expect_error(enrichment_network_vis(tissue="SKM-GN", cluster="8w_F1_M1", assays="silly"))
})

test_that("invalid cluster", {
  expect_error(enrichment_network_vis(tissue="SKM-GN", cluster="silly"))
})

test_that("no results", {
  df = data.frame(matrix(ncol=ncol(input), nrow=0))
  colnames(df) = colnames(input)
  expect_message(enrichment_network_vis(df))
})

test_that("duplicate pathways", {
  input = data.frame(data.table::data.table(MotrpacRatTraining6moData::GRAPH_PW_ENRICH)[
    tissue %in% c("SKM-GN","HEART") &
      ome %in% c("PROT","TRNSCRPT") & 
      grepl(":8w", cluster)])
  input$cluster = NULL
  expect_error(enrichment_network_vis(input))
})

test_that("returns HTML", {
  out = enrichment_network_vis(tissue="SKM-GN", 
                               cluster="8w_F1_M1",
                               return_html = TRUE, 
                               overwrite_html = TRUE, 
                               out_html = "/tmp/test.HTML")
  expect_equal(out, "/tmp/test.HTML")
})
