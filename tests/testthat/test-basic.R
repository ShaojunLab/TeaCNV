# make test data
data <- data.frame(
  segName = c("chr1_827077_120942799", "chr1_148951389_248907553", 
              "chr2_19449_91766890", "chr2_95076467_108240173"),
  SegMean = c(1.426278, 1.421164, 1.451020, 1.368492),
  ratio_map = c(1.47, 1.47, 1.47, 1.47),
  stringsAsFactors = FALSE
)

data <- matrix(
  c(
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 0,
    2, 3, 4, 5, 6, 7, 8, 9, 0, 1,
    3, 4, 5, 6, 7, 8, 9, 0, 1, 2,
    4, 5, 6, 7, 8, 9, 0, 1, 2, 3,
    5, 6, 7, 8, 9, 0, 1, 2, 3, 4,
    6, 7, 8, 9, 0, 1, 2, 3, 4, 5,
    7, 8, 9, 0, 1, 2, 3, 4, 5, 6,
    8, 9, 0, 1, 2, 3, 4, 5, 6, 7,
    9, 0, 1, 2, 3, 4, 5, 6, 7, 8,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 0,
    2, 3, 4, 5, 6, 7, 8, 9, 0, 1,
    3, 4, 5, 6, 7, 8, 9, 0, 1, 2,
    4, 5, 6, 7, 8, 9, 0, 1, 2, 3,
    5, 6, 7, 8, 9, 0, 1, 2, 3, 4,
    6, 7, 8, 9, 0, 1, 2, 3, 4, 5,
    7, 8, 9, 0, 1, 2, 3, 4, 5, 6,
    8, 9, 0, 1, 2, 3, 4, 5, 6, 7,
    9, 0, 1, 2, 3, 4, 5, 6, 7, 8
  ),
  nrow = 20, ncol = 10, byrow = TRUE
)
rownames(data) <- paste0("Row_", 1:20)
colnames(data) <- paste0("Col_", 1:10)
df <- data.frame(
  ColumnName = colnames(data),
  Group = rep(c("Group1", "Group2"), 10),
  stringsAsFactors = FALSE
)

data_out <- data.frame(
  Group1 = c(20, 25, 20, 25, 20, 25, 20, 25, 20, 25, 
             20, 25, 20, 25, 20, 25, 20, 25, 20, 25),
  Group2 = c(25, 20, 25, 20, 25, 20, 25, 20, 25, 20, 
             25, 20, 25, 20, 25, 20, 25, 20, 25, 20),
  row.names = paste0("Row_", 1:20)
)


test_that("TeaCNV function returns expected output", {
  result <- TeaCNV::pseudo_bulk_v2(data,df)
  expect_equal(result, data_out) 
})
