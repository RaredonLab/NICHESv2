test_that("build_edge_meta: adds Sending.* and Receiving.* columns for each meta col", {
  cell.meta <- data.table::data.table(
    cell      = c("c1","c2","c3"),
    cell.type = c("A","B","A"),
    x         = c(0.1, 0.2, 0.3),
    y         = c(1.1, 1.2, 1.3)
  )
  data.table::setkey(cell.meta, cell)

  edge.list <- data.table::data.table(
    edge          = c("c1|c2","c2|c3"),
    SendingCell   = c("c1","c2"),
    ReceivingCell = c("c2","c3")
  )

  em <- build_edge_meta(edge.list, cell.meta)

  expect_true(all(c("Sending.cell.type","Receiving.cell.type",
                    "Sending.x","Receiving.x",
                    "Sending.y","Receiving.y") %in% names(em)))
  expect_equal(em$Sending.cell.type,   c("A","B"))
  expect_equal(em$Receiving.cell.type, c("B","A"))
  expect_equal(em$Sending.x,   c(0.1, 0.2))
  expect_equal(em$Receiving.x, c(0.2, 0.3))
})

test_that("build_edge_meta: Joint.* column built only for categorical cols", {
  cell.meta <- data.table::data.table(
    cell      = c("c1","c2"),
    cell.type = c("A","B"),     # categorical → Joint expected
    x         = c(0.5, 0.7)     # numeric → no Joint
  )
  data.table::setkey(cell.meta, cell)

  edge.list <- data.table::data.table(
    edge = "c1|c2",
    SendingCell = "c1",
    ReceivingCell = "c2"
  )

  em <- build_edge_meta(edge.list, cell.meta)
  expect_true("Joint.cell.type" %in% names(em))
  expect_false("Joint.x" %in% names(em))
  expect_equal(em$Joint.cell.type, "A|B")
})

test_that("build_edge_meta: empty cell.meta (only 'cell' col) returns edge.list unchanged", {
  cell.meta <- data.table::data.table(cell = c("c1","c2"))
  data.table::setkey(cell.meta, cell)

  edge.list <- data.table::data.table(
    edge = "c1|c2", SendingCell = "c1", ReceivingCell = "c2"
  )
  em <- build_edge_meta(edge.list, cell.meta)
  expect_equal(sort(names(em)), sort(c("edge","SendingCell","ReceivingCell")))
  expect_equal(nrow(em), 1L)
})

test_that("build_edge_meta: result is keyed on 'edge'", {
  cell.meta <- data.table::data.table(cell = c("c1","c2"), x = c(1, 2))
  data.table::setkey(cell.meta, cell)
  edge.list <- data.table::data.table(
    edge = c("c1|c2","c2|c1"),
    SendingCell = c("c1","c2"),
    ReceivingCell = c("c2","c1")
  )
  em <- build_edge_meta(edge.list, cell.meta)
  expect_equal(data.table::key(em), "edge")
})

test_that("build_edge_meta: factor columns also generate Joint.* columns", {
  cell.meta <- data.table::data.table(
    cell  = c("c1","c2"),
    group = factor(c("X","Y"))
  )
  data.table::setkey(cell.meta, cell)
  edge.list <- data.table::data.table(
    edge = "c1|c2", SendingCell = "c1", ReceivingCell = "c2"
  )
  em <- build_edge_meta(edge.list, cell.meta)
  expect_true("Joint.group" %in% names(em))
  expect_equal(em$Joint.group, "X|Y")
})
