devtools::test(filter = "compute_LRM_Score") # 9 PASSED
devtools::test(filter = "smooth_expression") # 10 PASSED
devtools::test(filter = "extract_neighborhoods") # 36 PASSED
devtools::test(filter = "add_neighborhoods") # 22 PASSED
devtools::test(filter = "merge_NICHESObjects") # 29 PASSED
devtools::test(filter = "Explore_Radius") # 29 PASSED
devtools::test(filter = "extract_NICHESInputs_Seurat") # 18 PASSED
devtools::test(filter = "extract_NICHESInputs_AnnData") # 24 PASSED
devtools::test(filter = "create_Edgelist_Spatial") # 15 PASSED
devtools::test(filter = "create_Edgelist_PseudoBulk") # 18 PASSED
devtools::test(filter = "create_Edgelist_SamplingByType") # 30 PASSED
devtools::test(filter = "compute_CellToCell") # 20 PASSED
devtools::test(filter = "compute_PseudoBulk_Connectivity") # 18 PASSED
devtools::test(filter = "create_NICHESObject") # 20 PASSED
devtools::test(filter = "load_LRM_database") # 11 PASSED
devtools::test(filter = "add_cell_meta") # 6 PASSED
devtools::test(filter = "add_segmentation") # 12 PASSED
devtools::test(filter = "Aggregate_CellToCell_Edges") # 16 PASSED
devtools::test(filter = "CellToCell_InSitu") # 22 PASSED, 1 SKIP - intentionally
devtools::test(filter = "NICHESObject_methods") # 28 PASSED

devtools::test()

devtools::document()
devtools::check()

check_result <- capture.output(devtools::check())
writeLines(check_result, "./working-scripts-NW/check_log.txt")


devtools::load_all()
lapply(c("ncomms8866_orig","ncomms8866_human","ncomms8866_mouse","ncomms8866_pig","ncomms8866_rat"),
       function(nm) { d <- get(nm); cat(nm, ":", nrow(d), "rows\n"); print(names(d)); cat("\n") })