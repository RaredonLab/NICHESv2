# create_NICHESObject.R
grep -c "edge.filter.mode" R/create_NICHESObject.R        # expect >= 4
grep -c "remove.autocrine" R/create_NICHESObject.R        # expect >= 3
grep -c "add_neighborhoods" R/create_NICHESObject.R       # expect >= 2
grep -c "return(" R/create_NICHESObject.R                 # expect 0 (no early returns remain)

# add_neighborhoods.R
grep -c "overwrite" R/add_neighborhoods.R                 # expect >= 5
grep -c "remove.autocrine.*attr" R/add_neighborhoods.R    # expect 1
grep -c "stale" R/add_neighborhoods.R                     # expect >= 2

# build_NICHESv2_example.R
grep -c "neighborhoods)" data-raw/build_NICHESv2_example.R  # expect 1

# test-create_NICHESObject.R
grep -c "Section 6" tests/testthat/test-create_NICHESObject.R  # expect 1