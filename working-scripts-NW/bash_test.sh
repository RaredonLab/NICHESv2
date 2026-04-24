# Run from package root after applying all edits
FILE="R/NICHESObject_class.R"

grep -c "15-slot specification"                    $FILE  # expect 1
grep -c "15 constituent slots"                     $FILE  # expect 1
grep -c "aggregations           = list()"          $FILE  # expect 2  (signature + structure)
grep -c "aggregations           = aggregations"    $FILE  # expect 1  (structure() call)
grep -c "names(x\$aggregations)"                  $FILE  # expect 1  (print method)
grep -c "names(object\$aggregations)"              $FILE  # expect 1  (summary method)
grep -c "Aggregation slots dropped on subset"      $FILE  # expect 1  ([method)