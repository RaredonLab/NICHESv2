FILE="R/create_NICHESObject.R"
grep -c "edge.filter.mode" $FILE   # expect 2 (comment + call)
grep -c "remove.autocrine" $FILE   # expect 2 (comment + call)
grep -c "@param edge"      $FILE   # expect 0
grep -c "@param remove"    $FILE   # expect 0