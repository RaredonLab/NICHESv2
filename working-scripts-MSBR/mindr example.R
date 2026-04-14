library('mindr')
# Example 1: From Markdown to other outputs ####

## Source document ####
input <- system.file("examples/mindr-md.Rmd", package = "mindr")

## file.show(input) # Open the input file with the default program, if any
input_txt <- readLines(input, encoding = "UTF-8")

## Convert to mind map text, markdown outline, R script, and HTML widget ####
mm_output <- mm(input_txt, output_type = c("mindmap", "markdown", "R", "widget"))
mm_output

## Save the output texts as files ####

### mind map ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".mm")
writeLines(mm_output$mindmap, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### markdown outline ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".md")
writeLines(mm_output$markdown, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### R script ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".R")
writeLines(mm_output$r, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### Widget ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".html")
htmlwidgets::saveWidget(mm_output$widget, file = output)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

## Generate directory according to the source document ####
temp_dir <- file.path(tempdir(), "mindr")
mm_output <- mm(input_txt, output_type = "dir", root = "mindr", md_list = TRUE,
                md_braces = TRUE, md_bookdown = TRUE, dir_to = temp_dir)
# system2('open', temp_dir) # Open the generated directory unlink(temp_dir,
# recursive = TRUE) # remove the generated directory

## More arguments ####
mm_output <- mm(input_txt, output_type = c("mindmap", "markdown", "R", "widget"),
                root = "mindr", md_list = TRUE, md_braces = TRUE, md_bookdown = TRUE)
mm_output

# Example 2: From mind map to other outputs ####

## Source document ####
input <- system.file("examples/mindr-mm.mm", package = "mindr")

## file.show(input) # Open the input file with the default program, if any
input_txt <- readLines(input, encoding = "UTF-8")

## Convert markdown outline, R script, and HTML widget ####
mm_output <- mm(input_txt, output_type = c("markdown", "R", "widget"))
mm_output

## Save the output texts as files ####

### markdown outline ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".md")
writeLines(mm_output$markdown, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### R script ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".R")
writeLines(mm_output$r, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### Widget ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".html")
htmlwidgets::saveWidget(mm_output$widget, file = output)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

## Generate directory according to the source document ####
temp_dir <- file.path(tempdir(), "mindr")
mm_output <- mm(input_txt, output_type = "dir", root = "mindr", dir_to = temp_dir)
# system2('open', temp_dir) # Open the generatecd directory unlink(temp_dir,
# recursive = TRUE) # remove the generated directory

# Example 3: From R script to other outputs ####

## Source document ####
input <- system.file("examples/mindr-r.R", package = "mindr")

## file.show(input) # Open the input file with the default program, if any
input_txt <- readLines(input, encoding = "UTF-8")

## Convert to mind map text, markdown text, and HTML widget ####
mm_output <- mm(input_txt, output_type = c("mindmap", "markdown", "widget"))
mm_output

## Save the output texts as files ####

### mind map ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".mm")
writeLines(mm_output$mindmap, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### R markdown ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".Rmd")
writeLines(mm_output$markdown, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### Widget ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".html")
htmlwidgets::saveWidget(mm_output$widget, file = output)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

## Generate directory according to the source document ####
temp_dir <- file.path(tempdir(), "mindr")
mm_output <- mm(input_txt, output_type = "dir", root = "mindr", dir_to = temp_dir)
# system2('open', temp_dir) # Open the generated directory unlink(temp_dir,
# recursive = TRUE) # remove the generated directory

# Example 4: From directory to other outputs ####

## Source directory ####
input <- system.file(package = "mindr")

## Convert to mind map text, markdown outline, R script, and HTML widget ####
mm_output <- mm(input, output_type = c("mindmap", "markdown", "R", "widget"))
mm_output

## Save the output texts as files ####

### mind map ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".mm")
writeLines(mm_output$mindmap, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### markdown outline ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".md")
writeLines(mm_output$markdown, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### R script ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".R")
writeLines(mm_output$r, output, useBytes = TRUE)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

### Widget ####
output <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".html")
htmlwidgets::saveWidget(mm_output$widget, file = output)
# file.show(output) # Open the output file with the default program, if any
message("Input:  ", input, "\nOutput: ", output)
# file.remove(output) # remove the output file

## Clone the source directory ####
temp_dir <- file.path(tempdir(), "mindr")
mm_output <- mm(input, output_type = "dir", dir_to = temp_dir)
# system2('open', temp_dir) # Open the generated directory unlink(temp_dir,
# recursive = TRUE) # remove the generated directory

