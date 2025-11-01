library('mindr')
library(data.tree)

## Source document
input_txt <- readLines("/Users/msbr/GitHub/NICHESMethods/markdown/RuleSetOrganization with Patterns.md", encoding = "unknown")

## Convert to mind map text, markdown outline, R script, and HTML widget ####
mm_output <- mindr::mm(input_txt, output_type = c("mindmap", "markdown", "R", "widget"))
mm_output

## Extract headings using regex
headings <- str_extract_all(mm_output$markdown, "^#+ .*$")
headings <- unlist(headings)
headings

## Create data frame structure
df <- data.frame(Basis=character(),
                 Category = character(),
                 Class = character(),
                 Family = character(),
                 Pattern = character())

## 1. For loop construction to build data frame
for(i in 1:length(headings)){
  level <- str_count(headings[i], "#") # Count the number of hashtags
  label <- str_trim(str_remove(headings[i], "^#+ ")) # Remove the hashes
  if(level==1){basis.stash <- label}
  if(level==2){category.stash <- label}
  if(level==3){class.stash <- label}
  if(level==4){family.stash <- label}
  if(level==5){pattern.stash <- label}
  if(level==5){
    df[i,] <- c(basis.stash,category.stash,class.stash,family.stash,pattern.stash)
  }
}
# 2. Remove the NAs induced by the for-loop construction above
df <- df[!is.na(df$Basis),]

## Inspect result
View(df)
