library(devtools)

## update docs, checks DESCRIPTION, creates vignettes
devtools::check(manual=TRUE)
# without running the examples
##devtools::check(cleanup = FALSE,args = c('--no-examples'),manual = TRUE,path = getwd())

## Build vignettes and manual
devtools::build_vignettes()   ## populates doc folder
devtools::build_manual(path="../doc")
