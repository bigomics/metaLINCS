library(devtools)

## Build vignettes and manual
#usethis::use_vignette()
devtools::document()          ## populates man folder with Rd
devtools::build_vignettes()   ## build vignette
devtools::build_manual(path="../doc")    ## populates doc folder

## update docs, checks DESCRIPTION, creates vignettes
devtools::check(manual=TRUE)
# without running the examples
##devtools::check(cleanup = FALSE,args = c('--no-examples'),manual = TRUE,path = getwd())

