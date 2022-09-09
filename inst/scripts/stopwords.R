# save tm::stopwords() in a data object to avoid another package dependency
library(tm)
library(devtools)
STOPWORDS = tm::stopwords()
usethis::use_data(STOPWORDS)
