library(devtools)

## 1. make package directory
# create_package()
#
# 2. make corresponding git
# use_git()
#
# 3. Write function
# strsplit1 <- function(x, split) { strsplit(x, split = split)[[1]] }
#
# 4. Save function
# use_r("can.be.numeric")
#check_numeric  
#categorize_tertiles
#constru_single_continuous
#constru_single_Mclust 
#constru_single
#constru
#
# 5. Load function
# load_all()
# install()
#
#
# 6. Dependencies
#use_package("parallel")
#use_package("survival")
#use_package("survminer")
#use_package("mclust")
#
# 8. License
# use_mit_license()
# . Edit description
# export(strsplit1)
# use_readme_rmd()
#
# . readme
# use_readme_rmd()
#
# . build namespace
document()
# . check
# check()

#install.packages(".", repos = NULL, type="source")
#library("CONSTRU")

# 
library(markdown)
markdownToHTML('README.md', 'README.html')

