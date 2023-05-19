if(!require("devtools")) install.packages("devtools")
if(!require("roxygen2")) install.packages("roxygen2")
if(!require("testthat")) install.packages("testthat")
usethis::create_package("./Rpackage/HBindex")
