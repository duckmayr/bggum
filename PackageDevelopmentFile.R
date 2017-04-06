
## Package Development File Template ##
library(devtools)
library(roxygen2)
current.code <- as.package('ggumR')
load_all(current.code)
document(current.code)
# check(current.code)
# install(pkg=current.code, local=TRUE)
# build(current.code, path=getwd())
