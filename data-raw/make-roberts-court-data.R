## I do not like strings as factors
options(stringsAsFactors = FALSE)
## Read in the raw SCDB data. I will not keep this in the repo,
## so it must be downloaded from http://supremecourtdatabase.org/
scdb <- read.csv("SCDB_2019_01_justiceCentered_Citation.csv")
## We only need a few columns (users can connect to the rest of the variables
## since they know the caseId if they're curious about something)
columns_to_keep <- c("caseId", "lexisCite", "term", "justiceName", "vote")
## We only want Roberts Court cases that are not evenly divided or
## per curiam decisions without oral argument
rows_to_keep <- scdb$chief == "Roberts" & !scdb$decisionType %in% c(2, 5)
## Now we can subset the data and write it out
roberts_court_scdb_subset <- scdb[rows_to_keep, columns_to_keep]
## Notice this assumes your working directory is bggum/data-raw/
write.csv(roberts_court_scdb_subset,
          file = "../vignettes/roberts_court.csv",
          row.names = FALSE)
