## Read in the Roberts Court data (working directory should be bggum/data-raw/)
options(stringsAsFactors = FALSE)
roberts_court <- read.csv("../vignettes/roberts_court.csv")
## Spread it into a response matrix after recoding the votes
roberts_court$vote <- ifelse(roberts_court$vote %in% c(2, 7), 0, 1)
library(dplyr)
library(tidyr)
library(bggum)
responses <- roberts_court %>%
    select(-caseId, -term) %>%
    spread(lexisCite, vote)
rownames(responses) <- responses$justiceName
responses <- as.matrix(responses[ , -1])
unanimous <- apply(responses, 2, function(x) length(unique(na.omit(x))) == 1)
responses <- responses[ , !unanimous]
## Tune hyperparameters
library(bggum)
set.seed(123)
proposal_sds <- tune_proposals(responses, 5000)
sapply(proposal_sds, mean)
set.seed(456)
temps <- tune_temperatures(responses, n_temps = 6, proposal_sds = proposal_sds)
round(temps, 2)
temps <- temps[1:6]
## Run two chains in parallel
library(parallel)
cl <- makeCluster(2, type = "FORK", outfile = "bggum-log2.txt")
clusterSetRNGStream(cl = cl, iseed = 789)
chains <- parLapplyLB(cl = cl, X = 1:2, fun = function(x) {
    ggumMC3(data = responses,
            sample_iterations = 50000,
            burn_iterations = 5000,
            proposal_sds = proposal_sds,
            temps = temps)
})
stopCluster(cl)
## Post-process
constraint <- which(rownames(responses) == "RBGinsburg")
processed_chains <- lapply(chains, post_process, constraint = constraint,
                           expected_sign = "-")
## Summarize posterior
posterior_summary <- summary(processed_chains)
## Save the posterior summary
saveRDS(posterior_summary, file = "../vignettes/posterior_summary.rds")
## Check convergence
library(coda)
convergence_stats <- gelman.diag(processed_chains)
write.csv(convergence_stats$psrf, file = "../vignettes/convergence.csv",
          row.names = FALSE)
## Write out Roberts' draws
roberts <- which(rownames(responses) == "JGRoberts")
iters <- nrow(chains[[1]])
idx <- seq(floor(iters / 1000), iters, floor(iters / 1000))
chain1_raw <- chains[[1]][idx, roberts]
chain2_raw <- chains[[2]][idx, roberts]
chain1_pp  <- processed_chains[[1]][idx, roberts]
chain2_pp  <- processed_chains[[2]][idx, roberts]
roberts_draws <- cbind(chain1_raw, chain2_raw, chain1_pp, chain2_pp)
write.csv(roberts_draws, file = "../vignettes/roberts_draws.csv",
          row.names = FALSE)

