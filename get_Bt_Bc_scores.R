
library(data.table)

args <- commandArgs(TRUE)

fname <- as.character(args[1])

f <- fread(fname, sep = "," , header = FALSE)

setkey(f, V1, V2)

match_idx <- cbind(unique(f[, "V1"]), unique(f[, "V1"]))

match_vals <- f[match_idx, "V5"]

write(match_vals, file = "match_vals.txt", sep = "\n")

