library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

results <- read.table(args[1], sep  = "\t", header =T, row.names=NULL)

sum_results <- results %>%
  group_by(Sample, Fusion) %>%
  summarise(JunctionReads = sum(JunctionReads),
            SpanningFragments = sum(SpanningFragments),
            UniqueStarts = sum(UniqueStarts))

filt_results <- sum_results %>%
  filter(JunctionReads >= 9) %>% 
  filter(UniqueStarts >= 4)

print(paste0("writing: ", args[1], " summed"))
write.csv(sum_results, file=paste0(args[1], "_summed.csv"), quote = F, row.names=F)
print(paste0("writing: ", args[1], " filtered"))
write.csv(filt_results, file=paste0(args[1], "_filtered.csv"), quote = F, row.names=F)
