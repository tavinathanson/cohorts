library('getopt')
library('deconstructSigs')

spec = matrix(c(
  'patient_input_dir', 'p'
))

result <- list()
counts <- list()

for(sample in ids){
  print(sample)
  data <- read.csv(sprintf("~/gc-bladder/post-bqsr_not-filtered_inputs/%s.input", sample), colClasses=c("Sample"="character"));
  
  data$Sample <- as.character(data$Sample)
  sigs.input <- mut.to.sigs.input(mut.ref = data, 
                                  sample.id = "Sample", 
                                  chr = "Chr", 
                                  pos = "Position", 
                                  ref = "Ref", 
                                  alt = "Alt")
  output <- whichSignatures(tumor.ref = sigs.input, 
                            signatures.ref = signatures.cosmic, 
                            sample.id = sample, contexts.needed = TRUE);
  
  # Add the sample name to a column named sample
  output$weights['Sample'] <- sample
  
  # Save the weights dataframe with a key/index of the `sample`
  result[[sample]] <- output$weights
  counts[[sample]] <- sigs.input
}

finalCounts = reshape::merge_all(counts)
counts_csv = cat(output_root, "_counts.csv")
sigs_csv = cat(output_root, "_signatures.csv")
write.csv(file = '~/gc-bladder/post-bqsr_not-filtered_outputs_counts.csv', finalCounts)
# Merge all the weight vectors into a single dataframe 
final = reshape::merge_all(result)
write.csv(file = '~/gc-bladder/post-bqsr_not-filtered_outputs_signatures.csv', final)