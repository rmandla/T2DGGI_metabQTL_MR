library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)
library(MRPRESSO)

args = commandArgs(trailingOnly=TRUE)

# read the outcome GWAS
protname <- args[1]
dataset = args[2]
exposure_type = args[3]
outcome_type = args[4]
N_exposure = as.integer(args[5])
N_outcome = as.double(args[6])
output_dir = args[7]
filter = args[8]
filter_name = ''

direction = paste0(exposure_type,'_to_',outcome_type)
outcome_path = paste0(output_dir,'_',protname,'_',dataset,'_',outcome_type,'_harmonized.txt')
exp_path = paste0(output_dir,'_',protname,'_',dataset,'_',exposure_type,'_harmonized.txt')
outcome_GWAS <- vroom(outcome_path)

#exposure
#exp_path = paste0(protname,'_',dataset,'_chrALL_clumped_withmeta.txt')

exp_dat <- read_exposure_data(filename = exp_path,
                              sep='\t',
                              snp_col='rsid',
                              beta_col = 'Effect',
                              se_col = 'SE',
                              effect_allele_col = 'eAllele',
                              other_allele_col = 'oAllele',
                              eaf_col = 'AF',
                              pval_col = 'P')

# read the outcome
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$rsid, type="outcome", snp_col="rsid", beta_col="Effect", se_col = "SE", eaf_col = "AF",
                                 effect_allele_col = "eAllele", other_allele_col = "oAllele", pval_col = "P", chr_col = "chr", pos_col = "pos")
# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exp_dat, outcome_dat=formatted_outcome)
# exclude rare variants
exp_dat_outcome %<>% filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))
exp_dat_outcome$samplesize.exposure <- N_exposure
exp_dat_outcome$samplesize.outcome <- N_outcome

if ( filter=='fstat') {
exp_dat_outcome$Fstat <- (exp_dat_outcome$beta.exposure)^2/(exp_dat_outcome$se.exposure)^2
exp_dat_outcome = exp_dat_outcome[exp_dat_outcome$Fstat>10,]
print('Running with Fstat filtering')
filter_name = '_fstat'
} else if ( filter=='steig') {
exp_dat_outcome <- steiger_filtering(exp_dat_outcome)
exp_dat_outcome <- subset(exp_dat_outcome, steiger_dir & steiger_pval<0.05)
filter_name = '_steig'
print('Running with Steiger filtering')
} else if ( filter=='both') {
exp_dat_outcome$Fstat <- (exp_dat_outcome$beta.exposure)^2/(exp_dat_outcome$se.exposure)^2
exp_dat_outcome = exp_dat_outcome[exp_dat_outcome$Fstat>10,]
exp_dat_outcome <- steiger_filtering(exp_dat_outcome)
exp_dat_outcome <- subset(exp_dat_outcome, steiger_dir & steiger_pval<0.05)
print('Running with Fstat and Steiger filtering')
filter_name = '_fstat_steig'
} else {
print('Running without filtering')
}

if (nrow(exp_dat_outcome)>1) {
# mr result
mr_results <- mr(exp_dat_outcome)
mr_name <- paste0(output_dir, '_',protname, ".",dataset,'.',direction,'.',filter_name,".mr.txt")
write.table(mr_results, file=mr_name, sep = '\t', quote = F)

# odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, '_',protname, ".",dataset,'.',direction,'.',filter_name,".or.txt")
write.table(OR, file=OR_name, sep = '\t', quote = F, row.names = F)

# scatter plot
pdf_name <- paste0(output_dir, '_',protname, ".",dataset,'.',direction,'.',filter_name,".scatter.pdf")
pdf(pdf_name, width = 10)

# Try plotting and ensure it doesn't throw an error
tryCatch({
  print("Generating scatter plot...")  # Debugging message
  plot_result <- mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0, 1) + theme_minimal()
  print("Scatter plot generated successfully.")  # Debugging message

  # Print the plot to the file
  print(plot_result)

}, error = function(e) {
  cat("Error generating scatter plot: ", conditionMessage(e), "\n")
})

# Ensure the pdf device is properly closed
dev.off()

# horizontal pleiotropy
tryCatch({
  pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
  pleio_name <- paste0(output_dir,'_',protname, ".",dataset,'.',direction,'.',filter_name,".pleio.txt")
  write.table(pleio_res, file=pleio_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# hetero test
tryCatch({
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$isquared <- abs(100*(hetero_res$Q - hetero_res$Q_df)/hetero_res$Q)  # I2 = 100%×(Q - df)/Q
  hetero_name <- paste0(output_dir,'_',protname, ".",dataset,'.',direction,'.',filter_name,".hetero.txt")
  write.table(hetero_res, file=hetero_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# steiger
print('RUNNING STEIGER')
steiger <- directionality_test(exp_dat_outcome)
steiger_name <- paste0(output_dir, '_',protname, ".",dataset,'.',direction,'.',filter_name,".steiger.txt")
write.table(steiger, file=steiger_name, sep = '\t', quote = F, row.names = F)

# MrPresso
tryCatch({
  mrpresso_res <- mr_presso(BetaExposure = "beta.exposure",
                            BetaOutcome = "beta.outcome",
                            SdOutcome = "se.outcome",
                            SdExposure = "se.exposure",
                            OUTLIERtest = TRUE,
                            DISTORTIONtest = TRUE,
                            data = exp_dat_outcome,
                            NbDistribution = 1000,
                            SignifThreshold = 0.05)
  mrpresso_name <- paste0(output_dir, '_', protname, '.',dataset,'.',direction,'.',filter_name, ".mrpresso.txt")

  # add global_rss, global_pval, distortion_indices, distortion_coef, distortion_pval
  mrpresso_df <- as.data.frame(mrpresso_res$`Main MR results`)
  mrpresso_df %<>% mutate(global_rss = mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs,
    global_pval =  mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue)
    #distortion_indices = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,
    #distortion_coef = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
    #distortion_pval = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue)

  write_tsv(mrpresso_df, file =  mrpresso_name)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
