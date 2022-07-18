library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_excel_functions.R")

# SNV.txt.gz
snv <- read_and_format_excel(
  input_path=file.path(work_dir, 'aan5951_tables1.xlsx'),
  sheetname='1C WES MAF (N=35)'
)
gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
write.table( snv , file=gz , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
close(gz)

# EXPR.txt.gz
expr <- read_and_format_excel(
  input_path=file.path(work_dir, 'aan5951_tables8.xlsx'),
  sheetname='temp3.csv'
)
gz <- gzfile(file.path(work_dir, 'EXPR.txt.gz'), "w")
write.table( expr , file=gz , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
close(gz)

# CLIN.txt
clin <- read_and_format_excel(
  input_path=file.path(work_dir, 'aan5951_tables1.xlsx'),
  sheetname='1B Discovery Clinical (N=35)'
)

clin_added <- read_and_format_excel(
  input_path=file.path(work_dir, 'aan5951_tables2.xlsx'),
  sheetname='S2B Val Clinical'
)

additional_patients <- colnames(expr)[str_detect(colnames(expr), 'PD')]
additional_patients <- str_replace_all(additional_patients, 'T', '')
clin_added_patients <- unlist(lapply(additional_patients, function(additional_patient){
  clin_added$patient_id[str_detect(clin_added$patient_id, additional_patient)]
}))
clin_added <- clin_added[clin_added$patient_id %in% clin_added_patients, ]
colnames(clin_added)[colnames(clin_added) == 'best_recist'] <- 'best_RECIST'
clin_added <- clin_added[, colnames(clin_added)[colnames(clin_added) %in% colnames(clin)]]
missing_cols <- colnames(clin)[!colnames(clin) %in% colnames(clin_added)]
for(col in missing_cols){
  clin_added[col] <- NA
}
clin_added <- clin_added[, colnames(clin)]
clin_added$patient_id <- str_replace_all(clin_added$patient_id, 'RCC-', '')
clin_added$patient_id <- paste0(clin_added$patient_id, 'T')
clin <- rbind(clin, clin_added)
write.table( clin , file=file.path(work_dir, 'CLIN.txt') , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
