args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t")

clin = cbind( clin[ , c( "patient_id" , "sex" , "age" , "best_RECIST" , "os_days" , "os_censor" , "pfs_days" , "pfs_censor" ) ] , "Kidney" , "PD-1/PD-L1" , NA , NA , NA , NA , NA , NA )
colnames(clin) = c( "patient" , "sex" , "age" , "recist" , "t.os"  ,"os","t.pfs", "pfs" , "primary" , "drug_type" , "histo" , "stage" , "dna" , "rna" , "response.other.info" , "response" )

clin$sex = ifelse(clin$sex %in% "FEMALE" , "F" , "M")

clin$t.os = clin$t.os/30.5
clin$t.pfs = clin$t.pfs/30.5
clin$response = Get_Response( data=clin )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
clin$rna[ clin$patient %in% case[ case$expr %in% 1 , ]$patient ] = "tpm"
clin$dna[ clin$patient %in% case[ case$snv %in% 1 , ]$patient ] = "wes"

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

