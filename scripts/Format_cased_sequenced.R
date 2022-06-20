library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
rownames(clin) = clin$patient_id

rna = as.matrix( fread( file.path(input_dir, "EXPR.txt.gz") , stringsAsFactors=FALSE  , sep="\t" , dec=',') )
rna = unique( sapply( colnames(rna)[-1] , function(x){ unlist( strsplit( x , "_T" , fixed=TRUE ) )[1] } ) )

patient = sort( unique( c( clin$patient_id , rna ) ) )

case = as.data.frame( cbind( patient , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) ) )
colnames(case) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = patient

case$snv = as.numeric( as.character( case$snv ) )
case$cna = as.numeric( as.character( case$cna ) )
case$expr = as.numeric( as.character( case$expr ) )

for( i in 1:nrow(case)){
	if( clin[ rownames(case)[i] ,]$wes %in% 1 ){
		case$snv[i] = 1
	}
	if( rownames(case)[i] %in% rna ){
		case$expr[i] = 1
	}
}

write.table( case , file=file.path(output_dir, "cased_sequenced.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )


