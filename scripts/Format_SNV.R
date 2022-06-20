library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

snv = as.data.frame( fread( file.path(input_dir, "SNV.txt.gz") , stringsAsFactors=FALSE , sep="\t" ))

data = cbind( snv[ , c("Start_position" , "Tumor_Sample_Barcode" , "Hugo_Symbol", "Variant_Classification" ) ] ,
				ifelse( snv[ , "Tumor_Seq_Allele1" ] %in% "-" , "" , snv[ , "Tumor_Seq_Allele1" ] ) , 
				ifelse( snv[ , "Tumor_Seq_Allele2" ] %in% "-" , "" , snv[ , "Tumor_Seq_Allele2" ] ) , 
				sapply( snv[ , "Chromosome" ] , function(x){ paste( "chr" , x , sep="" ) } )
			)

colnames(data) = c( "Pos" , "Sample" , "Gene" , "Effect" , "Ref" , "Alt" , "Chr" )

data$Ref = ifelse( data$Ref %in% "-" , "" , data$Ref )
data$Alt = ifelse( data$Alt %in% "-" , "" , data$Alt )

data = cbind ( data , apply( data[ , c( "Ref", "Alt" ) ] , 1 , function(x){ ifelse( nchar(x[1]) != nchar(x[2]) , "INDEL", "SNV") } ) )
colnames(data) = c( "Pos" , "Sample" , "Gene" , "Effect" , "Ref" , "Alt" , "Chr" , "MutType"  )


data$Sample = sapply( data$Sample , function(x){ unlist( strsplit( x , "-Tumor" , fixed=TRUE))[1] } )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
data = data[ data$Sample %in% case[ case$snv %in% 1 , ]$patient , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]


write.table( data , file= file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
