library(data.table)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr = data.frame( fread( file.path(input_dir, "EXPR.txt.gz") , stringsAsFactors=FALSE  , sep="\t" , dec=',') )

# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# 
# genes <- sapply( expr[,1] , function(x){ unlist( strsplit( x , "." , fixed=TRUE ))[1] })
# names(genes) = expr[,1]
# 
# ensembl <- as.matrix( getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
#               values = genes, mart= mart) )
# rownames(ensembl) = ensembl[ , "ensembl_gene_id"]
# 
# ensembl = ensembl[ intersect( genes , ensembl[ , "ensembl_gene_id" ] ) , ]
# genes = genes[ genes %in% intersect( genes , ensembl[ , "ensembl_gene_id" ] ) ]
# 
# expr = expr[ expr[,1] %in% names(genes) , ]
# rownames(expr)  = ensembl[ sapply( expr[,1] , function(x){ unlist( strsplit( x , "." , fixed=TRUE ))[1] }) , "hgnc_symbol" ]
# expr = expr[,-1]
# expr = expr[!( rownames(expr) %in% "" ), ]
# 
# rid = rownames(expr)
# cid = colnames(expr)
# expr = apply(apply(expr,2,as.character),2,as.numeric)
# colnames(expr) = cid
# rownames(expr) = rid

gene_id  <- expr$gene_id
expr = data.frame(apply(apply(expr,2,as.character),2,as.numeric))
expr$gene_id <- gene_id
#############################################################################
#############################################################################
## Remove duplicate genes
expr_uniq <- expr[!(expr$gene_id %in% expr[duplicated(expr$gene_id),]$gene_id), ]
rownames(expr_uniq) <- expr_uniq$gene_id
expr_uniq <- expr_uniq[, colnames(expr_uniq)[colnames(expr_uniq) != 'gene_id']]

expr_dup <- expr[expr$gene_id %in% expr[duplicated(expr$gene_id),], ]

# expr_uniq <- expr[!(rownames(expr)%in%rownames(expr[duplicated(rownames(expr)),])),]
# expr_dup <- expr[(rownames(expr)%in%rownames(expr[duplicated(rownames(expr)),])),]

if(length(expr_dup) > 0){
  expr_dup <- expr_dup[order(rownames(expr_dup)),]
}
id <- unique(expr_dup$gene_id)

expr_dup.rm <- NULL
# names <- NULL
for(j in 1:length(id)){
	tmp <- expr_dup[which(expr_dup$gene_id %in% id[j]), ]
	tmp <- tmp[1, ]
	rownames(tmp) <- tmp$gene_id
	tmp <- tmp[, -1]
	# tmp.sum <- apply(tmp,1,function(x){sum(as.numeric(as.character(x)),na.rm=T)})
	# tmp <- tmp[which(tmp.sum %in% max(tmp.sum,na.rm=T)), ]z
	if( !is.null(dim(tmp)) ){
	  expr_dup.rm <- rbind(expr_dup.rm,tmp) 
	  # names <- c(names,names(tmp.sum)[1])
	}   
}
expr <- rbind(expr_uniq,expr_dup.rm)
# rownames(expr) <- c(rownames(expr_uniq),names)
expr = expr[sort(rownames(expr)),]
#############################################################################
#############################################################################

colnames(expr) = sapply( colnames(expr) , function(x){ unlist( strsplit( x , "_T" , fixed=TRUE))[1] } )


case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
expr = log2( expr[ , case[ case$expr %in% 1 , ]$patient ] + 1 )

write.table( expr , file=file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
