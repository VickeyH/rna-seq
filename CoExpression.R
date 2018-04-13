# Script Name: CoExpression.R 
# Function:  Correlation Analysis
# Version: v1.0
# Usage:
#       source("A_B_Correlation.R")
#       A_B_Correlation( r_cutoff, p_cutoff, exp, gene_filter, 
#                        sample_filter,out)
# Author: yujingzhang@capitalbio.com
# Data: 2014-06-26

# Sub Function
get_sig_cor <- function( x, r_cutoff = 0.99,p_cutoff = 0.05,genelist,out,filename)  {

	r <- cor(as.matrix(x), method = "pearson")
	p <- matrix(NA, ncol(x), ncol(x))
	rownames(p) <- colnames(x)
	colnames(p) <- colnames(x)

	for (i in 1:ncol(x)) {
		for (j in 1:ncol(x)) {
			p[i,j] <- cor.test(x[,i], x[,j], method="pearson")$p.value
		}
	}

	if (r_cutoff >= 0) {
		sig <- abs(r) >= r_cutoff & p <= p_cutoff
	}else if (r_cutoff < 0) {
	    sig <- r < r_cutoff & p <= p_cutoff
	}

	sig[ lower.tri(sig, diag=T) ] <- FALSE

	ind <-  which( sig, arr.ind=F)
	xyind <-  which( sig, arr.ind=T)

	sig_cor <- data.frame(Source = rownames(r)[xyind[,1]],
	                      Target = colnames(r)[xyind[,2]],
						  Correlation = r[ind],
						  P.value = p[ind] )
	if(is.data.frame(genelist)){
		Source.GeneSymbol <- sapply(as.character(sig_cor$Source), function(x)
					            ifelse(is.na(genelist[x,"GeneSymbol"]) |
							           genelist[x,"GeneSymbol"] == "---",
						 	           x, as.character(genelist[x,"GeneSymbol"])) )
		Target.GeneSymbol <- sapply(as.character(sig_cor$Target), function(x)
					            ifelse(is.na(genelist[x,"GeneSymbol"]) |
							           genelist[x,"GeneSymbol"] == "---",
							           x, as.character(genelist[x,"GeneSymbol"])) )

		if (is.na(file.info(out)[2])) dir.create(out)
		write.table( cbind(sig_cor, Source.GeneSymbol, Target.GeneSymbol),
	             filename,col.names=T, row.names=F, sep="\t", quote=F)

	}else{
		if (is.na(file.info(out)[2])) dir.create(out)
		write.table( sig_cor,filename,col.names=T, row.names=F, sep="\t", quote=F)
	}
	

}


# Main Function
CoExpression <- function(r_cutoff,
                            p_cutoff,
                            exp_sgl,
			    gene_filter,
			    sample_filter,
			    out,
			    name1,
			    name2) {
	
	r_cutoff <- as.numeric(r_cutoff)
	p_cutoff <- as.numeric(p_cutoff)

	A <- read.table( exp_sgl, comment.char = "",quote="",
	                 header=T, row.names=1, sep="\t", check.names=F)
	
	# only use genes in gene_filter file
	genelist <- dir( pattern= gene_filter )
	if ( length(genelist)==1 ) {
		genelist <- read.table(gene_filter, check.names=F, 
		                       comment.char= "", header=T, sep="\t")
		if(grep("ProbeName",colnames(genelist))){
		    rownames(genelist) <- genelist[,"ProbeName"]
		}
		else if(grep("tracking_id",colnames(genelist))){
		    rownames(genelist) <- genelist[,"tracking_id"]
		}
		A <- A[ rownames(A) %in% rownames(genelist), ]
	}

	# only use samples in sample_filter file
	samplelist <- dir( pattern= sample_filter )
	if ( length(samplelist)==1 ) {
		samplelist <- read.table( sample_filter , check.names=F,
		                         comment.char = "", header=T, sep="\t")
		grp <- samplelist[,"group"]
		if(length(unique(grp))==1){
			EXP <- A[ , colnames(A) %in% samplelist[,"sampleID"] ]
			EXP_sig_cor <- get_sig_cor( t(EXP), r_cutoff = r_cutoff,
	                            p_cutoff = p_cutoff,genelist,out,"Co-Expression.xls" )
		}else{
			loc1 <- grp==1
			loc2 <- grp==2
			TEXP <- A[ , colnames(A) %in% samplelist[loc1,"sampleID"] ]
			
			NEXP <- A[ , colnames(A) %in% samplelist[loc2,"sampleID"] ]
			# compute significantly correlated case A-B pairs
			get_sig_cor( t(TEXP), r_cutoff = r_cutoff,
	                           p_cutoff = p_cutoff,genelist,out,paste(name1,"Co-Expression.xls",sep=".") )
	
	
	
			# compute significantly correlated control A-B pairs
			get_sig_cor( t(NEXP), r_cutoff = r_cutoff,
	                            p_cutoff = p_cutoff,genelist,out,paste(name2,"Co-Expression.xls",sep=".") )

		}
	}else{
		
		 get_sig_cor( t(A), r_cutoff = r_cutoff,
	                     p_cutoff = p_cutoff,genelist,out,"Co-Expression.xls" )
	}

	
	
}
