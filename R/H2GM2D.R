#' HTGM2Ddriver
#'
#' @import minimalistGODB
#' @import GoMiner
#' @import HTGM
#' @import grDevices
#' @import stats
#' @importFrom gplots heatmap.2
#' @import jaccard
#'
#' @description driver to invoke GoMiner and HTGM2D, and compare the results
#'
#' @param dir character string full path name to the directory acting as result repository
#' @param geneList character vector of user-supplied genes of interest
#' @param ontologies character vector of 2 ontologies e.g. c("biological_process","cellular_component")
#' @param GOGOA3 return value of subsetGOGOA() 
#' @param enrichThresh  numerical acceptance threshold for enrichment passed to GoMiner
#' @param countThresh numerical acceptance threshold for gene count passed to GoMiner
#' @param fdrThresh numerical acceptance threshold for fdr passed to GoMiner
#' @param nrand integer number of randomizations passed to GoMiner
#'
#' @examples
#' \dontrun{
#' # this example takes too long to run, and
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # This example is given in full detail in the package vignette.
#' # You can generate GOGOA3.RData using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO

#' # load("~/GODB_RDATA/GOGOA3.RData")
#' geneList<-cluster52
#' ontologies<-c("biological_process","cellular_component")
#' dir<-tempdir()
#' HTGM2Ddriver(dir,geneList,ontologies,GOGOA3,enrichThresh=2,countThresh=5,fdrThresh=0.10,nrand=100)
#' }
#' 
#' @return returns no value, but saves hyperlinked SVG heatmap files to a results directory
#' 
#' @export
HTGM2Ddriver<-
  function(dir,geneList,ontologies,GOGOA3,enrichThresh=2,countThresh=5,fdrThresh=0.10,nrand=100) {
    stamp<-gsub(":","_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
    title<-"HTGM2Dresults"
    subd<-sprintf("%s/%s_%s",dir,title,stamp)
    dir.create(subd)
    
    # first run regular GoMiner as basis for comparison later
    l<-list()
    for(ontology in ontologies)
      l[[ontology]]<-GoMiner(title=ontology,subd,geneList,GOGOA3,ontology,enrichThresh,countThresh,fdrThresh,nrand)
    #x_l<-l
    #save(x_l,file="data/x_l.RData")
    
    # next run HTGM2D
    mat<-HTGM2D(subd,geneList,ontologies,GOGOA3)
    #x_mat<-mat
    #save(x_mat,file="data/x_mat.RData")
    
    # finally compare regular GoMiner and HTGM2D results
    # hope is that HTGM2D retrieves some significant categories that are not retrieved by regular GoMiner
    compareGoMinerHTGM2D(subd,mat,l,ontologies)
  }

#' compareGoMinerHTGM2D
#'
#' @description compare the results of GoMiner and HTGM2D
#'
#' @param subd character string full path name to the output subdirectory
#' @param mat return value of Jaccard()
#' @param l of return values of GoMiner()
#' @param ontologies character vector of 2 ontologies e.g. c("biological_process","cellular_component")
#'
#' @examples
#' ontologies<-c("biological_process","cellular_component")
#' #load("data/x_l.Rdata")
#' #load("data/x_mat.Rdata")
#' subd<-tempdir()
#' compareGoMinerHTGM2D(subd,x_mat,x_l,ontologies)
#' 
#' @return returns no value, but saves files that list category difference between GoMiner and HTG2D
#' 
#' @export
compareGoMinerHTGM2D<-
  function(subd,mat,l,ontologies) {

    sink(sprintf("%s/%s",subd,"compareGoMinerHTGM2D.txt"))
    print(ontologies[1])
    sd1<-setdiff(unique(mat[,"catr"]),l[[ontologies[1]]]$thresh[,"Row.names"])
    print("ONLY IN HTGM2D:")
    print(sd1)
    
    sd2<-setdiff(l[[ontologies[1]]]$thresh[,"Row.names"],unique(mat[,"catr"]))
    print("ONLY IN GoMiner:")
    print(sd2)
    
    overlap<-intersect(unique(mat[,"catr"]),l[[ontologies[1]]]$thresh[,"Row.names"])
    print("OVERLAP:")
    print(overlap)
      
    print(ontologies[2])
    sd1<-setdiff(unique(mat[,"catc"]),l[[ontologies[2]]]$thresh[,"Row.names"])
    print("ONLY IN HTGM2D:")
    print(sd1)
    
    sd2<-setdiff(l[[ontologies[2]]]$thresh[,"Row.names"],unique(mat[,"catc"]))
    print("ONLY IN GoMiner:")
    print(sd2)
    
    overlap<-intersect(unique(mat[,"catc"]),l[[ontologies[2]]]$thresh[,"Row.names"])
    print("OVERLAP:")
    print(overlap)
    sink()
  }

#' HTGM2D
#' 
#' @description run 2D version of GoMiner
#' 
#' @param dir character string full path name to the directory acting as result repository
#' @param geneList character vector of user-supplied genes of interest
#' @param ontologies character vector of 2 ontologies e.g. c("biological_process","cellular_component")
#' @param GOGOA3 return value of subsetGOGOA()  
#' 
#' @examples
#' \dontrun{
#' # this example takes too long to run, and
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # This example is given in full detail in the package vignette.
#' # You can generate GOGOA3.RData using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO

#' # load("~/GODB_RDATA/GOGOA3.RData")
#' subd<-tempdir()
#' geneList<-cluster52
#' ontologies<-c("biological_process","cellular_component")
#' mat<-HTGM2D(subd,geneList,ontologies,GOGOA3)
#' }
#' 
#' @return returns the return value of Jaccard()
#' 
#' @export
HTGM2D<-
  function(dir,geneList,ontologies,GOGOA3) {
    m1<-catGenes(geneList,GOGOA3,ontologies[1])
    #x_m1<-m1
    #save(x_m1,file="data/x_m1.RData")
    m2<-catGenes(geneList,GOGOA3,ontologies[2])
    #x_m2<-m2
    #save(x_m2,file="data/x_m2.RData")
    
    mat<-Jaccard(dir,m1,m2)
    #x_jmat<-mat
    #save(x_jmat,file="data/x_jmat.RData")
    jHeatMap<-JaccardHeatMap(dir,mat)
    return(mat)
  }

#' catGenes
#' 
#' @description match up genes in gene list with categories in GOGOA3 database
#'
#' @param geneList character vector of user-supplied genes of interest
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ontology c("molecular_function","cellular_component","biological_process")
#' 
#' @examples
#' #load("data/GOGOA3small.RData")
#' geneList<-cluster52
#' m1<-catGenes(geneList,GOGOA3small,"biological_process")
#' 
#' @return returns a matrix of 1's and 0's indicating the presence or absence of gene-category pairs
#' 
#' @export
catGenes<-
	function(geneList,GOGOA3,ontology) {
		DB<-GOGOA3$ontologies[[ontology]]
		CATS<-unique(DB[,"GO_NAME"])
		
		m<-matrix(0,nrow=length(CATS),ncol=length(geneList))
		rownames(m)<-CATS
		colnames(m)<-geneList

		for(cat in CATS) {
			w<-which(DB[,"GO_NAME"]==cat)
			
			genes<-intersect(DB[w,"HGNC"],geneList)
			m[cat,genes]<-1		
		}
		
		return(m)
	}

#' JaccardHeatMap
#' 
#' @description use the Jaccard metric to construct 2D heat map
#'
#' @param dir character string containing path name of output directory
#' @param mat return value of Jaccard()
#'
#' @examples
#' #load("data/x_jmat.RData")
#' dir<-tempdir()
#' jHeatMap<-JaccardHeatMap(dir,x_jmat)
#' 
#' @return returns a Jaccard matrix of cat1 vs cat2 FDR, and also saves hyperlinked SVG heatmap files to a results directory
#' 
#' @export
JaccardHeatMap<-
	function(dir,mat) {
		catrs<-unique(mat[,"catr"])
		catcs<-unique(mat[,"catc"])
		back<-max(as.numeric(mat[,"p"]))
		jmat<-matrix(back,nrow=length(catrs),ncol=length(catcs))
		rownames(jmat)<-catrs
		colnames(jmat)<-catcs
				
		for(r in 1:nrow(mat))
			jmat[mat[r,"catr"],mat[r,"catc"]]<-as.numeric(mat[r,"p"])
		
    svgWidth<-(8.0 + 0.5*ncol(mat)) * 0.526
    svgHeight<-(8.0 + 0.5*nrow(mat)) * 0.526

    filename<-sprintf("%s/%s",dir,"HTGM2D.svg")
	  svg(filename,width=svgWidth,height=svgHeight)
    # trick - use row for 'key' to get more space for long category names, but suppress key
    hm<-heatmap.2(jmat,col = heat.colors(n=100,rev=FALSE),trace="none",lmat=rbind( c(0, 3),
        c(2,1), c(0,4) ),lhei=c(1,4,15),lwid=c(1,50),key=FALSE,margins = c(1, 30))
    dev.off()
    
    ff<-hyperlinks(filename,rownames(jmat[hm$rowInd,]),colnames(jmat[,hm$colInd]))
    
	return(jmat)
	}

#' Jaccard
#' 
#' @description create the heat map data that is needed as input to JaccardHeatMap()
#'
#' @param dir character string full pathname to the directory acting as result repository
#' @param m1 return value of catGenes
#' @param m2 return value of catGenes
#' @param thresh1 integer acceptance threshold for the number of genes in a cat
#' @param thresh2 integer acceptance threshold for the number of common genes in 2 cats
#' @param B integer a total bootstrap iteration
#'
#' @examples
#' #load("data/x_m1.RData")
#' #load("data/x_m2.RData")
#' mat<-Jaccard(dir=tempdir(),x_m1,x_m2)
#' 
#' @return returns a numerical matrix containing number of genes and associated p value in the intersection of 2 categories
#' 
#' @export
Jaccard<-
	function(dir,m1,m2,thresh1=2,thresh2=3,B=100) {
		l<-list()
		
		subdir<-sprintf("%s/%s",dir,"hyperGenes")
    	if(!dir.exists(subdir))
      		dir.create(subdir)
		
		cn<-c("catr","catc","lint","nr","nc","p")
		mat<-matrix(ncol=length(cn))
		rownames(mat)<-NA
		colnames(mat)<-cn
		total<-length(rownames(m1))*length(rownames(m2))
		count<-0
		for(cat1 in rownames(m1)) {
			if(sum(m1[cat1,])<=thresh1)
				next
			n<-0
			for(cat2 in rownames(m2)) {
				count<-count+1
				if(sum(m2[cat2,])<=thresh1)
					next			
				w<-which(m1[cat1,]&m2[cat2,]==TRUE)
				if(length(w)<=thresh2)
					next
					
					# names(w) are the desired common gene names
					label<-sprintf("%s__%s.txt",substr(cat1,1,10),substr(cat2,1,10))
					writeLines(sort(names(w)),sprintf("%s/%s",subdir,label))				
									
				x<-jaccard.test(m1[cat1,],m2[cat2,],method="bootstrap",fix='0',B=B,verbose=FALSE)
				if(x$pvalue<=.1) {
				#	n<-n+1
				#	if(n==10) {
				#		n<-0
				#		print(c(count,total,count/total))
				#	}
					
					v<-c(cat1,cat2,length(w),sum(m1[cat1,]),sum(m2[cat2,]),x$pvalue)
					rn<-rownames(mat)
					mat<-rbind(mat,v)
					rownames(mat)<-c(rn,sprintf("%s_____%s",cat1,cat2))
				}
			}
		}
		return(mat[-1,]) # first row of mat is just a bunch of NA's
	}
