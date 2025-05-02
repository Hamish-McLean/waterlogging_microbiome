
list_biome.load <-
function(
	countData,
	colData,
	taxData,
	phylipData=NULL,
	tax_order=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),
	tax_conf=0.65,
	RHB=NULL
){
	if(missing(countData)|missing(colData)|missing(taxData)) {
		if(missing(countData)) missing_data("countData")
		if(missing(colData)) missing_data("colData")
		if(missing(taxData)) missing_data("taxData")
		stop("Missing values")
	}

	# load otu count table
	countData <- fread(countData,header=T,sep="\t",row.names=1, comment.char = "")

	# load sample metadata
	colData   <- fread(colData,header=T,sep="\t",row.names=1)

	# load taxonomy data
	taxData   <- fread(taxData,header=F,sep=",",row.names=1)

	# reorder columns
	taxData<-taxData[,tax_order]

	# add best "rank" at 0.65 confidence and tidy-up the table
	taxData<-phyloTaxaTidy(taxData,tax_conf)

	# save data into a list
	ubiom <- list(
		countData=countData,
		colData=colData,
		taxData=taxData,
		RHB=RHB
	)

	# get unifrac dist
	if(!missing(phylipData)) ubiom$phylipData <- fread.phylip(phylipData)

	return(ubiom)

	missing_data<-function(x) {
		print(paste(x, "must be specified"))
	}
}

as.number <-
function(f,convert=as.numeric) {
	convert(levels(f))[f]
}

batchEffectDes <-
function(obj,batch) {

	rld <- log2(counts(obj,normalize=T)+1)
	mypca <- prcomp(t(rld))

	myformula <- as.formula(paste("mypca$x~",batch,sep=""))

	pc.res <- resid(aov(myformula,obj@colData))

	mu <- colMeans(t(rld))

	Xhat <- pc.res %*% t(mypca$rotation)
	Xhat <- t(scale(Xhat, center = -mu, scale = FALSE))
	Xhat <- (2^Xhat)-1

	Xhat[Xhat<0] <- 0
	Xhat <- round(Xhat,6)

	return(Xhat)
}
calcCorrelog<- function(
	pca,
	obj,
	pc,
	na.add,
	condition=c("Y","N"),
	condColumn=1,
	returnInp,
	returnCD,
	useMeans=F
){
	cond<-condition[1]

	pc.x <- scores(pca)[rownames(scores(pca))%in%rownames(sample_data(obj)[sample_data(obj)[[condColumn]]==cond]),]
	col.x <- sample_data(obj)[sample_data(obj)[[condColumn]]==cond,]

	pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
	if(useMeans) {
		pc.reshape <- dcast(pc.dt,Distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
	} else {
		pc.reshape <-pc.dt
	}
	names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
	pc.reshape<-pc.reshape[order(pc.reshape$Distance),]
	dvec <- pc.reshape$Distance

	if (!missing(na.add)) {
		inp <- pc.reshape[[pc]]
		dvec<- unlist(sapply(1:length(inp),function(i) if(i%in%na.add){return(c(mean(c(dvec[i],dvec[(i-1)])),dvec[i]))}else{return(dvec[i])}))
		inp1 <- sapply(1:length(inp),function(i) if(i%in%na.add){return(c(NA,inp[i]))}else{return(inp[i])})
	}else {
		inp1<-pc.reshape[[pc]]
	}
	
	cond<-condition[2]

	pc.x <- scores(pca)[rownames(scores(pca))%in%rownames(sample_data(obj)[sample_data(obj)[[condColumn]]==cond]),]
	col.x <-  sample_data(obj)[sample_data(obj)[[condColumn]]==cond,]
	pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
	if(useMeans) {
		pc.reshape <- dcast(pc.dt,Distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
	} else {
		pc.reshape <-pc.dt
	}	

	names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])	

	if (!missing(na.add)) {
		inp <- pc.reshape[[pc]]
		inp2 <- sapply(1:length(inp),function(i) if(i%in%na.add){return(c(NA,inp[i]))}else{return(inp[i])})
	}else {
		inp2<-pc.reshape[[pc]]
	}
	
	if(returnInp) {
		return(cbind(unlist(inp1),unlist(inp2),dvec))
	}

	ct1 <- correr1(unlist(inp1),returnCD)
	ca1 <- correr1(unlist(inp2),returnCD)

	d<-as.data.frame(cbind(ct1,ca1,dvec[1:(length(dvec)-2)]))
	d
}# collapses samples with unequal libraries to mean and recalculates size factors 

collapseReplicates2 <- 
function (object, groupby, run, renameCols = TRUE,simple=F,method=mean)
{
    if (!is.factor(groupby)) groupby <- factor(groupby)
    groupby <- droplevels(groupby)
    stopifnot(length(groupby) == ncol(object))
    sp <- split(seq(along = groupby), groupby)
	if(simple) {sizeFactors(object) <- 1}
    sizefactors <- sapply(sp, function(i) prod(sizeFactors(object)[i,drop=F]))

    countdata <- sapply(sp, function(i) 
		rowMeans(counts(object,normalize=T)[,i, drop = FALSE]*prod(sizeFactors(object)[i,drop=F]))
	)

    mode(countdata) <- "integer"
    colsToKeep <- sapply(sp, `[`, 1)
    collapsed <- object[, colsToKeep]
    dimnames(countdata) <- dimnames(collapsed)
    assay(collapsed) <- countdata
    if (!missing(run)) {
        stopifnot(length(groupby) == length(run))
        colData(collapsed)$runsCollapsed <- sapply(sp, function(i) paste(run[i],
            collapse = ","))
    }
    if (renameCols) {
        colnames(collapsed) <- levels(groupby)
    }
    sizeFactors(collapsed) <- sizefactors

    #stopifnot(sum(as.numeric(assay(object))) == sum(as.numeric(assay(collapsed))))
    collapsed
}
combineTaxa2 <-
function(
	taxData, 
	rank="species", 
	# cut-off for matching ranks
	confidence=0.95,
	returnFull=F
) {
	
	# data tables and dplyr are going to be used 
	require(plyr)
	require(data.table)
	
	# new filter 
	levelled <- taxaConfVec(taxData,confidence,which(colnames(taxData)==rank))
	taxData <- data.table(taxData[grep("[^\\(].[^\\)]$",levelled),],keep.rownames="OTU")

	# group any remaining OTUs by the rank value, will return the OTUs as a character list
	taxData <- ddply(taxData,rank,summarize,OTUS=list(as.character(OTU)))
	
	# return taxa associated with more than one OTU
	if(!returnFull) {
		return(taxData[lapply(taxData[,2],function(x) length(unlist(x)))>1,])
	}

	taxData
}

combineTaxa <-
function(
	# path to the taxonomy file
	path, 
	# regex to specify filter for species (or other rank) column ( (s),(g),(f) and etc.)
	rank="species", 
	# cut-off for matching ranks
	confidence=0.95, 
	# column order MUST match k,p,c,o,f,g,s,k_conf,p_conf,c,conf,o_conf,f_conf,g_conf,s_conf
	# use column_order to reorder the columns to the required use -99 to keep the same order(or some other "large" number)  
	column_order=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),
	# whether the taxonomy file contains a header
	header=F, 
	# taxonomy file seperator
	sep=",",
	 # column with row names, set to NULL if no row names 
	row.names=1,
	# pass any further arguments to read.table
	...
	
) {
	
	# data tables and dplyr are going to be used 
	require(plyr)
	require(data.table)
	
	# read taxonomy file into a data frame
	taxData <- read.table(path,header=header,sep=sep,row.names=row.names,...)

	# reorder columns
	taxData<-taxData[,column_order]

	# add best "rank" at confidence and tidy-up the table
	taxData<-phyloTaxaTidy(taxData,confidence)

	# create a regex string out of the rank value
	rank_reg <- paste0("\\(",tolower(substr(rank,1,1)),"\\)")

	# filter the data for ranks at at least the given confidence
	taxData <- data.table(taxData[grep(rank_reg,taxData$rank),],keep.rownames="OTU")

	# group any remaining OTUs by the rank value, will return the OTUs as a character list
	taxData <- ddply(taxData,~rank,summarize,OTUS=list(as.character(OTU)))
	
	# return taxa associated with more than one OTU
	taxData[lapply(taxData[,2],function(x) length(unlist(x)))>1,]

}

combCounts <- 
function(
	# results from combineTaxa function
	combinedTaxa,
	# associated count data table
	countData
) {
	
	start<-nrow(countData)+1
	countData<-rbind(countData,t(sapply(combinedTaxa[,2],function(x) colSums(rbind(countData[rownames(countData)%in%unlist(x),],0)))))
	end <- nrow(countData)
	rownames(countData)[start:end] <- lapply(combinedTaxa[,2],function(x) paste(x[[1]][1],length(unlist(x)),sep="_"))
	countData <- countData[!row.names(countData)%in%unlist(combinedTaxa[,2]),]
  	countData[complete.cases(countData),]
}


combTaxa <- 
function (
	# results from combineTaxa function
	combinedTaxa,
	# the taxonomy table which was used for combineTaxa
	taxData
) { 
	keep <- taxData[unlist(lapply(combinedTaxa[,2],"[[",1)),]
	rownames(keep)<-lapply(combinedTaxa[,2],function(x) paste(x[[1]][1],length(unlist(x)),sep="_"))
	taxData<- taxData[!row.names(taxData)%in%unlist(combinedTaxa[,2]),]
	rbind(taxData,keep)
}
combine_biom <- function(locX,locY) {
	biom1 <- read.table(locX,header=T,sep="\t", comment.char="")	
	biom2 <- read.table(locY,header=T,sep="\t", comment.char="")
	biom <- merge(biom1,biom2,by.x="X.OTU.ID",by.y="X.OTU.ID",all=T)
	biom[is.na(biom)] <- 0
	return(biom)	
}########################################################################
#
# This is a modified (for speed) version of the CRAN package cooccur
# only hyper is currently available for prob calculation
#
########################################################################

coprob2 <-
function(max_inc,j,min_inc,nsite){
    require(gmp)
    as.matrix(round(chooseZ(max_inc,j) * chooseZ(nsite - max_inc, min_inc - j),0) / round(chooseZ(nsite,min_inc),0))
}

cooccur2 <- 
function (mat, type = "spp_site", thresh = TRUE, spp_names = FALSE,
    true_rand_classifier = 0.1, prob = "hyper", site_mask = NULL,
    only_effects = FALSE, eff_standard = TRUE, eff_matrix = FALSE)
{
    require(cooccur)
    if (type == "spp_site") {
        spp_site_mat <- mat
    }
    if (type == "site_spp") {
        spp_site_mat <- t(mat)
    }
    if (spp_names == TRUE) {
        spp_key <- data.frame(num = 1:nrow(spp_site_mat), spp = row.names(spp_site_mat))
    }
    if (!is.null(site_mask)) {
        if (nrow(site_mask) == nrow(spp_site_mat) & ncol(site_mask) ==
            ncol(spp_site_mat)) {
            N_matrix <- create.N.matrix(site_mask)
        }
        else {
            stop("Incorrect dimensions for site_mask, aborting.")
        }
    }
    else {
        site_mask <- matrix(data = 1, nrow = nrow(spp_site_mat),ncol = ncol(spp_site_mat))
        N_matrix <- matrix(data = ncol(spp_site_mat), nrow = nrow(spp_site_mat),ncol = nrow(spp_site_mat))
    }
    spp_site_mat[spp_site_mat > 0] <- 1
    tsites <- ncol(spp_site_mat)
    nspp <- nrow(spp_site_mat)
    spp_pairs <- choose(nspp, 2)
    incidence <- prob_occur <- obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow = spp_pairs,ncol = 3)
    incidence <- prob_occur <- matrix(nrow = nrow(N_matrix), ncol = ncol(N_matrix))

### hacked bit - need updating to add back only_effects=T (easy) and combinations probablity (not so easy, my current method requires huge amounts of memory) 

	mat_mask <- as.matrix(spp_site_mat*site_mask)
	incidence <- t(apply(mat_mask,1,function(m) site_mask%*%m))
	pairs <- t(apply(mat_mask,1,function(m) mat_mask%*%m))
	diag(incidence) <- NA
	prob_occur <- incidence/N_matrix

	obs_cooccur <- pairs
	prob_cooccur <- prob_occur*t(prob_occur)
	exp_cooccur <- prob_cooccur*N_matrix

	if (thresh) {
		n_pairs <- sum(prob_cooccur>=0,na.rm=T)/2
		t_table <- exp_cooccur>=1
		prob_cooccur <- prob_cooccur*t_table
		obs_cooccur <- obs_cooccur*t_table
		exp_cooccur <- exp_cooccur*t_table
		n_omitted <- n_pairs - sum(t_table,na.rm=T)
	}

	sp1_inc=incidence
	sp2_inc=t(incidence)
	max_inc <- pmax(sp1_inc, sp2_inc)
	min_inc <- pmin(sp1_inc, sp2_inc)
	nsite <- N_matrix
	psite <- nsite + 1

	arr <- array(c(min_inc,nsite,max_inc,(sp1_inc + sp2_inc)),c(nrow(nsite),ncol(nsite),4))

	effect_func <- function(a=arr){
		X <- apply(a,1:2, function(x) {
			x[is.na(x)]<-0
			i <- x[1]-((x[2]-x[4])*(-1)+abs((x[2]-x[4])*(-1)))/2
			ii <- (i+abs(i))/2
			#ii[is.na(ii)]<-0
			y<-rep(1,ii)
			return(y)
		})
	}
	
	hyper_func <- function(a=arr){
		X <- apply(a, 1:2 , function(x) {
			x[is.na(x)]<-0
			y<-phyper(0:x[1],x[1],x[2]-x[1],x[3])
			y<-c(y[1],(y[-1]-y[-length(y)]))
			return(y)
		})
		return(X)
	}

	hyper_func_quick <- function(a=arr){
		X <- apply(a, 1:2 , function(x) {
			x[is.na(x)]<-0
			y<-phyper(0:x[1],x[1],x[2]-x[1],x[3])
			y<-c(y[1],(y[-1]-y[-length(y)]))
			return(y)
		})
		return(X)
	}

	hyper_abundance_func <- function(a=arr){
		X <- apply(a, 1:2 , function(x) {
			x[is.na(x)]<-0
			y<-phyper(0:x[1],x[1],x[3],x[3])
			y<-c(y[1],(y[-1]-y[-length(y)]))
			return(y)
		})
		return(X)
	}

	comb_func_himem  <- function(a=arr){
		a[is.na(a)] <- 0
		y <- sapply(0:max(a[,,2]),function(i) coprob2(a[,,3],i,a[,,1],a[,,2]))
		y <- do.call(cbind, y)
		start <- (a[,,4]-a[,,2]+abs(a[,,4]-a[,,2]))/2
		end <- a[,,1]
		#i <- a[,,1]-((a[,,2]-a[,,4])*(-1)+abs((a[,,2]-a[,,4])*(-1)))/2
		#len <- (i+abs(i))/2
		x <- matrix(nrow=nrow(y),ncol=ncol(y))
		a1 <- ifelse(col(x)>=(as.vector(start)+1),1,NA)
		a2 <- ifelse(col(x)<=(as.vector(end)+1),1,NA)
		return(array(as.numeric(y)*a1*a2,c(10,10,11)))
	}

	comb_func  <- function(a=arr){
		X <- apply(arr,1:2, function(x) {
			x[is.na(x)]<-0
			y<-as.numeric(coprob2(x[3],0:x[2],x[1],x[2]))
			start <- (x[4]-x[2]+abs(x[4]-x[2]))/2			
			#i <- x[1]-((x[2]-x[4])*(-1)+abs((x[2]-x[4])*(-1)))/2
			#len <- (i+abs(i))/2
			#ii[is.na(ii)]<-0
			return(y[start:(x[1]+1)])
		})
		return(X)
	}
#return(comb_func())
	if (only_effects) {
		prob_share_site <- effect_func()
	} else {
		if (prob == "hyper") {
			prob_share_site <- hyper_func()
		} else if (prob == "comb") {		
			prob_share_site <- comb_func()
		} else if (prob == "comb_hi") {		
			prob_share_site <- comb_func_himem()	
		} else {
			print("Unsupported probability model specified\n. Returning effect sizes only")
			prob_share_site<-effect_func()
			only_effect=T
		}
    	}

	prob_share_site<- prob_share_site[which(lower.tri(prob_share_site))]
	obs_cooccur<- obs_cooccur[which(lower.tri(obs_cooccur))]
	prob_cooccur<- prob_cooccur[which(lower.tri(prob_cooccur))]
	exp_cooccur<- exp_cooccur[which(lower.tri(exp_cooccur))]
	t_table <- t_table[which(lower.tri(t_table))]

	sp <- matrix(rep(seq(1,nspp),nspp),nrow=nspp,ncol=nspp)
	sp1 <- t(sp)[which(lower.tri(sp[-nrow(sp),-ncol(sp)],diag=T))]
	sp2 <- sp[which(lower.tri(sp,diag=F))]

	sp <- matrix(rep(rownames(sp1_inc),ncol(sp1_inc)),nrow=ncol(sp1_inc),ncol=ncol(sp1_inc))
	sp1_name <- t(sp)[which(lower.tri(sp[-nrow(sp),-ncol(sp)],diag=T))]
	sp2_name <- sp[which(lower.tri(sp,diag=F))]

	sp1_inc<- sp1_inc[which(lower.tri(sp1_inc))]
	sp2_inc<- sp2_inc[which(lower.tri(sp2_inc))]

	p_lt <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[1:(obs_cooccur[i]+1)]))
	p_gt <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[(obs_cooccur[i]+1):length(unlist(prob_share_site[i]))]))
	p_exactly_obs <- sapply(seq(1,length(prob_share_site)),function(i) sum(unlist(prob_share_site[i])[(obs_cooccur[i]+1)]))

	p_lt <- round(p_lt, 5)
	p_gt <- round(p_gt, 5)
	p_exactly_obs <- round(p_exactly_obs, 5)

	prob_cooccur <- round(prob_cooccur, 3)
	exp_cooccur <- round(exp_cooccur, 1)

	output<-data.frame(
		sp1=sp1,
		sp2=sp2,
		sp1_inc=sp2_inc,
		sp2_inc=sp1_inc,
		obs_cooccur=obs_cooccur,
		prob_cooccur=prob_cooccur,
		exp_cooccur=exp_cooccur,
		p_lt=p_lt,
		p_gt=p_gt,
		sp1_name=sp1_name,
		sp2_name=sp2_name
	)
	if (thresh) {
		# could do this earlier and save a couple of minutes execution time
		output <- output[t_table,]
	}
####

    true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >=
        0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <=
        (tsites * true_rand_classifier)), ]))
    output_list <- list(call = match.call(), results = output,
        positive = nrow(output[output$p_gt < 0.05, ]), negative = nrow(output[output$p_lt <
            0.05, ]), co_occurrences = (nrow(output[output$p_gt <
            0.05 | output$p_lt < 0.05, ])), pairs = nrow(output),
        random = true_rand, unclassifiable = nrow(output) - (true_rand +
            nrow(output[output$p_gt < 0.05, ]) + nrow(output[output$p_lt <
            0.05, ])), sites = N_matrix, species = nspp, percent_sig = (((nrow(output[output$p_gt <
            0.05 | output$p_lt < 0.05, ])))/(nrow(output))) *
            100, true_rand_classifier = true_rand_classifier)
    if (spp_names == TRUE) {
        output_list$spp_key <- spp_key
        output_list$spp.names = row.names(spp_site_mat)
    }
    else {
        output_list$spp.names = c(1:nrow(spp_site_mat))
    }
    if (thresh == TRUE) {
        output_list$omitted <- n_omitted
        output_list$pot_pairs <- n_pairs
    }
    class(output_list) <- "cooccur"
    if (only_effects == F) {
        output_list
    }
    else {
        effect.sizes(mod = output_list, standardized = eff_standard,
            matrix = eff_matrix)
    }
}
correr1 <- function(
	x,
	returnData=F
) {
	y <- x
	count<-1
	mycorr <- NULL
	while (length(x) >2) {
		if(returnData) {
			mycorr[[count]] <- cbind(y[1:(length(x)-1)],x[2:length(x)])
		} else {
			mycorr[count] <- cor(y[1:(length(x)-1)],x[2:length(x)],use="pairwise.complete.obs")
		}
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}
correr1_dat <- function(x) {
	y <- x
	mycorr <- NULL
	count<-1
	while (length(x) >2) {
		mycorr[[count]] <- cbind(y[1:(length(x)-1)],x[2:length(x)])
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}
correr2 <- function(
	X,
	breaks
){
	y <- seq(0,max(X[,2])+1,breaks)
	vec <- X[,2] 
	test <- nearest.vec(y,vec)	
	test[9] <- NA
	x <- X[X[,2]==test,1]
	mycorr <- numeric(0)
	count<-1
	while (length(x) >2) {
		mycorr[count] <- cor(y[1:(length(x)-1)],x[2:length(x)],use="pairwise.complete.obs")
		count <- count +1
		x<-x[-1]
	}
	return(mycorr)
}

countTaxa <- function(
	obj,
	taxon="phylum"
){
	# sums number of unique entries at given level 
	suppressPackageStartupMessages(require(dplyr))
	#suppressPackageStartupMessages(require(data.table))
	data.frame(Taxa=as.data.frame(as.matrix(obj))[[taxon]])

	dtx <-data.frame(Taxa=as.data.frame(as.matrix(obj))[[taxon]])
	dtx %>% group_by(Taxa) %>% summarise(count=length(Taxa)) #dply method
	
#	return(setDT(dtx)[, .N, keyby=Taxa]) # data table method - doesn't work in function??
	
}

countTaxa2 <- 
function(
        obj,
        taxon="phylum"
){
        # sums number of unique entries at given level
        suppressPackageStartupMessages(require(data.table))
        data.table(Taxa=as.data.frame(as.matrix(obj))[[taxon]])[, .N, keyby=Taxa]
}# functions for use with data table

#sets NA to 0 by reference
unsetNA = function(DT) {
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

#sets 0 to NA by ref
setNA = function(DT) {
  for (j in seq_len(ncol(DT)))
    set(DT,which(DT[[j]]==0),j,NA)
}
  

des_to_pca <- function(obj) {
	X<-prcomp(t(assay(varianceStabilizingTransformation(obj))))
	X$percentVar<-X$sdev^2/sum(X$sdev^2)
	X
}

des_to_phylo <- function(
	obj
){
	suppressPackageStartupMessages(require(phyloseq))
	suppressPackageStartupMessages(require(DESeq2))
	X<-
	phyloseq(
		otu_table(assay(obj),taxa_are_rows=T),
	 	sample_data(as.data.frame(colData(obj)))
	 )
}

dt_to_df <- function(
	DT,
	row_names=1
){
	DF <- as.data.frame(DT)
	row.names(DF) <- DF[,row_names]
	DF <- DF[,-row_names]
	DF
} 
fltTaxon <- function(
	obj,
	taxon="phylum",
	out="phylo"
){
# same as phyloseq tax_glom (drops NA columns returned by tax_glom), but works on S3 biom data (i.e. ubiom)
# perhaps tax_glom is good with big datasets as fltTaxon is miles faster - aggregate will get pretty slow for large datasets
	if(class(obj)[[1]]=="phyloseq") {
		obj <- phylo_to_ubiom(obj)
	}
	n <- which(colnames(obj[[2]])==taxon)
	x <- aggregate(obj[[1]],by=obj[[2]][,1:n],sum)
	ls.biom <- list(x[,(n+1):ncol(x)],x[,1:n],obj[[3]])
	names(ls.biom) <- c("countData","taxonomy","colData")
	if(out=="phylo") {
		return(ubiom_to_phylo(ls.biom))
	}
	return(ls.biom)	
}

fread.phylip <- 
function(path,...){
	dt <- fread(path,...)
	m <- as.matrix(dt[,-1])
	rownames(m) <- colnames(m) <- dt[[1]]
	return(m)
}# function sometimes useful for replacing calcfactors
geoMeans <- function(d,dummy) {
	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}	
	gm = apply(counts(d), 1, gm_mean)
	sizeFactors(estimateSizeFactors(d, geoMeans =gm))
}

geoSet <- function(d,dummy) {
	tryCatch ({
		sizeFactors(estimateSizeFactors(d))
	}, error = function(e) {
		cat("Sizefactors: Using geoMeans","\n")
		return(geoMeans(d))
	})
	
}#' @title Remove a layer from a compiled ggplot2 object.
#' @export
#' @description Removes specified layers from a ggplot object.
#' @param p ggplot2 plot object
#' @param geom character string of the name of the layer to remove
#' @param idx numeric of which index of geom to remove
#' @examples
#' p=ggplot(iris,aes(x =Sepal.Length,y=Sepal.Width))
#' p=p+geom_point(aes(colour=Species))+geom_line()
#' p
#' pnew=p%>%remove_geom('point',1)
#' pnew
remove_geom=function(p, geom,idx,labels=NULL){
  layers=p$layers
  layer.type=lapply(p$layers,function(x) class(x$geom))
  a.rm=which(grepl(paste0("(?i)", geom), layer.type))
  
  if(length(a.rm)==0) stop(paste0("There are less no ",geom," layers available in the plot to remove"), call. = FALSE)
  
  if(idx>length(a.rm)) stop(paste0("There are less than ",idx," ",geom," layers available in the plot to remove"), call. = FALSE)
  
  if(length(a.rm)>=idx) a.rm=a.rm[idx]
  layers <- layers[-a.rm]
  p$layers <- layers
  lapply(labels,function(l) p$labels[[l]] <<- NULL)
  p
}

ggplot_legend <- 
function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

get_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

ggplot_legend <- 
function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

theme_blank <- function(base_size = 11, base_family = "")
{
	theme_bw(base_size = base_size, base_family = base_family) %+replace% 
	theme(
		panel.border = element_blank(), 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank()
	)
}

theme_classic_thin <- function(base_size = 11, base_family = "") 
{
	theme_blank(base_size = base_size, base_family = base_family) %+replace% 
	theme(	
		axis.line.x = element_line(size=0.3,colour = "black"),
		axis.line.y = element_line(size=0.3,colour = "black"),
		axis.text = element_text(colour = "black")
	)
}


theme_facet_blank <- function(base_size = 11, base_family = "",angle=-90,t=2,r=0,b=0,l=0,hjust=0,vjust=0)
{
	theme_classic_thin(base_size = base_size, base_family = base_family) %+replace% 
	theme(
		panel.border = element_rect(colour = "black", fill=NA, size=0.5),
		axis.text.x = element_text(angle = angle, margin = margin(t=t,r=r,b=b,l=l),hjust=hjust,vjust=vjust)
	)
}



#  theme(
#	legend.title = element_blank(),
#	legend.position="bottom", 
#   	legend.direction="horizontal",
#   	legend.key = element_rect(colour = NA),
#	text=element_text(size=18)
#  )
import_ubiom <- function (
	locX,
	locY,
	locZ
){
	options(stringsAsFactors = FALSE)
	countData <- read.table(locX,header=T,sep="\t", comment.char="")
	rownames(countData ) <- countData [,1]
	countData <- countData [,-1]
	taxonomy <- read.csv(locY,header=F)
	taxonomy <- taxonomy [,c(1,2,4,6,8,10,12,14)]
	rownames(taxonomy) <- taxonomy[,1]
	taxonomy <- taxonomy[,-1]
	colnames(taxonomy) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
	colData <- read.table(locZ,sep="\t",header=T)
	rownames(colData) <- colData [,1]
	colData <- colData[,-1,drop=FALSE]
	countData <- countData[,rownames(colData)]
	ls.biom <- list(countData,colData, taxonomy)
	names(ls.biom) <- c("countData","colData","taxonomy")
	return(ls.biom)
}

loadData <-
function(
	countData,
	colData,
	taxData,
	phylipData,
	tax_order=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),
	tax_conf=0.65,
	RHB=NULL
){
	if(missing(countData)|missing(colData)|missing(taxData)) {
		if(missing(countData)) missing_data("countData")
		if(missing(colData)) missing_data("colData")
		if(missing(taxData)) missing_data("taxData")
		stop("Missing values")
	}

	# load otu count table
	countData <- read.table(countData,header=T,sep="\t",row.names=1, comment.char = "")

	# load sample metadata
	colData <- read.table(colData,header=T,sep="\t",row.names=1)

	# load taxonomy data
	taxData <- read.table(taxData,header=F,sep=",",row.names=1)

	# reorder columns
	taxData<-taxData[,tax_order]

	# add best "rank" at 0.65 confidence and tidy-up the table
	taxData<-phyloTaxaTidy(taxData,tax_conf)

	# save data into a list
	ubiom <- list(
		countData=countData,
		colData=colData,
		taxData=taxData,
		RHB=RHB
	)

	# get unifrac dist
	if(!missing(phylipData)) ubiom$phylipData <- fread.phylip(phylipData)

	return(ubiom)

	missing_data<-function(x) {
		print(paste(x, "must be specified"))
	}
}# adapted version of phyloTaxaTidy to mash names by conf
phyloTidy <- function(taxData,conf,setCount=T) {
  X <- taxData[,1:7]
  Y <- taxData[,9:15]
  Y[,1] <-1.00
  X[apply(Y,2,function(y)y<conf)] <- "unknown"
  n <- c("kingdom","phylum","class","order","family","genus","species")
  X <-apply(X,2,function(x) {
    x <- sub("\\(.*","",x)
    x <- sub("_SH[0-9]+\\..*","",x)
    x <- gsub("_"," ",x)
    #x <- sub("_"," ",x)
    x
  })
  
  td <- as.data.table(X,keep.rownames = "OTU")
  cols <- names(td)[-1]
  num <- td[,.N,by=cols]
  num2 <- invisible(apply(num,1,function(x) as.data.table(t(matrix(x[-8],nrow=7,ncol=as.numeric(x[8]))))))
  invisible(lapply(num2,setnames,cols))
  invisible(lapply(num2,function(x)x[,counts:=1:nrow(x)]))
  td2 <- do.call("rbind",num2)
  td[,(cols):=lapply(.SD,as.character),.SDcols=cols]
  setorderv(td, cols[-1])
  setorderv(td2,cols[-1])
  
  td$counts <- td2$counts
  td[,(cols):=lapply(.SD,as.factor),.SDcols=cols]
  taxData <- as.data.frame(td)
  rownames(taxData) <- taxData$OTU
  taxData <- taxData[,-1]
  if(setCount) {
    taxData$species <- as.factor(paste(taxData$species,taxData$counts,sep=" "))
    taxData <- taxData[,-8]
  }
  taxData
}

funq <- function(lfc,p){lfc[c(which(p>0.1),which(is.na(p)))] <- 0;lfc}

qf <- function(x)sub("Fungi;unknown.*","REMOVE",x)

# get the full taxonomy for each taxon_id
# get_full_taxa <- function(taxon_id,obj) {
#   t_var   <- c(obj$taxon_indexes()[[taxon_id]],obj$supertaxa()[[taxon_id]])
#   t_names <- sapply(t_var,function(x) {obj$taxa[[x]]$get_name()})
#   paste(rev(t_names),collapse = ";")
# }

get_full_taxa <- function(taxon_id,t1=t1,t2=t2) {
  t_names <- t2[c(taxon_id,unlist(t1[taxon_id]))]
  paste(rev(t_names),collapse = ";")
}

MCDESeq <- function (OBJ,formula,Rank,colData) {
  
  o <- OBJ$clone(deep=T)
  formula <- formula(formula)
  # set rank for collapsing
  #  Rank = "genus"
  # collapse OTU abundances at each taxonomic rank at and above Rank 
  o$data$tax_abund <- o %>% 
    taxa::filter_taxa(taxon_ranks == Rank,supertaxa=T) %>%
    calc_taxon_abund(data = "otu_table",cols = colData$sample_id)
  
  # metacoder differential analysis - this is bobbins  
  #FUN_ENDO$data$diff_table <- compare_groups(FUN_ENDO, dataset = "otu_table",
  #                                           cols = colData$sample_id, # What columns of sample data to use
  #                                           groups = colData$status) # What category each sample is assigned to
  
  # get the otu table with taxon_ids  - rejig to get taxon_id as rownames
  
  
  otu_table <- as.data.table(o$data$otu_table) #metacoder:::get_taxmap_table(OBJ, "otu_table")
  # merge duplicate "species"
  numeric_cols <- which(sapply(otu_table, is.integer))
  otu_table <- as.data.frame(otu_table[, lapply(.SD, sum), by = taxon_id, .SDcols = numeric_cols])
  rownames(otu_table) <- otu_table[,1]
  otu_table <- otu_table[,-1]
  
  
  # get Rank and above abundances - rejig to get taxon_id as rownames
  tax_abund <- as.data.frame(o$data$tax_abund) #metacoder:::get_taxmap_table(OBJ, "tax_abund")
  rownames(tax_abund) <- tax_abund[,1]
  tax_abund <- tax_abund[,-1]
  
  # set character columns to factors
  colData <- as.data.frame(unclass(as.data.frame(colData)))
  
  # create new dds object from combined tables
  dds <- DESeqDataSetFromMatrix(rbind(otu_table,tax_abund),colData,~1)
  #cols <- colnames(colData(dds))[sapply(colData(dds),is.factor)]
  #dds$status <- as.factor(dds$status)
  design(dds) <- formula#~status
  dds <- DESeq(dds)
  return(dds)
  # res <- results(dds,alpha=alpha,contrast=contrast)
  # 
  # # make a data table from the resuults
  # res_merge <- as.data.table(as.data.frame(res),keep.rownames="taxon_id")
  # 
  # # add second log fold change column with fold changes for non sig taxa set to 0 (for colouring purposes on graphs)
  # res_merge[,log2FoldChange2:=funq(log2FoldChange,padj)]
  # 
  # # add full taxonomy to results (this is slow - good reason not to use S4 object model!)
  # 
  # # order the results
  # setkey(res_merge,taxon_id)
  # 
  # # these are required due to the restrictive interace of taxa/metacoder 
  # # reduces the time required to get the full taxonomy of each taxon id by about x10,000
  # t1 <- o$supertaxa()
  # t2 <- o$taxon_names()
  # res_merge[,taxonomy:=sapply(1:nrow(res_merge),get_full_taxa,t1,t2)]
  # 
  # # add results as diff_table (could be called anything - but for consistency with metacoder keep same names)
  # #
  # as_tibble(res_merge)
}

MCDESres <- function (OBJ,dds,Rank,colData,contrast,alpha=0.1) {
  
  o <- OBJ$clone(deep=T)
  # set rank for collapsing
  #  Rank = "genus"
  # collapse OTU abundances at each taxonomic rank at and above Rank 
  o$data$tax_abund <- o %>% 
    taxa::filter_taxa(taxon_ranks == Rank,supertaxa=T) %>%
    calc_taxon_abund(data = "otu_table",cols = colData$sample_id)
  
  # metacoder differential analysis - this is bobbins  
  #FUN_ENDO$data$diff_table <- compare_groups(FUN_ENDO, dataset = "otu_table",
  #                                           cols = colData$sample_id, # What columns of sample data to use
  #                                           groups = colData$status) # What category each sample is assigned to
  
  # get the otu table with taxon_ids  - rejig to get taxon_id as rownames
  otu_table <- as.data.table(o$data$otu_table) #metacoder:::get_taxmap_table(OBJ, "otu_table")
  # merge duplicate "species"
  numeric_cols <- which(sapply(otu_table, is.integer))
  otu_table <- as.data.frame(otu_table[, lapply(.SD, sum), by = taxon_id, .SDcols = numeric_cols])
  rownames(otu_table) <- otu_table[,1]
  otu_table <- otu_table[,-1]  
  # get Rank and above abundances - rejig to get taxon_id as rownames
  tax_abund <- as.data.frame(o$data$tax_abund) #metacoder:::get_taxmap_table(OBJ, "tax_abund")
  rownames(tax_abund) <- tax_abund[,1]
  tax_abund <- tax_abund[,c(-1)]
  
  # set character columns to factors
  colData <- as.data.frame(unclass(as.data.frame(colData)))
  
  # create new dds object from combined tables
  #dds <- DESeqDataSetFromMatrix(rbind(otu_table,tax_abund),colData,~1)
  #cols <- colnames(colData(dds))[sapply(colData(dds),is.factor)]
  #dds$status <- as.factor(dds$status)
  #design(dds) <- formula#~status
  #dds <- DESeq(dds)
  res <- results(dds,alpha=alpha,contrast=contrast)
  
  # make a data table from the resuults
  res_merge <- as.data.table(as.data.frame(res),keep.rownames="taxon_id")
  
  # add second log fold change column with fold changes for non sig taxa set to 0 (for colouring purposes on graphs)
  res_merge[,log2FoldChange2:=funq(log2FoldChange,padj)]
  
  # add full taxonomy to results (this is slow - good reason not to use S4 object model!)
  
  # order the results
  setkey(res_merge,taxon_id)
  
  # these are required due to the restrictive interace of taxa/metacoder 
  # reduces the time required to get the full taxonomy of each taxon id by about x10,000
  t1 <- o$supertaxa()
  t2 <- o$taxon_names()
  res_merge[,taxonomy:=sapply(1:nrow(res_merge),get_full_taxa,t1,t2)]
  
  # add results as diff_table (could be called anything - but for consistency with metacoder keep same names)
  #
  as_tibble(res_merge)
}


nearest.vec <- function(
	x, 
	vec,
	breaks=1
) {
    smallCandidate <- findInterval(x, vec, all.inside=TRUE)
    largeCandidate <- smallCandidate + 1
    nudge <- 2 * x > vec[smallCandidate] + vec[largeCandidate]
    empty <- abs(vec[smallCanditate]-break)
    return(smallCandidate + nudge)
}

normHTS <- function(
	obj,
	controlSamples)
{

	gm_mean = function(x, na.rm=TRUE){
		exp(sum(log(x[x > 0]), na.rm=na.rm) / sum(!is.na(x)))
	}

	uniqueRuns <- unique(sub("S[0-9]*","",row.names(sample_data(obj))))

	dds.all <- sapply(uniqueRuns,function(x) phylo_to_des(
		prune_samples(sub("S[0-9]*","",row.names(sample_data(obj)))==x,obj)))
	
	c.counts <- lapply(dds.all,function(x) counts(x,normalize=T)[,colnames(x)%in%controlSamples])

	cc1 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,1],error=function(e)NA))))
	cc2 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,2],error=function(e)NA))))
	cc3 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,3],error=function(e)NA))))
	cc4 <- as.data.frame(lapply(c.counts,function(x) return(tryCatch(x[,4],error=function(e)NA))))
	
	x2 <- t(data.frame(
		cc1=colSums(cc1)/colSums(cc1),
		cc2=colSums(cc2)/colSums(cc2),
		cc3=colSums(cc3)/colSums(cc3),
		cc4=colSums(cc4)/colSums(cc4)	
	))

	cc1[is.na(cc1)] <- apply(cc1,1,gm_mean)
	cc2[is.na(cc2)] <- apply(cc2,1,gm_mean)
	cc3[is.na(cc3)] <- apply(cc3,1,gm_mean)
	cc4[is.na(cc4)] <- apply(cc4,1,gm_mean)

	x1 <- rbind(
		estimateSizeFactorsForMatrix(cc1),
		estimateSizeFactorsForMatrix(cc2),
		estimateSizeFactorsForMatrix(cc3),
		estimateSizeFactorsForMatrix(cc4)
	) 
	
	x3 <- x1 * x2
	x4 <- apply(x3,2,function(x) median(x,na.rm=T))	
	x4[is.na(x4)] <- 1	
	
	xx <- unlist(lapply(dds.all,sizeFactors))

	xy <- xx*sapply(names(xx),function(x) x4[sub("\\..*","",x)])
	
	names(xy) <- sub(".*\\.","",names(xy))

	return(xy)

	#k <- which(is.na(x3), arr.ind=TRUE)
	#x3[k] <- apply(x3,1,gm_mean)[k[,1]]
	
	#estimateSizeFactorsForMatrix(x3)
}



norm_table <- 
function(phylo,calcFactors=function(o,...)estimateSizeFactorsForMatrix(o,...),...)
{

	sizeFactors <- calcFactors(otu_table(phylo),...)
	if(missing(sizeFactors)) {
		sizeFactors <- sample_data(phylo)$sizeFactors
	}
	
	 t(t(otu_table(phylo))/sizeFactors)
}

ordinate <- 
function (physeq, method = "DCA", distance = "bray", formula = NULL,
    ...)
{
    library(phyloseq)
    if (inherits(physeq, "formula")) {
        .Deprecated(msg = paste0("First argument, `physeq`, as formula is deprecated.\n",
            "There is now an explicit `formula` argument.\n",
            "Please revise method call accordingly."))
        formchar = as.character(physeq)
        if (length(formchar) < 3) {
            stop("Need both sides of formula in this deprecated syntax... Revisit ordinate() documentation / examples.")
        }
        physeq <- get(as.character(physeq)[2])
        newFormula = as.formula(paste0("~", formchar[length(formchar)]))
        return(ordinate(physeq, method = method, distance = distance,
            formula = newFormula, ...))
    }
    method_table <- c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS",
        "MDS", "PCoA")
    if (inherits(physeq, "character")) {
        if (physeq == "help") {
            cat("Available arguments to methods:\n")
            print(c(method_table))
            cat("Please be exact, partial-matching not supported.\n")
            cat("Can alternatively provide a custom distance.\n")
            cat("See:\n help(\"distance\") \n")
            return()
        }
        else if (physeq == "list") {
            return(c(method_table))
        }
        else {
            cat("physeq needs to be a phyloseq-class object, \n")
            cat("or a character string matching \"help\" or \"list\". \n")
        }
    }
    if (!inherits(physeq, "phyloseq") & !inherits(physeq, "otu_table")) {
        stop("Expected a phyloseq object or otu_table object.")
    }
    if (method == "DCA") {
        return(decorana(veganifyOTU(physeq), ...))
    }
    if (method %in% c("CCA", "RDA")) {
        return(cca.phyloseq(physeq, formula, method, ...))
    }
    if (method == "CAP") {
        return(capscale.phyloseq(physeq, formula, distance, ...))
    }
    if (method == "DPCoA") {
        return(DPCoA(physeq, ...))
    }
    if (inherits(distance, "dist")) {
        ps.dist <- distance
    }
    else if (class(distance) == "character") {
        vegdist_methods <- c("manhattan", "euclidean", "canberra",
            "bray", "kulczynski", "jaccard", "gower", "altGower",
            "morisita", "horn", "mountford", "raup", "binomial",
            "chao")
        if (method == "NMDS" & distance %in% vegdist_methods) {
            return(metaMDS(veganifyOTU(physeq), distance, ...))
        }
        ps.dist <- distance(physeq, distance, ...)
    }
    if (method %in% c("PCoA", "MDS")) {
        return(pcoa(ps.dist))
    }
    if (method == "NMDS") {
        return(metaMDS(ps.dist))
    }
}
phyloTaxaTidy <- function(obj,...) {
	if(ncol(obj)==7){
		colnames(obj) <- c(
			"kingdom","phylum","class", "order", "family","genus", "species"
		)
	} else {
		colnames(obj) <- c(
			"kingdom","phylum","class", "order", "family","genus", "species",
			"k_conf", "p_conf","c_conf","o_conf","f_conf","g_conf","s_conf"
		)
		
		obj[,8:14] <- as.numeric(unlist(obj[,8:14]))
		obj[obj[,8:14]==100,8:14] <- 0
		#obj[obj[,8:14]==-1,8:14] <- 1
		obj <- taxaConf(obj,...)
		obj <- obj[,c(1:7,15,8:14)]
	}
	#obj <- sub("_+"," ",obj)
	
	obj[,1:7] <- t(apply(obj[,1:7],1,taxonomyTidy))
	return(obj)
}
phylo_to_des <- function(
	obj,
	design=~1,
	fit=F,
	obj2=NA,	
	calcFactors=function(d,o)
	{
		sizeFactors(estimateSizeFactors(d))
	},
	...
){
	suppressPackageStartupMessages(require(phyloseq))
	suppressPackageStartupMessages(require(DESeq2))
	dds <-  phyloseq_to_deseq2(obj,design)
	sizeFactors(dds) <- calcFactors(dds,obj2)
    	if (fit) {
    	 	return(DESeq(dds,...))
    	} else {
    		return(dds)
    	}
} 

phylo_to_ubiom <- function(
	obj
){
	suppressPackageStartupMessages(require(phyloseq))
	list(
		countData=as.data.frame(obj@otu_table@.Data),
		taxonomy=as.data.frame(obj@tax_table@.Data),
		colData=as.data.frame(suppressWarnings(as.matrix(sample_data(obj)))) # suppresses a warning from the matrix call about loss of S4 state
	)
}

plotCorrelog <- function(
	pca,
	obj,
	pc="PC1",
	cutoff=15,
	xlim=NULL,
	ylim=NULL,
	na.add,
	returnData=F,
	returnInp=F,
	returnCD=F,
	data,
	legend=T,
	lpos="right",
	cols=c("#edf8b1","#2c7fb8"),
	lineWidth=1.5,
	useMeans=F,
	textSize=12,
	design=c("Tree","Grass") 
) {

	if(returnInp|returnCD){
		returnData=T
	}
	
	if(!missing(data)){
		d<-data
	} else {
		d <- calcCorrelog(pca,obj,pc,na.add,design,1,returnInp,returnCD,useMeans)
	}
	
	if(returnData) {
		return(d)
	}
	d<-d[order(d$V3),]
	d<- d[1:(length(d$V3[d$V3<=cutoff])),]
	names(d) <- c("Tree","Grass","Distance")
	d2 <- melt(d,id="Distance")
	colnames(d2) <- c("Distance","Sample","Correlation")

	g <- ggplot(d2)
	g <- g + coord_fixed(ratio = 10, xlim = xlim, ylim = ylim, expand = TRUE)
	
	# g <- g + theme_classic()

#	g <- g + theme_bw()
#	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#		legend.position=lpos,legend.text=element_text(size=10),legend.title=element_text(size=10), legend.background = element_rect(fill=alpha('white', 0)))
# 	g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))

	if(!legend) {
		g <- g + theme_classic(base_size=textSize) %+replace% theme(legend.position="none")
	} else {
		g <- g + theme_classic(base_size=textSize) + 
			theme(
			legend.position=lpos,
			legend.text=element_text(size=10),
			legend.title=element_text(size=10),
			legend.background = element_rect(fill=alpha('white', 0))
		)
	}


	g <- g + geom_line(na.rm=T,aes(x=Distance, y=Correlation, colour=Sample),size=lineWidth)
	g <- g + scale_colour_manual(values=cols)
	g
}

plotCummulativeReads <- 
function(countData,cutoffs=c(0.8,0.9,0.99,0.999),returnData="dtt",plot=TRUE,bysample=F)
{
	
	require(ggplot2)
	require(data.table)
	
	# calculate row sums from normalised read counts
    if(bysample){
		DT <- data.table(apply(countData,2,sort,decreasing=T))
   	} else {
		DT <- data.table(CD=rowSums(countData),OTU=rownames(countData))
		setorder(DT,-CD)
	}


#	DT <- DT[order("CD",decreasing=T)]
	print(paste("Returning data table as",returnData))
 	assign(returnData, DT, envir=globalenv())
	if(!plot)return("Not plotting")

	# calculate cumulative sum of ordered OTU counts 
    #DT <- cumsum(DT)
    if(!bysample){DT<-DT[-2]}
	DT <- cumsum(DT)
	#DT$CD <- cumsum(DT[,"CD"])

	suppressWarnings(if(cutoffs&!bysample) {
		# get cutoff poisiotns
		mylines <- data.table(RHB=sapply(cutoffs,function(i) {nrow(DT) - sum(DT$CD >= max(DT$CD,na.rm=T)*i)}),Cutoffs=cutoffs)	
#mylines <- data.table(RHB=sapply(cutoffs,function(i) {nrow(DT) - sum(DT >= max(DT,na.rm=T)*i)}),Cutoffs=cutoffs)	
	}) 

	# log the count values
	DT$CD <- log10(DT[,"CD"])

	# create an empty ggplot object from the data table
	g <- ggplot(data=DT,aes_string(x=seq(1,nrow(DT)),y="CD"))

	# remove plot background and etc.
	g <- g + theme_classic_thin()

	# a colour pallete that can be descriminated by the colur blind 
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	# plot cumulative reads
	g <- g + geom_line(size=1.5) + scale_colour_manual(values=cbbPalette) 

	if(exists("mylines")) {
		# add cutoff lines
		g <- g + geom_vline(data=mylines,aes(xintercept=RHB),colour=cbbPalette[2:(length(cutoffs)+1)])
	
		# label cutoff lines
		g <- g + geom_text(data = mylines, aes(label=Cutoffs,x=RHB, y = -Inf),angle = 90, inherit.aes = F, hjust = -.5, vjust = -.5)
	}

	# add axis lables
	g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")

	# print the plot
	g
}


plotRarefaction <- function(X,cutOff) {
  
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")  
  DT <- as.data.table(X)
  DT <- DT[,lapply(.SD,sort,decreasing=T)]
  DT <- cumsum(DT)
  
  
  DT <- DT[,lapply(.SD,log10)]
  DT[,Count:=1:nrow(DT)]
  
  DT <- melt(DT,id.vars="Count")
  
  g <- ggplot(DT[Count<=cutOff,],aes(x=Count,y=value,colour=variable))    
  g + 
    geom_line(size = 1.5)  + 
    scale_colour_manual(values = cbPalette) + 
    theme_classic_thin() %+replace% theme(legend.position = "none") + 
    xlab("Number of OTUs") + ylab(expression(~Log[10]~" sequence count"))
}
plotHeatmap <- 
function (dist,
	coldata,
	axis.labels=F,
	textSize=11,
	textFont="Helvetica"

) {
	require(ggplot2)
	#require(reshape)
	
	if (!missing(coldata)) {
		dist<-merge(dist,coldata,by="row.names",sort=F)
	}

	h1 <- melt(dist)
	colnames(h1) <-c("X","Y","Unidist")
	g <- ggplot(h1,aes(x=X,y=Y,fill=Unidist))
	if(axis.labels) {
		g <- g + theme_classic(textSize,textFont) %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1,),axis.title=element_blank())
	} else 	{
		g <- g + theme_classic(textSize,textFont) %+replace% theme(axis.text=element_blank(),axis.ticks=element_blank())
	}	

	g <- g + geom_raster()

	g <- g + scale_x_discrete("Var1")
	g <- g + scale_y_discrete("Var2")

	#g <- g + geom_tile(aes(fill=Unidist),colour="black")
	# g <- g + scale_fill_gradient(low="#000033", high = "#66CCFF", na.value = "black") # blue scale

	g <- g + scale_fill_gradient(low="black", high = "lightgrey", na.value = "white") # greyscale
	
	return(g)

}



plotOTUs <-
function (
	countData,
	colData,
	design="time",
	colour="Treatment",
	plotsPerPage=5,
	facet=formula(~OTU),
	Ylab="log_counts",
	line="smooth",
	returnData=F
) {
	suppressPackageStartupMessages(require(viridis))

	d <- data.frame(t(countData),colData)

	d <- melt(d,
		measure.vars=colnames(d)[1:(ncol(t(countData)))],
		id.vars = colnames(d)[(ncol(t(countData))+1):ncol(d)],
		variable.name = "OTU", 
		value.name = Ylab)

	d[[Ylab]] <- d[[Ylab]]+abs(min(d[[Ylab]]))

	ymax <- max(d[[Ylab]])
	allVars <- unique(d$OTU)
	noVars <- length(allVars)
	plotSequence <- c(seq(0, noVars-1, by = plotsPerPage), noVars)
	
	if(returnData) { return(d)}
	
	sapply(seq(2,length(plotSequence)),function(i) {
		start <- plotSequence[i-1] + 1
		end <- plotSequence[i]
		tmp <- d[d$OTU %in% allVars[start:end],]
		cat(unique(tmp$OTU), "\n")
		g <- ggplot(data=tmp,aes_string(y=Ylab, x=design,colour=colour),ylim=c(0,ymax))
		g <- g + theme_classic_thin(base_size = 16) %+replace% theme(panel.border=element_rect(colour="black",size=0.25,fill=NA),legend.position="bottom")
		g <- g + scale_colour_viridis(discrete=TRUE)
		g <- g + facet_grid(facet,scales="free_x")
		g <- g + geom_point(size=2)
		if(line=="smooth"){
			g <- g + stat_smooth(method=locfit, formula=y~lp(x),se=F)
		}else {
			g <- g + geom_line()
		}
		print(g)
	})

}	#' Ordination plot
#'
#' plotOrd can help in making pretty PCA/Ordination plots.
#' 
#' This is a function for making simple ggplot ordination plots from two column data with 
#' additional metadata. Based on \code{DESeq2::plotPCA}, but with extra functionality for labeling
#' points, and etc. All the additional functionality can be applied post return of ggplot object.
#' 
#'
#' @param obj a dataframe containing xy coordinates, e.g. principal component
#'   scores.
#' @param colData dataframe containing sample metadata.
#' @param design column(s) of colData to discriminate sample type (colour).
#' @param shapes column(s) of colData to discriminate sample type (shape).
#' @param label column(s) of colData to use for sample labeling.
#' @param facet column(s) of colData. Adds a faceting column to the 
#'  returned ggplot object. To use call g + facet_wrap(~facet). 
#' @param plot Plot either [point] or [label] (Default = "point").
#' @param labelSize Text size for labels(Default = 4).
#' @param labelPosition Label position relative to point(Default = c(-1.5,0.5)).
#' @param sublabels a numeric vector of labels to remove (Default = F).
#' @param cluster Set to turn on clustering, value is stat_ellipse confidence.
#' @param continuous T/F whether design is a continuos variable (default FALSE).
#' @param colourScale Vector used for continuous colour scales (Low to High)
#'   (Default = c(low="red", high="yellow")) #greyscale low="#000000", high="#DCDCDC".
#' @param cbPalette Use a predefined colour blind palette.
#' 	Max eight factors allowable in design (Default = F).
#' @param pointSize The size of plot points (Default = 2).
#' @param textSize The base text size for the graph (Default = 11).
#' @param textFont The base text font for the graph (Default = "Helvetica").
#' @param xlims, ylims Numeric vectors of axis limits, 
#'  e.g. c(-8,8) (Default unlimited).
#' @param legend Position of legend. 
#'   Set to "none" to remove legend (Default = "right").
#' @param legendDesign Display legend for design (Default=True).
#' @param legendShape Display legend for shapes (Default=True).
#' @param title Title (Default is to not use a title).
#' @param xlabel, ylabel Set axis labels.
#' @param axes Columns of obj to use for plotting (Default = c(1,2)).
#' @param alpha Numeric value, "prettifies" points by adding an extra outer circle 
#'   with given alpha value.
#' @param exclude vector of points to exclude (Default show all points).
#' @param noPlot T/F if set do not plot return list of: 
#'   [1] selected obj axes and [2] aesthetics (Default FALSE).
#' @param ... additional parameters (unused).
#' @return A ggplot scatter plot of the axes taken from obj, 
#' colours as per design and shapes as per shapes (unless noPlot set to TRUE).
#' @examples
#' d <- data.frame(PCA1=runif(10,-8,8),PCA2=runif(10,-4,6))
#' m <- data.frame(Sample=seq(1,10),
#'	Condition=rep(c("H","D"),5),
#'	Site=c(rep("A",5),rep("B",5))) 
#' 
#' plotOrd(obj=d,colData=m,design="Condition")
#'
#' plotOrd(obj=d,colData=m,design="Condition",shapes="Site", alpha=0.75)
#'
#' plotOrd(obj=d,colData=m,design="Condition",xlims=c(-2,2), label="Sample")
#'
#' plotOrd(obj=d,colData=m,design="Condition",pointSize=3, alpha=0.75, textSize=16, cluster=0.75)  

plotOrd <- function (
	obj,
	colData,
	design=NULL, # column(s) of colData
	shapes=NULL, # column(s) of colData
	label=NULL,  # column(s) of colData
	facet=NULL,  # column(s) of colData. This doesn't add a layer to the graph just adds facet as a data column, to use call g + facet_wrap(~facet) 
	plot="point", # or "label"
	labelSize=4, # for text label to point
	labelPosition=c(-1.5,0.5), # for text label to point
	sublabels=F,
	cluster=NULL,
	continuous=F,
	colourScale=c(low="red", high="yellow"), #greyscale low="#000000", high="#DCDCDC"
	cbPalette=F,
	pointSize=2,
	textSize=11,
	textFont="Helvetica",
	xlims=NULL,
	ylims=NULL,
	legend="right",
	legendDesign=T,
	legendShape=T,
	title=NULL,
	xlabel,
	ylabel,
	axes=c(1,2),
	alpha=NULL,

	exclude=T, # sometimes it can be useful to include a vector of points to exclude
	noPlot=F, 
	...
) {

	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(viridis))
	
	if(missing(obj)|missing(colData)) return(print("Error : please specify both obj and colData"))

	ef <- function(X) {
		cat("WARNING: Incorrect columns specified \"",X,"\"",sep="")
		return(NULL)
	}

	design <- tryCatch(colnames(colData[,design,drop=F]),error=function(e)ef(design))
	shapes <- tryCatch(colnames(colData[,shapes,drop=F]),error=function(e)ef(shapes))
	label  <- tryCatch(colnames(colData[,label,drop=F]), error=function(e)ef(label))
	facet  <- tryCatch(colnames(colData[,facet,drop=F]), error=function(e)ef(facet))

	# check if label is set if using label as a plot type
	ll<-T
	if(tolower(plot)=="label") {
		ll<-F
		shapes<-NULL
		alpha<-NULL
		if (!length(label)) {
			print("No label column specified defaulting to first column of colData")
			label=colnames(colData)[1]
		}
 
	}

	obj     <- obj[!rownames(obj)%in%exclude,]
	colData <- colData[!rownames(colData)%in%exclude,]
	
	if(!is.null(title)){if(title=="debugging"){invisible(mapply(assign, ls(),mget(ls()), MoreArgs=list(envir = globalenv())));return(title)}}

	d <- as.data.frame(obj[,axes])

	if(!is.null(cluster)) {
	#	km <- kmeans(obj,...)
	#	colData$Cluster<-as.factor(km$cluster)
	} 

	x = colnames(d)[1]
	y = colnames(d)[2]

	aes_map <-aes_string(x=x,y=y)

	if (length(design)) {
		colour <- if (length(design) > 1) {
			factor(apply(as.data.frame(colData[, design,drop = FALSE]), 1, paste, collapse = " : "))
		}
		else {
			as.factor(colData[[design]])
		}
		d <- cbind(d,colour=colour)
		if(continuous) {
			aes_map <- modifyList(aes_map,aes(colour=as.number(colour)))
		} else {
			aes_map <- modifyList(aes_map,aes_string(colour="colour"))
		}
	}

	if (length(shapes)) {
		shape <- if (length(shapes) > 1) {
			factor(apply(as.data.frame(colData[, shapes,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[shapes]])
		}
		d <- cbind(d,shapes = shape)
		aes_map <- modifyList(aes_map,aes_string(shape="shapes"))
	}

	if(length(label)) {
		label <- if (length(label) > 1) {
			factor(apply(as.data.frame(colData[, label,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[label]])
		}
		d <- cbind(d,label = label)
		d$label[sublabels] <- NA
		aes_map <- modifyList(aes_map,aes_string(label="label"))
	}

	if(length(facet)) {
		facet <- if (length(facet) > 1) {
			factor(apply(as.data.frame(colData[, facet,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[facet]])
		}
		d <- cbind(d,facet = facet)
		aes_map <- modifyList(aes_map,aes_string(facet="facet"))
	}


	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	if(noPlot) return(list(data=d,aes=aes_map))

	### ggplot ###
	g <- ggplot(data=d,aes_map) 
	
	g <- g + if(tolower(plot)=="label"){
		geom_label(size=(labelSize))
	} else {
		geom_point(size=pointSize,na.rm = TRUE)
	}
	g <- g + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)

	g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.position=legend)

	if(!is.null(design)) {
		if(continuous) {
			g <- g + scale_colour_gradient(low=colourScale[1], high=colourScale[2],name=design,guide=legendDesign)
		} else {
			if(cbPalette) {
				g<-g+scale_colour_manual(values=cbbPalette,guide=legendDesign) + guides(colour=guide_legend(title=design))
			} else {
				g<-g+scale_colour_viridis(discrete=TRUE,guide=legendDesign) + guides(colour=guide_legend(title=design))
			}
		}
	}

	if (!is.null(shapes)) {
		g <- g + scale_shape_discrete(name=shapes)
		if(!legendShape) {g <- g + guides(shapes="none")}
	}

	if(!is.null(alpha)) g <-g+ geom_point(size=(pointSize+(pointSize*1.5)),alpha=alpha)

	if(!is.null(cluster)) {
			g<-g+stat_ellipse(geom="polygon", level=cluster, alpha=0.2)
	}

	if (!missing(xlabel)) {g <- g + xlab(xlabel)}
	if (!missing(ylabel)) {g <- g + ylab(ylabel)}

	if(length(label)&ll) { 
		g <- g + geom_text(size=(labelSize), vjust=labelPosition[1], hjust=labelPosition[2],check_overlap = TRUE)
	}

	return(g)
}
plotPCA <- function (	
	obj, 
	design = "condition",
	labelby,
	ntop = 500,
	pcx = 1,
	pcy = 2, 
	returnData = FALSE,
	cofix=F,
	trans=T,
	addLabels=F,
	filter=1,
	filterFun=function(o,f){return(o)},
	transform=function(o,design,...) 
	{	
		suppressPackageStartupMessages(require(DESeq2))
		dots <- list(...)
		if(!is.null(dots$calcFactors)) {
			calcFactors <- dots$calcFactors
			dots$calcFactors<-NULL
			dots$object <- phylo_to_des(o,as.formula(design),calcFactors=calcFactors)
			assay(do.call(varianceStabilizingTransformation, dots))
		} else {
			assay(varianceStabilizingTransformation(phylo_to_des(o,as.formula(design)),...))
		}
	},...
) {
	suppressPackageStartupMessages(require(genefilter))
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(DESeq2))

	if(trans) {
		obj@otu_table@.Data <- transform(obj,as.formula(paste0("~",design)),...)
	}
	
	obj <- filterFun(obj,filter)    	

 	if (returnData) {
 		#d <- pca$x
 		#attr(d, "percentVar") <- percentVar
 		pca <- prcomp(t(otu_table(obj)))
 		pca$percentVar <- pca$sdev^2/sum(pca$sdev^2)
 		return(pca)
	}
 
 	rv <- rowVars(otu_table(obj))
 
	colData <- sample_data(obj)
	#suppressWarnings(as.data.frame(as.matrix(obj@sam_data)))
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	pca <- prcomp(t(otu_table(obj)[select, ]))
	percentVar <- pca$sdev^2/sum(pca$sdev^2)
	
	if (!all(design %in% names(colData))) {
		stop("the argument 'design' should specify columns of colData")
	}
	design.df <- as.data.frame(colData[, design,drop = FALSE])
	group <- if (length(design) > 1) {
		factor(apply(design.df, 1, paste, collapse = " : "))
	}
	else {
		as.factor(sample_data(obj)[[design]])
	}
	
	if (!missing(labelby)) {
		shape <- as.factor(sample_data(obj)[[labelby]])
		d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,design.df,shape = shape)
		colnames(d)[grep("shape", colnames(d))] <- labelby
	} else {
		d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,design.df)
	}

	colnames(d)[grep("group", colnames(d))] <- design

	if(cofix) {
		d[,1] <- d[,1] * percentVar[pcx]
		d[,2] <- d[,2] * percentVar[pcy]
	}

	g <- ggplot()
	g <- g + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
	if (!missing(labelby)) {
	g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group, shape=shape),size=3)
	g <- g + scale_shape_discrete(name=labelby)
	} else {
	g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group),size=3)
	}
	g <- g + scale_colour_discrete(name=design)
	g <- g + xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance"))
	g <- g + ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
	return(g)
}

plotPCAWithLabels <- function (
	object, 
	intgroup = "condition", 
	ntop = 500,
	pcx = 1,
	pcy = 2, 
	returnData = FALSE
){
    suppressPackageStartupMessages(require(genefilter))
    suppressPackageStartupMessages(require(ggplot2))
    
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(object@colData))) {
        stop("the argument 'intgroup' should specify columns of colData")
    }
    intgroup.df <- as.data.frame(object@colData[, intgroup,
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,
        intgroup.df, name = object$label)
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }

    ggplot() +
    geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group),size=3) +
    geom_text(data=d, mapping=aes(x=PC1, y=PC2, label=name,colour=group), size=3, vjust=2, hjust=0.5) +
    xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance")) +
    ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
}

plotTaxa <- function(
	obj=mybiom, 	# obj (phloseq) a phyloseq object which must include taxonomy and sample data (or alternatively an S3 list)
	taxon="phylum", 	# taxon (str) is the taxonomic level of interest
	design="all", 		# condition (str) describes how the samples should be grouped (must be column of sample data)
	proportional=T,	# proportional (bool) whether the graph should use proportional or absolute values
	cutoff=1, 	# cutoff (double) for proportional graphs. 
	topn=0, 		# topn (int)taxons to display (by total reads) for non-prortional graphs. T
	others=T, 	# combine values less than cutoff/topn into group "other"
	ordered=F, 	# order by value (max to min)
	stack="stack",	# stacked graph, set to "dodge" for side by side
	type=2, 		# type: (1) by sample (2) by taxonomy 
	fixed=F, 		# fixed is a ggplot parameter to apply coord_fixed(ratio = 0.1)
	ncol=1, 		# ncol is a ggplot paramter to use n columns for the legend
	trans=T,		# set to False if using original/prior transformed counts (useful for larger OTU tables)
	bw=F,		# set graph to black and white 
	ylab="",
	ylims=NULL,	# 
	margins=c(0.2,0.2,0.2,1.5),
	textSize=14,
	returnData=F,
	legend=T,
	NEW=T,
    conf=0.65,
	transform=function(o,design,...) 
	{	
		suppressPackageStartupMessages(require(DESeq2))
		dots <- list(...)
		
		if(!is.null(dots$calcFactors)) {
			calcFactors <- dots$calcFactors
			dots$calcFactors<-NULL
			if(length(dots)>=1) {
				 assay(varianceStabilizingTransformation(ubiom_to_des(o,design=design,calcFactors=calcFactors),unlist(dots)))
			} else {
				 assay(varianceStabilizingTransformation(ubiom_to_des(o,desing=design,calcFactors=calcFactors)))
			}
		} else {
			if(NEW) {
				assay(varianceStabilizingTransformation(o$dds))
			} else {
				assay(varianceStabilizingTransformation(ubiom_to_des(o,design=as.formula(design)),...))
			}
		}
	}, # data transformation function 
	... # arguments to pass to transform function (obviously they could just be set in the function, but this looks neater)
) {
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(scales))

	if(isS4(obj)) {
		obj <- phylo_to_ubiom(obj)
	} 


	if(trans) {
		temp <- design
		idx <- grep(design,colnames(obj[[3]]))
		if(length(unique(obj[[3]][idx]))<=1) {
			design<-1
		}
		obj[[1]] <- as.data.frame(transform(obj,as.formula(paste0("~",design)),...))
		#obj[[1]] <- as.data.frame(transform(obj,as.formula("~1"),...))
		design<-temp
	}
	#obj[[1]][obj[[1]]] <- obj[[1]][obj[[1]]]+abs(min(obj[[1]][obj[[1]]]))

	obj[[1]][obj[[1]]<0] <- 0


	if(NEW) {
		obj <- list(countData=obj$countData,taxData=obj$taxData,colData=obj$colData)
	}

	

	taxa_sum <- sumTaxa(obj,taxon=taxon,design=design,conf=conf)
	taxa_sum$taxon[grep("\\(",taxa_sum$taxon)] <- taxa_sum$taxon[sub("\\(.*"," incertae sedis",taxa_sum$taxon)]

	if(!topn) {
		obj[[3]]$MLUflop <- 1 #assigns the MLU flop digit
		tx <- sumTaxa(obj,taxon=taxon,"MLUflop")
		tx[,-1] <- prop.table(as.matrix(tx[,-1]),2)*100
		txk <- tx[tx[,2]>=cutoff,1]
	} else {
		taxa_sum[,ncol(taxa_sum)+1]<- 0
		taxa_sum <- taxa_sum[order(rowSums(taxa_sum[,-1]),decreasing=T),]
		taxa_sum <- taxa_sum[,-ncol(taxa_sum)]	
		txk <- taxa_sum[1:topn,1]
	}
	
	if(proportional) {
		taxa_sum[,-1] <- prop.table(as.matrix(taxa_sum[,-1]),2)*100
	}

	taxa_cut <- taxa_sum[taxa_sum[,1]%in%txk,]
	taxa_cut <- taxa_cut[order(taxa_cut[,1],decreasing=T),]
	if(others) {
		taxa_cut <- rbind(taxa_cut,setNames(data.frame(x="others" ,t(colSums(taxa_sum[!taxa_sum[,1]%in%txk,-1]))),names(taxa_cut)))
	}
	taxa_cut <- na.omit(taxa_cut)
	taxa_cut[,1] <- as.factor(taxa_cut[,1])
	if(ordered) {
		taxa_cut[,ncol(taxa_cut)+1] <- 0
		taxa_cut[,1] <- reorder(taxa_cut[,1],-rowSums(taxa_cut[,-1]))
		taxa_cut <- taxa_cut[,-ncol(taxa_cut)]
	}

	if(returnData) {
		return(taxa_cut)
	}

	md2 <- melt(taxa_cut,id=colnames(taxa_cut)[1])

	md2$variable <- factor(md2$variable, levels=levels(md2$variable)[order(levels(md2$variable))]  )

	md2$value <- as.numeric(md2$value)

	if (type==1) {
		g <- ggplot(md2,aes_string(x=md2[,2],y=md2[,3],fill=taxon))
	} else {
		colnames(md2) <- c("taxa",design,"value")
		g <- ggplot(md2,aes_string(x=as.factor(md2[,1]),y=md2[,3],fill=design))
	}
	

	if(bw) {
		g<-g+geom_bar(stat="identity",colour="white",fill="black",position="stack")
	} else {
		g <- g + geom_bar(stat="identity",colour="white",position="stack")		
	}

	g <- g  + xlab("")

	if (fixed) {
		g <- g  + coord_fixed(ratio = 0.1)
	} 

	scaleFUN<-function(x) sprintf("%.0f", x)

	g <- g + scale_y_continuous(expand = c(0, 0),labels = scaleFUN,limits=ylims)
	g <- g + ylab(ylab)
	g <- g + guides(fill=guide_legend(ncol=ncol))
	g <- g + theme_blank()
	g <- g + theme(
		axis.text.x = element_text(angle = 45, hjust = 1,size=textSize),
		plot.margin=unit(margins,"cm"), 
	    axis.line.y = element_line(colour = "black",size=1),
		axis.ticks.x=element_blank(),
		text=element_text(size=textSize),
		axis.title.y=element_text(size=(textSize-2))

	)
	if (type==1) {
		g <- g + theme(legend.text=element_text(face="italic"))
	} else if (type==2) {
		#g <- g + theme(axis.text.x=element_text(face="italic"))
	}
	if(!legend) {
		g<- g+guides(fill=FALSE)
	}
	
	return(g)
}


plot_ma <- function
(
	fitObj,
	xlim=c(-6,6),
	textSize=16,
	textFont="Helvetica",
	pointSize=3,
	legend=F,
	LOG=2,
	crush=T
)
{
	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	d <- as.data.table(fitObj)
	colnames(d) <- c("log2FoldChange","baseMean","padj")
	d$group<-"Not sig"
	d$group[d$padj<=0.05]<- "p <= 0.05" 
	d$group[abs(d$log2FoldChange)>1] <- "FC >= 2"	
	d$group[(d$padj<=0.05&(abs(d$log2FoldChange)>1))] <- "p <= 0.05 &\nFC >= 2" 
	d$shape<-16
	d$group<-as.factor(d$group)
	d$group<-factor(d$group,levels(d$group)[c(2,3,1,4)])

	if(crush){
		d$shape[d$log2FoldChange<xlim[1]] <- 25
		d$log2FoldChange[d$log2FoldChange<xlim[1]] <-xlim[1]
		d$shape[d$log2FoldChange>xlim[2]] <- 24
		d$log2FoldChange[d$log2FoldChange>xlim[2]] <-xlim[2]
	}

	g <- ggplot(data=d,aes(x=log2FoldChange,y=log(baseMean,LOG),colour=group,shape=shape))
	
	if(!legend) {
		g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.position="none")
	} else {
		g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.title=element_blank())
	}

	g <- g + scale_shape_identity() 
	g <- g + geom_point(size=pointSize)
	g <- g + scale_colour_manual(values=cbbPalette)
	g <- g + xlab(expression("Log"[2]*" fold change"))
	g <- g + ylab(expression("Log"[2]*" mean expression"))
	g <- g + xlim(xlim)
	g <- g + expand_limits(x = xlim[1], y = 5)
	g <- g + coord_flip()
	return(g)
}
##phyloseq alpha diversity function modified to allow return of data

plot_richness <-
function (physeq, x = "samples", color = NULL, shape = NULL, size=2, textSize=14,
    title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL,
    sortby = NULL,
    limits,
    returnData=F
    )
{
	library(phyloseq)
    erDF = estimate_richness(physeq, split = TRUE, measures = measures)
    measures = colnames(erDF)
    ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
    measures = measures[!measures %in% ses]
    if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
        DF <- data.frame(erDF, sample_data(physeq))
    }
    else {
        DF <- data.frame(erDF)
    }
    if (!"samples" %in% colnames(DF)) {
        DF$samples <- sample_names(physeq)
    }
    if (!is.null(x)) {
        if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
            x <- "samples"
        }
    }
    else {
        x <- "samples"
    }
    if(returnData) {
        return(DF)
    }

#	DF[,3:4] <- round(DF[,3:4],2)

    if(!missing(limits)) {
	  DF <- DF[DF[,1]>=limits[1],]
	  DF<-DF[!rowSums(DF[,1:(ncol(DF)-ncol(sample_data(physeq))-1)]>limits[2]),]
    }

    mdf = melt(DF, measure.vars = measures)
    mdf$se <- NA_integer_
    if (length(ses) > 0) {
        selabs = ses
        names(selabs) <- substr(selabs, 4, 100)
        substr(names(selabs), 1, 1) <- toupper(substr(names(selabs),
            1, 1))
        mdf$wse <- sapply(as.character(mdf$variable), function(i,
            selabs) {
            selabs[i]
        }, selabs)
        for (i in 1:nrow(mdf)) {
            if (!is.na(mdf[i, "wse"])) {
                mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
            }
        }
        mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
    }
    if (!is.null(measures)) {
        if (any(measures %in% as.character(mdf$variable))) {
            mdf <- mdf[as.character(mdf$variable) %in% measures,
                ]
        }
        else {
            warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
        }
    }
    if (!is.null(shsi)) {
        warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
    }
    if (!is.null(sortby)) {
        if (!all(sortby %in% levels(mdf$variable))) {
            warning("`sortby` argument not among `measures`. Ignored.")
        }
        if (!is.discrete(mdf[, x])) {
            warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
        }
        if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[,
            x])) {
            wh.sortby = which(mdf$variable %in% sortby)
            mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby,
                "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE,
                simplify = TRUE))))
        }
    }

    richness_map = aes_string(x = x, y = "value", colour = color,
        shape = shape)
    p = ggplot(mdf, richness_map) + geom_point(na.rm = TRUE,position = position_dodge(width = 0.75),size=size)
    if (any(!is.na(mdf[, "se"]))) {
        p = p + geom_errorbar(aes(ymax = value + se, ymin = value -
            se), width = 0.2,position = position_dodge(width = 0.75))
    }
    p = p +theme_bw(base_size = textSize, base_family = "Helvetica")
#axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
    p = p + theme(
      panel.grid.major.x = element_line(size = .5, color = "lightgrey"),
      panel.grid.major.y =element_blank())
    p = p + ylab("Alpha Diversity Measure")
    p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}


# cut down version of above - useful for projects without phyloseq objects
plot_alpha <-
function (	
  countData,
  colData,
  design=NULL,
  colour=NULL,
  shape=NULL,
  returnData=F,
  limits,
  measures=NULL,
  discrete=T,
  legend="right",
  pointSize=2,
  cbPalette=F,
  type="dot", # or box	
  ...
){
	suppressPackageStartupMessages(require(viridis))
	suppressPackageStartupMessages(require(vegan))

	simpleCap <- function(x) {
        	s <- strsplit(x, " ")[[1]]
        	paste(toupper(substring(s, 1,1)), substring(s, 2),
        	sep="", collapse=" ")
    	}

	retain<-measures
	# function to convet input to integer values with ceiling for values < 1
	alpha_counts <- function(X) {
		X[X==0]<-NA
		X[X<1] <- 1
		X[is.na(X)] <- 0
		return(round(X,0))
	}

	# convert counts to integers and transpose
	OTU <- t(alpha_counts(countData))

	# calculate diversity measures
	all_alpha <- data.table(
		t(vegan::estimateR(OTU)),
		shannon = vegan::diversity(OTU, index = "shannon"),
		simpson = vegan::diversity(OTU, index = "simpson"),
		keep.rownames="Samples"
	)


	if(returnData) {return(all_alpha)}

	# set colData to a data table (no ides why this is to data frame first...)
	colData <- as.data.table(as.data.frame(colData),keep.rownames="Samples")

	# get column names of all_alpha
	measures = colnames(all_alpha[,-1])

	# get standard error columns (Chao1 and ACE)
	ses = colnames(all_alpha)[grep("^se\\.", colnames(all_alpha))]

	# temp holder for se labels
	selabs = ses

	# rename se labels to the same as the corresponding measures column
	names(selabs) <- sub("se","S",selabs)

	# set measures to (measures - se) columns
	measures = measures[!measures %in% ses]

	all_alpha <- as.data.table(inner_join(all_alpha,colData))

	id.vars <- colnames(colData)

	# melt the data table by Sample
	mdf <- melt(all_alpha,meaures.vars=meaures,id.vars=id.vars,variable.factor = TRUE)

	X<-all_alpha[,c(id.vars,selabs),with=F]
	names(X) <- c(id.vars,names(selabs))
	X <- melt(X,id.vars=id.vars,value="se",variable.factor = TRUE)

	mdf <- left_join(mdf,X)

	# remove se rows
	mdf <- mdf[as.character(mdf$variable) %in% measures,]
	
	# capitalise indices
	mdf$variable <- sub("S\\.","",mdf$variable)
	mdf$variable <- sapply(mdf$variable,simpleCap)

	# refactor variable
	mdf$variable <- sub("Obs","Observed",mdf$variable)
	mdf$variable <- as.factor(mdf$variable)
	mdf$value <- as.numeric(mdf$value)

	# retain indices were interested in (o.k. this would be better done before calculating them)
	if (!is.null(retain)) {
	    retain<-sapply(tolower(retain),simpleCap)
    	mdf <- mdf[as.character(mdf$variable) %in% retain,]
    	mdf <- droplevels(mdf)
    	mdf$variable <- factor(mdf$variable, levels = retain)
    }
	
	arguments <- list(...)
	if("debugging"%in%names(arguments))if(arguments$debugging) return(mdf)

	if(!missing(limits)) {
		mdf <- mdf[!(mdf$variable==limits[3]&(mdf$value<as.numeric(limits[1])|mdf$value>as.numeric(limits[2]))),]
	}

	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	aes_map <-aes_string(y="value",variable="variable")
	if(length(design)) aes_map <- modifyList(aes_map,aes_string(x=design))
	if(length(colour)) aes_map <- modifyList(aes_map,aes_string(colour=colour))
	if(length(shape)) aes_map <- modifyList(aes_map,aes_string(shape=shape))

	# create ggplot object
	g <- ggplot(data=mdf,aes_map)

	# add a theme
	g <- g + theme_classic_thin() %+replace% theme(
		panel.border = element_rect(colour = "black", fill=NA, size=0.5),
		axis.text.x = element_text(angle = -90, vjust = 0.5,hjust = 0),
		legend.position=legend
	)

	if(type=="box") {
	  g <- g + geom_boxplot( position = position_dodge(width = 1))
	} else {
	  # add points
	  g <- g + geom_point(na.rm = TRUE,position = position_dodge(width = 0.5),size=pointSize)

	  # add error bars
	  #g <- g + geom_errorbar(aes(ymax = value + se, ymin = value -  se), width = 0.5,position = position_dodge(width = 0.5))
	}
	# add y label
	g <- g + ylab("Alpha Diversity Measure")

	# change colours to viridis
	if(cbPalette) {
		g<-g+scale_colour_manual(values=cbbPalette) #+ guides(colour=guide_legend(title=design))
	} else {
		g<-g+scale_colour_viridis(discrete=TRUE) #+ guides(colour=guide_legend(title=design))
	}

	# Add heading to each graph
	g <- g + facet_wrap(~variable, nrow = 1,scales="free_y")

	g

}

# modified phyloseq ordination plot (allows rescaling by percentage variation in the axes)
plot_ordination <-
function (
  physeq,
  ordination,
  type = "samples",
  axes = 1:2,
  color = NULL,
  shape = NULL,
  label = NULL,
  title = NULL,
  justDF = FALSE,
  continuous=F,
  colourScale=c(low="red", high="yellow"),
  cbPalette=F,
  xlims=NULL,
  ylims=NULL,
  alpha=NULL,
  rescale=T,
  ...
)
{
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(viridis))

    if (length(type) > 1) {
        warning("`type` can only be a single option,\\n            but more than one provided. Using only the first.")
        type <- type[[1]]
    }
    if (length(color) > 1) {
        warning("The `color` variable argument should have length equal to 1.",
            "Taking first value.")
        color = color[[1]][1]
    }
    if (length(shape) > 1) {
        warning("The `shape` variable argument should have length equal to 1.",
            "Taking first value.")
        shape = shape[[1]][1]
    }
    if (length(label) > 1) {
        warning("The `label` variable argument should have length equal to 1.",
            "Taking first value.")
        label = label[[1]][1]
    }
    official_types = c("sites", "species", "biplot", "split",
        "scree")
    if (!inherits(physeq, "phyloseq")) {
        if (inherits(physeq, "character")) {
            if (physeq == "list") {
                return(official_types)
            }
        }
        warning("Full functionality requires `physeq` be phyloseq-class ",
            "with multiple components.")
    }
    type = gsub("^.*site[s]*.*$", "sites", type, ignore.case = TRUE)
    type = gsub("^.*sample[s]*.*$", "sites", type, ignore.case = TRUE)
    type = gsub("^.*species.*$", "species", type, ignore.case = TRUE)
    type = gsub("^.*taxa.*$", "species", type, ignore.case = TRUE)
    type = gsub("^.*OTU[s]*.*$", "species", type, ignore.case = TRUE)
    type = gsub("^.*biplot[s]*.*$", "biplot", type, ignore.case = TRUE)
    type = gsub("^.*split[s]*.*$", "split", type, ignore.case = TRUE)
    type = gsub("^.*scree[s]*.*$", "scree", type, ignore.case = TRUE)
    if (!type %in% official_types) {
        warning("type argument not supported. `type` set to 'samples'.\\n",
            "See `plot_ordination('list')`")
        type <- "sites"
    }
    if (type %in% c("scree")) {
        return(plot_scree(ordination, title = title))
    }
    is_empty = function(x) {
        length(x) < 2 | suppressWarnings(all(is.na(x)))
    }
    specDF = siteDF = NULL
    trash1 = try({
        siteDF <- scores(ordination, choices = axes, display = "sites",
            physeq = physeq)
    }, silent = TRUE)
    trash2 = try({
        specDF <- scores(ordination, choices = axes, display = "species",
            physeq = physeq)
    }, silent = TRUE)
    siteSampIntx = length(intersect(rownames(siteDF), sample_names(physeq)))
    siteTaxaIntx = length(intersect(rownames(siteDF), taxa_names(physeq)))
    specSampIntx = length(intersect(rownames(specDF), sample_names(physeq)))
    specTaxaIntx = length(intersect(rownames(specDF), taxa_names(physeq)))
    if (siteSampIntx < specSampIntx & specTaxaIntx < siteTaxaIntx) {
        co = specDF
        specDF <- siteDF
        siteDF <- co
        rm(co)
    }
    else {
        if (siteSampIntx < specSampIntx) {
            siteDF <- specDF
            specDF <- NULL
        }
        if (specTaxaIntx < siteTaxaIntx) {
            specDF <- siteDF
            siteDF <- NULL
        }
    }
    if (is_empty(siteDF) & is_empty(specDF)) {
        warning("Could not obtain coordinates from the provided `ordination`. \\n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
        return(NULL)
    }
    if (is_empty(specDF) & type != "sites") {
        message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
        specDF <- data.frame(wascores(siteDF, w = veganifyOTU(physeq)),
            stringsAsFactors = FALSE)
    }
    if (is_empty(siteDF) & type != "species") {
        message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
        siteDF <- data.frame(wascores(specDF, w = t(veganifyOTU(physeq))),
            stringsAsFactors = FALSE)
    }
    specTaxaIntx <- siteSampIntx <- NULL
    siteSampIntx <- length(intersect(rownames(siteDF), sample_names(physeq)))
    specTaxaIntx <- length(intersect(rownames(specDF), taxa_names(physeq)))
    if (siteSampIntx < 1L & !is_empty(siteDF)) {
        warning("`Ordination site/sample coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
        siteDF <- NULL
    }
    if (specTaxaIntx < 1L & !is_empty(specDF)) {
        warning("`Ordination species/OTU/taxa coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
        specDF <- NULL
    }
    if (is_empty(siteDF) & is_empty(specDF)) {
        warning("Could not obtain coordinates from the provided `ordination`. \\n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
        return(NULL)
    }
    if (type %in% c("biplot", "split") & (is_empty(siteDF) |
        is_empty(specDF))) {
        if (is_empty(siteDF)) {
            warning("Could not access/evaluate site/sample coordinates. Switching type to 'species'")
            type <- "species"
        }
        if (is_empty(specDF)) {
            warning("Could not access/evaluate species/taxa/OTU coordinates. Switching type to 'sites'")
            type <- "sites"
        }
    }

	if (length(extract_eigenvalue(ordination)[axes]) > 0) {
        eigvec = extract_eigenvalue(ordination)
        fracvar = eigvec[axes]/sum(eigvec)
        percvar = round(100 * fracvar, 1)
    } else {
		percvar = 1
	}

    if (type != "species") {
		if(rescale) siteDF <- t(t(siteDF)*percvar)
        sdf = NULL
        sdf = data.frame(access(physeq, slot = "sam_data"), stringsAsFactors = FALSE)
        if (!is_empty(sdf) & !is_empty(siteDF)) {
            siteDF <- cbind(siteDF, sdf[rownames(siteDF), ])
        }
    }
    if (type != "sites") {
        tdf = NULL
        tdf = data.frame(access(physeq, slot = "tax_table"),
            stringsAsFactors = FALSE)
        if (!is_empty(tdf) & !is_empty(specDF)) {
            specDF = cbind(specDF, tdf[rownames(specDF), ])
        }
    }
    if (!inherits(siteDF, "data.frame")) {
        siteDF <- as.data.frame(siteDF, stringsAsFactors = FALSE)
    }
    if (!inherits(specDF, "data.frame")) {
        specDF <- as.data.frame(specDF, stringsAsFactors = FALSE)
    }
    DF = NULL
    DF <- switch(EXPR = type, sites = siteDF, species = specDF,
        {
            specDF$id.type <- "Taxa"
            siteDF$id.type <- "Samples"
            colnames(specDF)[1:2] <- colnames(siteDF)[1:2]
            DF = merge(specDF, siteDF, all = TRUE)
            if (!is.null(shape)) {
                DF <- rp.joint.fill(DF, shape, "Samples")
            }
            if (!is.null(shape)) {
                DF <- rp.joint.fill(DF, shape, "Taxa")
            }
            if (!is.null(color)) {
                DF <- rp.joint.fill(DF, color, "Samples")
            }
            if (!is.null(color)) {
                DF <- rp.joint.fill(DF, color, "Taxa")
            }
            DF
        })
    if (justDF) {
        return(DF)
    }
    if (!is.null(color)) {
        if (!color %in% names(DF)) {
            warning("Color variable was not found in the available data you provided.",
                "No color mapped.")
            color <- NULL
        }
    }



    if (!is.null(shape)) {
        if (!shape %in% names(DF)) {
            warning("Shape variable was not found in the available data you provided.",
                "No shape mapped.")
            shape <- NULL
        }
    }
    if (!is.null(label)) {
        if (!label %in% names(DF)) {
            warning("Label variable was not found in the available data you provided.",
                "No label mapped.")
            label <- NULL
        }
    }
    x = colnames(DF)[1]
    y = colnames(DF)[2]
    if (ncol(DF) <= 2) {
        message("No available covariate data to map on the points for this plot `type`")
        ord_map = aes_string(x = x, y = y)
    }
    else if (type %in% c("sites", "species", "split")) {
        ord_map = aes_string(x = x, y = y, color = color, shape = shape,
            na.rm = TRUE)
    }
    else if (type == "biplot") {
        if (is.null(color)) {
            ord_map = aes_string(x = x, y = y, size = "id.type",
                color = "id.type", shape = shape, na.rm = TRUE)
        }
        else {
            ord_map = aes_string(x = x, y = y, size = "id.type",
                color = color, shape = shape, na.rm = TRUE)
        }
    }
# return(ord_map)
    p <- ggplot(DF, ord_map) + geom_point(na.rm = TRUE)
	p <- p + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)

	if(continuous) {
		p <- p + scale_colour_gradient(low=colourScale[1], high=colourScale[2])
	} else {
		if(cbPalette) {
			p<-p+scale_colour_manual(values=cbbPalette)#	+ guides(colour=guide_legend(title=design))
		} else {
			p<-p+scale_colour_viridis(discrete=TRUE) # + guides(colour=guide_legend(title=design))
		}
	}

    if (type == "split") {
        p <- p + facet_wrap(~id.type, nrow = 1)
    }
    if (type == "biplot") {
        if (is.null(color)) {
            p <- update_labels(p, list(colour = "Ordination Type"))
        }
        p <- p + scale_size_manual("type", values = c(Samples = 5,
            Taxa = 2))
    }
    if (!is.null(label)) {
        label_map <- aes_string(x = x, y = y, label = label,
            na.rm = TRUE)
        p = p + geom_text(label_map, data = rm.na.phyloseq(DF,
            label), size = 2, vjust = 1.5, na.rm = TRUE)
    }
    if (!is.null(title)) {
        p = p + ggtitle(title)
    }
    if (length(extract_eigenvalue(ordination)[axes]) > 0 & !rescale) {
        strivar = as(c(p$label$x, p$label$y), "character")
        strivar = paste0(strivar, "   [", percvar, "%]")
        p = p + xlab(strivar[1]) + ylab(strivar[2])
    }
    return(p)
}

prcomp2 <-
function(obj) {
	X <- prcomp(t(obj))
	X$percentVar<-X$sdev^2/sum(X$sdev^2)
	X
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}# subsamples OTU data
# set depth to rarefy data, or use some a vector for a proportional scheme (e.g., if you want to combine different numbers of biological replicate per condition)
# sampling can be done with or without replacement

subsample <- function(
  obj,
  depth=1000,
  replace=T,
  sdr=.Random.seed
){
  if(length(depth)==1) {
    depth <- sapply(1:ncol(obj),function(i) min(depth,sum(obj[,i])))
  }

  l <- lapply(1:ncol(obj),function(i) {
    set.seed(sdr)
    if(replace) {
      otu <- sample(rownames(obj),depth[i],replace=T,prob=(obj[,i]/sum(obj[,i])))
	  } else {
      otu <- sample(unlist(sapply(seq_along(obj[,i]),function(x) {rep(rownames(obj)[x],obj[x,i])})))
      otu <- sample(otu,depth[i])
    }
    table(otu)
  })
  cts <- Reduce(function(...) merge(..., all=T,by = "otu"), l)
  rownames(cts) <- cts[,1]
  cts <- cts[-1]
  cts[is.na(cts)] <- 0
  colnames(cts) <- colnames(obj)
  cts
}

sumTaxa <- function(
	obj,
	taxon="phylum",
	design="all",
	proportional=F,
	conf=0.65,
	meanDiff=F,
	weighted=T
){
# sums by sample data 
	suppressPackageStartupMessages(require(plyr))
	suppressPackageStartupMessages(require(dplyr))
	suppressPackageStartupMessages(require(tibble))
	suppressPackageStartupMessages(require(reshape2))

	obj[[2]][,taxon] <- taxaConfVec(obj[[2]][,-8],conf=conf,level=which(colnames(obj[[2]])==taxon))
	dtx <- left_join(rownames_to_column(obj[[1]]),rownames_to_column(obj[[2]][,taxon,drop=F]))
	rownames(dtx) <- dtx[,1]
   	md <- melt(dtx,id=c(taxon,"rowname"))

	obj[[3]]$all <- "all"
	if(length(design)>1) {
		obj[[3]][[paste(design,collapse=":")]] <- factor(apply(as.data.frame(obj[[3]][, design,drop = FALSE]), 1, paste, collapse = " : "))
		design <- paste(design,collapse=":")		
	}

	if(meanDiff&length(levels(as.factor(obj[[3]][[design[1]]])))!=2) {
		warning("Design has more than 2 levels mean difference may be difficult to interpret")	
	}

	if (!meanDiff) {
		md <- md[,setdiff(names(md), "rowname")]
	}

	md$variable <- mapvalues(md$variable,from=rownames(obj[[3]]), to=as.character(obj[[3]][,design]))
	nd <- dcast(md,...~variable,sum)

	if(meanDiff) {
		temp <- nd[,c(-1,-2)]
		cols <- ncol(temp)
		nnd <- as.data.frame(lapply(seq(1,(cols-1)),function(i) lapply(seq((i+1),cols),function(j) 
		{x=as.data.frame(abs(temp[,i]-temp[,j]));
		 colnames(x)[1] <-paste(colnames(temp)[i],colnames(temp)[j],sep="|");
		 if(weighted){x<-x*((temp[,i]+temp[,j])/sum(temp[,c(i,j)]))}
		 return(x)}))) 
		nd <- aggregate(nnd,by=list(nd[[taxon]]),sum)
	}

	colnames(nd)[1] <- taxon
	if(proportional) {
		nd[-1] <-  prop.table(as.matrix(nd[,-1]),2)*100
	}
	return(nd)
}

sumTaxaAdvanced <-
function (
	obj,# list(countdata,taxData,colData)
	taxon="phylum", 	# taxon (str) is the taxonomic level of interest
	conf=0.65,
	design="all", 		# condition (str) describes how the samples should be grouped (must be column of sample data)
	proportional=T,	# proportional (bool) whether the graph should use proportional or absolute values
	cutoff=1, 	# cutoff (double) for proportional graphs. 
	topn=0, 		# topn (int)taxons to display (by total reads) for non-prortional graphs. T
	others=T, 	# combine values less than cutoff/topn into group "other"
	ordered=F 	# order by value (max to min)
) {
	
	obj$taxData <- 
	taxa_sum <- sumTaxa(obj,taxon=taxon,design=design,conf=conf)
	taxa_sum$taxon[grep("\\(",taxa_sum$taxon)] <- taxa_sum$taxon[sub("\\(.*"," incertae sedis",taxa_sum$taxon)]

	if(!topn) {
		obj[[3]]$MLUflop <- 1 #assigns the MLU flop digit
		tx <- sumTaxa(obj,taxon=taxon,"MLUflop")
		tx[,-1] <- prop.table(as.matrix(tx[,-1]),2)*100
		txk <- tx[tx[,2]>=cutoff,1]
	} else {
		taxa_sum[,ncol(taxa_sum)+1]<- 0
		taxa_sum <- taxa_sum[order(rowSums(taxa_sum[,-1]),decreasing=T),]
		taxa_sum <- taxa_sum[,-ncol(taxa_sum)]	
		txk <- taxa_sum[1:topn,1]
	}
	
	if(proportional) {
		taxa_sum[,-1] <- prop.table(as.matrix(taxa_sum[,-1]),2)*100
	}

	taxa_cut <- taxa_sum[taxa_sum[,1]%in%txk,]
	taxa_cut <- taxa_cut[order(taxa_cut[,1],decreasing=F),]
	if(others) {
		taxa_cut <- rbind(taxa_cut,setNames(data.frame(x="others" ,t(colSums(as.data.frame(taxa_sum[!taxa_sum[,1]%in%txk,-1])))),names(taxa_cut)))
	}
	taxa_cut <- na.omit(taxa_cut)
	taxa_cut[,1] <- as.factor(taxa_cut[,1])
	if(ordered) {
		taxa_cut <- taxa_cut[order(rowSums(as.data.frame(taxa_cut[,-1])),decreasing=T),]
	}
	
	return(taxa_cut)

}

taxaConf <- 
  function (obj, conf = 0.65, level = 7) 
  {
    rank <- apply(obj, 1, function(x) {
      l <- (level + 7)
      y <- abs(as.numeric(x[8:l]))
      i <- last(which(y >= conf))
      i[is.na(i)] <- 1 # Change as the which function returned "integer(0) not na.
      #i[isTRUE(identical(i,integer(0)))] <- 1 #changed back - seems to be working now with updated Rlang
			# Handle the case where i is empty
      if (length(i) == 0) {
        i <- 1
      }
      if (i == 1) {
        s = "(k)"
      }
      else if (i == 2) {
        s = "(p)"
      }
      else if (i == 3) {
        s = "(c)"
      }
      else if (i == 4) {
        s = "(o)"
      }
      else if (i == 5) {
        s = "(f)"
      }
      else if (i == 6) {
        s = "(g)"
      }
      else {
        s = "(s)"
      }
      ret <- paste(x[i], s, sep = "")
      return(ret)
    })
    X <- suppressMessages(as.data.frame(as.matrix(obj)))
    X$rank <- as.character(rank)
    return(X)
  }

taxaConfVec <-
function (obj,conf=0.65,level=7){
	rank <- apply(obj,1,function(x){
		l <- (level+7)
		y <-abs(as.numeric(x[8:l]))
		i <- last(which(y>=conf))
		i[is.na(i)] <-1 
		##edit
  	    s=""
	    if(i!=level) {
		 if(i==1){s="(k)"
		 } else if(i==2) {s="(p)"
		 } else if(i==3) {s="(c)"
		 } else if(i==4) {s="(o)"	
		 } else if(i==5) {s="(f)"
		 } else if(i==6) {s="(g)"
		 } else {s="(s)"}
        }
		ret <- paste(x[i],s,sep="")
 		return(ret)
	})
    return(as.character(rank))
}

taxonomyTidy <- function(x) {
	if (x[2]=="unknown") {x[2] <- paste(x[1],"(k)",sep="")}
	if (x[3]=="unknown") {if(any(grep('\\(',x[2]))) {x[3]<-x[2]}else{x[3]<-paste(x[2],"(p)",sep="")}}
	if (x[4]=="unknown") {if(any(grep('\\(',x[3]))) {x[4]<-x[3]}else{x[4]<-paste(x[3],"(c)",sep="")}}
	if (x[5]=="unknown") {if(any(grep('\\(',x[4]))) {x[5]<-x[4]}else{x[5]<-paste(x[4],"(o)",sep="")}}
	if (x[6]=="unknown") {if(any(grep('\\(',x[5]))) {x[6]<-x[5]}else{x[6]<-paste(x[5],"(f)",sep="")}}
	if (x[7]=="unknown") {if(any(grep('\\(',x[6]))) {x[7]<-x[6]}else{x[7]<-paste(x[6],"(g)",sep="")}}
	return(x)
}

taxa_to_funguild <- function(taxData,confidence=.6) {
  taxData <- phyloTaxaTidy(taxData,confidence)
  taxData$V <- gsub("(.*\\()([a-zA-Z])(\\))","\\2",taxData$rank)
  col_keep <- c(names(taxData)[1:7],"V")
  as.data.frame(t(apply(taxData[,col_keep],1,function(x) {
    x[-1:-match(x[8],substr(col_keep[1:7],1,1))] <- "unidentified"
    x[-8]
  }))) 
}


ubiom_to_des <- function(
	obj, 
	design=~1,
	fit=F,
	filter,
	calcFactors=function(d)
	{
		sizeFactors(estimateSizeFactors(d))
	},
	...
){
	suppressPackageStartupMessages(require(DESeq2))

	invisible(mapply(assign, names(obj), obj,MoreArgs=list(envir = environment())))
	
	colData <- colData[colnames(countData),,drop = FALSE]

	if(!missing(filter)) {
		filter <- eval(filter)
		colData <- droplevels(colData[filter,])
		countData <- countData[,filter]
	}

	dds <- 	suppressWarnings(DESeqDataSetFromMatrix(countData,colData,design))

	sizeFactors(dds) <- calcFactors(dds)

    if (fit) {
     	return(DESeq(dds,...))
    } else {
    	return(dds)
    }
} 
ubiom_to_phylo <- function(
	obj
){
	suppressPackageStartupMessages(require(phyloseq))
	phyloseq(
		otu_table(obj[[1]],taxa_are_rows=T),
	 	tax_table(as.matrix(obj[[3]])),
	 	sample_data(obj[[2]])
	 )
}
ubiome_to_taxmap <- function(
	obj#,
	#filter="rowSums(otu_table[,c(-1,-2)])>=0"
){
	suppressPackageStartupMessages(require(metacoder))
	# tibblise data
	otu_table <- as.tibble(obj[[1]],rownames="otu_id")
	sample_data  <- as.tibble(obj[[2]],rownames="sample_id")
	tax_data <- as.tibble(obj[[3]],rownames="otu_id")
	parsed_tax <- apply(tax_data,1,function(x) {
		x <- sub(".*\\(.*",NA,x[2:8])
		x <- x[!is.na(x)]
		x <- sub("_SH.*","",x)
		x <- gsub("_"," ",x)
		xx <- hierarchy(x)
		lapply(seq_along(xx$taxa),function(i) {
			xx$taxa[[i]]$rank$name <<- names(x)[i]
		})
		xx
	})
	output <- taxmap(.list = parsed_tax, named_by_rank = T)
	# set the taxon_id to the rank (rank is the lowest defined rank with a given confidence)
	t1 <- output$taxon_names()
	t2 <- sub("\\(.*","",tax_data$rank)
	t2 <- sub("_SH.*","",t2)
	t2 <- gsub("_"," ",t2)	
	t3 <- sapply(t2,function(x) names(t1[t1==x])[1])
	otu_table <- as.tibble(cbind(taxon_id=t3,otu_table,stringsAsFactors=F))
	#otu_counts <- otu_table[eval(filter),]
	output$data <- list(
		otu_table   = otu_table#,
		#otu_counts  = otu_counts,
		#sample_data = sample_data
	)
	output	
}

calc_taxon_abund <-
function (obj, data, cols = NULL, groups = NULL, out_names = NULL)
{
  do_it <- function(count_table, cols = cols, groups = groups) {
        my_print("Summing per-taxon counts from ", length(cols),
            " columns ", ifelse(length(unique(groups)) == length(unique(cols)),
                "", paste0("in ", length(unique(groups)), " groups ")),
            "for ", length(obj$taxon_ids()), " taxa")
        obs_indexes <- obj$obs(data)
        output <- lapply(split(cols, groups), function(col_index) {
            col_subset <- count_table[, col_index]
            vapply(obs_indexes, function(i) {
            	if(length(i)>0) {
                	sum(col_subset[i, ])
                } else numeric(1)
            }, numeric(1))
        })
        output <- as.data.frame(output, stringsAsFactors = FALSE)
        return(output)
    }
    output <- do_calc_on_num_cols(obj, data, cols = cols, groups = groups,
        other_cols = NULL, out_names = out_names, func = do_it)
    output <- cbind(data.frame(taxon_id = obj$taxon_ids(), stringsAsFactors = FALSE),
        output)
    dplyr::as.tbl(output)
}		     
