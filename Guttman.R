# Instructions: make sure that this script and input files are in the
# same directory and R's directory is set to that.
# Read in the four functions wrap(), Gm(), impute(), and p() by selecting
# all the code below the block of four lines using read.table().
# For the last function to work the package vegan needs to be installed.
# Select the dataset you want to analyze by running one of the four lines below.
# Or you can read in a file of your own having a similar format.
# Writing Gm(g, outfile="yes") gives you a Guttman coefficient (GC) and a scale in 
# a file called GSout.txt. If you are not interested in the scale simply write Gm(g).
# Writing wrap(g,10000) gives you an imputed matrix created by N=10000 attempts
# to replace NA's in a way such that the GC of the imputed matrix is as close 
# as possible to that of the original matrix, and a p-value based on 9999 permuations 
# of the imputed matrix output to the console. 

g <- read.table(file="1-2.txt",header=TRUE)
g <- read.table(file="1-3.txt",header=TRUE)
g <- read.table(file="1-3LOC.txt",header=TRUE)
g <- read.table(file="3-1.txt",header=TRUE)

wrap <- function(g,N) {
	ximp <- impute(g,N)
	pval <- p(ximp)
	return(pval)
}

Gm <- function(x, outfile="no") {
# transpose the matrix so verbs are columns
# reorder the matrix gathering 1's in the top left
# find optimal cut-off for minimizing errors of inclusion and exclusion
	x <- t(x)
	rows <- length(x[,1]); cols <- length(x[1,])
	horsum <- function(y) { sum(x[y,],na.rm=T) / (cols - length(which(is.na(x[y,1:cols])))) }
	versum <- function(z) { sum(x[,z],na.rm=T) / (rows - length(which(is.na(x[1:rows,z])))) }
	hsum <- c(); vsum <- c()
	for (i in 1:rows) {hsum[i] <- horsum(i)}
	for (i in 1:cols) {vsum[i] <- versum(i)}
	hs <- rev(order(hsum)); vs <- rev(order(vsum))
	left <- x[,vs] # puts the columns with most 1's to left
	x <- left[hs,] # puts the rows with most 1's on top
# determine the ranking of each verb by greedily getting a maximal number of 1's
# with a minimal number of errors
# vp is the number of legitimate ones scored by a given individual
	err <- c(); vp <- c()
	for (i in 1:rows) {
		err.zero <- c(); err.one <- c()
		if ( length(grep(0, x[i,]))==0 | length(grep(1, x[i,]))==0 ) {
			err[i] <- 0
			vp[i] <- length(grep(1, x[i,]))
		} else {
			for (j in 1:(cols-1)) {
				err.zero[j] <- length(grep(0,x[i,1:j]))
				err.one[j] <- length(grep(1,x[i,(j+1):cols]))
			}
		err.zo <- err.zero + err.one
		w <- max(which(err.zo==min(err.zo)))
		err[i] <- err.zo[w]
# the following will give the rank based on the number of 1's
		vp[i] <- length(grep(1, x[i,1:w]))
# the following will give the rank based on the cut-off for 1's to the left
# so NA's to the left of the cutoff are counted as ones
#		vp[i] <- w
		}
	}
	GC <- 100-round(100*sum(err)/(length(x)-length(which(is.na(x)))),2)
	if (outfile=="yes") {
		GS <- rownames(x); alt <- colnames(x)[1]
		ranking <- c(max(rank(vp,ties.method="min"))+1-rank(vp,ties.method="min"))
		write.table(cbind(GS,ranking),file="GSout.txt",quote=F,sep="\t",row.names=F)
	}
	return(GC)
}

# Creates a matrix where NA's are replaced by imputed 1'and 0's
impute <- function(x, N) {
	nas <- sum(apply(x, 2, function(y) length(which(is.na(y)))))
	pn <- matrix(0, nrow=nas, ncol=2)  # stands for position of nas
	GC <- Gm(x)
	xcol <- x
	counter <- 0
	for (i in 1:length(xcol[,1])) {
		for (j in 1:length(xcol[1,])) {
			if ( is.na(xcol[i,j]) ) {
				counter <- counter + 1
				pn[counter,1] <- i; pn[counter,2] <- j
				xcol[i,j] <- 1
				GC1 <- Gm(xcol)
				xcol[i,j] <- 0
				GC0 <- Gm(xcol)
				ifelse ( abs(GC - GC1) <= abs(GC - GC0), xcol[i,j] <- 1, xcol[i,j] <- 0 )
			}
		}
	}
	GCcol <- Gm(xcol)
	xrow <- x
	for (k in 1:length(xrow[1,])) {
		for (l in 1:length(xrow[,1])) {
			if ( is.na(xrow[l,k]) ) {
				xrow[l,k] <- 1
				GC1 <- Gm(xrow)
				xrow[l,k] <- 0
				GC0 <- Gm(xrow)
				ifelse ( abs(GC - GC1) <= abs(GC - GC0), xrow[l,k] <- 1, xrow[l,k] <- 0 )
			}
		}
	}
	GCrow <- Gm(xrow)
	ifelse ( abs(GC - GCcol) <= abs(GC - GCrow), ximp <- xcol, ximp <- xrow )
	cat("GC:",GC,"\n")
	cat("GCrow first iteration:",GCrow,"\n")
	cat("GCcol first iteration:",GCcol,"\n")
# at this point ximp, the imputed matrix should be have a pretty close match of GC
# now the positions at which there were NAs are randomly accessed
# and changes are made to 0 and 1 until the GC is as the original one
	GCimp <- Gm(ximp)
	GCnew <- 0
	rounds <- 0
	for (i in 1:N) { 
		rounds <- rounds + 1
		GCimp <- Gm(ximp)
		pos <- sample(nas, 1)
		old.state <- ximp[pn[pos,1],pn[pos,2]]
		if ( old.state==0 ) {
			new.state <- 1
		}
		if ( old.state==1 ) {
			new.state <- 0
		}
		ximp[pn[pos,1],pn[pos,2]] <- new.state
		GCnew <- Gm(ximp)
		if ( abs(GC - GCnew) < abs(GC - GCimp) ) {
			ximp[pn[pos,1],pn[pos,2]] <- new.state
		}
		if ( abs(GC - GCnew) >= abs(GC - GCimp) ) {
			ximp[pn[pos,1],pn[pos,2]] <- old.state
		}
		cat("improvement attempt",rounds,"--currently this close to matching GC:", abs(GC-Gm(ximp)),"\n")
		if (GCnew == GC) {
			return(ximp)
		}
		if (rounds == N) {
			cat("Stopping, tried", N," to get a GC equal to the original without luck\n")
			cat("I am returning the best I could do\n")
			return(ximp)
			break
		}
	}	
}

#p-values for imputed matrix
require(vegan)
p <- function(ximp) {
	gc_imp <- Gm(ximp)
	cat("GC of the imputed matrix:", gc_imp, "\n")
	mat <- permatfull(ximp, mtype="prab", times=9999)
	mat$perm[[10000]] <- mat$orig
	GCs <- rep(0,10000)
	for (j in 1:10000) {
		GCs[j] <- Gm(mat$perm[[j]])
	}
	pval <- (length(GCs[GCs >= gc_imp]))/length(GCs)
	return(pval)
}



