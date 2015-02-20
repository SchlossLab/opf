read.lt <- function(d){

	n <- as.numeric(d[1]); d <- d[-1]
	m <- matrix(rep(NA, n*n), nrow=n)
	
	groups <- d[1]; d <- d[-1]

	for(i in 2:n){
		groups <- c(groups, d[1]); d <- d[-1]
		m[i,1:(i-1)] <- as.numeric(d[1:(i-1)])
		d <- d[-c(1:(i-1))]
	}
	
	rownames(m) <- groups
	colnames(m) <- groups
	
	return(m)
}

read.trimatrix <- function(infile) {
    #file is distance matrix from mothur ex. preg.eggnog.thetayc.1.lt.dist
    matrixcols <- readLines(infile, n=1)
    d <- read.table(file=infile, fill=TRUE, header=F, skip=1, row.names=1, col.names=paste("V", 1:matrixcols))
    d$x<-NA
    mat<-data.matrix(d)
    mat[is.na(mat)]<-0
    mat <- mat + t(mat)
    colnames(mat)<-rownames(mat)
    return(mat)
}



subject <- c("101.1"="101", "101.3"="101", "160.1"="160", "160.3"="160", "195.1"="195", "195.3"="195",
			 "226.1"="226", "226.3"="226", "234.1"="234", "234.3"="234", "237.1"="237", "237.3"="237", 
			 "264.1"="264", "264.3"="264", "303.1"="303", "303.3"="303", "312.1"="312", "312.3"="312",
			 "342.1"="342", "342.3"="342")


a <- read.table(file="preg.16s.0.03.subsample.tshared", header=T, skip=1)
rownames(a) <- a$Group
a <- a[-1,-1]

nrow(a)
#[1] 1176

sum(apply(a>0, 1, sum)==ncol(a))
#[1] 2
sum(apply(a>0, 1, sum)==ncol(a))/nrow(a)
#[1] 0.00170068

sum(apply(a>0, 1, sum)==1)
#[1] 823
sum(apply(a>0, 1, sum)==1)/nrow(a)
#[1] 0.6998299

d <- scan("preg.16s.thetayc.0.03.lt.ave.dist", what="")
lt <- read.lt(d)
mean(lt, na.rm=T)
#[1] 0.8566066
median(lt, na.rm=T)
#[1] 0.8866105
quantile(lt, na.rm=T, probs=c(0.025, 0.975))
#     2.5%     97.5% 
#0.5685698 0.9734358 





a <- read.table(file="preg.kegg.1.subsample.tshared", header=T, skip=1)
rownames(a) <- a$Group
a <- a[-1,-1]

nrow(a)
#[1] 4002

sum(apply(a>0, 1, sum)==ncol(a))
#[1] 2036
sum(apply(a>0, 1, sum)==ncol(a))/nrow(a)
#[1] 0.5087456

sum(apply(a>0, 1, sum)==1)
#[1] 70 		#probably not meaningful
sum(apply(a>0, 1, sum)==1)/nrow(a)
#[1] 0.01749125 #probably not meaningful


a.subset <- a[apply(a, 1, sum) > (ncol(a)*20),]
a.subset <- a.subset + runif(nrow(a.subset)*ncol(a.subset), 0,0.01)

indiv <- levels(factor(subject))
w.pvalues <- rep(0, nrow(a.subset))
estimates <- rep(0, nrow(a.subset))
t.pvalues <- rep(0, nrow(a.subset))

for(i in 1:nrow(a.subset)){
	test <- wilcox.test(x=as.numeric(a.subset[i,paste("X", indiv, ".1", sep="")]), y=as.numeric(a.subset[i,paste("X", indiv, ".3", sep="")]), paired=T, conf.int=T)
	w.pvalues[i] <- test$p.value
	estimates[i] <- test$estimate

	test <- t.test(x=as.numeric(a.subset[i,paste("X", indiv, ".1", sep="")]), y=as.numeric(a.subset[i,paste("X", indiv, ".3", sep="")]), paired=T)#, conf.int=T)
	t.pvalues[i] <- test$p.value
}

pdf("kegg.pdf")
plot(estimates, w.pvalues, ylim=rev(c(0.001, 1)), log="y", xlim=c(-500, 500))
dev.off()


length(w.pvalues[w.pvalues<0.01])
min(w.pvalues)
min.p <- which.min(w.pvalues)
min(p.adjust(w.pvalues, method="BH"))
t1 <- as.numeric(a.subset[min.p,paste("X", indiv, ".1", sep="")])
t3 <- as.numeric(a.subset[min.p,paste("X", indiv, ".3", sep="")])
diff <- t3-t1
median(diff)
#[1] -13
rownames(a[min.p,])
#[1] "K03644"  #lipoic acid synthetase



d <- scan("preg.kegg.thetayc.1.lt.ave.dist", what="")
lt <- read.lt(d)
mean(lt, na.rm=T)
#[1] 0.1199133
median(lt, na.rm=T)
#[1] 0.10388
quantile(lt, na.rm=T, probs=c(0.025, 0.975))
#     2.5%     97.5% 
#0.0412175 0.2773472 





a <- read.table(file="preg.opf.0.30.subsample.tshared", header=T, skip=1)
rownames(a) <- a$Group
a <- a[-1,-1]

nrow(a)
#[1] 659740

sum(apply(a>0, 1, sum)==ncol(a))
#[1] 7239
sum(apply(a>0, 1, sum)==ncol(a))/nrow(a)
#[1] 0.0109725

sum(apply(a>0, 1, sum)==1)
#[1] 77316
sum(apply(a>0, 1, sum)==1)/nrow(a)
#[1] 0.1171916

a.subset <- a[apply(a, 1, sum) > (ncol(a)*20),]
a.subset <- a.subset + runif(nrow(a.subset)*ncol(a.subset), 0,0.01)


indiv <- levels(factor(subject))
w.pvalues <- rep(0, nrow(a.subset))
estimates <- rep(0, nrow(a.subset))
t.pvalues <- rep(0, nrow(a.subset))

for(i in 1:nrow(a.subset)){
	test <- wilcox.test(x=as.numeric(a.subset[i,paste("X", indiv, ".1", sep="")]), y=as.numeric(a.subset[i,paste("X", indiv, ".3", sep="")]), paired=T, conf.int=T)
	w.pvalues[i] <- test$p.value
	estimates[i] <- test$estimate

	test <- t.test(x=as.numeric(a.subset[i,paste("X", indiv, ".1", sep="")]), y=as.numeric(a.subset[i,paste("X", indiv, ".3", sep="")]), paired=T)#, conf.int=T)
	t.pvalues[i] <- test$p.value
}

pdf("opf.pdf")
plot(estimates, w.pvalues, ylim=rev(range(c(0.001, 1))), log="y", xlim=c(-500, 500))
dev.off()

min(w.pvalues)
min.p <- which.min(w.pvalues)
min(p.adjust(w.pvalues, method="BH"))
t1 <- as.numeric(a.subset[min.p,paste("X", indiv, ".1", sep="")])
t3 <- as.numeric(a.subset[min.p,paste("X", indiv, ".3", sep="")])
diff <- t3-t1
median(diff)
rownames(a[min.p,])


d <- scan("preg.opf.thetayc.0.30.lt.ave.dist", what="")
lt <- read.lt(d)
mean(lt, na.rm=T)
#[1] 0.5152775
median(lt, na.rm=T)
#[1] 0.489975
quantile(lt, na.rm=T, probs=c(0.025, 0.975))
#     2.5%     97.5% 
#0.2331931 0.7871524 





kegg <- read.table(file="preg.kegg.tshared", header=T, skip=1)
rownames(kegg) <- kegg$Group
kegg <- kegg[-1,-1]


opf <- read.table(file="preg.opf.tshared", header=T, skip=1)
rownames(opf) <- opf$Group
opf <- opf[-1,-1]


kegg.count <- apply(kegg, 2, sum)
opf.count <- apply(opf, 2, sum)
fraction <- 1-kegg.count/opf.count



#Figure 2

d <- scan("preg.16s.thetayc.0.03.lt.ave.dist", what="")
otu.lt <- read.lt(d)

d <- scan("preg.kegg.thetayc.1.lt.ave.dist", what="")
kegg.lt <- read.lt(d)

d <- scan("preg.opf.thetayc.0.30.lt.ave.dist", what="")
opf.lt <- read.lt(d)

samples <- c("101.1"="101", "101.3"="101", "160.1"="160", "160.3"="160", "195.1"="195", "195.3"="195",
			 "226.1"="226", "226.3"="226", "234.1"="234", "234.3"="234", "237.1"="237", "237.3"="237", 
			 "264.1"="264", "264.3"="264", "303.1"="303", "303.3"="303", "312.1"="312", "312.3"="312",
			 "342.1"="342", "342.3"="342")
subjects <- levels(factor(samples))

otu <- vector()
opf <- vector()
kegg <- vector()

for(i in subjects){
	t1 <- paste(i, "1", sep=".")
	t3 <- paste(i, "3", sep=".")
	
	otu[i] <- otu.lt[t3, t1]
	opf[i] <- opf.lt[t3, t1]	
	kegg[i] <- kegg.lt[t3, t1]
}

pdf("Figure2.pdf")
par(mar=c(5, 5, 1, 1))
plot(1, type="n", xlim=c(0,3), ylim=c(0,1), xaxt="n", xlab="Method of clustering sequences", ylab=expression(paste(Theta[YC], " distance between first and third trimester", sep="")))
points(x=rep(0.5, length(kegg)), y=kegg, pch=19)
points(x=rep(1.5, length(otu)), y=otu, pch=19)
points(x=rep(2.5, length(opf)), y=opf, pch=19)
segments(x0=rep(0.5, length(kegg)), x1=rep(1.5, length(otu)), y0=kegg, y1=otu)
segments(x0=rep(2.5, length(opf)), x1=rep(1.5, length(otu)), y0=opf, y1=otu)
axis(1, at=c(0.5, 1.5, 2.5), label=c("KEGG", "OTU", "OPF"))
dev.off()
