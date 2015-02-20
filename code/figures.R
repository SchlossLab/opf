library(RSvgDevice)
library(ggplot2)
library(gplots)
library(reshape)


# Figure 1

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                           c( N    = length2(xx[,col], na.rm=na.rm),
                              mean = mean   (xx[,col], na.rm=na.rm),
                              sd   = sd     (xx[,col], na.rm=na.rm)
                              )
                          },
                    measurevar,
                    na.rm
             )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


t.dists <- read.csv(file='distance_comparison.csv', sep=',', header=T)
dfc <- summarySE(t.dists, measurevar="dist_bars", groupvars=c("type","method"))

devSVG(file='Fig1.svg')
ggplot(dfc, aes(x=method, y=dist_bars, group=type, fill=type)) + 
    geom_errorbar(aes(ymin=dist_bars-se, ymax=dist_bars+se), width=.1, position=position_dodge(.9)) +
    geom_bar(position=position_dodge(), stat='identity', colour='black') + theme_bw() + ylab('Distance') + xlab(' ') + scale_x_discrete(limits=c("unifrac", "unifrac_uw", "16s","opf", "kegg"), labels=c('Weighted Unifrac', 'Unweighted Unifrac', '16S', 'OPF', 'KEGG')) + theme(legend.justification=c(0,1), legend.position=c(0,1)) + scale_fill_manual(values=c("#000000", '#666666', "white"), 
                       name="Comparison",
                       breaks=c("all", "family", "twins"),
                       labels=c("Unrelated", "Inside family", "Twin to twin"))
dev.off()


# Figure 2

#t.pcoa <- read.table(file='twin_pcoa_opf_kegg_16s_vars.csv', sep=',', header=T)

#devSVG(file='Fig2.svg')
#ggplot(t.pcoa, aes(axis1, axis2, color=weight, shape=as.factor(group), group=as.factor(family))) + geom_point(size=4) + theme_bw() + xlab("Dimension 1") + ylab("Dimension 2") + geom_polygon(aes(mapping=group, alpha=1)) + facet_wrap(~method, scales='free', ncol=1) + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(color = 'black'))
#dev.off()


# Figure 3

t.bin <- read.table(file='twin_opf_bin.csv', sep=' ', header=T, row.names=1)
rsum <- rowSums(t.bin)
t.bin['sum'] <- NA
t.bin$sum <- rsum
twin.o.sums <- matrix(nrow=18, ncol=1)
 
for (i in 1:18){
	s <- sum(t.bin$sum==i)
	twin.o.sums[i,1] <- s
}
tk.bin <- read.table(file='twin_kegg_bin.csv', sep=' ', header=T, row.names=1)
rksum <- rowSums(tk.bin)
tk.bin['sum'] <- NA
tk.bin$sum <- rksum
twin.k.sums <- matrix(nrow=18, ncol=1)
 
for (i in 1:18){
	s <- sum(tk.bin$sum==i)
 	twin.k.sums[i,1] <- s
}
ts.bin <- read.table(file='twin_16s_bin.csv', sep=' ', header=T, row.names=1)
rssum <- rowSums(ts.bin)
ts.bin['sum'] <- NA
ts.bin$sum <- rssum
twin.s.sums <- matrix(nrow=18, ncol=1)
 
for (i in 1:18){
 	s <- sum(ts.bin$sum==i)
 	twin.s.sums[i,1] <- s
}

devSVG(file='Fig3.svg')
plot(twin.o.sums, log='y', ylim=c(0.6, 6000), xlab='# of samples', ylab='Shared OPFs/OTUs/KEGG Categories', pch=15, cex=4)
points(twin.s.sums, pch=16, cex=4)
points(twin.k.sums, pch=17, cex=4)
legend(1, 100, c('OPF', 'OTU', 'KEGG'), pch=c(15, 16, 17), , cex=4)
dev.off()


# Figure 4

preg.pcoa <- read.csv(file='pregnancy_pcoa_vars.csv', sep=',', header=T)

devSVG(file='Fig4.svg')
ggplot(preg.pcoa, aes(axis1, axis2, shape=trimester)) + geom_point(size=5) + theme_bw()  + facet_wrap(~method, ncol=1, scales='free')
dev.off()

# Figure 5

#m.pcoa <- read.csv(file='mice/mice_pcoa_opf_kegg_16s_vars.txt', sep='\t', header=T)

#devSVG(file='Fig5.svg')
#ggplot(m.pcoa, aes(axis1, axis2, color=time)) + geom_point(size=4) + theme_bw() + xlab("Dimension 1") + ylab("Dimension 2") + facet_wrap(~method, scales='free', ncol=1) + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(color = 'black'))
#dev.off()


# Figure 6

#cor.mice <- read.table(file='mice/mice_correlation_table_for_heatmap.csv', sep=',', header=T)
#cor.mice <- as.matrix(cor.mice[,2:34])
#colors <- c("#003399", "#0000CC", "#0033CC", "#0000FF", "#3333FF", "#3366FF", "#0099FF", "#66CCFF", "#FFFFFF", "#FF9999", "#FF6666", "#FF3333", "#CC3333", "#993333", "#990000", "#660000")

#heatmap.2(cor.mice, dendrogram="both", col=colors, trace="none", margins=c(10, 0), keysize=1, density.info="none", cexCol=0.8, cexRow=0.6)
