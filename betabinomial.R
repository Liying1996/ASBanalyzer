library(VGAM)

##### Read Parameters #####

file_name <- commandArgs(trailingOnly=TRUE)[1]
title <- commandArgs(trailingOnly=TRUE)[2]
fdr <- as.numeric(commandArgs(trailingOnly=TRUE)[3])


##### All functions ##### 

# Process Data
process_data <- function(file_name, minN)
{
  sample <-  read.csv(file_name, sep = '\t', header = F)
  sample <- sample[, 1:6]
  sample <- sample[sample[,5]+sample[,6] >= minN, ] 
  sample$lower <- apply(sample[,5:6],1,min)
  sample$total <- sample[,5] + sample[,6]
  colnames(sample) <- c("Chr","Pos","Ref","Alt","Ref_counts","Alt_counts","lower","total")
  sample$allelicRatio <- sample$Ref_counts/sample$total
  return(sample)
}


 nulldistrib <- function(minN,maxN,p,w,binSize,yuplimit,distrib="binomial",b=0)
{
  d.combined = matrix(0,sum(seq(minN+1,maxN+1)),2)
  ptr = 1

  for (i in minN:maxN){
    ## doing the distribution
    k=seq(0,i)

    if(distrib == "binomial")
    {
      d = dbinom(k,i,p) ## binomial probability of all possible outcome
    }
    else if(distrib == "betabinomial")
    {
      d = dbetabinom(k,i,p,b)
    }

    ## weight each with actual counts in empirical
    d.w = d*w[i,1]

    if(i == minN)
    {
      d.combined[ptr:length(k),1] = k/i
      d.combined[ptr:length(k),2] = d.w
      colnames(d.combined) = c('allelicRatio','wBinDist')
    }
    else
    {
      d.combined[ptr:(ptr+length(k)-1),1] = k/i  # ref allele ratio
      d.combined[ptr:(ptr+length(k)-1),2] = d.w  # similar to FP count in AS test
    }

    ptr = ptr + length(k)
  }

  ## sort the d.combined distribution of all the n's
  d.combined.sorted = d.combined[ order(d.combined[,1],d.combined[,2]), ]

  ## bin it according to empirical distribution
  bins=pretty(0:1,binSize)
  start=0
  end=0

  d.combined.sorted.binned = matrix(0,length(bins)-1,2)

  for (z in 2:length(bins)) ##skip 0    ###sum up FP in one bins
  {
    start=bins[z-1]
    end=bins[z]

    row=z-1
    d.combined.sorted.binned[row,1] = (end-start)/2 +start ## equi of a $mid in hist

    d.combined.sorted.binned[row,2] = sum(d.combined.sorted[(d.combined.sorted[,1]<=end & d.combined.sorted[,1]>start),2])
    ## empirical right closed, left open ?hist; right=TRUE
    ## (range] so no double counts
    ## but zero gets excluded!!
    if(row==1)
    {
      d.combined.sorted.binned[row,2] = sum(d.combined.sorted[(d.combined.sorted[,1]<=end & d.combined.sorted[,1]>=start),2])
    }

    ## empirical right closed, left open ?hist; right=TRUE
    ## (range] so no double counts
    ## but zero gets excluded!!
    if(row==1)
    {
      d.combined.sorted.binned[row,2] = sum(d.combined.sorted[(d.combined.sorted[,1]<=end & d.combined.sorted[,1]>=start),2])
    }

  }

  ## change "counts" into density
  d.combined.sorted.binned[,2] = d.combined.sorted.binned[,2]/sum(d.combined.sorted.binned[,2])  #get probability

  return(d.combined.sorted.binned)
}

cutoff <- function(x,y) sum(y<=x)

fp <- function(w,p,p.thresh,distrib="binomial",b=0)
{
  ## doing the distribution; as.integer converts table entities to integers
  a=lapply(as.integer(w[,1]),function(x) seq(0,x))
  if(distrib == "binomial")
  {
    b = lapply(a,function(x) apply(as.data.frame(2*pbinom(x,max(x),p)),1,function(x) min(x,1)))  ###get binomial p value & make sure it's smaller than 1
  }
  else if(distrib == "betabinomial")
  {
    b = lapply(a,function(x) apply(as.data.frame(2*pbetabinom(x,max(x),p,rho)),1,function(x) min(x,1)))
  }
  ## find which ones are below threshold u
  d = lapply(b,function(x) x<=p.thresh)
  ## weight them by actual counts
  e = mapply(function(x,y,z) x*y*z, d, b, w[,2])
  ## sum up TP for all totals
  f = sapply(e,max)
  g = sum(f)

  return(g)
}

bisect <- function(p,p.sim,p.choice,fdr,fdr.threshold,by,distrib="binomial",b=0,w,p.thresh)
{
  p.fdr.e = matrix(0,100,3)
  e.prev = 10
  flag = 3
  ctr = 1
  p.fdr.e[ctr,1] = p.choice
  p.fdr.e[ctr,2] = fdr
  p.fdr.e[ctr,3] = e.prev

  while(flag)
  {
    start = max(0,(p.choice - by/2))
    end = p.choice + by/2
    by = by/4

    if(start==0){ start = 2.5e-5 } ## do not make it 0

    range = seq(start,end,by)

    for (i in range)
    {
      tp = cutoff(i,p)

      if(distrib == "binomial")
      {
        fp = fp(w,p.thresh,i,"binomial")
      }
      else if(distrib == "betabinomial")
      {
        fp = fp(w,p.thresh,i,"betabinomial",b)
      }

      fdr.ind = fp/tp
      e.curr = fdr.threshold - fdr.ind
      ctr = ctr + 1

      p.fdr.e[ctr,1] = i
      p.fdr.e[ctr,2] = fdr.ind
      p.fdr.e[ctr,3] = e.curr
      e.prev = p.fdr.e[(ctr-1),3]
      p.choice = i

      if(e.curr < 0){ break }

    }

    if(signif(p.fdr.e[ctr-1,3],3) == signif(p.fdr.e[ctr,3],3)){ flag = 0 }
  }
  return(p.fdr.e)
}


find_AS <- function(sample,rho,fdr)
{
  p=0.5
  FDR.thresh = fdr
  sample$p.bin = apply(data.frame(2 * mapply(pbinom,sample$lower,sample$total,p)),1,function(x) min(x,1))
  sample$p.betabin = apply(data.frame(2 * mapply(pbetabinom,sample$lower,sample$total,p,rho)),1,function(x) min(x,1))

  total <- sample$total
  w = as.data.frame(table(total), stringsAsFactors=F)  #frequency table of total reads

  p.thresh = data.frame( c(seq(0,0.01,by=0.001), seq(0.01,0.1,by=0.01)[-1], seq(0.1,1,by=0.1)[-1]) )
  # tx to tg

  tp.bin = apply(p.thresh,1,cutoff,y=sample$p.bin)+1                #number of SNPs with P < p.thresh for a series of p.thresh
  fp.bin = apply(p.thresh,1,function(x) fp(w,p,x,"binomial"))
  fdr.bin = fp.bin / tp.bin
  binomial = as.data.frame(cbind(p.thresh,fp.bin,tp.bin,fdr.bin))
  colnames(binomial) = c("pval","FP.bin","TP.bin","FDR.bin")

  tp.betabin = apply(p.thresh,1,cutoff,y=sample$p.betabin)+1
  fp.betabin = apply(p.thresh,1,function(x) fp(w,p,x,"betabinomial",b))
  fdr.betabin = fp.betabin / tp.betabin
  betabinomial = as.data.frame(cbind(p.thresh,fp.betabin,tp.betabin,fdr.betabin))
  colnames(betabinomial) = c("pval","FP.betabin","TP.betabin","FDR.betabin")

  FDR.txt = data.frame(cbind(p.thresh,binomial[,2],binomial[,3],binomial[,4],
                             betabinomial[,2],betabinomial[,3], betabinomial[,4]))
  colnames(FDR.txt) <- c("pval","FP.bin","TP.bin","FDR.bin",
                         "FP.betabin","TP.betabin","FDR.betabin")

  step= 0.0001

  p.choice.bin = max(p.thresh[,1][fdr.bin<=FDR.thresh])
  fdr.choice.bin = max(fdr.bin[fdr.bin<=FDR.thresh])

  p.choice.betabin = max(p.thresh[,1][fdr.betabin<=FDR.thresh])
  fdr.choice.betabin = max(fdr.betabin[fdr.betabin<=FDR.thresh])


  if (fdr.choice.bin == 0) {

    if (FDR.txt[2,3]==1){
      print (paste0("no SNP pass binomial FDR threshold, threshold = ",FDR.thresh))
      p.choice.bin.2 = "NULL"
      sample$bin.FDR = "not significant"
    } else{
      p.choice.bin.1 = as.data.frame(bisect(sample$p.bin,binomial[,2],p.choice.bin,fdr.choice.bin,FDR.thresh,step,"binomial",b=0,w,p))

      p.choice.bin.1 = p.choice.bin.1[p.choice.bin.1[,3]>0,]
      p.choice.bin.2 = p.choice.bin.1[nrow(p.choice.bin.1),1]
      for (i in 1:dim(sample)[1]){
        if (sample[i,10] < p.choice.bin.2){
          sample[i,12] = "significant"
        } else {
          sample[i,12] = "not significant"
        }
      }
    }
  } else {
    p.choice.bin.1 = as.data.frame(bisect(sample$p.bin,binomial[,2],p.choice.bin,fdr.choice.bin,FDR.thresh,step,"binomial",b=0,w,p))
    p.choice.bin.1 = p.choice.bin.1[p.choice.bin.1[,3]>0,]
    p.choice.bin.2 = p.choice.bin.1[nrow(p.choice.bin.1),1]
    for (i in 1:dim(sample)[1]){
      if (sample[i,10] < p.choice.bin.2){
        sample[i,12] = "significant"
      } else {
        sample[i,12] = "not significant"
      }
    }
  }

  if (fdr.choice.betabin ==0){
    if (FDR.txt[2,6]==1){
      print (paste0("no SNP pass beta-binomial FDR threshold, threshold = ",FDR.thresh))
      p.choice.betabin.2 = "NULL"
      sample$betabin.FDR = "not significant"
    } else{
      p.choice.betabin.1 = as.data.frame(bisect(sample$p.betabin,betabinomial[,2],p.choice.betabin,fdr.choice.betabin,FDR.thresh,step,"betabinomial",b,w,p))
      p.choice.betabin.1 = p.choice.betabin.1[p.choice.betabin.1[,3]>0,]
      p.choice.betabin.2 = p.choice.betabin.1[nrow(p.choice.betabin.1),1]
      for (i in 1:dim(sample)[1]){
        if (sample[i,11] < p.choice.betabin.2){
          sample[i,13] = "significant"
        } else {
          sample[i,13] = "not significant"
        }
      }
    }
  } else{
    p.choice.betabin.1 = as.data.frame(bisect(sample$p.betabin,betabinomial[,2],p.choice.betabin,fdr.choice.betabin,FDR.thresh,step,"betabinomial",b,w,p))
    p.choice.betabin.1 = p.choice.betabin.1[p.choice.betabin.1[,3]>0,]
    p.choice.betabin.2 = p.choice.betabin.1[nrow(p.choice.betabin.1),1]
    for (i in 1:dim(sample)[1]){
      if (sample[i,11] < p.choice.betabin.2){
        sample[i,13] = "significant"
      } else {
        sample[i,13] = "not significant"
      }
    }
  }
  colnames(sample)[12] <- "bin.FDR"
  colnames(sample)[13] <- "betabin.FDR"

  sample <- sample[order(sample$p.betabin),]


  a <- c("FDR.threshold ","p.choice.bin.old ","p.choice.bin ","p.choice.betabin.old ", "p.choice.betabin ")
  b <- c(FDR.thresh,p.choice.bin,p.choice.bin.2,p.choice.betabin, p.choice.betabin.2)
  ab <- data.frame(cbind(a,b))

  results <- list(Sample = sample, FDR.table = FDR.txt, thresholds = ab)
  return(results)
}


##### Main #####
p=0.5              
minN=6 

data <- process_data(file_name, minN)

if(max(data$total) < 2500){ maxN=max(data$total) }else { maxN=2500 }
apropor = length(data$total[data$total <= 2500]) / nrow(data)
yuplimit=0.15
binSize=40
bins=pretty(0:1,binSize)
r.min = 0
r.max = 1

 

## graded weights for SSE calculation
r = seq(r.min,r.max,(r.max - r.min)/((length(bins) - 1)/2))
r = r[2:length(r)]
if((length(bins)-1)%%2 != 0){ w.grad=c(r,sort(r[1:(length(r)-1)],decreasing=TRUE))
}else { w.grad=c(r,sort(r[1:length(r)],decreasing=TRUE)) }

## empirical allelic Ratio
data.match=data[data$total <= maxN & data$total >= minN, ]
h = hist(data.match$allelicRatio, xlim=range(0,1),breaks=bins,right=TRUE)
h$density =h$counts/sum(h$counts)

empirical = h$counts/sum(h$counts)


# weight by empirical counts
t = as.data.frame(table(data$total), stringsAsFactors=F)
w = matrix(0,max(data$total),1)

for (jj in 1:nrow(t)){
  w[as.integer(t[jj,1]),1] = t[jj,2]
}

d.combined.sorted.binned = nulldistrib(minN,maxN,p,w,binSize,yuplimit,distrib="binomial")

r.sta = 0
r.end = 0.99
r.by  = 0.1
b.range = seq(r.sta,r.end,by=r.by)
labels = matrix(0,50,1)
ctr = 1
sse = sum((empirical-d.combined.sorted.binned[,2])^2)
b.choice = 0
b.and.sse = matrix(0,50,2)
colnames(b.and.sse) <- c('b','sse')

for (k in b.range){
  e.combined.sorted.binned = nulldistrib(minN,maxN,p,w,binSize,yuplimit,distrib="betabinomial",b=k)
  # plot(e.combined.sorted.binned[,1],e.combined.sorted.binned[,2])
  ## minimize sse for betabinomials
  if(b.choice==0){ b.and.sse[1,1]=b.choice; b.and.sse[1,2]=sse }
  sse.bbin = sum(w.grad*((empirical-e.combined.sorted.binned[,2])^2))
  b.and.sse[ctr+1,1] = k     #rho
  b.and.sse[ctr+1,2] = sse.bbin    #SSE
  labels[ctr] = paste("betabin,b=",signif(k,2),"; SSE=",signif(sse.bbin,2))

  if(sse.bbin < sse){ sse = sse.bbin; b.choice = k }
  else if(sse.bbin > sse){ break }

  ctr = ctr + 1
}

ctr.ori = ctr             
b.and.sse.ori = b.and.sse

b.chosen = b.choice
sse.chosen = sse
flag = 3
if(b.chosen >= 0.9){flag = 0; newctr = ctr}

while(flag){
  r.sta = max(0,(b.choice - r.by/2))
  r.end = b.choice + r.by/2
  r.by  = r.by/4
  b.range = seq(r.sta,r.end,by=r.by)
  labels = matrix(0,50,1)
  newctr = 1
  sse = b.and.sse[1,2]
  b.choice = 0

  for (k in b.range)
  {
    e.combined.sorted.binned = nulldistrib(minN,maxN,p,w,binSize,yuplimit,distrib="betabinomial",b=k)

    ## minimize sse for betabinomials
    sse.bbin = sum(w.grad*((empirical-e.combined.sorted.binned[,2])^2))
    b.and.sse[(ctr+2),1] = k
    b.and.sse[(ctr+2),2] = sse.bbin
    labels[newctr] = paste("betabin,b=",signif(k,3),"; SSE=",signif(sse.bbin,3))

    if(sse.bbin < sse){ sse = sse.bbin; b.choice = k }
    else if(sse.bbin > sse){ break }

    ctr = ctr + 1
    newctr = newctr + 1
  }
  # print(paste("b.chosen=",b.choice,", SSE.chosen=",sse))
  labels = labels[1:(newctr+1),]

  if(signif(b.and.sse[ctr+2,2],3) == signif(b.and.sse[ctr+1,2],3)){ flag = 0 }
}

rho = signif(b.choice,3)


results <- find_AS(data,rho,fdr)

result <- results$Sample
write.table(result,file=paste0(strsplit(file_name,"\\.")[[1]][1],"_AS_", fdr,".txt"),sep="\t",row.names=FALSE,quote=FALSE)



FDR.table <- results$FDR.table
thresh <- results$thresholds

# density of sig AS SNPs with different AR

sub = result[result$bin.FDR=='significant',]
h.bin = hist(sub$allelicRatio, xlim=range(0,1),breaks=bins,right=TRUE)
empirical.bin = h.bin$counts/sum(h$counts)
sub=result[result$betabin.FDR=='significant',]
if (dim(sub)[1]>0){
  h.betabin=hist(sub$allelicRatio, xlim=range(0,1),breaks=bins,right=TRUE)
  empirical.betabin = h.betabin$counts/sum(h$counts)
} else {
  empirical.betabin = rep(0,50)
}

jpeg(file = paste0(strsplit(file_name,"\\.")[[1]][1], '_', fdr, ".jpeg"), height = 600, width = 800)
par(cex.axis=1, cex.lab=1, cex.main=1.5)
barplot(empirical, main=paste0(title, " (rho = ", signif(b.choice,3), ")"),ylab='density', xlab='allelic Ratio',names.arg=h$mids, ylim=c(0,1.15*yuplimit), xaxt="n")
par(new=TRUE)
barplot(empirical.bin, ylab='density', xlab='allelic Ratio', col='red',
        names.arg=h$mids, ylim=c(0,1.15*yuplimit), xaxt="n")
par(new=TRUE)
barplot(empirical.betabin, ylab='density', xlab='allelic Ratio', col='blue',
                     names.arg=h$mids, ylim=c(0,1.15*yuplimit), xaxt="n")
par(new=TRUE)
plot(d.combined.sorted.binned,ylim=c(0,yuplimit),pch=16,type='b',col='red',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i", cex = 0.8, xaxt="n")
par(new=TRUE)
plot(e.combined.sorted.binned,ylim=c(0,yuplimit), pch=17,type='b',col='blue',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i", cex = 0.8, xaxt="n")
axis(1,seq(0.01, 0.97, by=0.06),labels=c('0', rep(' ', 7), 0.5, rep(' ', 7),'1')) 

legend(0.01,0.14,c("empirical","binomial","beta-binomial"),
       col=c("grey","red","blue"),
       cex=1, pt.cex=2,text.col = "black", pch = 15, bg = 'white')
dev.off()
