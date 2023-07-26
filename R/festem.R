rm.outlier <- function(x,percent){
  tmp <- (x>=quantile(x,percent))
  Q13 <- quantile(x[tmp],c(0.25,0.75))
  if (Q13[2]==Q13[1]) Q13[2] <- Q13[1]+1
  upper <- max(Q13[2]*3-2*Q13[1],1)
  lower <- max(Q13[1]*3-2*Q13[2],0)
  outlier.index <- tmp & (x>upper | x<lower)
  mean.noout <- mean(x[tmp & (!outlier.index)])
  x[outlier.index] <- round(mean.noout)
  cat(sum(outlier.index),"\n")
  x
}

sub.sample <- function(x){
  library.size <- x[1]
  x <- x[-1]
  sample.count <- function(a,library.size){
    rpois(1,a/library.size)
  }
  apply(matrix(x), 1, sample.count,library.size = library.size)
}

my.min <- function(x){
  if (sum(is.na(x))==length(x)){
    NA
  } else {
    min(x,na.rm = T)
  }
}

my.max <- function(x){
  if (sum(is.na(x))==length(x)){
    NA
  } else {
    max(x,na.rm = T)
  }
}

festem <- function(counts,cluster.labels,batch.id,
                   prior.weight = 0.05, prior.weight.filter = 0.9,
                   earlystop = 1e-4, outlier_cutoff = 0.90,
                   min.percent = 0.01,min.cell.num = 30,
                   seed = 321,num.threads = 2){
  library(parallel)
  library(edgeR)
  set.seed(seed)
  cluster.labels <- factor(cluster.labels)
  levels(cluster.labels) <- 1:nlevels(cluster.labels)
  batch.id <- factor(batch.id)
  
  cl <- makeCluster(getOption("cl.cores", num.threads))
  clusterSetRNGStream(cl, iseed = seed)
  em.result <- vector("list",nlevels(batch.id))
  em.result.f <- vector("list",nlevels(batch.id))
  for (B in 1:nlevels(batch.id)){
    counts.tmp <- counts[,batch.id==(levels(batch.id)[B])]
    cluster.labels.tmp <- cluster.labels[batch.id==(levels(batch.id)[B])]
    # Removing those clusters with no more than 10 cells in this batch
    cluster.labels.tmp <- factor(cluster.labels.tmp)
    counts.tmp <- counts.tmp[,summary(cluster.labels.tmp)[cluster.labels.tmp]>10]
    cluster.labels.tmp <- cluster.labels.tmp[summary(cluster.labels.tmp)[cluster.labels.tmp]>10]
    cluster.labels.tmp <- factor(cluster.labels.tmp)
    levels(cluster.labels.tmp) <- 1:nlevels(cluster.labels.tmp)
    
    ## Outlier
    counts.tmp <- t(parApply(cl,counts.tmp,1,rm.outlier,percent = 0.95))
    rownames(counts.tmp) <- rownames(counts)
    
    ## Sub-sampling
    library.size <- calcNormFactors(counts.tmp)
    counts.tmp <- parApply(cl,rbind(library.size,counts.tmp),2,sub.sample)
    counts.tmp <- t(parApply(cl,counts.tmp,1,rm.outlier,percent = outlier_cutoff))
    rownames(counts.tmp) <- rownames(counts)
    
    nonzeros.num <- function(x){sum(x!=0)}
    tmp <- apply(counts.tmp, 1, nonzeros.num)
    min.cell <- min(min.percent*ncol(counts.tmp),min.cell.num)
    #sum(tmp>=min.cell)
    counts.tmp <- counts.tmp[tmp>=min.cell,]
    
    alpha.label <- numeric(nlevels(cluster.labels.tmp)-1)
    for (g in 1:length(alpha.label)) {
      alpha.label[g] <- sum(cluster.labels.tmp==g)/ncol(counts.tmp)
    }
    
    time.tmp <- Sys.time()
    em.result[[B]] <- parApply(cl,counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                               group.num = nlevels(cluster.labels.tmp),prior.weight=prior.weight,earlystop = earlystop)
    print(paste0("Batch ",B," -- ","Time cost: ",difftime(Sys.time(),time.tmp,units = "secs")))
    
    time.tmp <- Sys.time()
    em.result.f[[B]] <- parApply(cl,counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                                 group.num = nlevels(cluster.labels.tmp),prior.weight=prior.weight.filter,earlystop = earlystop)
    print(paste0("Batch ",B," -- ","Time cost: ",difftime(Sys.time(),time.tmp,units = "secs")))
    
  }
  
  ## EM
  p.tmp <- matrix(NA,nrow = nrow(counts),ncol = length(em.result))
  rownames(p.tmp) <- rownames(counts)
  for (i in 1:length(em.result)){
    p.tmp[colnames(em.result[[i]]),i] <- em.result[[i]][1,]
  }
  p.tmp <- parApply(cl,p.tmp,1,my.min)
  p.tmp <- p.tmp*length(em.result)
  
  em.tmp <- matrix(NA,nrow = nrow(counts),ncol = length(em.result))
  rownames(em.tmp) <- rownames(counts)
  for (i in 1:length(em.result)){
    em.tmp[colnames(em.result.f[[i]]),i] <- em.result.f[[i]][2,]
  }

  em.tmp <- apply(em.tmp,1,my.max)
  stopCluster(cl)
  
  result.EM = data.frame(names = rownames(counts), p = p.adjust(p.tmp,"BH"), EM = em.tmp)
  tmp.na <- result.EM[is.na(result.EM$p),1]
  result.EM <- result.EM[!is.na(result.EM$p),]
  gene.names <- result.EM[result.EM$p<0.05 & result.EM$EM>0,]
  gene.names <- gene.names[order(-gene.names$EM),]
  
  tmp <- result.EM[result.EM$p>=0.05 & result.EM$EM>0,]
  tmp <- tmp[order(-tmp$EM),]
  gene.names <- rbind(gene.names,tmp)
  tmp <- result.EM[result.EM$EM<=0,]
  tmp <- tmp[order(tmp$p,-tmp$EM),]
  gene.names <- rbind(gene.names,tmp)
  gene.names <- gene.names[,1]
  gene.names <- c(gene.names,tmp.na)
  EM <- gene.names
  
  list(stat = result.EM,genelist = EM)
}