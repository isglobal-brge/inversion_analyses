getFreqs <- function(x, index.inv, ...){
  xx <- lapply(index.inv, getFreqs.i, x=x, ...)
  ans <- data.frame(pop=attr(xx[[1]], "pops"), xx)
  colnames(ans)[-1] <- names(x)[index.inv]
  out <- as_tibble(ans)
  out
}


getFreqs.i <- function(i, x, grouping, ...){
  ff <- function(x){
    ans <- summary(x)$allele.freq[,2]
    return(ans["I"])
  }
  
  if (missing(grouping))
    pop <- x$country
  else
    pop <- x[ , which(names(x)==grouping)]
  
  inv <- SNPassoc::snp(x[,i], ...)
  freqs <- aggregate(inv, list(pop=pop), ff)
  names(freqs)[2] <- names(x)[i]
  ans <- freqs[,2]
  attr(ans, "pops") <- freqs[,1]
  ans
}

getHWE <- function(x, i){
  ff <- function(x){
    summary(x)$HWE
  }
  
  pop <- x$country
  inv <- SNPassoc::snp(x[,i])
  hwe <- aggregate(inv, list(pop=pop), ff)
  names(hwe)[2] <- "hwe"
  rownames(hwe) <- hwe$pop
  hwe
}

getData <- function(x, diseases, i, invs, ii, pca=TRUE, ...){
  ii <- which(colnames(invs)==ii)
  freqs <- getFreqs(invs, ii, ...)
  mask <- x$cause == diseases[i]
  x.i <- x[mask,]
  x.1000G <- x.i[x.i$location%in%rownames(freqs),
                 c("location", "year", "val")]
  data <- merge(x.1000G, freqs, by.x="location", by.y="pop")
  data2 <- data[data$year=="2016", ]
  if (pca)
    ans <- merge(data2, res.pca$x[,1:2], by.x="location", by.y="row.names")
  else
    ans <- data2
  ans
}

removeIC <- function(x){
  ans <- as.numeric(sapply(strsplit(as.character(x), " \\["), "[[", 1))
  ans
}  


plotSuperpop <- function(data, x, y, xlab="Inversion Frequency (%)", 
                         ylab="Age-standardized rate", pos.leg="topright") {
  xx <- data[ , which(colnames(data)==x)] 
  yy <- data[ , which(colnames(data)==y)]
  plot(xx, yy, type="n", xlab=xlab, ylab=ylab)
  points(xx,  yy,  pch=21, fg="grey40", bg=data$superpop, 
#         cex=(1/data$n)*400)
cex=1.4)
  wordcloud::textplot(xx, yy, data$lab, new=FALSE, pos=4)
  
  mm <- mixMod(data, x, y)
  o <- order(xx)
  lines(xx[o], predict(mm, level=0)[o], lwd=3, lty=2, col="darkgreen")
  rr <- round(sjstats::r2(mm)$r2*100, 1)
  pval <- formatC(summary(mm)$tTable[2,5])
  
  legend(pos.leg, c(eval(substitute(expression(R^2 == rr), list(rr=rr))), 
                    paste("p-value (mixed-model)= ", pval, sep="")), bty="n")
}

plotCountry <- function(data, x, y, xlab="Inversion Frequency (%)", 
                        ylab="Age-standardized rate"){
  xx <- data[ , which(colnames(data)==x)] 
  yy <- data[ , which(colnames(data)==y)]
  plot(xx, yy, type="n", xlab=xlab, ylab=ylab)
  points(xx,  yy,  pch=21, fg="grey40", bg=data$superpop, cex=2)
  wordcloud::textplot(xx, yy, data$lab, new=FALSE, pos=4)
}

mixMod <- function(data, i, k){
  x <- data[, which(colnames(data)==i)]
  y <- data[, which(colnames(data)==k)]
  z <- data$superpop
  dd <- data.frame(x=x, y=y, z=z)
  d.s <- groupedData(y ~ x | z, data=dd)
  mod <- lme(y ~ x, data=d.s, random = ~ 1 )
  mod
}

ecologicalAssoc_lmm <- function(x, invs.names, diseases) {
  n <- length(diseases)
  out <- list()
  kk <- 1
  for (k in invs.names){
    ans <- matrix(NA, nrow=n, ncol=2)
    ii <- 1
    for (i in diseases) {
      sel <- c(which(colnames(x)==i),
               which(colnames(x)==k),
               which(colnames(x)=="superpop"))
      d <- x[, sel]
      names(d) <- c("y", "x", "z")
      d.s <- groupedData(y ~ x | z, d)
      mod <- try(lme(y ~ x, data=d.s, random = ~ 1 ), TRUE)
      if (!inherits(mod, "try-error")){
        ans[ii, 1] <- summary(mod)$tTable[2,1]
        ans[ii, 2] <- summary(mod)$tTable[2,5]
      }
      ii <- ii + 1
    }
    rownames(ans) <- diseases
    colnames(ans) <- c("effect", "p-value")
    out[[kk]] <- ans
    kk <- kk + 1
  }  
  names(out) <- invs.names
  out
}


ecologicalAssoc_lm <- function(data, invs.names, diseases) {
  n <- length(diseases)
  out <- list()
  kk <- 1
  for (k in invs.names){
    ans <- matrix(NA, nrow=2, ncol=n)
    ii <- 1
    for (i in diseases) {
      sel <- c(which(colnames(data)==i),
               which(colnames(data)==k))
      d <- data[, sel]
      names(d) <- c("y", "x")
      mod <- try(lm(y~x, data=d), TRUE)
      if (!inherits(mod, "try-error")){
        ans[1, ii] <- summary(mod)$coefficients[2,1]
        ans[2, ii] <- summary(mod)$coefficients[2,4]
      }
      ii <- ii + 1
    }
    colnames(ans) <- diseases
    rownames(ans) <- c("effect", "p-value")
    out[[kk]] <- ans
    kk <- kk + 1
  }  
  names(out) <- invs.names
  out
}

getOBvar <- function(x){
  xcat <- cut(x, c(-Inf, 18, 25, 29, 39, Inf), 
              labels = 1:5)
  ans <- rep(NA, length(x))
  ans[xcat%in%c(2,3)] <- 0
  ans[xcat%in%c(4,5)] <- 1
  ans
}

addPheno <- function(x1, x2, bmicat=FALSE){
  ans <-merge(x1, x2, 
              by.x="SID", by.y="ID_1")
  if (!bmicat){
    ans$BMICAT <- cut(ans$bmi, 
                      c(-Inf, 18, 25, 29, 39, Inf), labels = 1:5)
  }
  
  ans$ob <- rep(NA, nrow(ans))
  ans$ob[ans$casco==0 & ans$BMICAT%in%c(2)] <- 0
  ans$ob[ans$casco==0 & ans$BMICAT%in%c(4,5)] <- 1
  
  ans$ob.l <- as.numeric(ans$BMICAT)
  # ans$ob.l[ans$ob.l==1] <- NA
  
  ans$ob2 <- rep(NA, nrow(ans))
  ans$ob2[ans$BMICAT%in%c(2)] <- 0
  ans$ob2[ans$BMICAT%in%c(4,5)] <- 1
  
  ans$obDiab <- rep(NA, nrow(ans))
  ans$obDiab[ans$casco==1 & ans$BMICAT%in%c(4,5)] <- 1
  ans$obDiab[ans$casco==0 & ans$BMICAT%in%c(2)] <- 0
  
  ans
}


concBand <- function(N, alpha=0.01, ...)
{
  
  # Calculate 95% confidence intervals. The jth order statistic from a uniform(0,1) 
  # sample has a beta(j,n-j+1) distribution (Casella & Berger, 2002, 2nd edition, pg
  # 230, Duxbury)
  
  conf <- 1-(alpha/2)
  c95<- vector()
  c05<- vector()
  for(i in 1:N){
    c95[i] <- qbeta(conf,i,N-i+1)
    c05[i] <- qbeta(1-conf,i,N-i+1)
  }
  ans <- cbind(c05, c95)
  ans
}

shade <- function(x1, y1, x2, y2, color = "lightblue") {
  n <- length(x2)
  polygon(c(x1, x2[n:1]), c(y1, y2[n:1]), border = NA, 
          col = color)
}


qqplotBand <- function(x,  main=NULL, lab, pos,  cex1=1, cex2=1, cex3=1, ...)
{
  o <- -log10(sort(x,decreasing=F))
  e <- -log10( 1:length(o)/length(o) )
  names(o) <- names(x)
  names(e) <- names(x)
  plot(e,o, main=main, 
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)), type="n", 
       cex.lab=cex1, cex.axis=cex2)
  
  N <-length(x)
  band <- concBand(N, ...)
  
  shade(e, -log(band[,1],10), e, -log(band[,2],10))
  
  points(e, o, ...)
  
  lines(e, e, col="blue")
  
  if (!missing(lab)) {
    for (i in 1:length(lab))
      {
        text(e[lab[i]], o[lab[i]], lab[i], pos=pos[i], cex=cex3, font=3)
    }
  }
}

plotCircosInv <- function(x, tissue, ...){
  sel <- grep(tissue, x$Tissue)
  dd <- x[sel,]
  rownames(dd)<- 1:nrow(dd)
  dd$symbol <- sapply(dd$symbol, function(x) x[1])
  dd$sign <- ifelse(sign(dd$logFC)==1, "blue", "red")
  dd2 <- dd[,c("chromosome", "start", "symbol",
               "to.gr", "posInv", "symbol")] 
  dd2$chromosome <-gsub("chr", "", dd2$chromosome)
  names(dd2) <- c("chr1", "po1", "gene1", "chr2", "po2", "gen2")
  cols <- ifelse(dd$rearrangement, "darkorange", "purple")
  
  invsAgg <- dd[!duplicated(dd$Inversion), 
                c("to.gr", "posInv", "Inversion")]
  
  par(mar=c(2, 2, 2, 2))
  plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="")
  circos(R=300, type="chr", cir="hg19", col=TRUE, 
         print.chr.lab=TRUE, W=10, cex=2)
  circos(R=290, cir="hg19", W=40, 
         mapping=dd[,c("chromosome", "start", "symbol")], 
         type="label", side="in", cex=1)
  circos(R=300, cir="hg19", W=50, 
         mapping=invsAgg, col="black",
         type="label", side="out", cex=1.5) 
  circos(R=180, cir="hg19", W=30, 
         mapping=dd[,c("chromosome", "start", "logFC")], 
         type="b3", lwd=3, B=TRUE, cutoff = 0, 
         col=dd$sign) 
  circos(R=160, cir="hg19", W=50, mapping=dd2,  
         type="link", col=cols, lwd=1.6) 
  legend("bottomleft", c("Cis-effect", "Trans-effect"), 
         lwd=2, col=c("darkorange", "purple"), cex=1.5)
  legend("bottomright", c("Up-regulation", "Down-regulation"), 
         pch=16, col=c("blue", "red"), cex=1.5)
  title(tissue, cex.main=2)
}


plotInvExpr.old <- function(i, v, x, ...){
  gene <- fData(x)[i, "hgnc_symbol"]
  ee <- data.frame(expr = v$E[i,],
                   inv = x$inv)
  bp <- ggplot(ee, aes(x=inv, y=expr, group=inv)) + 
    geom_jitter(alpha = 0.3, color = "tomato", width=0.1) +
    geom_boxplot(alpha = 0) + xlab("Inversion 8p23") +
    ylab("Gene expression") + ggtitle(gene)
  bp
}

plotInvExpr <- function(i, x, model="additive", 
                        xlab="Inversion 8p23.1", ...){
  gene <- fData(x)[i, "hgnc_symbol"]
  if (model=="additive")
  ee <- data.frame(expr = exprs(x)[i,],
                   inv = x$inv)
  if (model=="recessive")
    ee <- data.frame(expr = exprs(x)[i,],
                     inv = recessive(x$inv))
  bp <- ggplot(ee, aes(x=inv, y=expr, group=inv)) + 
    geom_jitter(alpha = 0.3, color = "tomato", width=0.1) +
    geom_boxplot(alpha = 0) + xlab(xlab) +
    ylab("Gene expression") + ggtitle(gene)
  bp
}




getMinPval <- function(x, ...){
  x <- pvalues(x)[,c(3,4,6)]
  pval <- apply(x, 1, min, na.rm=TRUE)
  model <- apply(x, 1, function(x) 
    c("dominant", "recessive", "additive")[which.min(x[!is.na(x)])])
  out <- pval
  names(out) <- rownames(x)
  out
}

getMinPval2 <- function(x, ...){
  xx1 <- p.value(x, df=1)
  xx2 <- p.value(x, df=2)
  pval <- apply(cbind(xx1, xx2), 1, min, na.rm=TRUE)
  out <- pval
  names(out) <- names(x)
  out
  
}

plotGTeX <- function(x, tit){
  x$Symbol <- getGene(x$symbol)
  x <- x[!is.na(x$Symbol),]
  x$beta <- ifelse(x$logFC>0, "Positive", "Negative")
  ggplot(x, aes(x = Tissue, y = Symbol, size=-log10(adj.P.Val)))+     geom_point(aes(col = beta)) +
    scale_size_continuous(breaks=c(2,3,4,6), range=c(2,6) ) + 
    labs(colour = "Effect of I allele") + 
    scale_color_manual(values=c("red", "blue")) + 
    xlab("Tissue") + ylab("Gene") +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5,
                                     hjust = 1)) +
    ggtitle(tit)
  
}

getGene <- function(x) {
  temp <- sapply(x, function(x) strsplit(x, ";"))
  ans <- unlist(lapply(temp, "[", 1))
  ans
}

getAlleleRef <- function(x){
  ss <- strsplit(as.character(x), "/")
  ans <- lapply(ss, "[", 1)
  ans
}


getAlleleRefs <- function(x){
  sum.x <- summary(x, print=FALSE)
  refs <- unlist(getAlleleRef(sum.x$alleles))
  ans <- data.frame(inversion=rownames(sum.x),
                       ref = refs, stringsAsFactors = FALSE)
  ans
}


getBetasSign <- function(x){
  ans <- SNPassoc::odds(x)
  beta <- rep(NA, nrow(ans))
  names(beta) <- rownames(ans)
  if(!attr(x, "quantitative")){
    beta[ans$OR >= 1 & ans$`p-value.log-additive` <= 0.05] <- 1 
    beta[ans$OR < 1 & ans$`p-value.log-additive` <= 0.05] <- -1 
  }
  else{
    beta[ans$upper < 0 & ans$`p-value` <= 0.05] <- -1 
    beta[ans$lower > 0 & ans$`p-value` <= 0.05] <- 1 
  }
  beta 
}
  