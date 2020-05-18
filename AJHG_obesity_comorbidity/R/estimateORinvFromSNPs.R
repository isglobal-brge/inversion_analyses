# Haplotype frequencies given MAFs and Dprime between SNPs
hapFreq<-function(Dprime, pA, pB){
  Dmax<-min(c((1-pA)*pB,pA*(1-pB)))
  D<-Dmax*Dprime
  mm.haplo<-matrix(NA,2,2)
  mm.haplo[1,1]<-(1-pA)*(1-pB)+D
  mm.haplo[1,2]<-pA*(1-pB)-D
  mm.haplo[2,1]<-(1-pA)*pB-D
  mm.haplo[2,2]<-pA*pB+D
  mm.haplo
}

# Bivariate distribution of SNPs
genoFreq<-function(myHapFreqs){
  mm.geno<-matrix(NA,3,3)
  mm.geno[1,1]<-myHapFreqs[1,1]^2
  mm.geno[1,2]<-2*myHapFreqs[1,1]*myHapFreqs[1,2]
  mm.geno[1,3]<-myHapFreqs[1,2]^2
  mm.geno[2,1]<-2*myHapFreqs[1,1]*myHapFreqs[2,1]
  mm.geno[2,2]<-2*myHapFreqs[1,1]*myHapFreqs[2,2]+2*myHapFreqs[1,2]*myHapFreqs[2,1]
  mm.geno[2,3]<-2*myHapFreqs[1,2]*myHapFreqs[2,2]
  mm.geno[3,1]<-myHapFreqs[2,1]^2
  mm.geno[3,2]<-2*myHapFreqs[2,1]*myHapFreqs[2,2]
  mm.geno[3,3]<-myHapFreqs[2,2]^2
  mm.geno
}

# Conditional genotype frequencies
transFreq<-function(myGenoFreqs){
  myGenoFreqs/rowSums(myGenoFreqs)
}

# genotype frequencies given the MAF
genoProb<-function(maf){
  c((1-maf)^2,2*maf*(1-maf),maf^2)
}


# Penetrance: Prob(Desease | SNP)
genoPen<-function(ORaA,ORAA,pD,pA,OR=TRUE){
  genoProbA<-genoProb(pA)
  pa<-1-pA
  if (OR){
    ff<-function(x) x*genoProbA[1]+x*ORaA*genoProbA[2]/(1+x*(ORaA-1))+x*ORAA*genoProbA[3]/(1+x*(ORAA-1))-pD
    PrDgaa<-uniroot(ff,c(0,1))$root
    PrDgAa <-PrDgaa*ORaA/(1+PrDgaa*(ORaA-1))
    PrDgAA <-PrDgaa*ORAA/(1+PrDgaa*(ORAA-1))
  } else {
    denom <- ORAA * pA^2 + ORaA * 2 * pA * pa + pa^2
    PrDgaa <- pD/denom
    PrDgAa <- ORaA * PrDgaa
    PrDgAA <- ORAA * PrDgaa
  }
  return(c(PrDgaa,PrDgAa,PrDgAA))
}

# translate OR of genotyped SNP to causal SNP
transOR<-function(ORaA,ORAA,pD,pA,pB,Dprime,OR=TRUE){
  penA<-genoPen(ORaA,ORAA,pD,pA,OR=OR)
  myTransFreq<-transFreq(genoFreq(hapFreq(Dprime,pA,pB)))
  penB<-myTransFreq%*%penA
  if (OR){
    ORBb<-(penB[2]/(1-penB[2]))/(penB[1]/(1-penB[1]))
    ORBB<-(penB[3]/(1-penB[3]))/(penB[1]/(1-penB[1]))
    ans<-c(ORBb,ORBB)
  } else{ # Relative Risk
    ORBb<-penB[2]/penB[1]
    ORBB<-penB[3]/penB[1]
    ans<-c(ORBb,ORBB)
  }
  ans
}

# translate r2 to Dprime
r2toDprime <- function(r2, pA, pB){
  r <- sqrt(r2) # assume positive r
  D <- r*sqrt(pA*(1-pA)*pB*(1-pB))
  Dmax<-min(c((1-pA)*pB,pA*(1-pB)))
  Dprime <- D/Dmax
  return(Dprime)
}



## obtain Odds Ratio of causal SNP given the OR of genotyepd SNP

transOR2 <- function(or, pD, pB, pA, Dprime, OR=TRUE){
  
  # or <- 1.15 # odds ratio per allele of genotyped locus
  # pD <- 0.01 # desease prevalence (not very important)
  # pB <- 0.3 # genotyped SNP MAF
  # pA <- 0.1 # causal SNP (inversion) MAF
  # Dprime <- 0.3
  
  ORbB <- or
  ORBB <- or^2
  
  # solve
  f <- function(x) transOR(x[1], x[2], pD, pA, pB, Dprime, OR=TRUE) - c(ORbB, ORBB)
  sol <- multiroot(f = f, start = c(1, 1))$root
  return(sol)
  
}

