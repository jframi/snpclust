clust.cmeans<-function(snp,  X.col= "X.Fluor", Y.col="Y.Fluor", call.col="Call", Rmin=1 , centers=3, seuil=0.8,  cols=c("red","green","blue","black","gray")){

  X.numcol<-match(X.col,colnames(snp))
  Y.numcol<-match(Y.col,colnames(snp))

  #make a theta R version of X and Y
  snpTR <- xy2ThetaR(snp[,c(X.numcol,Y.numcol)])

  # Filter snp on Rmin value
  snp.filt<-snp[snpTR$R>Rmin,]

  snp.tocm<-snp.filt

  snp.filt.cm<-cmeans(snp.tocm[,c(X.numcol,Y.numcol)],centers=centers)
  mbships<-sapply(seq_along(snp.filt.cm$cluster), function(a) snp.filt.cm$membership[a,snp.filt.cm$cluster[a]])

  plot(snp.filt[,X.numcol],snp.filt[,Y.numcol],xlim=c(0,max(snp[,X.numcol])),ylim=c(0,max(snp[,Y.numcol])))

  #Plot filtered points as black
  points(snp[snpTR$R<Rmin,9:10],col="black")

  #Plot result of clustering taking seuil into account
  colors<-cols[as.numeric(snp.filt.cm$cluster)]
  colors[mbships<=seuil]<-"grey"
  points(snp.filt[,X.numcol],snp.filt[,Y.numcol],col=colors)
}


man.reclust<-function(snp, X.col= "X.Fluor", Y.col="Y.Fluor", call.col="Call", what, update.all=F, cols=c("red","green","blue","black","gray")){
  if (!what%in%levels(snp[,call.col])){
    stop(paste(what, "is not part of", call.col, "column.\n", "Please use one of ", paste(levels(snp[,call.col]),collapse = ", ")))
  }
  X.numcol<-match(X.col,colnames(snp))
  Y.numcol<-match(Y.col,colnames(snp))
  call.numcol<-match(call.col,colnames(snp))

  plot(snp[,X.numcol],snp[,Y.numcol], main=paste("Please draw a polygon for ", what, " classe", "\n click in negative X when finished"), cex.main=0.8)
  points(snp[,X.numcol],snp[,Y.numcol],col=cols[as.numeric(snp[,call.col])])
  polyg<-draw.polygon()
  pointsin<-point.in.polygon(snp[,X.numcol],snp[,Y.numcol],polyg$x,polyg$y)
  if (update.all){
    snp[snp[,call.numcol]==what,call.numcol]<-"Unknown"
    snp[pointsin==1,call.numcol]<-what
  } else{
    snp[pointsin==1,call.numcol]<-what
  }
  plot(snp[,X.numcol],snp[,Y.numcol])
  points(snp[,X.numcol],snp[,Y.numcol],col=cols[as.numeric(snp[,call.numcol])])
  polygon(polyg$x,polyg$y,border="red")
  return(snp)
}

clust.lga<-function(snp, ignore=NULL, X.col= "X.Fluor", Y.col="Y.Fluor",call.col="Call", Rmin=1 ,  k=3,seuil=0.8,  cols=c("red","green","blue","black","gray"), ...){

  X.numcol<-match(X.col,colnames(snp))
  Y.numcol<-match(Y.col,colnames(snp))

  #make a theta R version of X and Y
  snpTR <- xy2ThetaR(snp[,c(X.numcol,Y.numcol)])

  # Filter snp on Rmin value
  snp.filt<-snp[snpTR$R>Rmin,]

  if (!is.null(ignore)){
    snp.filt<-snp.filt[snp.filt[,call.col]!=ignore,]
  }
  snp.tocm<-snp.filt

  snp.filt.lga<-lga(snp.tocm[,c(X.numcol,Y.numcol)],k = k, ...)
  #plot(snp.filt[,c(X.numcol,Y.numcol)], xlim=c(0,max(snp[,X.numcol])),ylim=c(0,max(snp[,Y.numcol])))
  plot(snp.filt.lga, xlim=c(0,max(snp[,X.numcol])),ylim=c(0,max(snp[,Y.numcol])))
  points(snp[snpTR$R<Rmin,c(X.numcol,Y.numcol)],col="grey")
  #colors<-cols[as.numeric(snp.filt.lga$cluster)]
  #colors[mbships<=seuil]<-"grey"
  #points(snp.filt[,c(X.numcol,Y.numcol)],col=colors)
  #nclust <- attr(snp.filt.lga, "k")
  #for (i in 1:nclust)
  #  abline(lm(y~x, data=data.frame(x=snp.filt[snp.filt.lga$cluster==i,X.numcol],y=snp.filt[snp.filt.lga$cluster==i,Y.numcol])) , lty=2,col=(i + 1))

}
