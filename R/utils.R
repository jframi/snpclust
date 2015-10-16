xy2ThetaR<-function(x){
  r<-apply(x,1,function(a) sqrt(a[1]^2+a[2]^2))
  theta<-apply(x,1,function(a) atan(a[2]/a[1]))
  return(data.frame( Theta=theta,R=r))

}

draw.polygon<-function(){
  p<-locator(1)
  pret<-p
  while (p$x>0){
    p2<-locator(1)
    if (p2$x>0){
      segments(p$x,p$y,p2$x,p2$y, col="red")
      pret$x<-c(pret$x,p2$x)
      pret$y<-c(pret$y,p2$y)
    }
    p<-p2
  }
  segments(pret$x[1],pret$y[1],tail(pret$x,1),tail(pret$y,1), col="red")
  return(pret)
}
