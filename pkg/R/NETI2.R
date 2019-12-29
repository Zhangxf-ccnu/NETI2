
#' @title Joint reconstruction of multiple gene NETworks by simultaneously capturing Inter-tumor and Intra-tumor heterogeneity
#' @description The complete procedure for reconstructing multiple gene networks
#'  using NETI2. For details, refer to Supplementary Section S3.2.
#' @import Matrix mvtnorm foreach QUIC
#' @usage NETI2(X,purity,lambda, tau, delta)
#' @param X  A list (length = \eqn{K}) of data matrices (\eqn{n_k \times p}), where K is the number of cancer subtypes and \eqn{n_k} is the sample size of k-th cancer subtype.
#' @param purity A list (length = \eqn{K}) of the tumor purity information vectors (\eqn{n_k \times 1}), where K is the number of cancer subtypes
#' and \eqn{n_k} is the sample size of k-th cancer subtype.
#' @param lambda  The common tuning parameter for controlling the overall degree sparsity of the estiamted gene networks. For details, refer to Supplementary TableS6.
#' @param tau   The tuning parameter balances the network size between non-cancerous and cancerous networks. For details, refer to Supplementary TableS6.
#' @param delta  The tuning parameter balances the network size between subtype specific networks. For details, refer to Supplementary TableS6.
#' @details The function is used to jointly reconstruction of multiple gene networks by simultaneously capturing inter-tumor
#'    and intra-tumor heterogeneity. For each cancer subtype, the observed gene  expression levels of tumor samples are assumed to be a
#'  mixture of expressions from non-cancerous and cancerous cells. The gene  expression levels of cancerous cells are assumed to follow subtype
#'  specific multivariate normal distributions, and the gene expression levels of non-cancerous cells across all subtypes are assumed to follow the
#'  same multivariate normal distribution. The precision matrices of these multivariate normal distributions are used to build the non-cancerous and subtype
#'  specific cancerous networks.  Given the observed gene expression data and tumor purity information data,
#' we use a penalized likelihood approach to estimate the model parameters. We develop an efficient iterative procedure based on Expectation
#' Maximization (EM) algorithm to solve the optimization problem with latent (unobserved) variables.
#'Each iteration of the EM algorithm consists of two steps: E-step and M-step. For details, refer to Supplementary Section S3.2.



#' @return
#' \item{\code{theta.y}}{A list (length = \eqn{K}) of estimated precision matrices from cancerous cells of different  cancer subtyes.}
#' \item{\code{theta.z}}{A matrix of estimated  precision matrice from non-cancerous cells shared by all cancer subtyes.}
#' \item{\code{LL.temp}}{Log-likelihood of the data for different EM iterations.}
#'
#' @references  Jia-Juan Tu, Le Ou-Yang, Hong Yan, Xiao-Fei Zhang and Hong Qin (2019), Joint reconstruction of multiple gene networks by simultaneously capturing inter-tumor
#'      and intra-tumor heterogeneity
#' @author Jia-Juan Tu
#' @seealso { \code{\link{generate.data}}, \code{\link{TCGA.BRCA}}}
#' @export
#' @examples # Simulation data
#' data.x= generate.data(p = 100, n = 100, K = 4, network.type ="ER", umin = 0.5, umax = 1)
#' result = NETI2(data.x$X,data.x$purity, lambda = 0.6, tau = 0.5,delta = 0.5)
#'
#' # TCGA breast cancer data
#' data("TCGA.BRCA")
#' result = NETI2(TCGA.BRCA$X,TCGA.BRCA$purity, lambda = 0.6, tau = 0.4,delta = 0.2)

NETI2<- function(X,purity,lambda = 1.6, tau = 0.4,delta = 0.2){

  K=dim(X)[1]
  lambday=lambda*tau
  lambdaz=lambda*(1-tau)
  re = EM1_whole(X,purity)
  Epurity = re$Epurity

  paraMY_mu = re$paraM_ma[,1:K]
  paraMY_sig = re$paraM_ma[,(K+2):(2*K+1)]

  paraMZ_mu = re$paraM_ma[,(K+1)];
  paraMZ_sig = as.matrix(re$paraM_ma[,(2*K+2)])

  re_em2 = EMalg(X,Epurity,paraMY_mu,paraMY_sig,paraMZ_mu,paraMZ_sig,lambday, lambdaz,delta)


  result = list(theta.y=re_em2$thetay, theta.z=re_em2$thetaz, LL.temp= re_em2$LL.temp)

}


### The included sub-functions are listed as following:
#(1) The mian function to estimate the sample purity
#library("foreach")
EM1_whole <- function (X, purity){

  K=dim(X)[1]
  p=ncol(X[[1]])
  n=matrix(0,K,1)
  err_tol=1e-3

  for (ii in 1:K){  n[ii]=  nrow(X[[ii]])}
  sumn=sum(n)

  xv=matrix(list(), K,1)
  uy=matrix(0, K,1)
  vy=matrix(0, K,1)
  uz_temp=matrix(0, K,1)
  vz_temp=matrix(0, K,1)

  paraM <-foreach(pi= 1:p)  %do%
  {

    for (ik in 1:K)  { xv[[ik]] = as.matrix(X[[ik]][,pi]) }

    for (ik in 1:K)  {
      uy[ik]=mean(xv[[ik]])
      uz_temp[ik]=mean(xv[[ik]])
      vy[ik]=var(xv[[ik]])
      vz_temp[ik]=var(xv[[ik]]) }

    uz = mean(uz_temp)
    vz = mean(vz_temp)

    for(i in 1:10000)
    {
      uz_old=uz
      uy_old=uy
      vz_old=vz
      vy_old=vy

      for (ik in 1:K)  {

        temp<-Estep.1d(uy[ik],uz,vy[ik],vz,1,xv[[ik]],purity[[ik]])

        mstep_temp=Mstep.1d(temp$Ey,temp$Ey2,temp$Ez,temp$Ez2)

        uy[ik]= mstep_temp$uy
        uz_temp[ik]=mstep_temp$uz

        vy[ik]=mstep_temp$vy
        vz_temp[ik]=mstep_temp$vz

      }

      uz = mean(uz_temp)
      vz = mean(vz_temp)


      if(abs(uy - uy_old)< abs(uy_old)*err_tol && abs(vy- vy_old)< abs(vy_old)*err_tol &&
         abs(uz - uz_old)< abs(uz_old)*err_tol && abs(vz- vz_old)< abs(vz_old)*err_tol)

      {
        break()}

    }
    c(uy,uz, vy,vz)

  }

  paraM_ma<-do.call("rbind",paraM)


  Epurity<-purity
  result = list(Epurity=Epurity,  paraM_ma= paraM_ma)

}

##
Estep.1d<-function(uy,uz,vy,vz, c,X,A)
{
  P<-c*A
  Ux<-P*uy+(1-P)*uz
  Vx<-P^2*vy+(1-P)^2*vz
  Sxy<-P*vy
  Ey<-uy+(X-Ux)*Sxy/Vx
  Vy<-vy-Sxy^2/Vx

  Ey2<-Vy+Ey^2

  Ez<-(X-Ey*P)/(1-P)
  Ez2<-(X^2-2*P*X*Ey+P^2*Ey2)/(1-P)^2

  list(Ey=Ey,Ey2=Ey2, Ez=Ez, Ez2=Ez2)
}
##
Mstep.1d<-function(Ey,Ey2,Ez,Ez2)
{
  uy<-mean(Ey)
  uz<-mean(Ez)
  vy<-mean(Ey2)-uy^2
  vz<-mean(Ez2)-uz^2
  list(uy=uy,uz=uz,vy=vy,vz=vz)
}

#(2) The mian function to estimate the gene networks from different cancer subtypes

EMalg<-function (X,Epurity,paraMY_mu,paraMY_sig,paraMZ_mu,paraMZ_sig,lambda_y, lambda_z,delta){

   convergecutoff=1e-4
  K=dim(X)[1]
  p=dim(X[[1]])[2]
  n=matrix(0,K,1)

  for (ii in 1:K){  n[ii]=  nrow(X[[ii]])}
  sumn=sum(n)
  nbar=sumn/K
  n_bar=(nbar^delta)*(n^(1-delta))
  lambday=lambda_y*n_bar
  lambdaz=lambda_z* sumn

  p=dim(X[[1]])[2]
  ESigmay=matrix(list(), K,1)
  ESigmaz=matrix(list(), K,1)
  Ey=matrix(list(), K,1)
  Ez=matrix(list(), K,1)
  Uy=matrix(list(), K,1)
  Sigmay=matrix(list(), K,1)
  thetay=matrix(list(), K,1)
  Uz_temp=matrix(list(), K,1)
  LL.temp=c()


  for (i in 1: K){
    Uy[[i]]=paraMY_mu[,i]
    Sigmay[[i]]=diag(paraMY_sig[,i])  }
  Uz=paraMZ_mu
  Sigmaz=diag(paraMZ_sig[,1])

  for (i in 1:200)

  {

    for (ij in 1: K){

      ESigmay_temp<-Estep.noloop(Uy[[ij]],Uz,Sigmay[[ij]],Sigmaz,1,X[[ij]],Epurity[[ij]])

      ESigmay[[ij]]=ESigmay_temp$ESigmay
      ESigmaz[[ij]]=ESigmay_temp$Esigmaz

      Ey[[ij]]=ESigmay_temp$Ey
      Ez[[ij]]=ESigmay_temp$Ez

    }


    for (ik in 1:K){

      Uy[[ik]] = apply(Ey[[ik]],2,mean)

      temp=matrix(Uy[[ik]], nrow=n[ik] , ncol=p, byrow=T)
      Sy=t(Ey[[ik]]-temp)%*%(Ey[[ik]]-temp)/n[ik] + ESigmay[[ik]]
      Sy_p=matrix(nearPD(Sy)$mat@x,p,p)


      lambdaym=matrix(lambday[ik],p,p)
      diag(lambdaym)=0
      re_y<-QUIC(Sy_p, lambdaym/n[ik],msg = 0, maxIter = 500)
      Sigmay[[ik]]=re_y$W
      thetay[[ik]]=re_y$X
    }


    Uz_temp[[ik]] = apply(Ez[[ik]],2,sum)
    Uz=do.call("rbind", Uz_temp)
    Uz=apply(Uz,2,sum)/sumn

    Sz=0
    for (k in 1:K){
      temp=matrix(Uz, nrow=n[k] , ncol=p, byrow=T)
      Sz=Sz+t(Ez[[k]]-temp)%*%(Ez[[k]]-temp) +ESigmaz[[k]]*n[k]
    }

    Sz_p=matrix(nearPD(Sz/sumn)$mat@x,p,p)


    lambdazm=matrix(lambdaz,p,p)
    diag(lambdazm)=0
    re_z<-QUIC(Sz_p,lambdazm/sumn,msg = 0, maxIter = 500)
    Sigmaz=re_z$W
    thetaz=re_z$X;

    sumlog= cal_logL(Uy,Uz,Sigmay,Sigmaz,X,Epurity)
    LL.temp<-c(LL.temp,sumlog)


    if(i>1 && abs(LL.temp[i-1]-LL.temp[i])<convergecutoff*abs(LL.temp[i-1]))
    { break()}

    i<-i+1

  }

  result = list(thetay = thetay, thetaz = thetaz, LL.temp = LL.temp)
}

##
Estep.noloop<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,AA)
{

  P<-as.vector(c.scale*AA)  ### percentage estimate
  P.m=diag(P)
  p=ncol(X)
  n=nrow(X)



  Ey<-matrix(NA,n, p)
  Ez<-Ey
  ESigmay<-matrix(0,p,p)
  Esigmaz<-ESigmay

  temp=eigen(Sigmay, symmetric=TRUE)
  M=temp$values
  U=t(temp$vectors)
  ### temp1=t(U)%*%diag(M)%*%U
  ### max(abs(temp1-Sigmay))
  Sy.half=t(U)%*%diag(M^0.5)
  Sy.half.inv=diag(1/M^0.5)%*%U
  ### max(abs(Sy.half%*%t(Sy.half)-Sigmay))

  W=Sy.half.inv%*%Sigmaz%*%t(Sy.half.inv)
  temp1=eigen(W, symmetric=TRUE)
  V=t(temp1$vectors)
  D=temp1$values
  ### temp2=t(V)%*%diag(D)%*%V
  ### max(abs(temp2-W))

  A=Sy.half%*%t(V)
  B=V%*%Sy.half.inv ## p by p

  BX=B%*%t(X)%*%P.m ##p by n
  BUy=(B%*%Uy)%*%matrix(P^2, 1, n) ## p by n
  BUz=(B%*%Uz)%*%matrix(P*(1-P),1,n) ## p by n
  tempB=BX-BUy-BUz ## p by n

  tempP2=matrix(P^2, p,n, byrow=TRUE) ##p by n
  tempC=tempP2+matrix(D, p,1)%*%matrix((1-P)^2, 1,n) ## p by n
  tempC=1/tempC

  tempBC=tempB*tempC ##p by n
  tempABC=A%*%tempBC ## p by n

  Ey=t(matrix(Uy, p,n)+tempABC) ##  n by p

  Cbar= diag(apply(tempP2*tempC, 1, mean)) ##p by p diagonal
  ESigmay=Sigmay-A%*%Cbar%*%t(A)

  P.2m=diag(1/(1-P))
  Ez=P.2m%*%(X-P.m%*%Ey)
  tempP4=matrix(P^4/(1-P)^2, p,n, byrow=TRUE) ##p by n
  CbarZ= diag(apply(tempP4*tempC, 1, mean)) ##p by p
  Esigmaz=mean(P^2/(1-P)^2)*Sigmay-A%*%CbarZ%*%t(A) ## p by p


  ESigmay<-(ESigmay+t(ESigmay))/2  ### garantee its symmetry
  Esigmaz<-(Esigmaz+t(Esigmaz))/2  ### garantee its symmetry

  #end=proc.time()
  list(Ey=Ey,Ez=Ez,ESigmay=ESigmay,Esigmaz=Esigmaz)
}
##
cal_logL<-function(Uy,Uz,Sigmay,Sigmaz,X,P){
  K=dim(X)[1]

  sumlogl=0
  for (ik in 1:K){
    re_logL<-calculateLL(Uy[[ik]],Uz,Sigmay[[ik]],Sigmaz,X[[ik]],P[[ik]])
    sumlogl=sumlogl+re_logL

  }

  sumlogl=sumlogl
}

#####
calculateLL<-function(Uy,Uz,Sigmay,Sigmaz,X,P)
{

  sum(sapply(1:nrow(X),function(i)
  {
    Ux<-P[i]*Uy+(1-P[i])*Uz
    Sigmax<-P[i]^2*Sigmay+(1-P[i])^2*Sigmaz
    dmvnorm(X[i,], Ux, Sigmax, log = TRUE)
  }))

}
