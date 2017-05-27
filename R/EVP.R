#'Optimization Algorithm to solve stratification problem
#'
#'@param X        Values of the population associated with stratification variable.
#'@param L        Number of strata.
#'@param cv       Value with target cv associated with stratification variable, between 0 and 1.
#'@param nhmin    Smallest possible sample size in any stratum. The default is 1.
#'@param Nmin     Smallest possible population size in any stratum. The default is 2.
#'@param imax     Maximum number of iterations. The default is 150.
#'@param tmax     Maximum number neighbors. The default is 7.
#'@param pmax     Maximum number of solutions - pool. The default is 5.
#'@param notbest  Maximum number of iterations without improvement. The default is 25.
#'@param range_s  Shaking interval. The default is 30.
#'@param range_b  Search interval. The default is 20.
#'@param cpu_time Maximum cpu time. The default is 5000.
#'@param cores    Numerical amount of CPUs requested for the cluster. The default is 2.
#'@return \item{TYPE}{Type of solution. If the solution was obtained by the exact method
#'                    returns "GLOBAL OPTIMA" else the solution was obtained by the metaheuristic
#'                     method returns "LOCAL OPTIMA".}
#'@return \item{bh}{Strata boundaries.}
#'@return \item{Nh}{Population size in each stratum.}
#'@return \item{n}{Total sample size.}
#'@return \item{nh}{Sample size in each stratum.The last one is a take-all stratum, then sample size
#'                  equal to population size.}
#'@return \item{cv}{Estimated coefficient of variation for the estimator of total of the stratification variable.}
#'@return \item{cputime}{Time consumed by the algorithm in seconds.}
#'@references Brito, J.A.M, Silva, P.L.N.,Semaan, G.S. and Maculan, N. (2015).
#'            Integer Programming Formulations Applied to Optimal Allocation in Stratified Sampling.
#'            Survey Methodology, 41: 427-442.
#'
#'            Brito, J.A.M, Semaan, G.S., Fadel, A.C. and Brito, L.R.(2017).
#'            An optimization approach applied to the optimal stratification problem,
#'            Communications in Statistics - Simulation and Computation.
#'
#'	          Hansen, P.; Mladenovic, N.; Perez-Brito, D. (2001).
#'	          Variable Neighborhood Decomposition Search.
#'	          Journal of Heuristics 7.4, pp. 335-350.
#'
#'           Glover, F.; Kochenberger, G. A., eds. (2003).
#'	         Handbook of Metaheuristics. Springer.
#'
#'@author Breno Trotta de Oliveira(brenotrotta@yahoo.com.br), Leonardo Lima and Jose Brito.
#'
#'@description
##'An Optimization Algorithm Applied to Univariate Stratification Problem.
#'It's aims to delimit the population strata and defining the allocation of sample,
#'considering  the following objective: minimize the sample size given a fixed precision level.
#'Exhaustive enumeration method is applied in small problems, while in problems with greater
#'complexity the algorithm is based on metaheuristics Variable Neighborhood Decomposition
#'Search with Path-Relinking.
#'
#'@export
#'@import Rglpk
#'@import snowfall
#'@import stratification
#'@import sampling
#'@import MultAlloc
#'@importFrom stats var
#'@importFrom utils combn
#'@examples
#'data(Sweden)
#'P75<-Sweden[,3]
#'solution1<-EVP(P75,L=2,cv=0.1,nhmin=5,imax=50,cores=2)
#'solution2<-EVP(P75,L=6,cv=0.1,cores=2,imax=25)

EVP<-function(X,L,cv,nhmin=1,Nmin=2,imax=150,tmax=7,pmax=5,notbest=25,range_s=30,range_b=20,cpu_time=5000,cores=2)
{
  calcula_Vh_Nh<-function(bk,X,H)
  {Vh<-rep(0,H)
  Nh<-rep(0,H)
  bk<-c(-Inf,bk,Inf)
  for(i in 1:H)
  {Xh<-which(X>bk[i] & X<=bk[i+1])
  Vh[i]<-var(X[Xh])
  if (is.na(Vh[i])==TRUE) {Vh[i]<-Inf}
  Nh[i]<-length(Xh)
  }
  return(cbind(Nh,Vh))
  }

  enumera_estratos<-function(X,H,Nmin,nmin,cvt)
  {
    tempo<-proc.time()
    b<-sort(unique(X))
    bk<-combn(b,H-1)
    u<-t(apply(bk,2,function(bk) calcula_Vh_Nh(bk,X,H)))
    notfeasible<-which(apply(u[,1:H],1,min)<Nmin)
    if (length(notfeasible)>0)
    {u<-u[-notfeasible,]
    bk<-bk[,-notfeasible]
    }
    tx<-sum(X)
    sb<-snowfall::sfApply(u,1,function(x) BSSM_FC(x[1:H],x[((H+1):(2*H))],tx,cvt,nmin,certain=T))
    ix<-which.min(unlist(lapply(sb,function(x) x$n)))
    tempo<-(proc.time()-tempo)[3]
    sfStop()
    if (H==2) {bh=bk[ix]}
    else {bh=bk[,ix]}
    return(list(bh=bh,Nh=u[ix,1:H],n=unlist(sb[[ix]][1]),nh=unlist(sb[[ix]][2]),cv=unlist(sb[[ix]][3]) , tempo=tempo))
  }

  perturbacao=function(X,L,S0,aleat,Q,range_s,w,Nmin)
  {
    posicoes_iniciais<-match(S0,Q)
    raio=seq(-range_s,range_s,1);
    raio=raio[raio!=0];
    itera=0
    repeat{
      itera=itera+1
      if (itera > 1000) {S1=S0; break}
      passo = sample(raio,1);
      posicoes_novas <- posicoes_iniciais
      posicoes_novas[aleat==1] <- posicoes_novas[aleat==1]+passo
      posicoes_novas<-sort(posicoes_novas)
      if (posicoes_novas[1] <= 0) next
      if  (posicoes_novas[L-1] > w) next
      S1=Q[posicoes_novas]
      u<-calcula_Vh_Nh(S1,X,L)
      naoviavel<-which(min(u[,1])<Nmin)
      if (length(naoviavel)==0 & any(S1!=S0)==T) {break}
    }
    return(S1)
  }

  buscalocal=function(X,L,S0,S1,aleat,Q,range_b,w,Nmin,k)
  {
    posicoes_iniciais<-match(S1,Q)
    raio1=seq(-range_b,range_b,1);
    raio1=raio1[raio1!=0];
    itera=0
    repeat{
      itera=itera+1
      if (itera > 50000) {S2=S1;u<-calcula_Vh_Nh(S2,X,L); break}
      passo1 = sample(raio1,k);
      passos=aleat
      passos[aleat==1] <- passo1
      posicoes_novas <- posicoes_iniciais + passos
      posicoes_novas<-sort(posicoes_novas)
      if (posicoes_novas[1] <= 0) next
      if  (posicoes_novas[L-1] > w) next
      S2=Q[posicoes_novas]
      u<-calcula_Vh_Nh(S2,X,L)
      naoviavel<-which(min(u[,1])<Nmin)
      if (length(naoviavel)==0 & any(S2!=S0)==T ) {break}
    }
    return(list(b=S2,Nh=u[,1],Vh=u[,2]))
  }

  cria_pool=function(Solucao, matriz)
  {
    if (length(matriz)>0) {
      if (any(Solucao!=matriz[1,])==T) {A<-rbind(Solucao,matriz)}
      else return(matriz)
    }
    else {
      A=NULL
      A=matrix(Solucao,nrow=1,byrow=T)
    }
    if (dim(A)[1] >pmax) {A=A[1:pmax,]}
    row.names(A)=NULL
    return(A)
  }

  path_relink<-function(pool,H,B,X,Nmin)
  {
    reconexao<-function(sb,sg,H)
    {

      conexoes<-NULL
      for(i in 1:(H-1))
      {
        ci<-seq(sb[i]+1,sg[i]-1)
        sbs<-cbind(ci,matrix(rep(sb[-i],length(ci)),ncol=H-2,byrow=TRUE))
        conexoes<-rbind(conexoes,sbs)
        sb[i]<-sg[i]
      }
      return(conexoes)
    }

    if (is.null(dim(pool)) == FALSE & length(pool) > (H-1) )
    {
      bks<-apply(pool,2,function (x) match(x,B))
      combinacoes<-combn(dim(bks)[1],2)
      if (dim(bks)[1]== 2)
      {w<-reconexao(bks[combinacoes[1],],bks[combinacoes[2],],H)
      w2=w
      }
      else {
        w<-apply(combinacoes,2,function(x) reconexao(bks[x[1],],bks[x[2],],H))
        w2=NULL
        for (i in 1:length(w)) {
          for (j in 1:( length(w[i]) / (H-1) ) ) {w2<-rbind( w2, w[[i]][j,] )}
        }
      }
      v<-t(apply(w2,1,function(x) duplicated(x)))
      id<-which(apply(v,1,function(x) all(x==FALSE))==TRUE)
      wm<-w2[id,]
      wm<-t(apply(wm,1,function(y) sort(y)))

      u1<-apply(wm,1,function(x) paste(x,sep="",collapse=""))
      id2<-which(duplicated(u1)==FALSE)
      solucoes<-wm[id2,]

      cuts<-matrix(B[t(solucoes)],ncol=H-1,byrow=TRUE)

      u<-t(apply(cuts,1,function(z) calcula_Vh_Nh(z,X,H)))
      naoviavel<-which(apply(u[,1:H],1,min)<Nmin)
      if (length(naoviavel)>0)
      {u<-u[-naoviavel,]
      cuts<-cuts[-naoviavel,]
      }
      return(list(matriz_b=cuts,NHSH=u))
    }
    else
      return()
  }

  Q=sort(unique(X)) # vetor com todos os valores populacionais de PO deduplicados
  w=length(Q) #qtd de valores distintos
  snowfall::sfInit(parallel=TRUE,cpus=cores)


  # A ================================== METODO EXATO =======================================================
  if ( choose(w,L-1)<1000001 & length(X)<4001 )
  {
    S=enumera_estratos(X,L,Nmin,nhmin,cv)
    names(S$n)=NULL;names(S$nh)=NULL;names(S$cv)=NULL;
    snowfall::sfStop()
    return(list(TYPE="GLOBAL OPTIMA", bh=S$bh,Nh=S$Nh,n=S$n,nh=S$nh,cv=S$cv, cputime =S$tempo))
  }

  else
  {
    # B ==================================== METODO VNDS ======================================================

    tempo0<-proc.time()
    kmax=L-1
    POOL=NULL
    it_pr=seq(from=5,to=imax,by=20)
    repeat {
      repeat {
        S0=sort(Q[sample(1:w,size=kmax)])
        u<-calcula_Vh_Nh(S0,X,L)
        naoviavel<-which(min(u[,1])<Nmin)
        if (length(naoviavel)==0) {break}
      }

      calcula_ini<- BSSM_FC(u[,1],u[,2],sum(X),cv,nhmin,certain=T)
      if (calcula_ini$n < 100 | L==2) {break}
    }
    S=list(bh=S0,n=calcula_ini$n,nh=calcula_ini$n,cv=calcula_ini$cvs)

    itera=1
    n_acum=S$n

    repeat {
      k=1
      for (k in 1:kmax)
      {

        aleat=srswor(k,kmax)

        #B.02============= SHAKING

        S1=perturbacao(X,L,S$bh,aleat,Q,range_s,w,Nmin)

        z<-as.list(rep(k,tmax))
        St<-lapply(z,function(x) buscalocal(X,L,S$bh,S1,aleat,Q,range_b,w,Nmin,x))

        U=NULL
        for (t in 1:tmax){
          U<-rbind(U,unlist(St[[t]]))
        }


        calcula=lapply(St,function(x)  BSSM_FC(x$Nh,x$Vh,sum(X),cv,nhmin,certain=T))
        ix<-which.min(unlist(lapply(calcula,function(x) x$n)))
        S2=list(bh=U[ix,1:(L-1)],Nh=U[ix,L:((2*L)-1)],n=unlist(calcula[[ix]][1]),nh=unlist(calcula[[ix]][2]),cv=unlist(calcula[[ix]][3]))
        names(S2$bh)=NULL;names(S2$Nh)=NULL;names(S2$n)=NULL;names(S2$nh)=NULL;names(S2$cv)=NULL;


        if (S2$n <= S$n) {S=S2; break}
        else {S=S;k=k+1}

      }



      POOL=cria_pool(S$bh,POOL)


      for ( i in it_pr )
      {
        if (itera==i & L>2)
        {
          PR=path_relink(POOL,L,Q,X,Nmin)
          if (is.null(PR)==FALSE)
          {
            sb<-apply(PR$NHSH,1,function(x) MultAlloc::BSSM_FC(x[1:L],x[((L+1):(2*L))],sum(X),cv,nhmin,certain=T))
            ix<-which.min(unlist(lapply(sb,function(x) x$n)))

            S3=list(bh=PR$matriz_b[ix,], Nh=PR$NHSH[ix,1:L], n=unlist(sb[[ix]][1]), nh=unlist(sb[[ix]][2]), cv=unlist(sb[[ix]][3]))
            names(S3$bh)=NULL;names(S3$Nh)=NULL;names(S3$n)=NULL;names(S3$nh)=NULL;names(S3$cv)=NULL;

            if (S3$n <= S$n) {S=S3}
          }

        }
      }


      n_acum=rbind(n_acum,S$n)

      tempo1<-proc.time()-tempo0;
      if (itera >= notbest) {media=mean(n_acum[(length(n_acum)-notbest+1):length(n_acum)])}
      else media=0

      if (itera==imax |
          tempo1[3] > cpu_time |
          n_acum[length(n_acum)]==media
      )   {break}
      else itera=itera+1

    }
    snowfall::sfStop()
    return(list(
      TYPE="LOCAL OPTIMA",
      bh=S$bh,Nh=S$Nh,n=S$n,nh=S$nh,cv=S$cv,
      cputime =tempo1[3]
    ))
  }
}
