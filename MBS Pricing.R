price_MBS <- function(Principal,MR,term,r0,r_mu,sigma,Kappa){
  r=MR/12
  an=(1-(1/(1+r))^term)/(r) # standard $1 annuity monthly
  PMT=Principal/an # actual monthly pmt without default
  
  Term=seq(0,term,1)
  loan_t=SP=IP=c(0)
  loan_t[1]=100000
  
  # construct CIR rt matrix
  dt=1/12
  m=term+122
  rt=matrix(0,1000,m)
  rt[,1]=r0
  
  #anti-variate control
  set.seed(123456789)
  dw=matrix(rnorm(1000*m/2,0,1), ncol=m, nrow=1000/2)
  dw=rbind(dw, -dw)
    
  for(k in 2:(m)){
    # Monte carlo on interest rate CIR 
    rt[,k]= abs(rt[,(k-1)]+Kappa*(r_mu-rt[,(k-1)])*dt+sigma*sqrt(rt[,(k-1)])*sqrt(dt)*dw[,(k-1)])
    #rt[,k]=ifelse(rt[,k]<=0, rt[,k-1], rt[,k])
  }
  
  # construct CPR_t
  CPR=matrix(0,1000,term)
  RI=BU=SG=PP=c(0)
  SY=c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.1,1.18,1.22,1.23,0.98)
  price=c(0)
  for(i in 1:1000){
    Ct=ZCB=c(0)
    for(k in 1:(term)){
      RI[k]=0.28+0.14*atan(-8.57+430*(MR-sum(rt[i,k:(k+120)])/120))
      BU[k]=0.3+0.7*loan_t[k]/100000
      SG[k]=min(1,k/12/30)
      
      # assuming loan inception in Jan
      month=k%%12
      month=ifelse(month==0, 12, month) # assign month number

      CPR[i,k]=RI[k]*BU[k]*SG[k]*SY[month]
      SP[k]=PMT-loan_t[k]*r
      PP[k]=(loan_t[k]-SP[k])*(1-(1-CPR[i,k])^(1/12))
      
      Ct[k]=PMT+PP[k]
      loan_t[k+1]=loan_t[k]*(1+r)-Ct[k]
    }
    
    
    # to keep track of loan balance, stop when paid-off
    j=1
    loan=loan_t
    while(loan[j]>0){
      loan[j+1]=loan[j]*(1+r)-Ct[j]
      ZCB[j]=Ct[j]*exp(-sum(rt[i,2:j])*dt)
      j=j+1
    }
    
    ZCB[j]=loan[j-1]*exp(-sum(rt[i,2:(j+1)])*dt)
    price[i]=sum(ZCB)
  }
  
  return (mean(price))
}
