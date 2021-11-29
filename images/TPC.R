## BA unimodal model

rBAU=function(temp,par_br,uni=1) {
  #input: temp=temperature range
  #       par_br=parameter vector
  #parameters:
  r0=par_br[[1]]  #scaling coefficient
  m=par_br[[2]]   #body-mass
  b=par_br[[3]]   #allometric scaling coefficient
  tpk=par_br[[4]] #optimal temperature
  E=par_br[[5]]   #activation energy
  Ed=par_br[[6]]  #deactivaiton energy
  k=8.617*10^-5   #Boltzmann constant
  # uni : specify if unimodal or exponential BA model

  if (uni==0) {
    l=1
  } else {
    l=1/(1+exp(-1/(k*temp)*(Ed-(Ed/tpk+k*log(E/(Ed-E)))*temp))) #decline phase
  }
  return(BR=r0*m^b*exp(-E/(k*temp))*l)
}

## Temperature

temp_seq=seq(285,308,length.out=50)

## Body mass

mH=1.34*10^-2 # herbivore body-mass

topt=298
### Parameters temperature dependent

par_aPH=c(aPH0=2*10^11,mH,baPH=0.05,topt=298,E=0.65,E2=1.15)
aPH_seq=rBAU(temp=temp_seq,par_br=par_aPH)

svg("images/TPC.svg", width=5,height=5)
dev.new()
plot(temp_seq,aPH_seq,lwd=2,axes=F,ann = F, type='l')
abline(v=topt)
dev.off()
