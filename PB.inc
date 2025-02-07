      PARAMETER (maxel=10001,naxel=1200,maxism=1200,maxmon=1201)
      COMMON/VECT/fdmon(0:maxel),fdsol(0:maxel),
     *ebelam(0:maxel),
     *Vexm(0:maxel),Vexs(0:maxel),Phimm(0:maxel),
     *Phiss(0:maxel),fem(0:maxel),ehbclam(0:maxel),
     *ae1(0:maxel),ae2(0:maxel),convp(0:maxel),cdtot(0:maxel),
     *fdes(0:maxel),
     *fdnmon(0:maxel),fenm(0:maxel),ebelanm(0:maxel),ehbclanm(0:maxel)
      COMMON/VAR/dz,closew,pis,pif,pit,rrT,vkrrT,hvk,scalem,scales,
     *dzpie,emscale,eblemb,ehblcmb,eblsmb,bcdt,
     *rrjkdiff,threqz,rtwelve,pie,rthree,rdz,btrams,
     *sclosew,elblemb,elhblcmb,elblsmb,distp1,dmitt,eplag,
     *seplag,tscalem,rnval,csurf,donn,scalec,
     *q1,q2,q3,p1,p2,p3,r2,r1,r0,s2,s1,s0,b2,b1,b0,r2sq,r1sq,
     *r0sq,rq3,cdmbulk,rp3,rnsval,ph,pqh2,pqh,rpqh3,pqh3,
     *dm1,dm2,dm3,rdm3,ds1,ds2,ds3,rds3,alpha,pi,gfact,
     *palpha,salpha,
     *srnsval,VL,VR,h,dhs,dhs3,rdhs3,dhs2,
     *bl,bl2,bl3,rbl3,
     *AA1,AA2,BB1,BB2,c1,c2,Y,scalees
      COMMON/HELTAL/istart,istp1,islut,ism,inw,nfack,imitt,
     *ist,ifin,istp1s,isluts,isms,inws,kst,kfin,ismms,ismm,
     *nfin,ibl




