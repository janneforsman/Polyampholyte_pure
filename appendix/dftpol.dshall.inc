      PARAMETER(maxel=12001,naxel=441,maxmon=301,maxism=21,naxmon=801)
      COMMON/VECT/fdmon(0:maxel),cdmon(0:maxel),fdsol(0:maxel),
     *cdsol(0:maxel),eblam(0:maxel),
     *Vexm(0:maxel),Vexs(0:maxel),Phimm(0:maxel),Phims(0:maxel),
     *Phiss(0:maxel),fem(0:maxel),ehbclam(0:maxel),
     *ebelam(0:maxel),
     *ae1(0:maxel),ae2(0:maxel),convp(0:maxel),
     *convs(0:maxel),Ymix(0:maxel),ae0(0:maxel),
     *etbclam(0:maxel),ehbelam(0:maxel),
     *convsm(0:maxel),convss(0:maxel),convmm(0:maxel),convms(0:maxel),
     *Vexrm(0:maxel),aehbclam(0:maxel),vdonn(0:maxel),ttfem(0:maxel),
     *tfdmon(0:maxel),fdgr(0:maxel)
      COMMON/VAR/dz,closew,pis,pif,pit,vk,rrT,vkrrT,hvk,scalem,scales,
     *dzpie,AA1,AA2,BB1,BB2,c1,c2,Y,emscale,eblemb,ehblcmb,eblsmb,bcdt,
     *rrjkdiff,threqz,rtwelve,pie,rthree,rdz,btrams,
     *sclosew,q1,q2,q3,p1,p2,p3,r2,r1,r0,s2,s1,s0,b2,b1,b0,r2sq,r1sq,
     *r0sq,Yfact,veq,rnmon,rrnmon,rrcmon,rq3,cdmbulk,cdsbulk,
     *cdmlbulk,cdslbulk,elblemb,elhblcmb,elblsmb,distp1,dmitt,eplag,
     *seplag,rrnarm,rnarm,tscalem,Pcollt,h,rrejkdiff,Rudifac,donn,fgrbr
      COMMON/HELTAL/istart,istp1,islut,ism,inw,nfack,imitt,nmon,
     *ist,ifin,istp1s,isluts,isms,inws,kst,kfin,iefin,iest,
     *kefin,kest




