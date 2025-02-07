      program platem
      implicit double precision (a-h,o-z)
      include 'PB.inc'
      ifc = 38
      ife = 36
      ifnc = 34      
      ins = 49
      issolv = 56
      isd = 58
      idpot = 63
      ksurf = 65
      ikh = 67      
      pi = acos(-1.d0)
      onethird = 1.d0/3.d0
      twothirds = 2.d0/3.d0
      bk = 1.38066D-23
      avno = 6.02214D23
      epszero = 8.854D-12
      elch = 1.602D-19
      T = 298.d0
      epsilon = 78.3
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      rphi = fourpi*0.2d0 
      aphi = fourpi*0.5d0
      pis = pi/6.d0
      pit = pi/3.d0
      pif = pi/4.d0     
      es22 = -32.d0*pi/9.d0
      a1 = 1.d0
      a2 = 2.45696d0
      b1 = 1.d0
      b2 = 4.10386d0
      c1 = -1.d0
      c2 = -3.75503
      AA1 = 2.d0*c1-2.d0*a1-4.d0
      AA2 = 2.d0*c2-2.d0*a2-4.d0
      BB1 = 3.d0-b1+a1-3.d0*c1
      BB2 = 3.d0-b2+a2-3.d0*c2
      Y = (9.82605d0-9.d0*pi*0.25d0)/(9.d0*pi*0.25d0-4.d0*pi/3.d0)      
      ddtol = 0.0000001d0
      delecok = 1.D-9
      smet = 1.d-10
      unorm = elch*elch/(4.d0*pi*epszero*epsilon*smet)
      bkT = bk*T
      T = bkT/unorm
      rrT = 1.d0/T
      bjerrum = rrT
      rbjerrum = rrT
      write(*,*) 'unorm,bkT /J = ',unorm,bkT
      write(*,*) 'T,bjerrum = ',T,bjerrum
c     CLOSE TO THE WALLS, THE DENSITY IS ASSUMED TO BE ZERO
      open (ife,file='fesfil',form='formatted')         
      open (ins,file='dpcinp',form='formatted')
      open (idpot,file='donnan',form='formatted')
      open (ksurf,file='surfdens',form='formatted')      
      open (ikh,file='separation',form='formatted')                  
      rewind idpot      
      rewind ife
      rewind ins
      rewind ksurf
      rewind ikh
      read(ins,*) 
      read(ins,*) bdmtrams,bdes
      bdm = 0.d0
      read(ins,*) 
      read(ins,*) nval,nsval
      nval = 1
      nsval = 1
      read(ins,*) 
      read(ins,*) nmontrams
      nmon = 0
      if (mod(nmon,2).eq.0) then
      abspolval = 0.d0
      else
      abspolval = 1.d0
      endif
      write(*,*) 'abs(polval) = ',abspolval
      frdonn = 0.001d0
      read(ksurf,*) surfdens      
      read(ins,*) 
      read(ins,*) htrams
      read(ikh,*) h      
      read(ins,*)
      read(ins,*) dz
      read(ins,*)      
      read(ins,*) dmm
      read(ins,*)      
      read(ins,*) dms
      read(ins,*)      
      read(ins,*) kread
      read(ins,*)      
      read(ins,*) ioimaxm
      read(ins,*)
      read(ins,*) trams
      read(ins,*)
      read(ins,*) dm1,ds1,bl     
      write(*,*) 'dm1,ds1 = ',dm1,ds1
      dhs = dm1
      write(*,*) 'dhs,bl = ',dhs,bl
      rbl = 1.d0/bl
      bl2 = bl*bl
      bl3 = bl2*bl
      rbl3 = 1.d0/bl3
      dhs2 = dhs*dhs
      dhs3 = dhs2*dhs
      rdhs3 = 1.d0/dhs3         

      rnval = real(nval)
      rrnval = 1.d0/rnval
      rnsval = real(nsval)
      rrnsval = 1.d0/rnsval
      halfh = 0.5d0*h
      nhalfh = int(halfh/dz+0.01)
      cbdm = 0.d0
      bds = bdes
      bdnm = bdm
      write(*,*) 'surdfdens = ',surfdens 
      write(*,*) 'bdm,bdes = ',bdm,bdes
      write(*,*) 'bdnm,bds = ',bdnm,bds
      tdmm = 1.d0-dmm
      tdms = 1.d0-dms
c      rrT = 1.d0/T
      rdz = 1.d0/dz
c     twopidz = twopi*dz
      twopidz = 0.5d0*dz*rbl      
      irdz = int(rdz+0.001d0)
      nfack = int((h+dm1)/dz+0.01d0)
      ibl = int(bl/dz+0.01d0)      
      imitt = int((h+dm1)/dz+0.01d0)/2
      closew = 0.5d0*dm1
      sclosew = closew      
      istart = int(closew/dz+0.01d0)
      iblnw = ibl+int(closew/dz+0.01d0)

      csurf = -2.d0*pi*surfdens*rrT*(h+dm1)          
      write(*,*) 'csurf = ',csurf            

      nfin = 2*imitt
      istp1 = istart+1 
c     islut = nfack
      islut = int((h+0.5d0*dm1)/dz+0.01d0)            
      istp1s = istp1
      isluts = islut
      inw = istp1+ibl-1      
      write(*,*) 'h = ',h
      write(*,*) 'closew,sclosew = ',closew,sclosew
      write(*,*) 'ism,ibl = ',ism,ibl
      write(*,*) 'nfack,nfin = ',nfack,nfin
      write(*,*) 'istp1,islut = ',istp1,islut
      write(*,*) 'imitt = ',imitt

      dm2 = dm1*dm1
      dm3 = dm2*dm1
      rdm3 = 1.d0/dm3
      ds2 = ds1*ds1
      ds3 = ds2*ds1
      rds3 = 1.d0/ds3      
      phzmm = -(6.d0*pi*dm2/5.d0)
      rphi = rphi*dm2
      aphi = aphi*dm2
      phzss = phzmm
      phzms = phzmm

      es22ss = 0.d0
      es22mm = es22ss
      es22ms = es22ss
      Pblj = 0.5d0*es22mm*bdm**2+0.5d0*es22ss*bds**2+
     *bdm*bds*es22ms

      q1 = ds1
      p1 = dm1
      q2 = q1*q1
      q3 = q2*q1
      p2 = p1*p1
      p3 = p2*p1
      dm2 = dm1*dm1
      dm3 = dm2*dm1
      rp3 = 1.d0/p3
      rq3 = 1.d0/q3
      rdm3 = 1.d0/dm3
      write(*,*) 'dm1,q1,p1 = ',dm1,q1,p1            
      ism = int(dm1/dz+0.01d0)
      isms = int(ds1/dz+0.01d0)
      ismms = int(p1/dz+0.01d0)
      write(*,*) 'ism,isms,ismms = ',ism,isms,ismms

      bdt = 2.d0*bdes*dhs3
      xsib = 1.d0-pis*bdt
      rxsib = 1.d0/xsib
      rxsibsq = rxsib*rxsib
      aex1= -(c1+1.d0)*dlog(xsib)-
     *0.5d0*(AA1*pis*bdt+BB1*(pis*bdt)**2)*rxsibsq
      aex2= -(c2+1.d0)*dlog(xsib)-
     *0.5d0*(AA2*pis*bdt+BB2*(pis*bdt)**2)*rxsibsq
      Vdae1dV = -pis*bdt*rxsib*(c1+1.d0-
     *0.5d0*(AA1+2.d0*BB1*pis*bdt)*rxsib-
     *pis*bdt*rxsibsq*(AA1+BB1*pis*bdt))
      Vdae2dV = -pis*bdt*rxsib*(c2+1.d0-
     *0.5d0*(AA2+2.d0*BB2*pis*bdt)*rxsib-
     *pis*bdt*rxsibsq*(AA2+BB2*pis*bdt))
      
      rNsda1 = pis*bds*rxsib*((c1+1.d0)-
     *0.5d0*AA1*rxsib*(1.d0+2.d0*pis*bdt*rxsib)-
     *BB1*pis*rxsib*bdt*(1.d0+pis*bdt*rxsib))

      rNesda1 = pis*bdes*rxsib*((c1+1.d0)-
     *0.5d0*AA1*rxsib*(1.d0+2.d0*pis*bdt*rxsib)-
     *BB1*pis*rxsib*bdt*(1.d0+pis*bdt*rxsib))
      
      baex1 = aex1
      baex2 = aex2
      exchps = aex1+(rNsda1+rNesda1)*dhs3
      exPs = -bds*Vdae1dV-bdes*Vdae1dV
      exP = exPs

      chemps = dlog(bds)+exchps
      chempes = dlog(bdes)+exchps
      Pbhs = bds+bdes+exPs
      Pb = Pbhs

      hsbconvs = (rNsda1+rNesda1)*dhs3
      trams = hsbconvs      
      strams = hsbconvs+aex1
      estrams = hsbconvs+aex1            
      scales = chemps
      scalees = chempes      
      
      hsbeblam = dexp(-strams+scales)
      hsbeeblam = dexp(-strams+scalees)      
      write(*,*) 'HSBCONVS = ',hsbconvs
      write(*,*) 'strams = ',strams
      write(*,*) 'bds = bdes = ',bds,bdes      
      write(*,*) 'max no. of iterations = ',ioimaxm
      write(*,*) 'bjerrum,rbjerrum = ',bjerrum, rbjerrum
      write(*,*) 'h,dz = ',h,dz
      write(*,*) 'coion chemical potential = ',chemps
      write(*,*) 'counterion chemical potential = ',chempes      
      write(*,*) 'total bulk pressure = ',Pb      
      write(*,*) 'rrT = ',rrT
      write(*,*) 'nval,nsval = ',nval,nsval
      write(*,*) 'surfdens = ',surfdens
      write(*,*) 'istp1,islut,inw = ',istp1,islut,inw
      write(*,*) 'ism,ibl = ',ism,ibl
      write(*,*) 'dmm,dms (density mixing param. mon.,solv.) = ',dmm,dms
      write(*,*) 'bond length (bl): ',bl
      write(*,*) 'monomer AND solvent hs diameter (dhs): ',dhs
      write(*,*) 'NOTE: MONOMERS/SOLVENT PARTICLES HAVE THE SAME SIZE!'       

      z = 0.5d0*dz
      iz = 1
      zp = z-dz
c      rewind 23
      do jz = 1,maxel
      zp = zp+dz
      diffz = abs(z-zp)
      iii = iabs(jz-iz)
      phiz = -2.d0*pi*diffz*rrT*dz
      Phimm(iii) = phiz
c      write(23,*) iii,Phimm(iii),jz
      enddo

      do i = 1,istp1
      fdsol(i) = 0.d0
      fdes(i) = 0.d0
      enddo
      
      if (kread.eq.0) then
      donn = 0.d0
      do iz = istp1,imitt       
      fdsol(iz) = bds
      fdes(iz) = bdes      
      enddo
      else
      read(idpot,*) donn         
c     do iz = istp1,imitt
      do iz = 1,imitt                
      read(ife,*) trams,fdsol(iz),fdes(iz)           
      enddo
      endif

      sum = 0.d0
      do i = istp1,imitt
      sum = sum-fdsol(i)+fdes(i)
      enddo
      elec = surfdens+sum*dz
      write(*,*) 'surfdens,elec (st) = ',surfdens,elec          

      elecA = elec      
      write(*,*) 'elecA,donn = ',elecA,donn
c     if (dabs(elecA).lt.delecok) goto 985
      if (dabs(elecA).lt.delecok) then
      ddonn = 0.d0
      donnB = donn
      elecB = 0.d0
      goto 985
      else      
      donnA = donn
      if (dabs(donn).lt.0.001d0) then
      ddonn = 0.05d0
      else
      ddonn = frdonn*donn
      endif
      if (elecA.lt.0.d0) then
      donnB = donnA+ddonn
      else
      donnB = donnA-ddonn
      endif
      ddd = donnB
 4949 continue
      donnB = ddd
      write(*,*) 'donnB,ddd = ',donnB,ddd
      sumfds = 0.d0
      sumfdes = 0.d0
      do i = istp1s,imitt
      sumfdes = fdes(i)*dexp(-(donnB-donn))+sumfdes      
      sumfds = fdsol(i)*dexp((donnB-donn))+sumfds
      enddo
      elecB = (sumfdes-sumfds)*dz+surfdens
      ddd = donnA-elecA*(donnB-donnA)/(elecB-elecA)      
      write(*,*) donn,donnB,elecB
      if (dabs(elecB).lt.delecok) goto 985
      donnA = donnB
      elecA = elecB
      donnB = ddd
      goto 4949
      write(*,*)
      endif      
 985  continue
      write(*,*) 'donn,donnB,elecB = ',donn,donnB,elecB

      sum = 0.d0
      do i = istp1,imitt
      fdes(i) = fdes(i)*dexp(-(donnB-donn)) 
      fdsol(i) = fdsol(i)*dexp((donnB-donn))      
      sum = sum-fdsol(i)+fdes(i)      
      enddo
      celec = surfdens+sum*dz
      write(*,*) 'surfdens,celec (st) = ',surfdens,celec

      jz = imitt+1
      do iz = imitt+1,nfin      
      jz = jz-1
      fdsol(iz) = fdsol(jz)
      fdes(iz) = fdes(jz)      
      enddo
      
      ddmax = -10000.
      niter = 0
 100  continue
      niter = niter+1
c     write(*,*)
      if (mod(niter,100).eq.0) then
      sum = 0.d0
      do i = istp1,imitt         
      sum = sum-fdsol(i)+fdes(i)
      enddo
c     surfdens = -sum*dz
      celec = surfdens+sum*dz      
      write(*,*) 'niter,ddmax,celec =',niter,ddmax,celec
      endif
      if (niter.gt.ioimaxm) goto 200      
      CALL CALCCD
      CALL AVEC
      CALL CALCELAM

      ddmax = 0.d0
      sumfds = 0.d0
      sumfdes = 0.d0
c      rewind 32
c      z = closew-0.5d0*dz      
      do i = istp1,imitt
c      z = z+dz
      sumfds = eblam(i)+sumfds
      sumfdes = eeblam(i)+sumfdes
c      write(32,'(4e12.5)') z,2.d0*c(i,1),dumsum,dumsunm                  
      enddo
c      stop

      stfdes = sumfdes*dz      
      stfds = sumfds*dz
      elecA = stfdes-stfds+surfdens
c      write(*,*) 'elecA,donn = ',elecA,donn
      if (dabs(elecA).lt.delecok) goto 585
      donnA = donn
      if (dabs(donn).lt.0.001d0) then
      ddonn = 0.05d0
      else
      ddonn = frdonn*donn
      endif
      if (elecA.lt.0.d0) then
      donnB = donnA+ddonn
      else
      donnB = donnA-ddonn
      endif
      ddd = donnB
 5959 continue
      donnB = ddd
c      write(*,*) 'donnB,ddd = ',donnB,ddd
      sumfds = 0.d0
      sumfdes = 0.d0
      do i = istp1s,imitt
      sumfdes = eeblam(i)*dexp(-(donnB-donn))+sumfdes      
      sumfds = eblam(i)*dexp((donnB-donn))+sumfds
      enddo
      elecB = (sumfdes-sumfds)*dz+surfdens
      ddd = donnA-elecA*(donnB-donnA)/(elecB-elecA)      
c      write(*,*) donn,donnB,elecB
      if (dabs(elecB).lt.delecok) goto 585
c      write(*,*) 'elecA,elecB = ',elecA,elecB
c      write(*,*) donnA,donnB,ddd
c      iop = iop+1
c      write(*,*) 'donnA,donnB = ',donnA,donnB,ddonn
c      write(*,*) (donnB-donnA)/(elecB-elecA),ddd
c      write(*,*) 
c      if (iop.gt.12) stop      
      donnA = donnB
      elecA = elecB
      donnB = ddd
      goto 5959
 585  continue
      do i = istp1,imitt      
      tfds = eblam(i)*dexp((donnB-donn))
      ddiff = abs(tfds-fdsol(i))/tfds
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdsol(i) = fdsol(i)*dms+tdms*tfds
      tfdes = eeblam(i)*dexp(-(donnB-donn))
      ddiff = abs(tfdes-fdes(i))/tfdes
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdes(i) = fdes(i)*dms+tdms*tfdes
      enddo
      jz = imitt+1
      sum = 0.d0
      do iz = imitt+1,islut      
      jz = jz-1
      fdsol(iz) = fdsol(jz)
      fdes(iz) = fdes(jz)
      sum = sum-fdsol(iz)+fdes(iz)
      enddo
      celec = surfdens+sum*dz
      if (ddmax.lt.ddtol.and.dabs(celec).lt.delecok) goto 200
      donn = donnB
      goto 100
 200  continue

      sum = 0.d0
      do i = istp1,imitt         
      sum = sum-fdsol(i)+fdes(i)
      enddo
      elec = surfdens+sum*dz            
      write(*,*) 'elec,surfdens = ',elec,surfdens

      rewind ife
      appsdens = surfdens
      rewind 90
c      write(90,*) 0.d0,appsdens      
      z = -0.5d0*dz
      do i = 1,nfin                            
      z = z+dz
      write(ife,*) z,fdsol(i),fdes(i)
      appsdens = appsdens+(fdes(i)-fdsol(i))*dz
      write(90,*) z,appsdens      
      enddo
c      z = z+0.5d0*dz
c      appsdens = appsdens+surfdens
c      write(90,*) z,appsdens
      
      avfdm = 0.d0
      avfdnm = 0.d0      
      avfds = 0.d0
      avfdes = 0.d0
      sumFid = 0.d0
      csumFid = 0.d0
      sumuu = 0.d0
      sumuw = 0.d0
      sumuljm = 0.d0
      altsumF = 0.d0
      sumexFreen = 0.d0
      csumexFreen = 0.d0      
      z = closew-0.5d0*dz
c     do iz = istp1s,isluts
      do iz = istp1s,imitt      
      z = z+dz
      fds = fdsol(iz)
      tfdes = fdes(iz)
      avfds = fds+avfds
      avfdes = tfdes+avfdes      
      suu = 0.d0
      suw = 0.d0
      do jz = istp1s,isluts      
      iii = iabs(jz-iz)
      suu = (-fdsol(jz)+fdes(jz))*Phimm(iii)+suu
      enddo
      sumuu = 0.5d0*(-fds+tfdes)*suu+sumuu
      sumuljm = 0.d0
      sumuw = (-fds+tfdes)*csurf+sumuw                  
      
      bslamb = dlog(eblam(iz))-scales
      beslamb = dlog(eeblam(iz))-scalees      

      aabslamb = dlog(eblam(iz))-scales-donn      
      aabeslamb = dlog(eeblam(iz))-scalees+donn
      
      cdt = cdtot(iz)*dm3
      pcdt = pis*cdt
      xsi = (1.d0-pcdt)
      rxsi = 1.d0/xsi
      sqrxsi = rxsi*rxsi
      flog = dlog(xsi)
      aex1 = -(c1+1.d0)*flog-0.5d0*(AA1+BB1*pcdt)*pcdt*sqrxsi
      aex2 = -(c2+1.d0)*flog-0.5d0*(AA2+BB2*pcdt)*pcdt*sqrxsi
      exFreen = (fds+tfdes)*aex1
      sumexFreen = exFreen+sumexFreen
      csumexFreen = (fds+tfdes)*aex1+
     *csumexFreen
      sumFid = sumFid+
     *fds*(dlog(fds)-1.d0-chemps)+tfdes*(dlog(tfdes)-1.d0-chempes)
      csumFid = csumFid+
     *fds*bslamb-fds+tfdes*beslamb-tfdes
      altsumF = altsumF+
     *fds*aabslamb-fds+tfdes*aabeslamb-tfdes
      enddo
      write(*,*) 'z = ',z
      sumexFreen = 2.d0*sumexFreen*dz
      csumexFreen = 2.d0*csumexFreen*dz
      write(*,*) 'sumexFreen, csumexFreen = ',sumexFreen,csumexFreen
      sumF = 2.d0*sumFid*dz+sumexFreen
      csumF = 2.d0*csumFid*dz+sumexFreen
      altsumF = 2.d0*altsumF*dz      
      write(*,*) sumF,csumF
      sumuljm = 2.d0*sumuljm*dz      
      sumuu = 2.d0*sumuu*dz
      sumuw = 2.d0*sumuw*dz
      write(*,*) 'sumuu = ',sumuu
      write(*,*) 'sumuw = ',sumuw
      write(*,*) 'sumuljm = ',sumuljm
      sumF = sumF+sumuu+sumuw
      csumF = csumF+sumuu+sumuw
      altsumF = altsumF-sumuu
      write(*,*) sumF,csumF,altsumF      
      write(*,*) 'sumF = ',sumF
      write(*,*) 'csumF = ',csumF
      write(*,*) 'altsumF = ',altsumF
      omega = sumF
      comega = csumF
      aomega = altsumF
      
      write(*,*) 'grand pot. (excl. surf-surf):  ',omega
      write(*,*) omega
      write(*,*) comega
      write(*,*) 'alt. grand pot. (excl. surf-surf):  ',aomega
      write(*,*) aomega      

      uss = -2.d0*pi*surfdens*surfdens*(h+dm1)*rrT                
      write(*,*) 'u(surf-surf), uss  = ',uss
      totomega = omega+uss
      ctotomega = comega+uss
      atotomega = aomega+uss
      write(*,*)            
      write(*,*) 'total grand pot. (incl. surf-surf):  ',totomega
      write(*,*) totomega+Pbhs*(h+dm1)
      write(*,*) ctotomega+Pbhs*(h+dm1)
      write(*,*) atotomega-2.d0*surfdens*donn+Pbhs*(h+dm1)      
      write(*,*) donn
      write(*,*) surfdens
      write(*,*)      
      write(*,*) 'elec,donn (fin) = ',elec,donn
      write(*,*) 'surfdens = ',surfdens      
      rewind idpot
      write(idpot,*) donn

      fp1S = fdsol(istp1)
      fn1S = fdsol(istp1+2)
      c0Skv = fdsol(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fsm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'solvent contact density - quad. extr.: ',fsm
      fp1S = fdes(istp1)
      fn1S = fdes(istp1+2)
      c0Skv = fdes(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fesm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'co-solvent contact density - quad. extr.: ',fesm      
      Pidc = fsm+fesm
      write(*,*) 
      write(*,*) 'ideal part of contact pressure: '
      write(*,'(1e25.14)') Pidc
      write(*,*) Pidc-Pbhs
      write(*,*) Pidc-Pbhs-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) Pidc-Pbhs-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) 
      
      avfds = avfds*dz/h
      avfdes = avfdes*dz/h      
      write(*,*) 'avfds = ',avfds
      write(*,*) 'avfdes = ',avfdes      

      write(*,*) 'donn = ',donn
      write(*,*) donn
      fdonn = donn
      write(*,*) 'fdonn = donn = ',fdonn
      write(*,*) fdonn
       
c     z = zmitt
      z = (h+dm1)/2.d0      
      du = 0.d0
      zp = closew-0.5d0*dz
      write(*,*) 'zmitt,zp(init) = ',z,zp
      do j = istp1,islut      
      zp = zp+dz
      diffz = dabs(zp-z)
      phic = -2.d0*pi*diffz*rrT*dz      
      du = phic*(-fdsol(j)+fdes(j))+du      
      enddo
      write(*,*) 'zmitt,zp(final) = ',z,zp      
      pmid = du+fdonn+csurf
      write(*,*) 'potential at mid plane: ',pmid,du+donn,du
      write(*,*) pmid
      write(*,*) du+donn+csurf
      write(*,*) du+csurf

      z = 0.d0
      du = 0.d0
      zp = closew-0.5d0*dz      
      write(*,*) 'zlw,zp(init) = ',z,zp
      do j = istp1,islut            
      zp = zp+dz
      diffz = dabs(zp-z)
      phic = -2.d0*pi*diffz*rrT*dz
      du = phic*(-fdsol(j)+fdes(j))+du
      enddo
      write(*,*) 'zlw,zp(final) = ',z,zp            
      psurf = du+fdonn+csurf
      write(*,*) 'potential at surface: ',du+fdonn,du+donn,du
      write(*,*) psurf
      write(*,*) du+donn+csurf
      write(*,*) du+csurf

      write(*,*) 'net surface potential, psurf-pmid: '
      write(*,*) psurf-pmid

      write(*,*) 'grand potential-surfdens*donn:'
      write(*,*) totomega-surfdens*donn

      rewind 92
      z = -0.5d0*dz            
      write(*,*) 'z+0.5d0*dz = ',z+0.5d0*dz
      do i = 1,imitt
      z = z+dz
      zp = -0.5d0*dz
      du = 0.d0
      do j = 1,islut            
      zp = zp+dz
      diffz = dabs(zp-z)
      kkk = iabs(j-i)
      du = Phimm(kkk)*(-fdsol(j)+fdes(j))+du
      enddo
      write(92,*) z,du+csurf
      enddo
      write(*,*) 'ddmax,niter: ',ddmax,niter      
      
      STOP
      END

      subroutine CALCCD
      implicit double precision (a-h,o-z)
      include 'PB.inc'
      fp1S = (fdsol(istp1)+fdes(istp1))
      fn1S=(fdsol(istp1+2)+fdes(istp1+2))
      c0Skv=
     *(fdsol(istp1+1)+fdes(istp1+1))
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fwc = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwc.lt.0.d0) fwc = 0.d0
      z = sclosew-0.5d0*dz
      iii = min(istp1+ism-1,imitt)      
      do iz = istp1s,imitt
      z = z+dz
      zs = closew
      z1 = zs
      z2 = zs+0.5d0*dz
      f1 = fwc
      f2 = (fdsol(istp1)+fdes(istp1))
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = 2.d0*(f2-f1)*((z2-z)**4-(z1-z)**4)
      zs = z2-dz
      do jz = istp1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdsol(jz+1)+fdes(jz+1))
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdsol(iz+ism)+fdes(iz+ism))
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      cdtot(iz) = 0.75d0*sancint*rdm3
      enddo

      do iz = iii+1,imitt      
      z = z+dz
      zs = z-dm1
      z1 = zs
      z2 = zs+dz
      f1 = (fdsol(iz-ism)+fdes(iz-ism))
      f2 =
     *(fdsol(iz-ism+1)+fdes(iz-ism+1))
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do jz = iz-ism+1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdsol(jz+1)+fdes(jz+1))
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdsol(iz+ism)+fdes(iz+ism))
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      cdtot(iz) = 0.75d0*sancint*rdm3
      enddo
      return
      end


      subroutine AVEC
      implicit double precision (a-h,o-z)
      include 'PB.inc'
      do iz = istp1s,imitt
      cdt = cdtot(iz)*dm3
      pcdt = pis*cdt
      xsi = (1.d0-pcdt)
      rxsi = 1.d0/xsi
      sqrxsi = rxsi*rxsi
      flog = dlog(xsi)            
      ae1(iz) = -(c1+1.d0)*flog-0.5d0*(AA1+BB1*pcdt)*pcdt*sqrxsi      
      daex1 = rxsi*(c1+1.d0-0.5d0*AA1*rxsi*(1.d0+2.d0*pcdt*rxsi)-
     *BB1*pcdt*rxsi*(1.d0+pcdt*rxsi))
      convp(iz) = (fdsol(iz)+fdes(iz))*daex1*pis*dhs3
      enddo
      jz = imitt+1
      do iz = imitt+1,imitt+ibl      
      jz = jz-1      
      convp(iz) = convp(jz)
      ae1(iz) = ae1(jz)
      enddo
      return
      end

      subroutine CALCELAM
      implicit double precision (a-h,o-z)
      include 'PB.inc'
c      dimension trams(0:maxel)
      fp1S = convp(istp1s)
      c0Skv = convp(istp1s+1)
      fn1S = convp(istp1s+2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fwc = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      z = closew-0.5d0*dz
      iii = min(istp1s+ism-1,imitt)
      do iz = istp1,iii
      z = z+dz
      zs = sclosew
      z1 = zs
      z2 = zs+0.5d0*dz
      f1 = fwc
      f2 = fp1S
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = 2.d0*(f2-f1)*((z2-z)**4-(z1-z)**4)
      zs = z2-dz
      jjj = iz+ism      
      do jz = istp1s,jjj-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jjj)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      trams = sancint*0.75d0*rdhs3
      
      hz = h-z
      du = 0.d0
      dulj = 0.d0
      do jz = istp1s,isluts
      kkk = iabs(jz-iz)
      du = Phimm(kkk)*(-fdsol(jz)+fdes(jz))+du      
      enddo
      strams = trams+ae1(iz)
      eblam(iz) = dexp(-strams+scales+du+donn+csurf)
      eeblam(iz) = dexp(-strams+scalees-du-donn-csurf)            
      enddo

      do iz = iii+1,imitt
      z = z+dz
      zs = z-dm1
      z1 = zs
      z2 = zs+dz
      f1 = convp(iz-ism)
      f2 = convp(iz-ism+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      jjj = iz+ism      
      do jz = iz-ism+1,jjj-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jjj)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      trams = sancint*0.75d0*rdhs3

      hz = h-z
      du = 0.d0
      dulj = 0.d0
      do jz = istp1s,isluts
      kkk = iabs(jz-iz)
      du = Phimm(kkk)*(-fdsol(jz)+fdes(jz))+du      
      enddo
      strams = trams+ae1(iz)
      eblam(iz) = dexp(-strams+scales+du+donn+csurf)
      eeblam(iz) = dexp(-strams+scalees-du-donn-csurf)            
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
      eblam(iz) = eblam(jz)
      eeblam(iz) = eeblam(jz)      
      enddo
      
      return 
      end

      
