      program platem
      implicit double precision (a-h,o-z)
      include 'dftpol.dshall.inc'      
      dimension c(maxmon,0:maxel),fps(0:maxel),dVex2m(0:maxel),
     *etheta(0:maxel),ch(maxmon,0:maxel),dVexlw(0:maxel)
      write(*,*) 'Monomer 1 constrained to be within'        
      write(*,*) 'one sigma from z=closew, i.e the left wall'  
      write(*,*) '(flexible grafting)'    
      do 1 i = 0,maxel
      fdmon(i) = 0.d0
      fem(i) = 0.d0
      Vexm(i) = 0.d0
      dVex2m(i) = 0.d0
      dVexlw(i) = 0.d0
      fps(i) = 0.d0
      fdgr(i) = 0.d0
      do 2 j = 1,maxmon
      ch(j,i) = 0.d0
 2    c(j,i) = 0.d0
 1    continue
      ifc = 38
      ins = 49
      itau = 50
c      iefil = 54
c      issolv = 56
      isd = 61
      pi = acos(-1.d0)
      bk = 1.38066D-23
      avno = 6.02214D23
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      a1 = 1.d0
      a2 = 2.45696d0
      b1 = 1.d0
      b2 = 4.10386d0
c      c1 = 0.d0
      c1 = -1.d0
      c2 = -3.75503
      AA1 = 2.d0*c1-2.d0*a1-4.d0
      AA2 = 2.d0*c2-2.d0*a2-4.d0
      BB1 = 3.d0-b1+a1-3.d0*c1
      BB2 = 3.d0-b2+a2-3.d0*c2
      Y = (9.82605d0-9.d0*pi*0.25d0)/(9.d0*pi*0.25d0-4.d0*pi/3.d0)
      rphi = fourpi*0.2d0 
      aphi = fourpi*0.5d0
      pis = pi/6.d0
      pit = pi/3.d0
      pif = pi/4.d0
      ddtol = 0.0000001d0
c     CLOSE TO THE WALLS, THE DENSITY IS ASSUMED TO BE ZERO
      open (ifc,file='fcdfil',form='formatted')    
      open (ins,file='input.dftp',form='formatted')
      open (itau  ,file='taufil',form='formatted')
c      open (iefil,file='ekblamfil',form='formatted')
c      open (issolv,file='sigmasolv',form='formatted')
      open (isd,file='surfdens',form='formatted')
      rewind ifc
      rewind ins
      rewind itau
c      rewind iefil
c      rewind issolv
      rewind isd
      read(ins,*) bdm
      read(ins,*) bdstrams
      bds = 0.d0
      read(ins,*) nmon
      read(ins,*) trams
      read(ins,*) ioimaxm
      read(ins,*) ntrams
      read(ins,*) T
      read(ins,*) h
      read(ins,*) dz
      read(ins,*) closew
      read(ins,*) dmm
      read(ins,*) dms
      read(ins,*) kread
      read(ins,*) hwdist
      read(ins,*) awm
      read(ins,*) rwm
      read(ins,*) aws
      read(ins,*) rws
      read(ins,*) ufact
      ufact = 0.d0
      read(ins,*) ess
      read(ins,*) ems
      q1 = 1.d0
      sclosew = closew
c      read(issolv,*) q1
c      read(issolv,*) sclosew
      read(itau,*) tau
      read(itau,*) Af
      read(itau,*) zwc
      read(itau,*) zpmin
      read(isd,*) surfdens
      surfdens = dabs(surfdens)
      write(*,*) 'surfdens = ',surfdens
      q2 = q1*q1
      q3 = q2*q1
      q4 = q2*q2
      q10 = q4*q4*q2
      p1 = (1.d0+q1)*0.5d0
      p2 = p1*p1
      p3 = p2*p1
      p4 = p2*p2
      p10 = p4*p4*p2
      rq3 = 1.d0/q3
      write(*,*) 'q1,p1 = ',q1,p1
c     transformation to sigma(mon) units
      bds = bds*rq3
      tdmm = 1.d0-dmm
      tdms = 1.d0-dms
      rrT = 1.d0/T
      rdz = 1.d0/dz
      twopidz = twopi*dz
      irdz = int(rdz+0.001d0)
      nfack = int(h/dz+0.01d0)
      istart = int(closew/dz+0.01d0)
      istp1 = istart+1 
      islut = int((h-closew)/dz+0.01d0)      
      ists = int(sclosew/dz+0.01d0)
      istp1s = ists+1 
      isluts = int((h-sclosew)/dz+0.01d0)      
      ism = int(1.d0/dz+0.01d0)
      isms = int(q1/dz+0.01d0)
      inw = ism+int(closew/dz+0.01d0)
      inws = isms+int(sclosew/dz+0.01d0)
      pie = pi/8.d0
      dzpie = pie*dz
      rnmon = real(nmon)
      rrnmon = 1.d0/rnmon      
      rrcmon = 1.d0/(rnmon-2.d0)
      checknm = abs(rnmon*0.5d0-real(int(rnmon*0.5d0+0.0001)))
      bdpol = bdm/rnmon
      chempp = -(rnmon-1.d0)*dlog(fourpi)
      scalem = chempp/(2.d0*rnmon)
      emscale = 2.d0*scalem
      imitt = nfack/2
      write(*,*) 'GFD model (old Ym)'   
      write(*,*) 'no. of monomers/polymer = ',nmon
      write(*,*) 'max no. of iterations = ',ioimaxm
      write(*,*) 'temperature = ',T
      write(*,*) 'h,dz = ',h,dz
      write(*,*) 'closew,sclosew = ',closew,sclosew
      write(*,*) 'nfack = ',nfack
      write(*,*) 'istp1,islut,inw = ',istp1,islut,inw
      write(*,*) 'istp1s,isluts,inws = ',istp1s,isluts,inws
      write(*,*) 'ism,isms = ',ism,isms
      write(*,*) 'dmm,dms (density mixing param. mon.,solv.) = ',dmm,dms
      write(*,*) 'TRUNC. + SH. STEEP MORSE WALL POT.!'
      write(*,*)'w(z)=Af*((1.d0-exp(-2*(z-zpmin)))^2-1.)-w(zwc),(z<zwc)'
      write(*,*) 'Af = ',Af
      write(*,*) 'zpmin = ',zpmin
      write(*,*) 'zwc = ',zwc
      refVexm = Af*((1.d0-dexp(-2*(zwc-zpmin)))**2-1.d0)
      write(*,*) 'refVexm = ',refVexm
c      rewind 87
c      hwdist = zpmin-0.5d0*dlog(2.d0-dexp(-2*(zwc-zpmin)))
      write(*,*) 'hwdist = '
      write(*,'(1e25.14)') hwdist
      z = closew-0.5d0*dz
      do iz = istp1s,isluts
      z = z+dz
      zw = z+hwdist
      zh = h-z
      zhw = h+hwdist-z
      Vexm(iz) = 0.d0
      if (zw.lt.zwc) then
      Vexm(iz) = Af*((1.d0-dexp(-2.d0*(zw-zpmin)))**2-1.d0)-refVexm
      dVexlw(iz) = 4.d0*Af*(1.d0-
     *dexp(-2.d0*(zw-zpmin)))*dexp(-2.d0*(zw-zpmin))
      endif
      if (zhw.lt.zwc) then
      Vexm(iz) = Af*((1.d0-dexp(-2.d0*(zhw-zpmin)))**2-1.d0)-refVexm+
     *Vexm(iz)
      dVex2m(iz) = 4.d0*Af*(1.d0-
     *dexp(-2.d0*(zhw-zpmin)))*dexp(-2.d0*(zhw-zpmin))
      endif
c      write(87,'(2f25.14)') z,Vexm(iz)
      enddo
 4948 continue
c      stop

      if (kread.eq.0) then
      z = closew-0.5d0*dz
      do 12 iz = istp1,imitt
      z = z+dz
      if (iz.eq.istp1+ism) write(*,*) z,iz
      if (z.lt.rnmon) then
      fdmon(iz) = bdm*(rnmon-z)
      fem(iz) = fdmon(iz)*rrnmon*(rnmon-z)
      endif
      fdgr(iz) = 0.d0
      if (iz.lt.(istp1+ism)) fdgr(iz) = surfdens
 12   continue
      else
      do iz = istp1,imitt
      if (iz.lt.(istp1+ism)) then
      read(ifc,*) trams,fdmon(iz),fem(iz),fdgr(iz)
      else
      read(ifc,*) trams,fdmon(iz),fem(iz),trams
      fdgr(iz) = 0.d0
      endif
      enddo
      endif
      jz = imitt+1
      do iz = imitt+1,isluts
      jz = jz-1
      fdgr(iz) = fdgr(jz)
      fdmon(iz) = fdmon(jz)
      fem(iz) = fem(jz)
      enddo
      do i = istp1,islut
      etheta(i) = 0.d0
      enddo
      ddmax = 10000.
      niter = 0
 100  continue
      niter = niter+1    
      CALL CDMNEW
c     SOME USEFUL VECTORS ARE OBTAINED IN AVEC
      CALL AVEC
c     EVALUATION OF exp(-0.5*beta*lambda(iz)/2)=ehbclam(iz)
      CALL EBLMNEW
      do i = istp1,istp1+ism-1
      etheta(i) = ebelam(i)
c      etheta(i) = ehbclam(i)*ehbclam(i)
      enddo
      if (mod(niter,100).eq.0) write(*,*) 'ddmax,niter = ',ddmax,niter
      if (niter.gt.ioimaxm) then
      write(*,*) 'NITER.GT.IOIMAXM !',niter
      goto 200
      endif
      do  245 iz = istp1,inw
      sume = 0.d0
      rsume = 0.d0
      do 345 jz = istp1,iz+ism-1
      rsume = etheta(jz)+rsume
 345  sume = ebelam(jz)+sume
      trams = ehbclam(iz)*twopidz
      ch(nmon-1,iz) = (0.5d0*ebelam(iz+ism)+sume)*trams
 245  c(nmon-1,iz) = (0.5d0*etheta(iz+ism)+rsume)*trams
      do 445 iz = inw+1,islut-ism
      sume = 0.5d0*ebelam(iz-ism)
      rsume = 0.5d0*etheta(iz-ism)
      do 545 jz = iz-ism+1,iz+ism-1
      rsume = etheta(jz)+rsume
 545  sume = ebelam(jz)+sume
      trams = ehbclam(iz)*twopidz
      ch(nmon-1,iz) = (0.5d0*ebelam(iz+ism)+sume)*trams
 445  c(nmon-1,iz) = (0.5d0*etheta(iz+ism)+rsume)*trams
      do iz = islut-ism+1,islut
      sume = 0.5d0*ebelam(iz-ism)
      rsume = 0.5d0*etheta(iz-ism)
      do jz = iz-ism+1,islut
      rsume = etheta(jz)+rsume
      sume = ebelam(jz)+sume
      enddo
      ch(nmon-1,iz) = ehbclam(iz)*sume*twopidz
      c(nmon-1,iz) = etheta(iz)*rsume*twopidz
      enddo

      k = nmon-1
      do 745 mmm = 2,nmon-2 
      k = k-1
      do  845 iz = istp1,inw
      sume = 0.d0
      do 945 jz = istp1,iz+ism-1
 945  sume = ehbclam(jz)*c(k+1,jz)+sume
      c(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*c(k+1,iz+ism)+
     *sume)*twopidz
      sume = 0.d0
      do jz = istp1,iz+ism-1
      sume = ehbclam(jz)*ch(k+1,jz)+sume
      enddo
 845  ch(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*ch(k+1,iz+ism)+
     *sume)*twopidz
      do 1045 iz = inw+1,islut-ism
      sume = 0.5d0*ehbclam(iz-ism)*c(k+1,iz-ism)
      do 1145 jz = iz-ism+1,iz+ism-1
 1145 sume = ehbclam(jz)*c(k+1,jz)+sume
      c(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*c(k+1,iz+ism)+
     *sume)*twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(k+1,iz-ism)
      do jz = iz-ism+1,iz+ism-1
      sume = ehbclam(jz)*ch(k+1,jz)+sume
      enddo
 1045 ch(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*ch(k+1,iz+ism)+
     *sume)*twopidz
      do iz = islut-ism+1,islut
      sume = 0.5d0*ehbclam(iz-ism)*c(k+1,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*c(k+1,jz)+sume
      enddo
      c(k,iz) = ehbclam(iz)*sume*twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(k+1,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*ch(k+1,jz)+sume
      enddo
      ch(k,iz) = ehbclam(iz)*sume*twopidz
      enddo
 745  continue

      do  1245 iz = istp1,inw
      sume = 0.d0
      do 1345 jz = istp1,iz+ism-1
 1345 sume = ehbclam(jz)*c(2,jz)+sume
      c(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*c(2,iz+ism)+sume)*
     *twopidz
      sume = 0.d0
      do 8645 jz = istp1,iz+ism-1
 8645 sume = ehbclam(jz)*ch(2,jz)+sume
 1245 ch(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*ch(2,iz+ism)+sume)*
     *twopidz
c 1245 ch(1,iz) = 
c     *ehbclam(iz)*ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*ch(2,iz+ism)+sume)*
c     *twopidz
      do 1445 iz = inw+1,islut-ism
      sume = 0.5d0*ehbclam(iz-ism)*c(2,iz-ism)
      do 1545 jz = iz-ism+1,iz+ism-1
 1545 sume = ehbclam(jz)*c(2,jz)+sume
      c(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*c(2,iz+ism)+sume)*
     *twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(2,iz-ism)
      do 8745 jz = iz-ism+1,iz+ism-1
 8745 sume = ehbclam(jz)*ch(2,jz)+sume
 1445 ch(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*ch(2,iz+ism)+sume)*
     *twopidz
c 1445 ch(1,iz) = 
c     *ehbclam(iz)*ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*ch(2,iz+ism)+sume)*
c     *twopidz
      do iz = islut-ism+1,islut
      sume = 0.5d0*ehbclam(iz-ism)*c(2,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*c(2,jz)+sume
      enddo
      c(1,iz) = ebelam(iz)*sume*twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(2,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*ch(2,jz)+sume
      enddo
      ch(1,iz) = ebelam(iz)*sume*twopidz
c      ch(1,iz) = ehbclam(iz)*ehbclam(iz)*sume*twopidz
      enddo

      sumgrmon = 0.d0
      do i = istp1s,istp1s+ism-1
      sumgrmon = ch(1,i)+sumgrmon
      enddo
      sumgrmon = sumgrmon*dz
      fnorm = surfdens/sumgrmon     
       do i = istp1s,istp1s+ism-1
      fdgr(i) = ch(1,i)*fnorm
      enddo          

      if (ddmax.lt.ddtol) goto 200
      ddmax = 0.d0
      z = -0.5d0*dz
      do 9 i = istp1s,islut
      z = z+dz
      if (z.lt.rnmon) then
      dumsum = 0.d0 
      do 10 k = 2,nmon-1
 10   dumsum = c(k,i)*ch(nmon+1-k,i)+dumsum
      ttfem(i) = c(1,i)*fnorm
c      ttfem(i) = 2.d0*c(1,i)*fnorm
      dumsum = dumsum*fnorm
      if (i.lt.(istp1+ism)) then
      tfdmon(i) = dumsum+ttfem(i)+fdgr(i)
      else
      fdgr(i) = 0.d0
      tfdmon(i) = dumsum+ttfem(i)
      endif
      else
      tfdmon(i) = 0.d0
      ttfem(i) = 0.d0
      endif
 9    continue
      z = -0.5d0*dz
      k = islut+1
      do i = istp1,imitt
      k = k-1
      tfdm = tfdmon(i)+tfdmon(k)
      tfem = ttfem(i)+ttfem(k)
      z = z+dz
      if (z.lt.rnmon) then      
      ddiff = abs(tfdm-fdmon(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fem(i) = fem(i)*dmm+tdmm*tfem
      fdmon(i) = fdmon(i)*dmm+tdmm*tfdm
      else
      fdmon(i) = 0.d0
      fem(i) = 0.d0
      endif
      enddo
      jz = imitt+1
      do iz = imitt+1,isluts
      jz = jz-1
      fdgr(iz) = fdgr(jz)
      fdmon(iz) = fdmon(jz)
      fem(iz) = fem(jz)
      enddo
      goto 100
 200  continue

      write(*,*)
      write(*,*) 'fnorm = ',fnorm
      niz = 0
      sumfdm = 0.d0
      Freen = 0.d0
      z = 0.d0
      do 500 iz = istp1,imitt
      bclamb = 2.d0*(dlog(ehbclam(iz))-scalem)
      belamb = dlog(ebelam(iz))-emscale
      niz = niz+1
      fdm = fdmon(iz)   
      fde = fem(iz)+fdgr(iz)
      fdc = fdm-fde
      sumfdm = fdm+sumfdm
      cdt = cdmon(iz)
      pcdt = pis*cdt
      xsi = (1.d0-pcdt)
      rxsi = 1.d0/xsi
      sqrxsi = rxsi*rxsi
      flog = dlog(xsi)
      aex1 = -(c1+1.d0)*flog-0.5d0*(AA1+BB1*pcdt)*pcdt*sqrxsi
      aex2 = -(c2+1.d0)*flog-0.5d0*(AA2+BB2*pcdt)*pcdt*sqrxsi
      exFreen = fdc*Y*(aex2-aex1)+0.5d0*fde*aex2
 500  Freen = Freen+fdc*bclamb+fde*belamb+fdm*(Vexm(iz)-rrnmon)+
     *fdgr(iz)*dlog(fnorm)+exFreen

      Freen = 2.d0*Freen*dz
      write(*,*) '/beta*(grand pot.) = ',Freen
      write(*,*) '(grand pot.)/epsmm = ',Freen*T      
      rewind ifc
      Fs2m = 0.d0
      Fslw = 0.d0
      Pmacross = 0.d0
      avfdm = 0.d0
      avfem = 0.d0
      z = sclosew-0.5d0*dz
      do 61 iz = istp1s,isluts
      z = z+dz
      Fs2m = fdmon(iz)*dVex2m(iz)+Fs2m
      Fslw = fdmon(iz)*dVexlw(iz)+Fslw
      if (iz.lt.(imitt+1)) then
      Pmacross = fdmon(iz)*dVex2m(iz)+Pmacross
      endif
      write(ifc,'(1f12.5,3e25.14)') z,fdmon(iz),fem(iz),fdgr(iz)
      avfdm = fdmon(iz)+avfdm
 61   avfem = fem(iz)+avfem
      write(*,*)
      write(*,*) '/beta*(Freen+Pb*h) (Pb = 0) = '
      write(*,*) Freen
      write(*,*)
      fp1S = fdmon(isluts)
      fn1S = fdmon(isluts-2)
      c0Skv = fdmon(isluts-1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fwcmq = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwcmq.lt.0.d0) fwcmq = 0.d0
      fwcml = fp1S+0.5d0*(fp1S-c0Skv)
      if (fwcml.lt.0.d0) fwcmq = 0.d0
      Fs2m = fwcmq-Fs2m*dz

      fp1S = fdmon(istp1)
      fn1S = fdmon(istp1+2)
      c0Skv = fdmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fwcmq = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwcmq.lt.0.d0) fwcmq = 0.d0
      fwcml = fp1S+0.5d0*(fp1S-c0Skv)
      if (fwcml.lt.0.d0) fwcmq = 0.d0
      Fslw = fwcmq-Fslw*dz

      fp1S = fdmon(istp1+ism)
      fn1S = fdmon(istp1+ism+2)
      c0Skv = fdmon(istp1+ism+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fmdp = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'nm(ism+delta) - quad. extr. = ',fmdp
      fp1S = fdmon(istp1+ism-1)
      fn1S = fdmon(istp1+ism-3)
      c0Skv = fdmon(istp1+ism-2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fmdm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'nm(izslab-delta) - quad. extr. = ',fmdm

      fp1S = fdgr(istp1+ism-1)
      fn1S = fdgr(istp1+ism-3)
      c0Skv = fdgr(istp1+ism-2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fgrbr = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'nmgraft(izslab-delta) - quad. extr. = ',fgrbr

      write(*,*) 'correcting for graft interactions'      

      Pdonn = fmdp-fmdm
      write(*,*) 'Pdonn = ',Pdonn
      Pnlw = Fslw+Pdonn

      Pnet = Fs2m+Pdonn
      write(*,*) 
      write(*,*)'/beta*solvation force (right wall) = '
      write(*,*)  Fs2m-fgrbr
      write(*,*)  Pnet
      write(*,*)
      write(*,*)'/beta*solvation force (left wall) = '
      write(*,*)  Fslw-fgrbr
      write(*,*)  Pnlw
      write(*,*)
      Pmacross = -2.d0*Pmacross*dz
      write(*,*) 'monomer contact density (left wall) = ',fwcmq
      avfdm = avfdm/real((islut+1-istp1))
      avfem = avfem/real((islut+1-istp1))
      write(*,*) 'avfdm = ',avfdm
      fp1S = fdmon(imitt+1)
      fn1S = fdmon(imitt+3)
      c0Skv = fdmon(imitt+2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fmm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'monomer mid plane density - quad. extr. = ',fmm
      Pidmid = fmm
      write(*,*) 'ideal part of mid plane pressure:'
      write(*,'(1e25.14)') Pidmid
      write(*,*) 'Pmacross:'
      write(*,'(1e25.14)') Pmacross
      write(*,*) 
 9998 continue
      write(*,*) 'ddmax,niter = ',ddmax,niter
 9999 continue
      STOP
      END


      subroutine CDMNEW
      implicit double precision (a-h,o-z)
      include 'dftpol.dshall.inc'
      fp1S = fdmon(istp1)
      fn1S = fdmon(istp1+2)
      c0Skv = fdmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fwc = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwc.lt.0.d0) fwc = 0.d0
      z = sclosew-0.5d0*dz
      do 33 iz = istp1s,istp1+ism-1
      z = z+dz
      zs = closew
      z1 = zs
      z2 = zs+0.5d0*dz
      f1 = fwc
      f2 = fdmon(istp1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
c      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      sanci4 = 2.d0*(f2-f1)*((z2-z)**4-(z1-z)**4)
      zs = z2-dz
      do 34 jz = istp1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 34   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
 33   cdmon(iz) = 0.75d0*sancint

      do 36 iz = istp1+ism,imitt
      z = z+dz
      zs = z-1.d0
      z1 = zs
      z2 = zs+dz
      f1 = fdmon(iz-ism)
      f2 = fdmon(iz-ism+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-ism+1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
 36   cdmon(iz) = 0.75d0*sancint
      return
      end



      subroutine AVEC
      implicit double precision (a-h,o-z)
c      include 'dftpol.inc'
      include 'dftpol.dshall.inc'
      do 72 iz = istp1,imitt
      cdt = cdmon(iz)
      pcdt = pis*cdt
      xsi = (1.d0-pcdt)
      rxsi = 1.d0/xsi
      sqrxsi = rxsi*rxsi
      flog = dlog(xsi)      
      ae1(iz) = -(c1+1.d0)*flog-0.5d0*(AA1+BB1*pcdt)*pcdt*sqrxsi
      ae2(iz) = -(c2+1.d0)*flog-0.5d0*(AA2+BB2*pcdt)*pcdt*sqrxsi
      daex1 = rxsi*(c1+1.d0-0.5d0*AA1*rxsi*(1.d0+2.d0*pcdt*rxsi)-
     *BB1*pcdt*rxsi*(1.d0+pcdt*rxsi))
      daex2 = rxsi*(c2+1.d0-0.5d0*AA2*rxsi*(1.d0+2.d0*pcdt*rxsi)-
     *BB2*pcdt*rxsi*(1.d0+pcdt*rxsi))
cc 72   convp(iz) = (Y*(fdmon(iz)-fem(iz)-fdgr(iz))*(daex2-daex1)+
cc     *0.5d0*(fem(iz)+fdgr(iz))*daex2)*pis
c 72   convp(iz) = (Y*(fdmon(iz)-fem(iz))*(daex2-daex1)+
c     *fem(iz)*daex2)*pis
 72   convp(iz) = (Y*(fdmon(iz)-fem(iz)-fdgr(iz))*(daex2-daex1)+
     *0.5d0*(fem(iz)+fdgr(iz))*daex2)*pis
      jz = imitt+1
      do 502 iz = imitt+1,islut
      jz = jz-1
      ae1(iz) = ae1(jz)
      ae2(iz) = ae2(jz)
 502  convp(iz) = convp(jz)
      return
      end

      subroutine EBLMNEW
      implicit double precision (a-h,o-z)
c      include 'dftpol.inc'
      include 'dftpol.dshall.inc'
      fp1S = convp(istp1)
      c0Skv = convp(istp1+1)
      fn1S = convp(istp1+2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fwc = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      z = closew-0.5d0*dz
      do 33 iz = istp1,istp1+ism-1
      z = z+dz
      zs = closew
      z1 = zs
      z2 = zs+0.5d0*dz
      f1 = fwc
      f2 = fp1S
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = 2.d0*(f2-f1)*((z2-z)**4-(z1-z)**4)
      zs = z2-dz
      do 34 jz = istp1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      zp = z1
      f1 = f2
      f2 = convp(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 34   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      zp = z1
      f1 = f2
      f2 = convp(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      trams = sancint*0.75d0
      emtrams = trams+0.5d0*ae2(iz)
c      emtrams = trams+ae2(iz)
      cmtrams = trams+Y*(ae2(iz)-ae1(iz))
      ebelam(iz) = dexp(-(emtrams+Vexm(iz))+emscale)
 33   ehbclam(iz) = dexp(-0.5d0*(cmtrams+Vexm(iz))+scalem)

      do 36 iz = istp1+ism,imitt
      z = z+dz
      zs = z-1.d0
      z1 = zs
      z2 = zs+dz
      zp = zs
      f1 = convp(iz-ism)
      f2 = convp(iz-ism+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-ism+1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      zp = zs
      f1 = f2
      f2 = convp(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      zp = zs
      f1 = f2
      f2 = convp(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      trams = sancint*0.75d0
      emtrams = trams+0.5d0*ae2(iz)
c      emtrams = trams+ae2(iz)
      cmtrams = trams+Y*(ae2(iz)-ae1(iz))
      ebelam(iz) = dexp(-(emtrams+Vexm(iz))+emscale)
 36   ehbclam(iz) = dexp(-0.5d0*(cmtrams+Vexm(iz))+scalem)
      jz = imitt+1
      do 502 iz = imitt+1,islut
      jz = jz-1
      ebelam(iz) = ebelam(jz)
 502  ehbclam(iz) = ehbclam(jz)
      return 
      end

