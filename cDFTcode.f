      program platem
      implicit double precision (a-h,o-z)
      include 'PB.inc'
      dimension c(0:maxel,maxmon),
     *chA(0:maxel),chB(0:maxel),rehbclam(0:maxel),
     *rehbclanm(0:maxel),
     *ttfenm(0:maxel),ttfem(0:maxel),ttfdmon(0:maxel),
     *ttfdnmon(0:maxel)
c     ALTERNATING CHARGES!
      ifc = 38
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
      delecok = 1.D-11
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
      open (ifc,file='fcdfil',form='formatted')
      open (ifnc,file='fncdfil',form='formatted')      
      open (ins,file='dpcinp',form='formatted')
      open (idpot,file='donnan',form='formatted')
      open (ksurf,file='surfdens',form='formatted')      
      open (ikh,file='separation',form='formatted')                  
      rewind idpot      
      rewind ifc
      rewind ifnc      
      rewind ins
      rewind ksurf
      rewind ikh
      read(ins,*) 
      read(ins,*) bdm,bdestrams
      read(ins,*) 
      read(ins,*) nval,nsval
      nval = 1
      nsval = 1
      read(ins,*) 
      read(ins,*) nmon
c      if (mod(nmon,2).eq.0) then
c      abspolval = 0.d0
c      else
c      abspolval = 1.d0
c      endif
      if (mod(nmon,2).eq.0) then
      write(*,*) 'nmon must be odd',nmon
      stop   
      else
      abspolval = 1.d0
      endif      
      write(*,*) 'abs(polval) = ',abspolval
      read(idpot,*) donn
      frdonn = 0.01d0
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
      
c     read(ins,*) dm1,ds1,bl
      read(ins,*) dm1,dmw,bl
      ds1 = dm1
      
      write(*,*) 'dm1,dmw = ',dm1,dmw
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
      bdnm = bdm
      write(*,*) 'surdfdens = ',surfdens 
      write(*,*) 'bdm = ',bdm
      write(*,*) 'bdnm = ',bdnm
      tdmm = 1.d0-dmm
      tdms = 1.d0-dms
c      rrT = 1.d0/T
      rdz = 1.d0/dz
c     twopidz = twopi*dz
      twopidz = 0.5d0*dz*rbl      
      irdz = int(rdz+0.001d0)
      ibl = int(bl/dz+0.01d0)      

c     closew = 0.5d0*dm1
      closew = 0.5d0*dmw
      write(*,*) 'closew = 0.5*dmw = ',closew
c      nfack = int((h+dm1)/dz+0.01d0)
c      imitt = int((h+dm1)/dz+0.01d0)/2
      nfack = int((h+dmw)/dz+0.01d0)
      imitt = int((h+dmw)/dz+0.01d0)/2            
      
      sclosew = closew      
      istart = int(closew/dz+0.01d0)
      iblnw = ibl+int(closew/dz+0.01d0)

c     csurf = -2.d0*pi*surfdens*rrT*(h+dm1)
      csurf = -2.d0*pi*surfdens*rrT*(h+dmw)                
      write(*,*) 'csurf = ',csurf            

      nfin = 2*imitt
      istp1 = istart+1 
cc     islut = nfack
      
cc     islut = int((h+0.5d0*dm1)/dz+0.01d0)
      islut = int((h+0.5d0*dmw)/dz+0.01d0)
c      islut = istart+int(h/dz+0.01d0)      
      
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

      rnmon = dfloat(nmon)
      rrnmon = 1.d0/rnmon
      Yfact = (rnmon-2.d0)*Y
      bdpol = bdm/rnmon
      nhalfmon = nmon/2      
      write(*,*) 'nmon = ',nmon
      write(*,*) 'bdm,bdpol = ',bdm,bdpol
      
      es22ss = 0.d0
      es22mm = es22ss
      es22ms = es22ss
      Pblj = 0.5d0*es22mm*bdm**2

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

      bdt = 2.d0*bdm*dhs3
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
      
      rNpda1 = pis*bdm*rxsib*((c1+1.d0)-
     *0.5d0*AA1*rxsib*(1.d0+2.d0*pis*bdt*rxsib)-
     *BB1*pis*rxsib*bdt*(1.d0+pis*bdt*rxsib))
      rNpda2 = pis*bdm*rxsib*((c2+1.d0)-
     *0.5d0*AA2*rxsib*(1.d0+2.d0*pis*bdt*rxsib)-
     *BB2*pis*rxsib*bdt*(1.d0+pis*bdt*rxsib))

      baex1 = aex1
      baex2 = aex2

      exchpp = Yfact*(aex2-aex1)+aex2+
     *(2.d0*Yfact*(rNpda2-rNpda1)+2.d0*rNpda2)*dhs3
      exPp = 2.d0*bdpol*(Yfact*(-Vdae2dV+Vdae1dV)-Vdae2dV)
      exP = exPp

      chempp = dlog(bdpol)+exchpp
c     note bdm = bdnm => chempp = chempnp
      Pbhs = 2.d0*bdpol+exPp
      Pb = Pbhs

      scalem = chempp/(2.d0*rnmon)
      emscale = 2.d0*scalem
      write(*,*) 'bdm = bdnm = ',bdm,bdnm
      write(*,*) 'max no. of iterations = ',ioimaxm
      write(*,*) 'bjerrum,rbjerrum = ',bjerrum, rbjerrum
      write(*,*) 'h,dz = ',h,dz
      write(*,*) 'polyion chemical potential = ',chempp
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
      
      if (kread.eq.0) then
      donn = 0.d0
      do iz = istp1,imitt       
      fem(iz) = 2.d0*bdm/rnmon
      fdmon(iz) = bdm
      fenm(iz) = 2.d0*bdm/rnmon
      fdnmon(iz) = bdm      
      enddo
      else
c     do iz = istp1,imitt
      do iz = 1,imitt                
      read(ifc,*) trams,fdmon(iz),fem(iz)
      read(ifnc,*) trams,fdnmon(iz),fenm(iz)      
      enddo
      endif

      sum = 0.d0
      sumn = 0.d0
      do i = istp1,imitt
      sum = sum+fdmon(i)-fdnmon(i)
      sumn = sumn+fdnmon(i)
      enddo
      sum = sum*dz
      sumn = sumn*dz
      elec = surfdens+sum
      write(*,*) 'surfdens,elec (st) = ',surfdens,elec          

      if (dabs(elec).gt.delecok) then
      useful = 0.5d0*surfdens/sumn
      x1 = useful+(useful**2+(sum+sumn)/sumn)**0.5d0
c      x2 = useful-(useful**2+(sum+sumn)/sumn)**0.5d0         
c      write(*,*) x1,x2
      ddonn = dlog(x1)
      sumfdm = 0.d0
      sumfdnm = 0.d0
      do iz = istp1s,imitt
      fdmon(iz) = fdmon(iz)*dexp(-abspolval*ddonn)
      fdnmon(iz) = fdnmon(iz)*dexp(abspolval*ddonn)      
      sumfdm = sumfdm+fdmon(iz)
      sumfdnm = sumfdnm+fdnmon(iz)
      enddo
      stfdm = sumfdm*dz
      stfdnm = sumfdnm*dz      
      elecB = stfdm-stfdnm+surfdens
      donn = donn+ddonn
      write(*,*) 'ddonn = ',ddonn
      write(*,*) 'donn = ',donn,elecB
      endif
      sum = 0.d0
      do i = istp1,imitt
      sum = sum+fdmon(i)-fdnmon(i)
      enddo
      celec = surfdens+sum*dz
      write(*,*) 'surfdens,celec (st) = ',surfdens,celec

      do i = 1,istp1
      fem(i) = 0.d0
      fdmon(i) = 0.d0
      fenm(i) = 0.d0
      fdnmon(i) = 0.d0      
      enddo
      
      jz = imitt+1
      do iz = imitt+1,nfin      
      jz = jz-1
      fdmon(iz) = fdmon(jz)      
      fem(iz) = fem(jz)
      fdnmon(iz) = fdnmon(jz)      
      fenm(iz) = fenm(jz)            
      enddo
      
      ddmax = -10000.
      niter = 0
 100  continue
      niter = niter+1
c     write(*,*)
      if (mod(niter,100).eq.0) then
      sum = 0.d0
      do i = istp1,imitt         
      sum = sum+fdmon(i)-fdnmon(i)
      enddo
c     surfdens = -sum*dz
      celec = surfdens+sum*dz      
      write(*,*) 'niter,ddmax,celec =',niter,ddmax,celec
      endif
      if (niter.gt.ioimaxm) goto 200      
      CALL CALCCD
      CALL AVEC
      CALL CALCELAM

      ddmax = -10000.
      do i = imitt-ibl,islut
      rehbclam(i) = 1.d0/ehbclam(i)
      rehbclanm(i) = 1.d0/ehbclanm(i)
      enddo
      
c     polymers with negative ends:
      nAB = 1
      do  iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = ebelanm(jz)+sume
      enddo
      tuu = ehbclam(iz)*(0.5d0*ebelanm(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
      c(iz,nmon-1) = tuu
      enddo
  
      do iz = iblnw+1,imitt
      sume = 0.5d0*ebelanm(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = ebelanm(jz)+sume
      enddo
      tuu = ehbclam(iz)*(0.5d0*ebelanm(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
      c(iz,nmon-1) = tuu
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
      c(iz,nmon-1) = tuu*rehbclam(jz)
      enddo

      k = nmon-1
      do mmm = 2,nmon-2
      k = k-1
      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do  iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      tuu = ehbclanm(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      chB(iz) = tuu*ehbclanm(iz)
      c(iz,k) = tuu
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chA(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      tuu = ehbclanm(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      chB(iz) = tuu*ehbclanm(iz)
      c(iz,k) = tuu
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
      tuu = chB(jz)
      chB(iz) = tuu
      c(iz,k) = tuu*rehbclanm(jz)
      enddo
      
      else
      do iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      tuu = ehbclam(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
      c(iz,k) = tuu
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chB(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      tuu = ehbclam(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
      c(iz,k) = tuu
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
      c(iz,k) = tuu*rehbclam(jz)
      enddo
      endif

      enddo

      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      c(iz,1) = ebelanm(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chA(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      c(iz,1) = ebelanm(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      enddo
      else
      do iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      c(iz,1) = ebelam(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chB(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      c(iz,1) = ebelam(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      enddo
      endif

      ddmax = 0.d0
      do i = istp1,imitt
      dumsum = 0.d0
      dumsunm = 0.d0
      do k = 2,nmon-1
      if (mod(k,2).eq.0) then
      dumsum = c(i,k)*c(i,nmon+1-k)+dumsum
      else
      dumsunm = c(i,k)*c(i,nmon+1-k)+dumsunm
      endif   
      enddo
      ttfenm(i) = 2.d0*c(i,1)
      ttfdnmon(i) = dumsunm+ttfenm(i)
      ttfdmon(i) = dumsum
      enddo

c     polymers with positive ends:
      nAB = 1
      do  iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = ebelam(jz)+sume
      enddo
      tuu = ehbclanm(iz)*(0.5d0*ebelam(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclanm(iz)
      c(iz,nmon-1) = tuu
      enddo
  
      do iz = iblnw+1,imitt
      sume = 0.5d0*ebelam(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = ebelam(jz)+sume
      enddo
      tuu = ehbclanm(iz)*(0.5d0*ebelam(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclanm(iz)
      c(iz,nmon-1) = tuu
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
c     c(iz,nmon-1) = tuu*rehbclam(jz)
      c(iz,nmon-1) = tuu*rehbclanm(jz)      
      enddo

      k = nmon-1
      do mmm = 2,nmon-2
      k = k-1
      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do  iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      tuu = ehbclam(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      chB(iz) = tuu*ehbclam(iz)
      c(iz,k) = tuu
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chA(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      tuu = ehbclam(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      chB(iz) = tuu*ehbclam(iz)
      c(iz,k) = tuu
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
      tuu = chB(jz)
      chB(iz) = tuu
      c(iz,k) = tuu*rehbclam(jz)
      enddo
      
      else
      do iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      tuu = ehbclanm(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclanm(iz)
      c(iz,k) = tuu
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chB(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      tuu = ehbclanm(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      chA(iz) = tuu*ehbclanm(iz)
      c(iz,k) = tuu
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
      c(iz,k) = tuu*rehbclanm(jz)
      enddo
      endif

      enddo

      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      c(iz,1) = ebelam(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chA(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chA(jz)+sume
      enddo
      c(iz,1) = ebelam(iz)*(0.5d0*chA(iz+ibl)+sume)*twopidz
      enddo
      else
      do iz = istp1,iblnw
      sume = 0.d0
      do jz = istp1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      c(iz,1) = ebelanm(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      enddo
      do iz = iblnw+1,imitt
      sume = 0.5d0*chB(iz-ibl)
      do jz = iz-ibl+1,iz+ibl-1
      sume = chB(jz)+sume
      enddo
      c(iz,1) = ebelanm(iz)*(0.5d0*chB(iz+ibl)+sume)*twopidz
      enddo
      endif

      ddmax = 0.d0
      sumfdm = 0.d0
      sumfdnm = 0.d0
c      rewind 32
c      z = closew-0.5d0*dz
      do i = istp1,imitt
c      z = z+dz
      dumsum = 0.d0
      dumsunm = 0.d0
      do k = 2,nmon-1
      if (mod(k,2).eq.0) then
      dumsunm = c(i,k)*c(i,nmon+1-k)+dumsunm
      else
      dumsum = c(i,k)*c(i,nmon+1-k)+dumsum
      endif   
      enddo
      ttfem(i) = 2.d0*c(i,1)
      ttfdmon(i) = dumsum+ttfem(i)+ttfdmon(i)
      ttfdnmon(i) = dumsunm+ttfdnmon(i)
      sumfdm = ttfdmon(i)+sumfdm
      sumfdnm = ttfdnmon(i)+sumfdnm
c      write(32,'(4e12.5)') z,2.d0*c(i,1),dumsum,dumsunm            
      enddo
c      stop
      
      sump = sumfdm*dz
      sumn = sumfdnm*dz
      sum = sump-sumn
      elec = sum+surfdens
      if (dabs(elec).gt.delecok) then
      useful = 0.5d0*surfdens/sumn
      x1 = useful+(useful**2+(sum+sumn)/sumn)**0.5d0
      ddonn = dlog(x1)
      else
      ddonn = 0.d0
      endif   

      sumfdm = 0.d0
      sumfdnm = 0.d0      
      do i = istp1s,imitt
      sumfdm = ttfdmon(i)*dexp(-abspolval*ddonn)+sumfdm
      sumfdnm = ttfdnmon(i)*dexp(abspolval*ddonn)+sumfdnm      
      enddo
      elec = (sumfdm-sumfdnm)*dz+surfdens
c      write(*,*) elec,ddonn
c      stop
            
      do i = istp1,imitt      
      tfdm = ttfdmon(i)*dexp(-abspolval*ddonn)
      ddiff = abs(tfdm-fdmon(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdmon(i) = fdmon(i)*dmm+tdmm*tfdm

      tfdm = ttfdnmon(i)*dexp(abspolval*ddonn)
      ddiff = abs(tfdm-fdnmon(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdnmon(i) = fdnmon(i)*dmm+tdmm*tfdm      

      tfdm = ttfem(i)*dexp(-abspolval*ddonn)
      ddiff = abs(tfdm-fem(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fem(i) = fem(i)*dmm+tdmm*tfdm

      tfdm = ttfenm(i)*dexp(abspolval*ddonn)
      ddiff = abs(tfdm-fenm(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fenm(i) = fenm(i)*dmm+tdmm*tfdm                  
      enddo

      jz = imitt+1
      sum = 0.d0
      do iz = imitt+1,islut      
      jz = jz-1
      fdmon(iz) = fdmon(jz)
      fem(iz) = fem(jz)
      fdnmon(iz) = fdnmon(jz)
      fenm(iz) = fenm(jz)            
      sum = sum+fdmon(iz)-fdnmon(iz)
      enddo
      celec = surfdens+sum*dz
      donn = donn+ddonn                  
      if (ddmax.lt.ddtol.and.dabs(celec).lt.delecok) goto 200
      goto 100
 200  continue

      sum = 0.d0
      do i = istp1,imitt         
      sum = sum+fdmon(i)-fdnmon(i)
      enddo
      elec = surfdens+sum*dz
      write(*,*) 'elec,surfdens = ',elec,surfdens

      write(*,*) 'bond consistency check (pos. chains)'
      
      ddmax = 0.d0
      sumposmon_pospol = 0.d0
      sumnegmon_pospol = 0.d0
      sumendmon_pospol = 0.d0
      sumposcmon = 0.d0
      sumfem = 0.d0
      do i = istp1,imitt
      dumsum = 0.d0
      dumsunm = 0.d0
      do k = 2,nmon-1
      if (mod(k,2).eq.0) then
      dumsunm = c(i,k)*c(i,nmon+1-k)+dumsunm
      else
      dumsum = c(i,k)*c(i,nmon+1-k)+dumsum
      endif   
      enddo
      endmon = 2.d0*c(i,1)
      sumposmon_pospol = endmon+dumsum+sumposmon_pospol
      sumposcmon = dumsum+sumposcmon
      sumnegmon_pospol = dumsunm+sumnegmon_pospol
      sumendmon_pospol = endmon+sumendmon_pospol
      sumfem = fem(i)+sumfem
      enddo
      write(*,*) sumposmon_pospol,sumnegmon_pospol
      write(*,*) sumposcmon,sumendmon_pospol,sumposcmon+sumendmon_pospol
      write(*,*) sumendmon_pospol,sumfem
      write(*,*) 'sumposmon_pospol/sumnegmon_pospol = ',
     *sumposmon_pospol/sumnegmon_pospol
      write(*,*) '(rnmon+1.)/(rnmon-1.) = ',(rnmon+1.)/(rnmon-1.)
      write(*,*) 'sumfem/sumendmon_pospol = ',sumfem/sumendmon_pospol
      
      rewind ifc
      rewind ifnc
      sumfem = 0.d0
      sumfenm = 0.d0
      count = 0.d0
      appsdens = surfdens
      rewind 90
c      write(90,*) 0.d0,appsdens
      z = -0.5d0*dz
c     do i = 1,imitt
      do i = 1,nfin                            
      z = z+dz
      write(ifc,*) z,fdmon(i),fem(i)
      write(ifnc,*) z,fdnmon(i),fenm(i)
      sumfem = sumfem+fem(i)
      sumfenm = sumfenm+fenm(i)
      count = count+1.
      appsdens = appsdens+(fdmon(i)-fdnmon(i))*dz
      write(90,*) z,appsdens
      enddo
c      z = z+0.5d0*dz
c      appsdens = appsdens+surfdens
c      write(90,*) z,appsdens            
      write(*,*) 'avfem,avfenm = ',sumfem/count,sumfenm/count
      write(*,*) 'polymer conc. ratio, avfem/avfenm: ',sumfem/sumfenm
      
      avfdm = 0.d0
      avfdnm = 0.d0      
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
      fdm = fdmon(iz)
      fdnm = fdnmon(iz)      
      fds = fdsol(iz)
      tfdes = fdes(iz)
      tfem = fem(iz)
      tfenm = fenm(iz)      
      avfdm = fdm+avfdm
      avfdnm = fdnm+avfdnm      
      suu = 0.d0
      suw = 0.d0
      do jz = istp1s,isluts      
      iii = iabs(jz-iz)
      suu = (fdmon(jz)-fdnmon(jz))*Phimm(iii)+suu
      enddo
      sumuu = 0.5d0*(fdm-fdnm)*suu+sumuu
      sumuljm = 0.d0
      sumuw = (fdm-fdnm)*csurf+sumuw                  
      
      bclamb = 2.d0*(dlog(ehbclam(iz))-scalem)
      belamb = dlog(ebelam(iz))-emscale
      bclanmb = 2.d0*(dlog(ehbclanm(iz))-scalem)
      belanmb = dlog(ebelanm(iz))-emscale      

      aabclamb = 2.d0*(dlog(ehbclam(iz))-scalem)+donn
      aabelamb = dlog(ebelam(iz))-emscale+donn
      aabclanmb = 2.d0*(dlog(ehbclanm(iz))-scalem)-donn
      aabelanmb = dlog(ebelanm(iz))-emscale-donn      
      
      tfdc = fdm-tfem
      tfdnc = fdnm-tfenm      
      cdt = cdtot(iz)*dm3
      pcdt = pis*cdt
      xsi = (1.d0-pcdt)
      rxsi = 1.d0/xsi
      sqrxsi = rxsi*rxsi
      flog = dlog(xsi)
      aex1 = -(c1+1.d0)*flog-0.5d0*(AA1+BB1*pcdt)*pcdt*sqrxsi
      aex2 = -(c2+1.d0)*flog-0.5d0*(AA2+BB2*pcdt)*pcdt*sqrxsi
      exFreen = (fdm-tfem)*Y*(aex2-aex1)+0.5d0*tfem*aex2+
     *(fdnm-tfenm)*Y*(aex2-aex1)+0.5d0*tfenm*aex2
      sumexFreen = exFreen+sumexFreen
      csumexFreen = tfdc*Y*(aex2-aex1)+0.5d0*tfem*aex2+
     *tfdnc*Y*(aex2-aex1)+0.5d0*tfenm*aex2+
     *csumexFreen
      sumFid = sumFid+tfdc*bclamb+tfem*belamb-fdm*rrnmon+
     *tfdnc*bclanmb+tfenm*belanmb-fdnm*rrnmon
      csumFid = csumFid+tfdc*bclamb+tfem*belamb-fdm*rrnmon+
     *tfdnc*bclanmb+tfenm*belanmb-fdnm*rrnmon
      altsumF = altsumF+tfdc*aabclamb+tfem*aabelamb-fdm*rrnmon+
     *tfdnc*aabclanmb+tfenm*aabelanmb-fdnm*rrnmon
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
c      write(*,*) 'csumF = ',csumF      
c      omega = sumFid+sumuu+sumuw+sumuljm
c      comega = csumFid-sumuu-sumuljm
c      aomega = altsumFid-sumuu-sumuljm
      omega = sumF
      comega = csumF
      aomega = altsumF
      
      write(*,*) 'grand pot. (excl. surf-surf):  ',omega
      write(*,*) omega
      write(*,*) comega
      write(*,*) 'alt. grand pot. (excl. surf-surf):  ',aomega
      write(*,*) aomega      

c     uss = -2.d0*pi*surfdens*surfdens*(h+dm1)*rrT
      uss = -2.d0*pi*surfdens*surfdens*(h+dmw)*rrT                      
      write(*,*) 'u(surf-surf), uss  = ',uss
      totomega = omega+uss
      ctotomega = comega+uss
      atotomega = aomega+uss
      write(*,*)            
      write(*,*) 'total grand pot. (incl. surf-surf):  ',totomega
c      write(*,*) totomega+Pbhs*(h+dm1)
c      write(*,*) ctotomega+Pbhs*(h+dm1)
c      write(*,*) atotomega-2.d0*surfdens*donn+Pbhs*(h+dm1)
      write(*,*) totomega+Pbhs*(h+dmw)
      write(*,*) ctotomega+Pbhs*(h+dmw)
      write(*,*) atotomega-2.d0*surfdens*donn+Pbhs*(h+dmw)            
      write(*,*) donn
      write(*,*) surfdens
      write(*,*)      
      write(*,*) 'elec,donn (fin) = ',elec,donn
      write(*,*) 'surfdens = ',surfdens      
      rewind idpot
      write(idpot,*) donn

      fp1S = fdmon(istp1)
      fn1S = fdmon(istp1+2)
      c0Skv = fdmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fmm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'cationic monomer contact density - quad. extr.: ',fmm
      fp1S = fdnmon(istp1)
      fn1S = fdnmon(istp1+2)
      c0Skv = fdnmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fnm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'anionic monomer contact density - quad. extr.: ',fnm      
      Pidc = fmm+fnm
      write(*,*) 
      write(*,*) 'ideal part of contact pressure: '
      write(*,'(1e25.14)') Pidc
      write(*,*) Pidc-Pbhs
      write(*,*) Pidc-Pbhs-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) Pidc-Pbhs-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) 
      
      avfdm = avfdm*dz/h
      avfdnm = avfdnm*dz/h      
      write(*,*) 'avfdm = ',avfdm
      write(*,*) 'avfdnm = ',avfdnm      
      write(*,*) 'donn = ',donn
      write(*,*) donn
      fdonn = donn
      write(*,*) 'fdonn = donn = ',fdonn
      write(*,*) fdonn
       
c     z = zmitt
c     z = (h+dm1)/2.d0
      z = (h+dmw)/2.d0            
      du = 0.d0
      zp = closew-0.5d0*dz
      write(*,*) 'zmitt,zp(init) = ',z,zp
      do j = istp1,islut      
      zp = zp+dz
      diffz = dabs(zp-z)
      phic = -2.d0*pi*diffz*rrT*dz      
      du = phic*(fdmon(j)-fdnmon(j))+du      
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
      du = phic*(fdmon(j)-fdnmon(j))+du
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
c      write(*,*) atotomega+surfdens*donn            

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
      du = Phimm(kkk)*(fdmon(j)-fdnmon(j))+du
      enddo
      write(92,*) z,du+csurf
      enddo      
      write(*,*) 'ddmax,niter: ',ddmax,niter      
      
      STOP
      END

      subroutine CALCCD
      implicit double precision (a-h,o-z)
      include 'PB.inc'
      fp1S = (fdmon(istp1)+fdnmon(istp1))
      fn1S=(fdmon(istp1+2)+fdnmon(istp1+2))
      c0Skv=
     *(fdmon(istp1+1)+fdnmon(istp1+1))
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
      f2 = (fdmon(istp1)+fdnmon(istp1))
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = 2.d0*(f2-f1)*((z2-z)**4-(z1-z)**4)
      zs = z2-dz
      do jz = istp1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdmon(jz+1)+fdnmon(jz+1))
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdmon(iz+ism)+fdnmon(iz+ism))
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
      f1 = (fdmon(iz-ism)+fdnmon(iz-ism))
      f2 =
     *(fdmon(iz-ism+1)+fdnmon(iz-ism+1))
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do jz = iz-ism+1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdmon(jz+1)+fdnmon(jz+1))
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = (fdmon(iz+ism)+fdnmon(iz+ism))
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
      ae2(iz) = -(c2+1.d0)*flog-0.5d0*(AA2+BB2*pcdt)*pcdt*sqrxsi
      daex1 = rxsi*(c1+1.d0-0.5d0*AA1*rxsi*(1.d0+2.d0*pcdt*rxsi)-
     *BB1*pcdt*rxsi*(1.d0+pcdt*rxsi))
      daex2 = rxsi*(c2+1.d0-0.5d0*AA2*rxsi*(1.d0+2.d0*pcdt*rxsi)-
     *BB2*pcdt*rxsi*(1.d0+pcdt*rxsi))
      convp(iz) = (Y*(fdmon(iz)-fem(iz))*(daex2-daex1)+
     *0.5d0*fem(iz)*daex2+
     *Y*(fdnmon(iz)-fenm(iz))*(daex2-daex1)+
     *0.5d0*fenm(iz)*daex2)*pis*dhs3
      enddo
      jz = imitt+1
      do iz = imitt+1,islut      
      jz = jz-1      
      convp(iz) = convp(jz)
      ae1(iz) = ae1(jz)
      ae2(iz) = ae2(jz)
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
      du = Phimm(kkk)*(fdmon(jz)-fdnmon(jz))+du      
      enddo
      emtrams = trams+0.5d0*ae2(iz)
      cmtrams = trams+Y*(ae2(iz)-ae1(iz))
      strams = trams+ae1(iz)
      ebelam(iz) = dexp(-emtrams+emscale-(du+donn+csurf))
      ehbclam(iz) = dexp(-0.5d0*(cmtrams+du+donn+csurf)+scalem)
      ebelanm(iz) = dexp(-emtrams+emscale+du+donn+csurf)
      ehbclanm(iz) = dexp(-0.5d0*(cmtrams-du-donn-csurf)+scalem)      
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
      du = Phimm(kkk)*(fdmon(jz)-fdnmon(jz))+du      
      enddo
      emtrams = trams+0.5d0*ae2(iz)
      cmtrams = trams+Y*(ae2(iz)-ae1(iz))
      strams = trams+ae1(iz)
      ebelam(iz) = dexp(-emtrams+emscale-(du+donn+csurf))
      ehbclam(iz) = dexp(-0.5d0*(cmtrams+du+donn+csurf)+scalem)
      ebelanm(iz) = dexp(-emtrams+emscale+du+donn+csurf)
      ehbclanm(iz) = dexp(-0.5d0*(cmtrams-du-donn-csurf)+scalem)            
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
c     eeblam(iz) = eeblam(jz)
      ebelam(iz) = ebelam(jz)            
      ehbclam(iz) = ehbclam(jz)
      ebelanm(iz) = ebelanm(jz)
      ehbclanm(iz) = ehbclanm(jz)      
      enddo
      
      return 
      end

      
