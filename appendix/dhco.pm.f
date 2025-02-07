      program platem
      implicit double precision (a-h,o-z)
      dimension fdmon(10001),fem(10001)
      ins = 49
      iefil = 54
      ifc = 56
c      isig = 58
      iamu = 61
      open (ins,file='input.dftp',form='formatted')
      open (ifc,file='fcdfil',form='formatted')
      open (iefil,file='dhfil',form='formatted')
c      open (isig,file='sigmasolv',form='formatted')
c      open (iamu,file='amufil',form='formatted')
      rewind ins
      rewind iefil
c      rewind iamu
      read(iefil,*) dh
      read(ins,*) bdm
      read(ins,*) bds
      read(ins,*) nmon
      read(ins,*) trams
      read(ins,*) ioimaxm
      read(ins,*) ntrams
      read(ins,*) T
      read(ins,*) h
c      read(ins,*) refdz
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
      read(ins,*) ess
      read(ins,*) ems
c      read(ins,*) bdist
c      read(ins,*) esqw
c      read(ins,*) dsqw
c      read(ins,*) bl
c      rewind isig
c      read(isig,*) sigs
c      read(isig,*) sclosew
      sclosew = 0.d0
c      read(iamu,*) trams
c      read(iamu,*) dz
      rdz = 1.d0/dz
      irdz = int(rdz+0.001d0)
      nfack = int(h/dz+0.01d0)
      istart = int(closew/dz+0.01d0)
      istp1 = istart+1 
      islut = int((h-closew)/dz+0.01d0)      
      ists = int(sclosew/dz+0.01d0)
      istp1s = ists+1 
      isluts = int((h-sclosew)/dz+0.01d0)      
      imitt = nfack/2
      rewind ifc
      do 13 iz = istp1s,imitt
 13   read(ifc,*) trams,fdmon(iz),fem(iz)

      h = h+dh
      newnf = int(h/dz+0.01d0)
      newimitt = newnf/2
      newisluts = int((h-sclosew)/dz+0.01d0)

      write(*,*) 'imitt,newimitt = ',imitt,newimitt 
      write(*,*) 'isluts,newisluts = ',isluts,newisluts
      
      do 14 iz = imitt+1,newimitt
      fdmon(iz) = fdmon(imitt)
 14   fem(iz) = fem(imitt)
      jz = newimitt+1
      do 15 iz = newimitt+1,newisluts
      jz = jz-1
      fdmon(iz) = fdmon(jz)
 15   fem(iz) = fem(jz)
      rewind ifc
      z = sclosew-0.5d0*dz
      do 16 iz = istp1s,newisluts
      z = z+dz
 16   write(ifc,'(1f21.5,2e25.14)') z,fdmon(iz),fem(iz)


      rewind ins
c      write(ins,*) bdm
c      write(ins,*) bds
      read(ins,*) bdm
      read(ins,*) bds
      write(ins,*) nmon
      write(ins,*) trams
      write(ins,*) ioimaxm
      write(ins,*) ntrams
      write(ins,*) T
      write(ins,*) h
      write(ins,*) dz
      write(ins,*) closew
      write(ins,*) dmm
      write(ins,*) dms
      write(ins,*) kread
      write(ins,*) hwdist
      write(ins,*) awm
      write(ins,*) rwm
      write(ins,*) aws
      write(ins,*) rws
      write(ins,*) ufact
      write(ins,*) ess
      write(ins,*) ems
c      write(ins,*) bdist 
c      write(ins,*) esqw
c      write(ins,*) dsqw
c      write(ins,*) bl
      stop
      end


