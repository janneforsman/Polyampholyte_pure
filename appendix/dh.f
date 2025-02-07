      program platem
      implicit double precision (a-h,o-z)
      ins = 49
      iefil = 54
c      iamu = 61
      open (ins,file='input.dftp',form='formatted')
      open (iefil,file='dhfil',form='formatted')
c      open (iamu,file='amufil',form='formatted')
      rewind ins
      rewind iefil
c      rewind iamu
      read(iefil,*) dh
      read(ins,*) bdm
      read(ins,*) bds
c      read(ins,'(1f20.16)') bdm
c      read(ins,'(1f20.16)') bds
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
      read(ins,*) ess
      read(ins,*) ems
c      read(ins,*) bdist
c      read(ins,*) esqw
c      read(ins,*) dsqw
c      read(ins,*) bl
c      read(iamu,*) trams
c      read(iamu,*) dz
      rewind ins
      h = h+dh
c      write(ins,'(1f20.16)') bdm
c      write(ins,'(1f20.16)') bds
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
c      write(ins,*) aws
      write(ins,'(1f18.14)') aws
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
