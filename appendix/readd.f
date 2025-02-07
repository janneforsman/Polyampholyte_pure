       program calculateptvir
       implicit double precision (a-h,o-z)
       dimension h(100),Fnet(100),Pnet(100),Pid(100),gp(100),
     *Pacr(100),Pcoll(100),Pbr(100),Pn2(100),Pc2(100),
     *avbdm(100),avbds(100),Pnlw(100),Ffree(100),Fgraft(100)
       character string*80 
       iFc = 13
c       open (iFc,file='Fh.cfm2',form='formatted')
       open (iFc,file='Fh',form='formatted')
       rewind iFc
       idPc = 14
c       open (idPc,file='dFh.cfm2',form='formatted')
       open (idPc,file='dFh',form='formatted')
       rewind idPc
       iPc = 15
c       open (iPc,file='Ph.cfm2',form='formatted')
       open (iPc,file='Ph',form='formatted')
       rewind iPc
       ig = 16
       minig = ig
c       open (ig,file='h36',form='formatted')
c       ig = ig+1
c       open (ig,file='h38',form='formatted')
c       ig = ig+1
c       open (ig,file='h7',form='formatted')
c       ig = ig+1
       open (ig,file='h8',form='formatted')
       ig = ig+1
       open (ig,file='h9',form='formatted')
       ig = ig+1
       open (ig,file='h10',form='formatted')
       ig = ig+1
       open (ig,file='h11',form='formatted')
       ig = ig+1
       open (ig,file='h12',form='formatted')
       ig = ig+1       
       open (ig,file='h13',form='formatted')
       ig = ig+1
       open (ig,file='h14',form='formatted')
       ig = ig+1
       open (ig,file='h15',form='formatted')
       ig = ig+1
       open (ig,file='h16',form='formatted')
       ig = ig+1
       open (ig,file='h17',form='formatted')
       ig = ig+1
       open (ig,file='h18',form='formatted')
       ig = ig+1
       open (ig,file='h19',form='formatted')
       ig = ig+1
       open (ig,file='h20',form='formatted')
       ig = ig+1
       open (ig,file='h21',form='formatted')
       ig = ig+1
       open (ig,file='h22',form='formatted')
       ig = ig+1
       open (ig,file='h23',form='formatted')
       ig = ig+1
       open (ig,file='h24',form='formatted')
       ig = ig+1
       open (ig,file='h25',form='formatted')
       ig = ig+1
       open (ig,file='h26',form='formatted')
       maxig = ig
       write(*,*) 'hwdist?'
       read(*,*) hwdist
       ngl = maxig-minig+1
       do 1 ig = minig,maxig
       rewind ig
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(string(9:27),*) h(ig)
       h(ig) = h(ig)+2.d0*hwdist
       write(*,*) 'h(ig) = ',h(ig)
 52    read(ig,'(a)') string
       if (string(1:13).eq.' /beta*(Freen') then
       read(ig,'(a)') string
       read(string(1:80),*) Fnet(ig)
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(ig,'(a)') string
       read(string(1:80),*) Pnet(ig)
c       read(ig,'(a)') string
c       read(string(1:80),*) Pn2(ig)
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(string(1:80),*) Pc2(ig)
c       read(ig,'(a)') string
c       read(string(1:80),*) Pnlw(ig)
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(string(1:80),*) Pc2(ig)
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(ig,'(a)') string
c       read(string(1:80),*) Pnlw(ig)
c       write(*,*) 'Fnet(ig),Pnet(ig) = ',Fnet(ig),Pnet(ig)
c       read(ig,'(a)') string
c       read(string(18:29),*) avbdm(ig)
c       read(string(30:80),*) avbds(ig)
       else
       goto 52
       endif
 1     continue
       do ig = minig,maxig
       write(iPc,'(1f14.7,2e25.14)') h(ig),Pnet(ig)-Pnet(maxig),Pnet(ig)
c       write(*,*) h(ig),Fnet(ig),Ffree(ig),Fgraft(ig)
       write(*,'(1f12.5,3e21.12)')h(ig),Fnet(ig)-Fnet(maxig),Fnet(ig)
       write(iFc,'(1f12.5,3e21.12)')h(ig),Fnet(ig)-Fnet(maxig),Fnet(ig)
       enddo

       if (maxig.gt.minig) then
       Pdnm = -(Fnet(ig+1)-Fnet(ig))/(h(ig+1)-h(ig))
       do ig = minig,maxig-1
       Pdnet = -(Fnet(ig+1)-Fnet(ig))/(h(ig+1)-h(ig))
       s = h(ig)+0.5d0*(h(ig+1)-h(ig))
       write(idPc,'(1f14.7,2e25.14)') s,Pdnet
       enddo
       endif
       STOP
       END




