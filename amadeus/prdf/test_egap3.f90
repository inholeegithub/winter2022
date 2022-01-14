!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module eigenvalues
       implicit none
       private
       save
       integer ne,nk,nbandi
       real*8 ef,temper
       real*8, allocatable :: eiv(:,:),wgt(:),cblowest(:),vbhighest(:)

       public :: ne,nk,nbandi,eiv,wgt,ef,temper,cblowest,vbhighest
       end module eigenvalues
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine egp_test(otname,einame,egp1,egp2,lfault)
       implicit none
       logical lfault
       character*280 otname,einame
       real*8 egp1,egp2
       logical lfault1,lfault2,lfault3

       lfault=.false.
       call read_outcar1(otname,lfault1)
       call read_eigenval(einame,lfault2)
       call check_gap(egp1,egp2,lfault3)
       if(lfault1 .or. lfault2 .or. lfault3) lfault=.true.
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_outcar1(otname,lfault1)
       USE eigenvalues, ONLY : ef
       implicit none
       character*280 otname
       logical lfault1
       character*7 ctest7
       logical lfault
  
       lfault=.false.
       ef=1.d111
       open(81,file=trim(otname),form='formatted')
       do 
       read(81,*,err=911,end=999) ctest7
       if(ctest7 == 'E-fermi')then
       backspace(81)
       read(81,101,err=911,end=999) ef
                              endif
  101  format(10x,f9.4)
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(81)
!      write(6,*) ef,' ef from OUTCAR'
!
       if(lfault)then
       ef=1.d111
       write(6,*) 'there is a falut with OUTCAR'
                 endif
       lfault1=lfault
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_eigenval(einame,lfault1)
       USE eigenvalues, ONLY : ef,wgt,eiv,ne,nk,nbandi,temper
       implicit none
       character*280 einame
       logical lfault1
       integer j,i,ik
       real*8 tmpx,tmpy,tmpz,sol1,sol2,tol
       real*8, external ::  fdft,zeroin
       logical lfault

       lfault=.false.
       open(81,file=trim(einame),form='formatted')
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999) ne,nk,nbandi
       allocate(eiv(nbandi,nk))
       allocate(wgt(nk))
       do ik=1,nk
!      write(6,*) ik
       read(81,*,err=911,end=999)
       read(81,*,err=911,end=999) tmpx,tmpy,tmpz,wgt(ik)
       do i=1,nbandi
       read(81,*,err=911,end=999) j,eiv(i,ik)
       enddo
       enddo
       goto 999
  911  continue
       lfault=.true.
  999  continue
       close(81)
! just skip new calculation, just use the ef value from OUTCAR
!      write(6,*) sum(wgt)
!      1.d-5 eV = 3.1577464d0 / (2.d0*13.6058d0)  K
!      au2kelvin= 3.1577464d5
       temper=1.d-4
       sol1=ef-1.d0
       sol2=ef+1.d0
!      write(6,*) sol1,fdft(sol1)
!      write(6,*) sol2,fdft(sol2)
!      tol=1.d-16 ; ef=zeroin(sol1,sol2,fdft,tol)
!      write(6,*) ef,' ef from calculaton, EIGENVAL'
       lfault1=lfault
       write(6,*) 'in eigenval',nk
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       real*8 function fdft(x)
       USE eigenvalues, ONLY : ne,wgt,eiv,temper,nbandi,nk
       implicit none
       real*8 x
       integer ik,i
       real*8 xxr

       fdft=float(ne)
       do ik=1,nk
       do i=1,nbandi
       xxr=(eiv(i,ik)-x)/temper ; if(xxr >  100.d0) xxr=100.d0 ; if(xxr < -100.d0) xxr=-100.d0
       fdft=fdft-(2.d0*wgt(ik))/(1.d0+exp(xxr))
       enddo
       enddo
       return
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine check_gap(egp1,egp2,lfault1)
       USE eigenvalues, ONLY : temper,ef,wgt,eiv,nbandi,nk,cblowest,vbhighest
       implicit none
       logical lfault1
       real*8 egp1,egp2
       integer ik,i,i1,i2
       real*8 t1,t2,test,obj,tst1,tst2
       real*8 occtol
       real*8 smallest,dgap
       logical lfault
       real*8 xxr
       real*8, allocatable :: wrk7(:),wrk8(:)
       integer, allocatable :: iwrk7(:),iwrk8(:)

       lfault=.false.
       occtol=1.d-10
       allocate(cblowest(nk),vbhighest(nk))
!      write(6,*) ef,' ef'
       dgap=1.d119
       do ik=1,nk

       i2=nbandi
       do i=nbandi,1,-1
       if(eiv(i,ik) >= ef)then
       xxr=(eiv(i,ik)-ef)/temper ; if(xxr >  100.d0) xxr=100.d0 ; if(xxr < -100.d0) xxr=-100.d0
       t2=(2.d0*wgt(ik))/(1.d0+exp(xxr))
       if(t2 < occtol) i2=i
                          endif
       enddo
       i1=1
       do i=1,nbandi
       if(eiv(i,ik) <= ef)then
       xxr=(eiv(i,ik)-ef)/temper ; if(xxr >  100.d0) xxr=100.d0 ; if(xxr < -100.d0) xxr=-100.d0
       t1=(2.d0*wgt(ik))/(1.d0+exp(xxr))
       if(t1 > occtol) i1=i
                          endif
       enddo
       cblowest(ik)=eiv(i2,ik)
       vbhighest(ik)=eiv(i1,ik)
!      t1=(2.d0        )/(1.d0+exp((cblowest(ik)-ef)/temper))
!      t2=(2.d0        )/(1.d0+exp((vbhighest(ik)-ef)/temper))
!      write(6,*) t1,t2
!      write(6,*) cblowest(ik)-vbhighest(ik)
       test=cblowest(ik)-vbhighest(ik)
       if(dgap > test)then
       dgap=test
                      endif
       enddo
       allocate(wrk7(nk))
       allocate(iwrk7(nk))
       allocate(wrk8(nk))
       allocate(iwrk8(nk))
       wrk7=cblowest
       wrk8=vbhighest
       call sortnr(nk,wrk7,iwrk7)
       call sortnr(nk,wrk8,iwrk8)
       write(6,*) 'cb lowest'
       write(6,'(20f12.6)') (wrk7(iwrk7(i)),i=1,10)
       write(6,*) 'vb maxima'
       write(6,'(20f12.6)') (wrk8(iwrk8(nk-i+1)),i=1,10)
       write(6,'(20f12.6)') wrk7(iwrk7(1))-wrk8(iwrk8(nk)), wrk7(iwrk7(2))-wrk8(iwrk8(nk)), wrk7(iwrk7(3))-wrk8(iwrk8(nk)), wrk7(iwrk7(4))-wrk8(iwrk8(nk)), wrk7(iwrk7(5))-wrk8(iwrk8(nk)), wrk7(iwrk7(6))-wrk8(iwrk8(nk)), wrk7(iwrk7(7))-wrk8(iwrk8(nk)), wrk7(iwrk7(8))-wrk8(iwrk8(nk))
       smallest=minval(cblowest)-maxval(vbhighest)

       tst1=wrk7(iwrk7(1))-wrk8(iwrk8(nk))
       tst2=wrk7(iwrk7(2))-wrk8(iwrk8(nk))
       if(abs(tst1-dgap) <1.d-8)then
       write(6,*) 'direct gap',tst1,'indirect gap',tst2
       egp1=tst2
       egp2=tst1
       smallest=tst2
                                endif
       if(smallest < dgap )then
       write(6,*) 'indirect gap, smallest indirect gap ', smallest,dgap
       egp1=smallest
       egp2=dgap
                           endif
       if(smallest <  dgap)then
       write(6,'(a17,2x,2f18.8,3x,i5)') 'indirect band gap',dgap,smallest,nk
                           else
       write(6,'(a15,2x,2f18.8,3x,i5,1x,a11)') 'direct band gap',dgap,smallest,nk,'+++++++++++'
                           endif
       call check_dist(lfault)
       write(6,*) '&'
       xxr=2.0d0
       call check_dist_ef(lfault,xxr)
       xxr=2.5d0
       call check_dist_ef(lfault,xxr)
       xxr=3.0d0
       call check_dist_ef(lfault,xxr)
       if(allocated(eiv)) deallocate(eiv)
       if(allocated(wgt)) deallocate(wgt)
       if(allocated(cblowest)) deallocate(cblowest)
       if(allocated(vbhighest)) deallocate(vbhighest)
       deallocate(wrk7,wrk8)
       deallocate(iwrk7,iwrk8)
       lfault1=lfault
       end
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine check_dist(lfault1)
       USE eigenvalues, ONLY : eiv,nbandi,nk
       implicit none
       logical lfault1
       integer ik,i,ih
       real*8 ehist(0:10000)
       real*8 de,avg,tmp

       avg=0.d0
       tmp=0.d0
       do ik=1,nk
       do i=1,nbandi-1
       de=eiv(i+1,ik)-eiv(i,ik)
       avg=avg+de
       tmp=tmp+1.d0
       enddo
       enddo
       avg=avg/tmp
       write(6,*) 'Wigner -> Poisson'
       tmp=0.d0
       ehist=0.d0
       do ik=1,nk
       do i=1,nbandi-1
       de=eiv(i+1,ik)-eiv(i,ik)
       de=de/avg
       ih=de/0.01d0
       ehist(ih)=ehist(ih)+1.d0
       tmp=tmp+1.d0
       enddo
       enddo
       tmp=sum(ehist)*0.01d0
       ehist=ehist/tmp
       do i=1,1000
       write(6,*) float(i-1)*0.01d0, ehist(i)
       enddo

       end
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine check_dist_ef(lfault1,ewindow)
       USE eigenvalues, ONLY : eiv,nbandi,nk,ef
       implicit none
       logical lfault1
       real*8 ewindow
       integer ik,i,ih
       real*8 ehist(0:10000)
       real*8 de,avg,tmp

       write(6,'(a1,2x,f18.8,2x,a10,2x,f18.8)') '#',ef,'ef,ewindow',ewindow
       avg=0.d0
       tmp=0.d0
       do ik=1,nk
       do i=1,nbandi-1
       if(abs(eiv(i,ik)-ef) < ewindow)then
       de=eiv(i+1,ik)-eiv(i,ik)
       avg=avg+de
       tmp=tmp+1.d0
                                      endif
       enddo
       enddo
       avg=avg/tmp
!      write(6,*) 'Wigner -> Poisson'
       tmp=0.d0
       ehist=0.d0
       do ik=1,nk
       do i=1,nbandi-1
       if(abs(eiv(i,ik)-ef) < ewindow)then
       de=eiv(i+1,ik)-eiv(i,ik)
       de=de/avg
       ih=de/0.01d0
       ehist(ih)=ehist(ih)+1.d0
       tmp=tmp+1.d0
                                      endif
       enddo
       enddo
       tmp=sum(ehist)*0.01d0
       ehist=ehist/tmp
       do i=1,1000
       write(6,*) float(i-1)*0.01d0, ehist(i)
       enddo

       write(6,*) '&'
       end
!
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program estimation_egap
       USE eigenvalues, ONLY : ef,wgt,eiv,ne,nk,nbandi,temper
       implicit none
       real*8 egp1,egp2
       character*280 einame,otname
       logical lfault1
!
       otname='OUTCAR'
       einame='EIGENVAL'
       call egp_test(otname,einame,egp1,egp2,lfault1)
       stop
       end program estimation_egap
!
      subroutine sortnr(n,arrin,indx)
!     sorts an array by the heapsort method
!     w. h. preuss et al. numerical recipes
      implicit real*8 (a-h,o-z)
      dimension arrin(n),indx(n)
      do 11 j=1,n
        indx(j)=j
 11   continue
      l=n/2+1
      ir=n
 10   continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
 20     if(j.le.ir)then
          if(j.lt.ir)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
      end
