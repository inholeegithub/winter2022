!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program bs_plot
       implicit none
       logical lfault1,lfault2
       character*80 einame,otname
       real*8 ef0

       otname='OUTCAR'
       einame='EIGENVAL'
       call read_outcar2(otname,ef0,lfault1)
       call read_plot_eigenval(einame,ef0,lfault2)
!      call plot_dos(ef0)

       end program bs_plot
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_plot_eigenval(einame,ef0,lfault1)
       implicit none
       character*80 einame
       real*8 ef0
       logical lfault1
       integer nk,nbandi,ne
       integer nstn
       integer j,i,ik
       real*8 tmpx,tmpy,tmpz,sol1,sol2,tol
       real*8 xxr,y1,y2,teste,ef
       real*8 array(10000,2)
       logical lfault
       real*8, allocatable :: eiv(:,:),wgt(:)

       ef=ef0
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
       lfault1=lfault
       eiv=eiv-ef ; ef=0.d0
!
!      teste=5.52
       teste=ef
!      nstn=40
       open(16,file='KPOINTS',form='formatted')
       read(16,*)
       read(16,*) nstn
       close(16)
!
       open(17,file='fort.17',form='formatted')
       do i=1,nbandi
       xxr=0.d0
       y1=-1.d222
       y2= 1.d222
       do ik=1,nk
       write(6,'(2f18.8)') xxr,eiv(i,ik)
       if(y1 < eiv(i,ik)) y1=eiv(i,ik)
       if(y2 > eiv(i,ik)) y2=eiv(i,ik)
       xxr=xxr+0.01d0
       enddo
       write(6,*) '&'
       write(17,'(2f18.8,2x,f18.8)') y2,y1,y1-y2
       enddo
!
       xxr=0.d0
       write(6,*) xxr,minval(eiv)-1.d0
       write(6,*) xxr,maxval(eiv)+1.d0
       write(6,*) '&'
       do ik=1,nk
       xxr=xxr+0.01d0
       if(mod(ik,nstn)==0)then
       write(6,*) xxr,minval(eiv)-1.d0
       write(6,*) xxr,maxval(eiv)+1.d0
       write(6,*) '&'
                          endif
       enddo
       close(17)
!
       open(17,file='fort.17',form='formatted')
       do i=1,nbandi
       read(17,*) array(i,1),array(i,2)
       enddo
       close(17)
!
       open(18,file='fort.18',form='formatted')
       do i=1,nbandi-1
       if(array(i,2)< teste .and. array(i+1,1) > teste) write(18,*) array(i+1,1)-array(i,2)
       enddo
       close(18)
!
       deallocate(eiv,wgt)
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine read_outcar2(otname,ef,lfault1)
       implicit none
       character*80 otname
       logical lfault1
       real*8 ef
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
!234567890
       subroutine plot_dos(ef0)
       implicit none
       integer npt
       real*8 ef0
       real*8 tmp1,tmp2,tmp3,tmp4
       integer i
       real*8, allocatable :: ee(:),dos(:),sdos(:)

       open(11,file='DOSCAR',form='formatted')
       read(11,*)
       read(11,*)
       read(11,*)
       read(11,*)
       read(11,*)
       read(11,*) tmp1,tmp2,npt,tmp3,tmp4
       allocate(ee(npt),dos(npt),sdos(npt))
       do i=1,npt
       read(11,*) ee(i),dos(i),sdos(i)
       ee(i)=ee(i)-ef0
       enddo
       close(11)

       do i=1,npt
       write(6,*) ee(i),dos(i)
       enddo
!      write(6,*) '&'
!      do i=1,npt
!      write(6,*) ee(i),sdos(i)
!      enddo
       deallocate(ee,dos,sdos)
       end subroutine plot_dos
