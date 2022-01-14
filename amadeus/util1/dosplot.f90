!234567890
       implicit none
       real*8 tmp,tmq,efermi
       real*8, allocatable :: egrid(:),edos(:)
       integer i,npt
       character*7 ch7
       character*1 ch1
  
       call system('grep E-fermi OUTCAR > zdelete')
       open(11,file='zdelete',form='formatted')
       read(11,*) ch7,ch1,efermi
       close(11, status='delete')
!      write(6,*) efermi
       open(1,file='DOSCAR',form='formatted')
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*) tmp,tmq,npt
!      write(6,*) npt
       allocate(egrid(npt))
       allocate(edos(npt))
       do i=1,npt
       read(1,*) egrid(i),edos(i)
       enddo
       close(1)
       do i=1,npt
       write(6,*) egrid(i)-efermi,edos(i)
       enddo
       deallocate(egrid)
       deallocate(edos)
       stop
       end
