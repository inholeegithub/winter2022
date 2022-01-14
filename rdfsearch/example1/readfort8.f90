!      Usage: ifort -o readfort8.x  numeral.o sortnr.o  readfort8.f90
!             readfort.x
!             python plot_fort11.py   
!             xmgrace fort.11
!234567890
!      Written by In-Ho Lee, KRISS, September 30, 2020.
       program readfort8
       implicit none
       integer nspecies,nproc,isize
       integer i,j,k,natot,iconf,nconf
       real*8 rr,voltol,refvol
       real*8, allocatable :: a1(:,:),a2(:,:),a3(:,:)
       real*8, allocatable :: pos(:,:,:),rvalue(:)
       real*8, allocatable :: sigmamatrix(:,:)
       integer, allocatable :: nelements(:),ispgset(:),indx(:)
       character*2, allocatable :: symbl(:)
       character*280 fname,gname,cwd
       logical lexist

       nproc=8
       isize=4
       open(1,file='rdfsearch.in',form='formatted')
       read(1,*) nspecies
       allocate(symbl(nspecies)) ; allocate(nelements(nspecies))
       allocate(sigmamatrix(nspecies,nspecies))
       read(1,*) (symbl(i),i=1,nspecies)
       read(1,*) (nelements(i),i=1,nspecies)
       read(1,*) refvol,voltol
       do i=1,nspecies
       read(1,*) (sigmamatrix(i,j),j=1,nspecies)
       enddo
       read(1,*)
       read(1,'(a280)') cwd
       close(1)
       cwd=adjustl(cwd)
       do i=1,280
       if(cwd(i:i) == ' ')then
       j=i
                          exit
                          endif
       enddo
       do i=j,280
       cwd(i:i)=' '
       enddo
       cwd=trim(cwd)
       i=len_trim(cwd) ; if(cwd(i:i) /= '/') cwd=trim(cwd)//'/' ; cwd=trim(cwd)
!
       nproc=0
       do i=0,1000-1
       call xnumeral(i,gname,isize) ; gname=trim(gname)
       gname=trim(cwd)//trim(gname)//'/fort.8'
       inquire(file=trim(gname),exist=lexist)
       if(lexist) nproc=nproc+1
       enddo
       write(6,*) nproc,' nproc'
!
       nconf=0
       do i=0,nproc-1
       call xnumeral(i,gname,isize) ; gname=trim(gname)
       gname=trim(cwd)//trim(gname)//'/fort.8'
       open(8,file=trim(gname),form='formatted')
       do k=1,100000
       read(8,*,end=999) 
       read(8,*,end=999) 
       read(8,*,end=999) 
       read(8,*,end=999) 
       read(8,*,end=999) 
       read(8,*,end=999) 
       read(8,*,end=999) 
       read(8,*,end=999) 
       do j=1,sum(nelements)
       read(8,*,end=999) 
       enddo
       nconf=nconf+1
       enddo
  999  continue
       close(8)
       enddo
       natot=sum(nelements)
       allocate(rvalue(nconf)) ; allocate(ispgset(nconf))
       allocate(a1(3,nconf),a2(3,nconf),a3(3,nconf))
       allocate(pos(natot,3,nconf))
       iconf=1
       do i=0,nproc-1
       call xnumeral(i,gname,isize) ; gname=trim(gname)
       gname=trim(cwd)//trim(gname)//'/fort.8'
       open(8,file=trim(gname),form='formatted')
       do k=1,100000
       if(iconf <= nconf)then
       read(8,*,end=991) rvalue(iconf),ispgset(iconf)
       read(8,*,end=991) 
       read(8,*,end=991) a1(1,iconf),a1(2,iconf),a1(3,iconf)
       read(8,*,end=991) a2(1,iconf),a2(2,iconf),a2(3,iconf)
       read(8,*,end=991) a3(1,iconf),a3(2,iconf),a3(3,iconf)
       read(8,*,end=991) 
       read(8,*,end=991) 
       read(8,*,end=991) 
       do j=1,sum(nelements)
       read(8,*,end=991) pos(j,1,iconf),pos(j,2,iconf),pos(j,3,iconf)
       enddo
       iconf=iconf+1
                         else
                         exit
                         endif
       enddo
  991  continue
       close(8)
       enddo
       open(11,file='fort.11',form='formatted')
       allocate(indx(nconf))
       call sortnr(nconf,rvalue,indx)
       isize=6
       iconf=1
       do i=nconf,1,-1
       j=indx(i)
       call xnumeral(iconf,gname,isize) ; gname=trim(gname)
       gname='POSCAR_'//trim(gname)
       if(iconf <= 200)then
       open(1,file=trim(gname),form='formatted')
       write(1,'(f8.4,1x,i5)') rvalue(j),ispgset(j)
       write(1,*) '1.'
       write(1,'(3f22.12)') a1(1,j),a1(2,j),a1(3,j)
       write(1,'(3f22.12)') a2(1,j),a2(2,j),a2(3,j)
       write(1,'(3f22.12)') a3(1,j),a3(2,j),a3(3,j)
       write(1,'(10(1x,a2,1x))') (symbl(k),k=1,nspecies)
       write(1,'(10i5)') (nelements(k),k=1,nspecies)
       write(1,*) 'direct'
       do k=1,natot
       write(1,'(3f20.12)') pos(k,1,j),pos(k,2,j),pos(k,3,j)
       enddo
       close(1)
                       endif
       write(11,'(i5,f12.6)') ispgset(j),rvalue(j) 
       iconf=iconf+1
       enddo
       close(11)
       deallocate(indx)
       deallocate(rvalue,pos,ispgset)
       deallocate(symbl,nelements,sigmamatrix)
       end program readfort8
!234567890
