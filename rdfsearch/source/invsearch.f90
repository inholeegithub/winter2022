!234567890
!      Written by In-Ho Lee, KRISS, April 25 (2004)
       subroutine arb_grid(cdivid,npoint,x_start,x_finish,array_g)
       implicit none
       character(3) cdivid
       integer npoint
       real*8 x_start,x_finish,array_g(npoint)
       real*8 delta
       integer i

       if(npoint == 1)then
       array_g(1)=x_start
       return
       endif
       if(cdivid == 'lin' .or. cdivid =='LIN')then
       delta=(x_finish-x_start)/float(npoint-1)
       do i=1,npoint
       array_g(i)=x_start+delta*float(i-1)
       enddo
       elseif(cdivid == 'inv' .or. cdivid == 'INV')then
       delta=(1.d0/x_finish-1.d0/x_start)/float(npoint-1)
       do i=1,npoint
       array_g(i)=1.d0/(1.d0/x_start+delta*float(i-1))
       enddo
       elseif(cdivid == 'log' .or. cdivid == 'LOG')then
       delta=dlog(x_finish/x_start)/float(npoint-1)
       do i=1,npoint
       array_g(i)=exp(dlog(x_start)+delta*float(i-1))
       enddo
       else
       write(6,*) 'input error, check cdivid'
       stop
       endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       real*8 function m33det(a)
       implicit none
       real*8, dimension(3,3), intent(in)  :: a

       m33det =  a(1,1)*a(2,2)*a(3,3)  &
               - a(1,1)*a(2,3)*a(3,2)  &
               - a(1,2)*a(2,1)*a(3,3)  &
               + a(1,2)*a(2,3)*a(3,1)  &
               + a(1,3)*a(2,1)*a(3,2)  &
               - a(1,3)*a(2,2)*a(3,1)
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       subroutine tolaty0(xyz,a1,a2,a3,ktot)
       implicit none
       integer ktot
       real*8 xyz(ktot,3),a1(3),a2(3),a3(3)
       integer j
       real*8 b(3,3),devid,d1,d2,d3

       devid=a1(1)*a2(2)*a3(3)-a1(2)*a2(1)*a3(3)-a1(1)*a2(3)*a3(2)   &
            +a1(3)*a2(1)*a3(2)+a1(2)*a2(3)*a3(1)-a1(3)*a2(2)*a3(1)
       b(1,1)=-a2(3)*a3(2)+a2(2)*a3(3)
       b(2,1)= a1(3)*a3(2)-a1(2)*a3(3)
       b(3,1)=-a1(3)*a2(2)+a1(2)*a2(3)
       b(1,2)= a2(3)*a3(1)-a2(1)*a3(3)
       b(2,2)=-a1(3)*a3(1)+a1(1)*a3(3)
       b(3,2)= a1(3)*a2(1)-a1(1)*a2(3)
       b(1,3)=-a2(2)*a3(1)+a2(1)*a3(2)
       b(2,3)= a1(2)*a3(1)-a1(1)*a3(2)
       b(3,3)=-a1(2)*a2(1)+a1(1)*a2(2)
       b(:,:)=b(:,:)/devid
       do j=1,ktot
       d1=b(1,1)*xyz(j,1)+b(1,2)*xyz(j,2)+b(1,3)*xyz(j,3)
       d2=b(2,1)*xyz(j,1)+b(2,2)*xyz(j,2)+b(2,3)*xyz(j,3)
       d3=b(3,1)*xyz(j,1)+b(3,2)*xyz(j,2)+b(3,3)*xyz(j,3)
       xyz(j,1)=d1 ; xyz(j,2)=d2 ; xyz(j,3)=d3
       enddo
       do j=1,ktot
       xyz(j,1)=xyz(j,1)-anint(xyz(j,1))
       xyz(j,2)=xyz(j,2)-anint(xyz(j,2))
       xyz(j,3)=xyz(j,3)-anint(xyz(j,3))
       enddo
       do j=1,ktot
       if(xyz(j,1) < 0.d0) xyz(j,1)=xyz(j,1)+1.d0
       if(xyz(j,2) < 0.d0) xyz(j,2)=xyz(j,2)+1.d0
       if(xyz(j,3) < 0.d0) xyz(j,3)=xyz(j,3)+1.d0
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       subroutine getinv(natoms,nsigma,ninv,a1,a2,a3,dir,sktinv,sref,nrref,rref)
       implicit none
       integer natoms,nsigma,ninv,nrref,n1,n2,n3
       real*8 sktinv(ninv,nsigma,nrref,natoms),sref(nsigma),rref(nrref)
       real*8 a1(3),a2(3),a3(3),dir(natoms,3)
       real*8 skinv,vkinv,tk(3,3),amat(3,3),sigma
!      real*8 tkinv1,tkinv2,tkinv3
       real*8 ainv1,ainv2,ainv3,ainv4,ainv5,ainv6,ainv7
       real*8 xi,yi,zi,xj,yj,zj,trace1,trace2
       real*8 ck,pi,arg,brg,rij,x,y,z,tmp,deno
       integer iatom,jatom,ii,jj,kk,ksigma,iref
       real*8, external :: m33det

       tmp=2.d0*rref(nrref)
       n1=int(tmp/sqrt(dot_product(a1,a1)))+1
       n2=int(tmp/sqrt(dot_product(a2,a2)))+1
       n3=int(tmp/sqrt(dot_product(a3,a3)))+1
       sktinv(:,:,:,:)=0.d0
       pi=4.d0*atan(1.d0)
       do iref=1,nrref
       do ksigma=1,nsigma
       sigma=sref(ksigma)
       ck=1.d0/(2.d0*pi)**(0.5d0)/sigma
       do iatom=1,natoms
       xi=dir(iatom,1)*a1(1)+dir(iatom,2)*a2(1)+dir(iatom,3)*a3(1)
       yi=dir(iatom,1)*a1(2)+dir(iatom,2)*a2(2)+dir(iatom,3)*a3(2)
       zi=dir(iatom,1)*a1(3)+dir(iatom,2)*a2(3)+dir(iatom,3)*a3(3)
       deno=0.d0
       skinv=0.d0 ; vkinv=0.d0 ; tk(:,:)=0.d0
       do jatom=1,natoms
       do ii=-n1,n1
       do jj=-n2,n2
       do kk=-n3,n3
       xj=dir(jatom,1)*a1(1)+dir(jatom,2)*a2(1)+dir(jatom,3)*a3(1)
       yj=dir(jatom,1)*a1(2)+dir(jatom,2)*a2(2)+dir(jatom,3)*a3(2)
       zj=dir(jatom,1)*a1(3)+dir(jatom,2)*a2(3)+dir(jatom,3)*a3(3)
       xj=xj+ii*a1(1)+jj*a2(1)+kk*a3(1)
       yj=yj+ii*a1(2)+jj*a2(2)+kk*a3(2)
       zj=zj+ii*a1(3)+jj*a2(3)+kk*a3(3)
       x=xi-xj ; y=yi-yj ; z=zi-zj 
       rij=sqrt(x*x+y*y+z*z)
       if(rij < 1.d-3) cycle
       arg=-((rij-rref(iref))/sigma)**2/2.d0 ; if(arg < -30.d0) arg=-30.d0
       brg=ck*exp(arg) 
       x=x/rij
       y=y/rij
       z=z/rij
       skinv=skinv+brg
       vkinv=vkinv+sqrt(x*x+y*y+z*z)*brg
       tk(1,1)=tk(1,1)+brg*x*x ; tk(1,2)=tk(1,2)+brg*x*y
       tk(1,3)=tk(1,3)+brg*x*z ; tk(2,1)=tk(2,1)+brg*y*x
       tk(2,2)=tk(2,2)+brg*y*y ; tk(2,3)=tk(2,3)+brg*y*z
       tk(3,1)=tk(3,1)+brg*z*x ; tk(3,2)=tk(3,2)+brg*z*y
       tk(3,3)=tk(3,3)+brg*z*z
       deno=deno+1.d0
       enddo
       enddo
       enddo
       enddo
       skinv=skinv/deno
       vkinv=vkinv/deno
       tk=tk/deno
       amat=matmul(tk,tk)
       trace1=tk(1,1)+tk(2,2)+tk(3,3)
       trace2=amat(1,1)+amat(2,2)+amat(3,3)
       ainv1=trace1 ; ainv2=(trace1**2-trace2)/2.d0 ; ainv2=max(ainv2,0.d0)
       ainv3=m33det(tk) ; ainv4=ainv1**2-2.d0*ainv2 ; ainv4=max(ainv4,0.d0)
       ainv5=ainv1**3-3.d0*ainv1*ainv2+3.d0*ainv3
       ainv6=2.d0/9.d0*(ainv1**2-3.d0*ainv2) ; ainv6=max(ainv6,0.d0)
       ainv7=2.d0/54.d0*(-9.d0*ainv1*ainv2+27.d0*ainv3+2.d0*ainv1**3)
       if(ninv == 9)then
       sktinv(1,ksigma,iref,iatom)=sktinv(1,ksigma,iref,iatom)+skinv
       sktinv(2,ksigma,iref,iatom)=sktinv(2,ksigma,iref,iatom)+vkinv
       sktinv(3,ksigma,iref,iatom)=sktinv(3,ksigma,iref,iatom)+ainv1
       sktinv(4,ksigma,iref,iatom)=sktinv(4,ksigma,iref,iatom)+sqrt(ainv2)
       sktinv(5,ksigma,iref,iatom)=sktinv(5,ksigma,iref,iatom)+sqrt(ainv4)
       sktinv(6,ksigma,iref,iatom)=sktinv(6,ksigma,iref,iatom)+sqrt(ainv6)
       sktinv(7,ksigma,iref,iatom)=sktinv(7,ksigma,iref,iatom)+(ainv3)**(1.d0/3.d0)
       sktinv(8,ksigma,iref,iatom)=sktinv(8,ksigma,iref,iatom)+(ainv5)**(1.d0/3.d0)
       sktinv(9,ksigma,iref,iatom)=sktinv(9,ksigma,iref,iatom)+(ainv7)**(1.d0/3.d0)
                    endif
       enddo
       enddo
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       subroutine invfeature1(ninv0,nsigma0,nrref0,ndeg,nspecies,nelements,symbl,qqq,vork)
       implicit none
       integer ninv0,nsigma0,nrref0,ndeg
       integer nspecies,nelements(nspecies)
       character*2 symbl(nspecies)
       real*8 qqq(ndeg),vork(ninv0*nsigma0*nrref0)
       integer ninv,nsigma,nrref
       integer natoms
       real*8 scale0,a1(3),a2(3),a3(3),vtest,r6(6),amat(3,3)
       real*8, allocatable :: dir(:,:)
       integer i,j,k,i0,j0,iref,ish
       character*1 ch1
       CHARACTER(3) cdivid
       real*8 x_start,x_finish
       real*8, allocatable :: sktinv(:,:,:,:)
       real*8, allocatable :: sref(:)
       real*8, allocatable :: rref(:),work(:,:,:)

       ish=ndeg-6
       do i=1,6
       r6(i)=qqq(ish+i)
       enddo
       call latmat(r6,amat,1)
       a1(:)=amat(1,:) ; a2(:)=amat(2,:) ; a3(:)=amat(3,:)
       natoms=sum(nelements)
       allocate(dir(natoms,3)) 
       do i=1,natoms
       dir(i,1)=qqq(3*(i-1)+1)
       dir(i,2)=qqq(3*(i-1)+2)
       dir(i,3)=qqq(3*(i-1)+3)
       enddo
!      cdivid='inv'
!      cdivid='log'
       cdivid='lin'
       nrref=nrref0
       allocate(rref(nrref))
       x_start=0.00d0 ; x_finish=8.00d0
       x_start=0.50d0 ; x_finish=8.00d0
       call arb_grid(cdivid,nrref,x_start,x_finish,rref)
!      cdivid='log'
!      cdivid='inv'
       cdivid='lin'
       nsigma=nsigma0
       allocate(sref(nsigma))
!      x_start=0.200d0  ; x_finish=1.20d0
!      x_start=0.050d0  ; x_finish=0.50d0
       x_start=0.100d0  ; x_finish=0.50d0
       call arb_grid(cdivid,nsigma,x_start,x_finish,sref)
       ninv=9
       ninv=ninv0
       allocate(sktinv(ninv,nsigma,nrref,natoms))
       call getinv(natoms,nsigma,ninv,a1,a2,a3,dir,sktinv,sref,nrref,rref)
       allocate(work(ninv,nsigma,nrref))
       work(:,:,:)=0.d0
       do iref=1,nrref
       do j=1,nsigma
       do k=1,ninv
       i0=0
       do i=1,nspecies
       do j0=1,nelements(i)
       i0=i0+1
       work(k,j,iref)=work(k,j,iref)+sktinv(k,j,iref,i0)
       enddo
       enddo
       enddo
       enddo
       enddo
       work=work/dble(i0)
       j0=0
       do k=1,ninv
       do j=1,nsigma
       do iref=1,nrref
       j0=j0+1
       vork(j0)=work(k,j,iref)
       enddo
       enddo
       enddo
       deallocate(dir)
       deallocate(work)
       deallocate(sktinv)
       deallocate(sref)
       deallocate(rref)
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 21, 2020.
       subroutine getatarget(fname,ninv0,nsigma0,nrref0,vork)
       implicit none
       character*280 fname
       integer ninv0,nsigma0,nrref0
       real*8 vork(ninv0*nsigma0*nrref0)
       integer nspecies
       real*8 scale0,a1(3),a2(3),a3(3),vtest,tmp
       character*2, allocatable:: symbl(:)
       integer, allocatable :: nelements(:)
       real*8, allocatable :: dir(:,:)
       integer i,j,k,i0,j0,ish,natoms,iref
       integer, parameter   :: nlen=1000
       integer              :: nwords,ipos
       character (len=nlen) :: text
       integer ninv,nsigma,nrref
       character*1 ch1
       CHARACTER(3) cdivid
       real*8 x_start,x_finish
       real*8, allocatable :: sktinv(:,:,:,:)
       real*8, allocatable :: sref(:),rref(:),work(:,:,:)
       real*8, external :: atomicmass

       open(1,file=trim(fname),form='formatted')
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,'(a1000)') text
       close(1)
       ipos = 1 ; nwords = 0
       loop: do
       i = verify(text(ipos:), ' ') !-- Find next non-blank.
       if (i == 0) exit loop        !-- No word found.
       nwords = nwords + 1          !-- Found something.
       ipos = ipos + i - 1          !-- Move to start of the word.
       i = scan(text(ipos:), ' ')   !-- Find next blank.
       if (i == 0) exit loop        !-- No blank found.
       ipos = ipos + i - 1          !-- Move to the blank.
       end do loop
       nspecies=nwords
       allocate(symbl(nspecies)) ; allocate(nelements(nspecies))
       open(1,file=trim(fname),form='formatted')
       read(1,*)
       read(1,*) scale0
       read(1,*) a1(1),a1(2),a1(3)
       read(1,*) a2(1),a2(2),a2(3)
       read(1,*) a3(1),a3(2),a3(3)
       if(scale0 > 0.d0)then
       a1=a1*scale0 ; a2=a2*scale0 ; a3=a3*scale0
                        endif
       if(scale0 < 0.d0)then
       vtest=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
            +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
            +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       vtest=abs(vtest)
       vtest=abs(scale0)/vtest ; vtest=vtest**(1.d0/3.d0)
       a1=a1*vtest ; a2=a2*vtest ; a3=a3*vtest
                        endif
       read(1,*) (symbl(i),i=1,nspecies)
       do i=1,nspecies
       symbl(i)=adjustl(symbl(i))
       enddo
       read(1,*) (nelements(i),i=1,nspecies)
       natoms=sum(nelements)
       vtest=(a1(2)*a2(3)-a1(3)*a2(2))*a3(1) &
            +(a1(3)*a2(1)-a1(1)*a2(3))*a3(2) &
            +(a1(1)*a2(2)-a1(2)*a2(1))*a3(3)
       vtest=abs(vtest)
       tmp=0.d0
       do i=1,nspecies
       tmp=tmp+atomicmass(symbl(i))*nelements(i)
       enddo
       tmp=tmp/vtest*(1.660539d-24)/(1.d-24)
       write(6,*) trim(fname), vtest, tmp, ' g/cm^3'
       allocate(dir(natoms,3))
  100  continue
       read(1,*) ch1
       if(ch1 == 'D') ch1='d'
       if(ch1 == 'C') ch1='c'
       if(ch1 == 'K') ch1='c'
       if(ch1 == 'k') ch1='c'
       if(ch1 == 'S') ch1='s'
       if(ch1 == 's') goto 100
       do i=1,natoms
       read(1,*) dir(i,1),dir(i,2),dir(i,3)
       enddo
       close(1)
       if(ch1 == 'c')then
       dir=dir*scale0
       call tolaty0(dir,a1,a2,a3,natoms)
                     endif
!      cdivid='inv'
!      cdivid='log'
       cdivid='lin'
       nrref=nrref0
       allocate(rref(nrref))
!      x_start=0.00d0 ; x_finish=8.00d0
       x_start=0.50d0 ; x_finish=8.00d0
       call arb_grid(cdivid,nrref,x_start,x_finish,rref)
!      cdivid='log'
!      cdivid='inv'
       cdivid='lin'
       nsigma=nsigma0
       allocate(sref(nsigma))
!      x_start=0.200d0  ; x_finish=1.20d0
!      x_start=0.050d0  ; x_finish=0.50d0
       x_start=0.100d0  ; x_finish=0.50d0
       call arb_grid(cdivid,nsigma,x_start,x_finish,sref)
       ninv=ninv0
       allocate(sktinv(ninv,nsigma,nrref,natoms))
       call getinv(natoms,nsigma,ninv,a1,a2,a3,dir,sktinv,sref,nrref,rref)
       allocate(work(ninv,nsigma,nrref))
       work(:,:,:)=0.d0
       do iref=1,nrref
       do j=1,nsigma
       do k=1,ninv
       i0=0
       do i=1,nspecies
       do j0=1,nelements(i)
       i0=i0+1
       work(k,j,iref)=work(k,j,iref)+sktinv(k,j,iref,i0)
       enddo
       enddo
       enddo
       enddo
       enddo
       work=work/dble(i0)
       j0=0
       do k=1,ninv
       do j=1,nsigma
       do iref=1,nrref
       j0=j0+1
       vork(j0)=work(k,j,iref)
       enddo
       enddo
       enddo
       deallocate(symbl) ; deallocate(nelements) ; deallocate(dir)
       deallocate(sktinv,sref,rref,work)
       end
!234567890
!      Written by In-Ho Lee, KRISS, October 2 (2020)
       program invsearch
!      use mpi_f08
       implicit none
       include 'mpif.h'
       character*280 fname,gname,cwd
       integer ninv0,nsigma0,nrref0
       integer nspecies,ndeg,ntargets
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       real*8 refvol,voltol,r6(6),amat(3,3),a1(3),a2(3),a3(3)
       integer ii,i0,i,j,k0,jj,ndata,ispg,ij,kl,iter,ish
       real*8 rr,prob,zz,tmp,pi,rsp,tau,dd,zd,probd,probrs
       real*8, allocatable :: sigmamatrix(:,:)
       real*8, allocatable :: vorkp(:),vorkq(:,:),ppp(:)
       real*8, allocatable :: wksp1(:),wksp2(:)
       logical lexist
       logical lpbc,lflag,lvcs
       integer myid,nproc,ierr
       integer isize

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
       if(myid == 0)then
       if(nproc > 1) print *, 'Parallel search:', nproc," processes are alive"
       if(nproc ==1) print *, 'Parallel search:', nproc," process is alive"
                    endif
       call random_seed()
       call random_number(tmp)
       ij=tmp+1000+1+myid
       call random_number(tmp)
       kl=tmp+123+1
       call rmarin(ij,kl)
       if(myid == 0)then
       open(1,file='rdfsearch.in',form='formatted')
       read(1,*) nspecies
       allocate(symbl(nspecies))
       allocate(nelements(nspecies))
       allocate(sigmamatrix(nspecies,nspecies))
       read(1,*) (symbl(i),i=1,nspecies)
       do i=1,nspecies
       symbl(i)=adjustl(symbl(i))
       enddo
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
                    endif
       call MPI_BCAST(nspecies,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       if(myid /= 0)then
       allocate(symbl(nspecies))
       allocate(nelements(nspecies))
       allocate(sigmamatrix(nspecies,nspecies))
                    endif
       call MPI_BCAST(cwd,280,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       i=2*nspecies
       call MPI_BCAST(symbl,i,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       i=nspecies*nspecies
       call MPI_BCAST(sigmamatrix,i,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(nelements,nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(refvol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(voltol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       isize=4
       if(myid == 0)then
       do i=0,nproc-1
       call xnumeral(i,gname,isize) ; gname=trim(gname)
       gname='mkdir '//trim(gname)  ; gname=trim(gname) ; call system(gname)
       enddo
                    endif
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call system('sleep 0.1')
       call xnumeral(myid,gname,isize) ; gname=trim(gname)
       gname=trim(cwd)//trim(gname)//'/fort.8'
       write(6,*) myid,trim(gname)
       open(8,file=trim(gname),form='formatted')
       ndeg=sum(nelements)*3+6
       allocate(ppp(ndeg))
       ninv0=9 ; nsigma0=1 ; nrref0=200
       ndata=ninv0*nsigma0*nrref0
       allocate(wksp1(ndata),wksp2(ndata))
       ntargets=1
       if(myid == 0)then
       open(2,file='targets',form='formatted')
       read(2,*) ntargets
       ntargets=iabs(ntargets)
                    endif
       call MPI_BCAST(ntargets,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       allocate(vorkq(ndata,ntargets+15*15),vorkp(ndata))
       if(myid == 0)then
       if(.true.)then
       do i=1,ntargets
       read(2,*) fname
       gname='fingers/'//trim(fname)//'_finplt'
       write(6,*) trim(gname)
       open(17,file=trim(gname),form='formatted')
       jj=0
       do k0=1,9
       do j=1,200
       jj=jj+1
       read(17,*) tmp,vorkq(jj,i) 
!      write(6,*) tmp,vorkq(jj,i) 
       enddo
       read(17,*)
       enddo
       close(17)
       enddo
       write(6,*) 'fingerprint files are provided',ntargets
                 endif
       if(.true.)then
       isize=5
       ii=ntargets
       do i=1,15*15
       jj=i-1
       ii=ii+1
       call xnumeral(jj,gname,isize) ; gname=trim(gname)
       gname='replica_fingers/'//trim(gname)//'_finplt'
       write(6,*) trim(gname)
       open(17,file=trim(gname),form='formatted')
       jj=0
       do k0=1,9
       do j=1,200
       jj=jj+1
       read(17,*) tmp,vorkq(jj,ii) 
!      write(6,*) tmp,vorkq(jj,ii) 
       enddo
       read(17,*)
       enddo
       close(17)
       enddo
       write(6,*) 'replica fingerprint files are provided',15*15
                 endif
       close(2)
                    endif
       i=ndata*(ntargets+15*15)
       call MPI_BCAST(vorkq,i,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       lpbc=.true. ; lvcs=.true. ; lflag=.true.
       do iter=1,1000000
       ispg=0
!      inquire(file='specific_ispg',exist=lexist)
!      if(lexist)then
!      open(2,file='specific_ispg',form='formatted')
!      read(2,*) ispg
!      close(2)
!      if(ispg < 0 .or. ispg > 230) ispg=0
!                endif
       call gen_latt_site(ispg,ndeg,nspecies,nelements,symbl,sigmamatrix,voltol,refvol,ppp,lpbc,lvcs,lflag)
       call invfeature1(ninv0,nsigma0,nrref0,ndeg,nspecies,nelements,symbl,ppp,vorkp)
       rr=-1.d0 ; i0=1
       do i=1,ntargets+15*15
       call pearson(vorkq(1,i),vorkp,ndata,tmp)
!      call kendl1(vorkq(1,i),vorkp,ndata,tau,zz,prob)
!      tmp=tau
!      write(6,*) tau
!      call spear(vorkq(1,i),vorkp,ndata,wksp1,wksp2,dd,zd,probd,rsp,probrs)
!      tmp=rsp
!      write(6,*) rsp
       if(tmp > rr)then
       rr=tmp ; i0=i
                   endif
       enddo
       if(rr > 0.6d0)then
       ish=ndeg-6
       do i=1,6
       r6(i)=ppp(ish+i) 
       enddo
       call latmat(r6,amat,1)
       a1(:)=amat(1,:) ; a2(:)=amat(2,:) ; a3(:)=amat(3,:)
       write(8,'(f9.4,2x,i5)') rr,ispg
       write(8,*) '1.'
       write(8,'(3f22.12)') a1(1),a1(2),a1(3)
       write(8,'(3f22.12)') a2(1),a2(2),a2(3)
       write(8,'(3f22.12)') a3(1),a3(2),a3(3)
       write(8,'(10(1x,a2,1x))') (symbl(i),i=1,nspecies)
       write(8,'(10(1x,i5,1x))') (nelements(i),i=1,nspecies)
       write(8,*) 'direct'
       do i=1,sum(nelements)
       write(8,'(3f18.12)') ppp(3*(i-1)+1),ppp(3*(i-1)+2),ppp(3*(i-1)+3)
       enddo
       call flush(8)
       tmp=sum(abs(vorkp(:)-vorkq(:,i0)))
!      write(6,*) tmp,' sum abs diff'
       tmp=sqrt(sum((vorkp(:)-vorkq(:,i0))**2))
!      write(6,*) tmp, ' sqrt sum squared'
       tmp=sqrt(dot_product(vorkp(:),vorkp(:)))*sqrt(dot_product(vorkq(:,i0),vorkq(:,i0)))
       tmp=dot_product(vorkp(:),vorkq(:,i0))/tmp
       pi=4.d0*atan(1.d0)
!      write(6,*) tmp,acos(tmp),acos(tmp)*180.d0/pi
       write(6,*) rr,acos(tmp)*180.d0/pi
                     endif
       enddo
       close(8)
!      open(2,iostat=i,file='specific_ispg',status='old')
!      if(i == 0) close(2,status='delete')
       deallocate(ppp)
       deallocate(vorkp,vorkq)
       deallocate(sigmamatrix,nelements,symbl)
       deallocate(wksp1,wksp2)
       call MPI_FINALIZE(ierr)
       end program invsearch
!234567890
