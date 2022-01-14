!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module qlabmod
       implicit none
       private
       save
       integer natot,nspecies,ndeg,nstrc
       real*8 rmax
       integer, allocatable :: itype(:),itypeext(:),ndum(:,:)
       real*8, allocatable :: qext(:),qlab(:,:,:),qlabtest(:,:,:),qlabsave(:,:,:,:)
       real*8, allocatable :: sigmamatrix(:,:)
       complex*16, allocatable :: ctdangl(:,:,:,:)
       logical lpbc
       
       public :: qlab_init,qlab_final,get_qlab,qlab_cmp
       public :: nstrc,qlabsave,qlabtest,qlab
       contains
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine qlab_init(nspecies0,nelements0,sigmamatrix0,rmax0,lpbc0)
       implicit none
       integer nspecies0,nelements0(nspecies0)
       real*8 sigmamatrix0(nspecies0,nspecies0),rmax0
       logical lpbc0
       integer i,j,k

       rmax=rmax0
       lpbc=lpbc0
       nspecies=nspecies0
       allocate(ctdangl(nspecies,nspecies,0:10,-10:10))
       allocate(ndum(nspecies,nspecies))
       natot=sum(nelements0) ; ndeg=6+3*natot
       allocate(itype(natot))
       k=0
       do i=1,nspecies0
       do j=1,nelements0(i)
       k=k+1
       itype(k)=i
       enddo
       enddo
       if(k /= natot)then
       write(6,*) k,natot,' k,natot'
                     stop
                     endif
       allocate(sigmamatrix(nspecies,nspecies))
       sigmamatrix=sigmamatrix0
       allocate(qlab(nspecies,nspecies,0:10))
       allocate(qlabtest(nspecies,nspecies,0:10))
       allocate(qlabsave(nspecies,nspecies,0:10,nstrc))
       end subroutine qlab_init
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine qlab_final()
       implicit none

       deallocate(sigmamatrix)
       deallocate(itype) ; deallocate(ctdangl) ; deallocate(ndum)
       deallocate(qlab) ; deallocate(qlabsave) ; deallocate(qlabtest)
       end subroutine qlab_final
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine qlab_cmp(dista)
       implicit none
       real*8 dista
       real*8 tmp,tmq
       integer i,j,l

       tmq=0.d0
       dista=0.d0 
       do i=1,nspecies
       do j=1,nspecies
       tmp=0.d0
       do l=0,10,2
       tmp=tmp+(qlabtest(i,j,l))**2
       enddo
       tmq=tmq+1.d0
       dista=dista+tmp
       enddo
       enddo
       dista=dista/tmq
       dista=sqrt(dista)
       end subroutine qlab_cmp
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine get_qlab(qqq)
       implicit none
       real*8 qqq(ndeg)
       real*8 pi,x,y,z,r,d1,d2,d3,theta,phi,arg
       real*8 aa1(3),aa2(3),aa3(3),t6(6),cmatrix(3,3),a1(3),a2(3),a3(3)
       integer n1,n2,n3,i,j,k,l,m,nb,ish,iti,itj,ncext(3)
       complex*16 ylm,ctmp(0:10,-10:10),cl(0:10)

       ctdangl(:,:,:,:)=cmplx(0.d0,0.d0) ; ndum(:,:)=0
       if(      lpbc)then
       do j=1,natot
       qqq(3*(j-1)+1)=qqq(3*(j-1)+1)-anint(qqq(3*(j-1)+1))
       qqq(3*(j-1)+2)=qqq(3*(j-1)+2)-anint(qqq(3*(j-1)+2))
       qqq(3*(j-1)+3)=qqq(3*(j-1)+3)-anint(qqq(3*(j-1)+3))
       if(qqq(3*(j-1)+1) < 0.d0) qqq(3*(j-1)+1)=qqq(3*(j-1)+1)+1.d0
       if(qqq(3*(j-1)+2) < 0.d0) qqq(3*(j-1)+2)=qqq(3*(j-1)+2)+1.d0
       if(qqq(3*(j-1)+3) < 0.d0) qqq(3*(j-1)+3)=qqq(3*(j-1)+3)+1.d0
       enddo
       ish=ndeg-6
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
       call get_extension(a1,a2,a3,rmax,ncext)
       n1=ncext(1) ; n2=ncext(2) ; n3=ncext(3)
       nb=(n1*n2*n3)*natot
       allocate(itypeext(nb))
       allocate(qext(3*nb+6))
       aa1=a1*dble(n1) ; aa2=a2*dble(n2) ; aa3=a3*dble(n3)
       nb=0
       do m=1,natot
       do i=0,n1-1
       do j=0,n2-1
       do k=0,n3-1
       nb=nb+1
       itypeext(nb)=itype(m)
       qext(3*(nb-1)+1)=qqq(3*(m-1)+1)+dble(i)
       qext(3*(nb-1)+2)=qqq(3*(m-1)+2)+dble(j)
       qext(3*(nb-1)+3)=qqq(3*(m-1)+3)+dble(k)
       enddo
       enddo
       enddo
       enddo
       do i=1,nb
       qext(3*(i-1)+1)=qext(3*(i-1)+1)/dble(n1)
       qext(3*(i-1)+2)=qext(3*(i-1)+2)/dble(n2)
       qext(3*(i-1)+3)=qext(3*(i-1)+3)/dble(n3)
       enddo
       do i=1,nb
       iti=itypeext(i)
       do j=1,nb
       if(j == i) cycle
       itj=itypeext(j)
       d1=qext(3*(i-1)+1)-qext(3*(j-1)+1)
       d2=qext(3*(i-1)+2)-qext(3*(j-1)+2)
       d3=qext(3*(i-1)+3)-qext(3*(j-1)+3)
       d1=d1-anint(d1)
       d2=d2-anint(d2)
       d3=d3-anint(d3)
       x=d1*aa1(1)+d2*aa2(1)+d3*aa3(1)
       y=d1*aa1(2)+d2*aa2(2)+d3*aa3(2)
       z=d1*aa1(3)+d2*aa2(3)+d3*aa3(3)
       r=sqrt(x*x+y*y+z*z)
       if(r > 6.d0*sigmamatrix(iti,itj)) cycle
       call xyz2rtp(x,y,z,r,theta,phi)
       do l=0,10,2
       do m=-l,l
       call sphhar(l,m,theta,phi,ylm)
       arg=-(r-2.d0*sigmamatrix(iti,itj))/2.d0
       if(arg < -50.d0) arg=-50.d0 ; if(arg >  50.d0) arg= 50.d0
       ctdangl(iti,itj,l,m)=ctdangl(iti,itj,l,m)+ylm*exp(arg)
       enddo
       enddo
       ndum(iti,itj)=ndum(iti,itj)+1
       enddo
       enddo
                     endif
       if(.not. lpbc)then
       do i=1,natot
       iti=itype(i)
       do j=1,natot
       if(j == i) cycle
       itj=itype(j)
       x=qqq(3*(i-1)+1)-qqq(3*(j-1)+1)
       y=qqq(3*(i-1)+2)-qqq(3*(j-1)+2)
       z=qqq(3*(i-1)+3)-qqq(3*(j-1)+3)
       r=sqrt(x*x+y*y+z*z)
       if(r > 6.d0*sigmamatrix(iti,itj)) cycle
       call xyz2rtp(x,y,z,r,theta,phi)
       do l=0,10,2
       do m=-l,l
       call sphhar(l,m,theta,phi,ylm)
       arg=-(r-2.d0*sigmamatrix(iti,itj))/2.d0
       if(arg < -50.d0) arg=-50.d0 ; if(arg >  50.d0) arg= 50.d0
       ctdangl(iti,itj,l,m)=ctdangl(iti,itj,l,m)+ylm*exp(arg)
       enddo
       enddo
       ndum(iti,itj)=ndum(iti,itj)+1
       enddo
       enddo
                     endif
       do i=1,nspecies
       do j=1,nspecies
       if(ndum(i,j) > 0) ctdangl(i,j,:,:)=ctdangl(i,j,:,:)/dble(ndum(i,j))
       enddo
       enddo
!
       pi=4.d0*atan(1.d0)
       qlab(:,:,:)=0.d0
       do i=1,nspecies
       do j=1,nspecies
       ctmp(:,:)=ctdangl(i,j,:,:)
       do l=0,10,2
       cl=cmplx(0.d0,0.d0)
       do m=-l,l
       cl(l)=cl(l)+ctmp(l,m)*conjg(ctmp(l,m))
       enddo
       cl(l)=sqrt((4.d0*pi)*cl(l)/(2.d0*dble(l)+1.d0))
       qlab(i,j,l)=real(cl(l))
       enddo
       enddo
       enddo
       if(allocated(itypeext)) deallocate(itypeext)
       if(allocated(qext)) deallocate(qext)
       end subroutine get_qlab

       end module qlabmod
!234567890
!      Written by In-Ho Lee, KRISS, April 21, 2016.
!      awk '{if(NR ==1) {print "12.34", $0} else { print }}' POSCAR
!      ifort -o batch4cdiffqlab.x strings.f90 numeral.f batch4cdiffqlab.f90 sphhar.f90
!      ./batch4cdiffqlab.x 
!      awk '{if(NF <3) {print}}' g1>g2
!      gnuplot> set pm3d
!      gnuplot> set palette rgbformulae 33,13,10
!      gnuplot> splot "matqlab.dat" with pm3d
       program batch4cdiffqlab
       USE qlabmod, ONLY : nstrc,qlab,qlabtest,qlabsave
       USE qlabmod, ONLY : qlab_cmp,qlab_init,qlab_final,get_qlab
       implicit none
       integer nspecies1,ndeg1,isize,npop
       integer i,j
       real*8 rmax00,tmp,omega,pi
       real*8, allocatable :: zmat(:,:)
       character*2, allocatable :: symbl1(:)
       integer, allocatable :: nelements1(:)
       real*8, allocatable :: qqq1(:),sigmamatrix1(:,:)
       character*280 fname,filea
       character*280 cmd
       logical lpbc0

       open(1,file='csa.in',form='formatted')
       read(1,*) nspecies1
       allocate(symbl1(nspecies1)) ; allocate(nelements1(nspecies1))
       allocate(sigmamatrix1(nspecies1,nspecies1))
       read(1,*) (symbl1(i),i=1,nspecies1)
       do i=1,nspecies1
       symbl1(i)=trim(adjustl(symbl1(i)))
       enddo
       read(1,*) (nelements1(i),i=1,nspecies1)
!      write(6,*) nelements1
       read(1,*) omega
       read(1,*)
       read(1,*)
       read(1,*)
       do i=1,nspecies1
       read(1,*) (sigmamatrix1(i,j),j=1,nspecies1)
!      write(6,*) (sigmamatrix1(i,j),j=1,nspecies1)
       enddo
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*) npop
       close(1)
       pi=4.d0*atan(1.d0)
       tmp=((3.d0*omega)/(4.d0*pi))**(1.d0/3.d0)
       write(6,*) omega,tmp
       nstrc=npop
!      write(6,*) nstrc
       ndeg1=6+3*sum(nelements1)
       allocate(qqq1(ndeg1))
       lpbc0=.true.
       isize=4
       do i=1,nstrc
       call xnumeral(i,fname,isize) ; fname='POSCAR_'//trim(fname) ; fname=trim(fname)
       cmd='cp '//trim(fname)//' POSCAR_a'  ; cmd=trim(cmd) ; call system(cmd)
       cmd="echo 5.00 >z1; head -n 1 POSCAR_a >>z1; awk 'ORS=NR%2?FS:RS' z1>z2; tail -n +2 POSCAR_a >>z2 ; mv z2 POSCAR_a ; rm z1"
       call system(cmd)
       filea='POSCAR_a'
       call get_poscar(filea,nspecies1,nelements1,symbl1,ndeg1,qqq1,rmax00)
       if(i == 1) call qlab_init(nspecies1,nelements1,sigmamatrix1,rmax00,lpbc0)
       call get_qlab(qqq1) ; qlabsave(:,:,:,i)=qlab(:,:,:)
       enddo
       allocate(zmat(nstrc,nstrc)) ; zmat=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j > i)then
       qlabtest(:,:,:)=qlabsave(:,:,:,j)-qlabsave(:,:,:,i)
       call qlab_cmp(tmp)
       zmat(i,j)=tmp ; zmat(j,i)=zmat(i,j)
                endif
       enddo
       enddo
       call qlab_final()
       fname='matqlab.dat'
       open(11,file=trim(fname),form='formatted')
       do i=1,nstrc
       do j=1,nstrc
       write(11,'(2i5,1x,f20.8)') i,j,zmat(i,j)
       enddo
       write(11,*)
       enddo
       close(11)
       deallocate(symbl1) ; deallocate(nelements1) ; deallocate(qqq1)
       call stats(nstrc,zmat)
       deallocate(zmat)
       deallocate(sigmamatrix1)
       stop
       end program batch4cdiffqlab
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine stats(nstrc,rmat)
       implicit none
       integer nstrc
       real*8 rmat(nstrc,nstrc)
       integer i,j,ip,iq,ir,nr
       real*8 rr,r0,r1,dr,tmr,avg,sig,xnorm
       real*8 rmax
       real*8, allocatable :: histo(:)
       real*8, external :: stepft

       rmax=maxval(rmat)
       nr=501 ; r1=rmax*1.1d0 ; r0=0.0d0 ; dr=(r1-r0)/dble(nr-1)
       allocate(histo(nr)) ; histo=0.d0
       avg=0.d0
       xnorm=0.0d0
       do i=1,nstrc
       do j=1,nstrc
       if(j <= i) cycle
       tmr=rmat(i,j)
       avg=avg+tmr
       xnorm=xnorm+1.0d0
       ir=(1+int(tmr/dr))-2 ; ip=max(ir,2) ; ir=(1+int(tmr/dr))+2 ; iq=min(ir,nr-1)
       do ir=ip,iq
       rr=r0+dr*float(ir-1)
       histo(ir)=histo(ir)+stepft(tmr-rr)*stepft(rr+dr-tmr)
       enddo 
       enddo
       enddo
       avg=avg/xnorm
       sig=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j <= i) cycle
       tmr=rmat(i,j)
       sig=sig+(tmr-avg)**2
       enddo
       enddo
       sig=sig/xnorm
       sig=sqrt(sig)
       open(11,file='matqlab_hist.dat',form='formatted')
       write(11,'(a1,2x,2f22.10)') '#', avg,sig
       do ir=1,nr
       rr=r0+dr*float(ir-1)
       if(rr <= rmax) write(11,'(f16.8,f20.10)') rr,histo(ir)
       enddo
       close(11)
       deallocate(histo)
       end 
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine get_extension(a1,a2,a3,rmax00,ncext)
       implicit none
       real*8 a1(3),a2(3),a3(3),rmax00
       integer ncext(3)
       real*8 v(3),h(3)

       call cross3(a1,a2,v) ; v=v/sqrt(sum(v*v)) ; h(3)=abs(sum(v*a3))
       call cross3(a3,a1,v) ; v=v/sqrt(sum(v*v)) ; h(2)=abs(sum(v*a2))
       call cross3(a2,a3,v) ; v=v/sqrt(sum(v*v)) ; h(1)=abs(sum(v*a1))
       v=rmax00/h+0.5d0
       ncext=nint(v)
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine cross3(a,b,c)
       implicit none
       real*8 a(3),b(3),c(3)

       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       end
!234567890
!      Written by In-Ho Lee, KRISS, March 27, 2016
       subroutine tolaty(xyz,a1,a2,a3,ktot)
       implicit none
       integer ktot
       real*8 xyz(ktot,3),a1(3),a2(3),a3(3)
       real*8 b(3,3),devid
       integer j,i
       real*8 d1,d2,d3

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
       xyz(j,1)=d1
       xyz(j,2)=d2
       xyz(j,3)=d3
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
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine latmat(rlat,wmat,ksign)
       implicit none
       integer ksign
       real*8 rlat(6),wmat(3,3)
       real*8 ra,rb,rc,cosinea,cosineb,cosinec
       real*8 epslat,tmr
       integer i,j

       epslat=1.0d-6
       if(ksign == 1)then
       wmat=0.0d0
       wmat(1,1)=rlat(1)
       wmat(2,1)=rlat(2)*cos(rlat(6))
       wmat(2,2)=rlat(2)*sin(rlat(6))
       wmat(3,1)=rlat(3)*cos(rlat(5))
       wmat(3,2)=rlat(3)*cos(rlat(4))*sin(rlat(6)) &
        -((rlat(3)*cos(rlat(5))-rlat(3)*cos(rlat(4))*cos(rlat(6)))/tan(rlat(6)))
       tmr=rlat(3)**2-wmat(3,1)**2-wmat(3,2)**2  ; if(tmr <= 1.0d-12) tmr=0.0d0
       wmat(3,3)=sqrt(tmr)
       do i=1,3
       do j=1,3
       if(abs(wmat(i,j)) < epslat) wmat(i,j)=0.0d0
       enddo
       enddo
                     else
       rlat=0.0d0
       ra=sqrt(wmat(1,1)**2+wmat(1,2)**2+wmat(1,3)**2)
       rb=sqrt(wmat(2,1)**2+wmat(2,2)**2+wmat(2,3)**2)
       rc=sqrt(wmat(3,1)**2+wmat(3,2)**2+wmat(3,3)**2)
       cosinea=(wmat(2,1)*wmat(3,1)+wmat(2,2)*wmat(3,2)+wmat(2,3)*wmat(3,3))/rb/rc
       cosineb=(wmat(1,1)*wmat(3,1)+wmat(1,2)*wmat(3,2)+wmat(1,3)*wmat(3,3))/rc/ra
       cosinec=(wmat(1,1)*wmat(2,1)+wmat(1,2)*wmat(2,2)+wmat(1,3)*wmat(2,3))/ra/rb
       rlat(1)=ra ; rlat(2)=rb ; rlat(3)=rc
       rlat(4)=acos(cosinea) ; rlat(5)=acos(cosineb) ; rlat(6)=acos(cosinec)
                     endif
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       real*8 function stepft(x)
       implicit none
       real*8 x
      
       stepft=0.d0 
       if(x >= 0.d0) stepft=1.d0 
       end function stepft
!234567890
!      Written by In-Ho Lee, KRISS, April 21, 2016.
       subroutine get_poscar(fname,nspecies,nelements,symbl,ndeg,qqq,rmax)
       USE strings, ONLY : parse,value
       implicit none
       character*280 fname
       integer nspecies,ndeg,nelements(nspecies)
       character*2 symbl(nspecies)
       real*8 rmax,qqq(ndeg)
       real*8 scale0,a1(3),a2(3),a3(3),tmp,cmatrix(3,3),t6(6),vtest
       real*8, allocatable :: dir(:,:),car(:,:)
       integer i,j,ks,ish
       integer ios,nargs
       character*200 str1
       character*200 args(40)
       character*20 delims
       logical lfault

       rmax=10.5d0
!      write(6,*) rmax,' rmax'
       lfault=.false.
       open(15,file=trim(fname),form='formatted')
       read(15,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs > 0)then
       call value(args(1),rmax,ios)
       if(ios /= 0)then
       rmax=10.5d0
       write(6,*) 'default rmax',rmax
                   endif
                    endif
       read(15,*) scale0
       read(15,*) a1(1),a1(2),a1(3)
       read(15,*) a2(1),a2(2),a2(3)
       read(15,*) a3(1),a3(2),a3(3)
       if(scale0 < 0.d0)then
       cmatrix(1,:)=a1(:) ; cmatrix(2,:)=a2(:) ; cmatrix(3,:)=a3(:)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       vtest=abs(scale0)/vtest ; vtest=vtest**(1.d0/3.d0)
       cmatrix=cmatrix*vtest
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
                        endif
       if(scale0 > 0.d0)then
       a1=a1*scale0 ; a2=a2*scale0 ; a3=a3*scale0
                        endif
       cmatrix(1,:)=a1(:) ; cmatrix(2,:)=a2(:) ; cmatrix(3,:)=a3(:)
       read(15,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       do i=1,nspecies
       symbl(i)=trim(adjustl(args(i)))
       enddo
       read(15,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       do i=1,nspecies
       call value(args(i),nelements(i),ios)
!      write(6,*) nelements(i)
       enddo
       ks=sum(nelements)
       allocate(dir(ks,3)) 
       read(15,'(a200)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs > 0)then
       if(args(1) == 'DIR' .or. args(1) == 'dir' .or. args(1) == 'D' .or. args(1) == 'd' .or. &
          args(1) == 'direct' .or. args(1) == 'Direct' .or. args(1) == 'Dir')then
       do j=1,ks
       read(15,*) dir(j,1),dir(j,2),dir(j,3)
       dir(j,1)=dir(j,1)-anint(dir(j,1))
       dir(j,2)=dir(j,2)-anint(dir(j,2))
       dir(j,3)=dir(j,3)-anint(dir(j,3))
       if(dir(j,1) < 0.d0) dir(j,1)=dir(j,1)+1.d0
       if(dir(j,2) < 0.d0) dir(j,2)=dir(j,2)+1.d0
       if(dir(j,3) < 0.d0) dir(j,3)=dir(j,3)+1.d0
       enddo
                                                                             else
       allocate(car(ks,3))
       do j=1,ks
       read(15,*) car(j,1),car(j,2),car(j,3)
       enddo
       dir=car
       call tolaty(dir,a1,a2,a3,ks)
       if(allocated(car)) deallocate(car)
                                                                             endif
                    endif
       goto 999
  911  continue
       lfault=.true.
  999  continue
       close(15)
       do i=1,ks
       qqq(3*(i-1)+1)=dir(i,1) ; qqq(3*(i-1)+2)=dir(i,2) ; qqq(3*(i-1)+3)=dir(i,3)
       enddo
       call latmat(t6,cmatrix,0)
       ndeg=3*sum(nelements)+6
       ish=ndeg-6
       do i=1,6
       qqq(ish+i)=t6(i)
       enddo
       if(allocated(dir)) deallocate(dir)
       end
!234567890
