!234567890
!      Written by In-Ho Lee, KRISS, April 21, 2016.
!      awk '{if(NR ==1) {print "12.34", $0} else { print }}' POSCAR
!      ifort -o batch4cdiffold.x strings.f90 numeral.f batch4cdiffold.f90
!      ./batch4cdiffold.x 
!      awk '{if(NF <3) {print}}' g1>g2
!      gnuplot> set pm3d
!      gnuplot> set palette rgbformulae 33,13,10
!      gnuplot> splot "mat.dat" with pm3d
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module bldist
       implicit none
       private
       save
       integer natot,ndeg,natotext1,natotext2,i222,ndim,nstrc
       real*8 rmax
       logical lpbc
       real*8, allocatable :: qext(:),wrk44(:)
       real*8, allocatable :: blsave(:,:),bltest(:),blsrtd(:)
       integer, allocatable :: iwrk44(:)
       public :: bldist_init,bldist_final,bldist_cmp,get_blsrtd
       public :: blsave,bltest,blsrtd

       contains
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine bldist_init(rmax00,i2220,nelements0,nspecies0,lpbc0,nstrc0)
       implicit none
       integer nstrc0,i2220,nspecies0,nelements0(nspecies0)
       real*8 rmax00
       logical lpbc0
       integer j

       i222=i2220 ; rmax=rmax00 ; lpbc=lpbc0 ; nstrc=nstrc0
       natot=sum(nelements0) ; ndeg=6+3*natot
       j=200*natot
       ndim=(j*(j-1))/2
       allocate(blsave(ndim,nstrc))
       allocate(bltest(ndim)) ; allocate(blsrtd(ndim))
       allocate(wrk44(ndim)) ; allocate(iwrk44(ndim))
       end subroutine bldist_init
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine bldist_final()
       implicit none

       deallocate(wrk44) ; deallocate(iwrk44)
       deallocate(blsave) ; deallocate(bltest) ; deallocate(blsrtd)
       end subroutine bldist_final
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine bldist_cmp(dista)
       implicit none
       real*8 dista
       integer i

       dista=0.d0 
       do i=1,ndim
       dista=dista+abs(bltest(i))
       enddo
       dista=sqrt(dista)
       end subroutine bldist_cmp
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine get_blsrtd(qqq)
       implicit none
       real*8 qqq(ndeg)
       real*8 x,y,z,d1,d2,d3
       real*8 aa1(3),aa2(3),aa3(3),t6(6),cmatrix(3,3),a1(3),a2(3),a3(3)
       integer ish,ncext(3),n1,n2,n3,natotext,ij,i,j,k,m

       if(lpbc)then
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
       if(n1 <=0 .or. n2 <=0 .or. n3 <=0)then
       n1=2 ; n2=2 ; n3=2
                                         endif
       if(i222 ==1)then
       n1=2 ; n2=2 ; n3=2
                   endif
       j=(n1*n2*n3)*natot
       allocate(qext(3*j+6))
       aa1=a1*dble(n1) ; aa2=a2*dble(n2) ; aa3=a3*dble(n3)
       natotext=0
       do m=1,natot
       do i=0,n1-1
       do j=0,n2-1
       do k=0,n3-1
       natotext=natotext+1
       qext(3*(natotext-1)+1)=qqq(3*(m-1)+1)+dble(i)
       qext(3*(natotext-1)+2)=qqq(3*(m-1)+2)+dble(j)
       qext(3*(natotext-1)+3)=qqq(3*(m-1)+3)+dble(k)
       enddo
       enddo
       enddo
       enddo
       do i=1,natotext
       qext(3*(i-1)+1)=qext(3*(i-1)+1)/dble(n1)
       qext(3*(i-1)+2)=qext(3*(i-1)+2)/dble(n2)
       qext(3*(i-1)+3)=qext(3*(i-1)+3)/dble(n3)
       enddo
       ij=0
       do i=1,natotext-1
       do j=i+1,natotext
       d1=qext(3*(i-1)+1)-qext(3*(j-1)+1)
       d2=qext(3*(i-1)+2)-qext(3*(j-1)+2)
       d3=qext(3*(i-1)+3)-qext(3*(j-1)+3)
       d1=d1-anint(d1)
       d2=d2-anint(d2)
       d3=d3-anint(d3)
       x=d1*aa1(1)+d2*aa2(1)+d3*aa3(1)
       y=d1*aa1(2)+d2*aa2(2)+d3*aa3(2)
       z=d1*aa1(3)+d2*aa2(3)+d3*aa3(3)
       ij=ij+1
       if(ij <= ndim) blsrtd(ij)=sqrt(x*x+y*y+z*z)
       enddo
       enddo
               endif
       if(.not. lpbc)then
       natotext=natot
       ij=0
       do i=1,natot-1
       do j=i+1,natot
       x=qqq(3*(i-1)+1)-qqq(3*(j-1)+1)
       y=qqq(3*(i-1)+2)-qqq(3*(j-1)+2)
       z=qqq(3*(i-1)+3)-qqq(3*(j-1)+3)
       ij=ij+1
       if(ij <= ndim) blsrtd(ij)=sqrt(x*x+y*y+z*z)
       enddo
       enddo
                     endif
       ij=min(ij,ndim)
       do i=1,ij
       wrk44(i)=blsrtd(i)
       enddo
       blsrtd=0.d0
       call sortnr(ij,wrk44,iwrk44)
       do i=1,ij
       blsrtd(i)=wrk44(iwrk44(i))
       if(blsrtd(i) > rmax) blsrtd(i)=0.0d0
       enddo
       if(allocated(qext)) deallocate(qext)
       end subroutine get_blsrtd

       end module bldist
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program batch4cdiffold
       USE bldist, ONLY : bldist_init,bldist_final,bldist_cmp,get_blsrtd
       USE bldist, ONLY : blsave,bltest,blsrtd
       implicit none
       integer nspecies1,ndeg1,npop
       integer i,j,i2220,nstrc0,isize
       real*8 rmax00,omega,tmp,pi
       logical lpbc0
       real*8, allocatable :: zmat(:,:)
       character*2, allocatable :: symbl1(:)
       integer, allocatable :: nelements1(:)
       real*8, allocatable :: qqq1(:)
       character*280 fname,filea
       character*280 cmd

       open(1,file='csa.in',form='formatted')
       read(1,*) nspecies1
       allocate(symbl1(nspecies1)) ; allocate(nelements1(nspecies1))
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
       read(1,*)
       enddo
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*) npop
       close(1)
!      write(6,*) npop
       nstrc0=npop
       pi=4.d0*atan(1.d0)
       tmp=((3.d0*omega)/(4.d0*pi))**(1.d0/3.d0)
       write(6,*) omega,tmp
       allocate(zmat(npop,npop))
       ndeg1=6+3*sum(nelements1)
       allocate(qqq1(ndeg1))
!
       lpbc0=.true.
       i2220=1
       i2220=0
       isize=4
       zmat=0.0d0
       do i=1,npop
       call xnumeral(i,fname,isize) ; fname='POSCAR_'//trim(fname) ; fname=trim(fname)
       cmd='cp '//trim(fname)//' POSCAR_a'  ; cmd=trim(cmd) ; call system(cmd)
       cmd="echo 10.00 >z1; head -n 1 POSCAR_a >>z1; awk 'ORS=NR%2?FS:RS' z1>z2; tail -n +2 POSCAR_a >>z2 ; mv z2 POSCAR_a ; rm z1"
       call system(cmd)
       filea='POSCAR_a'
       call get_poscar(filea,nspecies1,nelements1,symbl1,ndeg1,qqq1,rmax00)
       if(i == 1) call bldist_init(rmax00,i2220,nelements1,nspecies1,lpbc0,nstrc0)
       call get_blsrtd(qqq1) ; blsave(:,i)=blsrtd(:)
       enddo
       do i=1,npop
       do j=1,npop
       if(j > i)then
       bltest(:)=blsave(:,j)-blsave(:,i)
       call bldist_cmp(tmp)
       zmat(i,j)=tmp ; zmat(j,i)=zmat(i,j)
                endif
       enddo
       enddo
       call bldist_final()
       fname='matold.dat' 
       open(11,file=trim(fname),form='formatted')
       do i=1,npop
       do j=1,npop
       write(11,'(2i5,1x,f20.8)') i,j,zmat(i,j)
       enddo
       write(11,*)
       enddo
       close(11)
       deallocate(symbl1) ; deallocate(nelements1) ; deallocate(qqq1)
       call stats(npop,fname)
       deallocate(zmat)
       stop
       end program batch4cdiffold
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine stats(npop,fname)
       implicit none
       integer npop
       character*280 fname
       integer i,j,ip,iq,ir,nr
       real*8 rr,r0,r1,dr,tmr,avg,sig,xnorm
       real*8 rmax
       real*8, allocatable :: rmat(:,:),histo(:)
       real*8, external :: stepft

       allocate(rmat(npop,npop))
       open(1,file=trim(fname),form='formatted')
       do i=1,npop
       do j=1,npop
       read(1,*) ip,iq,rmat(i,j)
       if(i /= ip)then
                  stop
                  endif
       if(j /= iq)then
                  stop
                  endif
       enddo
       read(1,*)
       enddo
       close(1)
       rmax=maxval(rmat)
       nr=501 ; r1=rmax*1.1d0 ; r0=0.0d0 ; dr=(r1-r0)/dble(nr-1)
       allocate(histo(nr)) ; histo=0.d0
       avg=0.d0
       xnorm=0.0d0
       do i=1,npop
       do j=1,npop
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
       do i=1,npop
       do j=1,npop
       if(j <= i) cycle
       tmr=rmat(i,j)
       sig=sig+(tmr-avg)**2
       enddo
       enddo
       sig=sig/xnorm
       sig=sqrt(sig)
       open(11,file='matold_hist.dat',form='formatted')
       write(11,'(a1,2x,2f22.10)') '#', avg,sig
       do ir=1,nr
       rr=r0+dr*float(ir-1)
       if(rr <= rmax) write(11,'(f12.6,f22.12)') rr,histo(ir)
       enddo
       close(11)
       deallocate(rmat)
       deallocate(histo)
       end 
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       real*8 function stepft(x)
       implicit none
       real*8 x
      
       stepft=0.0d0 
       if(x >= 0.0d0) stepft=1.0d0 
       end function stepft
!234567890
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
