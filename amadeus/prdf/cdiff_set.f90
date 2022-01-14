!234567890
!      Written by In-Ho Lee, KRISS, April 29, 2016.
       module elemdist
       implicit none
       private
       save
       integer ndeg,natot,nspecies,nstrc,ndim
       real*8 rmax,rc1,rc2
       logical lpbc
       integer, allocatable :: nelements(:)
       real*8, allocatable :: hist1(:),histsave1(:,:),histtest1(:),qext(:),engg(:)
       real*8, allocatable :: hist2(:),histsave2(:,:),histtest2(:)
       public :: elemdist_init,elemdist_final,get_hist,elemdist_cmp
       public :: engg
       contains
!234567890
!      Written by In-Ho Lee, KRISS, April 29, 2016.
       subroutine elemdist_init(rmax0,nspecies0,nelements0,rc10,rc20,nstrc0,lpbc0)
       implicit none
       integer nspecies0,nelements0(nspecies0),nstrc0
       real*8 rmax0,rc10,rc20
       logical lpbc0
       integer i

       lpbc=lpbc0 ; rmax=rmax0
       if(rmax <= 0.d0) rmax=10.d0
       rc1=rc10 ; rc2=rc20
       nstrc=nstrc0
       nspecies=nspecies0
       allocate(engg(nstrc))
       allocate(nelements(nspecies))
       do i=1,nspecies
       nelements(i)=nelements0(i)
       enddo
       natot=sum(nelements) ; ndeg=3*natot+6
       ndim=natot*8
       allocate(hist1(0:ndim),histsave1(0:ndim,nstrc),histtest1(0:ndim))
       allocate(hist2(0:ndim),histsave2(0:ndim,nstrc),histtest2(0:ndim))
       end subroutine elemdist_init
!234567890
!      Written by In-Ho Lee, KRISS, April 29, 2016.
       subroutine elemdist_final(kmd,avg,sig)
       implicit none
       integer kmd
       real*8 avg,sig
       integer i,j
       character*280 fname
       real*8, allocatable :: wmat(:,:)

       if(kmd >= 1)then
       allocate(wmat(nstrc,nstrc)) ; wmat=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j > i)then
       call elemdist_cmp(i,j,wmat(i,j)) ; wmat(j,i)=wmat(i,j)
                endif
       enddo
       enddo
       call gen_avgsig(avg,sig,nstrc,wmat)
       if(kmd > 1)then
       fname='mat1.dat'         ; call plotdiff(fname,nstrc,wmat)
       fname='mat1_hist.dat'    ; call stats(fname,nstrc,wmat)
                  endif
       if(kmd > 1)then
       fname='ed1.dat'          ; call plotdiffe(fname,nstrc,engg,wmat)
                  endif
       deallocate(wmat)
                   endif
       if(allocated(engg)) deallocate(engg)
       if(allocated(nelements)) deallocate(nelements)
       if(allocated(hist1)) deallocate(hist1)
       if(allocated(hist2)) deallocate(hist2)
       if(allocated(histsave1)) deallocate(histsave1)
       if(allocated(histsave2)) deallocate(histsave2)
       if(allocated(histtest1)) deallocate(histtest1)
       if(allocated(histtest2)) deallocate(histtest2)
       end subroutine elemdist_final
!234567890
!      Written by In-Ho Lee, KRISS, April 29, 2016.
       subroutine get_hist(iisv,qqq)
       implicit none
       integer iisv
       real*8 qqq(ndeg)
       integer kk1,kk2,ish,n1,n2,n3,natotext,m,i,j,k
       real*8 r,x,y,z,d1,d2,d3,t6(6),a1(3),a2(3),a3(3),cmatrix(3,3),aa1(3),aa2(3),aa3(3)

       hist1(:)=0.0d0 ; hist2(:)=0.0d0
       if(lpbc)then
       ish=ndeg-6
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       a1(:)=cmatrix(1,:) ; a2(:)=cmatrix(2,:) ; a3(:)=cmatrix(3,:)
!      call get_extension(a1,a2,a3,rmax,ncext)
!      ncext=ncext*2
!      n1=ncext(1) ; n2=ncext(2) ; n3=ncext(3)
       n1=2 ; n2=2 ; n3=2
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
       do i=1,natotext
       kk1=0 ; kk2=0
       do j=1,natotext
       if(j == i) cycle
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
       if( r <= rc1) then
       kk1=kk1+1
                     endif
       if( r > rc1 .and. r <= rc2) then
       kk2=kk2+1
                                   endif
       enddo
       if(kk1 > -1 .and. kk1 <= ndim) hist1(kk1)=hist1(kk1)+1.0d0
       if(kk2 > -1 .and. kk2 <= ndim) hist2(kk2)=hist2(kk2)+1.0d0
       enddo
       deallocate(qext)
               endif
       if(.not. lpbc)then
       do i=1,natot
       kk1=0 ; kk2=0
       do j=1,natot
       if(j == i) cycle
       x=qqq(3*(i-1)+1)-qqq(3*(j-1)+1)
       y=qqq(3*(i-1)+2)-qqq(3*(j-1)+2)
       z=qqq(3*(i-1)+3)-qqq(3*(j-1)+3)
       r=sqrt(x*x+y*y+z*z)
       if( r <= rc1) then
       kk1=kk1+1
                     endif
       if( r > rc1 .and. r <= rc2) then
       kk2=kk2+1
                                   endif
       enddo
       if(kk1 > -1 .and. kk1 <= ndim) hist1(kk1)=hist1(kk1)+1.0d0
       if(kk2 > -1 .and. kk2 <= ndim) hist2(kk2)=hist2(kk2)+1.0d0
       enddo
                     endif
       if(iisv /= 0)then
       histsave1(:,iisv)=hist1(:) 
       histsave2(:,iisv)=hist2(:)
                    endif
       end subroutine get_hist
!234567890
!      Written by In-Ho Lee, KRISS, April 29, 2016.
       subroutine elemdist_cmp(ii,jj,dista)
       implicit none
       integer ii,jj
       real*8 dista
       integer k
       
       if(ii == 0 .and. jj /= 0)then
       histtest1(:)=histsave1(:,jj)-hist1(:) ; histtest2(:)=histsave2(:,jj)-hist2(:) 
                                endif
       if(ii /= 0 .and. jj == 0)then
       histtest1(:)=hist1(:)-histsave1(:,ii) ; histtest2(:)=hist2(:)-histsave2(:,ii) 
                                endif
       if(ii /= 0 .and. jj /= 0)then
       histtest1(:)=histsave1(:,jj)-histsave1(:,ii) ; histtest2(:)=histsave2(:,jj)-histsave2(:,ii) 
                                endif

       dista=0.d0
       do k=0,natot
       dista=dista+dble(k)*(2.d0*abs(histtest1(k))+abs(histtest2(k)))
       enddo
       end subroutine elemdist_cmp
 
       end module elemdist
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       module bldist
       implicit none
       private
       save
       integer natot,ndeg,natotext1,natotext2,i222,ndim,nstrc
       real*8 rmax
       logical lpbc
       real*8, allocatable :: qext(:),wrk44(:),engg(:)
       real*8, allocatable :: blsave(:,:),bltest(:),blsrtd(:)
       integer, allocatable :: iwrk44(:)
       public :: bldist_init,bldist_final,get_blsrtd,bldist_cmp
       public :: engg
       contains
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine bldist_init(rmax0,i2220,nelements0,nspecies0,nstrc0,lpbc0)
       implicit none
       integer nstrc0,i2220,nspecies0,nelements0(nspecies0)
       real*8 rmax0
       logical lpbc0
       integer j

       i222=i2220 ; rmax=rmax0 ; lpbc=lpbc0 ; nstrc=nstrc0
       if(rmax <= 0.d0) rmax=10.d0
       natot=sum(nelements0) ; ndeg=6+3*natot
       j=20*natot
       ndim=(j*(j-1))/2
       allocate(blsave(ndim,nstrc))
       allocate(bltest(ndim)) ; allocate(blsrtd(ndim))
       allocate(wrk44(ndim)) ; allocate(iwrk44(ndim))
       allocate(engg(nstrc))
       end subroutine bldist_init
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine bldist_final(kmd,avg,sig)
       implicit none
       integer kmd
       real*8 avg,sig
       integer i,j
       character*280 fname
       real*8, allocatable :: wmat(:,:)

       if(kmd >= 1)then
       allocate(wmat(nstrc,nstrc)) ; wmat=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j > i)then
       call bldist_cmp(i,j,wmat(i,j)) ; wmat(j,i)=wmat(i,j)
                endif
       enddo
       enddo
       call gen_avgsig(avg,sig,nstrc,wmat)
       if(kmd > 1)then
       if(i222 == 1)then
       fname='mat2.dat'       ; call plotdiff(fname,nstrc,wmat)
       fname='mat2_hist.dat'  ; call stats(fname,nstrc,wmat)
                    endif
       if(i222 == 2)then
       fname='mat3.dat'       ; call plotdiff(fname,nstrc,wmat)
       fname='mat3_hist.dat'  ; call stats(fname,nstrc,wmat)
                    endif
                  endif
       if(kmd > 1)then
       if(i222 == 1)then
       fname='ed2.dat'          ; call plotdiffe(fname,nstrc,engg,wmat)
                    endif
       if(i222 == 2)then
       fname='ed3.dat'          ; call plotdiffe(fname,nstrc,engg,wmat)
                    endif
                  endif
       deallocate(wmat)
                   endif
       if(allocated(engg)) deallocate(engg)
       if(allocated(wrk44)) deallocate(wrk44) 
       if(allocated(iwrk44)) deallocate(iwrk44)
       if(allocated(blsave)) deallocate(blsave) 
       if(allocated(bltest)) deallocate(bltest) 
       if(allocated(blsrtd)) deallocate(blsrtd)
       end subroutine bldist_final
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine bldist_cmp(ii,jj,dista)
       implicit none
       integer ii,jj
       real*8 dista
       integer k

       if(ii /= 0 .and. jj /= 0)then
       bltest(:)=blsave(:,jj)-blsave(:,ii)             
                                endif
       if(ii == 0 .and. jj /= 0)then
       bltest(:)=blsave(:,jj)-blsrtd(:)             
                                endif
       if(ii /= 0 .and. jj == 0)then
       bltest(:)=blsrtd(:)-blsave(:,ii)             
                                endif

       dista=0.d0 
       do k=1,ndim
       dista=dista+abs(bltest(k))
       enddo
       dista=sqrt(dista)
       end subroutine bldist_cmp
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine get_blsrtd(iisv,qqq)
       implicit none
       integer iisv
       real*8 qqq(ndeg)
       real*8 r,x,y,z,d1,d2,d3
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
       ncext=ncext*2
       n1=ncext(1) ; n2=ncext(2) ; n3=ncext(3)
       if(i222 >  6)then
       n1=i222 ; n2=i222 ; n3=i222
                    endif
       if(i222 == 5)then
       n1=5 ; n2=5 ; n3=5
                    endif
       if(i222 == 4)then
       n1=4 ; n2=4 ; n3=4
                    endif
       if(i222 == 3)then
       n1=3 ; n2=3 ; n3=3
                    endif
       if(i222 == 2)then
       n1=2 ; n2=2 ; n3=2
                    endif
       if(i222 == 1)then
       n1=1 ; n2=1 ; n3=1
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
       r=sqrt(x*x+y*y+z*z)
       if(r < rmax)then
       ij=ij+1
       if(ij <= ndim) blsrtd(ij)=r
                   endif
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
       r=sqrt(x*x+y*y+z*z)
       if(r < rmax)then
       ij=ij+1
       if(ij <= ndim) blsrtd(ij)=r
                   endif
       enddo
       enddo
                     endif
       ij=min(ij,ndim)
       do i=1,ij
       wrk44(i)=blsrtd(i)
       enddo
       call sortnr(ij,wrk44,iwrk44)
       do i=1,ij
       blsrtd(i)=wrk44(iwrk44(i))
       enddo
       do i=ij+1,ndim
       blsrtd(i)=0.d0
       enddo
       if(allocated(qext)) deallocate(qext)
       if(iisv /= 0)then
       blsave(:,iisv)=blsrtd(:)
                    endif
       end subroutine get_blsrtd

       end module bldist
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
       real*8, allocatable :: sigmamatrix(:,:),engg(:)
       complex*16, allocatable :: ctdangl(:,:,:,:)
       logical lpbc
       public :: qlab_init,qlab_final,get_qlab,qlab_cmp
       public :: engg
       contains
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine qlab_init(rmax0,nspecies0,nelements0,sigmamatrix0,nstrc0,lpbc0)
       implicit none
       integer nspecies0,nelements0(nspecies0),nstrc0
       real*8 sigmamatrix0(nspecies0,nspecies0),rmax0
       logical lpbc0
       integer i,j,k

       nstrc=nstrc0 ; rmax=rmax0 ; lpbc=lpbc0
       if(rmax <= 0.d0) rmax=10.d0
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
       allocate(engg(nstrc))
       end subroutine qlab_init
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine qlab_final(kmd,avg,sig)
       implicit none
       integer kmd
       real*8 avg,sig
       integer i,j
       character*280 fname
       real*8, allocatable :: wmat(:,:)

       if(kmd >= 1)then
       allocate(wmat(nstrc,nstrc)) ; wmat=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j > i)then
       call qlab_cmp(i,j,wmat(i,j))  ; wmat(j,i)=wmat(i,j)
                endif
       enddo
       enddo
       call gen_avgsig(avg,sig,nstrc,wmat)
       if(kmd > 1)then
       fname='mat4.dat'      ; call plotdiff(fname,nstrc,wmat)
       fname='mat4_hist.dat' ; call stats(fname,nstrc,wmat)
                  endif
       if(kmd > 1)then
       fname='ed4.dat'          ; call plotdiffe(fname,nstrc,engg,wmat)
                  endif
       deallocate(wmat)
                   endif
       if(allocated(engg)) deallocate(engg)
       if(allocated(sigmamatrix)) deallocate(sigmamatrix)
       if(allocated(itype)) deallocate(itype) 
       if(allocated(ctdangl)) deallocate(ctdangl) 
       if(allocated(ndum)) deallocate(ndum)
       if(allocated(qlab)) deallocate(qlab) 
       if(allocated(qlabsave)) deallocate(qlabsave) 
       if(allocated(qlabtest)) deallocate(qlabtest)
       end subroutine qlab_final
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine qlab_cmp(ii,jj,dista)
       implicit none
       integer ii,jj
       real*8 dista
       real*8 tmp,tmq
       integer i,j,l

       if(ii /= 0 .and. jj /= 0)then
       qlabtest(:,:,:)=qlabsave(:,:,:,jj)-qlabsave(:,:,:,ii) 
                                endif
       if(ii == 0 .and. jj /= 0)then
       qlabtest(:,:,:)=qlabsave(:,:,:,jj)-qlab(:,:,:) 
                                endif
       if(ii /= 0 .and. jj == 0)then
       qlabtest(:,:,:)=qlab(:,:,:)-qlabsave(:,:,:,ii) 
                                endif

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
       subroutine get_qlab(iisv,qqq)
       implicit none
       integer iisv
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
       ncext=ncext*2
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
       if(iisv /= 0)then
       qlabsave(:,:,:,iisv)=qlab(:,:,:)
                    endif
       end subroutine get_qlab

       end module qlabmod
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       module prdf
       implicit none
       private
       save
       integer nr,ndim,nstrc,ndeg,natot,kosine
       real*8 r0,r1,dr,rmax
       real*8, allocatable :: prdf0(:,:),prdfsave(:,:,:),prdftest(:,:),engg(:)
       integer, allocatable :: irow(:,:)
       integer nspecies
       integer, allocatable :: nelements(:)
       character*2, allocatable :: symbl(:)
       logical lpbc
       public :: prdf_init,prdf_final,get_prdf,prdf_cmp,stepft
       public :: engg
       contains
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine prdf_init(rmax0,nspecies0,nelements0,symbl0,nstrc0,lpbc0)
       implicit none
       integer nstrc0,nspecies0,nelements0(nspecies0)
       real*8 rmax0
       character*2 symbl0(nspecies0)
       logical lpbc0
       integer i,j,jprint

       kosine=0
       if(nstrc0 < 0)then
       nstrc0=iabs(nstrc0)
       kosine=1
                     endif
       lpbc=lpbc0 ; nstrc=nstrc0
       nspecies=nspecies0
       allocate(nelements(nspecies)) ; allocate(symbl(nspecies))
       allocate(irow(nspecies,nspecies))
       nelements=nelements0
       natot=sum(nelements) ; ndeg=3*natot+6
       ndim=0
       do i=1,nspecies
       do j=1,nspecies
       ndim=ndim+1
       irow(i,j)=ndim
       enddo
       enddo
       do i=1,nspecies
       symbl(i)=trim(adjustl(symbl0(i)))
       nelements(i)=nelements0(i)
       enddo
       if(rmax0 > 10.d0)then
       write(6,*) rmax0,' has been changed to',10.d0
       rmax0=10.d0
                        endif
       rmax=rmax0
       if(rmax <= 0.d0) rmax=10.d0
!      nr=1001 ; r0=0.0d0 ; r1=rmax*2.0d0 ; dr=(r1-r0)/dble(nr-1)
       nr=401  ; r0=0.0d0 ; r1=20.0d0    ; dr=(r1-r0)/dble(nr-1)
       jprint=1
       jprint=0
       if(jprint == 1)then
       write(6,'(2f18.8,1x,a7)') rmax,dr,'rmax,dr'
       write(6,'(3f18.8)') r0,r1,dr
       write(6,'(i5)') nr
                      endif
       if(.not. allocated(prdf0)) allocate(prdf0(nr,ndim))
       if(.not. allocated(prdfsave)) allocate(prdfsave(nr,ndim,nstrc))
       if(.not. allocated(prdftest)) allocate(prdftest(nr,ndim))
       allocate(engg(nstrc))
       end subroutine prdf_init
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine prdf_final(kmd,avg,sig)
       implicit none
       integer kmd
       real*8 avg,sig
       integer i,j
       character*280 fname
       real*8, allocatable :: wmat(:,:)

       if(kmd >= 1)then
       allocate(wmat(nstrc,nstrc)) ; wmat=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j > i)then
       call prdf_cmp(i,j,wmat(i,j))  ; wmat(j,i)=wmat(i,j)
                endif
       enddo
       enddo
       call gen_avgsig(avg,sig,nstrc,wmat)
       if(kmd > 1)then
       if(kosine == 1)then
       fname='mat6.dat'          ; call plotdiff(fname,nstrc,wmat)
       fname='mat6_hist.dat'     ; call stats(fname,nstrc,wmat)
                      else
       fname='mat5.dat'          ; call plotdiff(fname,nstrc,wmat)
       fname='mat5_hist.dat'     ; call stats(fname,nstrc,wmat)
                      endif
                  endif
       if(kmd > 1)then
       if(kosine == 1)then
       fname='ed6.dat'          ; call plotdiffe(fname,nstrc,engg,wmat)
                      else
       fname='ed5.dat'          ; call plotdiffe(fname,nstrc,engg,wmat)
                      endif
                  endif
       deallocate(wmat)
                   endif
       if(allocated(engg)) deallocate(engg)
       if(allocated(nelements)) deallocate(nelements)
       if(allocated(symbl)) deallocate(symbl)
       if(allocated(irow)) deallocate(irow)
       if(allocated(prdf0)) deallocate(prdf0)
       if(allocated(prdfsave)) deallocate(prdfsave)
       if(allocated(prdftest)) deallocate(prdftest)
       end subroutine prdf_final
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine prdf_cmp(ii,jj,v)
       implicit none
       integer ii,jj
       real*8 v
       integer ir
       real*8 rr

       if(kosine == 1)then
       call prdf_cmp0(ii,jj,v)
       return
                      endif

       if(ii /= 0 .and. jj == 0)then 
       prdftest(:,:)=prdf0(:,:)-prdfsave(:,:,ii)
                                endif
       if(ii == 0 .and. jj /= 0)then 
       prdftest(:,:)=prdfsave(:,:,jj)-prdf0(:,:)
                                endif
       if(ii /= 0 .and. jj /= 0)then 
       prdftest(:,:)=prdfsave(:,:,jj)-prdfsave(:,:,ii) 
                                endif
       do ir=1,nr
       rr=r0+dr*float(ir-1)
       if(rr > rmax) prdftest(ir,:)=0.d0
       enddo
       call frobeniusnorm(nr,ndim,prdftest,v)
!      write(6,*) v
       end subroutine prdf_cmp
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine prdf_cmp0(ii,jj,v)
       implicit none
       integer ii,jj
       real*8 v
       integer ir,i1,i2,i3
       real*8 rr,tmp,tmq,tmr,tms

       tmp=0.0d0 ; tmq=0.0d0 ; tmr=0.0d0
       do ir=1,nr
       rr=r0+dr*float(ir-1)
       if(rr < rmax)then
       do i1=1,nspecies
       do i2=1,nspecies
       if(i2 > i1) cycle
       tms=dble(nelements(i1))*dble(nelements(i2))
       i3=irow(i1,i2)
       if(ii /= 0 .and. jj /= 0)then 
       tmp=tmp+(prdfsave(ir,i3,jj)-1.d0)*(prdfsave(ir,i3,ii)-1.d0)*tms
       tmq=tmq+(prdfsave(ir,i3,ii)-1.d0)*(prdfsave(ir,i3,ii)-1.d0)*tms
       tmr=tmr+(prdfsave(ir,i3,jj)-1.d0)*(prdfsave(ir,i3,jj)-1.d0)*tms
                                endif
       if(ii == 0 .and. jj /= 0)then 
       tmp=tmp+(prdfsave(ir,i3,jj)-1.d0)*(prdf0(ir,i3)-1.d0)*tms
       tmq=tmq+(prdf0(ir,i3)-1.d0)*(prdf0(ir,i3)-1.d0)*tms
       tmr=tmr+(prdfsave(ir,i3,jj)-1.d0)*(prdfsave(ir,i3,jj)-1.d0)*tms
                                endif
       if(ii /= 0 .and. jj == 0)then 
       tmp=tmp+(prdf0(ir,i3)-1.d0)*(prdfsave(ir,i3,ii)-1.d0)*tms
       tmq=tmq+(prdfsave(ir,i3,ii)-1.d0)*(prdfsave(ir,i3,ii)-1.d0)*tms
       tmr=tmr+(prdf0(ir,i3)-1.d0)*(prdf0(ir,i3)-1.d0)*tms
                                endif
       enddo
       enddo
                    endif
       enddo
       tmq=sqrt(tmq) ; tmr=sqrt(tmr) ; v=0.5d0*(1.d0-tmp/(tmq*tmr))
       end subroutine prdf_cmp0
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       real*8 function stepft(x)
       implicit none
       real*8 x
      
       stepft=0.d0 
       if(x >= 0.d0) stepft=1.d0 
       end function stepft
!234567890
!      Written by In-Ho Lee, KRISS, April 15, 2016.
       subroutine frobeniusnorm(m,n,a,v)
       implicit none
       integer m,n
       real*8 a(m,n),v
       integer i,j

       v=0.d0
       do i=1,m
       do j=1,n
       v=v+(abs(a(i,j)))**2
       enddo
       enddo
       v=sqrt(v)
       end subroutine frobeniusnorm
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine get_prdf(iisv,qqq)
       implicit none
       integer iisv
       real*8 qqq(ndeg)
       integer, allocatable :: itype(:),nelementsext(:)
       real*8, allocatable :: dir(:,:),dirext(:,:)
       integer i1,i2,i3,j1,j2,k1,k2,ir,ip,iq,ks,ksext,ish,iprint,kprint
       integer i,j,k,ie,je,ke,k0,kk0,ncext(3)
       real*8 rr,pi,tmr,vtest,vtest0,vec(3),wec(3)
       real*8 a1(3),a2(3),a3(3),aa1(3),aa2(3),aa3(3),t6(6),cmatrix(3,3)

       pi=4.0d0*atan(1.0d0)
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
       ncext=ncext*2
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest0=abs(vtest)
       ks=sum(nelements)
       kprint=1
       kprint=0
       if(kprint == 1)then
       write(6,'(6f16.8,1x,f12.5)') t6(1),t6(2),t6(3),t6(4)*180.d0/pi,t6(5)*180.d0/pi,t6(6)*180.d0/pi,vtest0
       write(6,*) ncext
                      endif
       kprint=1
       kprint=0
!
       allocate(dir(natot,3))
       do i=1,natot
       dir(i,1)=qqq(3*(i-1)+1) ; dir(i,2)=qqq(3*(i-1)+2) ; dir(i,3)=qqq(3*(i-1)+3)
       enddo
       aa1=a1*dble(ncext(1)) ; aa2=a2*dble(ncext(2)) ; aa3=a3*dble(ncext(3))
       k=sum(nelements)*(ncext(1)*ncext(2)*ncext(3))
       allocate(nelementsext(nspecies)) ; allocate(itype(k))
       nelementsext(:)=nelements(:)*(ncext(1)*ncext(2)*ncext(3))
       ksext=sum(nelementsext)
!      write(6,*)k,ksext
       allocate(dirext(ksext,3))
       cmatrix(1,:)=aa1(:) ; cmatrix(2,:)=aa2(:) ; cmatrix(3,:)=aa3(:)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest=abs(vtest)
       k=0 ; kk0=0
       do i=1,nspecies
       k0=0
       do j=1,nelements(i)
       k0=k0+1 ; kk0=kk0+1
       do ie=0,ncext(1)-1
       do je=0,ncext(2)-1
       do ke=0,ncext(3)-1
       k=k+1
       dirext(k,1)=(dble(ie)+dir(kk0,1))/dble(ncext(1))
       dirext(k,2)=(dble(je)+dir(kk0,2))/dble(ncext(2))
       dirext(k,3)=(dble(ke)+dir(kk0,3))/dble(ncext(3))
       itype(k)=i
       enddo
       enddo
       enddo
       enddo
       if(k0 /= nelements(i))then
       write(6,*) 'error, k0',i
                             stop
                             endif
       enddo
       if(k /= ksext)then
       write(6,*) 'error ksext',k,ksext
                     stop
                     endif
       deallocate(dir) 
       i3=0
       do i1=1,nspecies
       do i2=1,nspecies
       i3=i3+1
       if(ndim < i3)then
       write(6,*) 'check ndim',ndim
                    stop
                    endif
       irow(i1,i2)=i3
       prdf0(:,i3)=0.d0
       enddo
       enddo
       k1=0
       do i1=1,nspecies
       do j1=1,nelementsext(i1)
       k1=k1+1
       if(k1 > ksext) stop
       k2=0
       do i2=1,nspecies
       i3=irow(i1,i2)
       do j2=1,nelementsext(i2)
       k2=k2+1
       if(k2 > ksext) stop
       if(k1 == k2) cycle
       if(itype(k1) == i1 .and. itype(k2) == i2)then
       vec(:)=dirext(k1,:)-dirext(k2,:) 
       vec(1)=vec(1)-anint(vec(1))
       vec(2)=vec(2)-anint(vec(2))
       vec(3)=vec(3)-anint(vec(3))
       wec(:)=vec(1)*aa1(:)+vec(2)*aa2(:)+vec(3)*aa3(:)
       tmr=sqrt(dot_product(wec,wec))
       ir=(1+int(tmr/dr))-2 ; ip=max(ir,2) ; ir=(1+int(tmr/dr))+2 ; iq=min(ir,nr-1)
       do ir=ip,iq
       rr=r0+dr*float(ir-1)
       prdf0(ir,i3)=prdf0(ir,i3)+stepft(tmr-rr)*stepft(rr+dr-tmr)*dble(ksext)/(dble(ksext)/vtest) &
                              /(dble(nelementsext(i1)))/(dble(nelementsext(i2)))
       enddo
                                                endif
       enddo
       enddo 
       enddo
       enddo
       deallocate(dirext) ; deallocate(nelementsext) ; deallocate(itype)
               endif
       if(.not. lpbc)then
       ish=ndeg-6
       do i=1,6
       t6(i)=qqq(ish+i)
       enddo
       call latmat(t6,cmatrix,1)
       ks=sum(nelements)
       vtest=(cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1) &
            +(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2) &
            +(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       vtest0=abs(vtest)
       allocate(itype(ks))
       k=0
       do i=1,nspecies
       do j=1,nelements(i)
       k=k+1
       itype(k)=i
       enddo
       enddo
       i3=0
       do i1=1,nspecies
       do i2=1,nspecies
       i3=i3+1
       if(ndim < i3)then
       write(6,*) 'check ndim',ndim
                    stop
                    endif
       irow(i1,i2)=i3
       prdf0(:,i3)=0.d0
       enddo
       enddo
       k1=0
       do i1=1,nspecies
       do j1=1,nelements(i1)
       k1=k1+1
       if(k1 > ks) stop
       k2=0
       do i2=1,nspecies
       i3=irow(i1,i2)
       do j2=1,nelements(i2)
       k2=k2+1
       if(k2 > ks) stop
       if(k1 == k2) cycle
       if(itype(k1) == i1 .and. itype(k2) == i2)then
       wec(1)=qqq(3*(k1-1)+1)-qqq(3*(k2-1)+1) 
       wec(2)=qqq(3*(k1-1)+2)-qqq(3*(k2-1)+2) 
       wec(3)=qqq(3*(k1-1)+3)-qqq(3*(k2-1)+3) 
       tmr=sqrt(dot_product(wec,wec))
       ir=(1+int(tmr/dr))-2 ; ip=max(ir,2) ; ir=(1+int(tmr/dr))+2 ; iq=min(ir,nr-1)
       do ir=ip,iq
       rr=r0+dr*float(ir-1)
       prdf0(ir,i3)=prdf0(ir,i3)+stepft(tmr-rr)*stepft(rr+dr-tmr)*dble(ks)/(dble(ks)/vtest0) &
                              /(dble(nelements(i1)))/(dble(nelements(i2)))
       enddo
                                                endif
       enddo
       enddo 
       enddo
       enddo
       deallocate(itype)
                     endif
       iprint=1
       iprint=0
       do i1=1,nspecies
       do i2=1,nspecies
       i3=irow(i1,i2)
       if(iprint == 1)then
       write(6,'(a1,1x,2i4)') '#', i1,i2
                      endif
       do ir=2,nr
       rr=r0+dr*float(ir-1)
       prdf0(ir,i3)=prdf0(ir,i3)/(dr*4.0d0*pi*rr**2)
       if(iprint == 1)then
       if(rr <= rmax) write(6,'(f12.6,f22.12)') rr,prdf0(ir,i3)
                      endif
       enddo
       if(iprint == 1)then
       write(6,*) '&'
                      endif
       enddo
       enddo
       if(iisv /= 0)then
       prdfsave(:,:,iisv)=prdf0(:,:)
                    endif
       end subroutine get_prdf

       end module prdf
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine plotdiffe(fname,nstrc,engg,zmat)
       implicit none
       character*80 fname
       integer nstrc
       real*8 zmat(nstrc,nstrc),engg(nstrc)
       integer i,j

       i=1
       open(11,file=trim(fname),form='formatted')
       write(11,'(2f20.8)') zmat(i,i),engg(i)-engg(i)
       do j=2,nstrc
       write(11,'(2f20.8)') zmat(i,j),engg(j)-engg(i)
       enddo
       close(11)
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine plotdiff(fname,nstrc,zmat)
       implicit none
       character*80 fname
       integer nstrc
       real*8 zmat(nstrc,nstrc)
       integer i,j

       open(11,file=trim(fname),form='formatted')
       do i=1,nstrc
       do j=1,nstrc
       write(11,'(2i5,1x,f20.8)') i,j,zmat(i,j)
       enddo
       write(11,*)
       enddo
       close(11)
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine stats(fname,nstrc,rmat)
       USE prdf, ONLY : stepft
       implicit none
       character*80 fname
       integer nstrc
       real*8 rmat(nstrc,nstrc)
       integer i,j,ip,iq,ir,nr
       real*8 rr,r0,r1,dr,tmr,avg,sig,xnorm
       real*8 rmax
       real*8, allocatable :: histo(:)

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
       open(11,file=trim(fname),form='formatted')
       write(11,'(a1,2x,2f22.10)') '#', avg,sig
       do ir=1,nr
       rr=r0+dr*float(ir-1)
       if(rr <= rmax) write(11,'(f12.6,f22.12)') rr,histo(ir)
       enddo
       close(11)
       deallocate(histo)
       end 
!234567890
!      Written by In-Ho Lee, KRISS, April 13, 2016.
       subroutine gen_avgsig(avg,sig,nstrc,wmat)
       implicit none
       integer nstrc
       real*8 avg,sig
       real*8 wmat(nstrc,nstrc)
       real*8 xnorm,tmr
       integer i,j

       avg=0.d0 ; xnorm=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j <= i) cycle
       tmr=wmat(i,j)
       avg=avg+tmr
       xnorm=xnorm+1.d0
       enddo
       enddo
       avg=avg/xnorm ; sig=0.d0
       do i=1,nstrc
       do j=1,nstrc
       if(j <= i) cycle
       tmr=wmat(i,j)
       sig=sig+(tmr-avg)**2
       enddo
       enddo
       sig=sig/xnorm ; sig=sqrt(sig)
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
       if(ncext(1) < 1) ncext(1)=1
       if(ncext(2) < 1) ncext(2)=1
       if(ncext(3) < 1) ncext(3)=1
       if(ncext(1) > 10) ncext(1)=10
       if(ncext(2) > 10) ncext(2)=10
       if(ncext(3) > 10) ncext(3)=10
       end
!234567890
!      Written by In-Ho Lee, KRISS, April 21, 2016.
!      awk '{if(NR ==1) {print "12.34", $0} else { print }}' POSCAR
!      ifort -o batch4cdiff.x strings.f90 numeral.f batch4cdiff.f90
!      ./batch4cdiff.x 
!      awk '{if(NF <3) {print}}' g1>g2
!      gnuplot> set pm3d
!      gnuplot> set palette rgbformulae 33,13,10
!      gnuplot> splot "mat.dat" with pm3d
       program batch4cdiff
       USE prdf, ONLY : prdf_init,prdf_final,get_prdf,prdf_cmp
       USE bldist, ONLY : bldist_init,bldist_final,get_blsrtd,bldist_cmp
       USE qlabmod, ONLY : qlab_init,qlab_final,get_qlab,qlab_cmp
       USE elemdist, ONLY : elemdist_init,elemdist_final,get_hist,elemdist_cmp
       USE prdf, ONLY : engg_prdf => engg
       USE bldist, ONLY : engg_bldist => engg
       USE qlabmod, ONLY : engg_qlabmod => engg
       USE elemdist, ONLY : engg_elemdist => engg
       implicit none
       integer nspecies1,ndeg1,isize,npop
       integer ktmp,i,j,nstrc0,i2220,idiff
       real*8 rmax0,tmp,omega,pi,rc10,rc20,avg,sig
       real*8, allocatable :: sigmamatrix1(:,:)
       character*2, allocatable :: symbl1(:)
       integer, allocatable :: nelements1(:)
       real*8, allocatable :: qqq1(:)
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
       read(1,*) omega
       read(1,*)
       read(1,*)
       read(1,*)
       do i=1,nspecies1
       read(1,*) (sigmamatrix1(i,j),j=1,nspecies1)
       enddo
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*) npop
       close(1)
       pi=4.d0*atan(1.d0)
       tmp=((3.d0*omega)/(4.d0*pi))**(1.d0/3.d0)
       rc10=2.3d0
       rc20=2.9d0
!      write(6,*) omega,tmp
       nstrc0=npop
       i2220=0
       lpbc0=.false.
       lpbc0=.true.
       idiff=6
!      idiff=5
!      idiff=4
!      idiff=3
!      idiff=2
!      idiff=1
!
       i2220=0
       if(idiff == 1)then
       i2220=0
                     endif
       if(idiff == 2)then
       i2220=1
                     endif
       if(idiff == 3)then
       i2220=2
                     endif
       if(idiff == 4)then
       i2220=0
                     endif
       if(idiff == 5)then
       i2220=0
                     endif
!
!      write(6,*) nstrc0
       pi=4.d0*atan(1.d0)
       tmp=((omega*3.d0)/(4.d0*pi))**(1.d0/3.d0)
!      write(6,*) omega,tmp
       ndeg1=6+3*sum(nelements1)
       allocate(qqq1(ndeg1))
       isize=4
       do i=1,nstrc0
       call xnumeral(i,fname,isize) ; fname='POSCAR_'//trim(fname) ; fname=trim(fname)
       cmd='cp '//trim(fname)//' POSCAR_a'  ; cmd=trim(cmd) ; call system(cmd)
       cmd="echo 6.00 >z1; head -n 1 POSCAR_a >>z1; awk 'ORS=NR%2?FS:RS' z1>z2; tail -n +2 POSCAR_a >>z2 ; mv z2 POSCAR_a ; rm z1"
       call system(cmd)
       filea='POSCAR_a'
       call get_poscar(filea,nspecies1,nelements1,symbl1,ndeg1,qqq1,rmax0)
!
       if(idiff == 1)then
       if(i == 1) call elemdist_init(rmax0,nspecies1,nelements1,rc10,rc20,nstrc0,lpbc0)
       open(1,file=trim(fname),form='formatted')
       read(1,*) engg_elemdist(i)
       close(1)
       call get_hist(i,qqq1)
                     endif
       if(idiff == 2)then
       if(i == 1) call bldist_init(rmax0,i2220,nelements1,nspecies1,nstrc0,lpbc0)
       open(1,file=trim(fname),form='formatted')
       read(1,*) engg_bldist(i)
       close(1)
       call get_blsrtd(i,qqq1)
                     endif
       if(idiff == 3)then
       if(i == 1) call bldist_init(rmax0,i2220,nelements1,nspecies1,nstrc0,lpbc0)
       open(1,file=trim(fname),form='formatted')
       read(1,*) engg_bldist(i)
       close(1)
       call get_blsrtd(i,qqq1)
                     endif
       if(idiff == 4)then
       if(i == 1) call qlab_init(rmax0,nspecies1,nelements1,sigmamatrix1,nstrc0,lpbc0)
       open(1,file=trim(fname),form='formatted')
       read(1,*) engg_qlabmod(i)
       close(1)
       call get_qlab(i,qqq1)
                     endif
       if(idiff == 5)then
       if(i == 1) call prdf_init(rmax0,nspecies1,nelements1,symbl1,nstrc0,lpbc0)
       open(1,file=trim(fname),form='formatted')
       read(1,*) engg_prdf(i)
       close(1)
       call get_prdf(i,qqq1)
                     endif
       if(idiff == 6)then
       ktmp=-nstrc0
       if(i == 1) call prdf_init(rmax0,nspecies1,nelements1,symbl1,ktmp,lpbc0)
       open(1,file=trim(fname),form='formatted')
       read(1,*) engg_prdf(i)
       close(1)
       call get_prdf(i,qqq1)
                     endif
!
       enddo
!
       if(idiff == 1) call elemdist_final(2,avg,sig)
       if(idiff == 2) call bldist_final(2,avg,sig)
       if(idiff == 3) call bldist_final(2,avg,sig)
       if(idiff == 4) call qlab_final(2,avg,sig)
       if(idiff == 5) call prdf_final(2,avg,sig)
       if(idiff == 6) call prdf_final(2,avg,sig)
       write(6,*) avg,sig
!
       deallocate(sigmamatrix1)
       deallocate(symbl1) ; deallocate(nelements1) ; deallocate(qqq1)
       stop
       end program batch4cdiff
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
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       real*8 function gaussianft(x,s,o)
       implicit none
       real*8 x,s,o
       real*8 pi,tmp
 
       pi=4.d0*atan(1.d0)
       tmp=-0.5d0*((x-o)/s)**2 ; if(tmp < -50.d0) tmp=-50.d0 ; if(tmp >  50.d0) tmp= 50.d0
       gaussianft=exp(tmp)/(s*sqrt(pi*2.d0))
       end
!234567890
