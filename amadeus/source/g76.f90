!234567890     
       subroutine get_intensity(fname13,xrdobj,lfault)
       USE strings, ONLY : parse,value
       implicit none
       character*280 fname13
       real*8 xrdobj
       logical lfault
       character*280 gname14
       real*8 yntensity(0:180),refyntensity(0:180)
       integer j,i,ih,ik,il
       real*8 xintsty,tthe,dhkl
       character*280 astring,bstring,cstring
       integer ios,nargs
       character*280 str1
       character*280 args(40)
       character*20 delims
       logical lfault0
       logical lexist

       lfault0=.false.
       inquire(file=trim(fname13),exist=lexist)
       if(.not. lexist)then
       lfault=.true. ; xrdobj=1.d97
                       return
                       endif
!      astring='/home/ihlee/abcd/0001/xrd.txt'
       astring=fname13
       bstring='xrd.txt'
       i=len_trim(astring)-len_trim(bstring)-5
       cstring=astring(1:i)//'refxrd.dat'
!
       gname14=cstring
       inquire(file=trim(gname14),exist=lexist)
       if(.not. lexist)then
       lfault=.true. ; xrdobj=1.d97
                       return
                       endif
       do i=0,180
       yntensity(i)=0.d0 ; refyntensity(i)=0.d0
       enddo
       open(14,file=trim(gname14),form='formatted')
       do
       read(14,'(a280)',err=811,end=899) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(args(1)(1:1) == '#') cycle
       if(args(1)(1:1) == '!') cycle
       if(nargs == 7)then
       if(args(1) == '2theta') cycle
       if(args(1) == '---') cycle
       call value(args(1),i,ios)
       call value(args(2),tthe,ios)
       call value(args(3),dhkl,ios)
       call value(args(4),ih,ios)
       call value(args(5),ik,ios)
       call value(args(6),il,ios)
       call value(args(7),xintsty,ios)
       j=tthe
       refyntensity(j)=refyntensity(j)+xintsty
!      write(6,*) tthe,xintsty
                     endif
       if(nargs == 2)then
       call value(args(1),tthe,ios)
       call value(args(2),xintsty,ios)
       j=tthe
       refyntensity(j)=refyntensity(j)+xintsty
!      write(6,*) tthe,xintsty
                     endif
       enddo
  811  continue
       lfault0=.true.
  899  continue
       close(14)
       if(lfault)then
       xrdobj=1.d97
                 return
                 endif
       open(13,file=trim(fname13),form='formatted')
       do
       read(13,'(a280)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(args(1)(1:1) == '#') cycle
       if(args(1)(1:1) == '!') cycle
       if(nargs == 7)then
       if(args(1) == '2theta') cycle
       if(args(1) == '---') cycle
       call value(args(1),i,ios)
       call value(args(2),tthe,ios)
       call value(args(3),dhkl,ios)
       call value(args(4),ih,ios)
       call value(args(5),ik,ios)
       call value(args(6),il,ios)
       call value(args(7),xintsty,ios)
       j=tthe
       yntensity(j)=yntensity(j)+xintsty
!      write(6,*) tthe,xintsty
                     endif
       enddo
  911  continue
       lfault0=.true.
  999  continue
       close(13)
       lfault=lfault0
       if(lfault)then
       xrdobj=1.d97
                 return
                 endif
       xrdobj=0.d0
       do i=0,180
       xrdobj=xrdobj+(yntensity(i)-refyntensity(i))**2
       enddo
       end 
!234567890     
