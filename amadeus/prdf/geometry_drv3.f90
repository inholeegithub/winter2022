module precision

! Real kinds

integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

! Integer kinds

integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

!Complex kinds

integer, parameter :: kc4 = kr4                            ! single precision complex
integer, parameter :: kc8 = kr8                            ! double precision complex

end module precision

module strings

use precision

private :: value_dr,value_sr,value_di,value_si
private :: write_dr,write_sr,write_di,write_si
private :: writeq_dr,writeq_sr,writeq_di,writeq_si

interface value  ! Generic operator for converting a number string to a 
                 ! number. Calling syntax is 'call value(numstring,number,ios)' 
                 ! where 'numstring' is a number string and 'number' is a 
                 ! real number or an integer (single or double precision).         
   module procedure value_dr
   module procedure value_sr
   module procedure value_di
   module procedure value_si
end interface

interface writenum  ! Generic  interface for writing a number to a string. The 
                    ! number is left justified in the string. The calling syntax
                    ! is 'call writenum(number,string,format)' where 'number' is
                    ! a real number or an integer, 'string' is a character string
                    ! containing the result, and 'format' is the format desired, 
                    ! e.g., 'e15.6' or 'i5'.
   module procedure write_dr
   module procedure write_sr
   module procedure write_di
   module procedure write_si
end interface

interface writeq  ! Generic interface equating a name to a numerical value. The
                  ! calling syntax is 'call writeq(unit,name,value,format)' where
                  ! unit is the integer output unit number, 'name' is the variable
                  ! name, 'value' is the real or integer value of the variable, 
                  ! and 'format' is the format of the value. The result written to
                  ! the output unit has the form <name> = <value>.
   module procedure writeq_dr
   module procedure writeq_sr
   module procedure writeq_di
   module procedure writeq_si
end interface


!**********************************************************************

contains

!**********************************************************************

subroutine parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

character(len=*) :: str,delims
character(len=len_trim(str)) :: strsav
character(len=*),dimension(:) :: args

strsav=str
call compact(str)
na=size(args)
do i=1,na
  args(i)=' '
end do  
nargs=0
lenstr=len_trim(str)
if(lenstr==0) return
k=0

do
   if(len_trim(str) == 0) exit
   nargs=nargs+1
   call split(str,delims,args(nargs))
   call removebksl(args(nargs))
end do   
str=strsav

end subroutine parse

!**********************************************************************

subroutine compact(str)

! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str)):: outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
isp=0
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  
  select case(ich)
  
    case(9,32)     ! space or tab character
      if(isp==0) then
        k=k+1
        outstr(k:k)=' '
      end if
      isp=1
      
    case(33:)      ! not a space, quote, or control character
      k=k+1
      outstr(k:k)=ch
      isp=0
      
  end select
  
end do

str=adjustl(outstr)

end subroutine compact

!**********************************************************************

subroutine removesp(str)

! Removes spaces, tabs, and control characters in string str

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)    
    case(0:32)  ! space, tab, or control character
         cycle       
    case(33:)  
      k=k+1
      outstr(k:k)=ch
  end select
end do

str=adjustl(outstr)

end subroutine removesp

!**********************************************************************

subroutine value_dr(str,rnum,ios)

! Converts number string to a double precision real number

character(len=*)::str
real(kr8)::rnum
integer :: ios

ilen=len_trim(str)
ipos=scan(str,'Ee')
if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
   ios=3
   return
end if
read(str,*,iostat=ios) rnum

end subroutine value_dr

!**********************************************************************

subroutine value_sr(str,rnum,ios)

! Converts number string to a single precision real number

character(len=*)::str
real(kr4) :: rnum
real(kr8) :: rnumd 

call value_dr(str,rnumd,ios)
if( abs(rnumd) > huge(rnum) ) then
  ios=15
  return
end if
if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
rnum=rnumd

end subroutine value_sr

!**********************************************************************

subroutine value_di(str,inum,ios)

! Converts number string to a double precision integer value

character(len=*)::str
integer(ki8) :: inum
real(kr8) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki8)

end subroutine value_di

!**********************************************************************

subroutine value_si(str,inum,ios)

! Converts number string to a single precision integer value

character(len=*)::str
integer(ki4) :: inum
real(kr8) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki4)

end subroutine value_si

!**********************************************************************

subroutine shiftstr(str,n)
 
! Shifts characters in in the string 'str' n positions (positive values
! denote a right shift and negative values denote a left shift). Characters
! that are shifted off the end are lost. Positions opened up by the shift 
! are replaced by spaces.

character(len=*):: str

lenstr=len(str)
nabs=iabs(n)
if(nabs>=lenstr) then
  str=repeat(' ',lenstr)
  return
end if
if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
return

end subroutine shiftstr

!**********************************************************************

subroutine insertstr(str,strins,loc)

! Inserts the string 'strins' into the string 'str' at position 'loc'. 
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. Trailing spaces of 'strins' are 
! removed prior to insertion

character(len=*):: str,strins
character(len=len(str))::tempstr

lenstrins=len_trim(strins)
tempstr=str(loc:)
call shiftstr(tempstr,lenstrins)
tempstr(1:lenstrins)=strins(1:lenstrins)
str(loc:)=tempstr
return

end subroutine insertstr

!**********************************************************************

subroutine delsubstr(str,substr)

! Deletes first occurrence of substring 'substr' from string 'str' and
! shifts characters left to fill hole. Trailing spaces or blanks are
! not considered part of 'substr'.

character(len=*):: str,substr

lensubstr=len_trim(substr)
ipos=index(str,substr)
if(ipos==0) return
if(ipos == 1) then
   str=str(lensubstr+1:)
else
   str=str(:ipos-1)//str(ipos+lensubstr:)
end if   
return

end subroutine delsubstr

!**********************************************************************

subroutine delall(str,substr)

! Deletes all occurrences of substring 'substr' from string 'str' and
! shifts characters left to fill holes.

character(len=*):: str,substr

lensubstr=len_trim(substr)
do
   ipos=index(str,substr)
   if(ipos == 0) exit
   if(ipos == 1) then
      str=str(lensubstr+1:)
   else
      str=str(:ipos-1)//str(ipos+lensubstr:)
   end if
end do   
return

end subroutine delall

!**********************************************************************

function uppercase(str) result(ucstr)

! convert string to upper case

character (len=*):: str
character (len=len_trim(str)):: ucstr

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')     
iquote=0
ucstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('a') .and. iav <= iachar('z')) then
    ucstr(i:i)=achar(iav+ioffset)
  else
    ucstr(i:i)=str(i:i)
  end if
end do
return

end function uppercase

!**********************************************************************

function lowercase(str) result(lcstr)

! convert string to lower case

character (len=*):: str
character (len=len_trim(str)):: lcstr

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
lcstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('A') .and. iav <= iachar('Z')) then
    lcstr(i:i)=achar(iav-ioffset)
  else
    lcstr(i:i)=str(i:i)
  end if
end do
return

end function lowercase

!**********************************************************************

subroutine readline(nunitr,line,ios)

! Reads line from unit=nunitr, ignoring blank lines
! and deleting comments beginning with an exclamation point(!)

character (len=*):: line

do  
  read(nunitr,'(a)', iostat=ios) line      ! read input line
  if(ios /= 0) return
  line=adjustl(line)
  ipos=index(line,'!')
  if(ipos == 1) cycle
  if(ipos /= 0) line=line(:ipos-1)
  if(len_trim(line) /= 0) exit
end do
return

end subroutine readline

!**********************************************************************

subroutine match(str,ipos,imatch)

! Sets imatch to the position in string of the delimiter matching the delimiter
! in position ipos. Allowable delimiters are (), [], {}, <>.

character(len=*) :: str
character :: delim1,delim2,ch

lenstr=len_trim(str)
delim1=str(ipos:ipos)
select case(delim1)
   case('(')
      idelim2=iachar(delim1)+1
      istart=ipos+1
      iend=lenstr
      inc=1
   case(')')
      idelim2=iachar(delim1)-1
      istart=ipos-1
      iend=1
      inc=-1
   case('[','{','<')
      idelim2=iachar(delim1)+2
      istart=ipos+1
      iend=lenstr
      inc=1
   case(']','}','>')
      idelim2=iachar(delim1)-2
      istart=ipos-1
      iend=1
      inc=-1
   case default
      write(*,*) delim1,' is not a valid delimiter'
      return
end select
if(istart < 1 .or. istart > lenstr) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if
delim2=achar(idelim2) ! matching delimiter

isum=1
do i=istart,iend,inc
   ch=str(i:i)
   if(ch /= delim1 .and. ch /= delim2) cycle
   if(ch == delim1) isum=isum+1
   if(ch == delim2) isum=isum-1
   if(isum == 0) exit
end do
if(isum /= 0) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if   
imatch=i

return

end subroutine match

!**********************************************************************

subroutine write_dr(rnum,str,fmt)

! Writes double precision real number rnum to string str using format fmt

real(kr8) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_dr

!***********************************************************************

subroutine write_sr(rnum,str,fmt)

! Writes single precision real number rnum to string str using format fmt

real(kr4) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_sr

!***********************************************************************

subroutine write_di(inum,str,fmt)

! Writes double precision integer inum to string str using format fmt

integer(ki8) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_di

!***********************************************************************

subroutine write_si(inum,str,fmt)

! Writes single precision integer inum to string str using format fmt

integer(ki4) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_si

!***********************************************************************

subroutine trimzero(str)

! Deletes nonsignificant trailing zeroes from number string str. If number
! string ends in a decimal point, one trailing zero is added.

character(len=*) :: str
character :: ch
character(len=10) :: exp

ipos=scan(str,'eE')
if(ipos>0) then
   exp=str(ipos:)
   str=str(1:ipos-1)
endif
lstr=len_trim(str)
do i=lstr,1,-1
   ch=str(i:i)
   if(ch=='0') cycle          
   if(ch=='.') then
      str=str(1:i)//'0'
      if(ipos>0) str=trim(str)//trim(exp)
      exit
   endif
   str=str(1:i)
   exit
end do
if(ipos>0) str=trim(str)//trim(exp)

end subroutine trimzero

!**********************************************************************

subroutine writeq_dr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr8) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_dr

!**********************************************************************

subroutine writeq_sr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr4) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_sr

!**********************************************************************

subroutine writeq_di(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki8) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_di

!**********************************************************************

subroutine writeq_si(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki4) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_si

!**********************************************************************

function is_letter(ch) result(res)

! Returns .true. if ch is a letter and .false. otherwise

character :: ch
logical :: res

select case(ch)
case('A':'Z','a':'z')
  res=.true.
case default
  res=.false.
end select
return

end function is_letter

!**********************************************************************

function is_digit(ch) result(res)

! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

character :: ch
logical :: res

select case(ch)
case('0':'9')
  res=.true.
case default
  res=.false.
end select
return

end function is_digit

!**********************************************************************

subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

character(len=*) :: str,delims,before
character,optional :: sep
logical :: pres
character :: ch,cha

pres=present(sep)
str=adjustl(str)
call compact(str)
lenstr=len_trim(str)
if(lenstr == 0) return        ! string str is empty
k=0
ibsl=0                        ! backslash initially inactive
before=' '
do i=1,lenstr
   ch=str(i:i)
   if(ibsl == 1) then          ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
   end if
   if(ch == '\') then          ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
   end if
   ipos=index(delims,ch)         
   if(ipos == 0) then          ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
   end if
   if(ch /= ' ') then          ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
   cha=str(i+1:i+1)            ! character is a space delimiter
   iposa=index(delims,cha)
   if(iposa > 0) then          ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
   else
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
end do
if(i >= lenstr) str=''
str=adjustl(str)              ! remove initial spaces
return

end subroutine split

!**********************************************************************

subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0
ibsl=0                        ! backslash initially inactive

do i=1,lenstr
  ch=str(i:i)
  if(ibsl == 1) then          ! backslash active
   k=k+1
   outstr(k:k)=ch
   ibsl=0
   cycle
  end if
  if(ch == '\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl

!**********************************************************************

end module strings  



!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program geometry_test
       USE strings, ONLY : parse,value
       implicit none
       integer nspecies,natom
       integer nex1,nex2,nex3
       real*8 a1(3),a2(3),a3(3),scale0,d1,d2,d3,tmp,amatrix(3,3),s6(6),pi,cellvol0,vtest
       integer, allocatable :: nelements(:),ityp(:)
       real*8, allocatable :: bdcut(:,:)
       real*8, allocatable :: qqq(:)
       real*8, allocatable :: ppp(:)
       character*2, allocatable :: symbl(:)
       character*1 c1
       integer na,i,j
       logical lpbc
       character*200 str1
       integer ios,nargs
       character*200 args(40)
       character*20 delims
!
       read(5,*) 
       read(5,*) scale0
       read(5,*) a1(1),a1(2),a1(3)
       read(5,*) a2(1),a2(2),a2(3)
       read(5,*) a3(1),a3(2),a3(3)
       if(scale0 < 0.d0)then
       amatrix(1,:)=a1(:) ; amatrix(2,:)=a2(:) ; amatrix(3,:)=a3(:)
       vtest=(amatrix(1,2)*amatrix(2,3)-amatrix(1,3)*amatrix(2,2))*amatrix(3,1) &
            +(amatrix(1,3)*amatrix(2,1)-amatrix(1,1)*amatrix(2,3))*amatrix(3,2) &
            +(amatrix(1,1)*amatrix(2,2)-amatrix(1,2)*amatrix(2,1))*amatrix(3,3)
       vtest=abs(vtest)
       vtest=abs(scale0)/vtest ; vtest=vtest**(1.d0/3.d0)
       amatrix=amatrix*vtest
       a1(:)=amatrix(1,:) ; a2(:)=amatrix(2,:) ; a3(:)=amatrix(3,:)
                        endif
       if(scale0 > 0.d0)then
       a1=a1*scale0 ; a2=a2*scale0 ; a3=a3*scale0
                        endif
       amatrix(1,:)=a1(:) ; amatrix(2,:)=a2(:) ; amatrix(3,:)=a3(:)
       cellvol0=         (amatrix(1,2)*amatrix(2,3)-amatrix(1,3)*amatrix(2,2))*amatrix(3,1)
       cellvol0=cellvol0+(amatrix(1,3)*amatrix(2,1)-amatrix(1,1)*amatrix(2,3))*amatrix(3,2)
       cellvol0=cellvol0+(amatrix(1,1)*amatrix(2,2)-amatrix(1,2)*amatrix(2,1))*amatrix(3,3)
       cellvol0=abs(cellvol0)
       call latmat(s6,amatrix,0)
       pi=4.d0*atan(1.d0)
       write(6,'(6f18.8)') (s6(i),i=1,3),(s6(i)*180.d0/pi,i=4,6)
       write(6,'( f18.8)') cellvol0
!
!      read(5,*) (symbl(j),j=1,nspecies)
!      read(5,*) (nelements(j),j=1,nspecies)
       read(5,'(a200)') str1
       delims=' '
       call parse(str1,delims,args,nargs)
       nspecies=nargs
       allocate(symbl(nspecies))
       allocate(nelements(nspecies))
       allocate(bdcut(nspecies,nspecies))
       bdcut(:,:)=2.75d0
       bdcut(1,1)=2.75d0
       do j=1,nargs
       symbl(j)=trim(args(j))
       enddo
!      print*, nargs, (symbl(j),j=1,nargs)
       read(5,'(a200)') str1
       delims=' '
       call parse(str1,delims,args,nargs)
       do j=1,nargs
       call value(args(j),i,ios)
       nelements(j)=i
       enddo
!      print*, nargs, (nelements(j),j=1,nargs)
!
       natom=sum(nelements)
       allocate(ityp(natom))
       allocate(qqq(3*natom)) ; allocate(ppp(3*natom))
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       ityp(na)=i
       enddo
       enddo
       if(natom /= na) stop
       read(5,*) c1
       if(c1 == 'S' .or. c1 == 's') read(5,*) c1
       if(c1 =='d' .or. c1 =='D')then
       do j=1,natom
       read(5,*) d1,d2,d3
       tmp=0.d0
!      call random_number(tmp) ; tmp=(tmp-0.5)*1.d-2
       qqq(3*(j-1)+1)=d1+tmp
!      call random_number(tmp) ; tmp=(tmp-0.5)*1.d-2
       qqq(3*(j-1)+2)=d2+tmp
!      call random_number(tmp) ; tmp=(tmp-0.5)*1.d-2
       qqq(3*(j-1)+3)=d3+tmp
       enddo
!      write(6,*) natom,'natom'
       ppp=qqq
       call tocarx(natom,ppp,a1,a2,a3)
       write(6,*) 'Cartesian'
       do j=1,natom
       write(6,*) ppp(3*(j-1)+1),ppp(3*(j-1)+2),ppp(3*(j-1)+3)
       enddo
       open(91,file='fort91.xyz',form='formatted')
       write(91,*) natom
       write(91,*) 
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       write(91,'(a2,2x,3f18.8)') symbl(i), ppp(3*(na-1)+1),ppp(3*(na-1)+2),ppp(3*(na-1)+3)
       enddo
       enddo
       close(91)
                                 endif
       if(c1 =='c' .or. c1 =='C')then
       do j=1,natom
       read(5,*) ppp(3*(j-1)+1),ppp(3*(j-1)+2),ppp(3*(j-1)+3)
       enddo
       open(91,file='fort91.xyz',form='formatted')
       write(91,*) natom
       write(91,*) 
       na=0
       do i=1,nspecies
       do j=1,nelements(i)
       na=na+1
       write(91,'(a2,2x,3f18.8)') symbl(i), ppp(3*(na-1)+1),ppp(3*(na-1)+2),ppp(3*(na-1)+3)
       enddo
       enddo
       close(91)
       qqq=ppp
       call tolatx(natom,qqq,a1,a2,a3)
       write(6,*) 'Direct'
       do j=1,natom
       write(6,*) qqq(3*(j-1)+1),qqq(3*(j-1)+2),qqq(3*(j-1)+3)
       enddo
                                 endif
!      call bdangles(nspecies,nelements,ityp,symbl,natom,bdcut,qqq,a1,a2,a3)
       nex1=1 ; nex2=1 ; nex3=1
       nex1=2 ; nex2=2 ; nex3=2
       call ext_cell(nex1,nex2,nex3,a1,a2,a3,natom,qqq,symbl,nelements,nspecies,ityp,bdcut)

       deallocate(bdcut)
       deallocate(qqq) ; deallocate(ppp)
       deallocate(symbl) ; deallocate(nelements) ; deallocate(ityp)
       end program geometry_test
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gen_coordination(natom,rc1,ncoord,a1,a2,a3,qqq)
       implicit none
       integer natom,ncoord(natom)
       real*8 rc1,a1(3),a2(3),a3(3),qqq(3*natom)
       integer i,j,kk1
       real*8 x,y,z,r,d1,d2,d3

       do i=1,natom
       kk1=0
       do j=1,natom
       if(j == i) cycle
       d1=qqq(3*(i-1)+1)-qqq(3*(j-1)+1)
       d2=qqq(3*(i-1)+2)-qqq(3*(j-1)+2)
       d3=qqq(3*(i-1)+3)-qqq(3*(j-1)+3)
       d1=d1-nint(d1)
       d2=d2-nint(d2)
       d3=d3-nint(d3)
       x=d1*a1(1)+d2*a2(1)+d3*a3(1)
       y=d1*a1(2)+d2*a2(2)+d3*a3(2)
       z=d1*a1(3)+d2*a2(3)+d3*a3(3)
       r=sqrt(x*x+y*y+z*z)
       if( r <= rc1)then
       kk1=kk1+1
                    endif
       enddo
       ncoord(i)=kk1
       enddo
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine bdangles(nspecies,nelements,ityp,symbl,natom,bdcut,qqq,a1,a2,a3)
       implicit none
       integer nspecies,natom,nelements(nspecies),ityp(natom)
       real*8 bdcut(nspecies,nspecies),qqq(3*natom),a1(3),a2(3),a3(3)
       character*2 symbl(nspecies)
       integer ih,ia,ib,ic,id
       real*8 xnorm1,xnorm2,xnorm3,d1,d2,d3,bdlng,ang,diang,bab,bbc,bcd,rc1
       real*8 ynorm1,ynorm2,ynorm3
       real*8 bdl_hist(0:1001),bda_hist(-10001:10001),bdd_hist(-10001:10001)
       real*8 avg1,avg2,avg3,sig1,sig2,sig3
       real*8 zipr,tmp,tmq
       integer, allocatable :: ncoord(:)
       real*8, external :: geometry

       bdl_hist=0.d0
       bda_hist=0.d0
       bdd_hist=0.d0
!
       open(71,file='fort.71',form='formatted')
       write(6,*) 'bond length analysis'
       ynorm1=0.d0
       avg1=0.d0
       sig1=0.d0
       xnorm1=0.d0
       do ia=1,natom-1
       do ib=ia+1,natom
       bdlng=geometry(a1,a2,a3,natom,qqq,ia,ib,0,0) 
       ih=bdlng/0.10d0 ; bdl_hist(ih)=bdl_hist(ih)+1.d0
       xnorm1=xnorm1+1.d0
       if(bdlng > bdcut(ityp(ia),ityp(ib))) cycle
       avg1=avg1+bdlng
       sig1=sig1+bdlng**2
       ynorm1=ynorm1+1.d0
       write(6,'(f16.5,2x,2i6)') bdlng,ia,ib
       enddo
       enddo
       if(xnorm1 >0.0d0) bdl_hist=bdl_hist/xnorm1
       do ih=0,1000
       if(float(ih)*0.1d0 < 15.d0) write(71,*) float(ih)*0.1d0,bdl_hist(ih)
       enddo
       close(71)
! 
       open(72,file='fort.72',form='formatted')
       write(6,*) 'bond angle analysis'
       ynorm2=0.d0
       avg2=0.d0
       sig2=0.d0
       xnorm2=0.d0
       do ic=1,natom
       do ia=1,natom-1
       if(ic == ia) cycle
       do ib=ia+1,natom
       if(ic == ib) cycle
       bdlng=geometry(a1,a2,a3,natom,qqq,ib,ic,0,0) ; if(bdlng > bdcut(ityp(ib),ityp(ic))) cycle
       bdlng=geometry(a1,a2,a3,natom,qqq,ia,ic,0,0) ; if(bdlng > bdcut(ityp(ia),ityp(ic))) cycle
       ang=geometry(a1,a2,a3,natom,qqq,ia,ic,ib,0)
       ih=ang/1.d0 ; bda_hist(ih)=bda_hist(ih)+1.d0
       xnorm2=xnorm2+1.d0
       avg2=avg2+ang
       sig2=sig2+ang**2
       ynorm2=ynorm2+1.d0
       write(6,'(f16.5,2x,3i6, 2f16.5)') ang,ia,ic,ib, &
          geometry(a1,a2,a3,natom,qqq,ib,ic,0,0), geometry(a1,a2,a3,natom,qqq,ia,ic,0,0)
       enddo
       enddo
       enddo
       if(xnorm2 >0.0d0) bda_hist=bda_hist/xnorm2
       do ih=-1000,1000
       if(float(ih)*1.d0 > -182.d0 .and.  float(ih)*1.d0 < 182.d0) write(72,*) float(ih)*1.d0,bda_hist(ih)
       enddo
       close(72)
!
       open(73,file='fort.73',form='formatted')
       write(6,*) 'dihedral angle analysis'
       ynorm3=0.d0
       avg3=0.d0
       sig3=0.d0
       xnorm3=0.d0
       do ib=1,natom-1
       do ic=ib+1,natom
       bdlng=geometry(a1,a2,a3,natom,qqq,ib,ic,0,0) ; if(bdlng > bdcut(ityp(ib),ityp(ic))) cycle
         bbc=geometry(a1,a2,a3,natom,qqq,ib,ic,0,0)

       do ia=1,natom
       if(ia == ib) cycle
       if(ia == ic) cycle
       bdlng=geometry(a1,a2,a3,natom,qqq,ia,ib,0,0) ; if(bdlng > bdcut(ityp(ia),ityp(ib))) cycle
         bab=geometry(a1,a2,a3,natom,qqq,ia,ib,0,0)
       do id=1,natom
       if(id == ia) cycle
       if(id == ib) cycle
       if(id == ic) cycle
       bdlng=geometry(a1,a2,a3,natom,qqq,ic,id,0,0) ; if(bdlng > bdcut(ityp(ic),ityp(id))) cycle
         bcd=geometry(a1,a2,a3,natom,qqq,ic,id,0,0)
       diang=geometry(a1,a2,a3,natom,qqq,ia,ib,ic,id)
       write(6,'(f16.5,2x,4i6,2x,3f16.5)') diang,ia,ib,ic,id, bab,bbc,bcd
       ih=diang/1.d0 ; bdd_hist(ih)=bdd_hist(ih)+1.d0
       xnorm3=xnorm3+1.d0
       avg3=avg3+diang
       sig3=sig3+diang**2
       ynorm3=ynorm3+1.d0
       enddo
       enddo
       enddo
       enddo
       if(xnorm3 >0.0d0) bdd_hist=bdd_hist/xnorm3
       do ih=-1000,1000
       if(float(ih)*1.d0 > -182.d0 .and.  float(ih)*1.d0 < 182.d0) write(73,*) float(ih)*1.d0,bdd_hist(ih)
       enddo
       close(73)
       avg1=avg1/ynorm1
       avg2=avg2/ynorm2
       avg3=avg3/ynorm3
       sig1=sig1/ynorm1
       sig2=sig2/ynorm2
       sig3=sig3/ynorm3
       sig1=sqrt(abs(sig1-avg1**2))
       sig2=sqrt(abs(sig2-avg2**2))
       sig3=sqrt(abs(sig3-avg3**2))
       write(6,'(2f18.6,2x,a7)') avg1,sig1,'avg,sig'
       write(6,'(2f18.6,2x,a7)') avg2,sig2,'avg,sig'
       write(6,'(2f18.6,2x,a7)') avg3,sig3,'avg,sig'
!
       tmp=0.d0
       do ih=-180,180
       tmp=tmp+bda_hist(ih)**2
       enddo
       write(6,*) 'IPR angle'
       zipr=0.d0
       do ih=-180,180
       zipr=zipr+bda_hist(ih)**4
       enddo
       write(6,*) zipr/tmp


       tmq=0.d0
       do ih=-180,180
       tmq=tmq+bdd_hist(ih)**2
       enddo
       write(6,*) 'IPR dihedral'
       zipr=0.d0
       do ih=-180,180
       zipr=zipr+bdd_hist(ih)**4
       enddo
       write(6,*) zipr/tmq
!
       rc1=bdcut(1,1)
       allocate(ncoord(natom))
       call gen_coordination(natom,rc1,ncoord,a1,a2,a3,qqq)
       write(6,*) 'coordination number analysis with rc1',rc1
       do ia=1,natom
       write(6,*) ia, ncoord(ia)
       enddo
       deallocate(ncoord)
!
       end
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine ext_cell(nex1,nex2,nex3,a1,a2,a3,natom,qqq,symbl1,nelements1,nspecies,ityp1,bdcut)
       implicit none
       integer nex1,nex2,nex3
       integer natom,nspecies,nelements1(nspecies),ityp1(natom)
       character*2 symbl1(nspecies)
       real*8 a1(3),a2(3),a3(3),qqq(3*natom),bdcut(nspecies,nspecies)
       real*8, allocatable :: qext(:)
       integer, allocatable :: nelements(:),ityp(:)
       character*2, allocatable :: symbl(:)
       integer nb,i,j,k,m
       real*8 aa1(3),aa2(3),aa3(3),x,y,z,r,d1,d2,d3

       aa1=a1*nex1
       aa2=a2*nex2
       aa3=a3*nex3
       nb=natom* (nex1*nex2*nex3)
       allocate(qext(3*nb))
       allocate(symbl(nspecies))
       do i=1,nspecies
       symbl(i)=symbl1(i)
       enddo
       allocate(nelements(nspecies))
       nelements=nelements1* (nex1*nex2*nex3)
       allocate(ityp(nb))
       nb=0
       do m=1,natom
       do i=0,nex1-1
       do j=0,nex2-1
       do k=0,nex3-1
       nb=nb+1
       ityp(nb)=ityp1(m)
       qext(3*(nb-1)+1)=qqq(3*(m-1)+1)+i
       qext(3*(nb-1)+2)=qqq(3*(m-1)+2)+j
       qext(3*(nb-1)+3)=qqq(3*(m-1)+3)+k
       enddo
       enddo
       enddo
       enddo
       do i=1,nb
       qext(3*(i-1)+1)=qext(3*(i-1)+1)/float(nex1)
       qext(3*(i-1)+2)=qext(3*(i-1)+2)/float(nex2)
       qext(3*(i-1)+3)=qext(3*(i-1)+3)/float(nex3)
       enddo

       call bdangles(nspecies,nelements,ityp,symbl,nb,bdcut,qext,aa1,aa2,aa3)

!      call ggsofk(nb,qext,aa1,aa2,aa3)
       call gofr(nb,qext,aa1,aa2,aa3)

       deallocate(symbl,nelements,ityp)
       deallocate(qext)
       return
       end
!
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gofr(natom,qqq,a1,a2,a3)
       implicit none
       integer natom
       real*8 qqq(3*natom),a1(3),a2(3),a3(3)
       real*8 pi,dr,tmp,cmatrix(3,3),cellvol0
       real*8 hist(10000),gr(10000),const
       real*8 rl,ru,xx,rr,xk
       real*8, external :: geometry
       integer i,j,ih,npt

       pi=4.d0*atan(1.d0)
       cmatrix(1,:)=a1(:) ; cmatrix(2,:)=a2(:) ; cmatrix(3,:)=a3(:)
       cellvol0=         (cmatrix(1,2)*cmatrix(2,3)-cmatrix(1,3)*cmatrix(2,2))*cmatrix(3,1)
       cellvol0=cellvol0+(cmatrix(1,3)*cmatrix(2,1)-cmatrix(1,1)*cmatrix(2,3))*cmatrix(3,2)
       cellvol0=cellvol0+(cmatrix(1,1)*cmatrix(2,2)-cmatrix(1,2)*cmatrix(2,1))*cmatrix(3,3)
       cellvol0=abs(cellvol0)
       dr=0.1d0
       gr=0.d0 ; hist=0.d0
       do i=1,natom-1
       do j=i+1,natom
       tmp=geometry(a1,a2,a3,natom,qqq,i,j,0,0)
       ih=int(tmp/dr)+1
       hist(ih)=hist(ih)+2.d0
       enddo
       enddo
       const=(4.d0*pi/3.d0)*(float(natom)/cellvol0)
       do i=1,1000
       rl=float(i-1)*dr
       ru=rl+dr
       xx=const*(ru**3-rl**3)
       gr(i)=hist(i)/float(natom)/xx
       enddo
       npt=1000
       do i=1,npt
       if(float(i-1)*dr > 20.d0) exit
       write(6,*) float(i-1)*dr,gr(i)
       enddo


       xk=0.1d0
       tmp=0.d0
       do i=1,1000
       rr=(float(i)*dr)
       tmp=tmp+rr**2 *sin(xk*rr)/(xk*rr) *gr(i)
       enddo
       tmp=tmp*dr
       tmp=tmp*4.d0*pi*(natom/cellvol0)+1.d0

       return
       end
!
!
!     ###################################################
!     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
!     ##              All Rights Reserved              ##
!     ###################################################
!
!     ################################################################
!     ##                                                            ##
!     ##  function geometry  --  evaluate distance, angle, torsion  ##
!     ##                                                            ##
!     ################################################################
!
!
!     "geometry" finds the value of the interatomic distance, angle
!     or dihedral angle defined by two to four input atoms
!
!
      function geometry(a1,a2,a3,natom,qqq,ia,ib,ic,id)
      implicit none
      integer natom
      real*8 qqq(3*natom),a1(3),a2(3),a3(3)
      real*8 d1,d2,d3
      integer ia,ib,ic,id
      real*8 xab,yab,zab
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 rab2,rcb2,rabc
      real*8 rt2,ru2,rtru
      real*8 cosine,sign
      real*8 radian 
      parameter (radian=57.29577951308232088d0)
      real*8 geometry
!
!
!     set default in case atoms are coincident or colinear
!
      geometry = 0.0d0
!
!     compute the value of the distance in angstroms
!
      if (ic .eq. 0) then
         xab = qqq(3*(ia-1)+1)-qqq(3*(ib-1)+1) 
         yab = qqq(3*(ia-1)+2)-qqq(3*(ib-1)+2)
         zab = qqq(3*(ia-1)+3)-qqq(3*(ib-1)+3)

         d1=xab-nint(xab)
         d2=yab-nint(yab)
         d3=zab-nint(zab)
         xab=d1*a1(1)+d2*a2(1)+d3*a3(1)
         yab=d1*a1(2)+d2*a2(2)+d3*a3(2)
         zab=d1*a1(3)+d2*a2(3)+d3*a3(3)


         geometry = sqrt(xab*xab + yab*yab + zab*zab)
!
!     compute the value of the angle in degrees
!
      else if (id .eq. 0) then
         xab =  qqq(3*(ia-1)+1)-qqq(3*(ib-1)+1)
         yab =  qqq(3*(ia-1)+2)-qqq(3*(ib-1)+2)
         zab =  qqq(3*(ia-1)+3)-qqq(3*(ib-1)+3)
         xcb =  qqq(3*(ic-1)+1)-qqq(3*(ib-1)+1)
         ycb =  qqq(3*(ic-1)+2)-qqq(3*(ib-1)+2)
         zcb =  qqq(3*(ic-1)+3)-qqq(3*(ib-1)+3)

         d1=xab-nint(xab)
         d2=yab-nint(yab)
         d3=zab-nint(zab)
         xab=d1*a1(1)+d2*a2(1)+d3*a3(1)
         yab=d1*a1(2)+d2*a2(2)+d3*a3(2)
         zab=d1*a1(3)+d2*a2(3)+d3*a3(3)

         d1=xcb-nint(xcb)
         d2=ycb-nint(ycb)
         d3=zcb-nint(zcb)
         xcb=d1*a1(1)+d2*a2(1)+d3*a3(1)
         ycb=d1*a1(2)+d2*a2(2)+d3*a3(2)
         zcb=d1*a1(3)+d2*a2(3)+d3*a3(3)


         rab2 = xab*xab + yab*yab + zab*zab
         rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
         rabc = sqrt(rab2 * rcb2)
         if (rabc .ne. 0.0d0) then
            cosine = (xab*xcb + yab*ycb + zab*zcb) / rabc
            cosine = min(1.0d0,max(-1.0d0,cosine))
            geometry = radian * acos(cosine)
         end if
!
!     compute the value of the dihedral angle in degrees
!
      else
         xba =   qqq(3*(ib-1)+1)-qqq(3*(ia-1)+1)
         yba =   qqq(3*(ib-1)+2)-qqq(3*(ia-1)+2)
         zba =   qqq(3*(ib-1)+3)-qqq(3*(ia-1)+3)
         xcb =   qqq(3*(ic-1)+1)-qqq(3*(ib-1)+1)
         ycb =   qqq(3*(ic-1)+2)-qqq(3*(ib-1)+2)
         zcb =   qqq(3*(ic-1)+3)-qqq(3*(ib-1)+3)
         xdc =   qqq(3*(id-1)+1)-qqq(3*(ic-1)+1)
         ydc =   qqq(3*(id-1)+2)-qqq(3*(ic-1)+2)
         zdc =   qqq(3*(id-1)+3)-qqq(3*(ic-1)+3)

         d1=xba-nint(xba)
         d2=yba-nint(yba)
         d3=zba-nint(zba)
         xba=d1*a1(1)+d2*a2(1)+d3*a3(1)
         yba=d1*a1(2)+d2*a2(2)+d3*a3(2)
         zba=d1*a1(3)+d2*a2(3)+d3*a3(3)

         d1=xcb-nint(xcb)
         d2=ycb-nint(ycb)
         d3=zcb-nint(zcb)
         xcb=d1*a1(1)+d2*a2(1)+d3*a3(1)
         ycb=d1*a1(2)+d2*a2(2)+d3*a3(2)
         zcb=d1*a1(3)+d2*a2(3)+d3*a3(3)

         d1=xdc-nint(xdc)
         d2=ydc-nint(ydc)
         d3=zdc-nint(zdc)
         xdc=d1*a1(1)+d2*a2(1)+d3*a3(1)
         ydc=d1*a1(2)+d2*a2(2)+d3*a3(2)
         zdc=d1*a1(3)+d2*a2(3)+d3*a3(3)



         xt = yba*zcb - ycb*zba
         yt = xcb*zba - xba*zcb
         zt = xba*ycb - xcb*yba
         xu = ycb*zdc - ydc*zcb
         yu = xdc*zcb - xcb*zdc
         zu = xcb*ydc - xdc*ycb
         rt2 = xt*xt + yt*yt + zt*zt
         ru2 = xu*xu + yu*yu + zu*zu
         rtru = sqrt(rt2 * ru2)
         if (rtru .ne. 0.0d0) then
            cosine = (xt*xu + yt*yu + zt*zu) / rtru
            cosine = min(1.0d0,max(-1.0d0,cosine))
            geometry = radian * acos(cosine)
            sign = xba*xu + yba*yu + zba*zu
            if (sign .lt. 0.0d0)  geometry = -geometry
         end if
      end if
      return
      end
!
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tocarx(natom,qqq,a1,a2,a3)
       implicit none
       integer natom
       real*8 qqq(3*natom),a1(3),a2(3),a3(3)
       real*8 x,y,z
       integer j

       do j=1,natom
       x=(a1(1)*qqq(3*(j-1)+1)+a2(1)*qqq(3*(j-1)+2)+a3(1)*qqq(3*(j-1)+3))
       y=(a1(2)*qqq(3*(j-1)+1)+a2(2)*qqq(3*(j-1)+2)+a3(2)*qqq(3*(j-1)+3))
       z=(a1(3)*qqq(3*(j-1)+1)+a2(3)*qqq(3*(j-1)+2)+a3(3)*qqq(3*(j-1)+3))
       qqq(3*(j-1)+1)=x
       qqq(3*(j-1)+2)=y
       qqq(3*(j-1)+3)=z
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine tolatx(natom,qqq,a1,a2,a3)
       implicit none
       integer natom
       real*8 qqq(3*natom),a1(3),a2(3),a3(3)
       real*8 d1,d2,d3,b(3,3),devid
       integer j

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
       do j=1,natom
       d1=(b(1,1)*qqq(3*(j-1)+1)+b(1,2)*qqq(3*(j-1)+2)+b(1,3)*qqq(3*(j-1)+3))
       d2=(b(2,1)*qqq(3*(j-1)+1)+b(2,2)*qqq(3*(j-1)+2)+b(2,3)*qqq(3*(j-1)+3))
       d3=(b(3,1)*qqq(3*(j-1)+1)+b(3,2)*qqq(3*(j-1)+2)+b(3,3)*qqq(3*(j-1)+3))
       qqq(3*(j-1)+1)=d1
       qqq(3*(j-1)+2)=d2
       qqq(3*(j-1)+3)=d3
       enddo
       do j=1,natom
       qqq(3*(j-1)+1)=qqq(3*(j-1)+1)-anint(qqq(3*(j-1)+1))
       qqq(3*(j-1)+2)=qqq(3*(j-1)+2)-anint(qqq(3*(j-1)+2))
       qqq(3*(j-1)+3)=qqq(3*(j-1)+3)-anint(qqq(3*(j-1)+3))
       enddo
       do j=1,natom
       if(qqq(3*(j-1)+1) <0.d0) qqq(3*(j-1)+1)=qqq(3*(j-1)+1)+1.d0
       if(qqq(3*(j-1)+2) <0.d0) qqq(3*(j-1)+2)=qqq(3*(j-1)+2)+1.d0
       if(qqq(3*(j-1)+3) <0.d0) qqq(3*(j-1)+3)=qqq(3*(j-1)+3)+1.d0
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, January 28, 2013.
       subroutine latmat(t6,am,idir)
       implicit none
       integer idir
       real*8 t6(6),am(3,3)
       real*8 ra,rb,rc,cosinea,cosineb,cosinec,anglea,angleb,anglec

       if(idir == 1)then
       am=0.d0
       am(1,1)=t6(1)
       am(2,1)=t6(2)*cos(t6(6))
       am(2,2)=t6(2)*sin(t6(6))
       am(3,1)=t6(3)*cos(t6(5))
       am(3,2)=t6(3)*cos(t6(4))*sin(t6(6))-((t6(3)*cos(t6(5))-t6(3)*cos(t6(4))*cos(t6(6)))/tan(t6(6)))
       am(3,3)=sqrt(t6(3)**2 -am(3,1)**2 -am(3,2)**2)
                    else
       t6=0.d0
       ra=sqrt(am(1,1)**2+am(1,2)**2+am(1,3)**2)
       rb=sqrt(am(2,1)**2+am(2,2)**2+am(2,3)**2)
       rc=sqrt(am(3,1)**2+am(3,2)**2+am(3,3)**2)
       cosinea=(am(2,1)*am(3,1)+am(2,2)*am(3,2)+am(2,3)*am(3,3))/rb/rc
       cosineb=(am(1,1)*am(3,1)+am(1,2)*am(3,2)+am(1,3)*am(3,3))/ra/rc
       cosinec=(am(1,1)*am(2,1)+am(1,2)*am(2,2)+am(1,3)*am(2,3))/ra/rb
       anglea=acos(cosinea)
       angleb=acos(cosineb)
       anglec=acos(cosinec)
       t6(1)=ra
       t6(2)=rb
       t6(3)=rc
       t6(4)=anglea
       t6(5)=angleb
       t6(6)=anglec
                    endif
       end 
!234567890
       subroutine ggsofk(natom,qqq,a1,a2,a3)
       implicit none
       integer natom
       real*8 qqq(3*natom),a1(3),a2(3),a3(3)
       real*8 xk,sofk
       real*8 dk,t2,t1
       integer i,npt

       npt=1000
       t2=21.d0
       t1=1.d-6
       
       dk=(t2-t1)/float(npt-1)
       do i=1,npt
       xk=t1+dk*float(i-1)
       call gsofk(xk,sofk,natom,qqq,a1,a2,a3)
       write(6,*) xk,sofk
       enddo

       return
       end
!234567890
       subroutine gsofk(xk,sofk,natom,qqq,a1,a2,a3)
       implicit none
       integer natom
       real*8 xk,sofk,qqq(3*natom),a1(3),a2(3),a3(3)
       real*8 first,second,tmq,tmp,brg,arg,r,x,y,z,d1,d2,d3
       integer i,j
       real*8, external :: aofr
        
       brg=1.d8
       first=1.d0
       do i=1,natom
       do j=1,natom
       if(j == i) cycle
       d1=qqq(3*(i-1)+1)-qqq(3*(j-1)+1)
       d2=qqq(3*(i-1)+2)-qqq(3*(j-1)+2)
       d3=qqq(3*(i-1)+3)-qqq(3*(j-1)+3)
       d1=d1-anint(d1)
       d2=d2-anint(d2)
       d3=d3-anint(d3)
       x=d1*a1(1)+d2*a2(1)+d3*a3(1)
       y=d1*a1(2)+d2*a2(2)+d3*a3(2)
       z=d1*a1(3)+d2*a2(3)+d3*a3(3)
       r=sqrt(x*x+y*y+z*z)
       arg=xk*r
       tmp=aofr(r,brg)
       if(arg > 1.d-6)then 
       tmp=aofr(r,brg)*sin(arg)/arg
                      endif
       first=first-tmp
       enddo
       enddo
       first=first/float(natom)
       second=0.d0
       do i=1,natom
       do j=1,natom
       if(j == i) cycle
       d1=qqq(3*(i-1)+1)-qqq(3*(j-1)+1)
       d2=qqq(3*(i-1)+2)-qqq(3*(j-1)+2)
       d3=qqq(3*(i-1)+3)-qqq(3*(j-1)+3)
       d1=d1-anint(d1)
       d2=d2-anint(d2)
       d3=d3-anint(d3)
       x=d1*a1(1)+d2*a2(1)+d3*a3(1)
       y=d1*a1(2)+d2*a2(2)+d3*a3(2)
       z=d1*a1(3)+d2*a2(3)+d3*a3(3)
       r=sqrt(x*x+y*y+z*z)
       arg=xk*r
       tmq=aofr(r,brg)
       if(arg > 1.d-6)then 
       tmq=aofr(r,brg)*sin(arg)/arg
                      endif
       second=second+tmq
       enddo
       enddo
       second=second/float(natom)
       sofk=first+second
       return
       end
!234567890
       real*8 function aofr(x,brg)
       implicit none
       real*8 x,brg
       real*8 tmp,tmq
 
       aofr=0.d0
       if(x < 2.d0* brg)then
       tmq=(1.d0-x/(2.d0*brg))**2  *(1.d0+x/(4.d0*brg))
       tmp=1.d0-(3.d0/4.d0)*(x/brg)+(1.d0/16.d0) *(x/brg)**3
       if(abs(tmp-tmq) >1.d-8)then
       write(6,*) 'something went wrong'
                              stop
                              endif
       aofr=tmq
                        endif
       return
       end
!234567890
