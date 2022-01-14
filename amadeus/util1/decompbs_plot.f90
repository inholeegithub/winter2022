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
! Input   : INCAR, OUTCAR, PROCAR
! Output  : sband.dat, pband.dat, dband.dat for xmgrace
!234567890
!      Written by In-Ho Lee, KRISS, November 5, 2018.
       program procar_analysis
       USE strings, ONLY : parse,value,lowercase
       implicit none
       logical lfault_stdout,lfault,lfault1
       logical lsoc
       character*280 stdname,otname
       integer ios,nargs
       character*280 str1
       character*280 args(40)
       character*20 delims
       integer i,j,k,irk,iband,iatom,ktot
       real*8 vec(3),wec(10),test,tmp,tmq,ef
       integer nk,natot,nbands
       integer ispd
       real*8, allocatable :: rk(:,:),rp(:,:),wgt(:)
       real*8, allocatable :: eig(:,:),occ(:,:),decomp(:,:,:)
       integer, allocatable :: iid(:,:)
       logical latom

       lsoc=.false.
       lfault=.false.
       open(19,file='INCAR',form='formatted')
       do
       read(19,'(a280)',err=811,end=899) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs >= 3)then
       if(lowercase(args(1)) == 'lsorbit')then
       if(lowercase(args(3)) == '.true.')then
       lsoc=.true.
!      write(6,*) args(1)
!      write(6,*) args(2)
!      write(6,*) args(3)
                                         endif
                                          endif
                     endif
       enddo
  811  continue
       lfault=.true.
  899  continue
       close(19)
       latom=.false.
       otname='OUTCAR'
       call read_outcar2(otname,ef,lfault1)
       lfault=.false.
       iatom=0
       ktot=0
       open(18,file='PROCAR',form='formatted')
       do
       read(18,'(a280)',err=911,end=999) str1
       delims=' '
       call parse(str1,delims,args,nargs)
       if(nargs == 12)then
       if(args(3) == 'k-points:' .and. args(7) == 'bands:' .and. args(11) == 'ions:')then
       call value(args(4),nk,ios)
       call value(args(8),nbands,ios)
       call value(args(12),natot,ios)
!      write(6,*) nk,nbands,natot
       allocate(rk(3,nk))
       allocate(wgt(nk))
       allocate(rp(3,natot))
       allocate(eig(nbands,nk))
       allocate(occ(nbands,nk))
       allocate(decomp(10,nbands,nk))
       allocate(iid(nbands,nk))
                                                                                     endif
                      endif
       if(nargs >=  6)then
! kpoint vector reading has problem in general due to the format problem in VASP
       if(args(1) == 'k-point')then
       if(args(3) == ':')then
       call value(args(2),irk,ios)
       call value(args(4),vec(1),ios)
       call value(args(5),vec(2),ios)
       call value(args(6),vec(3),ios)
       call value(args(9),tmp,ios)
       rk(:,irk)=vec(:)
       wgt(irk)=tmp
                         endif
                               endif
                      endif
       if(nargs ==  8)then
       if(args(1) == 'band')then
       call value(args(2),iband,ios)
       call value(args(5),tmp,ios)
       eig(iband,irk)=tmp
       call value(args(8),tmq,ios)
       occ(iband,irk)=tmq
                            endif
                      endif
       if(nargs ==  11)then
       if(args(1) == 'ion' .and. args(2) == 's' .and. args(11) == 'tot')then
       ktot=0
       iatom=-1
       latom=.true.
       cycle
                                                                        endif
                       endif
       if(nargs ==  11)then
       if(latom)then
       if(args(1) == 'tot')then
       ktot=ktot+1
       if(ktot == 4 .and. lsoc) ktot=0
       if(ktot > 1 .and. lsoc) cycle
       do j=1,10
       call value(args(j+1),wec(j),ios)
       enddo
       if(ktot == 1 .and. .not. lsoc)then
       decomp(:,iband,irk)=wec(:)
       ktot=0 ; tmp=-1.d20 ; k=0
       test=abs(decomp(1,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=1
                      endif
       test=abs(decomp(2,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=2
                      endif
       test=abs(decomp(3,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=3
                      endif
       test=abs(decomp(4,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=4
                      endif
       test=abs(decomp(5,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=5
                      endif
       test=abs(decomp(6,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=6
                      endif
       test=abs(decomp(7,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=7
                      endif
       test=abs(decomp(8,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=8
                      endif
       test=abs(decomp(9,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=9
                      endif
       iid(iband,irk)=k
                                     endif
       if(ktot == 1 .and.       lsoc)then
       decomp(:,iband,irk)=wec(:)
       tmp=-1.d20 ; k=0
       test=abs(decomp(1,iband,irk))
       if(tmp <= test)then
       tmp=test ;  k=1
                      endif
       test=abs(decomp(2,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=2
                      endif
       test=abs(decomp(3,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=3
                      endif
       test=abs(decomp(4,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=4
                      endif
       test=abs(decomp(5,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=5
                      endif
       test=abs(decomp(6,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=6
                      endif
       test=abs(decomp(7,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=7
                      endif
       test=abs(decomp(8,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=8
                      endif
       test=abs(decomp(9,iband,irk))
       if(tmp <= test)then
       tmp=test ; k=9
                      endif
       iid(iband,irk)=k
                                     endif
                           else
       call value(args(1),iatom,ios)
       do j=1,10
       call value(args(j+1),wec(j),ios)
       enddo
                           endif
                endif
                       endif
       enddo
  911  continue
       lfault=.true.
  999  continue
       close(18)
       eig=eig-ef ; ef=0.d0
       ispd=1
       call spdbsplot(nbands,nk,eig,iid,ispd)
       ispd=2
       call spdbsplot(nbands,nk,eig,iid,ispd)
       ispd=3
       call spdbsplot(nbands,nk,eig,iid,ispd)

       deallocate(eig)
       deallocate(occ)
       deallocate(rk)
       deallocate(wgt)
       deallocate(rp)
       deallocate(decomp)
       deallocate(iid)
       stop
       end program procar_analysis
!234567890
!      Written by In-Ho Lee, KRISS, November 11, 2018.
       subroutine spdbsplot(nbands,nk,eig,iid,ispd)
       implicit none
       integer nbands,nk,ispd,iid(nbands,nk)
       real*8 eig(nbands,nk)
       integer i,irk
       real*8 xxr,y1,y2,dk

       dk=0.01d0
       if(ispd == 1)then
       open(7,file='sband.dat',form='formatted')
       do i=1,nbands
       y1=-1.d222 ; y2= 1.d222
       xxr=0.d0-dk
       do irk=1,nk
       xxr=xxr+dk
       if(iid(i,irk) == 1) &
       write(7,'(2f18.8)') xxr,eig(i,irk)
       if(y1 < eig(i,irk)) y1=eig(i,irk)
       if(y2 > eig(i,irk)) y2=eig(i,irk)
       enddo
       enddo
       close(7)
                    endif
       if(ispd == 2)then
       open(7,file='pband.dat',form='formatted')
       do i=1,nbands
       y1=-1.d222 ; y2= 1.d222
       xxr=0.d0-dk
       do irk=1,nk
       xxr=xxr+dk
       if(iid(i,irk) == 2 .or. iid(i,irk) == 3 .or. iid(i,irk) == 4) &
!      if(iid(i,irk) == 2) &
       write(7,'(2f18.8)') xxr,eig(i,irk)
       if(y1 < eig(i,irk)) y1=eig(i,irk)
       if(y2 > eig(i,irk)) y2=eig(i,irk)
       enddo
       enddo
       close(7)
                    endif
       if(ispd == 3)then
       open(7,file='dband.dat',form='formatted')
       do i=1,nbands
       y1=-1.d222 ; y2= 1.d222
       xxr=0.d0-dk
       do irk=1,nk
       xxr=xxr+dk
       if(iid(i,irk) == 5 .or. iid(i,irk) == 6 .or. iid(i,irk) == 7 .or.  iid(i,irk) == 8 .or. iid(i,irk) ==9) &
!      if(iid(i,irk) == 3) &
       write(7,'(2f18.8)') xxr,eig(i,irk)
       if(y1 < eig(i,irk)) y1=eig(i,irk)
       if(y2 > eig(i,irk)) y2=eig(i,irk)
       enddo
       enddo
       close(7)
                    endif
       end
!234567890
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
