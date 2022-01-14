!
!234567890
      subroutine rmarin(ij,kl)
!  This subroutine and the next function generate random numbers. See 
!  the comments for SA for more information. The only changes from the 
!  orginal code is that (1) the test to make sure that RMARIN runs first 
!  was taken out since SA assures that this is done (this test didn't 
!  compile under IBM's VS Fortran) and (2) typing ivec as integer was 
!  taken out since ivec isn't used. With these exceptions, all following 
!  lines are original. 

! This is the initialization routine for the random number generator 
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
      real u(97), c, cd, cm
      integer i97, j97
      common /raset1/ u, c, cd, cm, i97, j97
      if( ij .lt. 0  .or.  ij .gt. 31328  .or.  &
          kl .lt. 0  .or.  kl .gt. 30081 ) then
          print '(a)', ' The first random number seed must have a value  &
    & between 0 and 31328'
          print '(a)',' The second seed must have a value between 0 and  &
    & 30081'
            stop
      endif
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
  3      continue
         u(ii) = s
  2   continue
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      return
      end
!
!234567890
      function ranmar()
      real u(97), c, cd, cm
      integer i97, j97
      common /raset1/ u, c, cd, cm, i97, j97
         uni = u(i97) - u(j97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         u(i97) = uni
         i97 = i97 - 1
         if(i97 .eq. 0) i97 = 97
         j97 = j97 - 1
         if(j97 .eq. 0) j97 = 97
         c = c - cd
         if( c .lt. 0.0 ) c = c + cm
         uni = uni - c
         if( uni .lt. 0.0 ) uni = uni + 1.0
         ranmar = uni
      return
      end
