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
