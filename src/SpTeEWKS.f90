Subroutine SpTeEWKS(y, t, sx, sy, n, m, MAXm, lambda, h, LOO, muhat)

implicit none

integer :: n, m(n), MAXm, i, j, ii, jj, LOO

double precision :: lambda, h, kt, ks, tmpt, tmps, tau

double precision :: y(1:n,1:MAXm), t(n), sx(1:n,1:MAXm), sy(1:n,1:MAXm), &
                    muhat(1:n,MAXm), yNUM, yDEN

tau=10D0

do i=1,n

do j=1,m(i)
    
   muhat(i,j)=0D0; yNUM=0D0; yDEN=0D0
 
   do ii=1,i

      tmpt=(t(i)-t(ii))*dble(n)

      if(tmpt <= tau) then 

         kt=(1D0-lambda)**tmpt
      
         do jj=1,m(ii)
    
            tmps=dsqrt(dble(sx(ii,jj)-sx(i,j))**2D0+dble(sy(ii,jj)-sy(i,j))**2D0)/h

            ks=max(0.75D0*(1D0-tmps**2D0),0D0)

            if(LOO .NE. 0 .AND. tmpt .EQ. 0D0 .AND. tmps .EQ. 0D0) then 

              ks=0D0

            end if

            yNUM=yNUM+y(ii,jj)*kt*ks; yDEN=yDEN+kt*ks

         end do

      end if

   end do
             
   muhat(i,j)=yNUM/yDEN

end do

end do


end Subroutine SpTeEWKS   