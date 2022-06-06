Subroutine SpTeWME(res, t, sx, sy, n, m, MAXm, gt, gs, stE1, NUM1, stE2, NUM2, covhat)

implicit none

integer :: n, m(n), MAXm, NUM1, NUM2, i1, i2, ii, jj

double precision :: gt, gs, kt, ks, tmp, num, den

double precision :: res(1:n,1:MAXm), t(n), sx(1:n,1:MAXm), sy(1:n,1:MAXm), &

                    stE1(NUM1,3), stE2(NUM2,3), wres1(NUM1), wres2(NUM2), &
                    
                    covhat(NUM1,NUM2)


do i1=1,NUM1

   wres1(i1)=0D0

   num=0D0; den=0D0

   do ii=1,n

      tmp=(t(ii)-stE1(i1,3))/gt

      if(tmp >= -1D0 .AND. tmp <= 1D0) then

         kt=0.75D0*(1D0-tmp**2D0)

         do jj=1,m(ii)

            tmp=dsqrt(dble(sx(ii,jj)-stE1(i1,1))**2D0+dble(sy(ii,jj)-stE1(i1,2))**2D0)/gs
          
            if(tmp >= -1D0 .AND. tmp <= 1D0) then
                  
               ks=0.75D0*(1D0-tmp**2D0)
                   
            else
                      
               ks=0D0
                  
            end if   

            num=num+res(ii,jj)*kt*ks
 
            den=den+kt*ks

         end do

      end if

   end do

   wres1(i1)=num/den   

end do


do i2=1,NUM2

   wres2(i2)=0D0

   num=0D0; den=0D0

   do ii=1,n

      tmp=(t(ii)-stE2(i2,3))/gt

      if(tmp >= -1D0 .AND. tmp <= 1D0) then

         kt=0.75D0*(1D0-tmp**2D0)

         do jj=1,m(ii)

            tmp=dsqrt(dble(sx(ii,jj)-stE2(i2,1))**2D0+dble(sy(ii,jj)-stE2(i2,2))**2D0)/gs
          
            if(tmp >= -1D0 .AND. tmp <= 1D0) then
                  
               ks=0.75D0*(1D0-tmp**2D0)
                   
            else
                      
               ks=0D0
                  
            end if   

            num=num+res(ii,jj)*kt*ks
 
            den=den+kt*ks

         end do

      end if

   end do

   wres2(i2)=num/den   

end do



do i1=1,NUM1

do i2=1,NUM2

   covhat(i1,i2)=0D0

   if(stE1(i1,1)==stE2(i2,1) .AND. stE1(i1,2)==stE2(i2,2) .AND. stE1(i1,3)==stE2(i2,3)) then

      num=0D0; den=0D0

      do ii=1,n

         tmp=(t(ii)-stE1(i1,3))/gt

         if(tmp >= -1D0 .AND. tmp <= 1D0) then

            kt=0.75D0*(1D0-tmp**2D0)

            do jj=1,m(ii)

               tmp=dsqrt(dble(sx(ii,jj)-stE1(i1,1))**2D0+dble(sy(ii,jj)-stE1(i1,2))**2D0)/gs
          
               if(tmp >= -1D0 .AND. tmp <= 1D0) then
                  
                  ks=0.75D0*(1D0-tmp**2D0)
                   
               else
                      
                  ks=0D0
                  
               end if   

               num=num+kt*ks*res(ii,jj)**2D0
 
               den=den+kt*ks

            end do

         end if

      end do

      covhat(i1,i2)=num/den
   
   else

      covhat(i1,i2)=wres1(i1)*wres2(i2)

   end if


end do

end do


end Subroutine SpTeWME       