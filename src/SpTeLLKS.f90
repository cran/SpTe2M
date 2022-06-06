Subroutine SpTeLLKS(y, t, sx, sy, n, m, MAXm, ht, hs, stE, NUM, eps, muhat)

implicit none

integer :: n, m(n), MAXm, NUM, i, ii, jj

double precision :: ht, hs, eps, kt, ks, tmp, detB

double precision :: y(1:n,1:MAXm), t(n), sx(1:n,1:MAXm), sy(1:n,1:MAXm), &
                    stE(NUM,3), muhat(NUM), B(4,4), Ky(4), adjB(4)


do i=1,NUM
    
   muhat(i)=0D0

   Ky(1)=0D0; Ky(2)=0D0;  Ky(3)=0D0;  Ky(4)=0D0
             
   B(1,1)=0D0; B(1,2)=0D0; B(1,3)=0D0; B(1,4)=0D0
                        
   B(2,1)=0D0; B(2,2)=0D0; B(2,3)=0D0; B(2,4)=0D0
                                
   B(3,1)=0D0; B(3,2)=0D0; B(3,3)=0D0; B(3,4)=0D0

   B(4,1)=0D0; B(4,2)=0D0; B(4,3)=0D0; B(4,4)=0D0
 
   do ii=1,n

      tmp=(t(ii)-stE(i,3))/ht

      if(tmp >= -1D0 .AND. tmp <= 1D0) then 

         if(eps .NE. 0D0) then 
            
            if(tmp > -eps .AND. tmp < eps) then
                  
                kt=3D0*(1D0-eps**2D0)/(4D0-3D0*eps-eps**3D0)/eps*dabs(tmp)
                
            else
                      
                kt=3D0/(4D0-3D0*eps-eps**3D0)*(1D0-tmp**2D0)
                  
            end if

         else 

            kt=0.75D0*(1D0-tmp**2D0)

         end if
      
         do jj=1,m(ii)
    
            tmp=dsqrt(dble(sx(ii,jj)-stE(i,1))**2D0+dble(sy(ii,jj)-stE(i,2))**2D0)/hs

            if(eps .NE. 0D0) then 
          
               if(tmp > -eps .AND. tmp < eps) then
                  
                  ks=3D0*(1D0-eps**2D0)/(4D0-3D0*eps-eps**3D0)/eps*dabs(tmp)
 
               else if ((tmp > -1D0 .AND. tmp < -eps) .OR. (tmp < 1D0 .AND. tmp > eps)) then
                      
                  ks=3D0/(4D0-3D0*eps-eps**3D0)*(1D0-tmp**2D0)
                   
               else
                       
                  ks=0D0
                  
               end if  
 
            else

               ks=max(0D0,0.75D0*(1D0-tmp**2D0))

            end if


            Ky(1)=Ky(1)+y(ii,jj)*kt*ks

            Ky(2)=Ky(2)+(t(ii)-stE(i,3))*y(ii,jj)*kt*ks
                       
            Ky(3)=Ky(3)+(sx(ii,jj)-stE(i,1))*y(ii,jj)*kt*ks
                              
            Ky(4)=Ky(4)+(sy(ii,jj)-stE(i,2))*y(ii,jj)*kt*ks
 
                             
            B(1,1)=B(1,1)+kt*ks

            B(1,2)=B(1,2)+(t(ii)-stE(i,3))*kt*ks
                      
            B(1,3)=B(1,3)+(sx(ii,jj)-stE(i,1))*kt*ks
                          
            B(1,4)=B(1,4)+(sy(ii,jj)-stE(i,2))*kt*ks 
                    
            B(2,2)=B(2,2)+(t(ii)-stE(i,3))**2D0*kt*ks 
                    
            B(2,3)=B(2,3)+(t(ii)-stE(i,3))*(sx(ii,jj)-stE(i,1))*kt*ks
                     
            B(2,4)=B(2,4)+(t(ii)-stE(i,3))*(sy(ii,jj)-stE(i,2))*kt*ks

            B(3,3)=B(3,3)+(sx(ii,jj)-stE(i,1))**2D0*kt*ks
                      
            B(3,4)=B(3,4)+(sx(ii,jj)-stE(i,1))*(sy(ii,jj)-stE(i,2))*kt*ks 
                     
            B(4,4)=B(4,4)+(sy(ii,jj)-stE(i,2))**2D0*kt*ks

            B(2,1)=B(1,2)
                 
            B(3,1)=B(1,3)
                          
            B(3,2)=B(2,3)
                      
            B(4,1)=B(1,4)
                      
            B(4,2)=B(2,4)
                     
            B(4,3)=B(3,4)
 
         end do

      end if

   end do

   adjB(1)=B(2,2)*B(3,3)*B(4,4)+B(3,2)*B(4,3)*B(2,4)+B(4,2)*B(2,3)*B(3,4)-&
     
             B(4,2)*B(3,3)*B(2,4)-B(4,3)*B(3,4)*B(2,2)-B(4,4)*B(3,2)*B(2,3)
             
   adjB(2)=B(2,1)*B(3,3)*B(4,4)+B(3,1)*B(4,3)*B(2,4)+B(4,1)*B(2,3)*B(3,4)-&
                     
             B(4,1)*B(3,3)*B(2,4)-B(4,3)*B(3,4)*B(2,1)-B(4,4)*B(3,1)*B(2,3)
             
   adjB(2)=-adjB(2)
             
   adjB(3)=B(2,1)*B(3,2)*B(4,4)+B(3,1)*B(4,2)*B(2,4)+B(4,1)*B(2,2)*B(3,4)-&
                     
             B(4,1)*B(3,2)*B(2,4)-B(4,2)*B(3,4)*B(2,1)-B(4,4)*B(3,1)*B(2,2)
                   
   adjB(4)=B(2,1)*B(3,2)*B(4,3)+B(3,1)*B(4,2)*B(2,3)+B(4,1)*B(2,2)*B(3,3)-&
                     
             B(4,1)*B(3,2)*B(2,3)-B(4,2)*B(3,3)*B(2,1)-B(4,3)*B(3,1)*B(2,2)
           
   adjB(4)=-adjB(4)
             
   detB=B(1,1)*adjB(1)+B(1,2)*adjB(2)+B(1,3)*adjB(3)+B(1,4)*adjB(4)
             
   muhat(i)=(adjB(1)*Ky(1)+adjB(2)*Ky(2)+adjB(3)*Ky(3)+adjB(4)*Ky(4))/detB

end do


end Subroutine SpTeLLKS    