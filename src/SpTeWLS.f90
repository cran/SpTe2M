Subroutine SpTeWLS(y, t, sx, sy, n, m, MAXm, NUM0, ht, hs, ht0, hs0, gt0, gs0, stE, NUM, muhat)

implicit none

integer :: n, m(n), MAXm, NUM0, NUM, SNUM(n), i, j, ii, jj, k, kk, Nullty, Ifault, &
           INDEX(NUM0)

double precision :: ht, hs, ht0, hs0, gt0, gs0, kt1, ks1, kt2, ks2, tmp, detB, eps, wt

double precision :: y(1:n,1:MAXm), t(n), sx(1:n,1:MAXm), sy(1:n,1:MAXm), &
                    stE(NUM,3), muhat(NUM), B(4,4), Ky(4), adjB(4), muhat0(NUM0), &
                    stE0(NUM0,3), res(1:n,1:MAXm)

double precision :: stE1(NUM0,3),y1(NUM0)

real(8) :: W(NUM0),SigmaI(NUM0,NUM0), Sigma(NUM0,NUM0), COV(NUM0**2), COVI(NUM0**2)

eps=0D0


SNUM(1)=0; k=0

do i=1,n-1

   k=k+m(i)
   
   SNUM(i+1)=k

end do

k=0

do i=1,n
do j=1,m(i)

   k=k+1
   stE0(k,1)=sx(i,j)
   stE0(k,2)=sy(i,j)
   stE0(k,3)=t(i)

end do
end do

call SpTeLLKS(y,t,sx,sy,n,m,MAXm,ht0,hs0,stE0,NUM0,eps,muhat0)

k=0

do i=1,n
do j=1,m(i)

   k=k+1
   res(i,j)=y(i,j)-muhat0(k)

end do
end do

call SpTeWME(res,t,sx,sy,n,m,MAXm,gt0,gs0,stE0,NUM0,stE0,NUM0,Sigma)

do i=1,NUM
    
   muhat(i)=0D0

   Ky(1)=0D0; Ky(2)=0D0;  Ky(3)=0D0;  Ky(4)=0D0
             
   B(1,1)=0D0; B(1,2)=0D0; B(1,3)=0D0; B(1,4)=0D0
                        
   B(2,1)=0D0; B(2,2)=0D0; B(2,3)=0D0; B(2,4)=0D0
                                
   B(3,1)=0D0; B(3,2)=0D0; B(3,3)=0D0; B(3,4)=0D0

   B(4,1)=0D0; B(4,2)=0D0; B(4,3)=0D0; B(4,4)=0D0

   k=0
 
   do ii=1,n

      tmp=(t(ii)-stE(i,3))/ht

      if(tmp >= -1D0 .AND. tmp <= 1D0) then         
                 
         do jj=1,m(ii)
    
            tmp=dsqrt(dble(sx(ii,jj)-stE(i,1))**2D0+dble(sy(ii,jj)-stE(i,2))**2D0)/hs
          
            if(tmp >= -1D0 .AND. tmp <= 1D0) then
                  
               k=k+1; stE1(k,1)=sx(ii,jj); stE1(k,2)=sy(ii,jj)

               stE1(k,3)=t(ii); y1(k)=y(ii,jj)

               INDEX(k)=SNUM(ii)+jj
                  
            end if   

         end do

      end if

   end do

   kk=0

   do ii=1,k

   do jj=1,ii

      kk=kk+1; COV(kk)=Sigma(INDEX(ii),INDEX(jj))

   end do
 
   end do

   call syminv(COV(1:kk),k,COVI(1:kk),W(1:k),Nullty,Ifault)

   do ii=1,k

   do jj=1,ii

      SigmaI(ii,jj)=COVI(jj+ii*(ii-1)/2) 

   end do

   end do

   do ii=1,k-1

   do jj=ii+1,k

      SigmaI(ii,jj)=SigmaI(jj,ii) 

   end do

   end do

   
   do ii=1,k
 
      tmp=(stE1(ii,3)-stE(i,3))/ht

      kt1=dsqrt(0.75D0*(1D0-tmp**2D0))

      tmp=dsqrt(dble(stE1(ii,1)-stE(i,1))**2D0+dble(stE1(ii,2)-stE(i,2))**2D0)/hs
      
      ks1=dsqrt(0.75D0*(1D0-tmp**2D0)) 

      do jj=1,k

         tmp=(stE1(jj,3)-stE(i,3))/ht

         kt2=dsqrt(0.75D0*(1D0-tmp**2D0))

         tmp=dsqrt(dble(stE1(jj,1)-stE(i,1))**2D0+dble(stE1(jj,2)-stE(i,2))**2D0)/hs
       
         ks2=dsqrt(0.75D0*(1D0-tmp**2D0))   

         wt=kt1*ks1*kt2*ks2*SigmaI(ii,jj)

         Ky(1)=Ky(1)+y1(jj)*wt

         Ky(2)=Ky(2)+(stE1(ii,3)-stE(i,3))*y1(jj)*wt
                       
         Ky(3)=Ky(3)+(stE1(ii,1)-stE(i,1))*y1(jj)*wt
                              
         Ky(4)=Ky(4)+(stE1(ii,2)-stE(i,2))*y1(jj)*wt
 
                             
         B(1,1)=B(1,1)+wt

         B(1,2)=B(1,2)+(stE1(jj,3)-stE(i,3))*wt
                      
         B(1,3)=B(1,3)+(stE1(jj,1)-stE(i,1))*wt
                          
         B(1,4)=B(1,4)+(stE1(jj,2)-stE(i,2))*wt 
                    
         B(2,2)=B(2,2)+(stE1(ii,3)-stE(i,3))*(stE1(jj,3)-stE(i,3))*wt 
                    
         B(2,3)=B(2,3)+(stE1(ii,3)-stE(i,3))*(stE1(jj,1)-stE(i,1))*wt
                     
         B(2,4)=B(2,4)+(stE1(ii,3)-stE(i,3))*(stE1(jj,2)-stE(i,2))*wt

         B(3,3)=B(3,3)+(stE1(ii,1)-stE(i,1))*(stE1(jj,1)-stE(i,1))*wt
                      
         B(3,4)=B(3,4)+(stE1(ii,1)-stE(i,1))*(stE1(jj,2)-stE(i,2))*wt 
                     
         B(4,4)=B(4,4)+(stE1(ii,2)-stE(i,2))*(stE1(jj,2)-stE(i,2))*wt

         B(2,1)=B(1,2)
                 
         B(3,1)=B(1,3)
                          
         B(3,2)=B(2,3)
                      
         B(4,1)=B(1,4)
                      
         B(4,2)=B(2,4)
                     
         B(4,3)=B(3,4)
 
      end do

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
              
 
end Subroutine SpTeWLS

   