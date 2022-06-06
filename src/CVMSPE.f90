Subroutine CVMSPE(res, t, sx, sy, n, m, MAXm, NUM, gt, gs, Nbw, dt, ds, mspe)

implicit none

integer :: n, m(n), MAXm, NUM, SNUM(n), ibw, i, j, ii, jj, k, kk, Nullty, Ifault, &
           INDEX(NUM), Nbw

double precision :: gt(Nbw), gs(Nbw), dt, ds, tmp

double precision :: res(1:n,1:MAXm), t(n), sx(1:n,1:MAXm), sy(1:n,1:MAXm), &
                    stE(NUM,3), pred(1:n,1:MAXm), mspe(Nbw)

double precision :: stE1(NUM,3),res1(NUM)

real(8) :: W(NUM), SigmaI(NUM,NUM), Sigma(NUM,NUM), COV(NUM**2), COVI(NUM**2)


SNUM(1)=0; k=0

do i=1,n-1

   k=k+m(i)
   
   SNUM(i+1)=k

end do

k=0

do i=1,n
do j=1,m(i)

   k=k+1
   stE(k,1)=sx(i,j)
   stE(k,2)=sy(i,j)
   stE(k,3)=t(i)

end do
end do


do ibw=1,Nbw

   mspe(ibw)=0D0

   call SpTeWME(res,t,sx,sy,n,m,MAXm,gt(ibw),gs(ibw),stE,NUM,stE,NUM,Sigma)

   do i=1,n

   do j=1,m(i)

      pred(i,j)=0D0; k=0

      do ii=1,n

         tmp=dabs(t(ii)-t(i))

         if(tmp < dt) then         
                 
            do jj=1,m(ii)
    
               tmp=dsqrt(dble(sx(ii,jj)-sx(i,j))**2D0+dble(sy(ii,jj)-sy(i,j))**2D0)
             
               if(tmp < ds .AND. (dt .NE. 0D0 .OR. ds .NE. 0D0)) then

                  k=k+1; stE1(k,1)=sx(ii,jj); stE1(k,2)=sy(ii,jj)

                  stE1(k,3)=t(ii); res1(k)=res(ii,jj)

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
   
         do jj=1,k

            pred(i,j)=pred(i,j)+Sigma(INDEX(ii),SNUM(i)+j)*SigmaI(ii,jj)*res1(jj)

         end do

      end do

   mspe(ibw)=mspe(ibw)+(res(i,j)-pred(i,j))**2D0
  
   end do

   end do  

   mspe(ibw)=mspe(ibw)/dble(NUM)
    
end do 


end Subroutine CVMSPE
