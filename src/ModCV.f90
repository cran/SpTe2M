Subroutine ModCV(y, t, sx, sy, n, m, MAXm, NUM, ht, hs, Nbw, eps, mcv)

implicit none

integer :: n, m(n), MAXm, Nbw, i, j, ibw, k, NUM

double precision :: ht(Nbw), hs(Nbw), eps

double precision :: y(1:n,1:MAXm), t(n), sx(1:n,1:MAXm), sy(1:n,1:MAXm), &
                    muhat(NUM), stE(NUM,3), mcv(Nbw)


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

   mcv(ibw)=0D0

   call SpTeLLKS(y,t,sx,sy,n,m,MAXm,ht(ibw),hs(ibw),stE,NUM,eps,muhat)
 
   k=0

   do i=1,n
   do j=1,m(i)

      k=k+1
      mcv(ibw)=mcv(ibw)+(y(i,j)-muhat(k))**2D0

   end do
   end do

   mcv(ibw)=mcv(ibw)/NUM

end do


end Subroutine ModCV      