      subroutine spvel(index, i1, i2, dim2, il, jl, kl, x, y, z,
     $     xo, yo, zo, ilower, iupper, jlower, jupper,
     $     klower, kupper, tintvl, qt)
c     
c     This subroutine is used to calculate the grid velocity due to
c     the motion of moving grid.
c     
c     qt(i,1) = -xt, qt(i,2) = -yt, qt(i,3) = -zt
c     
      implicit none
c     
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     dim2,                ! maxium face number in 3 directions
     $     i1, i2,              ! node index in two directions
     $     il, jl, kl,          ! cell number in 3 directions
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     index

      double precision, intent(in):: 
     $     tintvl

      double precision, dimension(ilower:dim2,3), intent(out):: qt

      double precision, intent(inout),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z, xo, yo, zo

c     LOCAL VARIABLES
c     
      double precision::
     $     tinv, dxdt, dydt, dzdt

      integer::
     $     i, j, k,             ! cell iteration index
     $     ilp, jlp, klp        ! il+1, jl+1, kl+1
c     
c *** SUBROUTINE START ***
c     
c *** set up some parameters
c     
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c     
c *** calculate the grid moving velocity at each node point
c     
      tinv = -0.125d0/tintvl
c     
      IF  ( INDEX.EQ.1 )  THEN
c     
c     ... xi direction
c     
         j = i1
         k = i2
c     
         do i = 1,il
c     
            dxdt =  x(i,j,k)+ x(i,j+1,k)+ x(i,j,k+1)+ x(i,j+1,k+1)+
     >           x(i+1,j,k)+ x(i+1,j+1,k)+ x(i+1,j,k+1)+ x(i+1,j+1,k+1)
     >           -xo(i,j,k)-xo(i,j+1,k)-xo(i,j,k+1)-xo(i,j+1,k+1)-
     >           xo(i+1,j,k)-xo(i+1,j+1,k)-xo(i+1,j,k+1)-xo(i+1,j+1,k+1)

            dydt =  y(i,j,k)+ y(i,j+1,k)+ y(i,j,k+1)+ y(i,j+1,k+1)+
     >           y(i+1,j,k)+ y(i+1,j+1,k)+ y(i+1,j,k+1)+ y(i+1,j+1,k+1)
     >           -yo(i,j,k)-yo(i,j+1,k)-yo(i,j,k+1)-yo(i,j+1,k+1)-
     >           yo(i+1,j,k)-yo(i+1,j+1,k)-yo(i+1,j,k+1)-yo(i+1,j+1,k+1)

            dzdt =  z(i,j,k)+ z(i,j+1,k)+ z(i,j,k+1)+ z(i,j+1,k+1)+
     >           z(i+1,j,k)+ z(i+1,j+1,k)+ z(i+1,j,k+1)+ z(i+1,j+1,k+1)
     >           -zo(i,j,k)-zo(i,j+1,k)-zo(i,j,k+1)-zo(i,j+1,k+1)-
     >           zo(i+1,j,k)-zo(i+1,j+1,k)-zo(i+1,j,k+1)-zo(i+1,j+1,k+1)

            qt(i,1) = dxdt*tinv
            qt(i,2) = dydt*tinv
            qt(i,3) = dzdt*tinv

         end do
c     
      ELSE IF  ( INDEX.EQ.2 )  THEN
c     
c     ... eta direction
c     
         i = i1
         k = i2
c     
         do j = 1,jl
c     
            dxdt =  x(i,j,k)+ x(i,j+1,k)+ x(i,j,k+1)+ x(i,j+1,k+1)+
     >           x(i+1,j,k)+ x(i+1,j+1,k)+ x(i+1,j,k+1)+ x(i+1,j+1,k+1)
     >           -xo(i,j,k)-xo(i,j+1,k)-xo(i,j,k+1)-xo(i,j+1,k+1)-
     >           xo(i+1,j,k)-xo(i+1,j+1,k)-xo(i+1,j,k+1)-xo(i+1,j+1,k+1)
c     
            dydt =  y(i,j,k)+ y(i,j+1,k)+ y(i,j,k+1)+ y(i,j+1,k+1)+
     >           y(i+1,j,k)+ y(i+1,j+1,k)+ y(i+1,j,k+1)+ y(i+1,j+1,k+1)
     >           -yo(i,j,k)-yo(i,j+1,k)-yo(i,j,k+1)-yo(i,j+1,k+1)-
     >           yo(i+1,j,k)-yo(i+1,j+1,k)-yo(i+1,j,k+1)-yo(i+1,j+1,k+1)
c     
            dzdt =  z(i,j,k)+ z(i,j+1,k)+ z(i,j,k+1)+ z(i,j+1,k+1)+
     >           z(i+1,j,k)+ z(i+1,j+1,k)+ z(i+1,j,k+1)+ z(i+1,j+1,k+1)
     >           -zo(i,j,k)-zo(i,j+1,k)-zo(i,j,k+1)-zo(i,j+1,k+1)-
     >           zo(i+1,j,k)-zo(i+1,j+1,k)-zo(i+1,j,k+1)-zo(i+1,j+1,k+1)
c     
            qt(j,1) = dxdt*tinv
            qt(j,2) = dydt*tinv
            qt(j,3) = dzdt*tinv

         end do

      ELSE IF  ( INDEX.EQ.3 )  THEN
c     
c     ... zeta direction
c     
         i = i1
         j = i2
c     
         do k = 1,kl
c     
            dxdt =  x(i,j,k)+ x(i,j+1,k)+ x(i,j,k+1)+ x(i,j+1,k+1)+
     >           x(i+1,j,k)+ x(i+1,j+1,k)+ x(i+1,j,k+1)+ x(i+1,j+1,k+1)
     >           -xo(i,j,k)-xo(i,j+1,k)-xo(i,j,k+1)-xo(i,j+1,k+1)-
     >           xo(i+1,j,k)-xo(i+1,j+1,k)-xo(i+1,j,k+1)-xo(i+1,j+1,k+1)

            dydt =  y(i,j,k)+ y(i,j+1,k)+ y(i,j,k+1)+ y(i,j+1,k+1)+
     >           y(i+1,j,k)+ y(i+1,j+1,k)+ y(i+1,j,k+1)+ y(i+1,j+1,k+1)
     >           -yo(i,j,k)-yo(i,j+1,k)-yo(i,j,k+1)-yo(i,j+1,k+1)-
     >           yo(i+1,j,k)-yo(i+1,j+1,k)-yo(i+1,j,k+1)-yo(i+1,j+1,k+1)

            dzdt =  z(i,j,k)+ z(i,j+1,k)+ z(i,j,k+1)+ z(i,j+1,k+1)+
     >           z(i+1,j,k)+ z(i+1,j+1,k)+ z(i+1,j,k+1)+ z(i+1,j+1,k+1)
     >           -zo(i,j,k)-zo(i,j+1,k)-zo(i,j,k+1)-zo(i,j+1,k+1)-
     >           zo(i+1,j,k)-zo(i+1,j+1,k)-zo(i+1,j,k+1)-zo(i+1,j+1,k+1)

            qt(k,1) = dxdt*tinv
            qt(k,2) = dydt*tinv
            qt(k,3) = dzdt*tinv

         end do

      END IF
c     
      end
