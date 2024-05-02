      subroutine distance(x, y, z, ilower, iupper, jlower, jupper,
     $     klower, kupper, il, jl, kl, xx, yy, zz, dist)
c     
c     compute the distance from inner cell to wall, used for turbulence viscosity
c     
      implicit none
c     
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl

      double precision, intent(in), dimension(ilower+1:iupper,
     $     jlower+1:jupper, klower+1:kupper)::
     $     x, y, z

      double precision, intent(out)::
     $     dist(6,il,jl,kl)

      double precision, dimension(0:il+1,0:jl+1,0:kl+1), intent(out)::
     $     xx, yy, zz
c     
c     LOCAL VARIABLES
c     
      integer::
     $     i, j, k,
     $     ipos, jpos, kpos

      double precision::
     $     len,
     $     dx1, dy1, dz1, dx2, dy2, dz2,
     $     nx, ny, nz

      integer::
     $     method = 2           ! 1: vector dot product
                                ! 2: straight connection
c     
c     *** START SUBROUTINE
c     
c     cell centers
c
      xx(0,1:jl,1:kl) = .25d0*(x(1,1:jl,1:kl)+x(1,2:jl+1,1:kl)
     $     +x(1,1:jl,2:kl+1)+x(1,2:jl+1,2:kl+1))  
      yy(0,1:jl,1:kl) = .25d0*(y(1,1:jl,1:kl)+y(1,2:jl+1,1:kl)
     $     +y(1,1:jl,2:kl+1)+y(1,2:jl+1,2:kl+1))  
      zz(0,1:jl,1:kl) = .25d0*(z(1,1:jl,1:kl)+z(1,2:jl+1,1:kl)
     $     +z(1,1:jl,2:kl+1)+z(1,2:jl+1,2:kl+1))  

      xx(il+1,1:jl,1:kl)=.25d0*(x(il+1,1:jl,1:kl)+x(il+1,2:jl+1,1:kl)
     $     +x(il+1,1:jl,2:kl+1)+x(il+1,2:jl+1,2:kl+1))  
      yy(il+1,1:jl,1:kl)=.25d0*(y(il+1,1:jl,1:kl)+y(il+1,2:jl+1,1:kl)
     $     +y(il+1,1:jl,2:kl+1)+y(il+1,2:jl+1,2:kl+1))  
      zz(il+1,1:jl,1:kl)=.25d0*(z(il+1,1:jl,1:kl)+z(il+1,2:jl+1,1:kl)
     $     +z(il+1,1:jl,2:kl+1)+z(il+1,2:jl+1,2:kl+1))  

      xx(1:il,0,1:kl) = .25d0*(x(1:il,1,1:kl)+x(2:il+1,1,1:kl)
     $     +x(1:il,1,2:kl+1)+x(2:il+1,1,2:kl+1))  
      yy(1:il,0,1:kl) = .25d0*(y(1:il,1,1:kl)+y(2:il+1,1,1:kl)
     $     +y(1:il,1,2:kl+1)+y(2:il+1,1,2:kl+1))  
      zz(1:il,0,1:kl) = .25d0*(z(1:il,1,1:kl)+z(2:il+1,1,1:kl)
     $     +z(1:il,1,2:kl+1)+z(2:il+1,1,2:kl+1))  

      xx(1:il,jl+1,1:kl)=.25d0*(x(1:il,jl+1,1:kl)+x(2:il+1,jl+1,1:kl)
     $     +x(1:il,jl+1,2:kl+1)+x(2:il+1,jl+1,2:kl+1))  
      yy(1:il,jl+1,1:kl)=.25d0*(y(1:il,jl+1,1:kl)+y(2:il+1,jl+1,1:kl)
     $     +y(1:il,jl+1,2:kl+1)+y(2:il+1,jl+1,2:kl+1))  
      zz(1:il,jl+1,1:kl)=.25d0*(z(1:il,jl+1,1:kl)+z(2:il+1,jl+1,1:kl)
     $     +z(1:il,jl+1,2:kl+1)+z(2:il+1,jl+1,2:kl+1))  

      xx(1:il,1:jl,0) = .25d0*(x(1:il,1:jl,1)+x(2:il+1,1:jl,1)
     $     +x(1:il,2:jl+1,1)+x(2:il+1,2:jl+1,1))  
      yy(1:il,1:jl,0) = .25d0*(y(1:il,1:jl,1)+y(2:il+1,1:jl,1)
     $     +y(1:il,2:jl+1,1)+y(2:il+1,2:jl+1,1))  
      zz(1:il,1:jl,0) = .25d0*(z(1:il,1:jl,1)+z(2:il+1,1:jl,1)
     $     +z(1:il,2:jl+1,1)+z(2:il+1,2:jl+1,1))  

      xx(1:il,1:jl,kl+1)=.25d0*(x(1:il,1:jl,kl+1)+x(2:il+1,1:jl,kl+1)
     $     +x(1:il,2:jl+1,kl+1)+x(2:il+1,2:jl+1,kl+1))  
      yy(1:il,1:jl,kl+1)=.25d0*(y(1:il,1:jl,kl+1)+y(2:il+1,1:jl,kl+1)
     $     +y(1:il,2:jl+1,kl+1)+y(2:il+1,2:jl+1,kl+1))  
      zz(1:il,1:jl,kl+1)=.25d0*(z(1:il,1:jl,kl+1)+z(2:il+1,1:jl,kl+1)
     $     +z(1:il,2:jl+1,kl+1)+z(2:il+1,2:jl+1,kl+1))  

      xx(1:il,1:jl,1:kl)=.125d0*
     $     (x(1:il,1:jl,1:kl)+x(2:il+1,1:jl,1:kl)
     $     +x(1:il,2:jl+1,1:kl)+x(2:il+1,2:jl+1,1:kl)
     $     +x(1:il,1:jl,2:kl+1)+x(2:il+1,1:jl,2:kl+1)
     $     +x(1:il,2:jl+1,2:kl+1)+x(2:il+1,2:jl+1,2:kl+1))
      yy(1:il,1:jl,1:kl) = .125d0*
     $     (y(1:il,1:jl,1:kl)+y(2:il+1,1:jl,1:kl)
     $     +y(1:il,2:jl+1,1:kl)+y(2:il+1,2:jl+1,1:kl)
     $     +y(1:il,1:jl,2:kl+1)+y(2:il+1,1:jl,2:kl+1)
     $     +y(1:il,2:jl+1,2:kl+1)+y(2:il+1,2:jl+1,2:kl+1))
      zz(1:il,1:jl,1:kl) = .125d0*
     $     (z(1:il,1:jl,1:kl)+z(2:il+1,1:jl,1:kl)
     $     +z(1:il,2:jl+1,1:kl)+z(2:il+1,2:jl+1,1:kl)
     $     +z(1:il,1:jl,2:kl+1)+z(2:il+1,1:jl,2:kl+1)
     $     +z(1:il,2:jl+1,2:kl+1)+z(2:il+1,2:jl+1,2:kl+1))


      if(method .eq. 1) then

c     to xi-lower surface
         ipos = 0
         do j = 1, jl
            do k = 1, kl

               dx1 = x(ipos,j+1,k+1)-x(ipos,j,k)
               dx2 = x(ipos,j,k+1)-x(ipos,j+1,k)
               dy1 = y(ipos,j+1,k+1)-y(ipos,j,k)
               dy2 = y(ipos,j,k+1)-y(ipos,j+1,k)
               dz1 = z(ipos,j+1,k+1)-z(ipos,j,k)
               dz2 = z(ipos,j,k+1)-z(ipos,j+1,k)
               nx = dy1*dz2-dy2*dz1
               ny = dx2*dz1-dx1*dz2
               nz = dx1*dy2-dx2*dy1
               len = dsqrt(nx*nx+ny*ny+nz*nz)
               nx = nx/len
               ny = ny/len
               nz = nz/len
               
               do i = 1, il
                  dx1 = xx(i,j,k)-xx(0,j,k)
                  dy1 = yy(i,j,k)-yy(0,j,k)
                  dz1 = zz(i,j,k)-zz(0,j,k)
                  dist(1,i,j,k) = abs(dx1*nx+dy1*ny+dz1*nz)
               end do
               
            end do
         end do

c     to xi-upper surface
         ipos = il+1
         do j = 1, jl
            do k = 1, kl

               dx1 = x(ipos,j+1,k+1)-x(ipos,j,k)
               dx2 = x(ipos,j,k+1)-x(ipos,j+1,k)
               dy1 = y(ipos,j+1,k+1)-y(ipos,j,k)
               dy2 = y(ipos,j,k+1)-y(ipos,j+1,k)
               dz1 = z(ipos,j+1,k+1)-z(ipos,j,k)
               dz2 = z(ipos,j,k+1)-z(ipos,j+1,k)
               nx = dy1*dz2-dy2*dz1
               ny = dx2*dz1-dx1*dz2
               nz = dx1*dy2-dx2*dy1
               len = dsqrt(nx*nx+ny*ny+nz*nz)
               nx = nx/len
               ny = ny/len
               nz = nz/len
               
               do i = 1, il
                  dx1 = xx(il+1,j,k)-xx(i,j,k)
                  dy1 = yy(il+1,j,k)-yy(i,j,k)
                  dz1 = zz(il+1,j,k)-zz(i,j,k)
                  dist(2,i,j,k) = abs(dx1*nx+dy1*ny+dz1*nz)
               end do
               
            end do
         end do


c     to eta-lower surface
         jpos = 0
         do i = 1, il
            do k = 1, kl

               dx1 = x(i+1,jpos,k+1)-x(i,jpos,k)
               dx2 = x(i,jpos,k+1)-x(i+1,jpos,k)
               dy1 = y(i+1,jpos,k+1)-y(i,jpos,k)
               dy2 = y(i,jpos,k+1)-y(i+1,jpos,k)
               dz1 = z(i+1,jpos,k+1)-z(i,jpos,k)
               dz2 = z(i,jpos,k+1)-z(i+1,jpos,k)
               nx = dy1*dz2-dy2*dz1
               ny = dx2*dz1-dx1*dz2
               nz = dx1*dy2-dx2*dy1
               len = dsqrt(nx*nx+ny*ny+nz*nz)
               nx = nx/len
               ny = ny/len
               nz = nz/len
               
               do j = 1, jl
                  dx1 = xx(i,j,k)-xx(i,0,k)
                  dy1 = yy(i,j,k)-yy(i,0,k)
                  dz1 = zz(i,j,k)-zz(i,0,k)
                  dist(3,i,j,k) = abs(dx1*nx+dy1*ny+dz1*nz)
               end do
               
            end do
         end do

c     to eta-upper surface

         jpos = jl+1
         do i = 1, il
            do k = 1, kl

               dx1 = x(i+1,jpos,k+1)-x(i,jpos,k)
               dx2 = x(i,jpos,k+1)-x(i+1,jpos,k)
               dy1 = y(i+1,jpos,k+1)-y(i,jpos,k)
               dy2 = y(i,jpos,k+1)-y(i+1,jpos,k)
               dz1 = z(i+1,jpos,k+1)-z(i,jpos,k)
               dz2 = z(i,jpos,k+1)-z(i+1,jpos,k)
               nx = dy1*dz2-dy2*dz1
               ny = dx2*dz1-dx1*dz2
               nz = dx1*dy2-dx2*dy1
               len = dsqrt(nx*nx+ny*ny+nz*nz)
               nx = nx/len
               ny = ny/len
               nz = nz/len
               
               do j = 1, jl
                  dx1 = xx(i,jl+1,k)-xx(i,j,k)
                  dy1 = yy(i,jl+1,k)-yy(i,j,k)
                  dz1 = zz(i,jl+1,k)-zz(i,j,k)
                  dist(4,i,j,k) = abs(dx1*nx+dy1*ny+dz1*nz)
               end do
               
            end do
         end do

c     to zeta-lower surface

         kpos = 0
         do i = 1, il
            do j= 1, jl

               dx1 = x(i+1,j+1,kpos)-x(i,j,kpos)
               dx2 = x(i,j+1,kpos)-x(i+1,j,kpos)
               dy1 = y(i+1,j+1,kpos)-y(i,j,kpos)
               dy2 = y(i,j+1,kpos)-y(i+1,j,kpos)
               dz1 = z(i+1,j+1,kpos)-z(i,j,kpos)
               dz2 = z(i,j+1,kpos)-z(i+1,j,kpos)
               nx = dy1*dz2-dy2*dz1
               ny = dx2*dz1-dx1*dz2
               nz = dx1*dy2-dx2*dy1
               len = dsqrt(nx*nx+ny*ny+nz*nz)
               nx = nx/len
               ny = ny/len
               nz = nz/len
               
               do k = 1, kl
                  dx1 = xx(i,j,k)-xx(i,j,0)
                  dy1 = yy(i,j,k)-yy(i,j,0)
                  dz1 = zz(i,j,k)-zz(i,j,0)
                  dist(5,i,j,k) = abs(dx1*nx+dy1*ny+dz1*nz)
               end do
               
            end do
         end do

c     to zeta-upper surface

         kpos = kl+1
         do i = 1, il
            do j= 1, jl

               dx1 = x(i+1,j+1,kpos)-x(i,j,kpos)
               dx2 = x(i,j+1,kpos)-x(i+1,j,kpos)
               dy1 = y(i+1,j+1,kpos)-y(i,j,kpos)
               dy2 = y(i,j+1,kpos)-y(i+1,j,kpos)
               dz1 = z(i+1,j+1,kpos)-z(i,j,kpos)
               dz2 = z(i,j+1,kpos)-z(i+1,j,kpos)
               nx = dy1*dz2-dy2*dz1
               ny = dx2*dz1-dx1*dz2
               nz = dx1*dy2-dx2*dy1
               len = dsqrt(nx*nx+ny*ny+nz*nz)
               nx = nx/len
               ny = ny/len
               nz = nz/len
               
               do k = 1, kl
                  dx1 = xx(i,j,kl+1)-xx(i,j,k)
                  dy1 = yy(i,j,kl+1)-yy(i,j,k)
                  dz1 = zz(i,j,kl+1)-zz(i,j,k)
                  dist(6,i,j,k) = abs(dx1*nx+dy1*ny+dz1*nz)
               end do
               
            end do
         end do
      else
         do i = 1, il
            dist(1,i,1:jl,1:kl) = dsqrt(
     $           (xx(i,1:jl,1:kl)-xx(0,1:jl,1:kl))**2
     $           +(yy(i,1:jl,1:kl)-yy(0,1:jl,1:kl))**2
     $           +(zz(i,1:jl,1:kl)-zz(0,1:jl,1:kl))**2)
            dist(2,i,1:jl,1:kl) = dsqrt(
     $           (xx(i,1:jl,1:kl)-xx(il+1,1:jl,1:kl))**2
     $           +(yy(i,1:jl,1:kl)-yy(il+1,1:jl,1:kl))**2
     $           +(zz(i,1:jl,1:kl)-zz(il+1,1:jl,1:kl))**2)
         end do

         do j = 1, jl
            dist(3, 1:il,j,1:kl) = dsqrt(
     $           (xx(1:il,j,1:kl)-xx(1:il,0,1:kl))**2
     $           +(yy(1:il,j,1:kl)-yy(1:il,0,1:kl))**2
     $           +(zz(1:il,j,1:kl)-zz(1:il,0,1:kl))**2)
            dist(4, 1:il,j,1:kl) = dsqrt(
     $           (xx(1:il,j,1:kl)-xx(1:il,jl+1,1:kl))**2
     $           +(yy(1:il,j,1:kl)-yy(1:il,jl+1,1:kl))**2
     $           +(zz(1:il,j,1:kl)-zz(1:il,jl+1,1:kl))**2)
         end do

         do k = 1, kl
            dist(5, 1:il,1:jl,k) = dsqrt(
     $           (xx(1:il,1:jl,k)-xx(1:il,1:kl,0))**2
     $           +(yy(1:il,1:jl,k)-yy(1:il,1:kl,0))**2
     $           +(zz(1:il,1:jl,k)-zz(1:il,1:kl,0))**2)
            dist(6, 1:il,1:jl,k) = dsqrt(
     $           (xx(1:il,1:jl,k)-xx(1:il,1:kl,kl+1))**2
     $           +(yy(1:il,1:jl,k)-yy(1:il,1:kl,kl+1))**2
     $           +(zz(1:il,1:jl,k)-zz(1:il,1:kl,kl+1))**2)
         end do

      end if

      end subroutine distance
