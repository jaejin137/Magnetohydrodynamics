      subroutine metric_all(x, y, z, lx, ly, lz, mx, my, mz, nx, ny, nz,
     $     il, jl, kl, ilower, iupper, jlower, jupper, klower, kupper)
c
c     compute matric for all cell centers in domain
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     il, jl, kl,
     $     ilower, iupper, jlower, jupper, klower, kupper
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x,  y, z
c
      double precision, intent(out), dimension(il+1,jl+1,kl+1)::
     $     lx, ly, lz, mx, my, mz, nx, ny, nz

      integer::
     $     i, j, k, ip, jp, kp
c
c     LOCAL VARIABLES
c
      double precision::
     $     dx1, dx2, dy1, dy2, dz1, dz2
c
c *** START SUBROUTINE ***
c
      ! lx, ly, lz
      
      lx = 0.d0
      ly = 0.d0
      lz = 0.d0

      do i = 1, il+1
         do j = 1, jl
            do k = 1, kl

               jp = j + 1
               kp = k + 1

               dx1 = x(i,jp,kp) - x(i,j,k)
               dy1 = y(i,jp,kp) - y(i,j,k)
               dz1 = z(i,jp,kp) - z(i,j,k)
               dx2 = x(i,j,kp)  - x(i,jp,k)
               dy2 = y(i,j,kp)  - y(i,jp,k)
               dz2 = z(i,j,kp)  - z(i,jp,k)

               lx(i,j,k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
               ly(i,j,k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
               lz(i,j,k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            end do
         end do
      end do

      lx(:il,:,:) = .5d0 * (lx(:il,:,:) + lx(2:il+1,:,:))
      lx(il+1,:,:) = 0.d0
      ly(:il,:,:) = .5d0 * (ly(:il,:,:) + ly(2:il+1,:,:))
      ly(il+1,:,:) = 0.d0
      lz(:il,:,:) = .5d0 * (lz(:il,:,:) + lz(2:il+1,:,:))
      lz(il+1,:,:) = 0.d0

      ! mx, my, mz
      
      mx = 0.d0
      my = 0.d0
      mz = 0.d0

      do j = 1, jl+1
         do i = 1, il
            do k = 1, kl

               ip = i + 1
               kp = k + 1

               dx1 = x(ip,j,kp) - x(i,j,k)
               dy1 = y(ip,j,kp) - y(i,j,k)
               dz1 = z(ip,j,kp) - z(i,j,k)
               dx2 = x(ip,j,k)  - x(i,j,kp)
               dy2 = y(ip,j,k)  - y(i,j,kp)
               dz2 = z(ip,j,k)  - z(i,j,kp)

               mx (i,j,k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
               my (i,j,k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
               mz (i,j,k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            end do
         end do
      end do

      mx(:,:jl,:) = .5d0 * (mx(:,:jl,:) + mx(:,2:jl+1,:))
      mx(:,jl+1,:) = 0.d0
      my(:,:jl,:) = .5d0 * (my(:,:jl,:) + my(:,2:jl+1,:))
      my(:,jl+1,:) = 0.d0
      mz(:,:jl,:) = .5d0 * (mz(:,:jl,:) + mz(:,2:jl+1,:))
      mz(:,jl+1,:) = 0.d0

      ! nx, ny, nz
      
      nx = 0.d0
      ny = 0.d0
      nz = 0.d0

      do k = 1, kl+1
         do i = 1, il
            do j = 1, jl

               ip = i + 1
               jp = j + 1

               dx1 = x(ip,jp,k) - x(i,j,k)
               dy1 = y(ip,jp,k) - y(i,j,k)
               dz1 = z(ip,jp,k) - z(i,j,k)
               dx2 = x(i,jp,k)  - x(ip,j,k)
               dy2 = y(i,jp,k)  - y(ip,j,k)
               dz2 = z(i,jp,k)  - z(ip,j,k)

               nx (i,j,k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
               ny (i,j,k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
               nz (i,j,k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            end do
         end do
      end do

      nx(:,:,:kl) = .5d0 * (nx(:,:,:kl) + nx(:,:,2:kl+1))
      nx(:,:,kl+1) = 0.d0
      ny(:,:,:kl) = .5d0 * (ny(:,:,:kl) + ny(:,:,2:kl+1))
      ny(:,:,kl+1) = 0.d0
      nz(:,:,:kl) = .5d0 * (nz(:,:,:kl) + nz(:,:,2:kl+1))
      nz(:,:,kl+1) = 0.d0

      end subroutine metric_all
