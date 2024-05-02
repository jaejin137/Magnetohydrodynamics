      subroutine vorticity(il, jl, kl, nl, ilower, iupper,
     $         jlower, jupper, klower, kupper, x, y, z, q, vol,
     $         lx, ly, lz, mx, my, mz, nx, ny, nz, u, v, w, omega)
c
c     compute vorticity
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     il, jl, kl, nl

      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)

      double precision, intent(in)::
     $     vol(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)

      double precision, intent(in), dimension(ilower+1:iupper,
     $     jlower+1:jupper, klower+1:kupper)::
     $     x, y, z

      double precision, intent(in),
     $     dimension(il+1,jl+1,kl+1)::
     $     lx, ly, lz, mx, my, mz, nx, ny, nz

      double precision, intent(in),
     $     dimension(ilower:iupper, jlower:jupper, klower:kupper)::
     $     u, v, w             ! velocity

      double precision, dimension(il,jl,kl), intent(out)::
     $     omega                ! vorticity
c
c     LOCAL VARIABLES
c
      integer:: i, j, k
      double precision::
     $     dudxi, dudeta, dudzeta,
     $     dvdxi, dvdeta, dvdzeta,
     $     dwdxi, dwdeta, dwdzeta,
     $     dudx, dudy, dudz,
     $     dvdx, dvdy, dvdz,
     $     dwdx, dwdy, dwdz
c
c *** START SUBROUTINE
c
      do i = 1, il
        do j = 1, jl
          do k = 1, kl
            dudxi = .5d0*(u(i+1,j,k)-u(i-1,j,k))
            dvdxi = .5d0*(v(i+1,j,k)-v(i-1,j,k))
            dwdxi = .5d0*(w(i+1,j,k)-w(i-1,j,k))

            dudeta = .5d0*(u(i,j+1,k)-u(i,j-1,k))
            dvdeta = .5d0*(v(i,j+1,k)-v(i,j-1,k))
            dwdeta = .5d0*(w(i,j+1,k)-w(i,j-1,k))

            dudzeta = .5d0*(u(i,j,k+1)-u(i,j,k-1))
            dvdzeta = .5d0*(v(i,j,k+1)-v(i,j,k-1))
            dwdzeta = .5d0*(w(i,j,k+1)-w(i,j,k-1))

c           dudx = dudxi*lx(i,j,k) + dudeta*mx(i,j,k)
c    $           + dudzeta*nx(i,j,k)
            dudy = dudxi*ly(i,j,k) + dudeta*my(i,j,k)
     $           + dudzeta*ny(i,j,k)
            dudz = dudxi*lz(i,j,k) + dudeta*mz(i,j,k)
     $           + dudzeta*nz(i,j,k)

            dvdx = dvdxi*lx(i,j,k) + dvdeta*mx(i,j,k)
     $           + dvdzeta*nx(i,j,k)
c           dvdy = dvdxi*ly(i,j,k) + dvdeta*my(i,j,k)
c    $           + dvdzeta*ny(i,j,k)
            dvdz = dvdxi*lz(i,j,k) + dvdeta*mz(i,j,k)
     $           + dvdzeta*nz(i,j,k)

            dwdx = dwdxi*lx(i,j,k) + dwdeta*mx(i,j,k)
     $           + dwdzeta*nx(i,j,k)
            dwdy = dwdxi*ly(i,j,k) + dwdeta*my(i,j,k)
     $           + dwdzeta*ny(i,j,k)
c           dwdz = dwdxi*lz(i,j,k) + dwdeta*mz(i,j,k)
c    $           + dwdzeta*nz(i,j,k)
            omega(i,j,k) = dsqrt(
     $       (dudy-dvdx)**2+(dvdz-dwdy)**2+(dwdx-dudz)**2)/vol(i,j,k)
          end do
        end do
      end do
c
      return
      end 
