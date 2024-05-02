      subroutine wall_dist(il, jl, kl, ides, ipos,
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     bc_xie_lower, bc_xie_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zta_lower, bc_zta_upper,
     $     x, y, z, xw, yw, zw, cdes, dist)
c     
c     compute the shortest distance from inner cell to wall
c     
      implicit none
c     
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl, ides, ipos

      double precision, intent(in), dimension(ilower+1:iupper,
     $     jlower+1:jupper, klower+1:kupper)::
     $     x, y, z

      integer, intent(in)::
     $     bc_xie_lower(jl, kl), bc_xie_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zta_lower(il, jl), bc_zta_upper(il, jl)

      double precision, intent(in), dimension(ipos)::
     $     xw, yw, zw

      double precision, intent(in):: cdes

      double precision, intent(out)::
     $     dist(il,jl,kl)

      double precision, dimension(0:il+1,0:jl+1,0:kl+1)::
     $     xx, yy, zz
c     
c     LOCAL VARIABLES
c     
      integer::
     $     i, j, k,
     $     ii, jj, kk,
     $     i1, j1, k1

      double precision::
     $     dist0, dx, dy, dz, dd(4), ddi, ddj, ddk, ddd
c     
c *** START SUBROUTINE
c     
c *** coordinates of the cell centers
c
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
c
c *** calculate the nearest wall distence
c
      do i = 1, il
        do j = 1, jl
          do k = 1, kl
            dist(i,j,k) = 1.d+12
          end do
        end do
      end do
c
      do i = 1, il
        do j = 1, jl
          do k = 1, kl
            do ii = 1, ipos
              dx = xw(ii)-xx(i,j,k)
              dy = yw(ii)-yy(i,j,k)
              dz = zw(ii)-zz(i,j,k)
              dist0 = dsqrt(dx**2+dy**2+dz**2)
              if (dist0.lt.dist(i,j,k)) dist(i,j,k) = dist0
            end do
          end do
        end do
      end do
c
c *** modify distance for DES simulation
c
      if (ides.ne.0) then
        do k = 1,kl
          do j = 1,jl
            do i = 1,il
              i1 = i+1
              j1 = j+1
              k1 = k+1
              dx = x(i1,j,k)-x(i,j,k)
              dy = y(i1,j,k)-y(i,j,k)
              dz = z(i1,j,k)-z(i,j,k)
              dd(1) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i1,j1,k)-x(i,j1,k)
              dy = y(i1,j1,k)-y(i,j1,k)
              dz = z(i1,j1,k)-z(i,j1,k)
              dd(2) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i1,j,k1)-x(i,j,k1)
              dy = y(i1,j,k1)-y(i,j,k1)
              dz = z(i1,j,k1)-z(i,j,k1)
              dd(3) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i1,j1,k1)-x(i,j1,k1)
              dy = y(i1,j1,k1)-y(i,j1,k1)
              dz = z(i1,j1,k1)-z(i,j1,k1)
              dd(4) = dsqrt(dx**2+dy**2+dz**2)
              ddi = 0.25d0*(dd(1)+dd(2)+dd(3)+dd(4))
c
              dx = x(i,j1,k)-x(i,j,k)
              dy = y(i,j1,k)-y(i,j,k)
              dz = z(i,j1,k)-z(i,j,k)
              dd(1) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i1,j1,k)-x(i1,j,k)
              dy = y(i1,j1,k)-y(i1,j,k)
              dz = z(i1,j1,k)-z(i1,j,k)
              dd(2) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i,j1,k1)-x(i,j,k1)
              dy = y(i,j1,k1)-y(i,j,k1)
              dz = z(i,j1,k1)-z(i,j,k1)
              dd(3) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i1,j1,k1)-x(i1,j,k1)
              dy = y(i1,j1,k1)-y(i1,j,k1)
              dz = z(i1,j1,k1)-z(i1,j,k1)
              dd(4) = dsqrt(dx**2+dy**2+dz**2)
              ddj = 0.25d0*(dd(1)+dd(2)+dd(3)+dd(4))
c
              dx = x(i,j,k1)-x(i,j,k)
              dy = y(i,j,k1)-y(i,j,k)
              dz = z(i,j,k1)-z(i,j,k)
              dd(1) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i1,j,k1)-x(i1,j,k)
              dy = y(i1,j,k1)-y(i1,j,k)
              dz = z(i1,j,k1)-z(i1,j,k)
              dd(2) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i,j1,k1)-x(i,j1,k)
              dy = y(i,j1,k1)-y(i,j1,k)
              dz = z(i,j1,k1)-z(i,j1,k)
              dd(3) = dsqrt(dx**2+dy**2+dz**2)

              dx = x(i1,j1,k1)-x(i1,j1,k)
              dy = y(i1,j1,k1)-y(i1,j1,k)
              dz = z(i1,j1,k1)-z(i1,j1,k)
              dd(4) = dsqrt(dx**2+dy**2+dz**2)
              ddk = 0.25d0*(dd(1)+dd(2)+dd(3)+dd(4))
c
              ddd = dmax1(ddi,ddj,ddk)
              dist(i,j,k) = dmin1(dist(i,j,k), cdes*ddd)
            end do
          end do
        end do
      end if
c
      return
      end subroutine wall_dist
