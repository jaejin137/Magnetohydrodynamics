      subroutine volume(x, y, z, vol, ilower,
     $     iupper, jlower, jupper, klower, kupper, il, jl, kl,
     $     bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper)
c
c     compute volume
c
c     IMPLICIT STATEMENT
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl
c      
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x,  y, z
c
      integer, intent(in)::
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl)

      double precision, intent(out), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol

c     
c     LOCAL VARIABLES
c
      double precision, dimension(:), allocatable::
     $     lx, ly, lz,
     $     mx, my, mz,
     $     nx, ny, nz,
     $     dx1, dx2,
     $     dy1, dy2,
     $     dz1, dz2             ! distance intervals
c
      integer::
     $     i, j, k,             ! iteration index
     $     jp, kp, dim2,         ! j+1, k+1
     $     bl

      integer:: invalidvol

      integer, parameter::
     $     in_bnd = 7,
     $     per_bnd = 10,
     $     sym_bnd = 8
      
c
c *** SUBROUTINE START ***
c      
      invalidvol = 0
c
      bl=1-ilower
c
      dim2 = max(iupper, jupper, kupper)
      allocate(lx(ilower+1:dim2), ly(ilower+1:dim2), lz(ilower+1:dim2),
     $     mx(ilower+1:dim2), my(ilower+1:dim2), mz(ilower+1:dim2),
     $     nx(ilower+1:dim2), ny(ilower+1:dim2), nz(ilower+1:dim2),
     $     dx1(ilower+1:dim2), dx2(ilower+1:dim2),
     $     dy1(ilower+1:dim2), dy2(ilower+1:dim2),
     $     dz1(ilower+1:dim2), dz2(ilower+1:dim2))
      

      do k = klower+1,kupper-1
         do j = jlower+1,jupper-1
            jp = j + 1
            kp = k + 1

c     define lx, ly, lz

            do i = ilower+1,iupper-1
               dx1(i) = x(i,jp,kp) - x(i,j,k)
               dy1(i) = y(i,jp,kp) - y(i,j,k)
               dz1(i) = z(i,jp,kp) - z(i,j,k)
               dx2(i) = x(i,j,kp)  - x(i,jp,k)
               dy2(i) = y(i,j,kp)  - y(i,jp,k)
               dz2(i) = z(i,j,kp)  - z(i,jp,k)

               lx(i) = 0.5d0 * (dy1(i) * dz2(i) - dz1(i) * dy2(i))
               ly(i) = 0.5d0 * (dz1(i) * dx2(i) - dx1(i) * dz2(i))
               lz(i) = 0.5d0 * (dx1(i) * dy2(i) - dy1(i) * dx2(i))
            end do
c
c     define mx, my, mz for cell i
c     (note that the vectors mx(i), ... are oriented in the i-direction
c     instead of j-direction as done in subroutine metric)
c
            do i = ilower+1,iupper-1
               dx1(i) = x(i+1,j,kp) - x(i,j,k)
               dy1(i) = y(i+1,j,kp) - y(i,j,k)
               dz1(i) = z(i+1,j,kp) - z(i,j,k)
               dx2(i) = x(i+1,j,k)  - x(i,j,kp)
               dy2(i) = y(i+1,j,k)  - y(i,j,kp)
               dz2(i) = z(i+1,j,k)  - z(i,j,kp)

               mx(i) = 0.5d0 * (dy1(i) * dz2(i) - dz1(i) * dy2(i))
               my(i) = 0.5d0 * (dz1(i) * dx2(i) - dx1(i) * dz2(i))
               mz(i) = 0.5d0 * (dx1(i) * dy2(i) - dy1(i) * dx2(i))
            end do
c
c     define nx, ny, nz for cell i
c     (note that the vectors nx(i), ... are oriented in the i-direction
c     instead of k-direction as done in subroutine metric)
c
            do i = ilower+1, iupper-1
               dx1(i) = x(i+1,jp,k) - x(i,j,k)
               dy1(i) = y(i+1,jp,k) - y(i,j,k)
               dz1(i) = z(i+1,jp,k) - z(i,j,k)
               dx2(i) = x(i,jp,k)  - x(i+1,j,k)
               dy2(i) = y(i,jp,k)  - y(i+1,j,k)
               dz2(i) = z(i,jp,k)  - z(i+1,j,k)

               nx(i) = 0.5d0 * (dy1(i) * dz2(i) - dz1(i) * dy2(i))
               ny(i) = 0.5d0 * (dz1(i) * dx2(i) - dx1(i) * dz2(i))
               nz(i) = 0.5d0 * (dx1(i) * dy2(i) - dy1(i) * dx2(i))
            end do
c
c     define volume
c
            do i = ilower+1, iupper-1
               dx1(i) = x(i+1,jp,kp) - x(i,j,k)
               dy1(i) = y(i+1,jp,kp) - y(i,j,k)
               dz1(i) = z(i+1,jp,kp) - z(i,j,k)
               vol(i,j,k) = 0.33333333d+0 * 
     >              ( (lx(i) + mx(i) + nx(i)) * dx1(i) +
     >              (ly(i) + my(i) + ny(i)) * dy1(i) +
     >              (lz(i) + mz(i) + nz(i)) * dz1(i) )
            end do

         end do
      end do

      do i = 1, il
         do j = 1, jl
            do k = 1, kl
               if(vol(i,j,k).le.0.d0) then
                  invalidvol = 1
                  write(*, "('vol =< 0 at: ',3(i5,1x))") i, j, k
                  print *, i,j,k,':'
                  print *, x(i,j,k), y(i,j,k), z(i,j,k)
                  print *, i+1,j,k,':'
                  print *, x(i+1,j,k), y(i+1,j,k), z(i+1,j,k)
                  print *, i,j+1,k,':'
                  print *, x(i,j+1,k), y(i,j+1,k), z(i,j+1,k)
                  print *, i+1,j+1,k,':'
                  print *, x(i+1,j+1,k), y(i+1,j+1,k), z(i+1,j+1,k)
                  print *, i,j,k+1,':'
                  print *, x(i,j,k+1), y(i,j,k+1), z(i,j,k+1)
                  print *, i+1,j,k+1,':'
                  print *, x(i+1,j,k+1), y(i+1,j,k+1), z(i+1,j,k+1)
                  print *, i,j+1,k+1,':'
                  print *, x(i,j+1,k+1), y(i,j+1,k+1), z(i,j+1,k+1)
                  print *, i+1,j+1,k+1,':'
                  print *, x(i+1,j+1,k+1), y(i+1,j+1,k+1),z(i+1,j+1,k+1)
                  stop
               end if
            end do
         end do
      end do
c---------------------------
      do j = 1, jl
         do k = 1, kl
            if (bc_xi_lower(j,k).ne.in_bnd) then
               if (bc_xi_lower(j,k).eq.per_bnd) then
                  do i=ilower+1,0
                     vol(i,j,k) = vol(i+il,j,k)
                  end do
               else
                  do i=ilower+1,0
                     vol(i,j,k) = vol(1-i,j,k)
                  end do
               end if
            end if
         end do
      end do
      do j = 1, jl
         do k = 1, kl
            if (bc_xi_upper(j,k).ne.in_bnd) then
               if (bc_xi_upper(j,k).eq.per_bnd) then
                  do i = il+1, iupper-1
                    vol(i,j,k) = vol(i-il,j,k)
                  end do
               else
                  do i = il+1, iupper-1
                     vol(i,j,k) =vol(2*il-i+1,j,k)
                  end do
               end if
            end if
         end do
      end do
c--------      
      do i = 1, il
         do k = 1, kl
            if (bc_eta_lower(i,k).ne.in_bnd) then
               if (bc_eta_lower(i,k).eq.per_bnd) then
                  do j = jlower+1, 0
                     vol(i,j,k) = vol(i,j+jl,k)
                  end do
               else
                  do j = jlower+1, 0
                     vol(i,j,k) =vol(i,1-j,k) 
                  end do
               end if
            end if
         end do
      end do
      do i = 1, il
         do k = 1, kl
            if (bc_eta_upper(i,k).ne.in_bnd) then
               if (bc_eta_upper(i,k).eq.per_bnd) then
                  do j = jl+1, jupper-1
                     vol(i,j,k) = vol(i,j-jl,k)
                  end do
               else
                  do j = jl+1, jupper-1
                    vol(i,j,k) =vol(i,2*jl-j+1,k)
                  end do
               end if
            end if
         end do
      end do
c--------     
      do i = 1, il
         do j = 1, jl
            if (bc_zeta_lower(i,j).ne.in_bnd) then
               if (bc_zeta_lower(i,j).eq.per_bnd) then
                   do k = klower+1, 0
                      vol(i,j,k) = vol(i,j,k+kl)
                   end do
               else
                   do k = klower+1, 0
                      vol(i,j,k) =vol(i,j,1-k) 
                   end do
               end if
            end if
         end do
      end do
      do i = 1, il
         do j = 1, jl
            if (bc_zeta_upper(i,j).ne.in_bnd ) then
               if (bc_zeta_upper(i,j).eq.per_bnd) then
                  do k = kl+1, kupper-1
                     vol(i,j,k) = vol(i,j,k-kl)
                  end do
               else
                  do k = kl+1, kupper-1
                     vol(i,j,k) =vol(i,j,2*kl-k+1) 
                  end do
               end if
            end if
         end do
      end do
c----------------------------------------------------    

      deallocate(lx, ly, lz, mx, my, mz, nx, ny, nz,
     $     dx1, dx2, dy1, dy2, dz1, dz2)

      if(invalidvol.gt.0) then
         write(*,*) 'invalid cell exists, program stop'
         stop
      end if
      
      end
