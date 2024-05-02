      subroutine mvgrid_af(dim2, il, jl, kl, nl, x, y, z, xo, yo, zo,
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $     bc_zeta_lower, bc_zeta_upper,
     $     dxyz, xref, yref, jfx, beta, xref0, yref0)
c**********************************************************************
c     This subroutine is used to update the grid points at each time
c     step for forced pitching airfoil.
c     (xref0,yref0): The reference point of torque
c**********************************************************************
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     dim2,               ! maxium face number in 3 directions
     $     il, jl, kl,          ! cell number in 3 directions
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     nl

      double precision, intent(in):: 
     $     dxyz(3)            ! moving distance during each time step

      double precision, intent(inout),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z, xo, yo, zo

      integer, dimension(jl,kl), intent(in)::
     $     bc_xi_lower, bc_xi_upper

      integer, dimension(il,kl), intent(in)::
     $     bc_eta_lower, bc_eta_upper

      integer, dimension(il,jl), intent(in)::
     $     bc_zeta_lower, bc_zeta_upper
c
c     LOCAL VARIABLES
c
      double precision::
     $     radi, rado,         ! radii of inner & outer boundaries
     $     beta,               ! stretching parameter
     $     dsn(dim2), sn(dim2),! clustering output
     $     ddx, ddy, ck, cc, alf, xp, yp, 
     $     xref, yref, xref0, yref0
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     jfx,                 ! j index for non-deformed part
     $     ilp, jlp, klp        ! il+1, jl+1, kl+1

      integer, parameter::
     $     peri = 10
c
c *** SUBROUTINE START ***
c
c *** set up some parameters
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
      k = 1
c
c *** translating in vertical direction
c
      do j = 1,jlp
        do i = 1,ilp
          x(i,j,k) = xo(i,j,k)
          y(i,j,k) = yo(i,j,k) + dxyz(1)
        end do
      end do

      xref0 = xref
      yref0 = yref + dxyz(1)
c
c *** rotating around the reference point(xref0,yref0) of old mesh
c
      alf = dxyz(2)
      do j = 1,jlp
        do i = 1,ilp
          xp = x(i,j,k) - xref0
          yp = y(i,j,k) - yref0
          x(i,j,k) =  xp*dcos(alf) + yp*dsin(alf)
          y(i,j,k) = -xp*dsin(alf) + yp*dcos(alf)
        end do
      end do
c
c *** translating the airfoil back to the new position
c
      do j = 1,jlp
        do i = 1,ilp
          x(i,j,k) = x(i,j,k) + xref0
          y(i,j,k) = y(i,j,k) + yref0
        end do
      end do
c
c *** whole mesh
c
      do k=2,klp
        do j=1,jlp
          do i=1,ilp
            x(i,j,k)=x(i,j,1)
            y(i,j,k)=y(i,j,1)
          end do
        end do
      end do
c
      return
      end
