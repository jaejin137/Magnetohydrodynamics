      subroutine mvgrid_cyl(dim2, il, jl, kl, nl, x, y, z, xo, yo, zo,
     $          ilower, iupper, jlower, jupper, klower, kupper,
     $          dxyz, jfx, beta)
c
c     This subroutine is used to update the grid points at each time
c     step.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     dim2,               ! maxium face number in 3 directions
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! number of equations to solve
     $     ilower, iupper, jlower, jupper, klower, kupper
c
      double precision, intent(in):: 
     $     dxyz(3)             ! moving distance during each time step
c
      double precision, intent(inout),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z, xo, yo, zo
c
c     LOCAL VARIABLES
c
      double precision::
     $     radi, rado,         ! radii of inner & outer boundaries
     $     beta,               ! stretching parameter
     $     eta, deta,          ! eta & interval for eta
     $     chh, at, pup, cup, cdw, ddx, ddy, ymv, alf, ymvp,
     $     dyg(dim2), dymv
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     jfx,                 ! j index for non-deformed part
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
c *** moving grids
c
c ... update the position for non-deformed part around the cylinder
c
      do 10 k = 1,klp
      do 10 j = 1,jlp
      do 10 i = 1,ilp
        x(i,j,k) = xo(i,j,k) + dxyz(1)
        y(i,j,k) = yo(i,j,k) + dxyz(2)
   10 continue
c
      return
      end
