        subroutine magfld(dstep,ilower,iupper,jlower,jupper,klower,
     $    kupper,il,jl,kl,nl,x_new,y_new,z_new,r_new,u,v,w,e,
     $    L_ref,mhdref,B_ref,I_in,B_x,B_y,B_z)

c       This subroutine defines the magnetic field due to infinitely
c       long line current source.

        implicit none

c       *********************** VARIABLES ***************************

c       Interface Variables

        integer, intent(in)::
     $  dstep,
     $  ilower, iupper,
     $  jlower, jupper,
     $  klower, kupper,
     $  il, jl, kl, nl
        
        double precision, intent(in)::
     $  L_ref,                  ! characteristic length
     $  B_ref,                  ! reference value for magnetic induction
     $  I_in                    ! current through the conductor

        character, intent(in)::
     $  mhdref              ! selection of reference

        double precision, intent(in), dimension(ilower:iupper,
     $    jlower:jupper,klower:kupper)::
     $  u, v, w,                ! velocity in each coordinates
     $  e,                      ! fluid dynamical energy density
     $  x_new, y_new, z_new,    ! coordinates of cell center
     $  r_new                   ! radial distance of cell center

c        double precision, intent(out), dimension(ilower:iupper,
c     $    jlower:jupper,klower:kupper)::
c     $  B_x_star, B_y_star, B_z_star    ! nondimensionalized Bs
        
c       Local Variables

        integer::
     $  i, j, k

        double precision, intent(out), dimension(ilower:iupper,
     $    jlower:jupper,klower:kupper)::
     $  B_x, B_y, B_z           ! magnetic field

        double precision::
     $  d, R_ref

        double precision, parameter::
     $  pi = 3.14159265359,
     $  mu_e0 = pi*4.0d-07     ! permeability of vacuum


c       ********************** Start Subroutine ************************

        d = L_ref/10           ! thickness of flat plate
        R_ref = d/4            ! radius of conductor

c       Calculate the magnetic field at each points
        if (mhdref .eq. 'B') then ! Assign B-field based on B_ref
            do k=klower,kupper
              do j=jlower,jupper
                do i=ilower,iupper
                  B_x(i,j,k) = 0.0d0
                  B_y(i,j,k) = -B_ref
                  B_Z(i,j,k) = 0.0d0
                end do
              end do
            end do
          elseif (mhdref .eq. 'I') then     ! Based on input current
            do k=klower,kupper
              do j=jlower,jupper
                do i=ilower,iupper
                  B_x(i,j,k) = 0.0d0
                  B_y(i,j,k) = -mu_e0*I_in/L_ref
                  B_z(i,j,k) = 0.0d0
                end do
              end do
            end do
          else
            write(*,*) "The magnetic induction can't be calculated."
        end if
            
        return
        
        end subroutine magfld
