        subroutine magfld(dstep,ilower,iupper,jlower,jupper,klower,
     $    kupper,il,jl,kl,nl,x_new,y_new,z_new,r_new,u,v,w,e,L_ref,
     $    mhdref,B_ref,I_in,B_x,B_y,B_z)

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

        character,intent(in)::
     $  mhdref
        
        double precision, intent(in)::
     $  L_ref,
     $  B_ref,
     $  I_in                    ! current through the conductor

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

        double precision::
     $  d,                      ! thickness of flatplate
     $  R_ref,                  ! Radius of conductor
     $  I_in_tmp

        double precision, intent(out), dimension(ilower:iupper,
     $    jlower:jupper,klower:kupper)::
     $  B_x, B_y, B_z           ! magnetic field

        double precision, parameter::
     $  pi = 3.14159265359,
     $  mu_e0 = pi*4.0d-07     ! permeability of vacuum


c       ********************** Start Subroutine ************************

        d = L_ref/10
        R_ref = d/4

        if (dstep==1) then
          B_x = 0.0d0
          B_y = 0.0d0
          B_z = 0.0d0
        end if

c       Calculate the magnetic field at each points
        if (mhdref.eq.'B') then ! Assign B-field based on B_ref
          I_in_tmp = (2*pi*R_ref/mu_e0)*B_ref
          do k=klower,kupper
            do j=jlower,jupper
              do i=ilower,iupper
                if (r_new(i,j,k).le.R_ref) then
                  B_x(i,j,k) = -mu_e0*I_in_tmp*y_new(i,j,k)/
     $              (2*pi*R_ref**2)
                  B_y(i,j,k) = mu_e0*I_in_tmp*x_new(i,j,k)/
     $              (2*pi*R_ref**2)
                  B_z(i,j,k) = 0.0d0
                else
                  B_x(i,j,k) = -mu_e0*I_in_tmp*y_new(i,j,k)/
     $              (2*pi*r_new(i,j,k)**2)
                  B_y(i,j,k) = mu_e0*I_in_tmp*x_new(i,j,k)/
     $              (2*pi*r_new(i,j,k)**2)
                  B_z(i,j,k) = 0.0d0
                end if
              end do
            end do
          end do
        else if (mhdref.eq.'I') then    ! Based on input current
          do k=klower,kupper
            do j=jlower,jupper
              do i=ilower,iupper
                if (r_new(i,j,k).le.R_ref) then
                    B_x(i,j,k) = -mu_e0*I_in*y_new(i,j,k)/
     $                (2*pi*R_ref**2)
                    B_y(i,j,k) = mu_e0*I_in*x_new(i,j,k)/
     $                (2*pi*R_ref**2)
                    B_z(i,j,k) = 0.0d0
                  else
                    B_x(i,j,k) = -mu_e0*I_in*y_new(i,j,k)/
     $                (2*pi*r_new(i,j,k)**2)
                    B_y(i,j,k) = mu_e0*I_in*x_new(i,j,k)/
     $                (2*pi*r_new(i,j,k)**2)
                    B_z(i,j,k) = 0.0d0
                end if
              end do
            end do
          end do
        end if

        return
        
        end subroutine magfld
