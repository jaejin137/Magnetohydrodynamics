        subroutine eleccond(dstep,ilower,iupper,jlower,jupper,klower,
     $    kupper,il,jl,kl,nl,x_new,y_new,z_new,r_new,u,v,w,e,t,
     $    t_thrsh, sigma_e_ref, sigma_e)

c       This subroutine defines the distribution of the electrical
c       conductivity.

        implicit none

c       Interface Variables

        integer, intent(in)::
     $  dstep,
     $  ilower, iupper,
     $  jlower, jupper,
     $  klower, kupper,
     $  il, jl, kl, nl

        DOUBLE PRECISION, INTENT(IN)::
     $  t_thrsh,
     $  sigma_e_ref
        
        double precision, intent(in), dimension(ilower:iupper,
     $  jlower:jupper,klower:kupper)::
     $  x_new, y_new, z_new,    ! coordinates of cell center
     $  r_new,                  ! radial distance of cell center
     $  u, v, w,                ! velocity components
     $  e,                      ! fluid mechanical energy
     $  t                       ! static temperature

c       Local Variables

        double precision, dimension(ilower:iupper,jlower:jupper,
     $  klower:kupper)::
     $  sigma_e                 ! electrical conductivity

        integer::
     $  i, j, k

        DOUBLE PRECISION, PARAMETER::
     $  PI=3.141592
c        double precision, parameter::
c     $  sigma_e_inf = 1.0d0    ! elelctrical conductivity at infinity

c       ******************* START SUBROUTINE *************************

c       Define the distribution of the electrical conductivity

        do k=klower,kupper
          do j=jlower,jupper
            do i=ilower,iupper
              if (t(i,j,k).ge.t_thrsh) then
                sigma_e = sigma_e_ref
              else
                sigma_e = 0.d0
              end if
c                sigma_e_star = sigma_e/sigma_e_inf
            end do
          end do
        end do

        return
        
        end subroutine
