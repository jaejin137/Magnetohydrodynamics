        subroutine source_mhd_lowrey(rank, dstep, intersteps,
     $  ilower, iupper, jlower, jupper, klower, kupper, il, jl, kl,
     $  nl, x_grd, y_grd, z_grd, q, gamma, mach_inf, rho_inf, T_inf,
     $  c_sound, L_ref, ionmodel, T_onset, sigma_e_ref, magtype, magref,
     $  B_ref, I_in, offset, sigma_e_star, B_x_star, B_y_star,
     $  B_z_star, B_norm, src_mhd)

c        subroutine source_mhd_lowrey(dstep, ilower, iupper, jlower,
c     $  jupper, klower, kupper, il, jl, kl, nl, bc_xi_lower,
c     $  bc_xi_upper, bc_eta_lower, bc_eta_upper, bc_zeta_lower,
c     $  bc_zeta_upper, x, y, z, q, vol, mach_inf, src_mhd)

c       This subroutine is to compute MHD source and to return the result
c       as a seperate source term from the fluid dynamical sources.

        use newtypes

        IMPLICIT NONE

c       *********************** VARIABLES ****************************

c       Interface Variables

        character, intent(in)::
     $  ionmodel,
     $  magtype*10,
     $  magref

        integer, intent(in)::
     $  rank,
     $  dstep,
     $  intersteps,
     $  ilower, iupper,
     $  jlower, jupper,
     $  klower, kupper,
     $  il, jl, kl, nl

c        integer, intent(in)::
c     $  bc_xi_lower(jl,kl), bc_xi_upper(jl,kl),
c     $  bc_eta_lower(il,kl), bc_eta_upper(il,kl),
c     $  bc_zeta_lower(il,jl), bc_zeta_upper(il,jl)

        double precision, intent(in)::
     $  gamma,         ! ratio of specific heats
     $  mach_inf,      ! Mach number at infinity
     $  rho_inf,       ! density at infinity
     $  T_inf,         ! temperature at infinity
     $  c_sound,       ! speed of sound
     $  L_ref,         ! length of the system
     $  T_onset,       ! threshold temperature for onset of nonzero
     $  sigma_e_ref,   ! reference value for sigma_e
     $  B_ref,         ! reference value for magnetic induction
     $  I_in,          ! reference value for input current
     $  offset(3)      ! offsets for new origin for B
        
        double precision, intent(in), dimension(ilower+1:iupper,
     $    jlower+1:jupper,klower+1:kupper)::
     $  x_grd, y_grd, z_grd  ! coordinates of grid points

        double precision, intent(in)::
     $  q(ilower:iupper,jlower:jupper,klower:kupper,nl)

        double precision,  dimension(ilower:iupper,
     $    jlower:jupper,klower:kupper)::
     $  sigma_e,                ! electrical conductivity
     $  B_x, B_y, B_z           ! magnetic field

        double precision, intent(out), dimension(ilower:iupper,
     $    jlower:jupper,klower:kupper)::
     $  sigma_e_star,                ! nondimensional sigma_e
     $  B_x_star, B_y_star, B_z_star,! nondimensional magnetic field
     $  B_norm                       ! norm of vector B_star

        double precision, intent(out)::
     $  src_mhd(il,jl,kl,2:5)


c       Local Variables

        integer::
     $  i, j, k, n,
     $  i_lbnd, j_lbnd, k_lbnd, ! indices of the lower boundary
     $  i_ubnd, j_ubnd, k_ubnd  ! indices of the upper boundary   
!     $  mhdsim

!        character::
!     $  fdatain_mhd*6

c        double precision, dimension(il+1,jl+1,kl+1) ::
c     $  lx, ly, lz, mx, my, mz, nx, ny, nz

        double precision::
     $  src_mhd_ext(ilower+1:iupper,jlower+1:jupper,klower+1:kupper,2:5)

        double precision, dimension(ilower:iupper,jlower:jupper,
     $  klower:kupper)::
     $  x_cnt, y_cnt, z_cnt,    ! coordinates of cell center
     $  x_new, y_new, z_new,    ! new coordinates translated from the
                                ! original origin to calculate B-field
     $  r_cnt,                  ! radial distance of cell center
     $  r_new,                  ! radial distance of cell center from
                                ! new origin
     $  delx, dely, delz,       ! cell size in each coordinates
     $  rho,                    ! density
     $  u, v, w,                ! velocitiy
     $  e,                      ! fluid mechanical energy
     $  t,                      ! static temperature
     $  J_x_star, J_y_star, J_z_star, ! nondimensionalized current density
     $  J_norm,                 ! norm of J at each points
     $  f_mhd_x, f_mhd_y, f_mhd_z, ! Lorentz force density
     $  f_mhd_norm,             ! norm of f_mhd at each points
     $  e_mhd                   ! magnetic energy density

        double precision::
     $  x_offset,       ! offset of the origin of B in x direction
     $  y_offset,       ! offset of the origin of B in y direction
     $  z_offset,       ! offset of the origin of B in z direction
     $  U_inf,          ! speed at infinity
     $  t_thrsh,        ! nondimensional threshold temperature for MHD
     $  sigma_e_inf,    ! electrical conductivity at infinity
     $  Re_m_inf,       ! magnetic Reynolds number
     $  sigma_e_max,    ! maximum value of elec. conductivity
     $  B_max           ! maximum value of magnetic field

        character::
     $  foutput_mhd*10,
     $  journal_caller*20

!        integer, parameter::
!     $  fdatain = 10

        type(journaldata)::
     $    journal_info
        
        double precision, parameter::
     $  pi = 3.14159265359,
     $  mu_e0 = pi*4.0d-7,      ! magnetic permeability of vacuum
     $  mu_e0_star = 1.0d0

c       Initialization
        if (dstep==1) then
          sigma_e_star = 0.0d0
          B_x_star = 0.0d0
          B_y_star = 0.0d0
          B_z_star = 0.0d0
          B_norm = 0.0d0
          src_mhd = 0.0d0
          x_offset = offset(1)
          y_offset = offset(2)
          z_offset = offset(3)
        end if

        t_thrsh = T_onset/T_inf        
        sigma_e_inf = sigma_e_ref
          
c       Set dummy src_mhd components
c
c       Assign 0 to 1st component of src_mhd for all time steps
c        src_mhd(:,:,:,1) = 0.0d0
        
c       Assign 0 to 6th and 7th components of src_mhd for turb=1 case      
c        if(nl>5) then
c          do n=6,nl
c            src_mhd(:,:,:,n) = 0.0d0
c          end do
c        end if

c       Define cell centers and cell sizes in terms of grid points
c       except boundaries
        do k=klower+1,kupper-1
          do j=jlower+1,jupper-1
            do i=ilower+1,iupper-1
              delx(i,j,k) = x_grd(i+1,j,k) -x_grd(i,j,k)
              dely(i,j,k) = y_grd(i,j+1,k) -y_grd(i,j,k)
              delz(i,j,k) = z_grd(i,j,k+1) -z_grd(i,j,k)
              x_cnt(i,j,k) = (x_grd(i,j,k) +x_grd(i+1,j,k))/2
              x_new(i,j,k) = x_cnt(i,j,k) -x_offset
              y_cnt(i,j,k) = (y_grd(i,j,k) +y_grd(i,j+1,k))/2
              y_new(i,j,k) = y_cnt(i,j,k) -y_offset
              z_cnt(i,j,k) = (z_grd(i,j,k) +z_grd(i,j,k+1))/2
              z_new(i,j,k) = z_cnt(i,j,k) -z_offset
c              write(*,*) 'x_cnt=',x_cnt(i,j,k),'y_cnt=',y_cnt(i,j,k)
            end do
          end do
        end do
        
c       Define cell centers at boundaries in terms of grid points
        do k=klower+1,kupper-1
          do j=jlower+1,jupper-1
            x_cnt(ilower,j,k) = x_cnt(ilower+1,j,k) -delx(ilower+1,j,k)
            x_new(ilower,j,k) = x_cnt(ilower,j,k) -x_offset
            x_cnt(iupper,j,k) = x_cnt(iupper-1,j,k) +dely(iupper-1,j,k)
            x_new(iupper,j,k) = x_cnt(iupper,j,k) -x_offset
          end do
        end do
        
        do k=ilower+1,kupper-1
          do i=ilower+1,iupper-1        
            y_cnt(i,jlower,k) = y_cnt(i,jlower+1,k) -dely(i,jlower+1,k)
            y_new(i,jlower,k) = y_cnt(i,jlower,k) -y_offset
            y_cnt(i,jupper,k) = y_cnt(i,jupper-1,k) +dely(i,jupper-1,k)
            y_new(i,jupper,k) = y_cnt(i,jupper,k) -y_offset
          end do
        end do
        
        do j=jlower+1,jupper-1
          do i=ilower+1,iupper-1
            z_cnt(i,j,klower) = z_cnt(i,j,klower+1) -delz(i,j,klower+1)
            z_new(i,j,klower) = z_cnt(i,j,klower) -z_offset
            z_cnt(i,j,kupper) = z_cnt(i,j,kupper-1) +delz(i,j,kupper-1)
            z_new(i,j,kupper) = z_cnt(i,j,kupper) -z_offset
          end do
        end do

c       Define radial distances of each cell centers
        r_cnt = sqrt(x_cnt**2 +y_cnt**2 +z_cnt**2)
        r_new = sqrt(x_new**2 +y_new**2 +z_new**2)
        
c       Rename qs with fluid mechanical variables
        rho = q(:,:,:,1)
        u = q(:,:,:,2)/rho
        v = q(:,:,:,3)/rho
        w = q(:,:,:,4)/rho
        e = q(:,:,:,5)/rho

c       Calculate the static temperature
        t = gamma*(gamma-1.0d0)*mach_inf**2*(e-.5d0*(u**2+v**2+w**2))

c       Calculate magnetic Reynolds number
        U_inf = mach_inf*c_sound
        Re_m_inf = mu_e0*sigma_e_ref*U_inf*L_ref
c       Check if the electrical conductivity changes in time
c
c       if (true) then
c         update the electrical conductivity
        call eleccond(dstep, ilower, iupper, jlower, jupper, klower,
     $    kupper, il, jl, kl, nl, x_new, y_new, z_new, r_new, u, v, w,
     $    e, t, ionmodel, t_thrsh, sigma_e_ref, sigma_e)
c       else
c         recycle the previous electrical conductivity
c       end if

c       Nondimensionalize the electrical conductivity
        sigma_e_star = sigma_e/sigma_e_inf

        sigma_e_max = maxval(sigma_e)

c       Check if the magnetic field changes in time
c
c       if (true) then
c         update the magnetic field
        call magfld(dstep, ilower, iupper, jlower, jupper, klower,
     $    kupper, il, jl, kl, nl, x_new, y_new, z_new, r_new,
     $    u, v, w, e, L_ref, magtype, magref, B_ref, I_in,
     $    B_x, B_y, B_z)
c       else
c         recycle the previous magnetic field
c       end if

c       Nondimensionalize the magnetic field and calculate norm and max
        B_x_star = B_x/(U_inf*sqrt(mu_e0*rho_inf))
        B_y_star = B_y/(U_inf*sqrt(mu_e0*rho_inf))
        B_z_star = B_z/(U_inf*sqrt(mu_e0*rho_inf))
        B_norm = sqrt(B_x_star**2 +B_y_star**2 +B_z_star**2)
c        B_max = max(maxval(B_x),maxval(B_y),maxval(B_z))
        B_max = maxval(B_norm)*(U_inf*sqrt(mu_e0*rho_inf))

c       Calculate the current density and its norm
        J_x_star = sigma_e_star*(v*B_z_star -w*B_y_star)
        J_y_star = sigma_e_star*(w*B_x_star -u*B_z_star)
        J_z_star = sigma_e_star*(u*B_y_star -v*B_x_star)
        J_norm = sqrt(J_x_star**2 +J_y_star**2 +J_z_star**2)

c       Calculate the Lorentz force density and its norm
        f_mhd_x = J_y_star*B_z_star -J_z_star*B_y_star
        f_mhd_y = J_z_star*B_x_star -J_x_star*B_z_star
        f_mhd_z = J_x_star*B_y_star -J_y_star*B_x_star
        f_mhd_norm = sqrt(f_mhd_x**2 +f_mhd_y**2 +f_mhd_z**2)

c       Calculate the magnetic field energy

        select case(ionmodel)

        case('F')           ! Frozen ionization model
          e_mhd = (f_mhd_x*u +f_mhd_y*v +f_mhd_z*w +1/sigma_e_star*
     $      (J_x_star**2 +J_y_star**2 +J_z_star**2))

        case('E')           ! Equilibrium ionization model
          do k=1,kl
            do j=1,jl
              do i=1,il
                if (sigma_e_star(i,j,k).gt.1.d-10) then
                  e_mhd(i,j,k) = (f_mhd_x(i,j,k)*u(i,j,k)
     $              +f_mhd_y(i,j,k)*v(i,j,k) +f_mhd_z(i,j,k)*w(i,j,k)
     $              +1/sigma_e_star(i,j,k)*(J_x_star(i,j,k)**2
     $              +J_y_star(i,j,k)**2 +J_z_star(i,j,k)**2))
                else
                  e_mhd(i,j,k) = 0.d0
                end if
              end do
            end do
          end do

        end select

c       Assign computed values of f_mhd and e_mhd to src_mhd and write
c       mhd variables in mhd.dat for reference

        do k=1,kl
          do j=1,jl
            do i=1,il
              src_mhd(i,j,k,2) = Re_m_inf*f_mhd_x(i,j,k)
              src_mhd(i,j,k,3) = Re_m_inf*f_mhd_y(i,j,k)
              src_mhd(i,j,k,4) = Re_m_inf*f_mhd_z(i,j,k)
              src_mhd(i,j,k,5) = Re_m_inf*e_mhd(i,j,k)
            end do
          end do
        end do

c       Assign dummy values of the mhd sources to the ghost layers
c       beyond the actual lower and upper boundary for consistency in
c       dimension with coordinate variables when writing mhd output file.

        i_lbnd = 1
        j_lbnd = 1
        k_lbnd = 1
        i_ubnd = il
        j_ubnd = jl
        k_ubnd = kl

        do k=1,kl
          do j=1,jl
            do i=1,il
              src_mhd_ext(i,j,k,2)= src_mhd(i,j,k,2)
              src_mhd_ext(i,j,k,3)= src_mhd(i,j,k,3)
              src_mhd_ext(i,j,k,4)= src_mhd(i,j,k,4)
              src_mhd_ext(i,j,k,5)= src_mhd(i,j,k,5)
            end do
          end do
        end do

        do k=klower+1,kl-1
          do j=jlower+1,jl-1
            do i=ilower+1,il-1
              src_mhd_ext(i,j,k,2)= src_mhd(i_lbnd,j_lbnd,k_lbnd,2)
              src_mhd_ext(i,j,k,3)= src_mhd(i_lbnd,j_lbnd,k_lbnd,3)
              src_mhd_ext(i,j,k,4)= src_mhd(i_lbnd,j_lbnd,k_lbnd,4)
              src_mhd_ext(i,j,k,5)= src_mhd(i_lbnd,j_lbnd,k_lbnd,5)
            end do
          end do
        end do

        do k=kl+1,kupper
          do j=jl+1,jupper
            do i=il+1,iupper
              src_mhd_ext(i,j,k,2)= src_mhd(i_ubnd,j_ubnd,k_ubnd,2)
              src_mhd_ext(i,j,k,3)= src_mhd(i_ubnd,j_ubnd,k_ubnd,3)
              src_mhd_ext(i,j,k,4)= src_mhd(i_ubnd,j_ubnd,k_ubnd,4)
              src_mhd_ext(i,j,k,5)= src_mhd(i_ubnd,j_ubnd,k_ubnd,5)
            end do
          end do
        end do

        if (mod(dstep,intersteps).eq.0) then
          write(foutput_mhd,'("mhd",i3.3,".dat")') rank+1
          open(50,file=foutput_mhd,form='formatted',status='replace')

500       format('TITLE = "MHD Flow"',/,
     $ 'VARIABLES = x, y, z, u, v, w, sigma_e, B_x, B_y, B_z,
     $ f_mhd_x, f_mhd_y, f_mhd_z, e_mhd',/,
     $ 'ZONE T = MHD, I = ',i5,', J = ',i5,', K = ',i5,',
     $ F = Point')
        
510       format(1x,3f16.9,4x,3f16.9,4x,f16.9,4x,3f16.9,4x,3f16.9,4x,
     $ f16.9)

c        Update the record of MHD quantities

          write(50,500) il+2, jl+2, kl+2
          do k=0,kl+1
            do j=0,jl+1
              do i=0,il+1
                write(50,510) x_grd(i,j,k), y_grd(i,j,k), z_grd(i,j,k),
     $ u(i,j,k), v(i,j,k), w(i,j,k), sigma_e_star(i,j,k),
     $ B_x_star(i,j,k), B_y_star(i,j,k), B_z_star(i,j,k),
     $ src_mhd_ext(i,j,k,2), src_mhd_ext(i,j,k,3), src_mhd_ext(i,j,k,4),
     $ src_mhd_ext(i,j,k,5)
              end do
            end do
          end do          

        end if
      
        close(50)

c       Journaling

        if (rank.eq.0) then
          if (dstep.eq.1) then
            journal_caller = 'source_mhd_lowrey'
c            journal_info%tag = ''
c            journal_info%char1 = ''
c            journal_info%char2 = ''
c            journal_info%intg1 = 0
c            journal_info%intg2 = 0
            journal_info%doub1 = sigma_e_inf
            journal_info%doub2 = sigma_e_max
            journal_info%doub3 = B_max
            journal_info%doub4 = Re_m_inf
            journal_info%doub5 = 0.0d0
            journal_info%doub6 = 0.0d0
            call journal(journal_caller,journal_info)
          end if
        end if

        return
        
        end subroutine source_mhd_lowrey
