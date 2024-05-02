      subroutine moving_grid(il, jl, kl, nl, ilower, iupper, jlower,
     $     jupper, klower, kupper, dim2, dstep, nn,  
     $     np, rank, snum, 
     $     xsl, xzl, time, 
     $     bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $     bc_zeta_lower, bc_zeta_upper, xs, ys, x, y, z, xo, yo, zo,
     $     qt, q, vol, vbd, xyz, udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $     wdk, kc, cs, mus, ub, wh, walf, uinf, nfp, alf0,
     $     nmode, hx, hy, hz, wpw, bcwake, cv1, lref, aref,
     $     xctr, yctr, xref, yref, xtor, ytor, jfx, beta, xref0, yref0,
     $     control)
c
c     update grid related info
c
      use datain

      implicit none

      type (datain_type), intent(in)::control
c
c     comment 1:
c     in parallel cascade case, some of the blades are oscillating,
c     others are not, grid updation is no only enabled when moving>0,
c     'strtp' is used to help determine the necessity of computing 
c     and exchanging mesh info (mpi). strtp=0 is also used to make
c     sure moving_gird is ignored. this is set in initial.f
c
c     Interface Variables
c
      integer, intent(in)::
     $     il, jl, kl,          ! computation domain dimension
     $     nl,                  ! equation number
     $     ilower, iupper,      ! extended comp domain dimen
     $     jlower, jupper,
     $     klower, kupper,
     $     dim2,                ! maxium face number in 3 directions
     $     dstep,               ! physical time marching steps
     $     nn,                  ! pseudo time step index
     $     np,            
     $     rank,            
     $     snum,                ! osci mesh number
     $     nfp,
     $     nmode
      integer::
     $     moving,              ! moving grid switch
     $     strtp,               ! moving grid struture type
     $     turb                 ! turbulence model switch
c
      double precision, intent(in)::
     $     xzl(2,3,6),
     $     time, cv1,           ! physical time
     $     cs, mus, ub,         ! parameters for induced vibration of cylinder
     $     kc, wh, walf,        ! parameters for induced vibration of airfoil
     $     alf0, uinf           ! parameters for induced vibration of airfoil
      double precision::
     $     gamma,
     $     machinf,
     $     reynolds,
     $     tref,
     $     tintvl,              ! physical time interval
     $     ronum                ! rossby number
c
      logical::
     $     suther,              ! viscosity method switch
     $     inviscid
c
      integer, intent(in)::     ! boundary conditions
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl)
      integer, intent(in):: bcwake(3,2)
c
      double precision, intent(in)::
     $     xs(il+1,jl+1,snum),     ! oscillating mesh for strtp=3
     $     ys(il+1,jl+1,snum)
c
      double precision, intent(in)::
     $     wpw(6),
     $     hx(il+1,kl+1,6), hy(il+1,kl+1,6), hz(il+1,kl+1,6)
c
      double precision, intent(inout),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z,             ! coords of current step
     $     xo, yo, zo           ! coords of previous step
c
      double precision, intent(out)::
     $     qt(ilower:iupper,jlower:jupper,klower:kupper,3)
c
      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(inout):: ! mesh cell volume
     $     vol(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)
c
      double precision, intent(inout)::
     $     xsl(4,3),
     $     vbd(3,3),
     $     xyz(3,3)
c
      double precision, dimension(jlower:jupper,
     $     klower:kupper,4), intent(inout):: udi, vdi, wdi
c
      double precision, dimension(ilower:iupper,
     $     klower:kupper,4), intent(inout):: udj, vdj, wdj
c
      double precision, dimension(ilower:iupper,
     $     jlower:jupper,4), intent(inout):: udk, vdk, wdk
c
c     Local Variables
c
      double precision::
     $     lref, aref, cl, cd, ct
c
      double precision::
     $     pj(6),
     $     ux(il+1,kl+1), uy(il+1,kl+1), uz(il+1,kl+1)
c
      double precision::
     $  visturb(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)
c
      double precision::
     $     dxyz(3)              ! moving distance of cylinder per each time step
c
      double precision:: xctr, yctr, xtor, ytor, beta
      double precision:: xref, yref, xref0, yref0 
c
      integer::
     $     i, ibeg, iend, jfx
c--------------------------parameter transfer

      moving=control%moving
      strtp=control%strtp
      turb=control%turb

      gamma=control%gamma
      machinf=control%machinf
      reynolds=control%reynolds
      ronum=control%ronum
      tref=control%tref
      tintvl=control%tintvl

      inviscid=control%inviscid
      suther=control%suther

c--------------------------------------------------------

c
c --- Begin ---
c
c ... save mesh in previous physical time step
c
      if (nn.eq.1) call save_mesh(il, jl, kl, ilower, iupper,
     $     jlower, jupper, klower, kupper, x, y, z, xo, yo, zo,
     $     udi, vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, xsl,
     $     xzl, time, nmode)
c
c ... calculate forces in induced moving
c
      if (moving.eq.2) then
        if(turb.ne.0) call turb_new(bcwake, x, y, z, ilower, iupper,
     $      jlower, jupper, klower, kupper, il, jl, kl, nl,
     $      q, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $      bc_eta_upper, bc_zeta_lower, bc_zeta_upper, vol,
     $      visturb,
     $      np, rank, udi, vdi, wdi, udj,
     $      vdj, wdj, udk, vdk, wdk, cv1, qt, control)

        if ( strtp.eq.4 ) then
          call modforce(nmode, il, jl, kl, nl, x, y, z, q, gamma,
     $                  ilower, iupper, jlower, jupper, klower,
     $                  kupper, hx, hy, hz, pj)
        else
          call cldm(rank, np, il, jl, kl, nl, 
     $      ilower, iupper, jlower, jupper, klower, kupper,
     $      bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $      bc_eta_upper, bc_zeta_lower, bc_zeta_upper,
     $      x, y, z, q, visturb,
     $      udj, vdj, wdj, lref, aref, xref0, yref0, cl, cd, ct,
     $      control)
        end if
      end if
c
      call structure_para(nn, moving, time, vbd, xyz, dxyz,
     $     dstep, tintvl, cd, ct, cl, xsl, xzl, strtp, cs, mus, ub,
     $     kc, wh, walf, uinf, nfp, alf0, il, jl, kl, nmode, pj, wpw,
     $     hx, hy, hz, ux, uy, uz)
c
      select case (strtp)
      case (1)             !induced cylinder
c
        call mvgrid_cyl(dim2, il, jl, kl, nl, x, y, z, xo, yo, zo,
     $       ilower, iupper, jlower, jupper, klower, kupper,
     $       dxyz, jfx, beta)
c
      case (2)             !forced pitching airfoil
c
        call mvgrid_af(dim2, il, jl, kl, nl, x, y, z, xo, yo, zo,
     $        ilower, iupper, jlower, jupper, klower, kupper,
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper,
     $        dxyz, xref, yref, jfx, beta, xref0, yref0)
c
      case (3,5)           !induced airfoil
c
        call mvgrid_af(dim2, il, jl, kl, nl, x, y, z, xo, yo, zo,
     $        ilower, iupper, jlower, jupper, klower, kupper,
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper,
     $        dxyz, xref, yref, jfx, beta, xref0, yref0)
c
      case (4)             !induced wing
        call mvgridw(il, jl, kl, x, y, z, xo, yo, zo, ux, uy, uz,
     $               ilower, iupper, jlower, jupper, klower, kupper)
      end select
c
c ... grid moving vels
c
      call compute_qt(il, jl, kl, ilower, iupper, jlower,
     $        jupper, klower, kupper, x, y, z, xo, yo, zo, tintvl,
     $        dim2, qt, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper)
c
      end subroutine moving_grid
