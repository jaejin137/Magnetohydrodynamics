      subroutine integrate_ALL(nn, x,
     $     y, z, vol, q, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,  ntotal,
     $     time, il, jl, kl, nl, dim1, dim2, ilower, iupper, jlower,
     $     jupper, klower, kupper, resdulmax, loc, rank, np, 
     $     visturb, dt,
     $     resdulave,
     $     qi, qii,  qiii,
     $      dstep, ke, udi,
     $     vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, qt,
     $     resdl, sigma, nrbc_ex,
     $     iblnu, ipt, jpt, kpt, ic1, ic2, ic3, tko, cb1, cb2,
     $     cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4,
c     $     dist,  diver, control)      !original
     $     dist,  diver, control,       !jaejin
c=================================jaejin(begin)=========================
c       Arguments for MHD simulation
     $     intersteps, rho_inf, T_inf, c_sound, L_ref, mhdsim,
     $     ionmodel, T_onset, sigma_e_ref, magtype, magref, B_ref,
     $     I_in, offset, sigma_e_star, B_x_star, B_y_star, B_z_star,
     $     B_norm, source_mhd)
c=================================jaejin(end)===========================

c
c     integrates equations of motion for one time step
c     from time level "n" to time level "n+1"
c
c     IMPLICIT STATEMENT
c
      use datain

      implicit none

      type (datain_type), intent(in)::control
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     il, jl, kl, nl,
     $     dim1, dim2,
     $     nn
c
      integer::
     $     rhs_order,
     $     limiter,
     $     rhs_scheme,
     $     idimen,
     $     lhs_order,
     $     vis_order,
     $     iter_gs,
     $     dual_t,
     $     lhs_scheme, precondition,
     $     integrate,          ! integrate method
     $     gcl,
     $     moving,source,
     $     unidt, nrbc_ex,
     $     turb, main_dir
c     
c      integer::
      integer, intent(in)::
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl),
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     rank, np,
     $     dstep
c
      double precision,intent(in)::
     $     diver(ilower:iupper,jlower:jupper,klower:kupper)
c
      double precision, dimension(jlower:jupper,
     $     klower:kupper,4), intent(in):: udi, vdi, wdi
c
      double precision, dimension(ilower:iupper,
     $     klower:kupper,4), intent(in):: udj, vdj, wdj
c
      double precision, dimension(ilower:iupper,
     $     jlower:jupper,4), intent(in):: udk, vdk, wdk
c
      double precision, intent(in)::
     $     delta(5)
c
      logical::
     $     inviscid, suther
c
      double precision, intent(in)::
     $     sigma(2)
c
      double precision::
     $     machinf,
     $     epsfactor,
     $     kfactor,
     $     reynolds,
     $     cfl,
     $     gamma,
     $     prandtl,
     $     prt,
     $     tref,
     $     ptotal,
     $     ttotal,
     $     angl1, angl2,
     $     poutlet,
     $     tintvl,
     $     k_prec,
     $     rayleigh,
     $     froude, epsilon, ronum
c
      integer, intent(inout)::
     $     ntotal
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(in)::
     $  visturb(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)
c
      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol
c
      double precision, intent(inout),
     $     dimension(ilower:iupper, jlower:jupper, klower:kupper, nl)::
     $     q, qi, qii, qiii
c     
      double precision, intent(in),
     $     dimension(ilower:iupper,jlower:jupper,klower:kupper,3)::
     $     qt
c
      double precision, intent(in):: dist(il,jl,kl)
c
      double precision::
     $     time
c
      integer, intent(out)::
     $     loc(3)
c
      logical, intent(in)::
     $     ke                   ! energy considered in energy
c
      double precision, intent(out):: resdl
c
c     LOCAL VARIABLES
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     k_rk                 ! for RK method
c
      double precision, dimension(il+1, jl+1, kl+1, nl)::
     $     rflux, sflux, tflux,
     $     flux
c
      double precision::
     $     rk_f(il,jl,kl,nl,4),
     $     sr(il,jl,kl)
c
      double precision,
     $     dimension(ilower:iupper,jlower:jupper,klower:kupper,nl)::
     $     q_n
c
      double precision, dimension(il,jl,kl,nl)::
     $     rhs
c
      double precision, dimension(il,jl,kl)::
     $     residual(il, jl, kl)
c
      double precision, intent(out)::
     $     resdulmax,          ! field maximum L2-Norm residual
     $     resdulave
c
      double precision::
     $     dt(il,jl,kl),        ! time marching interval
     $     dtemp,
     $     gamma1,
     $     dual_time,
     $     pinlet,              ! pressure at inlet
     $     rsdbci, rsdbcj,      ! residuals for i & j bc eqns
     $     rhsbi(jl,kl,nl,2),
     $     rhsbj(il,kl,nl,2),
     $     cst_nri(nl,nl,jl,kl,5),
     $     ced_nri(nl,nl,jl,kl,5),
     $     cst_nrj(nl,nl,il,kl,5),
     $     ced_nrj(nl,nl,il,kl,5)
c
      logical::
     $     checktime = .False.,     ! if systime check necessary
     $     wallfunc
c
c     SA 1eq turbulent model contants
c     
      double precision, intent(in)::
     $  tko, cb1, cb2, cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4
      integer, intent(in):: iblnu, ipt, jpt, kpt, ic1, ic2, ic3
c
c     k-omega model contants
c     
      double precision, parameter::
     $     tkb = .09d0,
     $     twb = .075d0,
     $     tks = .5d0,
     $     tws = .5d0,
     $     alpha = .555556d0
c
      double precision::
     $     v2,constant
c
c-----------for precodition
      double precision::
     $     r,u,v,w,vel,
     $     v1,p,t,th,c,
     $     drhodp,drhodt,cp,
     $     urr,theta
      double precision, dimension(nl,nl)::
     $     p_matrix
c=================================jaejin(begin)=========================
c       MHD variables
      integer, intent(in)::
     $  mhdsim,
     $  intersteps

      character, intent(in)::
     $  ionmodel,
     $  magtype*10,
     $  magref

      double precision, intent(in),dimension(ilower:iupper,jlower:
     $  jupper,klower:kupper)::
     $  sigma_e_star,
     $  B_x_star, B_y_star, B_z_star, B_norm
      
      double precision, intent(in)::
     $  rho_inf,
     $  T_inf,
     $  c_sound,
     $  L_ref,
     $  T_onset,
     $  sigma_e_ref,
     $  B_ref,
     $  I_in,
     $  offset(3),
     $  source_mhd(il,jl,kl,2:5)

c=================================jaejin(end)===========================
                                                  

c--------------------------parameter transfer
     
      dual_t=control%dual_t
      gcl=control%gcl
      idimen=control%idimen
      integrate=control%integrate
      iter_gs=control%iter_gs
      lhs_order=control%lhs_order
      lhs_scheme=control%lhs_scheme
      limiter=control%limiter
      main_dir=control%main_dir    ! not used
      moving=control%moving
      nrbc_ex=control%nrbc_ex
      precondition=control%precondition
      rhs_order=control%rhs_order
      rhs_scheme=control%rhs_scheme
      source=control%source
      turb=control%turb
      unidt=control%unidt
      vis_order=control%vis_order  !not used


      inviscid=control%inviscid
      suther=control%suther

      angl1=control%angl1
      angl2=control%angl2
      cfl=control%cfl
      epsilon=control%epsilon
      epsfactor=control%epsfactor
      froude=control%froude
      gamma=control%gamma
      kfactor=control%kfactor
      k_prec=control%k_prec
      machinf=control%machinf
      poutlet=control%poutlet
      ptotal=control%ptotal
      prandtl=control%prandtl
      prt=control%prt
      rayleigh=control%rayleigh
      reynolds=control%reynolds
      ronum=control%ronum
      tref=control%tref
      tintvl=control%tintvl
      ttotal=control%ttotal
c--------------------------------------------------------

c
c *** SUBROUTINE START ***
c
      wallfunc = .False.
      if(turb.eq.2) wallfunc=.True.
      gamma1 = gamma - 1.d0
c
      if(unidt.ne.0.and.dual_t.ne.1) then
        dual_time = time - dt(1,1,1)
      else
        dual_time = time
      end if
c
      call rhside(control, dt, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta, 
     $     rhs, residual, q, rflux, sflux,
     $     tflux, x, y, z, vol, il, jl, kl, nl,
     $     dim2, ilower, iupper, jlower, jupper, klower, kupper, rank,
     $     visturb,  ntotal,
     $     qi, qii, dstep, wallfunc, tks, tws, tkb, twb,
     $     alpha, ke, sr, udi, vdi, wdi, udj, vdj, wdj,
     $     udk, vdk, wdk, qt, rhsbi, rhsbj,
     $     iblnu, ipt, jpt, kpt, ic1, ic2, ic3, tko, cb1, cb2,
     $     cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4,
c     $     dist,  qiii, diver) !original
     $     dist,  qiii, diver, !jaejin
c=================================jaejin(begin)=========================
c      Arguments for MHD simulation
     $     intersteps, rho_inf, T_inf, c_sound, L_ref, mhdsim,
     $     ionmodel, T_onset, sigma_e_ref, magtype, magref, B_ref,
     $     I_in, offset, sigma_e_star, B_x_star, B_y_star, B_z_star,
     $     B_norm, source_mhd)
c=================================jaejin(end)===========================
c      print *, "first rhside call finished"
c
      SELECT CASE (integrate)
c
c     AF method
c
      CASE (1)

         do i = 1, idimen
            select case (i)
            case (1)
               flux = rflux
            case (2)
               flux = sflux
            case (3)
               flux = tflux
            end select
            call invert(i, dt, bc_xi_lower, bc_xi_upper,
     $           bc_eta_lower, bc_eta_upper, bc_zeta_lower,
     $           bc_zeta_upper, delta,
     $           q, flux, gamma, tref, prandtl, rhs,
     $           x, y, z, vol, il, jl, kl, nl, dim2, dim1,
     $           ilower, iupper, jlower, jupper, klower, kupper, rank,
     $           visturb, wallfunc, ke, qt, 
     $           udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $           wdk, tko, qiii, diver, control)
         end do

      q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:) + rhs
c      print *, "CASE 1 passed"
c
c     R-K method
c
      CASE (2)

         q_n(1:il,1:jl,1:kl,:)=q(1:il,1:jl,1:kl,:)

         do k_rk = 1,2

            call updaterk(rank, k_rk, il, jl, kl, nl, q, q_n, rhs,
     $           ilower, iupper, jlower, jupper, klower, kupper) 

            rk_f(:,:,:,:,k_rk)=rhs

      call rhside(control, dt, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta, 
     $     rhs, residual, q, rflux, sflux,
     $     tflux, x, y, z, vol, il, jl, kl, nl,
     $     dim2, ilower, iupper, jlower, jupper, klower, kupper, rank,
     $     visturb,  ntotal,
     $     qi, qii, dstep, wallfunc, tks, tws, tkb, twb,
     $     alpha, ke, sr, udi, vdi, wdi, udj, vdj, wdj,
     $     udk, vdk, wdk, qt, rhsbi, rhsbj,
     $     iblnu, ipt, jpt, kpt, ic1, ic2, ic3, tko, cb1, cb2,
     $     cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4,
c     $     dist,  qiii, diver) !original
     $     dist,  qiii, diver, !jaejin
c=================================jaejin(begin)=========================
c      Arguments for MHD simulation
     $     intersteps, rho_inf, T_inf, c_sound, L_ref, mhdsim,
     $     ionmodel, T_onset, sigma_e_ref, magtype, magref, B_ref,
     $     I_in, offset, sigma_e_star, B_x_star, B_y_star, B_z_star,
     $     B_norm, source_mhd)
c=================================jaejin(end)===========================

         end do

         call updaterk(rank, 3, il, jl, kl, nl, q, q_n, rhs,
     $        ilower, iupper, jlower, jupper, klower, kupper)

         rk_f(:,:,:,:,3)=rhs

       call rhside(control, dt, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta, 
     $     rhs, residual, q, rflux, sflux,
     $     tflux, x, y, z, vol, il, jl, kl, nl,
     $     dim2, ilower, iupper, jlower, jupper, klower, kupper, rank,
     $     visturb,  ntotal,
     $     qi, qii, dstep, wallfunc, tks, tws, tkb, twb,
     $     alpha, ke, sr, udi, vdi, wdi, udj, vdj, wdj,
     $     udk, vdk, wdk, qt, rhsbi, rhsbj,
     $     iblnu, ipt, jpt, kpt, ic1, ic2, ic3, tko, cb1, cb2,
     $     cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4,
c     $     dist,  qiii, diver) !original
     $     dist,  qiii, diver, !jaejin
c=================================jaejin(begin)=========================
c      Arguments for MHD simulation
     $     intersteps, rho_inf, T_inf, c_sound, L_ref, mhdsim,
     $     ionmodel, T_onset, sigma_e_ref, magtype, magref, B_ref,
     $     I_in, offset, sigma_e_star, B_x_star, B_y_star, B_z_star,
     $     B_norm, source_mhd)
c=================================jaejin(end)===========================

         rk_f(:,:,:,:,4)=rhs

         rhs = 0.16666667d0*(rk_f(:,:,:,:,1)+ 2.0*rk_f(:,:,:,:,2)
     $        +2.0*rk_f(:,:,:,:,3) +rk_f(:,:,:,:,4))

         q(1:il,1:jl,1:kl,:) = q_n(1:il,1:jl,1:kl,:)+rhs
c         print *, "CASE 2 passed"
c
c     Euler method
c
      CASE (3)
         if (precondition.ge.2) then
           if (precondition.eq.1) then
             constant=gamma*machinf*machinf
             do k=1,kl
               do j=1,jl
                 do i=1,il
                  r=q(i,j,k,1)
                  u=q(i,j,k,2)/r
                  v=q(i,j,k,3)/r
                  w=q(i,j,k,4)/r
                  vel=0.5*(u**2+v**2+w**2)
                  v1=sqrt(2.0*vel)
                  p=(gamma-1.0)*(q(i,j,k,5)-r*vel)
                  t=gamma*machinf**2*p/r
                  th=(q(i,j,k,5)+p)/r
                  c=sqrt((gamma-1.0)*(th-vel))
                  drhodp=gamma/c**2
                  drhodt=-r/t
                  cp=c**2/((gamma-1.0)*t)
                  cp=1.d0/((gamma-1.0)*machinf**2)
                  if (v1.gt.c) then
                     urr=c
                  elseif (v1.lt.c/100000.d0) then
                     urr=c/100000.d0
                  else
                     urr=v1
                  end if 
                  urr=dmin1(c,dmax1(v1,k_prec))
                  theta=(1.0/urr**2-drhodt/(r*cp))
                  call pre_matrix(nl,cp,theta,drhodt,th,r,
     $                            u,v,w,p_matrix)
                  call ivsnr(p_matrix,nl,1e-6)                  
                  qiii(i,j,k,:) = qiii(i,j,k,:)+
     $                       matmul(p_matrix(:,:),rhs(i,j,k,:))
c------------------
                  v2 = 0.5*(qiii(i,j,k,2)**2+qiii(i,j,k,3)**2+
     $                   qiii(i,j,k,4)**2)
                  q(i,j,k,1) = constant*qiii(i,j,k,1)/
     $                          qiii(i,j,k,5)
                  q(i,j,k,2) = qiii(i,j,k,2)*q(i,j,k,1)
                  q(i,j,k,3) = qiii(i,j,k,3)*q(i,j,k,1)
                  q(i,j,k,4) = qiii(i,j,k,4)*q(i,j,k,1)
                  q(i,j,k,5) = qiii(i,j,k,1)/(gamma-1.0)+
     $                          q(i,j,k,1)*v2
                 end do
                end do
              end do  
           end if
         else       
            q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:)+rhs
         end if
c         print *, "CASE 3 passed"
       
c     
c     Gauss-Sidell method
c
      CASE (4)
c         print *, "CASE 4 started"
c         print *, "Calling ursn..."
         call ursn(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux, 
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $        rhs, il, jl, kl, nl, dim1, dim2, ilower, iupper,
     $        jlower, jupper, klower, kupper, visturb, 
     $        dstep, checktime, wallfunc, ke, qt, sr, udi,
     $        vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, cst_nri, ced_nri,
     $        cst_nrj, ced_nrj, tko, qiii, 
     $        diver, control)
c         print *, "ursn call passed"
         if (precondition.ge.1) then
            qiii(1:il,1:jl,1:kl,:) = qiii(1:il,1:jl,1:kl,:)+
     $                               q_n(1:il,1:jl,1:kl,:)
            if (precondition.ge.2) then
               do k=1,kl
                 do j=1,jl
                  do i=1,il
                     p=qiii(i,j,k,1)
                     t=qiii(i,j,k,5)
                     call qstat(p,t,0.0,r)
                     q(i,j,k,1) = 1.0d0/r
                     q(i,j,k,2) = qiii(i,j,k,2)*q(i,j,k,1)
                     q(i,j,k,3) = qiii(i,j,k,3)*q(i,j,k,1)
                     q(i,j,k,4) = qiii(i,j,k,4)*q(i,j,k,1)
                     v2 = 0.5*(qiii(i,j,k,2)**2+qiii(i,j,k,3)**2+
     $                    qiii(i,j,k,4)**2)
                     call enthalpy(th,t,r,p)
                     q(i,j,k,5) =  q(i,j,k,1)*(th+v2)-p 
                   end do
                 end do
               end do             
            else
              constant=gamma*machinf*machinf
              do k=1,kl
                do j=1,jl
                  do i=1,il
                     v2 = 0.5*(qiii(i,j,k,2)**2+qiii(i,j,k,3)**2+
     $                    qiii(i,j,k,4)**2)
                     q(i,j,k,1) = constant*qiii(i,j,k,1)/
     $                            qiii(i,j,k,5)
                     q(i,j,k,2) = qiii(i,j,k,2)*q(i,j,k,1)
                     q(i,j,k,3) = qiii(i,j,k,3)*q(i,j,k,1)
                     q(i,j,k,4) = qiii(i,j,k,4)*q(i,j,k,1)
                     q(i,j,k,5) = qiii(i,j,k,1)/(gamma-1.0)+
     $                            q(i,j,k,1)*v2
                  end do
                end do
              end do       
            end if
         else
            q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:)+
     $                               q_n(1:il,1:jl,1:kl,:)
         end if
c         print *, "CASE 4 passed"
c
c     LU-SGS method
c
      CASE (5)

         call lu_sgs(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux,
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $        rhs, il, jl, kl, nl, dim1, dim2, ilower, iupper,
     $        jlower, jupper, klower, kupper, visturb, 
     $        dstep, checktime, wallfunc, ke, qt,  sr,  udi,
     $        vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, cst_nri, ced_nri,
     $        cst_nrj, ced_nrj,control)

         q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:)+q_n(1:il,1:jl,1:kl,:)
c         print *, "CASE 5 passed"
c
c     LU-Gauss-Sidell method
c
      CASE (11)

         call ursn_1(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux, 
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $        rhs, il, jl, kl, nl, dim1, dim2, ilower, iupper,
     $        jlower, jupper, klower, kupper, visturb, 
     $        dstep, checktime, wallfunc, ke, qt, sr, udi,
     $        vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, cst_nri, ced_nri,
     $        cst_nrj, ced_nrj, tko, qiii, 
     $        diver, control)

         q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:)+q_n(1:il,1:jl,1:kl,:)
c         print *, "CASE 11 passed"
c
c     LU-Gauss-Sidell method
c
      CASE (7)

         call ursn_2(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux, 
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $        rhs, il, jl, kl, nl, dim1, dim2, ilower, iupper,
     $        jlower, jupper, klower, kupper, visturb, 
     $        dstep, checktime, wallfunc, ke, qt, sr, udi,
     $        vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, cst_nri, ced_nri,
     $        cst_nrj, ced_nrj, tko, qiii, 
     $        diver, control)

         q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:)+q_n(1:il,1:jl,1:kl,:)
c         print *, "CASE 7 passed"
c
c     Gauss-Sidell method
c
      CASE (6)

         call lu_ursn(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux, 
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $        rhs, il, jl, kl, nl, dim1, dim2, ilower, iupper,
     $        jlower, jupper, klower, kupper, visturb, 
     $        dstep, checktime, wallfunc, ke, qt, sr, udi,
     $        vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, cst_nri, ced_nri,
     $        cst_nrj, ced_nrj, tko, qiii, 
     $        diver, control)


         q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:)+q_n(1:il,1:jl,1:kl,:)
c         print *, "CASE 6 passed"
c
c     Gauss-Sidell method
c
      CASE (12)

         call lu_ursn_gs(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux, 
     $        bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $        bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $        rhs, il, jl, kl, nl, dim1, dim2, ilower, iupper,
     $        jlower, jupper, klower, kupper, visturb, 
     $        dstep, checktime, wallfunc, ke, qt, sr, udi,
     $        vdi, wdi, udj, vdj, wdj, udk, vdk, wdk, cst_nri, ced_nri,
     $        cst_nrj, ced_nrj, tko, qiii, 
     $        diver, control)
         q(1:il,1:jl,1:kl,:) = q(1:il,1:jl,1:kl,:)+q_n(1:il,1:jl,1:kl,:)
c         print *, "CASE 12 passed"

      END SELECT

c      print *, "Integrating ..."
c
c     check abnormal p, t values
c
      call negative(il, jl, kl, nl, q, gamma,
     $     ilower, iupper, jlower, jupper, klower, kupper, ke,
     $     qt, ronum)
c
c     record maximum residual
c
      resdulmax = maxval(residual)
      loc = maxloc(residual)
      resdulave = sum(residual)/(il*jl*kl)
c
c     --- non-reflecting bc ---
c
      rsdbci = 0.d0
      rsdbcj = 0.d0
c
c ... in i direction
c
      if ( nrbc_ex.eq.1 ) then
        dtemp = ( 1.d0+0.5d0*gamma1*machinf*machinf )**3.5d0
        pinlet = ptotal / dtemp
c
        call integ_nri_GS(gamma, x, y, z, vol, q, q_n,
     $    bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $    bc_zeta_lower, bc_zeta_upper, iter_gs, il, jl, kl, nl, dim1,
     $    dim2, ilower, iupper, jlower, jupper, klower, kupper,
     $    rsdbci, dual_t, tintvl, qi, qii, dual_time, dt, pinlet,
     $    poutlet, machinf, ptotal, ttotal, sigma, rhsbi, cst_nri,
     $    ced_nri)
      end if
c
c ... in j direction
c
      if ( nrbc_ex.eq.2 ) then
        dtemp = ( 1.d0+0.5d0*gamma1*machinf*machinf )**3.5d0
        pinlet = ptotal / dtemp
c
        call integ_nrj_GS(gamma, x, y, z, vol, q, q_n,
     $    bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $    bc_zeta_lower, bc_zeta_upper, iter_gs, il, jl, kl, nl, dim1,
     $    dim2, ilower, iupper, jlower, jupper, klower, kupper,
     $    rsdbcj, dual_t, tintvl, qi, qii, dual_time, dt, pinlet,
     $    poutlet, rhsbj, cst_nrj, ced_nrj)
      end if
c
c ... in both i & j direction
c
      if ( nrbc_ex.eq.3 ) then
        dtemp = ( 1.d0+0.5d0*gamma1*machinf*machinf )**3.5d0
        pinlet = ptotal / dtemp
c
        call integ_nri_GS(gamma, x, y, z, vol, q, q_n,
     $    bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $    bc_zeta_lower, bc_zeta_upper, iter_gs, il, jl, kl, nl, dim1,
     $    dim2, ilower, iupper, jlower, jupper, klower, kupper,
     $    rsdbci, dual_t, tintvl, qi, qii, dual_time, dt, pinlet,
     $    poutlet, machinf, ptotal, ttotal, sigma, rhsbi, cst_nri,
     $    ced_nri)
c
        call integ_nrj_GS(gamma, x, y, z, vol, q, q_n,
     $    bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $    bc_zeta_lower, bc_zeta_upper, iter_gs, il, jl, kl, nl, dim1,
     $    dim2, ilower, iupper, jlower, jupper, klower, kupper,
     $    rsdbcj, dual_t, tintvl, qi, qii, dual_time, dt, pinlet,
     $    poutlet, rhsbj, cst_nrj, ced_nrj)
      end if
c
      resdl = rsdbci
      if ( rsdbcj.gt.resdl ) resdl = rsdbcj
c
      return
      end
