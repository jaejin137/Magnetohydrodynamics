      subroutine lhs_matrix(index, i1, i2, dt, x, y, z, vol, bclower,
     $   bcupper, delta, wallfunc, capa, capb, capc,
     $   q, flux, il, jl, kl, nl, dim2, dim1, ilower, iupper, jlower,
     $   jupper, klower, kupper, visturb, ke, udi, vdi, wdi, udj, vdj,
     $   wdj, udk, vdk, wdk, qtline, ci_nr, tko,  qiii, diverg,
     $   control)
c     
c     calculate the LHS matrices in xi-direction
c     
c     IMPLICIT STATEMENT
c     
      use datain

      implicit none

      type (datain_type), intent(in)::control
c     
c     INTERFACE VARIABLES

      integer, intent(in)::
     $     index,               ! 1-xi, 2-eta, 3-zeta
     $     i1, i2,              ! j k,  i k,   i j
     $     il, jl, kl, nl,
     $     dim1, dim2
      integer::
     $     blen,
     $     moving, precondition,
     $     main_dir

      logical, intent(in)::
     $     wallfunc, ke

      integer, intent(in)::
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     bclower, bcupper

      double precision, intent(in)::
     $     delta(5)

      double precision:: qconer(nl,2)

      integer::
     $     rhs_order,
     $     limiter,
     $     rhs_scheme,
     $     lhs_scheme,
     $     idimen,
     $     lhs_order

      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z

      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol

      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q, qiii

      double precision,intent(in)::
     $     diverg(ilower:iupper,jlower:jupper,klower:kupper)

      double precision, dimension(jlower:jupper,
     $     klower:kupper,4), intent(in):: udi, vdi, wdi

      double precision, dimension(ilower:iupper,
     $     klower:kupper,4), intent(in):: udj, vdj, wdj

      double precision, dimension(ilower:iupper,
     $     jlower:jupper,4), intent(in):: udk, vdk, wdk

      double precision, intent(in)::
     $  visturb(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)

      double precision, dimension(nl,nl,dim1), intent(out)::
     $     capa, capb, capc

      double precision, intent(in)::
     $     dt(il,jl,kl)        ! time marching interval
      double precision::
     $     gamma,
     $     prandtl,
     $     prt,
     $     tref,
     $     ptotal,
     $     ttotal,
     $     angl1,
     $     angl2,
     $     poutlet,
     $     k_prec,
     $     ronum

      double precision, intent(in)::
     $     flux(il+1, jl+1, kl+1, nl), ! e, f, g surface flux
     $     qtline(ilower:dim2,3)

      double precision::
     $     machinf, epsfactor, kfactor, reynolds

      logical::
     $     inviscid, suther
c
c     SA 1eq turbulent model contants
c
      double precision, intent(in):: tko

      double precision, dimension(nl,nl,3,2), intent(out)::
     $     ci_nr
c     
c     LOCAL VARIABLES
c     
      double precision, dimension(dim2)::
     $     capul, capur,
     $     capvl, capvr,
     $     capwl, capwr,        
     >     el, er,
     $     pl, pr,
     $     ql, qr,
     >     rhoinvl, rhoinvr,
     $     rhol, rhor,
     $     rhoul, rhour,
     $     rhovl, rhovr,
     $     rhowl, rhowr,
     $     rhoel, rhoer,
     $     rhotkl, rhotkr,
     $     rhotwl, rhotwr,
     $     sqrtinv,
     $     sqrtrhol, sqrtrhor,
     >     tenthalpyl, tenthalpyr,
     $     ul, ur,
     $     vl, vr,
     $     wl, wr,
     $     qtl, qtr,tl,tr

      double precision, dimension(ilower:dim2)::
     $     p, qq,
     $     rho, rhoe, rhou, rhov, rhow,
     $     rinv, t, u, v, w
      double precision, dimension(ilower:dim2)::
     $     diver

      double precision::
     $     fflux(dim2,nl)
      
      double precision, dimension(ilower+1:dim2)::
     $     lx, ly, lz,
     $     mx, my, mz,
     $     nx, ny, nz,
     $     lxhat, lyhat, lzhat,
     $     mxhat, myhat, mzhat,
     $     nxhat, nyhat, nzhat

      double precision, dimension(nl,nl,dim2)::
     $     a

      double precision  al(dim2,nl,nl),    ar(dim2,nl,nl),
     >     ll(dim2,nl,nl),    lr(dim2,nl,nl)
c     
      double precision::
     $     dtv(dim1),
     >     glx(dim2), gly(dim2), glz(dim2),
     >     molec(0:dim2),
     >     term1(dim2), term2(dim2), term3(dim2),
     >     thermal(dim2),
     >     viscous(dim2),
     $     volumeinv(ilower:dim2),
     $     vbds(3,4)

      double precision::
     $     coeff1, coeff2, coeff3, ! temporary coefficient
     $     c1,  gamma1

      integer::
     $     n,
     $     ii, jj,              ! iteration variable
     $     i,                   ! iteration variable
     $     imax                 ! surface number of current line

c     variables for combination

c      double precision, dimension(dim2)::
      double precision, dimension(ilower+1:dim2)::
     $     llx, lly, llz,
     $     mmx, mmy, mmz,
     $     nnx, nny, nnz,
     $     llxhat, llyhat, llzhat,
     $     mmxhat, mmyhat, mmzhat,
     $     nnxhat, nnyhat, nnzhat
     $                          ! metric in xi, eta, zeta direction
     $                          ! when index: 1,  2,  3

c      double precision, dimension(dim2)::
c     $     llx, lly, llz,
c     $     mmx, mmy, mmz,
c     $     nnx, nny, nnz,
c     $     llxhat, llyhat, llzhat,
c     $     mmxhat, mmyhat, mmzhat,
c     $     nnxhat, nnyhat, nnzhat,
c     $                          ! metric in xi, eta, zeta direction
c     $                          ! when index: 1,  2,  3
c     $     capvell, capvelr,    ! capu, capv, capw
c     $     alphak, alphaw,
c     $     volume(ilower+1:dim2),      ! viscosity on current line
c     $     vist(ilower+1:dim2),        ! viscosity on current line
c     $     qlocal(ilower:dim2,nl)

      double precision, dimension(dim2)::
     $     capvell, capvelr,    ! capu, capv, capw
     $     alphak, alphaw,
     $     volume(ilower:dim2),      ! viscosity on current line
     $     vist(ilower:dim2),        ! viscosity on current line
     $     qlocal(ilower:dim2,nl)

      integer::
     $     istart, istop

c..   VARIABLES FOR ZHA SCHEME
      
      double precision::
     $     dhat_l(dim2,nl,nl),
     $     dhat_r(dim2,nl,nl)

C..   VARIABLES for Van Leer scheme

      double precision::
     $     area_,               ! interface area
     $     vx1,vy1,vz1,         ! interface unit vector components
     $     u_l,v_l,w_l,         ! left velocity components in x,y,z dir
     $     d_l,                 ! left density
     $     e_l,                 ! left density      
     $     tk_l,                 ! left density      
     $     c_l,                 ! left speed of sound
     $     u_r,v_r,w_r,         ! right velocity components in x,y,z dir
     $     d_r,                 ! right density
     $     e_r,                 ! right density      
     $     tk_r,                 ! right density      
     $     c_r,                 ! right speed of sound
     $     aml(nl,nl),          ! left jacobian of van leer flux
     $     amr(nl,nl),          ! right jacobian of van leer flux
     $     qt_l, qt_r           ! GRID VELOCITIES AT LEFT AND RIGHT SIDES
c
c..   k-omega variables
      double precision, dimension(ilower:dim2)::
     $     rhotk, tk,
     $     rhotw, tw

      double precision, dimension(dim2)::
     $     tkl, tkr,
     $     twl, twr

      double precision, parameter::
     $     tks = .5d0,          ! sigma in k equation
     $     tws = .5d0           ! sigma in omega equation
c--------------for precondition---------------
      double precision::
     $     cl,cr,c,temp
c          urr

      double precision, dimension(0:dim2)::
     $     drhodpl,drhodpr,
     $     drhodtl,drhodtr,
     $     cpl,cpr

      double precision, dimension(nl,nl)::
     $     k_matrix
c
      double precision::
     $     ve, we, psi, theta, phil, phir
c     
c     *** SUBROUTINE START ***
c     
c     preliminary
c     
c--------------------------parameter transfer
      blen=control%blen
      idimen=control%idimen
      lhs_order=control%lhs_order
      lhs_scheme=control%lhs_scheme
      limiter=control%limiter
      main_dir=control%main_dir
      moving=control%moving
      precondition=control%precondition
      rhs_order=control%rhs_order
      rhs_scheme=control%rhs_scheme

      inviscid=control%inviscid
      suther=control%suther

      angl1=control%angl1
      angl2=control%angl2
      epsfactor=control%epsfactor
      gamma=control%gamma
      kfactor=control%kfactor
      k_prec=control%k_prec
      machinf=control%machinf
      poutlet=control%poutlet
      ptotal=control%ptotal
      prandtl=control%prandtl
      prt=control%prt
      reynolds=control%reynolds
      ronum=control%ronum
      tref=control%tref
      ttotal=control%ttotal
c--------------------------------------------------------
      gamma1 = gamma - 1.d0
      c1    = 2.d0 - gamma
      volume = 0.d0

c     compute metrics  (see subroutine rhs.f for explanation)

      select case (index)
      case (1)
         imax = il+1
      case (2)
         imax = jl+1
      case (3)
         imax = kl+1
      end select

      call metric(index,i1, i2, imax, il, jl, kl, x, y, z, lx, ly, lz,
     $     mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat, mxhat, myhat,
     $     mzhat, nxhat, nyhat, nzhat, ilower, iupper, jlower, jupper,
     $     klower, kupper, bclower, bcupper)

      select case (index)
      case (1)
         llx = lx
         lly = ly
         llz = lz
         mmx = mx
         mmy = my
         mmz = mz
         nnx = nx
         nny = ny
         nnz = nz
         llxhat = lxhat
         llyhat = lyhat
         llzhat = lzhat
         mmxhat = mxhat
         mmyhat = myhat
         mmzhat = mzhat
         nnxhat = nxhat
         nnyhat = nyhat
         nnzhat = nzhat
         istart = ilower
         istop = iupper
         volumeinv(ilower+1:iupper-1) =1.d0/vol(:,i1,i2)
         volume(ilower+1:iupper-1) =vol(:,i1,i2)
         vist(ilower+1:iupper-1) = visturb(:,i1,i2)
         dtv(:imax-1) = dt(:,i1,i2)
         fflux(:imax,:) = flux(:,i1,i2,:)
      case (2)
         llx = mx
         lly = my
         llz = mz
         mmx = nx
         mmy = ny
         mmz = nz
         nnx = lx
         nny = ly
         nnz = lz
         llxhat = mxhat
         llyhat = myhat
         llzhat = mzhat
         mmxhat = nxhat
         mmyhat = nyhat
         mmzhat = nzhat
         nnxhat = lxhat
         nnyhat = lyhat
         nnzhat = lzhat
         istart = jlower
         istop = jupper
         volumeinv(jlower+1:jupper-1) =1.d0/vol(i1,:,i2)
         volume(jlower+1:jupper-1) =vol(i1,:,i2)
         vist(jlower+1:jupper-1) = visturb(i1,:,i2)
         dtv(:imax-1) = dt(i1,:,i2)
         fflux(:imax,:) = flux(i1,:,i2,:)
      case (3)
         llx = nx
         lly = ny
         llz = nz
         mmx = lx
         mmy = ly
         mmz = lz
         nnx = mx
         nny = my
         nnz = mz
         llxhat = nxhat
         llyhat = nyhat
         llzhat = nzhat
         mmxhat = lxhat
         mmyhat = lyhat
         mmzhat = lzhat
         nnxhat = mxhat
         nnyhat = myhat
         nnzhat = mzhat
         istart = klower
         istop = kupper
         volumeinv(klower+1:kupper-1) =1.d0/vol(i1,i2,:)
         volume(klower+1:kupper-1) =vol(i1,i2,:)
         vist(klower+1:kupper-1) = visturb(i1,i2,:)
         dtv(:imax-1) = dt(i1,i2,:)
         fflux(:imax,:) = flux(i1,i2,:,:)
      end select

c     
c     collect variables as vector
c     
      select case (index)
      case (1)
         if (precondition.ge.2) then
            qlocal(ilower:iupper,:) = qiii(:,i1,i2,:)
            diver(ilower:iupper) = diverg(:,i1,i2)
         else
            qlocal(ilower:iupper,:) = q(:,i1,i2,:)
         end if
         if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            vbds(1,:) = udi(i1,i2,:)
            vbds(2,:) = vdi(i1,i2,:)
            vbds(3,:) = wdi(i1,i2,:)
         end if
      case (2)
         if (precondition.ge.2) then
            qlocal(jlower:jupper,:) = qiii(i1,:,i2,:)
            diver(jlower:jupper) = diverg(i1,:,i2)
         else
            qlocal(jlower:jupper,:) = q(i1,:,i2,:)
         end if
         if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            vbds(1,:) = udj(i1,i2,:)
            vbds(2,:) = vdj(i1,i2,:)
            vbds(3,:) = wdj(i1,i2,:)
         end if
      case (3)
         if (precondition.ge.2) then
            qlocal(klower:kupper,:) = qiii(i1,i2,:,:)
            diver(klower:kupper) = diverg(i1,i2,:)
         else
            qlocal(klower:kupper,:) = q(i1,i2,:,:)
         end if
         if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            vbds(1,:) = udk(i1,i2,:)
            vbds(2,:) = vdk(i1,i2,:)
            vbds(3,:) = wdk(i1,i2,:)
         end if
      end select

      do i = istart, istop
         rho(i) = qlocal(i,1)
         rhou(i) = qlocal(i,2)
         rhov(i) = qlocal(i,3)
         rhow(i) = qlocal(i,4)
         rhoe(i) = qlocal(i,5)
         if(nl.eq.6) then
            rhotk(i) = qlocal(i,6)
         else if(nl.eq.7) then
            rhotk(i) = qlocal(i,6)
            rhotw(i) = qlocal(i,7)
         end if
      end do
c     
c     boundary conditions
c     
      if (precondition.ge.2) then
         select case (index)
         case (1)
           qconer(:,1) = qiii(0,i1,i2,:)
           qconer(:,2) = qiii(il+1,i1,i2,:)
         case (2)
           qconer(:,1) = qiii(i1,0,i2,:)
           qconer(:,2) = qiii(i1,jl+1,i2,:)
         case (3)
           qconer(:,1) = qiii(i1,i2,0,:)
           qconer(:,2) = qiii(i1,i2,kl+1,:)
         end select
      else
         select case (index)
         case (1)
           qconer(:,1) = q(0,i1,i2,:)
           qconer(:,2) = q(il+1,i1,i2,:)
         case (2)
           qconer(:,1) = q(i1,0,i2,:)
           qconer(:,2) = q(i1,jl+1,i2,:)
         case (3)
           qconer(:,1) = q(i1,i2,0,:)
           qconer(:,2) = q(i1,i2,kl+1,:)
         end select
      end if
c
c      call boundary(index, dim2, i1, i2, il, jl, kl, nl,
c     $     ilower, iupper, jlower, jupper, klower, bclower, bcupper,
c     $     q, gamma, machinf, kupper, ptotal, ttotal, angl1, angl2,
c     $     poutlet, rho, rhoe, rhou, rhov, rhow, rhotk, rhotw, ke,
c     $     moving, vbds, x, y, z, blen, tintvl, qconer,
c     $     precondition,qiii,diver,index_lower,index_upper)
c
      if (precondition.ge.2) then
         do i = istart, istop
            u(i) = rhou(i)
            v(i) = rhov(i)
            w(i) = rhow(i)
            qq(i) = 0.5d0 * (u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
            p(i) = rho(i)
            t(i) = rhoe(i)
            call qstat(rho(i),rhoe(i),diver(i),rinv(i))
         end do
      else
         do i = istart, istop
            rinv(i) = 1.d0 / rho(i)
            u(i) = rinv(i) * rhou(i)
            v(i) = rinv(i) * rhov(i)
            w(i) = rinv(i) * rhow(i)
            qq(i) = 0.5d0 * (u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
            if (nl.eq.6) then
              tk(i) = rhotk(i) * rinv(i)
              if ( dabs(ronum).gt.1e-9 ) then
                ve = qtline(i,2)
                we = qtline(i,3)
                psi = ve*ve+we*we
                theta = qq(i)-0.5d0*psi
                p(i) = gamma1*( rhoe(i)-rho(i)*theta )
              else
                p(i) = gamma1*( rhoe(i)-rho(i)*qq(i) )
              end if
            else if(nl.eq.7.and.ke) then
               tk(i) = rhotk(i) * rinv(i)
               tw(i) = rhotw(i) * rinv(i)
               p(i) = gamma1 * (rhoe(i)-rhotk(i)-rho(i)*qq(i))
            else
               if ( dabs(ronum).gt.1e-9 ) then
                 ve = qtline(i,2)
                 we = qtline(i,3)
                 psi = ve*ve+we*we
                 theta = qq(i)-0.5d0*psi
                 p(i) = gamma1*( rhoe(i)-rho(i)*theta )
               else
                 p(i) = gamma1*( rhoe(i)-rho(i)*qq(i) )
               end if
            end if
            t(i) = p(i) * gamma * machinf * machinf / rho(i)
         end do
      end if
c     
c     determine left and right states for flow variables
c     
      call reconstruct(index,imax, lhs_order,dim2, il, jl, kl, nl,
     $     lx, ly, lz,
     $     mx, my, mz, nx, ny, nz, rho, rhou, rhov, rhow, rhoe,
     $     capul, capur, capvl, capvr, capwl, capwr, el, er, pl, pr,
     $     ql, qr, ul, ur, vl, vr, wl, wr, rhoinvl, rhoinvr, rhol, rhor,
     $     rhoul, rhour, rhovl, rhovr, rhowl, rhowr, rhoel, rhoer,
     $     rhotkl, rhotkr, rhotwl, rhotwr,
     $     sqrtinv, sqrtrhol, sqrtrhor, tenthalpyl, tenthalpyr,
     $     bclower, bcupper, rhotk, rhotw, tkl, tkr, twl, twr,
     $     .False., ke, qtline, qtl, qtr, ilower,
     $     tl, tr, diver, control)
      
      select case (index)
      case (1)
         capvell = capul
         capvelr = capur
      case (2)
         capvell = capvl
         capvelr = capvr
      case (3)
         capvell = capwl
         capvelr = capwr
      end select
c     
c     compute Jacobian matrices using left and right face values
c     
      al = 0.d0
      ar = 0.d0
c
c-------------------precondition--------------------
      if (precondition.ge.1) then
         if (precondition.ge.2) then
            do i=1,imax
                call qstat(pl(i),tl(i),diver(i),temp)
                rhol(i)=1.0d0/temp
                call cpt(tl(i),cpl(i))
                call drhdtpp(temp,tl(i),drhodtl(i))
                call drhdptt(temp,tl(i),drhodpl(i))
               
                call qstat(pr(i),tr(i),diver(i),temp)
                rhor(i)=1.0d0/temp
                call cpt(tr(i),cpr(i))
                call drhdtpp(temp,tr(i),drhodtr(i))
                call drhdptt(temp,tr(i),drhodpr(i)) 
c--------------------same as 1
            end do
         else
            do i=1,imax
               tl(i)=gamma*machinf**2*pl(i)/rhol(i)
               cl=sqrt((gamma-1.0)*(tenthalpyl(i)-ql(i)))
               drhodpl(i)=gamma/cl**2
               drhodtl(i)=-rhol(i)/tl(i)
               cpl(i)=cl**2/((gamma-1.0)*tl(i))
               cpl(i)=1.d0/((gamma-1.0)*machinf**2)

               tr(i)=gamma*machinf**2*pr(i)/rhor(i)
               cr=sqrt((gamma-1.0)*(tenthalpyr(i)-qr(i)))
               drhodpr(i)=gamma/cr**2
               drhodtr(i)=-rhor(i)/tr(i)
               cpr(i)=cr**2/((gamma-1.0)*tr(i))
               cpr(i)=1.d0/((gamma-1.0)*machinf**2)

            end do
         end if
c---------------------------------
      do i=1,imax
         al(i,1,1) = drhodpl(i)*capvell(i)
         ar(i,1,1) = drhodpr(i)*capvelr(i)
         al(i,1,2) = rhol(i)*llx(i)
         ar(i,1,2) = rhor(i)*llx(i)
         al(i,1,3) = rhol(i)*lly(i)
         ar(i,1,3) = rhor(i)*lly(i)
         al(i,1,4) = rhol(i)*llz(i)
         ar(i,1,4) = rhor(i)*llz(i)
         al(i,1,5) = drhodtl(i)*capvell(i)
         ar(i,1,5) = drhodtr(i)*capvelr(i)


         al(i,2,1) = drhodpl(i)*ul(i)*capvell(i) + llx(i)
         ar(i,2,1) = drhodpr(i)*ur(i)*capvelr(i) + llx(i)
         al(i,2,2) = rhol(i)*capvell(i) + rhol(i)*ul(i)*llx(i)
         ar(i,2,2) = rhor(i)*capvelr(i) + rhor(i)*ur(i)*llx(i)
         al(i,2,3) = rhol(i)*ul(i)*lly(i)
         ar(i,2,3) = rhor(i)*ur(i)*lly(i)
         al(i,2,4) = rhol(i)*ul(i)*llz(i)
         ar(i,2,4) = rhor(i)*ur(i)*llz(i)
         al(i,2,5) = drhodtl(i)*ul(i)*capvell(i)
         ar(i,2,5) = drhodtr(i)*ur(i)*capvelr(i)

         glx(i) = gamma1*llx(i)
         if(nl.eq.7.and.ke) then
            al(i,2,6) = -glx(i)
            ar(i,2,6) = -glx(i)
         end if

         al(i,3,1) = drhodpl(i)*vl(i)*capvell(i) + lly(i)
         ar(i,3,1) = drhodpr(i)*vr(i)*capvelr(i) + lly(i)
         al(i,3,2) = rhol(i)*vl(i)*llx(i)
         ar(i,3,2) = rhor(i)*vr(i)*llx(i)
         al(i,3,3) = rhol(i)*capvell(i) + rhol(i)*vl(i)*lly(i)
         ar(i,3,3) = rhor(i)*capvelr(i) + rhor(i)*vr(i)*lly(i)
         al(i,3,4) = rhol(i)*vl(i)*llz(i)
         ar(i,3,4) = rhor(i)*vr(i)*llz(i)
         al(i,3,5) = drhodtl(i)*vl(i)*capvell(i)
         ar(i,3,5) = drhodtr(i)*vr(i)*capvelr(i)

         gly(i) = gamma1*lly(i)
         if(nl.eq.7.and.ke) then
            al(i,3,6) = -gly(i)
            ar(i,3,6) = -gly(i)
         end if

         al(i,4,1) = drhodpl(i)*wl(i)*capvell(i) + llz(i)
         ar(i,4,1) = drhodpr(i)*wr(i)*capvelr(i) + llz(i)
         al(i,4,2) = rhol(i)*wl(i)*llx(i)
         ar(i,4,2) = rhor(i)*wr(i)*llx(i)
         al(i,4,3) = rhol(i)*wl(i)*lly(i)
         ar(i,4,3) = rhor(i)*wr(i)*lly(i)
         al(i,4,4) = rhol(i)*capvell(i) + rhol(i)*wl(i)*llz(i)
         ar(i,4,4) = rhor(i)*capvelr(i) + rhor(i)*wr(i)*llz(i)
         al(i,4,5) = drhodtl(i)*wl(i)*capvell(i)
         ar(i,4,5) = drhodtr(i)*wr(i)*capvelr(i)


         glz(i) = gamma1*llz(i)
         if(nl.eq.7.and.ke) then
            al(i,4,6) = -glz(i)
            ar(i,4,6) = -glz(i)
         end if


         al(i,5,1) = drhodpl(i)*tenthalpyl(i)*capvell(i)
         ar(i,5,1) = drhodpr(i)*tenthalpyr(i)*capvelr(i)
         al(i,5,2) = rhol(i)*ul(i)*capvell(i)+
     $               rhol(i)*tenthalpyl(i)*llx(i)
         ar(i,5,2) = rhor(i)*ur(i)*capvelr(i)+
     $               rhor(i)*tenthalpyr(i)*llx(i)
         al(i,5,3) = rhol(i)*vl(i)*capvell(i)+
     $               rhol(i)*tenthalpyl(i)*lly(i)
         ar(i,5,3) = rhor(i)*vr(i)*capvelr(i)+
     $               rhor(i)*tenthalpyr(i)*lly(i)
         al(i,5,4) = rhol(i)*wl(i)*capvell(i)+
     $               rhol(i)*tenthalpyl(i)*llz(i)
         ar(i,5,4) = rhor(i)*wr(i)*capvelr(i)+
     $               rhor(i)*tenthalpyr(i)*llz(i)
         al(i,5,5) = (drhodtl(i)*tenthalpyl(i)+
     $               rhol(i)*cpl(i))*capvell(i)
         ar(i,5,5) = (drhodtr(i)*tenthalpyr(i)+
     $               rhor(i)*cpr(i))*capvelr(i)



         if(nl.eq.7.and.ke) then
            al(i,5,1) = -gamma*capvell(i)*el(i) + 
     >           gamma1*capvell(i)*(2.d0*ql(i)+tkl(i))
            ar(i,5,1) = -gamma*capvelr(i)*er(i) + 
     >           gamma1*capvelr(i)*(2.d0*qr(i)+tkr(i))
         end if

         if(nl.eq.6) then
            al(i,6,1) = -tkl(i)*capvell(i)
            ar(i,6,1) = -tkr(i)*capvelr(i)
            al(i,6,2) = llx(i)*tkl(i)
            ar(i,6,2) = llx(i)*tkr(i)
            al(i,6,3) = lly(i)*tkl(i)
            ar(i,6,3) = lly(i)*tkr(i)
            al(i,6,4) = llz(i)*tkl(i)
            ar(i,6,4) = llz(i)*tkr(i)
            al(i,6,6) = capvell(i)
            ar(i,6,6) = capvelr(i)
         else if(nl.eq.7) then
            if(ke) then
               al(i,5,6) = -gamma1*capvell(i)
               ar(i,5,6) = -gamma1*capvelr(i)
            end if

            al(i,6,1) = -tkl(i)*capvell(i)
            ar(i,6,1) = -tkr(i)*capvelr(i)
            al(i,6,2) = llx(i)*tkl(i)
            ar(i,6,2) = llx(i)*tkr(i)
            al(i,6,3) = lly(i)*tkl(i)
            ar(i,6,3) = lly(i)*tkr(i)
            al(i,6,4) = llz(i)*tkl(i)
            ar(i,6,4) = llz(i)*tkr(i)
            al(i,6,6) = capvell(i)
            ar(i,6,6) = capvelr(i)

            al(i,7,1) = -twl(i)*capvell(i)
            ar(i,7,1) = -twr(i)*capvelr(i)
            al(i,7,2) = llx(i)*twl(i)
            ar(i,7,2) = llx(i)*twr(i)
            al(i,7,3) = lly(i)*twl(i)
            ar(i,7,3) = lly(i)*twr(i)
            al(i,7,4) = llz(i)*twl(i)
            ar(i,7,4) = llz(i)*twr(i)
            al(i,7,7) = capvell(i)
            ar(i,7,7) = capvelr(i)
         end if

         if(moving.ge.1) then
            if(nl.eq.7) then
               print *, '(lhs_matrix)',
     $              ' moving is not working well with kw yet'
               stop
            end if
            do n = 1, nl
               al(i,n,n) = al(i,n,n)+qtl(i)
               ar(i,n,n) = ar(i,n,n)+qtr(i)
            end do
         end if

      end do

      else
c-----------------------not precondition-------------
c     
      do i = 1,imax
         if ( dabs(ronum).gt.1e-9 ) then
           ve = qtline(i-1,2)
           we = qtline(i-1,3)
           psi = ve*ve+we*we
           phil   = ql(i)+0.5d0*psi
           ve = qtline(i,2)
           we = qtline(i,3)
           psi = ve*ve+we*we
           phir   = qr(i)+0.5d0*psi
         else
           phil   = ql(i)
           phir   = qr(i)
         end if
c
         al(i,1,2) = llx(i)
         ar(i,1,2) = llx(i)
         al(i,1,3) = lly(i)
         ar(i,1,3) = lly(i)
         al(i,1,4) = llz(i)
         ar(i,1,4) = llz(i)

         glx(i) = gamma1*llx(i)

         al(i,2,1) = -ul(i)*capvell(i) + glx(i)*phil
         ar(i,2,1) = -ur(i)*capvelr(i) + glx(i)*phir
         al(i,2,2) = capvell(i) + c1*llx(i)*ul(i)
         ar(i,2,2) = capvelr(i) + c1*llx(i)*ur(i)
         al(i,2,3) = ul(i)*lly(i) - glx(i)*vl(i)
         ar(i,2,3) = ur(i)*lly(i) - glx(i)*vr(i)
         al(i,2,4) = ul(i)*llz(i) - glx(i)*wl(i)
         ar(i,2,4) = ur(i)*llz(i) - glx(i)*wr(i)
         al(i,2,5) = glx(i)
         ar(i,2,5) = glx(i)
         if(nl.eq.7.and.ke) then
            al(i,2,6) = -glx(i)
            ar(i,2,6) = -glx(i)
         end if

         gly(i) = gamma1*lly(i)

         al(i,3,1) = -vl(i)*capvell(i) + gly(i)*phil
         ar(i,3,1) = -vr(i)*capvelr(i) + gly(i)*phir
         al(i,3,2) = vl(i)*llx(i) - gly(i)*ul(i)
         ar(i,3,2) = vr(i)*llx(i) - gly(i)*ur(i)
         al(i,3,3) = capvell(i) + c1*lly(i)*vl(i)
         ar(i,3,3) = capvelr(i) + c1*lly(i)*vr(i)
         al(i,3,4) = vl(i)*llz(i) - gly(i)*wl(i)
         ar(i,3,4) = vr(i)*llz(i) - gly(i)*wr(i)
         al(i,3,5) = gly(i)
         ar(i,3,5) = gly(i)
         if(nl.eq.7.and.ke) then
            al(i,3,6) = -gly(i)
            ar(i,3,6) = -gly(i)
         end if

         glz(i) = gamma1*llz(i)

         al(i,4,1) = -wl(i)*capvell(i) + glz(i)*phil
         ar(i,4,1) = -wr(i)*capvelr(i) + glz(i)*phir
         al(i,4,2) = wl(i)*llx(i) - glz(i)*ul(i)
         ar(i,4,2) = wr(i)*llx(i) - glz(i)*ur(i)
         al(i,4,3) = wl(i)*lly(i) - glz(i)*vl(i)
         ar(i,4,3) = wr(i)*lly(i) - glz(i)*vr(i)
         al(i,4,4) = capvell(i) + c1*llz(i)*wl(i)
         ar(i,4,4) = capvelr(i) + c1*llz(i)*wr(i)
         al(i,4,5) = glz(i)
         ar(i,4,5) = glz(i)
         if(nl.eq.7.and.ke) then
            al(i,4,6) = -glz(i)
            ar(i,4,6) = -glz(i)
         end if

         if(nl.eq.7.and.ke) then
            al(i,5,1) = -gamma*capvell(i)*el(i) + 
     >           gamma1*capvell(i)*(2.d0*ql(i)+tkl(i))
            ar(i,5,1) = -gamma*capvelr(i)*er(i) + 
     >           gamma1*capvelr(i)*(2.d0*qr(i)+tkr(i))
         else
            al(i,5,1) = -gamma*capvell(i)*el(i) + 
     >           2.d0*gamma1*capvell(i)*ql(i)
            ar(i,5,1) = -gamma*capvelr(i)*er(i) + 
     >           2.d0*gamma1*capvelr(i)*qr(i)
         end if
         al(i,5,2) = -gamma1*capvell(i)*ul(i) + llx(i)*tenthalpyl(i)
         ar(i,5,2) = -gamma1*capvelr(i)*ur(i) + llx(i)*tenthalpyr(i)
         al(i,5,3) = -gamma1*capvell(i)*vl(i) + lly(i)*tenthalpyl(i)
         ar(i,5,3) = -gamma1*capvelr(i)*vr(i) + lly(i)*tenthalpyr(i)
         al(i,5,4) = -gamma1*capvell(i)*wl(i) + llz(i)*tenthalpyl(i)
         ar(i,5,4) = -gamma1*capvelr(i)*wr(i) + llz(i)*tenthalpyr(i)
         al(i,5,5) = gamma*capvell(i)
         ar(i,5,5) = gamma*capvelr(i)

         if(nl.eq.6) then
            al(i,6,1) = -tkl(i)*capvell(i)
            ar(i,6,1) = -tkr(i)*capvelr(i)
            al(i,6,2) = llx(i)*tkl(i)
            ar(i,6,2) = llx(i)*tkr(i)
            al(i,6,3) = lly(i)*tkl(i)
            ar(i,6,3) = lly(i)*tkr(i)
            al(i,6,4) = llz(i)*tkl(i)
            ar(i,6,4) = llz(i)*tkr(i)
            al(i,6,6) = capvell(i)
            ar(i,6,6) = capvelr(i)
         else if(nl.eq.7) then
            if(ke) then
               al(i,5,6) = -gamma1*capvell(i)
               ar(i,5,6) = -gamma1*capvelr(i)
            end if

            al(i,6,1) = -tkl(i)*capvell(i)
            ar(i,6,1) = -tkr(i)*capvelr(i)
            al(i,6,2) = llx(i)*tkl(i)
            ar(i,6,2) = llx(i)*tkr(i)
            al(i,6,3) = lly(i)*tkl(i)
            ar(i,6,3) = lly(i)*tkr(i)
            al(i,6,4) = llz(i)*tkl(i)
            ar(i,6,4) = llz(i)*tkr(i)
            al(i,6,6) = capvell(i)
            ar(i,6,6) = capvelr(i)

            al(i,7,1) = -twl(i)*capvell(i)
            ar(i,7,1) = -twr(i)*capvelr(i)
            al(i,7,2) = llx(i)*twl(i)
            ar(i,7,2) = llx(i)*twr(i)
            al(i,7,3) = lly(i)*twl(i)
            ar(i,7,3) = lly(i)*twr(i)
            al(i,7,4) = llz(i)*twl(i)
            ar(i,7,4) = llz(i)*twr(i)
            al(i,7,7) = capvell(i)
            ar(i,7,7) = capvelr(i)
         end if

         if(moving.ge.1) then
            if(nl.eq.7) then
               print *, '(lhs_matrix)',
     $              ' moving is not working well with kw yet'
               stop
            end if
            do n = 1, nl
               al(i,n,n) = al(i,n,n)+qtl(i)
               ar(i,n,n) = ar(i,n,n)+qtr(i)
            end do
         end if

      end do
c
      end if
c-------------------end if of precondition------------------------

c     
      IF (lhs_scheme.eq.2) Then ! use Zha scheme
         
         call zha_matrix(index,imax, dim2, nl,
     >        delta, gamma, ul, ur, vl, vr, wl, wr,rhol, rhor, 
     >        pl, pr,rhoel, rhoer,ql,qr,capul, capur, capvl, capvr,
     $        capwl, capwr, lx, ly, lz,  mx, my, mz, nx, ny, nz,
     >        al,ar,dhat_l,dhat_r, ilower)

c..   calculate  A_hat_left and A_hat_right

         do i = 1, imax
            do ii= 1, nl
               do jj= 1, nl
                  al(i,ii,jj) = 0.5*(al(i,ii,jj)+dhat_l(i,ii,jj))
                  ar(i,ii,jj) = 0.5*(ar(i,ii,jj)-dhat_r(i,ii,jj))
               end do
            end do
         end do

      else if (lhs_scheme.eq.4) then ! use van leer  scheme

         do i = 1, imax
            
            area_ = dsqrt(llx(i)**2 + lly(i)**2 + llz(i)**2)
            vx1  = llx(i)/area_
            vy1  = lly(i)/area_
            vz1  = llz(i)/area_

            c_l = dsqrt(gamma*pl(i)/rhol(i))
            d_l  = rhol(i)
            u_l  = rhoul(i)/rhol(i)
            v_l  = rhovl(i)/rhol(i)
            w_l  = rhowl(i)/rhol(i)
            e_l  = rhoel(i)

            c_r = dsqrt(gamma*pr(i)/rhor(i))
            d_r  = rhor(i)
            u_r  = rhour(i)/rhor(i)
            v_r  = rhovr(i)/rhor(i)
            w_r  = rhowr(i)/rhor(i)
            e_r  = rhoer(i)
            if (nl.eq.6) then
              tk_l  = rhotkl(i)/rhol(i)
              tk_r  = rhotkr(i)/rhor(i)
            end if
c
            if ( moving.eq.1 ) then
              qt_l = qtl(i)/area_
              qt_r = qtr(i)/area_
            else
              qt_l = 0.d0
              qt_r = 0.d0
            end if
c
            call van_leer_matrix(gamma,nl,vx1,vy1,vz1,c_l,d_l,u_l,v_l,
     >           w_l,e_l,c_r,d_r,u_r,v_r,w_r,e_r,tk_l,tk_r,
     >           aml,amr,qt_l,qt_r) 

c..   calculate  a_hat_left and a_hat_right

            do ii= 1,nl
               do jj= 1,nl
                  al(i,ii,jj) = area_*aml(ii,jj)
                  ar(i,ii,jj) = area_*amr(ii,jj)
               end do
            end do
            
         END DO
         

      ELSE IF (lhs_scheme.eq.1) Then ! use Roe scheme
c     
c     Roe matrix
c     
         if (precondition.ge.1) then
            call roe_matrix_p(index,imax, nl,
     $        delta, sqrtinv, sqrtrhol, sqrtrhor,
     $        ul, ur, vl, vr, wl, wr, tenthalpyl, tenthalpyr,
     $        lx, ly, lz, lxhat, lyhat, lzhat, mx, my, mz, nx, ny, nz,
     $        mxhat, myhat, mzhat, nxhat, nyhat, nzhat, tkl, tkr,
     $        twl, twr, a, ke, qtl, qtr,ilower,
     $        rhol,rhor,tl,tr,pl,pr,vist,volumeinv,
     $        diver, control)
         else
           call roe_matrix(index, imax, nl, dim2, delta, sqrtinv,
     $     sqrtrhol, sqrtrhor, ul, ur, vl, vr, wl, wr, tenthalpyl,
     $     tenthalpyr, lx, ly, lz, lxhat, lyhat, lzhat, mx, my, mz, nx,
     $     ny, nz, mxhat, myhat, mzhat, nxhat, nyhat, nzhat, tkl, tkr,
     $     twl, twr, a, ke, qtl, qtr, ilower, qtline, control)
         end if


c     
c     compute A_hat_left and A_hat_right
c     
         do i = 1,imax
            al(i,:,:) = .5d0* (al(i,:,:) + a(:,:,i))
            ar(i,:,:) = .5d0* (ar(i,:,:) - a(:,:,i))
         end do

      end if

c     diffusion matrices
      ll = 0.d0
      lr = 0.d0
c----------------redefine the rho and rhoe for matrix_bnd
      if (precondition.ge.2) then
         do i=0,imax
            call qstat(p(i),t(i),diver(i),temp)
            rho(i)=1.0d0/temp
            call enthalpy(tenthalpyr(i),t(i),temp,p(i))
            tenthalpyr(i)=tenthalpyr(i)+qq(i)
  	    rhoe(i)=rho(i)*tenthalpyr(i)-p(i)
         end do
      end if
c------------------------
      if (.not.inviscid) then
         do i = 0, imax
            if  ( suther )  then
                molec(i) = dsqrt(t(i))**3 * 
     $                     (1.d0 + tref) / (t(i) + tref)
            else
                molec(i) = t(i)
            end if
         end do
c--------------------------
      if (precondition.ge.1) then
         if (precondition.ge.2) then
            do i=0,imax
               temp=1.0d0/rho(i)
               call cpt(t(i),cpr(i))
               call drhdtpp(temp,t(i),drhodtr(i))
               call drhdptt(temp,t(i),drhodpr(i))
c------------same as 1
             end do
          else
            do i=0,imax
               tenthalpyr(i)=(rhoe(i)+p(i))/rho(i)
               c=sqrt((gamma-1.0)*(tenthalpyr(i)-qq(i)))
               drhodpr(i)=gamma/c**2
               drhodtr(i)=-rho(i)/t(i)
               cpr(i)=c**2/((gamma-1.0)*t(i))
               cpr(i)=1.d0/((gamma-1.0)*machinf**2)
             end do
          end if
       end if
c---------------------------------

      do i = 1, imax-1
         volumeinv(i) = 1.d0/volume(i)
      end do

      i = 0
      if(bclower.eq.7) then
         volumeinv(i) = 1.d0/volume(i)
      else if(bclower.eq.10) then
         volumeinv(i) = 1.d0/volume(imax-1)
      else
         volumeinv(i) = volumeinv(i+1)
      end if

      i = imax
      if(bcupper.eq.7)then
         volumeinv(i) = 1.d0/volume(i)
      else if(bcupper.eq.10) then
         volumeinv(i) = 1.d0/volume(1)
      else
         volumeinv(i) = volumeinv(i-1)
      end if

      coeff1 = 0.5d0 / reynolds

      if (nl.eq.6) then
        do i = 1,imax
          viscous(i) = coeff1*(  ! alpha_tau
     $        (molec(i-1)+vist(i-1))*volumeinv(i-1)
     $        +(molec(i)+vist(i))*volumeinv(i))      
        end do
      else
        do i = 1,imax
          viscous(i) = coeff1*(  ! alpha_tau
     $        (molec(i-1)+vist(i-1)*reynolds)*volumeinv(i-1)
     $        +(molec(i)+vist(i)*reynolds)*volumeinv(i))      
        end do
      end if

      do i = 1,imax
        term1(i) = viscous(i)*(1.333333d0*llx(i)*llx(i) + 
     >        lly(i)*lly(i) + llz(i)*llz(i))  
        term2(i) = viscous(i)*0.333333d0*llx(i)*lly(i)  
        term3(i) = viscous(i)*0.333333d0*llx(i)*llz(i)
      end do

      do   i = 1,imax
         ll(i,2,1) =  term1(i)*u(i-1)*rinv(i-1) +
     >        term2(i)*v(i-1)*rinv(i-1) + 
     >        term3(i)*w(i-1)*rinv(i-1) 
         lr(i,2,1) = -(term1(i)*u(i)*rinv(i) +
     >        term2(i)*v(i)*rinv(i) +
     >        term3(i)*w(i)*rinv(i) )
         ll(i,2,2) = -term1(i)*rinv(i-1)
         lr(i,2,2) =  term1(i)*rinv(i)
         ll(i,2,3) = -term2(i)*rinv(i-1)
         lr(i,2,3) =  term2(i)*rinv(i)
         ll(i,2,4) = -term3(i)*rinv(i-1) 
         lr(i,2,4) =  term3(i)*rinv(i) 
         if(nl.eq.7.and.ke) then
            ll(i,2,6) = -llx(i)/3.d0 
            lr(i,2,6) = -llx(i)/3.d0 
         end if
      end do

      do i = 1,imax
         term1(i) = viscous(i)*0.333333d0*llx(i)*lly(i)
         term2(i) = viscous(i)*(llx(i)*llx(i) + 
     >        1.333333d0*lly(i)*lly(i) + llz(i)*llz(i))  
         term3(i) = viscous(i)*0.333333d0*lly(i)*llz(i)
      end do

      do i = 1,imax
         ll(i,3,1)  =  term1(i)*u(i-1)*rinv(i-1) +
     >        term2(i)*v(i-1)*rinv(i-1) + 
     >        term3(i)*w(i-1)*rinv(i-1) 
         lr(i,3,1)  = -(term1(i)*u(i)*rinv(i) +
     >        term2(i)*v(i)*rinv(i) +
     >        term3(i)*w(i)*rinv(i) )
         ll(i,3,2)  = -term1(i)*rinv(i-1)
         lr(i,3,2)  =  term1(i)*rinv(i)
         ll(i,3,3)  = -term2(i)*rinv(i-1)
         lr(i,3,3)  =  term2(i)*rinv(i)
         ll(i,3,4)  = -term3(i)*rinv(i-1) 
         lr(i,3,4)  =  term3(i)*rinv(i) 
         if(nl.eq.7.and.ke) then
            ll(i,3,6) = -lly(i)/3.d0
            lr(i,3,6) = -lly(i)/3.d0
         end if
      end do
c     
      do i = 1,imax
         term1(i) = viscous(i)*0.333333d0*llx(i)*llz(i)
         term2(i) = viscous(i)*0.333333d0*lly(i)*llz(i)
         term3(i) = viscous(i)*(llx(i)*llx(i) + lly(i)*lly(i) + 
     >        1.333333d0*llz(i)*llz(i))  
      end do

      do  i = 1,imax
         ll(i,4,1)  =  term1(i)*u(i-1)*rinv(i-1) +
     >        term2(i)*v(i-1)*rinv(i-1) +
     >        term3(i)*w(i-1)*rinv(i-1)
         lr(i,4,1)  = -(term1(i)*u(i)*rinv(i) +
     >        term2(i)*v(i)*rinv(i) +
     >        term3(i)*w(i)*rinv(i) )
         ll(i,4,2)  = -term1(i)*rinv(i-1)
         lr(i,4,2)  =  term1(i)*rinv(i)
         ll(i,4,3)  = -term2(i)*rinv(i-1)
         lr(i,4,3)  =  term2(i)*rinv(i)
         ll(i,4,4)  = -term3(i)*rinv(i-1)
         lr(i,4,4)  =  term3(i)*rinv(i)
         if(nl.eq.7.and.ke) then
            ll(i,4,6) = -llz(i)/3.d0
            lr(i,4,6) = -llz(i)/3.d0
         end if
      end do

      if (nl.eq.6) then
        do i = 1, imax
          alphak(i) = coeff1*((molec(i-1)+rhotk(i-1))*volumeinv(i-1)
     $               +(molec(i)+rhotk(i))*volumeinv(i))/tko
        end do
      else if (nl.eq.7) then
        do i = 1, imax
          alphak(i) = coeff1*((molec(i-1)+vist(i-1)*reynolds*tks)
     $        *volumeinv(i-1)
     $        +(molec(i)+vist(i)*reynolds*tks)
     $        * volumeinv(i))
          alphaw(i) = coeff1*((molec(i-1)+vist(i-1)*reynolds*tws)
     $        *volumeinv(i-1)
     $        +(molec(i)+vist(i)*reynolds*tws)
     $        * volumeinv(i))
        end do
      end if

      if (nl.eq.6) then
        do i = 1,imax
          viscous(i) = coeff1*(  
     $        (molec(i-1)/prandtl+vist(i-1)/prt)
     $        *volumeinv(i-1)
     $        +(molec(i)/prandtl+vist(i)/prt)
     $        * volumeinv(i))      
        end do
      else
        do i = 1,imax
          viscous(i) = coeff1*(  
     $        (molec(i-1)/prandtl+vist(i-1)*reynolds/prt)
     $        *volumeinv(i-1)
     $        +(molec(i)/prandtl+vist(i)*reynolds/prt)
     $        * volumeinv(i))      
        end do
      end if

      coeff2 = 1.d0/(gamma1*machinf*machinf) ! Cp
      if (precondition.ge.2) coeff2=cpr(i)   !cpr i.e. Cp
      coeff3 = gamma*gamma1*machinf*machinf 

      do i = 1,imax
         thermal(i) = coeff2*viscous(i) ! alpha_q
         term1(i) = llx(i)*llx(i) + lly(i)*lly(i) + llz(i)*llz(i)
         term2(i) = thermal(i)*term1(i)*coeff3
      end do

      do i = 1,imax
         if(nl.eq.7.and.ke) then
            ll(i,5,1) = term2(i)*rinv(i-1)*
     >           (rhoe(i-1)*rinv(i-1) - (u(i-1)**2 +
     >           v(i-1)**2 + w(i-1)**2)-tk(i-1))
         else
            ll(i,5,1) = term2(i)*rinv(i-1)*
     >           (rhoe(i-1)*rinv(i-1) - (u(i-1)**2 +
     >           v(i-1)**2 + w(i-1)**2))
         end if
         ll(i,5,1) = ll(i,5,1) + 
     >        0.5d0*((u(i)+u(i-1))*ll(i,2,1) -
     >        u(i-1)*rinv(i-1)*fflux(i,2) +
     >        (v(i)+v(i-1))*ll(i,3,1) - 
     >        v(i-1)*rinv(i-1)*fflux(i,3) +
     >        (w(i)+w(i-1))*ll(i,4,1) - 
     >        w(i-1)*rinv(i-1)*fflux(i,4) )
         if(nl.eq.7.and.ke) ll(i,5,1) = ll(i,5,1) +
     $        alphak(i)*tk(i-1)*rinv(i-1)*term1(i)
         if(nl.eq.7.and.ke) then
            lr(i,5,1) = -term2(i)*rinv(i)*
     >           (rhoe(i)*rinv(i) - (u(i)**2 +
     >           v(i)**2 + w(i)**2)-tk(i))
         else
            lr(i,5,1) = -term2(i)*rinv(i)*
     >           (rhoe(i)*rinv(i) - (u(i)**2 +
     >           v(i)**2 + w(i)**2))
         end if
         lr(i,5,1) = lr(i,5,1) + 
     >        0.5d0*((u(i)+u(i-1))*lr(i,2,1) -
     >        u(i)*rinv(i)*fflux(i,2) +
     >        (v(i)+v(i-1))*lr(i,3,1) - 
     >        v(i)*rinv(i)*fflux(i,3) +
     >        (w(i)+w(i-1))*lr(i,4,1) - 
     >        w(i)*rinv(i)*fflux(i,4) )
         if(nl.eq.7.and.ke) lr(i,5,1) = lr(i,5,1) -
     $        alphak(i)*tk(i)*rinv(i)*term1(i)
         ll(i,5,2) = term2(i)*u(i-1)*rinv(i-1) +
     >        0.5d0*((u(i)+u(i-1))*ll(i,2,2) +
     >        rinv(i-1)*fflux(i,2) +
     >        (v(i)+v(i-1))*ll(i,3,2) +
     >        (w(i)+w(i-1))*ll(i,4,2) )    
         lr(i,5,2) = -term2(i)*u(i)*rinv(i) +
     >        0.5d0*((u(i)+u(i-1))*lr(i,2,2) +
     >        rinv(i)*fflux(i,2) +
     >        (v(i)+v(i-1))*lr(i,3,2) +
     >        (w(i)+w(i-1))*lr(i,4,2) )    
         ll(i,5,3) = term2(i)*v(i-1)*rinv(i-1) +
     >        0.5d0*((u(i)+u(i-1))*ll(i,2,3) +
     >        (v(i)+v(i-1))*ll(i,3,3) +
     >        rinv(i-1)*fflux(i,3) +
     >        (w(i)+w(i-1))*ll(i,4,3) )    
         lr(i,5,3) = -term2(i)*v(i)*rinv(i) +
     >        0.5d0*((u(i)+u(i-1))*lr(i,2,3) +
     >        (v(i)+v(i-1))*lr(i,3,3) +
     >        rinv(i)*fflux(i,3) +
     >        (w(i)+w(i-1))*lr(i,4,3) )    
         ll(i,5,4) = term2(i)*w(i-1)*rinv(i-1) +
     >        0.5d0*((u(i)+u(i-1))*ll(i,2,4) +
     >        (v(i)+v(i-1))*ll(i,3,4) +
     >        (w(i)+w(i-1))*ll(i,4,4) +    
     >        rinv(i-1)*fflux(i,4) )
         lr(i,5,4) = -term2(i)*w(i)*rinv(i) +
     >        0.5d0*((u(i)+u(i-1))*lr(i,2,4) +
     >        (v(i)+v(i-1))*lr(i,3,4) +
     >        (w(i)+w(i-1))*lr(i,4,4) +
     >        rinv(i)*fflux(i,4) )
         ll(i,5,5) = -term2(i)*rinv(i-1) +
     >        0.5d0*((u(i)+u(i-1))*ll(i,2,5) +
     >        (v(i)+v(i-1))*ll(i,3,5) +
     >        (w(i)+w(i-1))*ll(i,4,5) )    
         lr(i,5,5) =  term2(i)*rinv(i) +
     >        0.5d0*((u(i)+u(i-1))*lr(i,2,5) +
     >        (v(i)+v(i-1))*lr(i,3,5) +
     >        (w(i)+w(i-1))*lr(i,4,5) )

         if(nl.eq.7) then
            if(ke) then
               ll(i,5,6) = .5d0*((u(i)+u(i-1))*ll(i,2,6)
     $              +(v(i)+v(i-1))*ll(i,3,6) +(w(i)+w(i-1))*ll(i,4,6))
               ll(i,5,6) = ll(i,5,6)
     $              +term2(i)*rinv(i-1) - alphak(i)*term1(i)*rinv(i-1)
               lr(i,5,6) = .5d0*((u(i)+u(i-1))*lr(i,2,6)
     $           +(v(i)+v(i-1))*lr(i,3,6) +(w(i)+w(i-1))*lr(i,4,6))
               lr(i,5,6) = ll(i,5,6)
     $           -term2(i)*rinv(i) + alphak(i)*term1(i)*rinv(i)
            end if

            ll(i,6,1) = alphak(i)*term1(i)*tk(i-1)*rinv(i-1)
            lr(i,6,1) = -alphak(i)*term1(i)*tk(i)*rinv(i)
            ll(i,6,6) = -alphak(i)*term1(i)*rinv(i-1)
            lr(i,6,6) = alphak(i)*term1(i)*rinv(i)
            
            ll(i,7,1) =  alphaw(i)*term1(i)*tw(i-1)*rinv(i-1)
            lr(i,7,1) = -alphaw(i)*term1(i)*tw(i)*  rinv(i)
            ll(i,7,7) = -alphaw(i)*term1(i)*rinv(i-1)
            lr(i,7,7) =  alphaw(i)*term1(i)*rinv(i)
         else if (nl.eq.6) then
           ll(i,6,1) = alphak(i)*term1(i)*tk(i-1)*rinv(i-1)
           lr(i,6,1) = -alphak(i)*term1(i)*tk(i)*rinv(i)
           ll(i,6,6) = -alphak(i)*term1(i)*rinv(i-1)
           lr(i,6,6) = alphak(i)*term1(i)*rinv(i)
         end if
c-------------------------------------
         if (precondition.ge.1) then
            call dwdq(nl,cpr(i-1),drhodpr(i-1),drhodtr(i-1),
     $                tenthalpyr(i-1),
     $                rho(i-1),u(i-1),v(i-1),w(i-1),k_matrix)
            ll(i,:,:)=matmul(ll(i,:,:),k_matrix(:,:))
            call dwdq(nl,cpr(i),drhodpr(i),drhodtr(i),tenthalpyr(i),
     $                rho(i),u(i),v(i),w(i),k_matrix)
            lr(i,:,:)=matmul(lr(i,:,:),k_matrix(:,:))
         end if

      end do

      end if
c---------------------end if of inviscid-----------------
c     
c     compute block matrices
c     
      capa = 0.d0
      capb = 0.d0
      capc = 0.d0

      do i = 1,imax
         dtv(i) = dtv(i) / volume(i)
      end do

      do ii = 1,nl
         do jj = 1,nl  
            do i = 1,imax-1
               capa(ii,jj,i) = dtv(i) * 
     >              ((al(i+1,ii,jj) - ll(i+1,ii,jj)) -
     >              (ar(i, ii,jj) - lr(i,  ii,jj)) )
               capb(ii,jj,i) = -dtv(i) * 
     >              ( al(i,ii,jj) - ll(i,ii,jj) )
               capc(ii,jj,i) =  dtv(i) * 
     >              (ar(i+1,ii,jj) - lr(i+1,ii,jj) )
            end do
         end do
      end do

c     ZERO OUT gammabc before defining it

      call matrix_bnd(index, i1, i2, il, jl, kl, x, y, z,
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     nl, dim1, dim2, imax-1, bclower, bcupper,
     $     qq, p, u, v, w, t, capa, capb, capc, index, 
     $     rho, rhoe, 
     $     rhotk, rhotw, wallfunc, ke, qconer,
     $     vbds, control)
c
c     viscous terms in xi direction for non-reflecting bc
c
      if ( index.eq.1 ) then
        do 1000  ii = 1,nl
        do 1000  jj = 1,nl
          dtv(1) = dt(1,i1,i2) / vol(1,i1,i2)
          ci_nr(ii,jj,1,1) =  dtv(1) * ( lr(1,ii,jj) - ll(2,ii,jj) )
          ci_nr(ii,jj,2,1) =  dtv(1) * ll(1,ii,jj)
          ci_nr(ii,jj,3,1) = -dtv(1) * lr(2,ii,jj)
1000    continue
c
        do 2000  ii = 1,nl
        do 2000  jj = 1,nl
          dtv(il) = dt(il,i1,i2) / vol(il,i1,i2)
          ci_nr(ii,jj,1,2) =  dtv(il) * ( lr(il,ii,jj) - ll(il+1,ii,jj))
          ci_nr(ii,jj,2,2) =  dtv(il) * ll(il,ii,jj)
          ci_nr(ii,jj,3,2) = -dtv(il) * lr(il+1,ii,jj)
2000    continue
      end if
c
c     viscous terms in ate direction for non-reflecting bc
c
      if ( index.eq.2 ) then
        do 3000  ii = 1,nl
        do 3000  jj = 1,nl
          dtv(1) = dt(i1,1,i2) / vol(i1,1,i2)
          ci_nr(ii,jj,1,1) =  dtv(1) * ( lr(1,ii,jj) - ll(2,ii,jj) )
          ci_nr(ii,jj,2,1) =  dtv(1) * ll(1,ii,jj)
          ci_nr(ii,jj,3,1) = -dtv(1) * lr(2,ii,jj)
3000    continue
c
        do 4000  ii = 1,nl
        do 4000  jj = 1,nl
          dtv(jl) = dt(i1,jl,i2) / vol(i1,jl,i2)
          ci_nr(ii,jj,1,2) =  dtv(jl) * ( lr(jl,ii,jj) - ll(jl+1,ii,jj))
          ci_nr(ii,jj,2,2) =  dtv(jl) * ll(jl,ii,jj)
          ci_nr(ii,jj,3,2) = -dtv(jl) * lr(jl+1,ii,jj)
4000    continue
      end if
c
c     viscous terms in zeta direction for non-reflecting bc
c
      if ( index.eq.3 ) then
        do 5000  ii = 1,nl
        do 5000  jj = 1,nl
          dtv(1) = dt(i1,i2,1) / vol(i1,i2,1)
          ci_nr(ii,jj,1,1) =  dtv(1) * ( lr(1,ii,jj) - ll(2,ii,jj) )
          ci_nr(ii,jj,2,1) =  dtv(1) * ll(1,ii,jj)
          ci_nr(ii,jj,3,1) = -dtv(1) * lr(2,ii,jj)
5000    continue
c
        do 6000  ii = 1,nl
        do 6000  jj = 1,nl
          dtv(kl) = dt(i1,i2,kl) / vol(i1,i2,kl)
          ci_nr(ii,jj,1,2) =  dtv(kl) * ( lr(kl,ii,jj) - ll(kl+1,ii,jj))
          ci_nr(ii,jj,2,2) =  dtv(kl) * ll(kl,ii,jj)
          ci_nr(ii,jj,3,2) = -dtv(kl) * lr(kl+1,ii,jj)
6000    continue
      end if
c
      return
      end
