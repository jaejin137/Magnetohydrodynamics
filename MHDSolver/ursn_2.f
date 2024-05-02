      subroutine ursn_2(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux,
     $     bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $     rhs, il, jl,
     $     kl, nl, dim1, dim2, ilower, iupper, jlower, jupper, klower,
     $     kupper, visturb, dstep, checktime, wallfunc, ke,
     $     qt, sr, udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $     wdk, cst_nri, ced_nri, cst_nrj, ced_nrj, tko, qiii, diver, 
     $     control)
c     
c     use Upwind Relaxation Sweeping(urs)  with Gauss-Seidel
c     to invert the implicit LHS, renamed from original urs.f
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
     $     il, jl, kl, nl, dim2, dim1,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     dstep

      integer::
     $     dual_t,
     $     moving, gcl

      double precision, intent(in)::
     $     rhs(il,jl,kl,nl),
     $     sr(il,jl,kl)

      double precision, dimension(jl,kl,4), intent(in)::
     $     udi, vdi, wdi

      double precision, dimension(il,kl,4), intent(in)::
     $     udj, vdj, wdj

      double precision, dimension(il,jl,4), intent(in)::
     $     udk, vdk, wdk

      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z

      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol

      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q,qiii

      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, 3)::
     $     qt

      double precision, dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1), intent(in)::
     $     visturb

      double precision, intent(in)::
     $     dt(il,jl,kl)
      double precision::
     $     gamma, prandtl, prt, tref, ptotal, ttotal,
     $     angl1, angl2, poutlet

      integer, intent(in)::
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl)
      
      double precision,intent(in)::
     $     diver(ilower:iupper,jlower:jupper,klower:kupper)

      double precision, intent(in)::
     $     delta(5)

      integer::
     $     rhs_order, limiter, rhs_scheme, idimen,
     $     lhs_order, iter_gs, lhs_scheme
      
      double precision::
     $     machinf, epsfactor, kfactor, reynolds,
     $     tintvl,
     $     ronum

      logical, intent(in)::
     $     checktime,           ! if systime check necessary
     $     wallfunc,
     $     ke
      logical::
     $     inviscid, suther

      double precision, dimension(il+1, jl+1, kl+1, nl), intent(in)::
     $     rflux, sflux, tflux

c
c     SA 1eq turbulent model contants
c
      double precision, intent(in):: tko

      double precision, intent(out), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q_n

      double precision, intent(out)::
     $     cst_nri(nl,nl,jl,kl,5), ced_nri(nl,nl,jl,kl,5),
     $     cst_nrj(nl,nl,il,kl,5), ced_nrj(nl,nl,il,kl,5)
c     
c     LOCAL VARIABLES
c     
      double precision, dimension(nl,nl,dim1)::
     $     capa, capb, capc,
     $     f(nl,dim1),
     $     soln(nl,dim1)

      double precision::
     $     qtline(-1:dim2,3)

      double precision, dimension(nl,nl,3,2)::
     $     ci_nr

c..   coefficient matrices

      double precision, dimension(nl,nl,il,jl,kl)::
     $     a_minus, a_plus,
     $     b_minus, b_plus,
     $     c_minus, c_plus,
     $     b_bar

      integer::
     $     n,
     $     ii, jj,
     $     i, j, k,
     $     nblocks,
     $     is,js,
     $     m,
     $     nj,
     $     bclower, bcupper

      double precision::
     $     cx1,cx2,
     $     cy1,cy2,
     $     cz1,cz2

      integer, parameter::
     $     bcperi = 10,
     $     bcwall = 3


      double precision QB(NL,NL),FIRST(NL),DQ(NL)     
      double precision, dimension(nl,nl,il,jl,kl)::
     $     IN_A
c--------------------------parameter transfer
      dual_t=control%dual_t
      gcl=control%gcl
      idimen=control%idimen
      iter_gs=control%iter_gs
      lhs_order=control%lhs_order
      lhs_scheme=control%lhs_scheme
      limiter=control%limiter
      moving=control%moving
      rhs_order=control%rhs_order
      rhs_scheme=control%rhs_scheme

      inviscid=control%inviscid
      suther=control%suther

      angl1=control%angl1
      angl2=control%angl2
      epsfactor=control%epsfactor
      gamma=control%gamma
      kfactor=control%kfactor
      machinf=control%machinf
      poutlet=control%poutlet
      ptotal=control%ptotal
      prandtl=control%prandtl
      prt=control%prt
      reynolds=control%reynolds
      ronum=control%ronum
      tref=control%tref
      tintvl=control%tintvl
      ttotal=control%ttotal
c--------------------------------------------------------
c ... the array q_n is used for dq in this subroutine
c
c *** SUBROUTINE START ***
c     
      q_n = 0.d0
c$$$      a_minus = 0.d0
c$$$      a_plus = 0.d0
c$$$      b_minus = 0.d0
c$$$      b_plus = 0.d0
c$$$      c_minus = 0.d0
c$$$      c_plus = 0.d0

c ... calculate the matrices in xi-direction

      do k = 1, kl
        do j = 1,jl
c
          bclower = bc_xi_lower(j,k)
          bcupper = bc_xi_upper(j,k)
c
          if(moving.ge.1.or.dabs(ronum).gt.1e-9)
     $           qtline(ilower:iupper,:) = qt(ilower:iupper,j,k,:)
c

          call lhs_matrix(1, j, k, dt, x, y, z, vol, bclower, bcupper,
     $           delta, wallfunc, capa,
     $           capb, capc, q, rflux, il, jl, kl, nl, dim2, dim1,
     $           ilower, iupper, jlower, jupper, klower, kupper,
     $           visturb, ke, udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $           wdk, qtline, ci_nr, tko, qiii, diver, control)
c
          a_minus(:,:,:,j,k) = capb(:,:,:il)
          a_plus(:,:,:,j,k)  = capc(:,:,:il)
          b_bar(:,:,:,j,k)   = capa(:,:,:il)

          do ii = 1,nl
             do jj = 1,nl
                cst_nri(ii,jj,j,k,1) = ci_nr(ii,jj,1,1) - capa(ii,jj, 1)
                ced_nri(ii,jj,j,k,1) = ci_nr(ii,jj,1,2) - capa(ii,jj,il)
             end do
          end do
c
        end do
      end do
c
c ... nrbc at j boundaries
c
      if ( bc_eta_lower(1,1).ge.14.and.bc_eta_lower(1,1).le.17) then
         do k  = 1,kl
            do i  = 1,il
               do ii = 1,nl
                  do jj = 1,nl
                     cst_nrj(ii,jj,i,k,2) = a_minus(ii,jj,i,1,k)
                     cst_nrj(ii,jj,i,k,3) =  a_plus(ii,jj,i,1,k)
                     cst_nrj(ii,jj,i,k,1) =   b_bar(ii,jj,i,1,k)
                  end do
               end do
            end do
         end do
      end if
c
      if ( bc_eta_upper(1,1).ge.14.and.bc_eta_upper(1,1).le.17) then
         do k  = 1,kl
            do i  = 1,il
               do ii = 1,nl
                  do jj = 1,nl
                     ced_nrj(ii,jj,i,k,2) = a_minus(ii,jj,i,jl,k)
                     ced_nrj(ii,jj,i,k,3) =  a_plus(ii,jj,i,jl,k)
                     ced_nrj(ii,jj,i,k,1) =   b_bar(ii,jj,i,jl,k)
                  end do
               end do
            end do
         end do
      end if

c
c ... calculate the matrices in eta-direction
c
      do k = 1, kl
        do i = 1,il
c
          bclower = bc_eta_lower(i,k)
          bcupper = bc_eta_upper(i,k)
          if(moving.ge.1.or.dabs(ronum).gt.1e-9)
     $           qtline(jlower:jupper,:) = qt(i,jlower:jupper,k,:)
c
          call lhs_matrix(2, i, k, dt, x, y, z, vol, bclower, bcupper,
     $           delta, wallfunc, capa,
     $           capb, capc, q, sflux, il, jl, kl, nl, dim2, dim1,
     $           ilower, iupper, jlower, jupper, klower, kupper,
     $           visturb, ke, udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $           wdk, qtline, ci_nr, tko,  qiii, diver, control)
c
          b_minus(:,:,i,:,k) = capb(:,:,:jl)
          b_plus(:,:,i,:,k)  = capc(:,:,:jl)
          b_bar(:,:,i,:,k)   = b_bar(:,:,i,:,k) + capa(:,:,:jl)
c
        do ii = 1,nl
           do jj = 1,nl
              cst_nrj(ii,jj,i,k,1) = cst_nrj(ii,jj,i,k,1)
     $             + ci_nr(ii,jj,1,1)
              ced_nrj(ii,jj,i,k,1) = ced_nrj(ii,jj,i,k,1)
     $             + ci_nr(ii,jj,1,2)
           end do
        end do
c
        end do
      end do
c
c ... nrbc at i boundaries
c
      if ( bc_xi_lower(1,1).ge.14.and.bc_xi_lower(1,1).le.17) then
         do k  = 1,kl
            do j  = 1,jl
               do ii = 1,nl
                  do jj = 1,nl
                     cst_nri(ii,jj,j,k,2) = b_minus(ii,jj,1,j,k)
                     cst_nri(ii,jj,j,k,3) =  b_plus(ii,jj,1,j,k)
                     cst_nri(ii,jj,j,k,1) =   b_bar(ii,jj,1,j,k) +
     >                    cst_nri(ii,jj,j,k,1)
                  end do
               end do
            end do
         end do
      end if
c
      if ( bc_xi_upper(1,1).ge.14.and.bc_xi_upper(1,1).le.17) then
         do k  = 1,kl
            do j  = 1,jl
               do ii = 1,nl
                  do jj = 1,nl
                     ced_nri(ii,jj,j,k,2) = b_minus(ii,jj,il,j,k)
                     ced_nri(ii,jj,j,k,3) =  b_plus(ii,jj,il,j,k)
                     ced_nri(ii,jj,j,k,1) =   b_bar(ii,jj,il,j,k) +
     >                    ced_nri(ii,jj,j,k,1)
                  end do
               end do
            end do
         end do
      end if
c
      if ( idimen.gt.2 ) then
c
c ... calculate the matrices in zeta-direction
c
        do j = 1,jl
          do i = 1,il
c
            bclower = bc_zeta_lower(i,j)
            bcupper = bc_zeta_upper(i,j)
c
            if(moving.ge.1.or.dabs(ronum).gt.1e-9)
     $              qtline(klower:kupper,:) = qt(i,j,klower:kupper,:)
c
            call lhs_matrix(3, i, j, dt, x, y, z, vol, bclower,
     $        bcupper, delta, wallfunc, 
     $        capa, capb, capc, q, tflux, il, jl, kl, nl,
     $        dim2, dim1, ilower, iupper, jlower, jupper, klower,
     $        kupper, visturb, ke, udi, vdi, wdi, udj, vdj, wdj,
     $        udk, vdk, wdk, qtline, ci_nr, tko, qiii, diver,
     $        control)
c
            c_minus(:,:,i,j,:) = capb(:,:,:kl)
            c_plus(:,:,i,j,:) = capc(:,:,:kl)
            b_bar(:,:,i,j,:) = b_bar(:,:,i,j,:) + capa(:,:,:kl)
c
          end do
        end do
c
c ... nrbc at i boundaries
c
        if ( bc_xi_lower(1,1).ge.14.and.bc_xi_lower(1,1).le.17 ) then
           do k  = 1,kl
              do j  = 1,jl
                 do ii = 1,nl
                    do jj = 1,nl
                       cst_nri(ii,jj,j,k,4) = c_minus(ii,jj,1,j,k)
                       cst_nri(ii,jj,j,k,5) =  c_plus(ii,jj,1,j,k)
                       cst_nri(ii,jj,j,k,1) =   b_bar(ii,jj,1,j,k) +
     >                      cst_nri(ii,jj,j,k,1)
                    end do
                 end do
              end do
           end do
        end if
c
        if ( bc_xi_upper(1,1).ge.14.and.bc_xi_upper(1,1).le.17) then
           do k  = 1,kl
              do j  = 1,jl
                 do ii = 1,nl
                    do jj = 1,nl
                       ced_nri(ii,jj,j,k,4) = c_minus(ii,jj,il,j,k)
                       ced_nri(ii,jj,j,k,5) =  c_plus(ii,jj,il,j,k)
                       ced_nri(ii,jj,j,k,1) =   b_bar(ii,jj,il,j,k) +
     >                      ced_nri(ii,jj,j,k,1)
                    end do
                 end do
              end do
           end do
        end if
c
c ... nrbc at j boundaries
c
        if ( bc_eta_lower(1,1).ge.14.and.bc_eta_lower(1,1).le.17) then
           do k  = 1,kl
              do i  = 1,il
                 do ii = 1,nl
                    do jj = 1,nl
                       cst_nrj(ii,jj,i,k,4) = c_minus(ii,jj,i,1,k)
                       cst_nrj(ii,jj,i,k,5) =  c_plus(ii,jj,i,1,k)
                       cst_nrj(ii,jj,i,k,1) =   b_bar(ii,jj,i,1,k) +
     >                      cst_nrj(ii,jj,i,k,1)
                    end do
                 end do
              end do
           end do
        end if
c
        if ( bc_eta_upper(1,1).ge.14.and.bc_eta_upper(1,1).le.17) then
           do k  = 1,kl
              do i  = 1,il
                 do ii = 1,nl
                    do jj = 1,nl
                       ced_nrj(ii,jj,i,k,4) = c_minus(ii,jj,i,jl,k)
                       ced_nrj(ii,jj,i,k,5) =  c_plus(ii,jj,i,jl,k)
                       ced_nrj(ii,jj,i,k,1) =   b_bar(ii,jj,i,jl,k) +
     >                      ced_nrj(ii,jj,i,k,1)
                    end do
                 end do
              end do
           end do
         
        end if
c
      end if
c
c ... k-omega wall function special treatment
c
      if(nl.eq.7) then
         do k = 1, kl
            do j = 1, jl
               do i = 1, il, il-1
                  if(i.eq.1.and.bc_xi_lower(j,k).eq.bcwall.or.
     $                 i.eq.il.and.bc_xi_upper(j,k).eq.bcwall) then
                     if(wallfunc) then
                        a_plus(6,:,i,j,k) = 0.d0
                        a_minus(6,:,i,j,k) = 0.d0
                        b_plus(6,:,i,j,k) = 0.d0
                        b_minus(6,:,i,j,k) = 0.d0
                        c_plus(6,:,i,j,k) = 0.d0
                        c_minus(6,:,i,j,k) = 0.d0
                        b_bar(6,:,i,j,k) = 0.d0
                     end if
                     a_plus(7,:,i,j,k) = 0.d0
                     a_minus(7,:,i,j,k) = 0.d0
                     b_plus(7,:,i,j,k) = 0.d0
                     b_minus(7,:,i,j,k) = 0.d0
                     c_plus(7,:,i,j,k) = 0.d0
                     c_minus(7,:,i,j,k) = 0.d0
                     b_bar(7,:,i,j,k) = 0.d0
                  end if
               end do
            end do
         end do
c
         do k = 1, kl
            do j = 1, jl, jl-1
               do i = 1, il
                  if(j.eq.1.and.bc_eta_lower(i,k).eq.bcwall.or.
     $                 j.eq.jl.and.bc_eta_upper(i,k).eq.bcwall) then
                     if(wallfunc) then
                        a_plus(6,:,i,j,k) = 0.d0
                        a_minus(6,:,i,j,k) = 0.d0
                        b_plus(6,:,i,j,k) = 0.d0
                        b_minus(6,:,i,j,k) = 0.d0
                        c_plus(6,:,i,j,k) = 0.d0
                        c_minus(6,:,i,j,k) = 0.d0
                        b_bar(6,:,i,j,k) = 0.d0
                     end if
                     a_plus(7,:,i,j,k) = 0.d0
                     a_minus(7,:,i,j,k) = 0.d0
                     b_plus(7,:,i,j,k) = 0.d0
                     b_minus(7,:,i,j,k) = 0.d0
                     c_plus(7,:,i,j,k) = 0.d0
                     c_minus(7,:,i,j,k) = 0.d0
                     b_bar(7,:,i,j,k) = 0.d0
                  end if
               end do
            end do
         end do
c
         if(idimen.gt.2) then
c
            do k = 1, kl, kl-1
               do j = 1, jl
                  do i = 1, il
                     if(k.eq.1.and.bc_zeta_lower(i,j).eq.bcwall.or.
     $                    k.eq.kl.and.bc_zeta_upper(i,j).eq.bcwall) then
                        if(wallfunc) then
                           a_plus(6,:,i,j,k) = 0.d0
                           a_minus(6,:,i,j,k) = 0.d0
                           b_plus(6,:,i,j,k) = 0.d0
                           b_minus(6,:,i,j,k) = 0.d0
                           c_plus(6,:,i,j,k) = 0.d0
                           c_minus(6,:,i,j,k) = 0.d0
                           b_bar(6,:,i,j,k) = 0.d0
                        end if
                        a_plus(7,:,i,j,k) = 0.d0
                        a_minus(7,:,i,j,k) = 0.d0
                        b_plus(7,:,i,j,k) = 0.d0
                        b_minus(7,:,i,j,k) = 0.d0
                        c_plus(7,:,i,j,k) = 0.d0
                        c_minus(7,:,i,j,k) = 0.d0
                        b_bar(7,:,i,j,k) = 0.d0
                     end if
                  end do
               end do
            end do

         end if

      end if
c
c ... add unit matrix
c     
      do ii = 1,nl
         b_bar(ii,ii,:,:,:) = b_bar(ii,ii,:,:,:)+1.d0
c
         if (moving.ge.1.and.gcl.eq.1 ) then
            b_bar(ii,ii,:,:,:) = b_bar(ii,ii,:,:,:)-dt*sr
         end if
c           
         if(dual_t.eq.1) then
            if(dstep.eq.1) then
               b_bar(ii,ii,:,:,:) = b_bar(ii,ii,:,:,:)
     $              + dt/tintvl
            else
               b_bar(ii,ii,:,:,:) = b_bar(ii,ii,:,:,:)
     $              + 1.5d0*dt/tintvl
            end if
         end if
      end do
c
      if ( dabs(ronum).gt.1e-9 ) then
        do 135 k  = 1, kl
        do 135 j  = 1, jl
        do 135 i  = 1, il
          cx1 = dt(i,j,k)*ronum
          b_bar(3,4,i,j,k) = b_bar(3,4,i,j,k)+cx1
          b_bar(4,3,i,j,k) = b_bar(4,3,i,j,k)-cx1
c
          cx2 = cx1/q(i,j,k,1)
          b_bar(3,1,i,j,k) = b_bar(3,1,i,j,k)+q(i,j,k,4)*cx2
          b_bar(4,1,i,j,k) = b_bar(4,1,i,j,k)-q(i,j,k,3)*cx2
135     continue
      end if
C--------------------------------------------------
           do k = 1, kl
               do i = 1, il
		  DO J=1,JL
	                QB(:,:)=b_bar(:,:,I,J,K)
			CALL IVSNR(QB,NL,.1D-06)
	                IN_A(:,:,I,J,K)=QB(:,:)
	            END DO
	         END DO
	      END DO

c..   start Gauss-Seidel Relaxation Sweeping
      
c..   eta-constant line
 
      if(jl.gt.1) then

         do m=1,iter_gs
            do k = 1, kl
               do i = 1, il
                  
                  if(mod(m,2).eq.1) then
                     is =i
                  else
                     is =il+1-i
                  end if
                  
                  capb(:,:,:jl) = b_minus(:,:,is,:,k) 
                  capc(:,:,:jl) = b_plus(:,:,is,:,k) 
                  capa(:,:,:jl) = b_bar(:,:,is,:,k)
c     
c     form right-hand-side vector f(n,j)
c     
c..   initilize zero
                  f = 0.d0

                  cz1=1.0
                  cz2=1.0
                  cx1=1.0
                  cx2=1.0

                  if(k.eq.1)  cz1=0.0
                  if(k.eq.kl) cz2=0.0
                  if(i.eq.1)  cx1=0.0
                  if(i.eq.il) cx2=0.0

                  do j = 1,jl
                     do n = 1,nl
                        do nj = 1,nl
                           f(n,j) = f(n,j) 
     >                          -cx1*
     $                          a_minus(n,nj,is,j,k)*q_n(is-1,j,k,nj)
     >                          -cx2*
     $                          a_plus(n,nj,is,j,k)*q_n(is+1,j,k,nj)
     >                          -cz1*
     $                          c_minus(n,nj,is,j,k)*q_n(is,j,k-1,nj)
     >                          -cz2*
     $                          c_plus(n,nj,is,j,k)*q_n(is,j,k+1,nj)
                        end do
                     end do
                  end do
                  
                  do n = 1, nl
                     f(n,:jl) = f(n,:jl) + rhs(is,:,k,n) 
                  end do

c     solve linear system
c     
                  nblocks = jl

                  if(bc_eta_lower(i,k).ne.bcperi
     $               .or.bc_eta_lower(i,k).ne.20) then
C------------------------->				 
		     DO J=1,JL
c	                FIRST(:)=0.
c	                IF (J.EQ.1) THEN
c	                   FIRST(:)=0.
c		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $					 CAPB(N,NJ,J)*Q_N(IS,J-1,K,NJ)+
     $					 CAPC(N,NJ,J)*Q_N(IS,J+1,K,NJ)

	                      END DO
	                   END DO
c	                END IF
	                FIRST(:)=F(:,J)-FIRST(:)
			DO N=1,NL
	                   DQ(N)=0.
	                   DO NJ=1,NL
	                      DQ(N)=DQ(N)+IN_A(N,NJ,IS,J,K)*FIRST(NJ)
	                   END DO
	                   Q_N(IS,J,K,N)=DQ(N)
	                END DO
	             ENDDO  !END OF J

C-----<-------------------------
		     DO J=JL,1,-1
c	                FIRST(:)=0.
c	                IF (J.EQ.JL) THEN
c	                   FIRST(:)=0.
c		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $					 CAPB(N,NJ,J)*Q_N(IS,J-1,K,NJ)+
     $					 CAPC(N,NJ,J)*Q_N(IS,J+1,K,NJ)
	                      END DO
	                   END DO
c	                END IF
	                FIRST(:)=F(:,J)-FIRST(:)
			DO N=1,NL
	                   DQ(N)=0.
	                   DO NJ=1,NL
	                      DQ(N)=DQ(N)+IN_A(N,NJ,IS,J,K)*FIRST(NJ)
	                   END DO
	                   Q_N(IS,J,K,N)=DQ(N)
	                END DO
	             ENDDO  !END OF J
C	                
C                     call block(nblocks, nl, dim1, capa, capb, capc,
C     $                    soln, f)
                  else
                     soln(:,:jl) = f
                     call pbtrip(capb(:,:,:jl), capa(:,:,:jl),
     $                    capc(:,:,:jl), soln(:,:jl), 1, jl, nl)

                  do n = 1,nl
                     q_n(is,1:jl,k,n) = soln(n,1:jl)
                  end do

                  end if

c
               end do
            end do
         end do

      end if

c..   xi-constant line

      if(il.gt.1) then

         do m=1,iter_gs
            do k = 1, kl
               do j = 1, jl
                  
                  if(mod(m,2).eq.1) js =j
                  if(mod(m,2).eq.0) js =jl+1-j

                  do i = 1,il
                     capb(:,:,i) = a_minus(:,:,i,js,k) 
                     capc(:,:,i) = a_plus(:,:,i,js,k) 
                     capa(:,:,i) = b_bar(:,:,i,js,k)
                  end do

c     
c     form right-hand-side vector f(n,j)
c     
                  do n = 1,nl
                     f(n,:il) = rhs(:il,js,k,n) 
                  end do

                  cy1=1.0
                  cy2=1.0
                  cz1=1.0
                  cz2=1.0
                  
                  if(j.eq.1)  cy1=0.0
                  if(j.eq.jl) cy2=0.0
                  if(k.eq.1)  cz1=0.0
                  if(k.eq.kl) cz2=0.0

                  do i = 1,il
                     do n = 1,nl
                        do nj = 1,nl
                           f(n,i) = f(n,i) 
     >                          -cy1*
     $                          b_minus(n,nj,i,js,k)*q_n(i,js-1,k,nj)
     >                          -cy2*
     $                          b_plus(n,nj,i,js,k)*q_n(i,js+1,k,nj)
     >                          -cz1*
     $                          c_minus(n,nj,i,js,k) *q_n(i,js,k-1,nj)
     >                          -cz2*
     $                          c_plus(n,nj,i,js,k)*q_n(i,js,k+1,nj)
                        end do
                     end do
                  end do

c     
c     solve linear system
c     
                  nblocks = il

                  if(bc_xi_lower(j,k).ne.bcperi
     $               .or.bc_xi_lower(j,k).ne.20) then
C------------------------->				 
		     DO I=1,IL
c	                FIRST(:)=0.
c	                IF (I.EQ.1) THEN
c	                   FIRST(:)=0.
c		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $					 CAPB(N,NJ,I)*Q_N(I-1,JS,K,NJ)+
     $					 CAPC(N,NJ,I)*Q_N(I+1,JS,K,NJ)
	                      END DO
	                   END DO
c	                END IF
	                FIRST(:)=F(:,I)-FIRST(:)
			DO N=1,NL
	                   DQ(N)=0.
	                   DO NJ=1,NL
	                      DQ(N)=DQ(N)+IN_A(N,NJ,I,JS,K)*FIRST(NJ)
	                   END DO
	                   Q_N(I,JS,K,N)=DQ(N)
	                END DO
	             ENDDO  !END OF I

C-----<-------------------------
		     DO I=IL,1,-1
c	                FIRST(:)=0.
c	                IF (I.EQ.IL) THEN
c	                   FIRST(:)=0.
c		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $					 CAPB(N,NJ,I)*Q_N(I-1,JS,K,NJ)+
     $					 CAPC(N,NJ,I)*Q_N(I+1,JS,K,NJ)
	                      END DO
	                   END DO
c	                END IF
	                FIRST(:)=F(:,I)-FIRST(:)
			DO N=1,NL
	                   DQ(N)=0.
	                   DO NJ=1,NL
	                      DQ(N)=DQ(N)+IN_A(N,NJ,I,JS,K)*FIRST(NJ)
	                   END DO
	                   Q_N(I,JS,K,N)=DQ(N)
	                END DO
	             ENDDO  !END OF I
C                     call block(nblocks, nl, dim1, capa, capb, capc,
C     $                    soln, f)
                  else
                     soln(:,:il) = f
                     call  pbtrip(capb(:,:,:il),capa(:,:,:il),
     $                    capc(:,:,:il), soln(:,:il),1, il, nl)

                  do   i = 1,il
                     do   n = 1,nl
                        q_n(i,js,k,n) = soln(n,i)
                     end do
                  end do

                  end if


c     
               end do
            end do
         end do

      end if    

c..   zeta-constant line

      if(kl.gt.1) then

         do m=1,iter_gs
            do j = 1, jl
               do i = 1, il
                  
                  if(mod(m+j,2).eq.1) is =i
                  if(mod(m+j,2).eq.0) is =il+1-i

                  do k = 1, kl
                     capb(:,:,k) = c_minus(:,:,is,j,k) 
                     capc(:,:,k) = c_plus(:,:,is,j,k) 
                     capa(:,:,k) = b_bar(:,:,is,j,k)
                  end do

c     
c     form right-hand-side vector f(n,j)
c     
c..   initilize zero

                  f = 0.d0
                  
                  cx1=1.0
                  cx2=1.0
                  cy1=1.0
                  cy2=1.0
                  
                  if(i.eq.1)  cx1=0.0
                  if(i.eq.il) cx2=0.0
                  if(j.eq.1)  cy1=0.0
                  if(j.eq.jl) cy2=0.0

                  do  k = 1, kl
                     do   n = 1,nl
                        do   nj = 1,nl
                           f(n,k) = f(n,k) 
     >                          -cx1*
     $                          a_minus(n,nj,is,j,k)*q_n(is-1,j,k,nj)
     >                          -cx2*
     $                          a_plus(n,nj,is,j,k)*q_n(is+1,j,k,nj)
     >                          -cy1*
     $                          b_minus(n,nj,is,j,k) *q_n(is,j-1,k,nj)
     >                          -cy2*
     $                          b_plus(n,nj,is,j,k)*q_n(is,j+1,k,nj)
                        end do
                     end do
                  end do

                  do   k = 1, kl
                     do   n = 1, nl
                        f(n,k) = f(n,k) + rhs(is,j,k,n) 
                     end do
                  end do
c     
c     solve linear system
c     
                  nblocks = kl

                  if(bc_zeta_lower(i,j).ne.bcperi
     $               .or.bc_zeta_lower(i,j).ne.20) then
C------------------------->				 
		     DO K=1,KL
c	                FIRST(:)=0.
c	                IF (K.EQ.1) THEN
c	                   FIRST(:)=0.
c		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $					 CAPB(N,NJ,K)*Q_N(IS,J,K-1,NJ)+
     $					 CAPC(N,NJ,K)*Q_N(IS,J,K+1,NJ)
	                      END DO
	                   END DO
c	                END IF
	                FIRST(:)=F(:,K)-FIRST(:)
			DO N=1,NL
	                   DQ(N)=0.
	                   DO NJ=1,NL
	                      DQ(N)=DQ(N)+IN_A(N,NJ,IS,J,K)*FIRST(NJ)
	                   END DO
	                   Q_N(IS,J,K,N)=DQ(N)
	                END DO
	             ENDDO  !END OF K

C-----<-------------------------
		     DO K=KL,1,-1
c	                FIRST(:)=0.
c	                IF (K.EQ.KL) THEN
c	                   FIRST(:)=0.
c		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $					 CAPB(N,NJ,K)*Q_N(IS,J,K-1,NJ)+
     $					 CAPC(N,NJ,K)*Q_N(IS,J,K+1,NJ)
	                      END DO
	                   END DO
c	                END IF
                        FIRST(:)=F(:,K)-FIRST(:)
			DO N=1,NL
	                   DQ(N)=0.
	                   DO NJ=1,NL
	                      DQ(N)=DQ(N)+IN_A(N,NJ,IS,J,K)*FIRST(NJ)
	                   END DO
	                   Q_N(IS,J,K,N)=DQ(N)
	                END DO
	             ENDDO  !END OF K
C                     call block(nblocks, nl, dim1, capa, capb, capc,
C     $                    soln, f)
                  else
                     soln(:,:kl) = f
                     call pbtrip(capb(:,:,:kl), capa(:,:,:kl),
     $                    capc(:,:,:kl), soln(:,:kl), 1, kl, nl)

                  do   k = 1,kl
                     do   n = 1,nl
                        q_n(is,j,k,n) = soln(n,k)
                     end do
                  end do

                  end if
c     

c     
               end do
            end do
         end do

      end if

      return
      end
