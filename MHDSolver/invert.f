      subroutine invert(index, dt, bc_xi_lower, bc_xi_upper,
     $     bc_eta_lower, bc_eta_upper, bc_zeta_lower, bc_zeta_upper,
     $     delta, q, flux,
     $     gamma, tref, prandtl, rhs, x, y, z, vol, il, jl, kl, nl,
     $     dim2, dim1, ilower, iupper, jlower, jupper, klower, kupper,
     $     rank, visturb, wallfunc, ke,qt,
     $     udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $     wdk, tko, qiii, diver,control)
c
c     defines and solves the linear system for xi, eta and zeta factorization
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
     $     il, jl, kl,          ! cell number in xi, eta, zeta
     $     nl,                  ! equation number
     $     dim2, dim1,          ! maximum surface, cell number
     $     bc_xi_lower(jl, kl),
     $     bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl),
     $     bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl),
     $     bc_zeta_upper(il, jl),
     $     ilower, iupper,      ! array memory limit in xi
     $     jlower, jupper,      ! array memory limit in eta
     $     klower, kupper,      ! array memory limit in zeta
     $     rank,                ! to deleted
     $     index                ! 1 xi, 2 eta, 3 zeta

      integer::
c     $     limiter,
c     $     rhs_scheme,
c     $     lhs_scheme,
     $     idimen,moving,
c     $     lhs_order,
c     $     rhs_order,
     $     dual_t              ! not applied in AF yet

      double precision,intent(in)::
     $     diver(ilower:iupper,jlower:jupper,klower:kupper)


      double precision, intent(in)::
     $     delta(5)
      
      logical, intent(in)::
     $     wallfunc,
     $     ke
      logical::
     $     inviscid,            ! invisid flow
     $     suther              ! sutherland method for viscosity

      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z              ! mesh info

      double precision, intent(in)::
     $     q(ilower:iupper,jlower:jupper,klower:kupper,nl),
     $     qiii(ilower:iupper,jlower:jupper,klower:kupper,nl),
     $     vol(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)
     $                          ! volume info

      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, 3)::
     $     qt

      double precision, intent(in)::
     $     flux(il+1, jl+1, kl+1, nl), ! viscid surface flux
     $     dt(il,jl,kl)        ! time marching interval
      double precision::
     $     gamma,
     $     prandtl,
     $     prt,
     $     tref,
     $     machinf,
     $     epsfactor,
     $     kfactor,
     $     reynolds,            ! reynolds number
     $     ptotal,
     $     ttotal,
     $     angl1, angl2,        ! inlet flow angles
     $     poutlet,              ! outlet static pressure
     $     ronum

      double precision, intent(in)::
     $  visturb(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)

      double precision, intent(inout)::
     $     rhs(il, jl, kl, nl)  ! right hand side

      double precision, dimension(jlower:jupper,
     $     klower:kupper,4), intent(in):: udi, vdi, wdi

      double precision, dimension(ilower:iupper,
     $     klower:kupper,4), intent(in):: udj, vdj, wdj

      double precision, dimension(ilower:iupper,
     $     jlower:jupper,4), intent(in):: udk, vdk, wdk

c
c     Local variables
c
      integer::
     $     n,                   ! equation index
     $     imax,                ! maximum cell number,
     $     i1max, i2max,        ! maximum cell number on another 2 dirs
     $     i1, i2,              ! interaction indices for i1max, i2max
     $     bclower, bcupper
c
c     SA 1eq turbulent model contants
c
      double precision, intent(in):: tko

      double precision, dimension(nl,nl,dim1)::
     $     capa, capb, capc,
     $     f(nl,dim1), soln(nl,dim1)

      double precision::
     $     qtline(ilower:dim2,3)

      double precision, dimension(nl,nl,3,2)::
     $     ci_nr
c
c     Subroutine Start
c
c--------------------------parameter transfer
      dual_t=control%dual_t
      idimen=control%idimen
      moving=control%moving

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
      ttotal=control%ttotal
c--------------------------------------------------------
      select case (index)
      case (1)
         imax = il
         i1max = jl
         i2max = kl
      case (2)
         imax = jl
         i1max = il
         i2max = kl
      case (3)
         imax = kl
         i1max = il
         i2max = jl
      end select

      do i1 = 1, i1max
         do i2 = 1, i2max

            select case (index)
            case (1)
               bclower = bc_xi_lower(i1,i2)
               bcupper = bc_xi_upper(i1,i2)
            case (2)
               bclower = bc_eta_lower(i1,i2)
               bcupper = bc_eta_upper(i1,i2)
            case (3)
               bclower = bc_zeta_lower(i1,i2)
               bcupper = bc_zeta_upper(i1,i2)
            end select
            
            if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
              if (index.eq.1) then
                qtline(ilower:iupper,:) = qt(ilower:iupper,i1,i2,:)
              elseif (index.eq.2) then
                qtline(jlower:jupper,:) = qt(i1,jlower:jupper,i2,:)
              else
                qtline(klower:kupper,:) = qt(i1,i2,klower:kupper,:)
              end if
            end if

            call lhs_matrix(index, i1, i2, dt, x, y, z, vol,
     $           bclower, bcupper, delta, wallfunc,
     $           capa, capb,
     $           capc, q, flux, il, jl, kl, nl, dim2, dim1, ilower,
     $           iupper, jlower, jupper, klower, kupper, visturb, ke,
     $           udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $           wdk, qtline, ci_nr, tko, qiii, diver, control)           
            
            do n = 1, nl
               capa(n,n,:) = 1.d0 + capa(n,n,:)
            end do

            do n = 1, nl
               select case (index)
               case (1)
                  f(n,:imax) = rhs(:,i1,i2,n)
               case (2)
                  f(n,:imax) = rhs(i1,:,i2,n)
               case (3)
                  f(n,:imax) = rhs(i1,i2,:,n)
               end select
            end do
            
            call block(imax, nl, dim1, capa, capb, capc, soln, f)

            do n = 1, nl
               select case (index)
               case (1)
                  rhs(:,i1,i2,n) = soln(n,:imax)
               case (2)
                  rhs(i1,:,i2,n) = soln(n,:imax)
               case (3)
                  rhs(i1,i2,:,n) = soln(n,:imax)
               end select
            end do

         end do
      end do

      return
      end
