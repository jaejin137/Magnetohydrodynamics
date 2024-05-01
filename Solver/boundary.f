      subroutine boundary(index, dim2, i1, i2, il, jl, kl, nl,
     $     ilower, iupper, jlower, jupper, klower, bclower, bcupper,
     $     q, kupper, rho, rhoe, rhou, rhov, rhow, rhotk, rhotw, ke,
     $     vbds, x, y, z, blen, qconer,
     $     qiii, diver, control)
c
c     set up boundary conditions for inviscid part (left side)
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
     $     il, jl, kl,
     $     nl, blen,
     $     dim2,
     $     index,
     $     i1, i2,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     bclower, bcupper

      integer::
     $     moving, precondition,
     $     main_dir
c     
      double precision::
     $     machinf, gamma, angl1, angl2, poutlet,
     $     tintvl, ronum
c     
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, dimension(3,4), intent(in)::
     $     vbds                  ! mesh end boun vels of the line
c
c      double precision, intent(inout)::  !  ????
c     $     ptotal, ttotal
      double precision::
     $     ptotal, ttotal
c     
      double precision, intent(inout)::
     $     q(ilower:iupper,jlower:jupper,klower:kupper,nl),
     $     qiii(ilower:iupper,jlower:jupper,klower:kupper,nl)
c     
      double precision, dimension(ilower:dim2),intent(inout)::
     $     rho,
     $     rhoe,
     $     rhou,
     $     rhov,
     $     rhow,diver
c
      double precision, intent(in):: qconer(nl,2)
c     
c     LOCAL VARIABLES
c     
      integer::
     $     bcindex,
     $     i, j, k,
     $     in(blen), out(blen),
     $     num,
     $     nlm, nlp,
     $     bndp,              ! periodical boundary point
     $     pdir,              ! periodic boundary direction 1: upper, -1: lower
     $     ii, iin, iout
c     
      double precision::
     $     gamma1,
     $     tout,
     $     k1, k2,
     $     pout, pin,
     $     uout, vout, wout, qqout,
     $     dxp, dyp, coef, x2db, y2db, z2db
c     
      double precision, dimension(nl)::
     $     qin, qout
c     
c     variables for bc9, p extrapolation
c
      integer::
     $     extra
c
      double precision::
     $     pratio, mach, kcoef, main_vel
c     
c ... variables for komega model
c
      double precision, dimension(ilower:dim2),intent(inout)::
     $     rhotk,
     $     rhotw
c
      double precision:: x0(3,2),x1(3),x2(3),b0
c
      logical, intent(in)::
     $     ke
c
      integer:: dl1, dl2
      integer:: ip, jp, kp, iii, itplk, kk
c
      double precision:: dp, xp, yp, zp, dx1, dx2, dy1, dy2, dz1, dz2
c     
      double precision:: temp, rinv
c
      double precision:: temp_bnd  !for fixed temperature, BCTYPE=19
c
      double precision::
     $     dag, rg, vg, wg, utmp, vtmp, ags, rginv,
     $     psi, theta, rhob, rhoub, rhovb, rhowb, rhoeb, as2,
     $     rd, ud, vd, wd, cd, ub, vb, wb, qq,
     $     ptb, ttb, angr, rad, ttl, ptl, p_bar,
     $     crm, co2, gammap, cb, tb, pb
c     
c     *** SUBROUTINE START ***
c     

      main_dir=control%main_dir
      moving=control%moving
      precondition=control%precondition

      angl1=control%angl1
      angl2=control%angl2
      gamma=control%gamma
      machinf=control%machinf
      poutlet=control%poutlet
      ptotal=control%ptotal
      ronum=control%ronum
      tintvl=control%tintvl
      ttotal=control%ttotal
c--------------------------------------------------------
      gamma1 = gamma-1.d0
      dag = 0.174532925d0
c
      select case (index)
      case(1)
        j = i1
        k = i2
        nlm = il
        nlp = il+1
        dl1 = jl
        dl2 = kl
      case(2)
        i = i1
        k = i2
        nlm = jl
        nlp = jl+1
        dl1 = il
        dl2 = kl
      case(3)
        i = i1
        j = i2 
        nlm = kl
        nlp = kl+1
        dl1 = il
        dl2 = jl
      end select
c
      do num = 1, 2
c
        select case (num)
        case(1)
          bcindex = bclower
          do ii = 1,blen
            in(ii) = ii
            out(ii) = 1-ii
          end do
          bndp = nlm
          pdir = -1
          temp_bnd = 4.0d0 !for fixed temperature
        case(2)
          bcindex = bcupper
          do ii = 1,blen
            in(ii) = nlm-ii+1
            out(ii) = nlp+ii-1
          end do
          bndp = 1
          pdir = 1
          temp_bnd = 1.0d0  !for fixed temperature
        end select
c
c--------------------------------
c
        if (precondition.ge.2) then
c
        select case (index)
        case(1)
          qin(:) = qiii(in(1),j,k,:)
        case(2)
          qin(:) = qiii(i,in(1),k,:)
        case(3)
          qin(:) = qiii(i,j,in(1),:)
        end select
        qout(:) = qconer(:,num)
c
        select case (bcindex)
        case(1)                ! zero gradient
          iin = in(1)
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = rho(iin)
            rhoe(iout) = rhoe(iin)
            rhou(iout) = rhou(iin)
            rhov(iout) = rhov(iin)
            rhow(iout) = rhow(iin)
            if(nl.eq.6) then
              rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
            end if
          end do
        case(2,15,17)                ! supersonic inflow
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = qout(1)
            rhou(iout) = qout(2)
            rhov(iout) = qout(3)
            rhow(iout) = qout(4)
            rhoe(iout) = qout(5)
            if(nl.eq.6) then
               rhotk(iout) = qout(6)
            else if(nl.eq.7) then
               rhotk(iout) = qout(6)
               rhotw(iout) = qout(7)
            end if
          end do
        case(3)                ! no slip adiabatic boundary
          if(moving.eq.0) then ! fixed bound
            do ii = 1,blen
              iin = in(ii)
              iout = out(ii)
              rho(iout) = rho(iin)
              rhoe(iout) = rhoe(iin)
              rhou(iout) = -rhou(iin)
              rhov(iout) = -rhov(iin)
              rhow(iout) = -rhow(iin)
              if(nl.eq.6) then
                rhotk(iout) = -rhotk(iin)
              else if(nl.eq.7) then
                rhotk(iout) = -rhotk(iin)
                rhotw(iout) = -rhotw(iin)
              end if
            end do
          else                ! moving bound, not available for kw model
             write(*,*)'from boundary,this case is not available'
            if ( num.eq.1 ) then
              dxp = ( x(i+1,1,k)-x(i,1,k) )*
     >              ( y(i,2,k)-y(i,1,k)+y(i+1,2,k)-y(i+1,1,k) )
              dyp = ( y(i+1,1,k)-y(i,1,k) )*
     >              ( x(i,2,k)-x(i,1,k)+x(i+1,2,k)-x(i+1,1,k) )
              coef = 2.d0/( dxp-dyp )
              dxp = -coef*( y(i+1,1,k)-y(i,1,k) )
              dyp =  coef*( x(i+1,1,k)-x(i,1,k) )
            else
              dxp = (x(i+1,jl+1,k)-x(i,jl+1,k) )*
     >              (y(i,jl+1,k)-y(i,jl,k)+y(i+1,jl+1,k)-y(i+1,jl,k))
              dyp = (y(i+1,jl+1,k)-y(i,jl+1,k) )*
     >              (x(i,jl+1,k)-x(i,jl,k)+x(i+1,jl+1,k)-x(i+1,jl,k))
              coef = 2.d0/( dxp-dyp )
              dxp = -coef*( y(i+1,jl+1,k)-y(i,jl+1,k) )
              dyp =  coef*( x(i+1,jl+1,k)-x(i,jl+1,k) )
            end if
c
            x2db = vbds(1,1) - vbds(1,3)
            y2db = vbds(2,1) - vbds(2,3)
c
            do ii = 1,blen
              iin = in(ii)
              iout = out(ii)
               rho(iout) = rho(iin)
               rhou(iout) = 2.d0*rho(iout)*vbds(1,1)-rhou(iin)
               rhov(iout) = 2.d0*rho(iout)*vbds(2,1)-rhov(iin)
               rhow(iout) =                    rhow(iin)
               uout = rhou(iout)/rho(iout)
               vout = rhov(iout)/rho(iout)
               wout = rhow(iout)/rho(iout)
               qqout = 0.5d0 * (uout**2 + vout**2 + wout**2)
               coef = rho(iout)/( dxp*dxp + dyp*dyp )/tintvl
               pin = gamma1*(qin(5) -.5d0*(rhou(iin)**2+rhov(iin)**2
     $              +rhow(iin)**2)/rho(iin)) ! zero pressure gradient
               pout = pin + coef*(dxp*x2db + dyp*y2db)
               rhoe(iout) = pout/gamma1 + rho(iout)*qqout
            end do
          end if
c
        case(19)                ! no slip isothermal boundary
          if(moving.eq.0) then ! 
            do ii = 1,blen
              iin = in(ii)
              iout = out(ii)
              rho(iout) = rho(iin)
              rhoe(iout) = temp_bnd
              rhou(iout) = -rhou(iin)
              rhov(iout) = -rhov(iin)
              rhow(iout) = -rhow(iin)
              if(nl.eq.6) then
                rhotk(iout) = -rhotk(iin)
              else if(nl.eq.7) then
                rhotk(iout) = -rhotk(iin)
                rhotw(iout) = -rhotw(iin)
              end if
            end do
          else                ! moving bound, not available for kw model
             write(*,*)'from boundary,this case is not available'
          end if
c
        case(4)                ! inlet BC with w = 0  (temporary BC)
          select case (main_dir)
          case(1)
            if(num.eq.1) then
              x0(1,1)=x(1,j+1,k)-x(1,j,k+1)
              x0(2,1)=y(1,j+1,k)-y(1,j,k+1)
              x0(3,1)=z(1,j+1,k)-z(1,j,k+1)
              x0(1,2)=x(1,j+1,k+1)-x(1,j,k)
              x0(2,2)=y(1,j+1,k+1)-y(1,j,k)
              x0(3,2)=z(1,j+1,k+1)-z(1,j,k)
            else
              x0(1,2)=x(il+1,j+1,k)-x(il+1,j,k+1)
              x0(2,2)=y(il+1,j+1,k)-y(il+1,j,k+1)
              x0(3,2)=z(il+1,j+1,k)-z(il+1,j,k+1)
              x0(1,1)=x(il+1,j+1,k+1)-x(il+1,j,k)
              x0(2,1)=y(il+1,j+1,k+1)-y(il+1,j,k)
              x0(3,1)=z(il+1,j+1,k+1)-z(il+1,j,k)
            end if
          case(2)
            if(num.eq.1) then
              x0(1,1)=x(i,1,k+1)-x(i+1,1,k)
              x0(2,1)=y(i,1,k+1)-y(i+1,1,k)
              x0(3,1)=z(i,1,k+1)-z(i+1,1,k)
              x0(1,2)=x(i+1,1,k+1)-x(i,1,k)
              x0(2,2)=y(i+1,1,k+1)-y(i,1,k)
              x0(3,2)=z(i+1,1,k+1)-z(i,1,k)
            else
              x0(1,2)=x(i,jl+1,k+1)-x(i+1,jl+1,k)
              x0(2,2)=y(i,jl+1,k+1)-y(i+1,jl+1,k)
              x0(3,2)=z(i,jl+1,k+1)-z(i+1,jl+1,k)
              x0(1,1)=x(i+1,jl+1,k+1)-x(i,jl+1,k)
              x0(2,1)=y(i+1,jl+1,k+1)-y(i,jl+1,k)
              x0(3,1)=z(i+1,jl+1,k+1)-z(i,jl+1,k)
            end if
          case(3)
            if(num.eq.1) then
              x0(1,1)=x(i+1,j,1)-x(i,j+1,1)
              x0(2,1)=y(i+1,j,1)-y(i,j+1,1)
              x0(3,1)=z(i+1,j,1)-z(i,j+1,1)
              x0(1,2)=x(i+1,j+1,1)-x(i,j,1)
              x0(2,2)=y(i+1,j+1,1)-y(i,j,1)
              x0(3,2)=z(i+1,j+1,1)-z(i,j,1)
            else
              x0(1,2)=x(i+1,j,kl+1)-x(i,j+1,kl+1)
              x0(2,2)=y(i+1,j,kl+1)-y(i,j+1,kl+1)
              x0(3,2)=z(i+1,j,kl+1)-z(i,j+1,kl+1)
              x0(1,1)=x(i+1,j+1,kl+1)-x(i,j,kl+1)
              x0(2,1)=y(i+1,j+1,kl+1)-y(i,j,kl+1)
              x0(3,1)=z(i+1,j+1,kl+1)-z(i,j,kl+1)
            end if
          end select

          x1(1)=x0(2,1)*x0(3,2)-x0(2,2)*x0(3,1)
          x1(2)=x0(1,2)*x0(3,1)-x0(1,1)*x0(3,2)
          x1(3)=x0(1,1)*x0(2,2)-x0(1,2)*x0(2,1)
          b0 = dsqrt(x1(1)**2+x1(2)**2+x1(3)**2)
          x2(1) = x1(1)/b0
          x2(2) = x1(2)/b0
          x2(3) = x1(3)/b0
c
          iin = in(1)
          b0 = dsqrt(rhou(iin)**2+rhov(iin)**2+rhow(iin)**2)
          uout = b0*x2(1)
          vout = b0*x2(2)
          wout = 0.d0
c
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = rho(iin) 
            rhou(iout) = uout
            rhov(iout) = vout
            rhow(iout) = wout
            rhoe(iout) = rhoe(iin)
          end do
c
        case(5)                ! subsonic outflow (fixed static pressure at freestream)
          iin = in(1)
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = poutlet
            rhou(iout) = rhou(iin)
            rhov(iout) = rhov(iin)
            rhow(iout) = rhow(iin)
            rhoe(iout) = rhoe(iin)
            if(nl.eq.6) then
              rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
              if(ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
            end if
          end do
        case(6)                ! subsonic inflow (zero gradient in p)
          do ii = 1,blen
            iout = out(ii)
            iin = iout-pdir
            rho(iout) = qin(1)
            rhou(iout) = qout(2)
            rhov(iout) = qout(3)
            rhow(iout) = qout(4)
            rhoe(iout) = qout(5)
            if(nl.eq.6) then
              rhotk(iout) = qout(6)
            else if(nl.eq.7) then
              rhotk(iout) = qout(6)
              rhotw(iout) = qout(7)
              if(ke) then
                pout = gamma1*(qin(5)-rhotk(iin) -.5d0*
     $            (rhou(iin)**2+rhov(iin)**2+rhow(iin)**2)/rho(iin))
                rhoe(iout) = pout/gamma1 + rhotk(iout)+.5d0*
     $            (rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)/rho(iout)
              end if
            end if
          end do
        case(7)                ! inner boundary for mpi communication
        case(100)              ! specified boundary

        case(8)       ! symmetry boundary
          do ii = 1,blen
            iin = in(ii)
            iout = out(ii)
            rho(iout) = rho(iin)
            rhoe(iout) = rhoe(iin)
            if(nl.eq.6) then
               rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
               rhotk(iout) = rhotk(iin)
               rhotw(iout) = rhotw(iin)
            end if
            select case (index)
            case(1)
              jp = j + 1
              kp = k + 1
              select case (num)
              case(1)
                dx1 = x(1,jp,kp) - x(1,j,k)
                dy1 = y(1,jp,kp) - y(1,j,k)
                dz1 = z(1,jp,kp) - z(1,j,k)
                dx2 = x(1,j,kp)  - x(1,jp,k)
                dy2 = y(1,j,kp)  - y(1,jp,k)
                dz2 = z(1,j,kp)  - z(1,jp,k)
              case(2)
                dx1 = x(il+1,jp,kp) - x(il+1,j,k)
                dy1 = y(il+1,jp,kp) - y(il+1,j,k)
                dz1 = z(il+1,jp,kp) - z(il+1,j,k)
                dx2 = x(il+1,j,kp)  - x(il+1,jp,k)
                dy2 = y(il+1,j,kp)  - y(il+1,jp,k)
                dz2 = z(il+1,j,kp)  - z(il+1,jp,k)
              end select
            case(2)
              ip = i + 1
              kp = k + 1
              select case (num)
              case(1)
                dx1 = x(ip,1,kp) - x(i,1,k)
                dy1 = y(ip,1,kp) - y(i,1,k)
                dz1 = z(ip,1,kp) - z(i,1,k)
                dx2 = x(ip,1,k)  - x(i,1,kp)
                dy2 = y(ip,1,k)  - y(i,1,kp)
                dz2 = z(ip,1,k)  - z(i,1,kp)
              case(2)
                dx1 = x(ip,jl+1,kp) - x(i,jl+1,k)
                dy1 = y(ip,jl+1,kp) - y(i,jl+1,k)
                dz1 = z(ip,jl+1,kp) - z(i,jl+1,k)
                dx2 = x(ip,jl+1,k)  - x(i,jl+1,kp)
                dy2 = y(ip,jl+1,k)  - y(i,jl+1,kp)
                dz2 = z(ip,jl+1,k)  - z(i,jl+1,kp)
              end select
            case(3)
              ip = i + 1
              jp = j + 1
              select case (num)
              case(1)
                dx1 = x(ip,jp,1) - x(i,j,1)
                dy1 = y(ip,jp,1) - y(i,j,1)
                dz1 = z(ip,jp,1) - z(i,j,1)
                dx2 = x(i,jp,1)  - x(ip,j,1)
                dy2 = y(i,jp,1)  - y(ip,j,1)
                dz2 = z(i,jp,1)  - z(ip,j,1)
              case(2)
                dx1 = x(ip,jp,kl+1) - x(i,j,kl+1)
                dy1 = y(ip,jp,kl+1) - y(i,j,kl+1)
                dz1 = z(ip,jp,kl+1) - z(i,j,kl+1)
                dx2 = x(i,jp,kl+1)  - x(ip,j,kl+1)
                dy2 = y(i,jp,kl+1)  - y(ip,j,kl+1)
                dz2 = z(i,jp,kl+1)  - z(ip,j,kl+1)
              end select
            end select
                                                                                
            xp = dabs(dy1 * dz2 - dz1 * dy2)
            yp = dabs(dz1 * dx2 - dx1 * dz2)
            zp = dabs(dx1 * dy2 - dy1 * dx2)
            dp = dmax1(xp,yp,zp)

            if (dp.eq.xp) then
c *** symmetry about x-plane (x = cont)
              rhou(iout) = -rhou(iin)
              rhov(iout) = rhov(iin)
              rhow(iout) = rhow(iin)
            else if (dp.eq.yp) then
c *** symmetry about y-plane (y = cont)
              rhou(iout) = rhou(iin)
              rhov(iout) = -rhov(iin)
              rhow(iout) = rhow(iin)
            else if (dp.eq.zp) then
c *** symmetry about z-plane (z = cont)
              rhou(iout) = rhou(iin)
              rhov(iout) = rhov(iin)
              rhow(iout) = -rhow(iin)
            end if
          end do
c
        case(9)         ! fixed total pressure and temperature at subsonic inlet
          extra = 0           ! 0: vel, 1: p extrapolated
          k1 = dtan(angl1)    ! angle between first and main direction
          k2 = dtan(angl2)    ! angle between second and main direction
                              ! main direction is one of u, v, and w
                              ! first and second direction is the first and
                              ! second direction in the order of u, v, and w
                              ! except the selected main direction
          if(extra.eq.0) then
            do ii = 1,blen
              iout = out(ii)
              iin = iout-pdir
              select case (main_dir)
              case(1)
                uout = rhou(iin)
                vout = uout*k1
                wout = uout*k2
              case(2)
                vout = rhov(iin)
                uout = vout*k1
                wout = vout*k2
              case(3)
                wout = rhow(iin)
                uout = wout*k1
                vout = wout*k2
              end select

              rhou(iout)=uout
              rhov(iout)=vout
              rhow(iout)=wout
              qqout = 0.5d0*(uout**2+vout**2+wout**2)
      	      call qstat(ptotal,ttotal,0.0d0,temp)
	      call htvp(ttotal,temp,ptotal,qqout
     $		 ,rhoe(out),rho(out),rinv)

              if(nl.eq.6) then
                rhotk(iout) = rho(iout)*qout(6)/qout(1)
              else if(nl.eq.7) then
                rhotk(iout) = rho(iout)*qout(6)/qout(1)
                rhotw(iout) = rho(iout)*qout(7)/qout(1)
                if(ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
              end if
            end do
          else if(extra.eq.1) then
            do ii = 1,blen
              iout = out(ii)
              iin = iout-pdir
              pout = gamma1*(qin(5) -.5d0*
     $             (rhou(iin)**2+rhov(iin)**2+rhow(iin)**2)/rho(iin))
              pratio = (ptotal/pout)**(gamma1/gamma) ! = ttotal/t
              mach = dsqrt(2.d0/gamma1*(pratio-1.d0)) ! local mach number
              tout = ttotal/pratio
              kcoef = dsqrt(1.d0+k1*k1+k2*k2)
              main_vel = (1.d0/kcoef)*mach*dsqrt(tout)/machinf
              select case (main_dir)
              case(1)
                uout = main_vel
                vout = main_vel*k1
                wout = main_vel*k2
              case(2)
                vout = main_vel
                uout = main_vel*k1
                wout = main_vel*k2
              case(3)
                wout = main_vel
                uout = main_vel*k1
                vout = main_vel*k2
              end select
              qqout = .5d0*(main_vel*kcoef)**2
              rho(out) = gamma*machinf*machinf*pout/tout
              rhou(iout) = rho(iout)*uout
              rhov(iout) = rho(iout)*vout
              rhow(iout) = rho(iout)*wout
              rhoe(iout) = rho(iout)*
     $             (tout/(gamma*gamma1*machinf*machinf)+qqout)
            end do
          end if
        case(10,20)               ! periodical boundary
          do ii = 1,blen
            iout = out(ii)
            iin = in(ii)
            rho(iout)=rho(iin)
            rhou(iout)=rhou(iin)
            rhov(iout)=rhov(iin)
            rhow(iout)=rhow(iin)
            rhoe(iout)=rhoe(iin)
            if(nl.eq.6) then
              rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
            end if
          end do
        case(11)                ! subsonic outflow (fixed static pressure)
          iin = in(1)
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = qout(1)
            rhou(iout) = rhou(iin)
            rhov(iout) = rhov(iin)
            rhow(iout) = rhow(iin)
            rhoe(iout) = rhoe(iin)
            if(nl.eq.6) then
              rhotk(iout) = -rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
              if(ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
            end if
          end do
        case(101,102,103,104,105,106,107,108,109,110)                ! no slip isothermal boundary
          if(moving.eq.0) then ! 
            do ii = 1,blen
              iin = in(ii)
              iout = out(ii)
              rho(iout) = rho(iin)
              rhoe(iout) = control%t(bcindex-100)
              rhou(iout) = -rhou(iin)
              rhov(iout) = -rhov(iin)
              rhow(iout) = -rhow(iin)
              if(nl.eq.6) then
                rhotk(iout) = -rhotk(iin)
              else if(nl.eq.7) then
                rhotk(iout) = -rhotk(iin)
                rhotw(iout) = -rhotw(iin)
              end if
            end do
          else                ! moving bound, not available for kw model
             write(*,*)'from boundary,this case is not available'
          end if
c
        end select
        select case (index)
        case(1)
          do ii = 1,blen
            iout = out(ii)
            qiii(iout,j,k,1) = rho(iout)
            qiii(iout,j,k,2) = rhou(iout)
            qiii(iout,j,k,3) = rhov(iout)
            qiii(iout,j,k,4) = rhow(iout)
            qiii(iout,j,k,5) = rhoe(iout)
            if (nl.eq.6) then
              qiii(iout,j,k,6) = rhotk(iout)
            end if

            call qstat(rho(iout),rhoe(iout),0.0d0,rinv)     
            q(iout,j,k,1) = 1.0d0/rinv
c------same as 1
c            q(iout,j,k,1) = gamma*machinf**2*rho(iout)/rhoe(iout)

            q(iout,j,k,2) = rhou(iout)/rinv
            q(iout,j,k,3) = rhov(iout)/rinv
            q(iout,j,k,4) = rhow(iout)/rinv
            qqout = 0.5d0*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
            call enthalpy(temp,rhoe(iout),rinv,rho(iout))            
            q(iout,j,k,5) = (temp+qqout)/rinv-rho(iout)

c-------same as 1
c            temp=gamma/(gamma-1.0)*rho(iout)/q(iout,j,k,1)
c            q(iout,j,k,5) = q(iout,j,k,1)*(temp+qqout)-rho(iout)
          end do

        case(2)
          do ii = 1,blen
            iout = out(ii)
            qiii(i,iout,k,1) = rho(iout)
            qiii(i,iout,k,2) = rhou(iout)
            qiii(i,iout,k,3) = rhov(iout)
            qiii(i,iout,k,4) = rhow(iout)
            qiii(i,iout,k,5) = rhoe(iout)
            if (nl.eq.6) then
              qiii(i,iout,k,6) = rhotk(iout)
            end if

            call qstat(rho(iout),rhoe(iout),0.0d0,rinv)     
            q(i,iout,k,1) = 1.0d0/rinv
c------same as 1
c            q(i,iout,k,1) = gamma*machinf**2*rho(iout)/rhoe(iout)

            q(i,iout,k,2) = rhou(iout)/rinv
            q(i,iout,k,3) = rhov(iout)/rinv
            q(i,iout,k,4) = rhow(iout)/rinv
            qqout = 0.5d0*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
            call enthalpy(temp,rhoe(iout),rinv,rho(iout))            
            q(i,iout,k,5) = (temp+qqout)/rinv-rho(iout)

c-------same as 1
c            temp=gamma/(gamma-1.0)*rho(iout)/q(i,iout,k,1)
c            q(i,iout,k,5) = q(i,iout,k,1)*(temp+qqout)-rho(iout)
          end do
        case(3)
          do ii = 1,blen
            iout = out(ii)
            qiii(i,j,iout,1) = rho(iout)
            qiii(i,j,iout,2) = rhou(iout)
            qiii(i,j,iout,3) = rhov(iout)
            qiii(i,j,iout,4) = rhow(iout)
            qiii(i,j,iout,5) = rhoe(iout)
            if (nl.eq.6) then
              qiii(i,j,iout,6) = rhotk(iout)
            end if

            call qstat(rho(iout),rhoe(iout),0.0d0,rinv)     
            q(i,j,iout,1) = 1.0d0/rinv
c------same as 1
c            q(i,j,iout,1) = gamma*machinf**2*rho(iout)/rhoe(iout)

            q(i,j,iout,2) = rhou(iout)/rinv
            q(i,j,iout,3) = rhov(iout)/rinv
            q(i,j,iout,4) = rhow(iout)/rinv
            qqout = 0.5d0*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
            call enthalpy(temp,rhoe(iout),rinv,rho(iout))            
            q(i,j,iout,5) = (temp+qqout)/rinv-rho(iout)

c-------same as 1
c            temp=gamma/(gamma-1.0)*rho(iout)/q(i,j,iout,1)
c            q(i,j,iout,5) = q(i,j,iout,1)*(temp+qqout)-rho(iout)
           end do
        end select
c
c--------------------------------
        else
c---------------old--------------
c--- precondition = 0 -----------
c
        select case (index)
        case(1)
          qin(:) = q(in(1),j,k,:)
        case(2)
          qin(:) = q(i,in(1),k,:)
        case(3)
          qin(:) = q(i,j,in(1),:)
        end select
        qout(:) = qconer(:,num)
c
        select case (bcindex)
        case(1)                ! zero gradient
          iin = in(1)
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = rho(iin)
            rhoe(iout) = rhoe(iin)
            rhou(iout) = rhou(iin)
            rhov(iout) = rhov(iin)
            rhow(iout) = rhow(iin)
            if(nl.eq.6) then
              rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
            end if
          end do
        case(2,15,17)                ! supersonic inflow
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = qout(1)
            rhou(iout) = qout(2)
            rhov(iout) = qout(3)
            rhow(iout) = qout(4)
            rhoe(iout) = qout(5)
            if(nl.eq.6) then
               rhotk(iout) = qout(6)
            else if(nl.eq.7) then
               rhotk(iout) = qout(6)
               rhotw(iout) = qout(7)
            end if
          end do
        case(12)               ! special periodical boundary for 3-d wing
          kk = i
          if ( i.lt.1  ) kk =  1
          if ( i.gt.il ) kk = il
          ip = il + 1 - kk
          do ii = 1,blen
            iin  = in(ii)
            iout = out(ii)
             rho(iout) = q(ip,iin,k,1)
            rhou(iout) = q(ip,iin,k,2)
            rhov(iout) = q(ip,iin,k,3)
            rhow(iout) = q(ip,iin,k,4)
            rhoe(iout) = q(ip,iin,k,5)
            if (nl.eq.6) then
              rhotk(iout) = q(ip,iin,k,6)
            else if (nl.eq.7) then
              rhotk(iout) = q(ip,iin,k,6)
              rhotw(iout) = q(ip,iin,k,7)
            end if
          end do
        case(13)               ! special bc for rotor tip clearance
          ip =  itplk(i,k)
          do ii = 1,blen
            iin  = in(ii)
            iout = out(ii)
             rho(iout) = q(ip,iin,k,1)
            rhou(iout) = q(ip,iin,k,2)
            rhov(iout) = q(ip,iin,k,3)
            rhow(iout) = q(ip,iin,k,4)
            rhoe(iout) = q(ip,iin,k,5)
            if (nl.eq.6) then
              rhotk(iout) = q(ip,iin,k,6)
            else if (nl.eq.7) then
              rhotk(iout) = q(ip,iin,k,6)
              rhotw(iout) = q(ip,iin,k,7)
            end if
          end do
        case(3)                ! no slip adiabatic boundary
          if (moving.eq.0) then   ! stationary wall or rotating coordinates
            if ( num.eq.1 ) then
              y2db = 2.d0*vbds(2,1)
              z2db = 2.d0*vbds(3,1)
            else
              y2db = 2.d0*vbds(2,2)
              z2db = 2.d0*vbds(3,2)
            end if
            iin = in(1)
            do ii = 1,blen
              iout = out(ii)
               rho(iout) = rho(iin)
              rhou(iout) = -rhou(iin)
              rhov(iout) = y2db*rho(iin)-rhov(iin)
              rhow(iout) = z2db*rho(iin)-rhow(iin)
              rhoe(iout) = rhoe(iin)
              if (nl.eq.6) then
                rhotk(iout) = -rhotk(iin)
              else if(nl.eq.7) then
                rhotk(iout) = -rhotk(iin)
                rhotw(iout) = -rhotw(iin)
              end if
            end do
          else                ! moving bound, not available for kw model
            if ( num.eq.1 ) then
              dxp = ( x(i+1,1,k)-x(i,1,k) )*
     >              ( y(i,2,k)-y(i,1,k)+y(i+1,2,k)-y(i+1,1,k) )
              dyp = ( y(i+1,1,k)-y(i,1,k) )*
     >              ( x(i,2,k)-x(i,1,k)+x(i+1,2,k)-x(i+1,1,k) )
              coef = 2.d0/( dxp-dyp )
              dxp = -coef*( y(i+1,1,k)-y(i,1,k) )
              dyp =  coef*( x(i+1,1,k)-x(i,1,k) )
            else
              dxp = (x(i+1,jl+1,k)-x(i,jl+1,k) )*
     >              (y(i,jl+1,k)-y(i,jl,k)+y(i+1,jl+1,k)-y(i+1,jl,k))
              dyp = (y(i+1,jl+1,k)-y(i,jl+1,k) )*
     >              (x(i,jl+1,k)-x(i,jl,k)+x(i+1,jl+1,k)-x(i+1,jl,k))
              coef = 2.d0/( dxp-dyp )
              dxp = -coef*( y(i+1,jl+1,k)-y(i,jl+1,k) )
              dyp =  coef*( x(i+1,jl+1,k)-x(i,jl+1,k) )
            end if
c
            x2db = vbds(1,1) - vbds(1,3)
            y2db = vbds(2,1) - vbds(2,3)
c
            do ii = 1,blen
              iin = in(ii)
              iout = out(ii)
               rho(iout) = rho(iin)
               rhou(iout) = 2.d0*rho(iout)*vbds(1,1)-rhou(iin)
               rhov(iout) = 2.d0*rho(iout)*vbds(2,1)-rhov(iin)
               rhow(iout) =                    rhow(iin)
               uout = rhou(iout)/rho(iout)
               vout = rhov(iout)/rho(iout)
               wout = rhow(iout)/rho(iout)
               qqout = 0.5d0 * (uout**2 + vout**2 + wout**2)
               coef = rho(iout)/( dxp*dxp + dyp*dyp )/tintvl
               pin = gamma1*(qin(5) -.5d0*(rhou(iin)**2+rhov(iin)**2
     $              +rhow(iin)**2)/rho(iin)) ! zero pressure gradient
               pout = pin + coef*(dxp*x2db + dyp*y2db)
               rhoe(iout) = pout/gamma1 + rho(iout)*qqout
            end do
          end if
        case(19)                ! no slip isothermal boundary
          if(moving.eq.0) then 
            do ii = 1,blen
              iin = in(ii)
              iout = out(ii)
              rho(iout) = rho(iin)
              pout = rho(iout)*temp_bnd/(gamma*machinf**2)
              rhou(iout) = -rhou(iin)
              rhov(iout) = -rhov(iin)
              rhow(iout) = -rhow(iin)
              qqout = 0.5*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
     $               /rho(iout)  !is multiplyed by rho(iout)
              rhoe(iout) = pout/gamma1+qqout

c---------------------------no slip isothermal ,
c                           zero normal pressure gradient boundary
              rhou(iout) = -rhou(iin)/rho(iin)
              rhov(iout) = -rhov(iin)/rho(iin)
              rhow(iout) = -rhow(iin)/rho(iin)
              qqout = 0.5*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
              pout = gamma1*(rhoe(iin)-rho(iin)*qqout)
              rho(iout) = gamma*machinf**2*pout/temp_bnd
              rhou(iout) = rhou(iout)*rho(iout)
              rhov(iout) = rhov(iout)*rho(iout)
              rhow(iout) = rhow(iout)*rho(iout)
              rhoe(iout) = pout/gamma1+qqout*rho(iout)
              if(nl.eq.6) then
                rhotk(iout) = -rhotk(iin)
              else if(nl.eq.7) then
                rhotk(iout) = -rhotk(iin)
                rhotw(iout) = -rhotw(iin)
              end if
            end do
          else    
            write(*,*)'not in use'
          end if
        case(23)     ! periodical boundary for o mesh
          iii = il+1-i
          if ( i.ge.102 ) dag = -dag
c         if ( i.ge.52 ) dag = -dag
          do ii = 1,blen
            iout = out(ii)
            iin  = in(ii)
            rg    = q(iii,iin,k,1)
            rginv = 1.0/rg
            vg    = rginv*q(iii,iin,k,3)
            wg    = rginv*q(iii,iin,k,4)
c
            rho(iout)  = rg
            rhou(iout) = q(iii,iin,k,2)
            rhoe(iout) = q(iii,iin,k,5)
            y2db = y(iii,iin+1,k)+y(iii,iin+1,k+1)+
     >             y(iii+1,iin+1,k)+y(iii+1,iin+1,k+1)+
     >             y(iii,iin,k)+y(iii,iin,k+1)+
     >             y(iii+1,iin,k)+y(iii+1,iin,k+1)
            z2db = z(iii,iin+1,k)+z(iii,iin+1,k+1)+
     >             z(iii+1,iin+1,k)+z(iii+1,iin+1,k+1)+
     >             z(iii,iin,k)+z(iii,iin,k+1)+
     >             z(iii+1,iin,k)+z(iii+1,iin,k+1)
            call angle(y2db, z2db, ags)
            utmp =  vg*dcos(ags)+wg*dsin(ags)
            vtmp = -vg*dsin(ags)+wg*dcos(ags)
            ags = ags+dag
            vout = utmp*dcos(ags)-vtmp*dsin(ags)
            wout = utmp*dsin(ags)+vtmp*dcos(ags)
            rhov(out) = rho(out)*vout
            rhow(out) = rho(out)*wout
c
            if (nl.eq.6) then
              rhotk(iout) = q(iii,iin,k,6)
            else if (nl.eq.7) then
              rhotk(iout) = q(iii,iin,k,6)
              rhotw(iout) = q(iii,iin,k,7)
            end if
          end do
        case(4)                ! inlet BC with w = 0  (temporary BC)
          select case (index)
          case(1)
            if(num.eq.1) then
              x0(1,1)=x(1,j+1,k)-x(1,j,k+1)
              x0(2,1)=y(1,j+1,k)-y(1,j,k+1)
              x0(3,1)=z(1,j+1,k)-z(1,j,k+1)
              x0(1,2)=x(1,j+1,k+1)-x(1,j,k)
              x0(2,2)=y(1,j+1,k+1)-y(1,j,k)
              x0(3,2)=z(1,j+1,k+1)-z(1,j,k)
            else
              x0(1,2)=x(il+1,j+1,k)-x(il+1,j,k+1)
              x0(2,2)=y(il+1,j+1,k)-y(il+1,j,k+1)
              x0(3,2)=z(il+1,j+1,k)-z(il+1,j,k+1)
              x0(1,1)=x(il+1,j+1,k+1)-x(il+1,j,k)
              x0(2,1)=y(il+1,j+1,k+1)-y(il+1,j,k)
              x0(3,1)=z(il+1,j+1,k+1)-z(il+1,j,k)
            end if
          case(2)
            if(num.eq.1) then
              x0(1,1)=x(i,1,k+1)-x(i+1,1,k)
              x0(2,1)=y(i,1,k+1)-y(i+1,1,k)
              x0(3,1)=z(i,1,k+1)-z(i+1,1,k)
              x0(1,2)=x(i+1,1,k+1)-x(i,1,k)
              x0(2,2)=y(i+1,1,k+1)-y(i,1,k)
              x0(3,2)=z(i+1,1,k+1)-z(i,1,k)
            else
              x0(1,2)=x(i,jl+1,k+1)-x(i+1,jl+1,k)
              x0(2,2)=y(i,jl+1,k+1)-y(i+1,jl+1,k)
              x0(3,2)=z(i,jl+1,k+1)-z(i+1,jl+1,k)
              x0(1,1)=x(i+1,jl+1,k+1)-x(i,jl+1,k)
              x0(2,1)=y(i+1,jl+1,k+1)-y(i,jl+1,k)
              x0(3,1)=z(i+1,jl+1,k+1)-z(i,jl+1,k)
            end if
          case(3)
            if(num.eq.1) then
              x0(1,1)=x(i+1,j,1)-x(i,j+1,1)
              x0(2,1)=y(i+1,j,1)-y(i,j+1,1)
              x0(3,1)=z(i+1,j,1)-z(i,j+1,1)
              x0(1,2)=x(i+1,j+1,1)-x(i,j,1)
              x0(2,2)=y(i+1,j+1,1)-y(i,j,1)
              x0(3,2)=z(i+1,j+1,1)-z(i,j,1)
            else
              x0(1,2)=x(i+1,j,kl+1)-x(i,j+1,kl+1)
              x0(2,2)=y(i+1,j,kl+1)-y(i,j+1,kl+1)
              x0(3,2)=z(i+1,j,kl+1)-z(i,j+1,kl+1)
              x0(1,1)=x(i+1,j+1,kl+1)-x(i,j,kl+1)
              x0(2,1)=y(i+1,j+1,kl+1)-y(i,j,kl+1)
              x0(3,1)=z(i+1,j+1,kl+1)-z(i,j,kl+1)
            end if
          end select

          x1(1)=x0(2,1)*x0(3,2)-x0(2,2)*x0(3,1)
          x1(2)=x0(1,2)*x0(3,1)-x0(1,1)*x0(3,2)
          x1(3)=x0(1,1)*x0(2,2)-x0(1,2)*x0(2,1)
          b0 = dsqrt(x1(1)**2+x1(2)**2+x1(3)**2)
          x2(1) = x1(1)/b0
          x2(2) = x1(2)/b0
          x2(3) = x1(3)/b0
c
          iin = in(1)
          b0 = dsqrt(rhou(iin)**2+rhov(iin)**2+rhow(iin)**2)/rho(iin)
          uout = b0*x2(1)
          vout = b0*x2(2)
          wout = 0.d0
          qqout = 0.5d0*(uout**2+vout**2+wout**2)
          pout = ptotal/(ttotal/
     $         (ttotal-(gamma-1.d0)*machinf**2*qqout))
     $         **(gamma/(gamma-1.d0))
c
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = gamma*machinf**2*pout/
     $         (ttotal-(gamma-1.d0)*machinf**2*qqout)
            rhou(iout) = rho(iout)*uout
            rhov(iout) = rho(iout)*vout
            rhow(iout) = rho(iout)*wout
            rhoe(iout) = rho(iout)*
     $         ((ttotal-(gamma-1.d0)*machinf**2*qqout)/
     $         (gamma*(gamma-1.d0)*machinf**2)+qqout)
          end do
c
        case(5)       ! subsonic outflow (fixed static pressure at freestream)
          iin = in(1)
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = rho(iin)
            rhou(iout) = rhou(iin)
            rhov(iout) = rhov(iin)
            rhow(iout) = rhow(iin)
            if ( dabs(ronum).gt.1e-9 ) then
              psi = vbds(2,2)*vbds(2,2)+vbds(3,2)*vbds(3,2)
              theta = 0.5d0*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
     >                /rho(iout)/rho(iout)-0.5d0*psi
              rhoe(iout) = poutlet/gamma1 + rho(iout)*theta
            else
              rhoe(iout) = poutlet/gamma1 + 0.5d0*
     $            (rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)/rho(iout)
            end if
            if(nl.eq.6) then
              rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
              if(ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
            end if
          end do
        case(22)       ! subsonic outflow (fixed static pressure at freestream)
          iin = in(1)
          rd =  rho(iin)
          ud = rhou(iin)/rd
          vd = rhov(iin)/rd
          wd = rhow(iin)/rd
          qq = 0.5d0*( ud*ud+vd*vd+wd*wd )
          psi = vbds(2,2)*vbds(2,2)+vbds(3,2)*vbds(3,2)
          theta = qq-0.5d0*psi
          pin = gamma1*( rhoe(iin)-rd*theta )
          as2 = gamma*pin/rd
          cd = dsqrt( as2 )
c
          call pbar(il, jl, kl, nl, x, y, z, q,
     >              ilower, iupper, jlower, jupper, klower, kupper,
     >              gamma, ronum, p_bar)
c
          pout = pin-p_bar+poutlet
          dp = pout-pin
          coef = rd*cd
          ub = ud-dp/coef
          vb = vd
          wb = wd
          qq = 0.5d0*( ub*ub+vb*vb+wb*wb )
          theta = qq-0.5d0*psi
          rhob  = rd+dp/as2
          rhoub = rhob*ub
          rhovb = rhob*vb
          rhowb = rhob*wb
          rhoeb = pout/gamma1+rhob*theta
c
          do ii = 1,blen
            iout = out(ii)
             rho(iout) =  rhob
            rhou(iout) = rhoub
            rhov(iout) = rhovb
            rhow(iout) = rhowb
            rhoe(iout) = rhoeb
c
            if(nl.eq.6) then
              rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
              if(ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
            end if
          end do
        case(6)       ! subsonic inflow (zero gradient in p)
          do ii = 1,blen
            iout = out(ii)
            iin = iout-pdir
            rho(iout) = qout(1)
            rhou(iout) = qout(2)
            rhov(iout) = qout(3)
            rhow(iout) = qout(4)
            pout = gamma1*(rhoe(iin)-.5d0*
     $         (rhou(iin)**2+rhov(iin)**2+rhow(iin)**2)/rho(iin))
            rhoe(iout) = pout/gamma1 + .5d0*
     $         (rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)/rho(iout)
            if(nl.eq.6) then
              rhotk(iout) = qout(6)
            else if(nl.eq.7) then
              rhotk(iout) = qout(6)
              rhotw(iout) = qout(7)
              if(ke) then
                pout = gamma1*(qin(5)-rhotk(iin) -.5d0*
     $            (rhou(iin)**2+rhov(iin)**2+rhow(iin)**2)/rho(iin))
                rhoe(iout) = pout/gamma1 + rhotk(iout)+.5d0*
     $            (rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)/rho(iout)
              end if
            end if
          end do
        case(21)       ! subsonic inflow for rotor
          ii = i
          kk = k
          if ( i.lt.1 )  ii =  1
          if ( i.gt.il ) ii = il
          if ( k.lt.1 )  kk =  1
          if ( k.gt.kl ) kk = kl
c
          call yzsfj(il, jl, kl, x, y, z, ilower, iupper,
     >               jlower, jupper, klower, kupper,
     >               ii, jl+1, kk, yp, zp)
          rad = dsqrt( yp*yp+zp*zp )
          ttl = ttb(ttotal,rad)
          ptl = ptb(ptotal,rad)
          call angle(yp, zp, ags)
c
          iin = in(1)
          rd = rho(iin)
          ud = rhou(iin)/rd
          vd = rhov(iin)/rd
          wd = rhow(iin)/rd
          qq = 0.5d0*( ud*ud+vd*vd+wd*wd )
          psi = vbds(2,2)*vbds(2,2)+vbds(3,2)*vbds(3,2)
          theta = qq-0.5d0*psi
          pin = gamma1*( rhoe(iin)-rd*theta )
          as2 = gamma*pin/rd
          cd = dsqrt( as2 )
          crm = ud-2.d0*cd/gamma1
          vd = vd+vbds(2,2)
          wd = wd+vbds(3,2)
          qq = 0.5d0*( ud*ud+vd*vd+wd*wd )
          co2 = as2+gamma1*qq
          gammap = gamma+1.d0
          coef = dsqrt( gammap*co2/gamma1-0.5d0*gamma1*crm*crm )
          cb = gamma1/gammap*( coef-crm )
c
          if ( cb.lt.0.d0 ) stop
          tb = ttl*( cb*cb/co2 )
          pb = ptl*( tb/ttl )**3.5
          rhob = gamma*pb/as2
c         rhob = pb*gamma*machinf*machinf/tb
          dx1 = angr(il, jl, kl, ii, kk, x, y, z,
     >               ilower, iupper, jlower, jupper, klower, kupper)
          k1 = dx1*dcos(ags)
          k2 = dx1*dsin(ags)
          coef = 1.d0+k1*k1+k2*k2
c
          ub = dsqrt( 2.d0*( ttl-tb )/gamma1/coef )/machinf
          vb = k1*ub
          wb = k2*ub
          vb = vb-vbds(2,2)
          wb = wb-vbds(3,2)
          qq = 0.5d0*( ub*ub+vb*vb+wb*wb )
          theta = qq-0.5d0*psi
          rhoub = rhob*ub
          rhovb = rhob*vb
          rhowb = rhob*wb
          rhoeb = pb/gamma1+rhob*theta
c
          do ii = 1, blen
            iout = out(ii)
             rho(iout) =  rhob
            rhou(iout) = rhoub
            rhov(iout) = rhovb
            rhow(iout) = rhowb
            rhoe(iout) = rhoeb
            if (nl.eq.6) then
              rhotk(iout) = qout(6)
            else if (nl.eq.7) then
              rhotk(iout) = qout(6)
              rhotw(iout) = qout(7)
              if (ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
            end if
          end do
        case(7)       ! inner boundary for mpi communication
        case(100)     ! specified boundary
        case(8)       ! symmetry boundary
          do ii = 1,blen
            iin = in(ii)
            iout = out(ii)
            rho(iout) = rho(iin)
            rhoe(iout) = rhoe(iin)
            if(nl.eq.6) then
               rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
               rhotk(iout) = rhotk(iin)
               rhotw(iout) = rhotw(iin)
            end if
            select case (index)
            case(1)
              jp = j + 1
              kp = k + 1
              select case (num)
              case(1)
                dx1 = x(1,jp,kp) - x(1,j,k)
                dy1 = y(1,jp,kp) - y(1,j,k)
                dz1 = z(1,jp,kp) - z(1,j,k)
                dx2 = x(1,j,kp)  - x(1,jp,k)
                dy2 = y(1,j,kp)  - y(1,jp,k)
                dz2 = z(1,j,kp)  - z(1,jp,k)
              case(2)
                dx1 = x(il+1,jp,kp) - x(il+1,j,k)
                dy1 = y(il+1,jp,kp) - y(il+1,j,k)
                dz1 = z(il+1,jp,kp) - z(il+1,j,k)
                dx2 = x(il+1,j,kp)  - x(il+1,jp,k)
                dy2 = y(il+1,j,kp)  - y(il+1,jp,k)
                dz2 = z(il+1,j,kp)  - z(il+1,jp,k)
              end select
            case(2)
              ip = i + 1
              kp = k + 1
              select case (num)
              case(1)
                dx1 = x(ip,1,kp) - x(i,1,k)
                dy1 = y(ip,1,kp) - y(i,1,k)
                dz1 = z(ip,1,kp) - z(i,1,k)
                dx2 = x(ip,1,k)  - x(i,1,kp)
                dy2 = y(ip,1,k)  - y(i,1,kp)
                dz2 = z(ip,1,k)  - z(i,1,kp)
              case(2)
                dx1 = x(ip,jl+1,kp) - x(i,jl+1,k)
                dy1 = y(ip,jl+1,kp) - y(i,jl+1,k)
                dz1 = z(ip,jl+1,kp) - z(i,jl+1,k)
                dx2 = x(ip,jl+1,k)  - x(i,jl+1,kp)
                dy2 = y(ip,jl+1,k)  - y(i,jl+1,kp)
                dz2 = z(ip,jl+1,k)  - z(i,jl+1,kp)
              end select
            case(3)
              ip = i + 1
              jp = j + 1
              select case (num)
              case(1)
                dx1 = x(ip,jp,1) - x(i,j,1)
                dy1 = y(ip,jp,1) - y(i,j,1)
                dz1 = z(ip,jp,1) - z(i,j,1)
                dx2 = x(i,jp,1)  - x(ip,j,1)
                dy2 = y(i,jp,1)  - y(ip,j,1)
                dz2 = z(i,jp,1)  - z(ip,j,1)
              case(2)
                dx1 = x(ip,jp,kl+1) - x(i,j,kl+1)
                dy1 = y(ip,jp,kl+1) - y(i,j,kl+1)
                dz1 = z(ip,jp,kl+1) - z(i,j,kl+1)
                dx2 = x(i,jp,kl+1)  - x(ip,j,kl+1)
                dy2 = y(i,jp,kl+1)  - y(ip,j,kl+1)
                dz2 = z(i,jp,kl+1)  - z(ip,j,kl+1)
              end select
            end select
                                                                                
            xp = dabs(dy1 * dz2 - dz1 * dy2)
            yp = dabs(dz1 * dx2 - dx1 * dz2)
            zp = dabs(dx1 * dy2 - dy1 * dx2)
            dp = dmax1(xp,yp,zp)

            if (dp.eq.xp) then
c *** symmetry about x-plane (x = cont)
              rhou(iout) = -rhou(iin)
              rhov(iout) = rhov(iin)
              rhow(iout) = rhow(iin)
            else if (dp.eq.yp) then
c *** symmetry about y-plane (y = cont)
              rhou(iout) = rhou(iin)
              rhov(iout) = -rhov(iin)
              rhow(iout) = rhow(iin)
            else if (dp.eq.zp) then
c *** symmetry about z-plane (z = cont)
              rhou(iout) = rhou(iin)
              rhov(iout) = rhov(iin)
              rhow(iout) = -rhow(iin)
            end if
          end do
        case(9)         ! fixed total pressure and temperature at subsonic inlet
          extra = 0           ! 0: vel, 1: p extrapolated
          k1 = dtan(angl1)    ! angle between first and main direction
          k2 = dtan(angl2)    ! angle between second and main direction
                              ! main direction is one of u, v, and w
                              ! first and second direction is the first and
                              ! second direction in the order of u, v, and w
                              ! except the selected main direction
          if(extra.eq.0) then
            do ii = 1,blen
              iout = out(ii)
              iin = iout-pdir
              select case (main_dir)
              case(1)
                uout = rhou(iin)/rho(iin)
                vout = uout*k1
                wout = uout*k2
              case(2)
                vout = rhov(iin)/rho(iin)
                uout = vout*k1
                wout = vout*k2
              case(3)
                wout = rhow(iin)/rho(iin)
                uout = wout*k1
                vout = wout*k2
              end select

              qqout = 0.5d0*(uout**2+vout**2+wout**2)
              pout = ptotal/(ttotal/
     $             (ttotal-(gamma-1.d0)*machinf**2*qqout))
     $             **(gamma/(gamma-1.d0))
              rho(iout) = gamma*machinf**2*pout/
     $             (ttotal-(gamma-1.d0)*machinf**2*qqout)
              rhou(iout) = rho(iout)*uout
              rhov(iout) = rho(iout)*vout
              rhow(iout) = rho(iout)*wout
              rhoe(iout) = rho(iout)*
     $             ((ttotal-(gamma-1.d0)*machinf**2*qqout)/
     $             (gamma*(gamma-1.d0)*machinf**2)+qqout)
              if(nl.eq.6) then
                rhotk(iout) = rho(iout)*qout(6)/qout(1)
              else if(nl.eq.7) then
                rhotk(iout) = rho(iout)*qout(6)/qout(1)
                rhotw(iout) = rho(iout)*qout(7)/qout(1)
                if(ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
              end if
            end do
          else if(extra.eq.1) then
            do ii = 1,blen
              iout = out(ii)
              iin = iout-pdir
              pout = gamma1*(qin(5) -.5d0*
     $             (rhou(iin)**2+rhov(iin)**2+rhow(iin)**2)/rho(iin))
              pratio = (ptotal/pout)**(gamma1/gamma) ! = ttotal/t
              mach = dsqrt(2.d0/gamma1*(pratio-1.d0)) ! local mach number
              tout = ttotal/pratio
              kcoef = dsqrt(1.d0+k1*k1+k2*k2)
              main_vel = (1.d0/kcoef)*mach*dsqrt(tout)/machinf
              select case (main_dir)
              case(1)
                uout = main_vel
                vout = main_vel*k1
                wout = main_vel*k2
              case(2)
                vout = main_vel
                uout = main_vel*k1
                wout = main_vel*k2
              case(3)
                wout = main_vel
                uout = main_vel*k1
                vout = main_vel*k2
              end select
              qqout = .5d0*(main_vel*kcoef)**2
              rho(out) = gamma*machinf*machinf*pout/tout
              rhou(iout) = rho(iout)*uout
              rhov(iout) = rho(iout)*vout
              rhow(iout) = rho(iout)*wout
              rhoe(iout) = rho(iout)*
     $             (tout/(gamma*gamma1*machinf*machinf)+qqout)
            end do
          end if
        case(10,20)               ! periodical boundary
          do ii = 1,blen
            iout = out(ii)
            iin = bndp+pdir*(ii-1)
            rho(iout)=rho(iin)
            rhou(iout)=rhou(iin)
            rhov(iout)=rhov(iin)
            rhow(iout)=rhow(iin)
            rhoe(iout)=rhoe(iin)
            if(nl.eq.6) then
              rhotk(iout) = rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
            end if
          end do
        case(11)                ! subsonic outflow (fixed static pressure)
          pout = gamma1*(qout(5) -.5d0*
     $           (qout(2)**2+qout(3)**2+qout(4)**2)/qout(1))
          iin = in(1)
          do ii = 1,blen
            iout = out(ii)
            rho(iout) = rho(iin)
            rhou(iout) = rhou(iin)
            rhov(iout) = rhov(iin)
            rhow(iout) = rhow(iin)
            rhoe(iout) = pout/gamma1 + .5d0*
     $        (rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)/rho(iout)
            if(nl.eq.6) then
              rhotk(iout) = -rhotk(iin)
            else if(nl.eq.7) then
              rhotk(iout) = rhotk(iin)
              rhotw(iout) = rhotw(iin)
              if(ke) rhoe(iout) = rhoe(iout)+rhotk(iout)
            end if
          end do
c
        case(14)            ! non-reflective bc for solid wall
c             
          do ii =1,blen
            iout = out(ii)
            iin = iout-pdir
            uout = -rhou(iin)/rho(iin)
            vout = -rhov(iin)/rho(iin)
            wout = -rhow(iin)/rho(iin)
            qqout = 0.5d0*( uout*uout+vout*vout+wout*wout )
            rhoe(iout) = qout(5)
            pin = gamma1*( rhoe(iin)-rho(iin)*qqout )
            coef = rho(iin)/pin
            dxp = 1.d0/gamma1 + coef*qqout
            pout = rhoe(iout)/dxp
            rho(iout) = coef*pout
            rhou(iout) = rho(iout)*uout
            rhov(iout) = rho(iout)*vout
            rhow(iout) = rho(iout)*wout
c             
            if ( index.eq.2 ) then
              q(i,iout,k,1) = rho(iout)
              q(i,iout,k,2) = rhou(iout)
              q(i,iout,k,3) = rhov(iout)
              q(i,iout,k,4) = rhow(iout)
            end if
          end do
c
        case(16)            ! non-reflective bc ( inflow )
c             
          k1 = dtan(angl1)
          k2 = dtan(angl2)
          do ii =1,blen
            iout = out(ii)
            iin = iout-pdir
            uout = rhou(iin)/rho(iin)
            vout = k1*uout
            wout = k2*uout
            qqout = 0.5d0*( uout*uout+vout*vout+wout*wout )
            rho(iout)  = rho(iin)
            rhoe(iout) = qout(5)
            pout = gamma1*( rhoe(iout)-rho(iout)*qqout )
            rhou(iout) = rho(iout)*uout
            rhov(iout) = rho(iout)*vout
            rhow(iout) = rho(iout)*wout
          end do

        case(101,102,103,104,105,106,107,108,109,110)                ! no slip isothermal boundary
          if(moving.eq.0) then 
            do ii = 1,blen
              iin = in(ii)
              iout = out(ii)
c              rho(iout) = rho(iin)
c              pout = rho(iout)*control%t(bcindex-100)/(gamma*machinf**2)
c              rhou(iout) = -rhou(iin)
c              rhov(iout) = -rhov(iin)
c              rhow(iout) = -rhow(iin)
c              qqout = 0.5*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
c     $               /rho(iout)  !is multiplyed by rho(iout)
c              rhoe(iout) = pout/gamma1+qqout

c---------------------------no slip isothermal ,
c                           zero normal pressure gradient boundary
              rhou(iout) = -rhou(iin)/rho(iin)
              rhov(iout) = -rhov(iin)/rho(iin)
              rhow(iout) = -rhow(iin)/rho(iin)
              qqout = 0.5*(rhou(iout)**2+rhov(iout)**2+rhow(iout)**2)
              pout = gamma1*(rhoe(iin)-rho(iin)*qqout)
              rho(iout) = gamma*machinf**2*pout/control%t(bcindex-100)
              rhou(iout) = rhou(iout)*rho(iout)
              rhov(iout) = rhov(iout)*rho(iout)
              rhow(iout) = rhow(iout)*rho(iout)
              rhoe(iout) = pout/gamma1+qqout*rho(iout)
              if(nl.eq.6) then
                rhotk(iout) = -rhotk(iin)
              else if(nl.eq.7) then
                rhotk(iout) = -rhotk(iin)
                rhotw(iout) = -rhotw(iin)
              end if
            end do
          else    
            write(*,*)'not in use'
          end if
c
        end select
        select case (index)
        case(1)
          do ii = 1,blen
            iout = out(ii)
            q(iout,j,k,1) = rho(iout)
            q(iout,j,k,2) = rhou(iout)
            q(iout,j,k,3) = rhov(iout)
            q(iout,j,k,4) = rhow(iout)
            q(iout,j,k,5) = rhoe(iout)
            if (nl.eq.6) then
              q(iout,j,k,6) = rhotk(iout)
            end if
          end do
        case(2)
          do ii = 1,blen
            iout = out(ii)
            q(i,iout,k,1) = rho(iout)
            q(i,iout,k,2) = rhou(iout)
            q(i,iout,k,3) = rhov(iout)
            q(i,iout,k,4) = rhow(iout)
            q(i,iout,k,5) = rhoe(iout)
            if (nl.eq.6) then
              q(i,iout,k,6) = rhotk(iout)
            end if
          end do
        case(3)
          do ii = 1,blen
            iout = out(ii)
            q(i,j,iout,1) = rho(iout)
            q(i,j,iout,2) = rhou(iout)
            q(i,j,iout,3) = rhov(iout)
            q(i,j,iout,4) = rhow(iout)
            q(i,j,iout,5) = rhoe(iout)
            if (nl.eq.6) then
              q(i,j,iout,6) = rhotk(iout)
            end if
          end do
        end select
        end if
c---------------end of precondition-------------
      end do
c
      return
      end

