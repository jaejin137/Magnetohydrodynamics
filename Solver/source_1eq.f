      subroutine source_1eq(blnu, il, jl, kl, nl,
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     x, y, z, q, bc_xie_lower, bc_xie_upper,
     $     bc_eta_lower, bc_eta_upper, bc_zta_lower, bc_zta_upper,
     $     vol, dist,
     $     iblnu, ipt, jpt, kpt, ic1, ic2, ic3, tko, cb1, cb2,
     $     cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4, source,
     $     control)
c
c     compute source term for SA_1eq model
c
      use datain

      implicit none

      type (datain_type), intent(in)::control
c
c     include header file
c
      include "/usr/local/include/mpif.h"
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     blnu, il, jl, kl, nl

      integer, intent(in)::
     $     bc_xie_lower(jl, kl), bc_xie_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zta_lower(il, jl), bc_zta_upper(il, jl)

      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)

      double precision, intent(in), dimension(ilower+1:iupper,
     $     jlower+1:jupper, klower+1:kupper)::
     $     x, y, z,
     $     vol(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)

      double precision, dimension(il,jl,kl), intent(out)::
     $     source               ! source terms of SA 1-eq model

      double precision::
     $     reynolds, tref, gamma, machinf, gamma1

      logical:: suther
c
c     SA 1eq turbulent model contants
c
      double precision, intent(in)::
     $  tko, cb1, cb2, cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4
      integer, intent(in):: iblnu, ipt, jpt, kpt, ic1, ic2, ic3

      double precision, intent(in):: dist(il,jl,kl)
c
c     LOCAL VARIABLES
c
      double precision,
     $     dimension(il+1,jl+1,kl+1)::
     $     lx, ly, lz, mx, my, mz, nx, ny, nz

      double precision, dimension(ilower:iupper, jlower:jupper,
     $     klower:kupper)::
     $     rho,                 ! density
     $     u, v, w,             ! velocity
     $     tk                   ! k, omega

      double precision::
     $     term1,term2,term3,term4,term5
      double precision::
     $     drdxie,drdeta,drdzta,drdx,drdy,drdz,
     $     dtkdxie,dtkdeta,dtkdzta,dtkdx,dtkdy,dtkdz,
     $     cap_x,fv1,fv2,ft1,ft2,g,r,fw,
     $     delt_u,vislam,t,qq,s_bar,gt,omegat,
     $     dx, dy, dz, dd, dtr, xptr(5)

      double precision, dimension(il,jl,kl)::
     $     omega                ! vorticity

      integer::
     $     ierr, i, j, k, ii, jj, kk
c
c *** START SUBROUTINE
c
c-----------------------------------------
      suther=control%suther

      gamma=control%gamma
      machinf=control%machinf
      reynolds=control%reynolds
      tref=control%tref
c-----------------------------------------
      gamma1=gamma-1.0
c *** compute metric
c
      call metric_all(x, y, z, lx, ly, lz, mx, my, mz, nx, ny, nz,
     $     il, jl, kl, ilower, iupper, jlower, jupper, klower, kupper)
c
c *** compute variable values
c
      do i = ilower, iupper
        do j = jlower, jupper
          do k = klower, kupper
            rho(i,j,k) = q(i,j,k,1)
            u(i,j,k) = q(i,j,k,2)/rho(i,j,k)
            v(i,j,k) = q(i,j,k,3)/rho(i,j,k)
            w(i,j,k) = q(i,j,k,4)/rho(i,j,k)
            tk(i,j,k) = q(i,j,k,6)/rho(i,j,k)
          end do
        end do
      end do

      call vorticity(il, jl, kl, nl, ilower, iupper,
     $     jlower, jupper, klower, kupper, x, y, z, q, vol,
     $     lx, ly, lz, mx, my, mz, nx, ny, nz, u, v, w, omega)
c
c *** trip point
c
      if (blnu.eq.iblnu) then
        xptr(1) = x(ipt,jpt,kpt)
        xptr(2) = y(ipt,jpt,kpt)
        xptr(3) = z(ipt,jpt,kpt)
        xptr(4) = omega(ipt,jpt,kpt)
        dx = x(ipt+ic1,jpt+ic2,kpt+ic3) - x(ipt,jpt,kpt)
        dy = y(ipt+ic1,jpt+ic2,kpt+ic3) - y(ipt,jpt,kpt)
        dz = z(ipt+ic1,jpt+ic2,kpt+ic3) - z(ipt,jpt,kpt)
        xptr(5) = dsqrt(dx**2+dy**2+dz**2)
        call mpi_bcast(xptr, 5, mpi_double_precision, iblnu-1,
     $        mpi_comm_world, ierr)
      else
        call mpi_bcast(xptr, 5, mpi_double_precision, iblnu-1,
     $        mpi_comm_world, ierr)
      end if
      omegat = xptr(4)
      dd = xptr(5)
c
      do i = 1, il
        do j = 1, jl
          do k = 1, kl
            qq = .5d0*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
            t = gamma * gamma1 * machinf * machinf *
     $          (q(i,j,k,5) - rho(i,j,k) * qq) / rho(i,j,k)
            if ( suther ) then
              vislam = dsqrt(t)**3*(1.d0+tref)/(t + tref)
            else
              vislam = t
            end if

            cap_x = q(i,j,k,6)/vislam
            fv1 = cap_x**3/(cap_x**3+cv1**3)
            fv2 = 1.d0-cap_x/(1.d0+cap_x*fv1)
c
c *** trip function
c
            delt_u = dsqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
            dx = x(i,j,k)-xptr(1)
            dy = y(i,j,k)-xptr(2)
            dz = z(i,j,k)-xptr(3)
            dtr = dsqrt(dx**2+dy**2+dz**2)
            gt = min(0.1d0,delt_u/(omegat*dd))
            ft1 = ct1*gt*dexp(-ct2*(omegat/delt_u)**2*
     $            (dist(i,j,k)**2+(dtr*gt)**2))
c
            ft2 = ct3*exp(-ct4*cap_x**2)
            s_bar = omega(i,j,k)+fv2*tk(i,j,k)/
     $              (reynolds*(cap_k*dist(i,j,k))**2)

            term1 = cb1*q(i,j,k,6)*(1.d0-ft2)*s_bar

            dtkdxie = .5d0*(tk(i+1,j,k)-tk(i-1,j,k))
            dtkdeta = .5d0*(tk(i,j+1,k)-tk(i,j-1,k))
            dtkdzta = .5d0*(tk(i,j,k+1)-tk(i,j,k-1))
            dtkdx = (dtkdxie*lx(i,j,k) + dtkdeta*mx(i,j,k)
     $             +dtkdzta*nx(i,j,k))/vol(i,j,k)
            dtkdy = (dtkdxie*ly(i,j,k) + dtkdeta*my(i,j,k)
     $             +dtkdzta*ny(i,j,k))/vol(i,j,k)
            dtkdz = dtkdxie*lz(i,j,k) + dtkdeta*mz(i,j,k)
     $             +dtkdzta*nz(i,j,k)/vol(i,j,k)
            
            drdxie = .5d0*(rho(i+1,j,k)-rho(i-1,j,k))
            drdeta = .5d0*(rho(i,j+1,k)-rho(i,j-1,k))
            drdzta = .5d0*(rho(i,j,k+1)-rho(i,j,k-1))
            drdx = (drdxie*lx(i,j,k) + drdeta*mx(i,j,k)
     $            +drdzta*nx(i,j,k))/vol(i,j,k)
            drdy = (drdxie*ly(i,j,k) + drdeta*my(i,j,k)
     $            +drdzta*ny(i,j,k))/vol(i,j,k)
            drdz = (drdxie*lz(i,j,k) + drdeta*mz(i,j,k)
     $            +drdzta*nz(i,j,k))/vol(i,j,k)

            term2 = -(vislam/rho(i,j,k)+tk(i,j,k))
     $             *(dtkdx*drdx+dtkdy*drdy+dtkdz*drdz)/tko

            term3 = cb2*rho(i,j,k)*(dtkdx**2+dtkdy**2+dtkdz**2)/tko

            r = tk(i,j,k)/
     $          (reynolds*s_bar*(cap_k*dist(i,j,k))**2)
            if (r.gt.10.d0) r = 10.d0
            g = r+cw2*(r**6-r)
            fw = g*((1.d0+cw3**6)/(g**6+cw3**6))**(1.d0/6.d0)

            term4 = -rho(i,j,k)*(cw1*fw-cb1*ft2/cap_k**2)
     $             *(tk(i,j,k)/dist(i,j,k))**2

            term5 = reynolds*rho(i,j,k)*ft1*delt_u**2
            
            source(i,j,k) = term1+
     $                       (term2+term3+term4)/reynolds+term5
               
          end do
        end do
      end do

      end subroutine source_1eq
