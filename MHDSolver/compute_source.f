      subroutine compute_source(x, y, z, ilower, iupper, jlower, jupper,
     $     klower, kupper, il, jl, kl, nl, q, bc_xi_lower, bc_xi_upper,
     $     bc_eta_lower, bc_eta_upper, bc_zeta_lower, bc_zeta_upper,
     $     vol, visturb, tkb, twb, alpha, source, ke)
c
c     compute turbulence viscosity for k-w model
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     il, jl, kl, nl

      integer, intent(in)::
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl)

      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)

      double precision, intent(in), dimension(ilower+1:iupper,
     $     jlower+1:jupper, klower+1:kupper)::
     $     x, y, z,
     $     vol(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)


      double precision, dimension(0:il+1,0:jl+1,0:kl+1), intent(in)::
     $     visturb

      double precision, dimension(il,jl,kl,6:7), intent(out)::
     $     source               ! source terms of k-omega model

      double precision, intent(in)::
     $     alpha,               ! alpha in omega equation
     $     twb,                 ! in omega equation
     $     tkb                  ! in k equation

      logical, intent(in)::
     $     ke
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
     $     tk, tw               ! k, omega

      double precision::
     $     dudxi=0.d0, dudeta=0.d0, dudzeta=0.d0,
     $     dvdxi=0.d0, dvdeta=0.d0, dvdzeta=0.d0,
     $     dwdxi=0.d0, dwdeta=0.d0, dwdzeta=0.d0,
     $     dudx=0.d0, dudy=0.d0, dudz=0.d0,
     $     dvdx=0.d0, dvdy=0.d0, dvdz=0.d0,
     $     dwdx=0.d0, dwdy=0.d0, dwdz=0.d0,
     $     txx, txy, txz,
     $     tyy, tyz, tzz,
     $     pk

      integer::
     $     i, j, k
      
      logical::
     $     ilow = .False.,
     $     iupp = .False.,          ! indicate if inner bound exist
     $     jlow = .False.,
     $     jupp = .False.,
     $     klow = .False.,
     $     kupp = .False.

      integer, parameter::
     $     bcinnr = 7,          ! inner boundary
     $     bcperi = 10          ! periodical boundary
c
c *** START SUBROUTINE
c
      ! compute metric
      call metric_all(x, y, z, lx, ly, lz, mx, my, mz, nx, ny, nz,
     $     il, jl, kl, ilower, iupper, jlower, jupper, klower, kupper)

      ! compute variable values
      rho = q(:,:,:,1)
      u = q(:,:,:,2)/rho
      v = q(:,:,:,3)/rho
      w = q(:,:,:,4)/rho
      tk = q(:,:,:,6)/rho
      tw = q(:,:,:,7)/rho

      do i = 1, il
         do j = 1, jl
            do k = 1, kl

               ilow = bc_xi_lower(j,k).eq.bcinnr
     $              .or.bc_xi_lower(j,k).eq.bcperi
     $              .or.bc_xi_lower(j,k).eq.20
               iupp = bc_xi_upper(j,k).eq.bcinnr
     $              .or.bc_xi_upper(j,k).eq.bcperi
     $              .or.bc_xi_upper(j,k).eq.20
               jlow = bc_eta_lower(i,k).eq.bcinnr
     $              .or.bc_eta_lower(i,k).eq.bcperi
     $              .or.bc_eta_lower(i,k).eq.20
               jupp = bc_eta_upper(i,k).eq.bcinnr
     $              .or.bc_eta_upper(i,k).eq.bcperi
     $              .or.bc_eta_upper(i,k).eq.20
               klow = bc_zeta_lower(i,j).eq.bcinnr
     $              .or.bc_zeta_lower(i,j).eq.bcperi
     $              .or.bc_zeta_lower(i,j).eq.20
               kupp = bc_zeta_upper(i,j).eq.bcinnr
     $              .or.bc_zeta_upper(i,j).eq.bcperi
     $              .or.bc_zeta_upper(i,j).eq.20

               if(il.gt.1) then
                  if(i.eq.1.and.(.not.ilow)) then
                     dudxi = u(i+1,j,k)-u(i,j,k)
                     dvdxi = v(i+1,j,k)-v(i,j,k)
                     dwdxi = w(i+1,j,k)-w(i,j,k)
                  else if(i.eq.il.and.(.not.iupp)) then
                     dudxi = u(i,j,k)-u(i-1,j,k)
                     dvdxi = v(i,j,k)-v(i-1,j,k)
                     dwdxi = w(i,j,k)-w(i-1,j,k)
                  else
                     dudxi = .5d0*(u(i+1,j,k)-u(i-1,j,k))
                     dvdxi = .5d0*(v(i+1,j,k)-v(i-1,j,k))
                     dwdxi = .5d0*(w(i+1,j,k)-w(i-1,j,k))
                  end if
               end if

               if(jl.gt.1) then
                  if(j.eq.1.and.(.not.jlow)) then
                     dudeta = u(i,j+1,k)-u(i,j,k)
                     dvdeta = v(i,j+1,k)-v(i,j,k)
                     dwdeta = w(i,j+1,k)-w(i,j,k)
                  else if(j.eq.jl.and.(.not.jupp)) then
                     dudeta = u(i,j,k)-u(i,j-1,k)
                     dvdeta = v(i,j,k)-v(i,j-1,k)
                     dwdeta = w(i,j,k)-w(i,j-1,k)
                  else
                     dudeta = .5d0*(u(i,j+1,k)-u(i,j-1,k))
                     dvdeta = .5d0*(v(i,j+1,k)-v(i,j-1,k))
                     dwdeta = .5d0*(w(i,j+1,k)-w(i,j-1,k))
                  end if
               end if

               if(kl.gt.1) then
                  if(k.eq.1.and.(.not.klow)) then
                     dudzeta = u(i,j,k+1)-u(i,j,k)
                     dvdzeta = v(i,j,k+1)-v(i,j,k)
                     dwdzeta = w(i,j,k+1)-w(i,j,k)
                  else if(k.eq.kl.and.(.not.kupp)) then
                     dudzeta = u(i,j,k)-u(i,j,k-1)
                     dvdzeta = v(i,j,k)-v(i,j,k-1)
                     dwdzeta = w(i,j,k)-w(i,j,k-1)
                  else
                     dudzeta = .5d0*(u(i,j,k+1)-u(i,j,k-1))
                     dvdzeta = .5d0*(v(i,j,k+1)-v(i,j,k-1))
                     dwdzeta = .5d0*(w(i,j,k+1)-w(i,j,k-1))
                  end if
               end if

               dudx = (dudxi*lx(i,j,k) + dudeta*mx(i,j,k)
     $              + dudzeta*nx(i,j,k))/vol(i,j,k)
               dudy = (dudxi*ly(i,j,k) + dudeta*my(i,j,k)
     $              + dudzeta*ny(i,j,k))/vol(i,j,k)
               dudz = (dudxi*lz(i,j,k) + dudeta*mz(i,j,k)
     $              + dudzeta*nz(i,j,k))/vol(i,j,k)

               dvdx = (dvdxi*lx(i,j,k) + dvdeta*mx(i,j,k)
     $              + dvdzeta*nx(i,j,k))/vol(i,j,k)
               dvdy = (dvdxi*ly(i,j,k) + dvdeta*my(i,j,k)
     $              + dvdzeta*ny(i,j,k))/vol(i,j,k)
               dvdz = (dvdxi*lz(i,j,k) + dvdeta*mz(i,j,k)
     $              + dvdzeta*nz(i,j,k))/vol(i,j,k)

               dwdx = (dwdxi*lx(i,j,k) + dwdeta*mx(i,j,k)
     $              + dwdzeta*nx(i,j,k))/vol(i,j,k)
               dwdy = (dwdxi*ly(i,j,k) + dwdeta*my(i,j,k)
     $              + dwdzeta*ny(i,j,k))/vol(i,j,k)
               dwdz = (dwdxi*lz(i,j,k) + dwdeta*mz(i,j,k)
     $              + dwdzeta*nz(i,j,k))/vol(i,j,k)

               if(ke) then
                  txx = visturb(i,j,k)*(4.d0*dudx-2.d0*dvdy-2.d0*dwdz)
     $                 /3.d0
     $                 -2.d0*rho(i,j,k)*tk(i,j,k)/3.d0
                  tyy = visturb(i,j,k)*(4.d0*dvdy-2.d0*dudx-2.d0*dwdz)
     $                 /3.d0
     $                 -2.d0*rho(i,j,k)*tk(i,j,k)/3.d0
                  tzz = visturb(i,j,k)*(4.d0*dwdz-2.d0*dudx-2.d0*dvdy)
     $                 /3.d0
     $                 -2.d0*rho(i,j,k)*tk(i,j,k)/3.d0
               else
                  txx = visturb(i,j,k)*(4.d0*dudx-2.d0*dvdy-2.d0*dwdz)
     $                 /3.d0
                  tyy = visturb(i,j,k)*(4.d0*dvdy-2.d0*dudx-2.d0*dwdz)
     $                 /3.d0
                  tzz = visturb(i,j,k)*(4.d0*dwdz-2.d0*dudx-2.d0*dvdy)
     $                 /3.d0
               end if

               txy = visturb(i,j,k)*(dudy+dvdx)
               txz = visturb(i,j,k)*(dudz+dwdx)
               tyz = visturb(i,j,k)*(dvdz+dwdy)

               pk = txx*dudx+tyy*dvdy+tzz*dwdz
     $              +txy*(dudy+dvdx)+txz*(dudz+dwdx)+tyz*(dvdz+dwdy)

               source(i,j,k,6) = pk-tkb*rho(i,j,k)*tk(i,j,k)*tw(i,j,k)
               source(i,j,k,7) = tw(i,j,k)*(alpha*pk/tk(i,j,k)-
     $              twb*rho(i,j,k)*tw(i,j,k))
               
            end do
         end do
      end do

      end subroutine compute_source
