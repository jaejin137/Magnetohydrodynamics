      subroutine wall_func(reynolds, machinf, suther, tref, gamma, q, x,
     $     y, z, dir, nl, lx, ly, lz, mx, my, mz, nx, ny, nz, tkw, tww,
     $     yplus, wallfunc, tkb, twb, check, ke)
c
c     wall treatment calculation and its influence on rhs matrix
c     more treatment about lhs matrix is applied in matrix_bnd
c
      implicit none
c
c     INTERFACE VARIABLES

      integer, intent(in)::
     $     nl,                  ! equation number
     $     dir                  ! wall location
     $                          ! 1, xi_lower, 2, xi_upper,
     $                          ! 3, eta_lower, 4, eta_upper,
     $                          ! 5, zeta_lower, 6, zeta_upper,

      double precision, intent(in)::
     $     reynolds,            ! reynolds number
     $     machinf,             ! reference mach number
     $     tref,
     $     gamma,
     $     q(nl),
     $     lx, ly, lz,          ! center metrics
     $     mx, my, mz,
     $     nx, ny, nz,
     $     twb,
     $     tkb


      double precision, dimension(2,2,2), intent(in)::
     $     x, y, z              ! wall cell geometry

      logical, intent(in)::
     $     suther,              ! viscosity method
     $     wallfunc,
     $     ke

      double precision, intent(out)::
     $     tkw,                 ! k value
     $     tww,                 ! omega value
     $     yplus
c
c     LOCAL VARIABLES

      double precision::
     $     niu,                 ! dynamic viscosity
     $     delta,               ! iteration difference
     $     utau,                ! friction velocity
     $     dist,                ! distance from cell center to wall
     $     vel,                 ! velocity parallel to wall
     $     rho,
     $     u, v, w,             ! velocity components
     $     e,                   ! total energy
     $     t,                   ! temperature
     $     tk,                  ! kinetic turbulent energy
     $     x0, y0, z0,          ! cell center location
     $     distoplane,
     $     leninv,
     $     llx, lly, llz,
     $     mmx, mmy, mmz,
     $     nnx, nny, nnz,
     $     resqrt,              ! sqrt(reynolds)
     $     uplus

      double precision, dimension(2,2)::
     $     surfx, surfy, surfz

      double precision, parameter:: ! k omega constants
     $     kappa = .41d0,
     $     b = 5.d0,
     $     eps = 1.d-6

      logical, intent(in)::
     $     check
c
c     SUBROUTINE BEGIN
c
c..   center location

      x0 = sum(x)/size(x)
      y0 = sum(y)/size(y)
      z0 = sum(z)/size(z)

c..   wall surface

      select case (dir)
      case (1)
         surfx = x(1,:,:)
         surfy = y(1,:,:)
         surfz = z(1,:,:)
      case (2)
         surfx = x(2,:,:)
         surfy = y(2,:,:)
         surfz = z(2,:,:)
      case (3)
         surfx = x(:,1,:)
         surfy = y(:,1,:)
         surfz = z(:,1,:)
      case (4)
         surfx = x(:,2,:)
         surfy = y(:,2,:)
         surfz = z(:,2,:)
      case (5)
         surfx = x(:,:,1)
         surfy = y(:,:,1)
         surfz = z(:,:,1)
      case (6)
         surfx = x(:,:,2)
         surfy = y(:,:,2)
         surfz = z(:,:,2)
      end select

c..   distance from center to wall

      dist = distoplane(x0, y0, z0, surfx, surfy, surfz)

c..   variables

      rho = q(1)
      u = q(2)/q(1)
      v = q(3)/q(1)
      w = q(4)/q(1)

      e = q(5)/q(1)
      tk = q(6)/q(1)
      if(ke) then
         t =gamma*(gamma-1.d0)*machinf*machinf*(e-tk-.5d0*(u*u+v*v+w*w))
      else
         t = gamma*(gamma-1.d0)*machinf*machinf*(e-.5d0*(u*u+v*v+w*w))
      end if
c..   dynamic viscosity

      if ( suther ) then
         niu = (dsqrt(t)**3*(1.d0+tref)/(t + tref))/rho
      else
         niu = t/rho
      end if

c..   unit center metrics

      leninv = 1.d0/dsqrt(lx**2+ly**2+lz**2)
      llx = lx*leninv
      lly = ly*leninv
      llz = lz*leninv
      leninv = 1.d0/dsqrt(mx**2+my**2+mz**2)
      mmx = mx*leninv
      mmy = my*leninv
      mmz = mz*leninv
      leninv = 1.d0/dsqrt(nx**2+ny**2+nz**2)
      nnx = nx*leninv
      nny = ny*leninv
      nnz = nz*leninv

c..   velocity parallel to wall

      select case (dir)
      case (1,2)
         vel = dsqrt((u*mmx+v*mmy+w*mmz)**2+(u*nnx+v*nny+w*nnz)**2)
      case (3,4)
         vel = dsqrt((u*llx+v*lly+w*llz)**2+(u*nnx+v*nny+w*nnz)**2)
      case (5,6)
         vel = dsqrt((u*llx+v*lly+w*llz)**2+(u*mmx+v*mmy+w*mmz)**2)
      end select

      if(wallfunc) then
c..   solve utau using Newton-Raphson method

         utau = 1.d-10          ! friction velocity initial
         delta = 1.d0           ! iteration difference of utau
c     
         resqrt = dsqrt(reynolds)

         do while(dabs(delta/utau).gt.eps)
            delta = (vel*utau*resqrt -
     $           utau*utau/kappa*log(utau*dist*resqrt/niu)-b*utau*utau)
     $           /(vel*resqrt+utau/kappa)
            utau = utau+delta

         end do

         tkw = utau*utau/dsqrt(tkb)/reynolds ! k
         tww = dsqrt(tkw)/(tkb**.25d0*kappa*dist) ! omega

         uplus = vel*resqrt/utau
         yplus = dist*utau*resqrt/niu

      else
         tww = 6.d0*niu/(twb*dist*dist*reynolds)
      end if

      end subroutine wall_func

      
