      subroutine metric(index, i1, i2, imax, il, jl, kl, x, y, z,
     $     lx, ly, lz, mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat,
     $     mxhat, myhat, mzhat, nxhat, nyhat, nzhat, ilower, iupper,
     $     jlower, jupper, klower, kupper, bclower, bcupper)
c     
c     DEFINITION
c
c     compute metrics on xi, eta or zeta faces
c
c     index   quantity computed               definition
c
c       1     lx, ly, lz             area-weighted normal to xi-face
c             lxhat, lyhat, lzhat    unit normal to xi-face
c             mxhat, myhat, mzhat    ARBITRARY unit normal to (lx,ly,lz)
c             nxhat, nyhat, nzhat    ARBITRARY unit normal to (lx,ly,lz)
c
c      -1     lx, ly, lz             area-weighted normal to xi-face
c             mx, my, mz             AVERAGE value of (mx,my,mz) on xi-face
c             nx, ny, nz             AVERAGE value of (nx,ny,nz) on xi-face
c
c
c       2     mx, my, mz             area-weighted normal to eta-face
c             mxhat, myhat, mzhat    unit normal to eta-face
c             lxhat, lyhat, lzhat    ARBITRARY unit normal to (mx,my,mz)
c             nxhat, nyhat, nzhat    ARBITRARY unit normal to (mx,my,mz)
c
c      -2     mx, my, mz             area-weighted normal to eta-face
c             lx, ly, lz             AVERAGE value of (lx,ly,lz) on eta-face
c             nx, ny, nz             AVERAGE value of (nx,ny,nz) on eta-face
c
c     
c       3     nx, ny, nz             area-weighted normal to zeta-face
c             nxhat, nyhat, nzhat    unit normal to zeta-face
c             lxhat, lyhat, lzhat    ARBITRARY unit normal to (nx,ny,nz)
c             mxhat, myhat, mzhat    ARBITRARY unit normal to (nx,ny,nz)
c
c      -3     nx, ny, nz             area-weighted normal to zeta-face
c             lx, ly, lz             AVERAGE value of (lx,ly,lz) on zeta-face
c             mx, my, mz             AVERAGE value of (mx,my,mz) on zeta-face
c     
c
c     DISCUSSION
c
c     arbitrary unit normals
c
c       one the arbitrary unit normal on the xi-face, namely 
c       (mxhat,myhat,mzhat), is obtained from (lxhat,lyhat,lzhat) 
c       by cross product with an arbitrary random vector.  the
c       second arbitrary unit normal on the xi-face, namely
c       (nxhat,nyhat,nzhat), is obtained by cross product of
c       (lxhat,lyhat,lzhat) and (mxhat,myhat,mzhat).  thus, these
c       three unit vectors form an orthogonal set.
c
c       it is important to emphasize, however, that on the xi-face,
c       (mxhat,myhat,mzhat) is in no way related to the unit normal
c       on adjacent eta-faces, and similarly for (nxhat,nyhat,nzhat).
c
c     average unit normals
c
c       on an xi-face, the average value of (mx,my,mz) is
c       defined as the arithmetic average of the four values of
c       the unit normals on the eta-faces which share intersect
c       the xi-face.  thus, the average value of (mx,my,mz)
c       at point P is the average of the values of (mx,my,mz)
c       at points a1, b1, c1 and d1.
c
c           eta
c
c           ^
c           |                         "i" denotes centroid (i,j,k)
c
c           + ---b1-- + ---d1-- +
c           |         |         |
c           |   "i"   P         |    
c           |         |         |
c           + ---a1-- + ---c1-- +     ---> xi
c
c
c       similarly, on an eta-face, the average value of (lx,ly,lz)
c       at point P is the average of the values of (lx,ly,lz) at
c       points a2, b2, c2, d2.
c
c
c           eta

c
c           ^
c           |
c
c           + ------- +
c           |         |
c           c2        d2              "i" denotes centroid (i,j,k)
c           |         |
c           + ---P--- + 
c           |         |        
c           a2  "i"   b2
c           |         |        
c           + ---b--- +      ---> xi
c
c
c
c
c.............  extra option to calculate the mutually orthogonal vectors .........

c      A different method is incorporated into this subroutine to calculate the[u 
c      three mutually orthogonal unit vectors (AIAA-90-0129, Zha & Liu). No random 
c      vector is needed and therefore the unreliability of the random number 
c      generator is removed. A control parameter, kvector, is added.

c      when   kvector  = 0, use random vector method
c      when   kvector  = 1, use present method

c      Zha, 28/07/94

c.................................................................................
c
c     DEFINITION OF VARIABLES
c
c       lxm(i),lym(i),lzm(i)      (lx,ly,lz) at face i-1/2 for cell (i,j,k)
c
c       lxp(i),lyp(i),lzp(i)      (lx,ly,lz) at face i+1/2 for cell (i,j,k)
c
c       mxm(i),mym(i),mzm(i)      (mx,my,mz) at face j-1/2 for cell (i,j,k)
c
c       mxp(i),myp(i),mzp(i)      (mx,my,mz) at face j+1/2 for cell (i,j,k)
c
c       nxm(i),nym(i),nzm(i)      (nx,ny,nz) at face k-1/2 for cell (i,j,k)
c
c       nxp(i),nyp(i),nzp(i)      (nx,ny,nz) at face k+1/2 for cell (i,j,k)
c
c     IMPLICIT STATEMENT
c     
      implicit none
c     
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     imax,                ! maxium surface number
     $     index,               ! face type index
     $     i1, i2,              ! input face poistion index
     $     il, jl, kl,          ! cell number in 3 directions
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     bclower, bcupper

      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x,  y, z

c     imax-ilower = iupper, jupper, kupper
      double precision, dimension(ilower+1:imax-ilower), intent(out)::
     $     lxhat, lyhat, lzhat,
     $     mxhat, myhat, mzhat,
     $     nxhat, nyhat, nzhat,
     $     lx, ly, lz,
     $     mx, my, mz,
     $     nx, ny, nz
c     
c     LOCAL VARIABLES
c     
      double precision, dimension(ilower+1:imax-ilower)::
     $     lxm, lxp, lym,
     >     lyp, lzm, lzp, mxm, mxp, mym, myp, mzm, mzp,
     >     nxm, nxp, nym, nyp, nzm, nzp
c     
      integer::
     $     i, j, k,             ! cell iteration index
     $     ip, jp, kp           ! i+1, j+1, k+1
c     
      double precision::
     $     rprimeinv,
     $     rprime,
     $     dx1, dx2,
     $     dy1, dy2,
     $     dz1, dz2

      integer, parameter::
     $     in_bnd = 7,
     $     per_bnd = 10,
     $     sym_bnd = 8
c     
c *** SUBROUTINE START ***
c     
      SELECT CASE (index)

      CASE (1)                  ! xi surfaces

         j  = i1 
         k  = i2
         jp = j + 1
         kp = k + 1
         
         do i = 1, imax
            dx1 = x(i,jp,kp) - x(i,j,k)
            dy1 = y(i,jp,kp) - y(i,j,k)
            dz1 = z(i,jp,kp) - z(i,j,k)
            dx2 = x(i,j,kp)  - x(i,jp,k)
            dy2 = y(i,j,kp)  - y(i,jp,k)
            dz2 = z(i,j,kp)  - z(i,jp,k)

            lx(i) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            ly(i) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            lz(i) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = lx(i)*lx(i) + ly(i)*ly(i) + lz(i)*lz(i)
            dx2 = dsqrt(dx1)
            dy1 = 1.d0 / dx2

            lxhat(i) = lx(i) * dy1
            lyhat(i) = ly(i) * dy1
            lzhat(i) = lz(i) * dy1
            
            rprime    = dsqrt( lxhat(i)**2 + lyhat(i)**2 )
            if (rprime.ne.0.d0) then
               rprimeinv = 1.d0/rprime
             
               mxhat(i) = - lyhat(i) * rprimeinv
               myhat(i) =   lxhat(i) * rprimeinv
               mzhat(i) =   0.d0
            
               nxhat(i) =   lxhat(i) * lzhat(i) * rprimeinv
               nyhat(i) =   lyhat(i) * lzhat(i) * rprimeinv
               nzhat(i) = - rprime  
            else
               rprime    = dsqrt( lxhat(i)**2 + lzhat(i)**2 )
               rprimeinv = 1.d0/rprime
             
               mxhat(i) = - lzhat(i) * rprimeinv
               myhat(i) =   0.d0
               mzhat(i) =   lxhat(i) * rprimeinv
            
               nxhat(i) =   lxhat(i) * lyhat(i) * rprimeinv
               nyhat(i) = - rprime  
               nzhat(i) =   lyhat(i) * lzhat(i) * rprimeinv
            end if
         end do

      CASE (2)                  ! eta surface

         i  = i1 
         k  = i2
         ip = i + 1
         kp = k + 1

         do j = 1, imax
            dx1 = x(ip,j,kp) - x(i,j,k)
            dy1 = y(ip,j,kp) - y(i,j,k)
            dz1 = z(ip,j,kp) - z(i,j,k)
            dx2 = x(ip,j,k)  - x(i,j,kp)
            dy2 = y(ip,j,k)  - y(i,j,kp)
            dz2 = z(ip,j,k)  - z(i,j,kp)

            mx(j) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            my(j) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            mz(j) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
 
            dx1 = mx(j)*mx(j) + my(j)*my(j) + mz(j)*mz(j)
            dx2 = dsqrt(dx1)
            dy1 = 1.d0 / dx2

            mxhat(j) = mx(j) * dy1
            myhat(j) = my(j) * dy1
            mzhat(j) = mz(j) * dy1
            
            rprime    = dsqrt( myhat(j)**2 + mzhat(j)**2 )
            if (rprime.ne.0.d0) then
               rprimeinv = 1.d0/rprime
            
               nxhat(j) =   0.0
               nyhat(j) = - mzhat(j) * rprimeinv
               nzhat(j) =   myhat(j) * rprimeinv
            
               lxhat(j) = - rprime    
               lyhat(j) =   mxhat(j) * myhat(j) * rprimeinv
               lzhat(j) =   mxhat(j) * mzhat(j) * rprimeinv
            else
               rprime    = dsqrt( myhat(j)**2 + mxhat(j)**2 )
               rprimeinv = 1.d0/rprime
            
               nxhat(j) =   myhat(j) * rprimeinv
               nyhat(j) = - mxhat(j) * rprimeinv
               nzhat(j) =   0.0
            
               lxhat(j) =   mxhat(j) * mzhat(j) * rprimeinv
               lyhat(j) =   mzhat(j) * myhat(j) * rprimeinv
               lzhat(j) = - rprime    
            end if
         end do   
         
      CASE (3)                  ! zeta surface

         i  = i1 
         j  = i2
         ip = i + 1
         jp = j + 1

         do k = 1, imax
            dx1 = x(ip,jp,k) - x(i,j,k)
            dy1 = y(ip,jp,k) - y(i,j,k)
            dz1 = z(ip,jp,k) - z(i,j,k)
            dx2 = x(i,jp,k)  - x(ip,j,k)
            dy2 = y(i,jp,k)  - y(ip,j,k)
            dz2 = z(i,jp,k)  - z(ip,j,k)

            nx(k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            ny(k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            nz(k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = nx(k)*nx(k) + ny(k)*ny(k) + nz(k)*nz(k)
            dx2 = dsqrt(dx1)
            dy1 = 1.d0 / dx2

            nxhat(k) = nx(k) * dy1
            nyhat(k) = ny(k) * dy1
            nzhat(k) = nz(k) * dy1
            
            rprime    = dsqrt( nxhat(k)**2 + nzhat(k)**2 )
            if (rprime.ne.0.d0) then
               rprimeinv = 1.d0/rprime

               lxhat(k) =   nzhat(k) * rprimeinv
               lyhat(k) =   0.d0
               lzhat(k) = - nxhat(k) * rprimeinv
            
               mxhat(k) =   nxhat(k) * nyhat(k) * rprimeinv
               myhat(k) = - rprime
               mzhat(k) =   nzhat(k) * nyhat(k) * rprimeinv
            else
               rprime    = dsqrt( nyhat(k)**2 + nzhat(k)**2 )
               rprimeinv = 1.d0/rprime

               lxhat(k) =   0.d0
               lyhat(k) =   nzhat(k) * rprimeinv
               lzhat(k) = - nyhat(k) * rprimeinv
            
               mxhat(k) = - rprime
               myhat(k) =   nyhat(k) * nxhat(k) * rprimeinv
               mzhat(k) =   nzhat(k) * nxhat(k) * rprimeinv
            end if
         end do   

      CASE (-1)                 ! xi on viscous term

         j  = i1
         k  = i2
         jp = j+1 
         kp = k+1 

c     lx,ly,lz  on xi-face

         do i = ilower+1, iupper
            dx1 = x(i,jp,kp) - x(i,j,k)
            dy1 = y(i,jp,kp) - y(i,j,k)
            dz1 = z(i,jp,kp) - z(i,j,k)
            dx2 = x(i,j,kp)  - x(i,jp,k)
            dy2 = y(i,j,kp)  - y(i,jp,k)
            dz2 = z(i,j,kp)  - z(i,jp,k)

            lx(i) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            ly(i) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            lz(i) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

c     mx,my,mz  on xi-face

         do i = ilower+1,iupper-1
            dx1 = x(i+1,jp,kp) - x(i,jp,k)
            dy1 = y(i+1,jp,kp) - y(i,jp,k)
            dz1 = z(i+1,jp,kp) - z(i,jp,k)
            dx2 = x(i+1,jp,k)  - x(i,jp,kp)
            dy2 = y(i+1,jp,k)  - y(i,jp,kp)
            dz2 = z(i+1,jp,k)  - z(i,jp,kp)

            mxp(i) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            myp(i) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            mzp(i) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = x(i+1,j,kp) - x(i,j,k)
            dy1 = y(i+1,j,kp) - y(i,j,k)
            dz1 = z(i+1,j,kp) - z(i,j,k)
            dx2 = x(i+1,j,k)  - x(i,j,kp)
            dy2 = y(i+1,j,k)  - y(i,j,kp)
            dz2 = z(i+1,j,k)  - z(i,j,kp)

            mxm(i) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            mym(i) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            mzm(i) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

         do i = ilower+2,iupper-1
            mx(i) = 0.25d0*(mxm(i)+mxm(i-1)+mxp(i)+mxp(i-1))
            my(i) = 0.25d0*(mym(i)+mym(i-1)+myp(i)+myp(i-1))
            mz(i) = 0.25d0*(mzm(i)+mzm(i-1)+mzp(i)+mzp(i-1))
         end do

         if (bclower.ne.in_bnd.and.bclower.ne.per_bnd
     $       .and.bclower.ne.sym_bnd) then
           i = 1
           mx(i)   = 0.5d0 * (mxm(i) + mxp(i))
           my(i)   = 0.5d0 * (mym(i) + myp(i))
           mz(i)   = 0.5d0 * (mzm(i) + mzp(i))
           do i = ilower+2,0
             mx(i) = mx(1)
             my(i) = my(1)
             mz(i) = mz(1)
           end do
         end if
         
         if (bcupper.ne.in_bnd.and.bcupper.ne.per_bnd
     $       .and.bcupper.ne.sym_bnd) then
           i = imax
           mx(i) = 0.5d0 * (mxm(i-1) + mxp(i-1))
           my(i) = 0.5d0 * (mym(i-1) + myp(i-1))
           mz(i) = 0.5d0 * (mzm(i-1) + mzp(i-1))
           do i = imax+1,iupper-1
             mx(i) = mx(imax)
             my(i) = my(imax)
             mz(i) = mz(imax)
           end do
	 end if

c     nx,ny,nz  on xi-face

         do i = ilower+1,iupper-1
            dx1 = x(i+1,jp,kp) - x(i,j,kp)
            dy1 = y(i+1,jp,kp) - y(i,j,kp)
            dz1 = z(i+1,jp,kp) - z(i,j,kp)
            dx2 = x(i,jp,kp)  - x(i+1,j,kp)
            dy2 = y(i,jp,kp)  - y(i+1,j,kp)
            dz2 = z(i,jp,kp)  - z(i+1,j,kp)

            nxp(i) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            nyp(i) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            nzp(i) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = x(i+1,jp,k) - x(i,j,k)
            dy1 = y(i+1,jp,k) - y(i,j,k)
            dz1 = z(i+1,jp,k) - z(i,j,k)
            dx2 = x(i,jp,k)  - x(i+1,j,k)
            dy2 = y(i,jp,k)  - y(i+1,j,k)
            dz2 = z(i,jp,k)  - z(i+1,j,k)

            nxm(i) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            nym(i) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            nzm(i) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

         do i = ilower+2, iupper-1
            nx(i) = 0.25d0*(nxm(i)+nxm(i-1)+nxp(i)+nxp(i-1))
            ny(i) = 0.25d0*(nym(i)+nym(i-1)+nyp(i)+nyp(i-1))
            nz(i) = 0.25d0*(nzm(i)+nzm(i-1)+nzp(i)+nzp(i-1))
         end do

         if (bclower.ne.in_bnd.and.bclower.ne.per_bnd
     $       .and.bclower.ne.sym_bnd) then
           i = 1
           nx(i)   = 0.5d0 * (nxm(i) + nxp(i))
           ny(i)   = 0.5d0 * (nym(i) + nyp(i))
           nz(i)   = 0.5d0 * (nzm(i) + nzp(i))
	   do i = ilower+2, 0
             nx(i)   = nx(1)
             ny(i)   = ny(1)
             nz(i)   = nz(1)
	   end do
         end if

         if (bcupper.ne.in_bnd.and.bcupper.ne.per_bnd
     $       .and.bcupper.ne.sym_bnd) then
           i = imax
           nx(i) = 0.5d0 * (nxm(i-1) + nxp(i-1))
           ny(i) = 0.5d0 * (nym(i-1) + nyp(i-1))
           nz(i) = 0.5d0 * (nzm(i-1) + nzp(i-1))
           do i = imax+1, iupper-1
             nx(i) = nx(imax)
             ny(i) = ny(imax)
             nz(i) = nz(imax)
	   end do
         end if

      case (-2)                 ! eta surface on viscous terms

         i  = i1 
         k  = i2
         ip = i + 1
         kp = k + 1

c     mx,my,mz  on eta-face

         do  j = jlower+1, jupper
            dx1 = x(ip,j,kp) - x(i,j,k)
            dy1 = y(ip,j,kp) - y(i,j,k)
            dz1 = z(ip,j,kp) - z(i,j,k)
            dx2 = x(ip,j,k)  - x(i,j,kp)
            dy2 = y(ip,j,k)  - y(i,j,kp)
            dz2 = z(ip,j,k)  - z(i,j,kp)
            
            mx(j) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            my(j) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            mz(j) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

c     lx,ly,lz  on eta-face

         do j = jlower+1,jupper-1
            dx1 = x(ip,j+1,kp) - x(ip,j,k)
            dy1 = y(ip,j+1,kp) - y(ip,j,k)
            dz1 = z(ip,j+1,kp) - z(ip,j,k)
            dx2 = x(ip,j,kp)  - x(ip,j+1,k)
            dy2 = y(ip,j,kp)  - y(ip,j+1,k)
            dz2 = z(ip,j,kp)  - z(ip,j+1,k)

            lxp(j) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            lyp(j) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            lzp(j) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = x(i,j+1,kp) - x(i,j,k)
            dy1 = y(i,j+1,kp) - y(i,j,k)
            dz1 = z(i,j+1,kp) - z(i,j,k)
            dx2 = x(i,j,kp)  - x(i,j+1,k)
            dy2 = y(i,j,kp)  - y(i,j+1,k)
            dz2 = z(i,j,kp)  - z(i,j+1,k)

            lxm(j) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            lym(j) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            lzm(j) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do
c     
         do j = jlower+2,jupper-1
           lx(j) = 0.25d0*(lxm(j)+lxm(j-1)+lxp(j)+lxp(j-1))
           ly(j) = 0.25d0*(lym(j)+lym(j-1)+lyp(j)+lyp(j-1))
           lz(j) = 0.25d0*(lzm(j)+lzm(j-1)+lzp(j)+lzp(j-1))
	 end do

         if (bclower.ne.in_bnd.and.bclower.ne.per_bnd
     $       .and.bclower.ne.sym_bnd) then
           j = 1
           lx(j)   = 0.5d0 * (lxm(j) + lxp(j))
           ly(j)   = 0.5d0 * (lym(j) + lyp(j))
           lz(j)   = 0.5d0 * (lzm(j) + lzp(j))
	   do j = jlower+2,0
             lx(j)   = lx(1)
             ly(j)   = ly(1)
             lz(j)   = lz(1)
           end do
         end if

         if (bcupper.ne.in_bnd.and.bcupper.ne.per_bnd
     $       .and.bcupper.ne.sym_bnd) then
           j = imax
           lx(j) = 0.5d0 * (lxm(j-1) + lxp(j-1))
           ly(j) = 0.5d0 * (lym(j-1) + lyp(j-1))
           lz(j) = 0.5d0 * (lzm(j-1) + lzp(j-1))
	   do j = imax+1, jupper-1
             lx(j)   = lx(imax)
             ly(j)   = ly(imax)
             lz(j)   = lz(imax)
           end do
         end if

c     nx,ny,nz  on eta-face

         do j = jlower+1,jupper-1
            dx1 = x(ip,j+1,kp) - x(i,j,kp)
            dy1 = y(ip,j+1,kp) - y(i,j,kp)
            dz1 = z(ip,j+1,kp) - z(i,j,kp)
            dx2 = x(i,j+1,kp)  - x(ip,j,kp)
            dy2 = y(i,j+1,kp)  - y(ip,j,kp)
            dz2 = z(i,j+1,kp)  - z(ip,j,kp)

            nxp(j) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            nyp(j) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            nzp(j) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = x(ip,j+1,k) - x(i,j,k)
            dy1 = y(ip,j+1,k) - y(i,j,k)
            dz1 = z(ip,j+1,k) - z(i,j,k)
            dx2 = x(i,j+1,k)  - x(ip,j,k)
            dy2 = y(i,j+1,k)  - y(ip,j,k)
            dz2 = z(i,j+1,k)  - z(ip,j,k)

            nxm(j) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            nym(j) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            nzm(j) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

         do j = jlower+2, jupper-1
            nx(j) = 0.25d0*(nxm(j)+nxm(j-1)+nxp(j)+nxp(j-1))
            ny(j) = 0.25d0*(nym(j)+nym(j-1)+nyp(j)+nyp(j-1))
            nz(j) = 0.25d0*(nzm(j)+nzm(j-1)+nzp(j)+nzp(j-1))
         end do

         if (bclower.ne.in_bnd.and.bclower.ne.per_bnd
     $       .and.bclower.ne.sym_bnd) then
           j = 1
           nx(j)   = 0.5d0 * (nxm(j) + nxp(j))
           ny(j)   = 0.5d0 * (nym(j) + nyp(j))
           nz(j)   = 0.5d0 * (nzm(j) + nzp(j))
	   do j = jlower+2, 0
             nx(j)   = nx(1)
             ny(j)   = ny(1)
             nz(j)   = nz(1)
           end do
         end if

         if(bcupper.ne.in_bnd.and.bcupper.ne.per_bnd
     $      .and.bcupper.ne.sym_bnd) then
           j = imax
           nx(j) = 0.5d0 * (nxm(j-1) + nxp(j-1))
           ny(j) = 0.5d0 * (nym(j-1) + nyp(j-1))
           nz(j) = 0.5d0 * (nzm(j-1) + nzp(j-1))
	   do j = imax+1, jupper-1
             nx(j)   = nx(imax)
             ny(j)   = ny(imax)
             nz(j)   = nz(imax)
           end do
         end if

      CASE (-3)                 ! zeta surface on viscous terms

         i  = i1 
         j  = i2
         ip = i + 1
         jp = j + 1

c     nx,ny,nz  on zeta-face

         do k = klower+1, kupper
            dx1 = x(ip,jp,k) - x(i,j,k)
            dy1 = y(ip,jp,k) - y(i,j,k)
            dz1 = z(ip,jp,k) - z(i,j,k)
            dx2 = x(i,jp,k)  - x(ip,j,k)
            dy2 = y(i,jp,k)  - y(ip,j,k)
            dz2 = z(i,jp,k)  - z(ip,j,k)

            nx(k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            ny(k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            nz(k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

c     lx,ly,lz  on zeta-face

         do k = klower+1,kupper-1
            dx1 = x(ip,jp,k+1) - x(ip,j,k)
            dy1 = y(ip,jp,k+1) - y(ip,j,k)
            dz1 = z(ip,jp,k+1) - z(ip,j,k)
            dx2 = x(ip,j,k+1)  - x(ip,jp,k)
            dy2 = y(ip,j,k+1)  - y(ip,jp,k)
            dz2 = z(ip,j,k+1)  - z(ip,jp,k)

            lxp(k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            lyp(k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            lzp(k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = x(i,jp,k+1) - x(i,j,k)
            dy1 = y(i,jp,k+1) - y(i,j,k)
            dz1 = z(i,jp,k+1) - z(i,j,k)
            dx2 = x(i,j,k+1)  - x(i,jp,k)
            dy2 = y(i,j,k+1)  - y(i,jp,k)
            dz2 = z(i,j,k+1)  - z(i,jp,k)

            lxm(k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            lym(k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            lzm(k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

         do k = klower+2,kupper-1
           lx(k) = 0.25d0*(lxm(k)+lxm(k-1)+lxp(k)+lxp(k-1))
           ly(k) = 0.25d0*(lym(k)+lym(k-1)+lyp(k)+lyp(k-1))
           lz(k) = 0.25d0*(lzm(k)+lzm(k-1)+lzp(k)+lzp(k-1))
	 end do

         if (bclower.ne.in_bnd.and.bclower.ne.per_bnd
     $       .and.bclower.ne.sym_bnd) then
           k = 1
           lx(k) = 0.5d0 * (lxm(k) + lxp(k))
           ly(k) = 0.5d0 * (lym(k) + lyp(k))
           lz(k) = 0.5d0 * (lzm(k) + lzp(k))
	   do k = klower+2, 0
             lx(k) = lx(1)
             ly(k) = ly(1)
             lz(k) = lz(1)
           end do
         end if

         if (bcupper.ne.in_bnd.and.bcupper.ne.per_bnd
     $       .and.bcupper.ne.sym_bnd) then
           k = imax
           lx(k) = 0.5d0 * (lxm(k-1) + lxp(k-1))
           ly(k) = 0.5d0 * (lym(k-1) + lyp(k-1))
           lz(k) = 0.5d0 * (lzm(k-1) + lzp(k-1))
	   do k = imax+1,kupper-1
             lx(k) = lx(imax)
             ly(k) = ly(imax)
             lz(k) = lz(imax)
           end do
         end if
                  
c     mx,my,mz  on zeta-face
         
         do k = klower+1,kupper-1
            dx1 = x(ip,jp,k+1) - x(i,jp,k)
            dy1 = y(ip,jp,k+1) - y(i,jp,k)
            dz1 = z(ip,jp,k+1) - z(i,jp,k)
            dx2 = x(ip,jp,k)  - x(i,jp,k+1)
            dy2 = y(ip,jp,k)  - y(i,jp,k+1)
            dz2 = z(ip,jp,k)  - z(i,jp,k+1)

            mxp(k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            myp(k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            mzp(k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)

            dx1 = x(ip,j,k+1) - x(i,j,k)
            dy1 = y(ip,j,k+1) - y(i,j,k)
            dz1 = z(ip,j,k+1) - z(i,j,k)
            dx2 = x(ip,j,k)  - x(i,j,k+1)
            dy2 = y(ip,j,k)  - y(i,j,k+1)
            dz2 = z(ip,j,k)  - z(i,j,k+1)

            mxm(k) = 0.5d0 * (dy1 * dz2 - dz1 * dy2)
            mym(k) = 0.5d0 * (dz1 * dx2 - dx1 * dz2)
            mzm(k) = 0.5d0 * (dx1 * dy2 - dy1 * dx2)
         end do

         do k = klower+2,kupper-1
           mx(k) = 0.25d0*(mxm(k)+mxm(k-1)+mxp(k)+mxp(k-1))
           my(k) = 0.25d0*(mym(k)+mym(k-1)+myp(k)+myp(k-1))
           mz(k) = 0.25d0*(mzm(k)+mzm(k-1)+mzp(k)+mzp(k-1))
	 end do

         if (bclower.ne.in_bnd.and.bclower.ne.per_bnd
     $       .and.bclower.ne.sym_bnd) then
           k = 1
           mx(k) = 0.5d0 * (mxm(k) + mxp(k))
           my(k) = 0.5d0 * (mym(k) + myp(k))
           mz(k) = 0.5d0 * (mzm(k) + mzp(k))
	   do k = klower+2, 0
             mx(k) = mx(1)
             my(k) = my(1)
             mz(k) = mz(1)
           end do
         end if

         if (bcupper.ne.in_bnd.and.bcupper.ne.per_bnd
     $       .and.bcupper.ne.sym_bnd) then
           k = imax
           mx(k) = 0.5d0 * (mxm(k-1) + mxp(k-1))
           my(k) = 0.5d0 * (mym(k-1) + myp(k-1))
           mz(k) = 0.5d0 * (mzm(k-1) + mzp(k-1))
	   do k = imax+1, kupper
             mx(k)   = mx(imax)
             my(k)   = my(imax)
             mz(k)   = mz(imax)
           end do
         end if

      END SELECT

      return
      end
