

      subroutine diverg(x,y,z,vol,ilower,iupper,jlower,jupper,
     $     klower,kupper,il,jl,kl,nl,q,diver) 
c
c     compute turbulence viscosity
c
      implicit none

      integer,intent(in)::
     $     il,jl,kl,nl,
     $     ilower,iupper,jlower,jupper,klower,kupper
c
c     INTERFACE VARIABLES
c
      double precision,intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x,y,z
      double precision,intent(in)::
     $     vol(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)

      double precision,intent(out)::
     $   diver(ilower:iupper,jlower:jupper,klower:kupper)

      double precision,intent(in)::
     $     q(ilower:iupper,jlower:jupper,klower:kupper,nl)
c
c     LOCAL VARIABLES
c
      double precision,
     $     dimension(il+1,jl+1,kl+1)::
     $     lx,ly,lz,mx,my,mz,nx,ny,nz
      double precision,dimension(ilower:iupper,jlower:jupper,
     $     klower:kupper)::
     $     u,v,w

      double precision::
     $     dudxi,dudeta,dudzeta,dvdxi,dvdeta,dvdzeta,dwdxi,dwdeta,
     $     dwdzeta

      integer::
     $     i,j,k

c
c *** START SUBROUTINE
c
!compute metric
      call metric_all(x, y, z, lx, ly, lz, mx, my, mz, nx, ny, nz,
     $     il, jl, kl, ilower, iupper, jlower, jupper, klower, kupper)

! compute variable values
	diver=0.0d0
         do  k = klower,kupper
            do  j = jlower,jupper
               do  i = ilower,iupper
                   u(i,j,k)=q(i,j,k,2)
                   v(i,j,k)=q(i,j,k,3)
                   w(i,j,k)=q(i,j,k,4)
               end do
            end do
         end do

      do i = 1, il
         do j = 1, jl
           do k = 1, kl
              dudxi = u(i+1,j,k)-u(i-1,j,k)
              dvdxi = v(i+1,j,k)-v(i-1,j,k)
              dwdxi = w(i+1,j,k)-w(i-1,j,k)

              dudeta = u(i,j+1,k)-u(i,j-1,k)
              dvdeta = v(i,j+1,k)-v(i,j-1,k)
              dwdeta = w(i,j+1,k)-w(i,j-1,k)

              dudzeta = u(i,j,k+1)-u(i,j,k-1)
              dvdzeta = v(i,j,k+1)-v(i,j,k-1)
              dwdzeta = w(i,j,k+1)-w(i,j,k-1)
                     
              diver(i,j,k)=dudxi*lx(i,j,k)+dudeta*mx(i,j,k)
     $             +dudzeta*nx(i,j,k)+dvdxi*ly(i,j,k)
     $             +dvdeta*my(i,j,k)+dvdzeta*ny(i,j,k)
     $             +dwdxi*lz(i,j,k)+dwdeta*mz(i,j,k)
     $             +dwdzeta*nz(i,j,k)
            end do
         end do
      end do


      return
      end 

