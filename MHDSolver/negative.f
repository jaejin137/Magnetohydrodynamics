c======================================================================c
c
      subroutine negative(il, jl, kl, nl, q, gamma,
     $     ilower, iupper, jlower, jupper, klower, kupper, ke,
     $     qt, ronum)
c
c     checks for negative density, pressure and temperature
c
c     IMPLICIT STATEMENT
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     il, jl, kl,
     $     nl, 
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper
c
      double precision, intent(inout)::
     $     q(ilower:iupper,jlower:jupper,klower:kupper,nl)
c
      double precision, intent(in),
     $     dimension(ilower:iupper,jlower:jupper,klower:kupper,3)::
     $     qt
c
      double precision, intent(in)::
     $     gamma, ronum
c
      logical, intent(in)::
     $     ke
c
c     LOCAL VARIABLES
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     nr, np, n            ! number for cells with
     $                          ! rho and pressure negative value
c
      double precision::
     $     gamma1, psi, qq, theta,
     $     p                    ! pressure
c
c --- SUBROUTINE START ---
c      
      gamma1 = gamma-1.d0
c
      nr = 0
      np = 0
      do k = 1,kl
        do j = 1,jl
          do i = 1,il
            if ( q(i,j,k,1).lt.0.d0 )  then
              nr = nr + 1
              print *, 'rho<0: ', i, j, k, q(i,j,k,1)
            end if
            if ( dabs(ronum).gt.1e-9 ) then
              psi = qt(i,j,k,2)*qt(i,j,k,2)+qt(i,j,k,3)*qt(i,j,k,3)
              qq  = 0.5d0*( q(i,j,k,2)*q(i,j,k,2)+q(i,j,k,3)*q(i,j,k,3)+
     >              q(i,j,k,4)*q(i,j,k,4) )/q(i,j,k,1)/q(i,j,k,1)
              theta = qq-0.5d0*psi
              p = gamma1*( q(i,j,k,5)-q(i,j,k,1)*theta )
            else
              p = gamma1*( q(i,j,k,5)-0.5d0*(q(i,j,k,2)**2 +
     >            q(i,j,k,3)**2 + q(i,j,k,4)**2)/q(i,j,k,1) )
            end if
            if ( p.lt.0.d0 ) then
              np = np + 1
              print *, 'p<0: ', i, j, k, p
              write(*,*) ( q(i,j,k,n),n=1,5 )
            end if
          end do
        end do
      end do
c
      if(np.gt.0) print *, np, 'negative pressure points found'
      if(nr.gt.0) print *, nr, 'negative density points found'
c           
      if((np.gt.0).or.(nr.gt.0)) stop
c
      return
      end 
c
c======================================================================c
