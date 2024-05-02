      subroutine strafnlr(dtow, tintvl, cl, ct, mus,
     $                    walf, wh, xsl, resdl, uinf, alf0)
c
c ********************************************************************
c     This subroutine is used to solve the governing equations for
c     two-degree-of-freedom model: M*(dS/dt)+K*S = Q.
c ********************************************************************
c
      implicit none
c
c     INTERFACE VARIABLES
c
      double precision, intent(in)::
     $     cl, ct,              ! lift & torque coefficients
     $     tintvl,              ! real time interval
     $     dtow,                ! iteration time interval
     $     walf,                ! undamped natural torsional freq.
     $     wh,                  ! undamped natural bending freq.
     $     uinf,                ! free stream velocity
     $     alf0,                ! Initial angle of attack without wind
     $     mus                  ! mass ratio
c
      double precision, intent(inout)::
     $     xsl(4,3)            ! structual variables at future, current,
                               ! previous
c
      double precision, intent(out)::
     $     resdl               ! residual, L2 norm
c
c     LOCAL VARIABLES
c
      double precision::
     $     xalf,                ! distance of elastic axis to center of mass
     $     ralf,                ! radius of gyration about the elastic axis
     $     ra2,                 ! square of ralf
     $     dltalf,              ! Lehr's pitch-damping coefficient
     $     dlth,                ! Lehr's plunge-damping coefficient
     $     ca,                  ! coefficients appear in gvrn eqn
     $     dsl(4),              ! increment of result
     $     tvls,                ! time interval for structural eqns
     $     lref                 ! reference length
c
      double precision::
     $     am(4,4),             ! M matrix
     $     ak(4,4),             ! K matrix
     $     bm(4)                ! Q vector
c
      integer::
     $     i, j                 ! local variables
c
c *** SUBROUTINE START
c
      xalf = 0.0484d0
      ralf = 0.197d0
      dltalf = 0.0041d0
      dlth = 0.0073d0
c
      lref = 0.3d0
      tvls = tintvl
c
      do i=1,3
        xsl(3,i) = xsl(3,i) - alf0
      end do
c
      ra2 = ralf**2
      ca = 0.63661977d0/mus
c
c ... M matrix
c
      am = 0.d0
      am(1,1) = 1.d0
      am(2,2) = 1.d0
      am(2,4) = xalf
      am(3,3) = 1.d0
      am(4,2) = am(2,4)
      am(4,4) = ra2
c
c ... K matrix
c
      ak = 0.d0
      ak(1,2) = -1.d0
      ak(2,1) = wh**2
      ak(2,2) = 2.d0*dlth*wh
      ak(3,4) = -1.d0
      ak(4,3) = ra2*walf**2
      ak(4,4) = 2.d0*ra2*dltalf*walf
c
c ... Q vector
c
      bm(1) = 0.d0
      bm(2) = -ca*cl
      bm(3) = 0.d0
      bm(4) = ca*ct
c
c ... b_hat matrix
c
      bm(1) = bm(1)                    - ak(1,2)*xsl(2,1)
      bm(2) = bm(2) - ak(2,1)*xsl(1,1) - ak(2,2)*xsl(2,1)
      bm(3) = bm(3)                    - ak(3,4)*xsl(4,1)
      bm(4) = bm(4) - ak(4,3)*xsl(3,1) - ak(4,4)*xsl(4,1)
c
      ca = 0.5d0/tvls
      do i=1,4
        dsl(i) = ca*(3.d0*xsl(i,1)-4.d0*xsl(i,2)+xsl(i,3))
      end do
c
      bm(1) = bm(1) - dsl(1)
      bm(2) = bm(2) - am(2,2)*dsl(2) - am(2,4)*dsl(4)
      bm(3) = bm(3) - dsl(3)
      bm(4) = bm(4) - am(4,2)*dsl(2) - am(4,4)*dsl(4)
c
      bm(1) = bm(1)*dtow
      bm(2) = bm(2)*dtow
      bm(3) = bm(3)*dtow
      bm(4) = bm(4)*dtow
c
c ... a_hat matrix
c
      ca = 1.5d0/tvls
      do j = 1, 4
        do i = 1, 4
          am(i,j) = dtow*(ca*am(i,j)+ak(i,j))
        end do
      end do
c
      do i = 1, 4
        am(i,i) = 1.d0 + am(i,i)
      end do
c     
c ... solving the system by gauss_seidel iteration
c     
      call gauss_seidel(am, bm, dsl)
c
c ... updates the result vecter
c
      do i = 1, 4
        xsl(i,1) = xsl(i,1) + dsl(i)
      end do
c
      do i=1,3
        xsl(3,i) = xsl(3,i) + alf0
      end do
c
c ... calculates the residual, L2 norm
c
      resdl = dsqrt( dsl(1)*dsl(1) + dsl(2)*dsl(2) +
     $               dsl(3)*dsl(3) + dsl(4)*dsl(4) )
c
      return
      end
