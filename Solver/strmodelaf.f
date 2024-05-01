C======================================================================C
clist()
      subroutine strmodelaf(dtow, tintvl, cd, cl, ct, xsl, resdl,
     $                      uinf, walf)
c
c     This subroutine is used to solve the governing equations for
c     two-degree-of-freedom model.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      double precision, intent(in)::
     $     cd, cl, ct,          ! drag, lift & torque coefficients
     $     tintvl,              ! real time interval
     $     dtow,                ! iteration time interval
     $     uinf,                ! free stream velocity
     $     walf                 ! uncoupled natural frequency of typical
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
     $     xalf,                ! static unbalance
     $     ra2,                 ! squared radius of gyration
     $     ust,                 ! reduced velocity
     $     mus,                 ! airfoil to fluid mass ratio
     $     wh,                  ! uncoupled natural frequency of typical
     $                          ! section in torsion & bending, respectively
     $     cct, ccl,            ! coefficients appear in gvrn eqn
     $     dsl(4)               ! increment of result
c
      double precision::
     $     am(4,4),             ! k matrix
     $     mm(4,4),             ! m matrix
     $     bm(4),               ! b matrix
     $     tvls,                ! time interval for structural eqns
     $     lref,                ! reference length
     $     diag
c
      integer::
     $     i, j                 ! local variables
c
c *** SUBROUTINE START ***
c
c ... coefficients in governing equtions
c
      xalf = 1.8d0
      ra2  = 3.48d0
      lref = 0.5d0
      mus  = 60.d0
      wh = walf
c
      ust = 2.d0*uinf/walf/lref
      ccl = 0.318309892d0*ust*ust*cl/mus
      cct = 0.636619783d0*ust*ust*ct/mus
c
      tvls = walf*lref*tintvl/uinf
c
c ... k matrix
c
      am = 0.d0
      am(1,2) = -1.d0
      am(2,1) =  (wh/walf)**2
      am(3,4) = -1.d0
      am(4,3) =  ra2
c
c ... m matrix
c
      mm = 0.d0
      mm(1,1) = 1.d0
      mm(2,2) = 1.d0
      mm(2,4) = xalf
      mm(3,3) = 1.d0
      mm(4,2) = xalf
      mm(4,4) = ra2
c
c ... q matrix
c
      bm(1) = 0.d0
      bm(2) = ccl
      bm(3) = 0.d0
      bm(4) = cct
c
c ... b_hat matrix
c
      bm(1) = bm(1) - am(1,2)*xsl(2,1)
      bm(2) = bm(2) - am(2,1)*xsl(1,1)
      bm(3) = bm(3) - am(3,4)*xsl(4,1)
      bm(4) = bm(4) - am(4,3)*xsl(3,1)
c
      diag = 0.5d0/tvls
      bm(1) = bm(1) - 
     >     diag*(3.d0*xsl(1,1)-4.d0*xsl(1,2)+xsl(1,3))*mm(1,1)
c
      bm(2) = bm(2) -
     >     diag*(3.d0*xsl(2,1)-4.d0*xsl(2,2)+xsl(2,3))*mm(2,2)-
     >     diag*(3.d0*xsl(4,1)-4.d0*xsl(4,2)+xsl(4,3))*mm(2,4)
c
      bm(3) = bm(3) -
     >     diag*(3.d0*xsl(3,1)-4.d0*xsl(3,2)+xsl(3,3))*mm(3,3)
c
      bm(4) = bm(4) - 
     >     diag*(3.d0*xsl(2,1)-4.d0*xsl(2,2)+xsl(2,3))*mm(4,2)-
     >     diag*(3.d0*xsl(4,1)-4.d0*xsl(4,2)+xsl(4,3))*mm(4,4)
c
      bm(1) = bm(1)*dtow
      bm(2) = bm(2)*dtow
      bm(3) = bm(3)*dtow
      bm(4) = bm(4)*dtow
c
c ... a_hat matrix
c
      diag = 1.5d0*dtow/tvls
      do 10 j = 1, 4
      do 10 i = 1, 4
        am(i,j) = dtow * am(i,j) + diag * mm(i,j)
   10 continue
c
      diag = 1.d0
      do 30 i = 1, 4
        am(i,i) = diag + am(i,i)
   30 continue
c
c ... solving the system by gauss_seidel iteration
c
      call gauss_seidel(am, bm, dsl)
c
c ... updates the result vecter
c
      do 60 i = 1, 4
        xsl(i,1) = xsl(i,1) + dsl(i)
   60 continue
c
c ... calculates the residual, L2 norm
c
      resdl = dsqrt( dsl(1)*dsl(1) + dsl(2)*dsl(2) +
     >               dsl(3)*dsl(3) + dsl(4)*dsl(4) )
c
      return
      end
c
C======================================================================C
