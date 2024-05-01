      subroutine gauss_seidel(am, bm, xm)
c
c     This subroutine is used to solve linear system by gauss seidel
c     iteration.
c
      implicit none
c
c *** INTERFACE VARIABLES
c
      double precision, intent(in)::
     $     am(4,4),           ! a matrix of the system
     $     bm(4)              ! b matrix of the system
c
      double precision, intent(out)::
     $     xm(4)              ! solution vecter
c
c *** LOCAL VARIABLES
c
      integer::
     $     nn,                ! size of the linear system
     $     max,               ! max number of iteration
     $     iter,              ! number of iteration
     $     i, j               ! local variables ndex
c
      double precision::
     $     error,             ! criterion of accuracy
     $     dif, difmx,        ! error for each variable & max error
     $     temp, sum          ! local variables
c
c *** SUBROUTINE START ***
c
c ... set up some parameters
c
      nn = 4
      max = 1000
      iter = 0
      error = 1.e-12
      xm = 0.d0
c
c ... iteration starts
c
   10 difmx = 0.d0
      iter = iter + 1
      if ( iter.gt.max ) stop 1
c
      do 90 i = 1, nn
c
        temp = xm(i)
        sum = 0.d0
c
        do 30 j = 1, nn
          if ( j.eq.i ) go to 30
          sum = sum + xm(j)*am(i,j)
   30   continue
c
        xm(i) = ( bm(i) - sum )/am(i,i)
        dif = dabs( temp - xm(i) )
        if ( dif.gt.difmx ) difmx = dif
c
   90 continue
c
      if ( difmx.gt.error ) go to 10

c     write(*,*)'gauss-seidel ',i,difmx
c
      return
      end
