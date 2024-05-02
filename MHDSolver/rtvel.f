c======================================================================c
clist()
      subroutine rtvel(il, jl, kl, x, y, z, ilower, iupper,
     $                 jlower, jupper, klower, kupper, qt, ronum)
c
c     This subroutine is used to calculate the grid rotating velocity
c     around the x-axis.
c
      implicit none
c
c     qt(i,1) = 0.d0, qt(i,2) = -ronum*zctr, qt(i,3) = ronum*yctr
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl          ! cell number in 3 directions
c
      double precision, intent(in):: 
     $     ronum               ! 1/Rossby number
c
      double precision, intent(out)::
     $     qt(ilower:iupper,jlower:jupper,klower:kupper,3)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
c     LOCAL VARIABLES
c
      double precision::
     $     yctr, zctr
c
      integer::
     $     n,
     $     i, j, k              ! cell iteration index
c
c --- SUBROUTINE START ---
c
c *** calculate the relative velocity at each node point ***
c
      do 10 k = 1, kl
      do 10 j = 1, jl
      do 10 i = 1, il
        call yzctr(il, jl, kl, x, y, z, ilower, iupper,
     >             jlower, jupper, klower, kupper,
     >             i, j, k, yctr, zctr)
c
        qt(i,j,k,1) =  0.d0
        qt(i,j,k,2) = -ronum*zctr
        qt(i,j,k,3) =  ronum*yctr
   10 continue
c
c *** boundary nodes ***
c
      do 70 n = 1,  3
      do 70 k = 1, kl
      do 70 j = 1, jl
        do 30 i = ilower, 0
          qt(i,j,k,n) = qt(1,j,k,n)
   30   continue
        do 50 i = il+1, iupper
          qt(i,j,k,n) = qt(il,j,k,n)
   50   continue
   70 continue
c
      do 130 n = 1,  3
      do 130 k = 1, kl
      do 130 i = 1, il
        do 90 j = jlower, 0
          qt(i,j,k,n) = qt(i,1,k,n)
   90   continue
        do 110 j = jl+1, jupper
          qt(i,j,k,n) = qt(i,jl,k,n)
  110   continue
  130 continue
c
      do 190 n = 1,  3
      do 190 j = 1, jl
      do 190 i = 1, il
        do 150 k = klower, 0
          qt(i,j,k,n) = qt(i,j,1,n)
  150   continue
        do 170 k = kl+1, kupper
          qt(i,j,k,n) = qt(i,j,kl,n)
  170   continue
  190 continue
c
      return
      end
c
c======================================================================c
