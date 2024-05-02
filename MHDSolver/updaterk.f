      subroutine updaterk(rank, irk, il, jl, kl, nl, q, q_n, rhs,
     $     ilower, iupper, jlower, jupper, klower, kupper)
c
c     update q using RK method
c
c     IMPLICIT STATEMENT
c
      implicit none
c
c     include header file
c
      include "/usr/local/include/mpif.h"
c
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     il, jl, kl,
     $     nl,
     $     irk,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     rank

      double precision, dimension(il,jl,kl,nl),intent(in)::
     $     rhs

      double precision, intent(inout), 
     $     dimension(ilower:iupper,jlower:jupper,klower:kupper,nl)::
     $     q, q_n
c
c     LOCAL VARIABLES
c
      integer i, j, k, n        ! iteration index
      
      double precision rk_coeft
c
c *** SUBROUTINE START ***
c
      if(irk.le.2) rk_coeft = 0.5d0
      if(irk.gt.2) rk_coeft = 1.0d0

      do k = 1,kl
        do j = 1,jl
          do i = 1,il
            do n = 1,nl
              q(i,j,k,n)=q_n(i,j,k,n)+rk_coeft*rhs(i,j,k,n)
            end do
          end do
        end do
      end do
      return
      end
