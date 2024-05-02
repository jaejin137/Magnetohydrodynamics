C======================================================================C
clist(6)
      subroutine urs_nri(dt, x, y, z, q, q_n, q_bi, gamma, iter_gs,
     $     rhsbi, bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $     bc_zeta_lower, bc_zeta_upper, il, jl, kl, nl,
     $     dim1, dim2, ilower, iupper, jlower, jupper, klower, kupper,
     $     vol, dual_t, tintvl, dual_time, pinlet, poutlet, machinf,
     $     ptotal, ttotal, mx2, sigma, cst_nri, ced_nri)
c
c     use Upwind Relaxation Sweeping(urs) with Gauss-Seidel
c     to invert the implicit LHS
c
c     IMPLICIT STATEMENT
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     il, jl, kl, nl, dim2, dim1, iter_gs,
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     dual_t
c
      double precision, intent(in)::
     $     rhsbi(jl,kl,nl,2), ptotal, machinf, ttotal, mx2, sigma(2)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol
c
      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q
c
      double precision, intent(in)::
     $     q_n(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in)::
     $     gamma, dt(il,jl,kl), tintvl, dual_time, pinlet,
     $     poutlet
c
      integer, intent(in)::
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl)
c
      double precision, intent(inout)::
     $     cst_nri(nl,nl,jl,kl,5), ced_nri(nl,nl,jl,kl,5)
c
      double precision, intent(out),
     $     dimension(jlower:jupper, klower:kupper, nl, 2)::
     $     q_bi
c
c     LOCAL VARIABLES
c
      double precision, dimension(nl,nl,dim1)::
     $     capa, capb, capc, f(nl,dim1), soln(nl,dim1)
c
      integer:: ilp, jlp, klp
c
      double precision:: am(nl,nl)
c
      integer:: n, ii, jj, i, j, k, nblocks,
     $          extra
c
c ... the array q_bi is used for dQ in this subroutine
c
c *** SUBROUTINE START ***
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c
c *** outflow boundary ( left-hand-side )
c
      if ( bc_xi_upper(jl, kl).eq.15 ) then
c
        i = il
        do 100 j = 1,jl
        do 100 k = 1,kl
          call sl_nri(j, k, il, jl, kl, nl, x, y, z,
     $              ilower, iupper, jlower, jupper, klower, kupper,
     $              gamma, q, am, poutlet, mx2, sigma)
c
          do 80 jj = 1,nl
          do 80 ii = 1,nl
            ced_nri(ii,jj,j,k,1) = ced_nri(ii,jj,j,k,1) +
     >                           dt(i,j,k) * am(ii,jj) / vol(i,j,k)
   80     continue
c
          do 90 ii = 1,nl
            ced_nri(ii,ii,j,k,1) = ced_nri(ii,ii,j,k,1) + 1.d0
c
            if ( dual_t.eq.1 ) then
              if ( dual_time.lt.1e-10 ) then
                ced_nri(ii,ii,j,k,1) = ced_nri(ii,ii,j,k,1)
     >                             + dt(i,j,k)/tintvl
              else
                ced_nri(ii,ii,j,k,1) = ced_nri(ii,ii,j,k,1)
     >                             + 1.5d0*dt(i,j,k)/tintvl
              end if
            end if
   90     continue
  100   continue
c
c ... start Gauss-Seidel Relaxation Sweeping
c ... eta-constant line
c
        do 200 k = 1, kl
c       
          do 110 j  = 1,jl
          do 110 ii = 1,nl
          do 110 jj = 1,nl
            capb(ii,jj,j) = ced_nri(ii,jj,j,k,2) 
            capc(ii,jj,j) = ced_nri(ii,jj,j,k,3) 
            capa(ii,jj,j) = ced_nri(ii,jj,j,k,1)
  110     continue
c
c ... form right-hand-side vector f(n,j)
c
c ... initilize f(n,j)
c
          do 170 j = 1,jl
          do 170 n = 1,nl
            f(n,j) = rhsbi(j,k,n,2) 
  170     continue
c
c ... solve linear system
c
          nblocks = jl
          call block(nblocks, nl, dim1, capa, capb, capc,
     >               soln, f)
c         
          do 190 j = 1,jl
          do 190 n = 1,nl
            q_bi(j,k,n,2) = soln(n,j)
  190     continue
c
  200   continue
c
      end if
c
c *** inflow boundary ( left-hand-side )
c
      if ( bc_xi_lower(1,1).eq.16.or.bc_xi_lower(1,1).eq.17 ) then
c
        extra = 0 
        if ( bc_xi_lower(1,1).eq.17 ) extra = 1
c
        i = 1
        do 230 j = 1,jl
        do 230 k = 1,kl
          call isl_nri(j, k, il, jl, kl, nl, x, y, z,
     $            ilower, iupper, jlower, jupper, klower, kupper,
     $            gamma, q, am, pinlet, machinf, ptotal, ttotal,
     $            mx2, sigma, extra)
c
          do 210 jj = 1,nl
          do 210 ii = 1,nl
            cst_nri(ii,jj,j,k,1) = cst_nri(ii,jj,j,k,1) +
     >                           dt(i,j,k) * am(ii,jj) / vol(i,j,k)
  210     continue
c
            do 220 ii = 1,nl
              cst_nri(ii,ii,j,k,1) = cst_nri(ii,ii,j,k,1) + 1.d0
c
            if ( dual_t.eq.1 ) then
              if ( dual_time.lt.1e-10 ) then
                cst_nri(ii,ii,j,k,1) = cst_nri(ii,ii,j,k,1)
     >                             + dt(i,j,k)/tintvl
              else
                cst_nri(ii,ii,j,k,1) = cst_nri(ii,ii,j,k,1)
     >                             + 1.5d0*dt(i,j,k)/tintvl
              end if
            end if
  220     continue
  230   continue
c
c ... start Gauss-Seidel Relaxation Sweeping
c ... eta-constant line
c
        do 280 k = 1, kl
c       
          do 240 j  = 1,jl
            capb(1,1,j) = cst_nri(5,5,j,k,2) 
            capc(1,1,j) = cst_nri(5,5,j,k,3) 
            capa(1,1,j) = cst_nri(5,5,j,k,1)
  240     continue
c
c ... form right-hand-side vector f(n,j)
c
c ... initilize f(n,j)
c
          do 260 j = 1,jl
            f(1,j) = rhsbi(j,k,5,1)
  260     continue
c
c ... solve linear system
c
          nblocks = jl
          call block1(nblocks, nl, dim1, 1, capa, capb, capc,
     >                soln, f)
c         
          do 270 j = 1,jl
            q_bi(j,k,5,1) = soln(1,j)
  270     continue
c
  280   continue
c
      end if
c
      return
      end
c
C======================================================================C
