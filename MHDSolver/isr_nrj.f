C======================================================================C
c This part is to enforced nrbc in j direction at minimum index side.  c
C======================================================================C
clist(1)
      subroutine isr_nrj(iq, kq, il, jl, kl, nl, x, y, z,
     $                  ilower, iupper, jlower, jupper, klower, kupper,
     $                  bc_eta_lower, bc_eta_upper,
     $                  gamma, q, bm, poutlet, mx2)
c
c     This subroutine is used to calculate the source term added to
c     the rhs of the nrflbc equations in j direction.
c
      implicit none
c
c *** INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     iq, kq,             ! indexes for current point to calculate
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl                  ! equation number
c
      integer, intent(in)::
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl)
c
      double precision, intent(in)::
     $     gamma, poutlet, mx2,
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(out)::
     $     bm(5)               ! right-hand-side matrix 
c
c *** LOCAL VARIABLES
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     ictr,                ! the index for the variables used for solving the bc
     $     l,                   ! index for most right volumes
     $     ip, jp, kp           ! temporary index
c
      double precision::
     $     gamma1, rr, uu, vv, ww, qq, pp,
     $     xix, xiy, xiz, jcb, bxi,
     $     xxi, xet, xzt, yxi, yet, yzt, zxi, zet, zzt,
     $     cf0
c
      double precision::
     $     pv(5,3),             ! primitive variables at three volumes
     $     jcbn(3),             ! jacobian at three volumes
     $     pvpxi(5),            ! derivative ppv/pxi
     $     alf,                 ! constant used in difference approximation
     $     xli(5),              ! italic Ls in xi direction
     $     plpv(5,5),           ! pl/pv
     $     dxi(5)               ! vector d
c
c *** SUBROUTINE START ***
c
c ... set up some parameters
c
      gamma1 = gamma - 1.d0
c
      ictr = 1
c
c ... the loop starts
c
      k = kq
      i = iq
      kp = k + 1
      ip = i + 1
c
c ... jacobian and primitive variables at outlet
c
      cf0 = 64.d0
c
      do 10 j = 0, 2
        l = j + 1
c
c ... jacobian ( center of volume )
c
        if ( l.ge.2 ) then
          jp = j + 1
          xxi = x(ip,j,k) -x(i,j,k) +x(ip,jp,k) -x(i,jp,k)+
     >          x(ip,j,kp)-x(i,j,kp)+x(ip,jp,kp)-x(i,jp,kp)
          xet = x(i,jp,k) -x(i,j,k) +x(ip,jp,k) -x(ip,j,k)+
     >          x(i,jp,kp)-x(i,j,kp)+x(ip,jp,kp)-x(ip,j,kp)
          xzt = x(i,j,kp) -x(i,j,k) +x(ip,j,kp) -x(ip,j,k)+
     >          x(i,jp,kp)-x(i,jp,k)+x(ip,jp,kp)-x(ip,jp,k)
c
          yxi = y(ip,j,k) -y(i,j,k) +y(ip,jp,k) -y(i,jp,k)+
     >          y(ip,j,kp)-y(i,j,kp)+y(ip,jp,kp)-y(i,jp,kp)
          yet = y(i,jp,k) -y(i,j,k) +y(ip,jp,k) -y(ip,j,k)+
     >          y(i,jp,kp)-y(i,j,kp)+y(ip,jp,kp)-y(ip,j,kp)
          yzt = y(i,j,kp) -y(i,j,k) +y(ip,j,kp) -y(ip,j,k)+
     >          y(i,jp,kp)-y(i,jp,k)+y(ip,jp,kp)-y(ip,jp,k)
c
          zxi = z(ip,j,k) -z(i,j,k) +z(ip,jp,k) -z(i,jp,k)+
     >          z(ip,j,kp)-z(i,j,kp)+z(ip,jp,kp)-z(i,jp,kp)
          zet = z(i,jp,k) -z(i,j,k) +z(ip,jp,k) -z(ip,j,k)+
     >          z(i,jp,kp)-z(i,j,kp)+z(ip,jp,kp)-z(ip,j,kp)
          zzt = z(i,j,kp) -z(i,j,k) +z(ip,j,kp) -z(ip,j,k)+
     >          z(i,jp,kp)-z(i,jp,k)+z(ip,jp,kp)-z(ip,jp,k)
c
          xix = xxi*( yet*zzt-yzt*zet )
          xiy = xet*( yxi*zzt-yzt*zxi )
          xiz = xzt*( yxi*zet-yet*zxi )
c
          qq = xix - xiy + xiz
          jcbn(l) = cf0/qq
        end if
c
c ... variables
c
        rr = 1.d0 / q(i,j,k,1)
        pv(1,l) =      q(i,j,k,1)
        pv(2,l) = rr * q(i,j,k,2)
        pv(3,l) = rr * q(i,j,k,3)
        pv(4,l) = rr * q(i,j,k,4)
        qq = 0.5d0 * ( pv(2,l)*pv(2,l) + pv(3,l)*pv(3,l) +
     >                 pv(4,l)*pv(4,l) )
        pv(5,l) = gamma1 * ( q(i,j,k,5) - pv(1,l) * qq )
   10 continue
      jcbn(1) = jcbn(2)
c
c ... compute metrics ( at min j boundary )
c
      j = 1
      jp = j + 1
      xet = x(i,jp,k) -x(i,j,k) +x(ip,jp,k) -x(ip,j,k)+
     >      x(i,jp,kp)-x(i,j,kp)+x(ip,jp,kp)-x(ip,j,kp)
      yet = y(i,jp,k) -y(i,j,k) +y(ip,jp,k) -y(ip,j,k)+
     >      y(i,jp,kp)-y(i,j,kp)+y(ip,jp,kp)-y(ip,j,kp)
      zet = z(i,jp,k) -z(i,j,k) +z(ip,jp,k) -z(ip,j,k)+
     >      z(i,jp,kp)-z(i,j,kp)+z(ip,jp,kp)-z(ip,j,kp)
c
      j = 1
      xxi = x(ip,j,k)-x(i,j,k)+x(ip,j,kp)-x(i,j,kp)
      yxi = y(ip,j,k)-y(i,j,k)+y(ip,j,kp)-y(i,j,kp)
      zxi = z(ip,j,k)-z(i,j,k)+z(ip,j,kp)-z(i,j,kp)
c
      xzt = x(i,j,kp)-x(i,j,k)+x(ip,j,kp)-x(ip,j,k)
      yzt = y(i,j,kp)-y(i,j,k)+y(ip,j,kp)-y(ip,j,k)
      zzt = z(i,j,kp)-z(i,j,k)+z(ip,j,kp)-z(ip,j,k)
c
      xix = xxi*( yet*zzt-yzt*zet )
      xiy = xet*( yxi*zzt-yzt*zxi )
      xiz = xzt*( yxi*zet-yet*zxi )
c
      cf0 = 16.d0
      qq = xix - xiy + xiz
      jcb = cf0/qq
c
      cf0 = 0.25d0
      qq = yzt*zxi-yxi*zzt
      xix = cf0*jcb*qq
c
      qq = xxi*zzt-xzt*zxi
      xiy = cf0*jcb*qq
c
      qq = xzt*yxi-xxi*yzt
      xiz = cf0*jcb*qq
c
      qq  = xix*xix + xiy*xiy + xiz*xiz
      bxi = dsqrt( qq )
c
c ... the primitive variables at i = ictr
c
      rr = pv(1,ictr)
      uu = pv(2,ictr)
      vv = pv(3,ictr)
      ww = pv(4,ictr)
      pp = pv(5,ictr)
      qq = 0.5d0 * ( uu*uu+vv*vv+ww*ww )
c
c ... calculate derivative ppv/pxi
c
      do 20 i = 1, 5
        pvpxi(i) = pv(i,2)/jcbn(2)-pv(i,1)/jcbn(1)
   20 continue
c
c ... compute xli, plpv, & dxi
c
      if ( bc_eta_lower(1,1).eq.14 ) then
        call gdwall(gamma, xix, xiy, xiz, jcb, bxi, pv,
     >              pvpxi, alf, xli, plpv, dxi, 1)
      end if
c
      if ( bc_eta_lower(1,1).eq.15 ) then
        call gdsdl(gamma, xix, xiy, xiz, jcb, bxi, pv,
     >             pvpxi, alf, poutlet, mx2, xli, plpv, dxi)
      end if
c
c ... d matrix
c
      bm(1) = dxi(1)
      bm(2) = uu * dxi(1) + rr * dxi(2)
      bm(3) = vv * dxi(1) + rr * dxi(3)
      bm(4) = ww * dxi(1) + rr * dxi(4)
      bm(5) = qq * dxi(1) + rr * ( uu * dxi(2) +
     >        vv * dxi(3) + ww * dxi(4) ) + dxi(5) / gamma1
c
      return
      end
c
C======================================================================C
C======================================================================C
clist(2)
      subroutine isl_nrj(iq, kq, il, jl, kl, nl, x, y, z,
     $                  ilower, iupper, jlower, jupper, klower, kupper,
     $                  bc_eta_lower, bc_eta_upper,
     $                  gamma, q, am, poutlet, mx2)
c
c     This subroutine is used to calculate the coefficient term added
c     to the lhs of the nrflbc equations in j direction.
c
      implicit none
c
c *** INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     iq, kq,             ! indexes for current point to calculate
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl                  ! equation number
c
      integer, intent(in)::
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl)
c
      double precision, intent(in)::
     $     gamma, poutlet, mx2,
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(out)::
     $     am(5,5)             ! left-hand-side matrix
c
c *** LOCAL VARIABLES
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     ictr,                ! the index for the variables used for solving the bc
     $     l,                   ! index for most right volumes
     $     ip, jp, kp           ! temporary index
c
      double precision::
     $     gamma1, rr, uu, vv, ww, qq, pp, cc,
     $     xix, xiy, xiz, jcb, bxi,
     $     xxi, xet, xzt, yxi, yet, yzt, zxi, zet, zzt,
     $     cf0
c
      double precision::
     $     pv(5,3),             ! primitive variables at three volumes
     $     jcbn(3),             ! jacobian at three volumes
     $     pvpxi(5),            ! derivative ppv/pxi
     $     alf,                 ! constant used in difference approximation
     $     xli(5),              ! italic Ls in xi direction
     $     plpv(5,5),           ! pl/pv
     $     dxi(5),              ! vector d
     $     pdpv(5,5)            ! pd/pv
c
c *** SUBROUTINE START ***
c
c ... set up some parameters
c
      gamma1 = gamma - 1.d0
c
      ictr = 1
c
c ... the loop starts
c
      k = kq
      i = iq
      kp = k + 1
      ip = i + 1
c
c ... jacobian and primitive variables at outlet
c
      cf0 = 64.d0
c
      do 10 j = 0, 2
        l = j + 1
c
c ... jacobian ( center of volume )
c
        if ( l.ge.2 ) then
          jp = j + 1
          xxi = x(ip,j,k) -x(i,j,k) +x(ip,jp,k) -x(i,jp,k)+
     >          x(ip,j,kp)-x(i,j,kp)+x(ip,jp,kp)-x(i,jp,kp)
          xet = x(i,jp,k) -x(i,j,k) +x(ip,jp,k) -x(ip,j,k)+
     >          x(i,jp,kp)-x(i,j,kp)+x(ip,jp,kp)-x(ip,j,kp)
          xzt = x(i,j,kp) -x(i,j,k) +x(ip,j,kp) -x(ip,j,k)+
     >          x(i,jp,kp)-x(i,jp,k)+x(ip,jp,kp)-x(ip,jp,k)
c
          yxi = y(ip,j,k) -y(i,j,k) +y(ip,jp,k) -y(i,jp,k)+
     >          y(ip,j,kp)-y(i,j,kp)+y(ip,jp,kp)-y(i,jp,kp)
          yet = y(i,jp,k) -y(i,j,k) +y(ip,jp,k) -y(ip,j,k)+
     >          y(i,jp,kp)-y(i,j,kp)+y(ip,jp,kp)-y(ip,j,kp)
          yzt = y(i,j,kp) -y(i,j,k) +y(ip,j,kp) -y(ip,j,k)+
     >          y(i,jp,kp)-y(i,jp,k)+y(ip,jp,kp)-y(ip,jp,k)
c
          zxi = z(ip,j,k) -z(i,j,k) +z(ip,jp,k) -z(i,jp,k)+
     >          z(ip,j,kp)-z(i,j,kp)+z(ip,jp,kp)-z(i,jp,kp)
          zet = z(i,jp,k) -z(i,j,k) +z(ip,jp,k) -z(ip,j,k)+
     >          z(i,jp,kp)-z(i,j,kp)+z(ip,jp,kp)-z(ip,j,kp)
          zzt = z(i,j,kp) -z(i,j,k) +z(ip,j,kp) -z(ip,j,k)+
     >          z(i,jp,kp)-z(i,jp,k)+z(ip,jp,kp)-z(ip,jp,k)
c
          xix = xxi*( yet*zzt-yzt*zet )
          xiy = xet*( yxi*zzt-yzt*zxi )
          xiz = xzt*( yxi*zet-yet*zxi )
c
          qq = xix - xiy + xiz
          jcbn(l) = cf0/qq
        end if
c
c ... variables
c
        rr = 1.d0 / q(i,j,k,1)
        pv(1,l) =      q(i,j,k,1)
        pv(2,l) = rr * q(i,j,k,2)
        pv(3,l) = rr * q(i,j,k,3)
        pv(4,l) = rr * q(i,j,k,4)
        qq = 0.5d0 * ( pv(2,l)*pv(2,l) + pv(3,l)*pv(3,l) +
     >                 pv(4,l)*pv(4,l) )
        pv(5,l) = gamma1 * ( q(i,j,k,5) - pv(1,l) * qq )
   10 continue
      jcbn(1) = jcbn(2)
c
c ... compute metrics ( at min j boundary )
c
      j = 1
      jp = j + 1
      xet = x(i,jp,k) -x(i,j,k) +x(ip,jp,k) -x(ip,j,k)+
     >      x(i,jp,kp)-x(i,j,kp)+x(ip,jp,kp)-x(ip,j,kp)
      yet = y(i,jp,k) -y(i,j,k) +y(ip,jp,k) -y(ip,j,k)+
     >      y(i,jp,kp)-y(i,j,kp)+y(ip,jp,kp)-y(ip,j,kp)
      zet = z(i,jp,k) -z(i,j,k) +z(ip,jp,k) -z(ip,j,k)+
     >      z(i,jp,kp)-z(i,j,kp)+z(ip,jp,kp)-z(ip,j,kp)
c
      j = 1
      xxi = x(ip,j,k)-x(i,j,k)+x(ip,j,kp)-x(i,j,kp)
      yxi = y(ip,j,k)-y(i,j,k)+y(ip,j,kp)-y(i,j,kp)
      zxi = z(ip,j,k)-z(i,j,k)+z(ip,j,kp)-z(i,j,kp)
c
      xzt = x(i,j,kp)-x(i,j,k)+x(ip,j,kp)-x(ip,j,k)
      yzt = y(i,j,kp)-y(i,j,k)+y(ip,j,kp)-y(ip,j,k)
      zzt = z(i,j,kp)-z(i,j,k)+z(ip,j,kp)-z(ip,j,k)
c
      xix = xxi*( yet*zzt-yzt*zet )
      xiy = xet*( yxi*zzt-yzt*zxi )
      xiz = xzt*( yxi*zet-yet*zxi )
c
      cf0 = 16.d0
      qq = xix - xiy + xiz
      jcb = cf0/qq
c
      cf0 = 0.25d0
      qq = yzt*zxi-yxi*zzt
      xix = cf0*jcb*qq
c
      qq = xxi*zzt-xzt*zxi
      xiy = cf0*jcb*qq
c
      qq = xzt*yxi-xxi*yzt
      xiz = cf0*jcb*qq
c
      qq  = xix*xix + xiy*xiy + xiz*xiz
      bxi = dsqrt( qq )
c
c ... the primitive variables at i = ictr
c
      rr = pv(1,ictr)
      uu = pv(2,ictr)
      vv = pv(3,ictr)
      ww = pv(4,ictr)
      pp = pv(5,ictr)
      qq = 0.5d0 * ( uu*uu+vv*vv+ww*ww )
      cc = dsqrt(gamma*pp/rr)
c
c ... calculate derivative ppv/pxi
c
      alf = -1.d0/jcbn(1)
c
      do 20 i = 1, 5
        pvpxi(i) = pv(i,2)/jcbn(2)-pv(i,1)/jcbn(1)
   20 continue
c
c ... compute xli, plpv, & dxi
c
      if ( bc_eta_lower(1,1).eq.14 ) then
        call gdwall(gamma, xix, xiy, xiz, jcb, bxi, pv,
     >              pvpxi, alf, xli, plpv, dxi, 1)
      end if
c
      if ( bc_eta_lower(1,1).eq.15 ) then
        call gdsdl(gamma, xix, xiy, xiz, jcb, bxi, pv,
     >             pvpxi, alf, poutlet, mx2, xli, plpv, dxi)
      end if
c
c ... matrix pd/pv
c
      cf0 = 1.d0/cc/1.4142136d0
      pdpv(1,1) = cf0*( xli(4)+xli(5)+rr*( plpv(4,1)+plpv(5,1) ) )
      cf0 = rr/cc/1.4142136d0
      pdpv(1,2) = cf0*( plpv(4,2)+plpv(5,2) )
      pdpv(1,3) = cf0*( plpv(4,3)+plpv(5,3) )
      pdpv(1,4) = cf0*( plpv(4,4)+plpv(5,4) )
      pdpv(1,5) = cf0*( plpv(4,5)+plpv(5,5) )
c
      cf0 = xix/bxi/1.4142136d0
      pdpv(2,1) = cf0*( plpv(4,1)-plpv(5,1) )
      pdpv(2,2) = cf0*( plpv(4,2)-plpv(5,2) )
      pdpv(2,3) = cf0*( plpv(4,3)-plpv(5,3) )
      pdpv(2,4) = cf0*( plpv(4,4)-plpv(5,4) )
      pdpv(2,5) = cf0*( plpv(4,5)-plpv(5,5) )
c
      cf0 = xiy/bxi/1.4142136d0
      pdpv(3,1) = cf0*( plpv(4,1)-plpv(5,1) )
      pdpv(3,2) = cf0*( plpv(4,2)-plpv(5,2) )
      pdpv(3,3) = cf0*( plpv(4,3)-plpv(5,3) )
      pdpv(3,4) = cf0*( plpv(4,4)-plpv(5,4) )
      pdpv(3,5) = cf0*( plpv(4,5)-plpv(5,5) )
c
      cf0 = xiz/bxi/1.4142136d0
      pdpv(4,1) = cf0*( plpv(4,1)-plpv(5,1) )
      pdpv(4,2) = cf0*( plpv(4,2)-plpv(5,2) )
      pdpv(4,3) = cf0*( plpv(4,3)-plpv(5,3) )
      pdpv(4,4) = cf0*( plpv(4,4)-plpv(5,4) )
      pdpv(4,5) = cf0*( plpv(4,5)-plpv(5,5) )
c
      cf0 = cc/1.4142136d0
      pdpv(5,1) = cf0*( xli(4)+xli(5)+rr*( plpv(4,1)+plpv(5,1) ) )
      cf0 = rr*cc/1.4142136d0
      pdpv(5,2) = cf0*( plpv(4,2)+plpv(5,2) )
      pdpv(5,3) = cf0*( plpv(4,3)+plpv(5,3) )
      pdpv(5,4) = cf0*( plpv(4,4)+plpv(5,4) )
      pdpv(5,5) = cf0*( plpv(4,5)+plpv(5,5) )
c
c ... matrix theta = pd/pQ
c
      am(1,1) = pdpv(1,1)
      am(1,2) = pdpv(1,2)
      am(1,3) = pdpv(1,3)
      am(1,4) = pdpv(1,4)
      am(1,5) = pdpv(1,5)
c
      am(2,1) = uu*(pdpv(1,1)-dxi(1)/rr)+dxi(2)+rr*pdpv(2,1)
      am(2,2) = dxi(1)/rr+uu*pdpv(1,2)+rr*pdpv(2,2)
      am(2,3) = uu*pdpv(1,3)+rr*pdpv(2,3)
      am(2,4) = uu*pdpv(1,4)+rr*pdpv(2,4)
      am(2,5) = uu*pdpv(1,5)+rr*pdpv(2,5)
c
      am(3,1) = vv*(pdpv(1,1)-dxi(1)/rr)+dxi(3)+rr*pdpv(3,1)
      am(3,2) = vv*pdpv(1,2)+rr*pdpv(3,2)
      am(3,3) = dxi(1)/rr+vv*pdpv(1,3)+rr*pdpv(3,3)
      am(3,4) = vv*pdpv(1,4)+rr*pdpv(3,4)
      am(3,5) = vv*pdpv(1,5)+rr*pdpv(3,5)
c
      am(4,1) = ww*(pdpv(1,1)-dxi(1)/rr)+dxi(4)+rr*pdpv(4,1)
      am(4,2) = ww*pdpv(1,2)+rr*pdpv(4,2)
      am(4,3) = ww*pdpv(1,3)+rr*pdpv(4,3)
      am(4,4) = dxi(1)/rr+ww*pdpv(1,4)+rr*pdpv(4,4)
      am(4,5) = ww*pdpv(1,5)+rr*pdpv(4,5)
c
      am(5,1) = qq*( pdpv(1,1)-2.d0*dxi(1)/rr ) +
     >          rr*( uu*pdpv(2,1)+vv*pdpv(3,1)+ww*pdpv(4,1) ) +
     >          pdpv(5,1)/gamma1
      am(5,2) = uu*dxi(1)/rr + qq*pdpv(1,2) + dxi(2) +
     >          rr*( uu*pdpv(2,2)+vv*pdpv(3,2)+ww*pdpv(4,2) ) +
     >          pdpv(5,2)/gamma1
      am(5,3) = vv*dxi(1)/rr + qq*pdpv(1,3) + dxi(3) +
     >          rr*( uu*pdpv(2,3)+vv*pdpv(3,3)+ww*pdpv(4,3) ) +
     >          pdpv(5,3)/gamma1
      am(5,4) = ww*dxi(1)/rr + qq*pdpv(1,4) + dxi(4) +
     >          rr*( uu*pdpv(2,4)+vv*pdpv(3,4)+ww*pdpv(4,4) ) +
     >          pdpv(5,4)/gamma1
      am(5,5) = qq*pdpv(1,5) +
     >          rr*( uu*pdpv(2,5)+vv*pdpv(3,5)+ww*pdpv(4,5) ) +
     >          pdpv(5,5)/gamma1
c
      return
      end
c
C======================================================================C
C======================================================================C
clist(3)
      subroutine gdsdl(gamma, xix, xiy, xiz, jcb, bxi, pv, pvpxi,
     >                 alf, poutlet, mx2, xli, plpv, dxi)
c
c This subroutine is used to calculate vector d for sideflow at lower
c boundary. Absolutely no-reflecting boundary.
c
c *** INTERFACE VARIABLES
c
      double precision, intent(in)::
     $     gamma, xix, xiy, xiz, jcb, bxi, pv(5,3), pvpxi(5),
     $     alf, poutlet, mx2
c
      double precision, intent(out)::
     $     xli(5),              ! italic Ls in xi direction
     $     plpv(5,5),           ! pl/pv
     $     dxi(5)               ! d vector
c
c *** LOCAL VARIABLES
c
      double precision::
     $     gamma1, rr, uu, vv, ww, pp, c2, cc, capu,
     $     capc, dcpu, pupv(5), pfpv(5), cfx, amrp, lmd,
     $     aac, bbc,
     $     cf0, cf1, cf2        ! temporary variables
c
c *** SUBROUTINE START ***
c
c ... set up some parameters
c
      gamma1 = gamma - 1.d0
      dcpu = xix*pvpxi(2) + xiy*pvpxi(3) + xiz*pvpxi(4)
c
c *** min index boundary
c
      cf2 = xix*pv(2,2)+xiy*pv(3,2)+xiz*pv(4,2)
c
      if ( cf2.lt.0.d0 ) then
c
c --- inside ---
c
        rr = pv(1,2)
        uu = pv(2,2)
        vv = pv(3,2)
        ww = pv(4,2)
        pp = pv(5,2)
        c2 = gamma * pp / rr
        cc = dsqrt( c2 )
        capu = xix*uu + xiy*vv + xiz*ww
        capc = cc*bxi
        amrp = alf-pvpxi(1)/rr
c
c ... calculate italic Ls
c
        xli(1) = capu*( xix*pvpxi(1)+xiz*pvpxi(3)-xiy*pvpxi(4)-
     >                  xix*pvpxi(5)/c2 )/bxi
        xli(2) = capu*( xiy*pvpxi(1)-xiz*pvpxi(2)+xix*pvpxi(4)-
     >                  xiy*pvpxi(5)/c2 )/bxi
        xli(3) = capu*( xiz*pvpxi(1)+xiy*pvpxi(2)-xix*pvpxi(3)-
     >                  xiz*pvpxi(5)/c2 )/bxi
        lmd = capu-capc
        xli(5) = lmd*( pvpxi(5)/rr/cc-dcpu/bxi )/1.4142136d0
c
c ... compute pl/pv
c
        pupv(1) = -capu/rr
        pupv(2) = xix/rr
        pupv(3) = xiy/rr
        pupv(4) = xiz/rr
        pupv(5) = 0.d0
c
c ... plpv(1,i)
c
        pfpv(1) = xix*alf-xiz*( rr*pvpxi(3)+vv*pvpxi(1) )/rr/
     >            rr+( alf-2.d0*pvpxi(1)/rr )*(xiy*ww-xiz*vv)/
     >            rr+xiy*( rr*pvpxi(4)+ww*pvpxi(1) )/rr/rr
        pfpv(2) =  0.d0
        pfpv(3) =  xiz*amrp/rr
        pfpv(4) = -xiy*amrp/rr
        pfpv(5) =  0.d0
c
        cf0 = pvpxi(5)/c2
        cfx = xix*pvpxi(1)+xiz*pvpxi(3)-xiy*pvpxi(4)-xix*cf0
        plpv(1,1) = ( pupv(1)*cfx+capu*pfpv(1) )/bxi
        plpv(1,2) = ( pupv(2)*cfx+capu*pfpv(2) )/bxi
        plpv(1,3) = ( pupv(3)*cfx+capu*pfpv(3) )/bxi
        plpv(1,4) = ( pupv(4)*cfx+capu*pfpv(4) )/bxi
        plpv(1,5) = ( pupv(5)*cfx+capu*pfpv(5) )/bxi
c
c ... plpv(2,i)
c
        pfpv(1) = xiy*alf+xiz*( rr*pvpxi(2)+uu*pvpxi(1) )/rr/
     >            rr+( alf-2.d0*pvpxi(1)/rr )*(xiz*uu-xix*ww)/
     >            rr-xix*( rr*pvpxi(4)+ww*pvpxi(1) )/rr/rr
        pfpv(2) = -xiz*amrp/rr
        pfpv(3) =  0.d0
        pfpv(4) =  xix*amrp/rr
        pfpv(5) =  0.d0
c      
        cfx = xiy*pvpxi(1)-xiz*pvpxi(2)+xix*pvpxi(4)-xiy*cf0
        plpv(2,1) = ( pupv(1)*cfx+capu*pfpv(1) )/bxi
        plpv(2,2) = ( pupv(2)*cfx+capu*pfpv(2) )/bxi
        plpv(2,3) = ( pupv(3)*cfx+capu*pfpv(3) )/bxi
        plpv(2,4) = ( pupv(4)*cfx+capu*pfpv(4) )/bxi
        plpv(2,5) = ( pupv(5)*cfx+capu*pfpv(5) )/bxi
c
c ... plpv(3,i)
c
        pfpv(1) = xiz*alf-xiy*( rr*pvpxi(2)+uu*pvpxi(1) )/rr/
     >            rr+( alf-2.d0*pvpxi(1)/rr )*(xix*vv-xiy*uu)/
     >            rr+xix*( rr*pvpxi(3)+vv*pvpxi(1) )/rr/rr
        pfpv(2) =  xiy*amrp/rr
        pfpv(3) = -xix*amrp/rr
        pfpv(4) =  0.d0
        pfpv(5) =  0.d0
c
        cfx = xiz*pvpxi(1)+xiy*pvpxi(2)-xix*pvpxi(3)-xiz*cf0
        plpv(3,1) = ( pupv(1)*cfx+capu*pfpv(1) )/bxi
        plpv(3,2) = ( pupv(2)*cfx+capu*pfpv(2) )/bxi
        plpv(3,3) = ( pupv(3)*cfx+capu*pfpv(3) )/bxi
        plpv(3,4) = ( pupv(4)*cfx+capu*pfpv(4) )/bxi
        plpv(3,5) = ( pupv(5)*cfx+capu*pfpv(5) )/bxi
c
c ... plpv(5,i)
c
        cf0 = bxi*pvpxi(5)/rr/cc
        pfpv(1) =  dcpu/rr+capu*amrp/rr-cf0/rr
        pfpv(2) = -xix*amrp/rr
        pfpv(3) = -xiy*amrp/rr
        pfpv(4) = -xiz*amrp/rr
        pfpv(5) =  0.d0
c
        cf1 = 1.d0/bxi/1.4142136d0
        cfx = -dcpu+cf0
        plpv(5,1) = ( pupv(1)*cfx+lmd*pfpv(1) )*cf1
        plpv(5,2) = ( pupv(2)*cfx+lmd*pfpv(2) )*cf1
        plpv(5,3) = ( pupv(3)*cfx+lmd*pfpv(3) )*cf1
        plpv(5,4) = ( pupv(4)*cfx+lmd*pfpv(4) )*cf1
        plpv(5,5) = ( pupv(5)*cfx+lmd*pfpv(5) )*cf1
c
c --- outside ---
c
        rr = pv(1,1)
        uu = pv(2,1)
        vv = pv(3,1)
        ww = pv(4,1)
        pp = pv(5,1)
        c2 = gamma * pp / rr
        cc = dsqrt( c2 )
c
c ... calculate italic Ls
c
        xli(4) = 0.d0
c
c ... plpv(4,i)
c
        plpv(4,1) = 0.d0
        plpv(4,2) = 0.d0
        plpv(4,3) = 0.d0
        plpv(4,4) = 0.d0
        plpv(4,5) = 0.d0
c
        aac = rr/cc
        bbc = rr*cc
c
      else
c
c --- inside ---
c
        rr = pv(1,2)
        uu = pv(2,2)
        vv = pv(3,2)
        ww = pv(4,2)
        pp = pv(5,2)
        c2 = gamma * pp / rr
        cc = dsqrt( c2 )
        capu = xix*uu + xiy*vv + xiz*ww
        capc = cc*bxi
        amrp = alf-pvpxi(1)/rr
c
c ... calculate italic Ls
c
        lmd = capu-capc
        xli(5) = lmd*( pvpxi(5)/rr/cc-dcpu/bxi )/1.4142136d0
c
c ... compute pl/pv
c
        pupv(1) = -capu/rr
        pupv(2) = xix/rr
        pupv(3) = xiy/rr
        pupv(4) = xiz/rr
        pupv(5) = 0.d0
c
c ... plpv(5,i)
c
        cf0 = bxi*pvpxi(5)/rr/cc
        pfpv(1) =  dcpu/rr+capu*amrp/rr-cf0/rr
        pfpv(2) = -xix*amrp/rr
        pfpv(3) = -xiy*amrp/rr
        pfpv(4) = -xiz*amrp/rr
        pfpv(5) =  0.d0
c
        cf1 = 1.d0/bxi/1.4142136d0
        cfx = -dcpu+cf0
        plpv(5,1) = ( pupv(1)*cfx+lmd*pfpv(1) )*cf1
        plpv(5,2) = ( pupv(2)*cfx+lmd*pfpv(2) )*cf1
        plpv(5,3) = ( pupv(3)*cfx+lmd*pfpv(3) )*cf1
        plpv(5,4) = ( pupv(4)*cfx+lmd*pfpv(4) )*cf1
        plpv(5,5) = ( pupv(5)*cfx+lmd*pfpv(5) )*cf1
c
c --- outside ---
c
        rr = pv(1,1)
        uu = pv(2,1)
        vv = pv(3,1)
        ww = pv(4,1)
        pp = pv(5,1)
        c2 = gamma * pp / rr
        cc = dsqrt( c2 )
c
c ... calculate italic Ls
c
        xli(1) = 0.d0
        xli(2) = 0.d0
        xli(3) = 0.d0
        xli(4) = 0.d0
c
c ... plpv(1,i)
c
        plpv(1,1) = 0.d0
        plpv(1,2) = 0.d0
        plpv(1,3) = 0.d0
        plpv(1,4) = 0.d0
        plpv(1,5) = 0.d0
c
c ... plpv(2,i)
c
        plpv(2,1) = 0.d0
        plpv(2,2) = 0.d0
        plpv(2,3) = 0.d0
        plpv(2,4) = 0.d0
        plpv(2,5) = 0.d0
c
c ... plpv(3,i)
c
        plpv(3,1) = 0.d0
        plpv(3,2) = 0.d0
        plpv(3,3) = 0.d0
        plpv(3,4) = 0.d0
        plpv(3,5) = 0.d0
c
c ... plpv(4,i)
c
        plpv(4,1) = 0.d0
        plpv(4,2) = 0.d0
        plpv(4,3) = 0.d0
        plpv(4,4) = 0.d0
        plpv(4,5) = 0.d0
c
        aac = rr/cc
        bbc = rr*cc
c
      end if
c
c ... compute dxi
c
      cf1 = 1.d0/bxi
      cf2 = 1.d0/1.4142136d0
c
      dxi(1) = cf1*( xix*xli(1)+xiy*xli(2)+xiz*xli(3) )+
     >         cf2*aac*( xli(4)+xli(5) )
c
      cf0 = cf2*xix
      dxi(2) = cf1*( -xiz*xli(2)+xiy*xli(3)+cf0*( xli(4)-xli(5) ) )
c
      cf0 = cf2*xiy
      dxi(3) = cf1*(  xiz*xli(1)-xix*xli(3)+cf0*( xli(4)-xli(5) ) )
c
      cf0 = cf2*xiz
      dxi(4) = cf1*( -xiy*xli(1)+xix*xli(2)+cf0*( xli(4)-xli(5) ) )
c
      dxi(5) = cf2*bbc*( xli(4)+xli(5) )
c
      return
      end
c
C======================================================================C
