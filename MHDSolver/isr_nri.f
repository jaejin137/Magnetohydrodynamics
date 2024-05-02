C======================================================================C
c This part is to enforced nrbc in i direction at minimum index side.  c
C======================================================================C
clist(1)
      subroutine isr_nri(jq, kq, il, jl, kl, nl, x, y, z,
     $               ilower, iupper, jlower, jupper, klower, kupper,
     $               gamma, q, bm, pinlet, machinf, ptotal, ttotal,
     $               mx2, sigma, extra)
c
c     This subroutine is used to calculate the source term added to
c     the rhs of the nrflbc equations in i direction.
c
      implicit none
c
c *** INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     jq, kq,             ! indexes for current point to calculate
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! equation number
     $     extra               ! = 0 internal flow; = 1 open field
c
      double precision, intent(in)::
     $     gamma, pinlet, machinf, ptotal, ttotal, mx2, sigma(2),
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
     $     gamma1, rr, uu, vv, ww, qq, pp, psh,
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
      j = jq
      kp = k + 1
      jp = j + 1
c
c ... jacobian and primitive variables at outlet
c
      cf0 = 64.d0
c
      do 10 i = 0, 2
        l = i + 1
c
c ... jacobian ( center of volume )
c
        if ( l.ge.2 ) then
          ip = i + 1
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
        pv(4,l) = 0.d0
        qq = 0.5d0 * ( pv(2,l)*pv(2,l) + pv(3,l)*pv(3,l) +
     >                 pv(4,l)*pv(4,l) )
        pv(5,l) = gamma1 * ( q(i,j,k,5) - pv(1,l) * qq )
   10 continue
      jcbn(1) = jcbn(2)
c
c ... compute metrics ( at mim i boundary )
c
      i = 1
      ip = i + 1
      xxi = x(ip,j,k) -x(i,j,k) +x(ip,jp,k) -x(i,jp,k)+
     >      x(ip,j,kp)-x(i,j,kp)+x(ip,jp,kp)-x(i,jp,kp)
      yxi = y(ip,j,k) -y(i,j,k) +y(ip,jp,k) -y(i,jp,k)+
     >      y(ip,j,kp)-y(i,j,kp)+y(ip,jp,kp)-y(i,jp,kp)
      zxi = z(ip,j,k) -z(i,j,k) +z(ip,jp,k) -z(i,jp,k)+
     >      z(ip,j,kp)-z(i,j,kp)+z(ip,jp,kp)-z(i,jp,kp)
c
      i = 1
      xet = x(i,jp,k)-x(i,j,k)+x(i,jp,kp)-x(i,j,kp)
      yet = y(i,jp,k)-y(i,j,k)+y(i,jp,kp)-y(i,j,kp)
      zet = z(i,jp,k)-z(i,j,k)+z(i,jp,kp)-z(i,j,kp)
c
      xzt = x(i,j,kp)-x(i,j,k)+x(i,jp,kp)-x(i,jp,k)
      yzt = y(i,j,kp)-y(i,j,k)+y(i,jp,kp)-y(i,jp,k)
      zzt = z(i,j,kp)-z(i,j,k)+z(i,jp,kp)-z(i,jp,k)
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
      qq = yet*zzt-yzt*zet
      xix = cf0*jcb*qq
c
      qq = xzt*zet-xet*zzt
      xiy = cf0*jcb*qq
c
      qq = xet*yzt-xzt*yet
      xiz = cf0*jcb*qq
c
      qq  = xix*xix + xiy*xiy + xiz*xiz
      bxi  = dsqrt( qq )
c
c ... the primitive variables at i = ictr
c
      rr = pv(1,ictr)
      uu = pv(2,ictr)
      vv = pv(3,ictr)
      ww = pv(4,ictr)
      ww = 0.d0
      pp = pv(5,ictr)
      qq = 0.5d0 * ( uu*uu+vv*vv+ww*ww )
c
c ... calculate derivative ppv/pxi
c
      alf = -1.d0/jcbn(1)
      do 20 i = 1, 5
        pvpxi(i) = pv(i,2)/jcbn(2)-pv(i,1)/jcbn(1)
   20 continue
c
c ... compute xli, plpv, & dxi
c
      if ( extra.eq.0 ) then
        psh = ptotal*( 1.d0-gamma1*qq*machinf*machinf/ttotal )
     >        **3.5d0
      else
        psh = pinlet
      end if
c
      call gdinl(gamma, xix, xiy, xiz, jcb, bxi, pv, pvpxi,
     >           alf, psh, mx2, sigma, xli, plpv, dxi)
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
      subroutine isl_nri(jq, kq, il, jl, kl, nl, x, y, z,
     $               ilower, iupper, jlower, jupper, klower, kupper,
     $               gamma, q, am, pinlet, machinf, ptotal, ttotal,
     $               mx2, sigma, extra)
c
c     This subroutine is used to calculate the coefficient term added
c     to the lhs of the nrflbc equations in i direction.
c
      implicit none
c
c *** INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     jq, kq,             ! indexes for current point to calculate
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! equation number
     $     extra               ! = 0 internal flow; = 1 open field
c
      double precision, intent(in)::
     $     gamma, pinlet, machinf, ptotal, ttotal, mx2, sigma(2),
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
     $     gamma1, rr, uu, vv, ww, qq, pp, cc, psh,
     $     xix, xiy, xiz, jcb, bxi,
     $     xxi, xet, xzt, yxi, yet, yzt, zxi, zet, zzt,
     $     cf0, cf1, cf2
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
      j = jq
      kp = k + 1
      jp = j + 1
c
c ... jacobian and primitive variables at outlet
c
      cf0 = 64.d0
c
      do 10 i = 0, 2
        l = i + 1
c
c ... jacobian ( center of volume )
c
        if ( l.ge.2 ) then
          ip = i + 1
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
        pv(4,l) = 0.d0
        qq = 0.5d0 * ( pv(2,l)*pv(2,l) + pv(3,l)*pv(3,l) +
     >                 pv(4,l)*pv(4,l) )
        pv(5,l) = gamma1 * ( q(i,j,k,5) - pv(1,l) * qq )
   10 continue
      jcbn(1) = jcbn(2)
c
c ... compute metrics ( at inlet boundary )
c
      i = 1
      ip = i + 1
      xxi = x(ip,j,k) -x(i,j,k) +x(ip,jp,k) -x(i,jp,k)+
     >      x(ip,j,kp)-x(i,j,kp)+x(ip,jp,kp)-x(i,jp,kp)
      yxi = y(ip,j,k) -y(i,j,k) +y(ip,jp,k) -y(i,jp,k)+
     >      y(ip,j,kp)-y(i,j,kp)+y(ip,jp,kp)-y(i,jp,kp)
      zxi = z(ip,j,k) -z(i,j,k) +z(ip,jp,k) -z(i,jp,k)+
     >      z(ip,j,kp)-z(i,j,kp)+z(ip,jp,kp)-z(i,jp,kp)
c
      i = 1
      xet = x(i,jp,k)-x(i,j,k)+x(i,jp,kp)-x(i,j,kp)
      yet = y(i,jp,k)-y(i,j,k)+y(i,jp,kp)-y(i,j,kp)
      zet = z(i,jp,k)-z(i,j,k)+z(i,jp,kp)-z(i,j,kp)
c
      xzt = x(i,j,kp)-x(i,j,k)+x(i,jp,kp)-x(i,jp,k)
      yzt = y(i,j,kp)-y(i,j,k)+y(i,jp,kp)-y(i,jp,k)
      zzt = z(i,j,kp)-z(i,j,k)+z(i,jp,kp)-z(i,jp,k)
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
      qq = yet*zzt-yzt*zet
      xix = cf0*jcb*qq
c
      qq = xzt*zet-xet*zzt
      xiy = cf0*jcb*qq
c
      qq = xet*yzt-xzt*yet
      xiz = cf0*jcb*qq
c
      qq  = xix*xix + xiy*xiy + xiz*xiz
      bxi  = dsqrt( qq )
c
c ... the primitive variables at i = ictr
c
      rr = pv(1,ictr)
      uu = pv(2,ictr)
      vv = pv(3,ictr)
      ww = pv(4,ictr)
      ww = 0.d0
      pp = pv(5,ictr)
      qq = 0.5d0 * ( uu*uu+vv*vv+ww*ww )
      cc = dsqrt(gamma*pp/rr)
c
c ... calculate derivative ppv/pxi
c
      alf = -1.d0/jcbn(1)
      do 20 i = 1, 5
        pvpxi(i) = pv(i,2)/jcbn(2)-pv(i,1)/jcbn(1)
   20 continue
c
c ... compute xli, plpv, & dxi
c
      if ( extra.eq.0 ) then
        psh = ptotal*( 1.d0-gamma1*qq*machinf*machinf/ttotal )
     >        **3.5d0
      else
        psh = pinlet
      end if
c
      call gdinl(gamma, xix, xiy, xiz, jcb, bxi, pv, pvpxi,
     >           alf, psh, mx2, sigma, xli, plpv, dxi)
c
c ... matrix pd/pv
c
c ... pdpv(1,i)
c
      cf1 = ( xix*plpv(1,1)+xiy*plpv(2,1)+xiz*plpv(3,1) )/bxi
      cf2 = rr*( plpv(4,1)+plpv(5,1) )
      cf0 = 1.d0/cc/1.4142136d0
      pdpv(1,1) = cf1+cf0*( xli(4)+xli(5)+cf2 )
c
      cf0 = rr/cc/1.4142136d0
      cf1 = ( xix*plpv(1,2)+xiy*plpv(2,2)+xiz*plpv(3,2) )/bxi
      cf2 = cf0*( plpv(4,2)+plpv(5,2) )
      pdpv(1,2) = cf1+cf2
c
      cf1 = ( xix*plpv(1,3)+xiy*plpv(2,3)+xiz*plpv(3,3) )/bxi
      cf2 = cf0*( plpv(4,3)+plpv(5,3) )
      pdpv(1,3) = cf1+cf2
c
      cf1 = ( xix*plpv(1,4)+xiy*plpv(2,4)+xiz*plpv(3,4) )/bxi
      cf2 = cf0*( plpv(4,4)+plpv(5,4) )
      pdpv(1,4) = cf1+cf2
c
      cf1 = ( xix*plpv(1,5)+xiy*plpv(2,5)+xiz*plpv(3,5) )/bxi
      cf2 = cf0*( plpv(4,5)+plpv(5,5) )
      pdpv(1,5) = cf1+cf2
c
c ... pdpv(2,i)
c
      cf0 = xix/1.4142136d0
      cf1 = 1.d0/bxi
      cf2 = cf0*( plpv(4,1)-plpv(5,1) )
      pdpv(2,1) = cf1*(-xiz*plpv(2,1)+xiy*plpv(3,1)+cf2)
c
      cf2 = cf0*( plpv(4,2)-plpv(5,2) )
      pdpv(2,2) = cf1*(-xiz*plpv(2,2)+xiy*plpv(3,2)+cf2)
c
      cf2 = cf0*( plpv(4,3)-plpv(5,3) )
      pdpv(2,3) = cf1*(-xiz*plpv(2,3)+xiy*plpv(3,3)+cf2)
c
      cf2 = cf0*( plpv(4,4)-plpv(5,4) )
      pdpv(2,4) = cf1*(-xiz*plpv(2,4)+xiy*plpv(3,4)+cf2)
c
      cf2 = cf0*( plpv(4,5)-plpv(5,5) )
      pdpv(2,5) = cf1*(-xiz*plpv(2,5)+xiy*plpv(3,5)+cf2)
c
c ... pdpv(3,i)
c
      cf0 = xiy/1.4142136d0
      cf1 = 1.d0/bxi
      cf2 = cf0*( plpv(4,1)-plpv(5,1) )
      pdpv(3,1) = cf1*(xiz*plpv(1,1)-xix*plpv(3,1)+cf2)
c
      cf2 = cf0*( plpv(4,2)-plpv(5,2) )
      pdpv(3,2) = cf1*(xiz*plpv(1,2)-xix*plpv(3,2)+cf2)
c
      cf2 = cf0*( plpv(4,3)-plpv(5,3) )
      pdpv(3,3) = cf1*(xiz*plpv(1,3)-xix*plpv(3,3)+cf2)
c
      cf2 = cf0*( plpv(4,4)-plpv(5,4) )
      pdpv(3,4) = cf1*(xiz*plpv(1,4)-xix*plpv(3,4)+cf2)
c
      cf2 = cf0*( plpv(4,5)-plpv(5,5) )
      pdpv(3,5) = cf1*(xiz*plpv(1,5)-xix*plpv(3,5)+cf2)
c
c ... pdpv(4,i)
c
      cf0 = xiz/1.4142136d0
      cf1 = 1.d0/bxi
      cf2 = cf0*( plpv(4,1)-plpv(5,1) )
      pdpv(4,1) = cf1*(-xiy*plpv(1,1)+xix*plpv(2,1)+cf2)
c
      cf2 = cf0*( plpv(4,2)-plpv(5,2) )
      pdpv(4,2) = cf1*(-xiy*plpv(1,2)+xix*plpv(2,2)+cf2)
c
      cf2 = cf0*( plpv(4,3)-plpv(5,3) )
      pdpv(4,3) = cf1*(-xiy*plpv(1,3)+xix*plpv(2,3)+cf2)
c
      cf2 = cf0*( plpv(4,4)-plpv(5,4) )
      pdpv(4,4) = cf1*(-xiy*plpv(1,4)+xix*plpv(2,4)+cf2)
c
      cf2 = cf0*( plpv(4,5)-plpv(5,5) )
      pdpv(4,5) = cf1*(-xiy*plpv(1,5)+xix*plpv(2,5)+cf2)
c
c ... pdpv(5,i)
c
      cf2 = rr*( plpv(4,1)+plpv(5,1) )
      cf0 = cc/1.4142136d0
      pdpv(5,1) = cf0*( xli(4)+xli(5)+cf2 )
c
      cf0 = rr*cc/1.4142136d0
      pdpv(5,2) = cf0*( plpv(4,2)+plpv(5,2) )
c
      pdpv(5,3) = cf0*( plpv(4,3)+plpv(5,3) )
c
      pdpv(5,4) = cf0*( plpv(4,4)+plpv(5,4) )
c
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
      subroutine gdinl(gamma, xix, xiy, xiz, jcb, bxi, pv, pvpxi,
     >                 alf, pinlet, mx2, sigma, xli, plpv, dxi)
c
c This subroutine is used to calculate vector d for inflow at min index
c boundary.
c
      implicit none
c
c *** INTERFACE VARIABLES
c
      double precision, intent(in)::
     $     gamma, xix, xiy, xiz, jcb, bxi, pv(5,3), pvpxi(5),
     $     alf, pinlet, mx2, sigma(2)
c
      double precision, intent(out)::
     $     xli(5),              ! italic Ls in xi direction
     $     plpv(5,5),           ! pl/pv
     $     dxi(5)               ! d vector
c
c *** LOCAL VARIABLES
c
      double precision::
     $     gamma1, rr, uu, vv, ww, pp, c2, cc, dpc, capu,
     $     capc, dcpu, pupv(5), pfpv(5), cfx, amrp, lmd,
     $     aac, bbc, ck, qq, pppv(5),
     $     cf0, cf1, cf2        ! temporary variables
c
c *** SUBROUTINE START ***
c
c ... set up some parameters
c
      gamma1 = gamma - 1.d0
      dcpu = xix*pvpxi(2) + xiy*pvpxi(3) + xiz*pvpxi(4)
c
c --- outside ---
c
      rr = pv(1,1)
      uu = pv(2,1)
      vv = pv(3,1)
      ww = pv(4,1)
      pp = pv(5,1)
      qq = uu*uu+vv*vv+ww*ww
      c2 = gamma * pp / rr
      cc = dsqrt( c2 )
      dpc = pp-pinlet
c
c ... calculate italic Ls
c
      cf0 = sigma(1)*( 1.d0+mx2 )
      ck = cf0/rr/jcb/1.4142136d0
      xli(4) = ck*dpc
c
c ... compute pl/pv
c
c ... pp/pv
c
      pppv(1) =  0.5d0*gamma1*qq
      pppv(2) = -gamma1*uu
      pppv(3) = -gamma1*vv
      pppv(4) = -gamma1*ww
      pppv(5) =  gamma1
c
c ... plpv(4,i)
c
      plpv(4,1) = ck*pppv(1)
      plpv(4,2) = ck*pppv(2)
      plpv(4,3) = ck*pppv(3)
      plpv(4,4) = ck*pppv(4)
      plpv(4,5) = ck*pppv(5)
c
      aac = rr/cc
      bbc = rr*cc
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
      xli(4) = xli(4)+xli(5)
      cf0 = rr/cc/bxi/1.4142136d0
      xli(1) = -cf0*xix*( xli(4)+xli(5) )
      xli(2) = -cf0*xiy*( xli(4)+xli(5) )
      xli(3) = -cf0*xiz*( xli(4)+xli(5) )
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
c ... plpv(1,i)
c
      cf0 = xix*rr/cc/bxi/1.4142136d0
      plpv(1,1) = -cf0*( 2.d0*xli(5)/rr+2.d0*plpv(5,1)+plpv(4,1) )
      plpv(1,2) = -cf0*( 2.d0*plpv(5,2)+plpv(4,2) )
      plpv(1,3) = -cf0*( 2.d0*plpv(5,3)+plpv(4,3) )
      plpv(1,4) = -cf0*( 2.d0*plpv(5,4)+plpv(4,4) )
      plpv(1,5) = -cf0*( 2.d0*plpv(5,5)+plpv(4,5) )
c
c ... plpv(2,i)
c
      cf0 = xiy*rr/cc/bxi/1.4142136d0
      plpv(2,1) = -cf0*( 2.d0*xli(5)/rr+2.d0*plpv(5,1)+plpv(4,1) )
      plpv(2,2) = -cf0*( 2.d0*plpv(5,2)+plpv(4,2) )
      plpv(2,3) = -cf0*( 2.d0*plpv(5,3)+plpv(4,3) )
      plpv(2,4) = -cf0*( 2.d0*plpv(5,4)+plpv(4,4) )
      plpv(2,5) = -cf0*( 2.d0*plpv(5,5)+plpv(4,5) )
c
c ... plpv(3,i)
c
      cf0 = xiz*rr/cc/bxi/1.4142136d0
      plpv(3,1) = -cf0*( 2.d0*xli(5)/rr+2.d0*plpv(5,1)+plpv(4,1) )
      plpv(3,2) = -cf0*( 2.d0*plpv(5,2)+plpv(4,2) )
      plpv(3,3) = -cf0*( 2.d0*plpv(5,3)+plpv(4,3) )
      plpv(3,4) = -cf0*( 2.d0*plpv(5,4)+plpv(4,4) )
      plpv(3,5) = -cf0*( 2.d0*plpv(5,5)+plpv(4,5) )
c
c ... plpv(4,i)
c
      plpv(4,1) = plpv(4,1)+plpv(5,1)
      plpv(4,2) = plpv(4,2)+plpv(5,2)
      plpv(4,3) = plpv(4,3)+plpv(5,3)
      plpv(4,4) = plpv(4,4)+plpv(5,4)
      plpv(4,5) = plpv(4,5)+plpv(5,5)
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
