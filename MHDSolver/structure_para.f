      subroutine structure_para(nn, moving, time, vbd, xyz, dxyz, 
     $     dstep, tintvl, cd, ct, cl, xsl, xzl, strtp, cs, mus, ub, 
     $     kc, wh, walf, uinf, nfp, alf0, il, jl, kl, nmode, pj, wpw,
     $     hx, hy, hz, ux, uy, uz)
c     
c     compute structure parameters, the updated parameters
c     saved in xsl(:,1)
c     
      implicit none
c
c     interface variables
c
      integer, intent(in)::
     $     nn,
     $     moving,
     $     dstep,
     $     strtp,
     $     nfp,
     $     il, jl, kl,
     $     nmode
c
      double precision, intent(in)::
     $     time,
     $     tintvl,
     $     cd,
     $     ct,
     $     cl,
     $     cs, mus, ub,
     $     kc, walf, wh, uinf, alf0
c
      double precision, intent(inout)::
     $     xyz(3,3), vbd(3,3), xsl(4,3), xzl(2,3,6)
c
      double precision, intent(in)::
     $     pj(6),
     $     wpw(6),
     $     hx(il+1,kl+1,6), hy(il+1,kl+1,6), hz(il+1,kl+1,6)
c
      double precision, intent(out)::
     $     dxyz(3),
     $     ux(il+1,kl+1), uy(il+1,kl+1), uz(il+1,kl+1)
c
c     local variables
c
      double precision::
     $     amp,
     $     aop,
     $     kc1,
     $     tstr,
     $     dtow,
     $     resdl,
     $     ust, kcw, csw,
     $     dxyz1(3), cf1
c
      integer::
     $     i, k, ii, jj, ktp2,
     $     nfpw,
     $     istm, iedm, kstm, kedm ! dimension for middle section of the wing
c     
c     subroutine begin
c     
c ... convert xyz, vbd to xsl at 3 time steps
c
      select case (strtp)
c
      case (1)    !induced cylinder
        dtow = 5.0d0
        call strmodel(dtow, tintvl, cd, cl, xsl, cs, mus, ub,
     $                resdl)
c
        dxyz(1) = xsl(1,1)-xsl(1,2)
        dxyz(2) = xsl(3,1)-xsl(3,2)
c
        xyz(1,1) = xsl(1,1)
        xyz(2,1) = xsl(3,1)
c
        vbd(1,1) = xsl(2,1)
        vbd(2,1) = xsl(4,1)
      case (2)    !forced pitching airfoil
        amp = 0.00d0*0.017453293d0
        aop = 1.01d0*0.017453293d0
        kc1 = 0.404d0
        tstr = kc1*tintvl
c
        xsl(1,1) = 0.d0
        xsl(2,1) = 0.d0
        xsl(3,1) = amp+aop*dsin(kc1*time)
        xsl(4,1) = ( xsl(3,1)-xsl(3,2) )/tstr
c
        dxyz(1) = xsl(1,1)-xsl(1,2)
        dxyz(2) = xsl(3,1)-xsl(3,2)
        dxyz(3) = 0.d0
c
        xyz(1,1) = xsl(1,1)
        xyz(2,1) = xsl(3,1)
c
        vbd(1,1) = xsl(2,1)
        vbd(2,1) = xsl(4,1)
      case (3)    !induced airfoil
        amp = 0.017453293d0
c
c     --- case elaf_11:Mach = 0.925, V = 3.875 ---
c       kc = 0.06663197154765449d0
c       nfp = 629
c     --- case elaf_2 :Mach = 0.925, V = 4.000, dt = 0.3 ---
c       kc = 0.06454972243679028d0
c       nfp = 973
c     --- case elaf_8 :Mach = 0.925, V = 5.150, dt = 0.5 ---
c       kc = 0.05013570674702157d0
c       nfp = 752
c     --- case elaf_9 :Mach = 0.925, V = 5.125, dt = 0.5 ---
c       kc = 0.05038027117017777d0
c       nfp = 749
c
        tstr = kc*tintvl
        if ( dstep.le.nfp ) then
          xsl(1,1) = 0.d0
          xsl(2,1) = 0.d0
          xsl(3,1) = amp*dsin(kc*time)
          xsl(4,1) = ( xsl(3,1)-xsl(3,2) )/tstr
          resdl = 0.d0
        else
          dtow = 50.d0
          call strmodelaf(dtow, tintvl, cd, -cl, ct, xsl, resdl,
     >                    uinf, walf) 
        end if
        dxyz(1) = -0.5d0*(xsl(1,1)-xsl(1,2))
        dxyz(2) =         xsl(3,1)-xsl(3,2)
        dxyz(3) = 0.d0
c
        xyz(1,1) = -0.5d0*xsl(1,1)
        xyz(2,1) =        xsl(3,1)
c
        vbd(1,1) = xsl(2,1)
        vbd(2,1) = xsl(4,1)
      case (4)    !induced wing
        nfpw = 306
        amp = 0.001d0
        kcw = 0.815400088398293d0
        ust = 2.729080851553200E-002
        csw = 0.d0
        tstr = kcw*tintvl
        if ( dstep.le.nfpw ) then
          xzl(1,1,1) = amp*dsin( 0.205514735d0*time )
          xzl(2,1,1) = 0.205514735d0*amp*
     >                 dcos( 0.205514735d0*time )
c         xzl(2,1,1) = ( xzl(1,1,1)-xzl(1,2,1) )/tstr
        else
          dtow = 50.d0
          call strctr(nmode, dtow, tstr, pj, xzl, csw, wpw, ust,
     >                resdl)
        end if
        write(6,*) 'dstep = ',dstep,xzl(1,1,1),xzl(2,1,1)
c
        ktp2 = 25
        istm = ( il+2 )/2
        iedm = il+1
        kstm =    1
        kedm = ktp2
        do 140 k = kstm, kedm
        do 140 i = istm, iedm
          dxyz1 = 0.d0
          do 130 jj = 1, nmode
            cf1 = xzl(1,1,jj)-xzl(1,2,jj)
            dxyz1(1) = dxyz1(1) + cf1*hx(i,k,jj)
            dxyz1(2) = dxyz1(2) + cf1*hy(i,k,jj)
            dxyz1(3) = dxyz1(3) + cf1*hz(i,k,jj)
  130     continue
          ux(i,k) = dxyz1(1)
          uy(i,k) = dxyz1(2)
          uz(i,k) = dxyz1(3)
c
          ii = 2*istm - i
          ux(ii,k) = ux(i,k)
          uy(ii,k) = uy(i,k)
          uz(ii,k) = uz(i,k)
  140   continue
      case (5)    !induced NLR7301 airfoil
        amp = 0.017453293d0
c
        tstr = kc*tintvl
        if ( dstep.le.nfp ) then
          xsl(1,1) = 0.d0
          xsl(2,1) = 0.d0
          xsl(3,1) = alf0+amp*dsin(kc*time)
          xsl(4,1) = ( xsl(3,1)-xsl(3,2) )/tstr
          resdl = 0.d0
        else
          dtow = 50.d0
          call strafnlr(dtow, tintvl, cl, ct, mus,
     $                  walf, wh, xsl, resdl, uinf, alf0)
        end if
        dxyz(1) = xsl(1,2)-xsl(1,1)
        dxyz(2) = xsl(3,1)-xsl(3,2)
        dxyz(3) = 0.d0
c
        xyz(1,1) = -xsl(1,1)
        xyz(2,1) = xsl(3,1)*180.d0/3.1415926d0
c
        vbd(1,1) = -xsl(2,1)
        vbd(2,1) = xsl(4,1)

      end select
c
      return
      end
