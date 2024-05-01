c======================================================================c
clist()
      subroutine wall_rtv(il, jl, kl, ilower, iupper, jlower, jupper,
     $     klower, kupper, x, y, z, udi, vdi, wdi, udj, vdj, wdj,
     $     udk, vdk, wdk, ronum)
c
c     compute wall boundary velocity for rotor flow (relative flow field)
c     including rotating hub from x = -0.264 (cm) to x = 4.521 (cm)
c
      implicit none
c
c     interface variables
c
      integer, intent(in)::
     $     il, jl, kl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper
c
      double precision, intent(in)::
     $     ronum
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,
     $     klower+1:kupper):: x, y, z
c
      double precision, dimension(jlower:jupper,
     $     klower:kupper,4), intent(out):: udi, vdi, wdi
c
      double precision, dimension(ilower:iupper,
     $     klower:kupper,4), intent(out):: udj, vdj, wdj
c
      double precision, dimension(ilower:iupper,
     $     jlower:jupper,4), intent(out):: udk, vdk, wdk
c
c     local variables
c
      integer::
     $     ilp, jlp, klp, i, j, k, ii
c
      double precision::
     $     xrst, xred, coef, xct, yct, zct,
     $     x1, x2, y1, y2, xgv, ygt
c
c --- subroutine start ---
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c
      coef = 0.25d0*ronum
csp
      xrst = -0.0264d0
      xred =  0.4521d0
csp
      udi = 0.d0
      vdi = 0.d0
      wdi = 0.d0
c
      udj = 0.d0
      vdj = 0.d0
      wdj = 0.d0
c
      udk = 0.d0
      vdk = 0.d0
      wdk = 0.d0
c
c *** i surface ***
c
      i = 1
      do 10 j = 1, jl
      do 10 k = 1, kl
        udi(j,k,1) = 0.d0 
        vdi(j,k,1) = 0.d0
        wdi(j,k,1) = 0.d0
   10 continue
c
      i = ilp
      do 30 j = 1, jl
      do 30 k = 1, kl
        udi(j,k,2) = 0.d0 
        vdi(j,k,2) = 0.d0
        wdi(j,k,2) = 0.d0
   30 continue
c
c *** j surface ***
c
      j = 1
      do 50 k = 1, kl
      do 50 i = 1, il
        udj(i,k,1) = 0.d0 
        vdj(i,k,1) = 0.d0
        wdj(i,k,1) = 0.d0
   50 continue
c
      j = jlp
      do 70 k = 1, kl
      do 70 i = 1, il
        call yzsfj(il, jl, kl, x, y, z, ilower, iupper,
     >             jlower, jupper, klower, kupper,
     >             i, j, k, yct, zct)
        udj(i,k,2) = 0.d0 
        vdj(i,k,2) =-ronum*zct
        wdj(i,k,2) = ronum*yct
   70 continue
c
c *** k surface ***
c
      k = 1
      do 90 j = 1, jl
      do 90 i = 1, il
        udk(i,j,1) = 0.d0 
        xct = 0.25d0*( x(i,j,k)+x(i+1,j,k)+x(i,j+1,k)+x(i+1,j+1,k) )
        if ( xct.gt.xrst.and.xct.lt.xred ) then
          vdk(i,j,1) = 0.d0
          wdk(i,j,1) = 0.d0
        else
c
c         vdk(i,j,1) = coef*( z(i,j,k)+z(i+1,j,k)+z(i,j+1,k)+
c    >                        z(i+1,j+1,k) )
c         wdk(i,j,1) =-coef*( y(i,j,k)+y(i+1,j,k)+y(i,j+1,k)+
c    >                        y(i+1,j+1,k) )
c
          call yzsfk(il, jl, kl, x, y, z, ilower, iupper,
     >               jlower, jupper, klower, kupper,
     >               i, j, k, yct, zct)
          vdk(i,j,1) = ronum*zct
          wdk(i,j,1) =-ronum*yct
        end if
   90 continue
c
      k = klp
      do 110 j = 1, jl
      do 110 i = 1, il
        udk(i,j,2) = 0.d0 
c
c       vdk(i,j,2) = coef*(z(i,j,k)+z(i+1,j,k)+z(i,j+1,k)+z(i+1,j+1,k))
c       wdk(i,j,2) =-coef*(y(i,j,k)+y(i+1,j,k)+y(i,j+1,k)+y(i+1,j+1,k))
c
        call yzsfk(il, jl, kl, x, y, z, ilower, iupper,
     >             jlower, jupper, klower, kupper,
     >             i, j, k, yct, zct)
        vdk(i,j,2) = ronum*zct
        wdk(i,j,2) =-ronum*yct
  110 continue
c
c *** boundary nodes ***
c
      do 730 ii = 1, 2
c
c *** i surfaces ***
c
c ... four sides ...
c
        do 170 k = 1, kl
          x1 = 1.d0
          x2 = 2.d0
          do 130 j = 0, jlower, -1
            udi(j,k,ii) = 0.d0
c
            xgv = dfloat(j)
            y1 = vdi(1,k,ii)
            y2 = vdi(2,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdi(j,k,ii) = ygt
c
            y1 = wdi(1,k,ii)
            y2 = wdi(2,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdi(j,k,ii) = ygt
  130     continue
c
          x1 = dfloat(jl)
          x2 = dfloat(jl-1)
          do 150 j = jlp, jupper
            udi(j,k,ii) = 0.d0
c
            xgv = dfloat(j)
            y1 = vdi(jl,k,ii)
            y2 = vdi(jl-1,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdi(j,k,ii) = ygt
c
            y1 = wdi(jl,k,ii)
            y2 = wdi(jl-1,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdi(j,k,ii) = ygt
  150     continue
  170   continue
c
        do 230 j = 1, jl
          x1 = 1.d0
          x2 = 2.d0
          do 190 k = 0, klower, -1
            udi(j,k,ii) = 0.d0
c
            xgv = dfloat(k)
            y1 = vdi(j,1,ii)
            y2 = vdi(j,2,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdi(j,k,ii) = ygt
c
            y1 = wdi(j,1,ii)
            y2 = wdi(j,2,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdi(j,k,ii) = ygt
  190     continue
c
          x1 = dfloat(kl)
          x2 = dfloat(kl-1)
          do 210 k = klp, kupper
            udi(j,k,ii) = 0.d0
c
            xgv = dfloat(k)
            y1 = vdi(j,kl,ii)
            y2 = vdi(j,kl-1,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdi(j,k,ii) = ygt
c
            y1 = wdi(j,kl,ii)
            y2 = wdi(j,kl-1,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdi(j,k,ii) = ygt
  210     continue
  230   continue
c
c ... four corners ...
c
        do 250 k = klower, 0
        do 250 j = jlower, 0
          udi(j,k,ii) = 0.d0
          vdi(j,k,ii) = 0.5d0*( vdi(1,k,ii)+vdi(j,1,ii) )
          wdi(j,k,ii) = 0.5d0*( wdi(1,k,ii)+wdi(j,1,ii) )
  250   continue
c
        do 270 k = klp, kupper
        do 270 j = jlower, 0
          udi(j,k,ii) = 0.d0
          vdi(j,k,ii) = 0.5d0*( vdi(1,k,ii)+vdi(j,kl,ii) )
          wdi(j,k,ii) = 0.5d0*( wdi(1,k,ii)+wdi(j,kl,ii) )
  270   continue
c
        do 290 k = klower, 0
        do 290 j = jlp, jupper
          udi(j,k,ii) = 0.d0
          vdi(j,k,ii) = 0.5d0*( vdi(jl,k,ii)+vdi(j,1,ii) )
          wdi(j,k,ii) = 0.5d0*( wdi(jl,k,ii)+wdi(j,1,ii) )
  290   continue
c
        do 310 k = klp, kupper
        do 310 j = jlp, jupper
          udi(j,k,ii) = 0.d0
          vdi(j,k,ii) = 0.5d0*( vdi(jl,k,ii)+vdi(j,kl,ii) )
          wdi(j,k,ii) = 0.5d0*( wdi(jl,k,ii)+wdi(j,kl,ii) )
  310   continue
c
c *** j surfaces ***
c
c ... four sides ...
c
        do 370 k = 1, kl
          x1 = 1.d0
          x2 = 2.d0
          do 330 i = 0, ilower, -1
            udj(i,k,ii) = 0.d0
c
            xgv = dfloat(i)
            y1 = vdj(1,k,ii)
            y2 = vdj(2,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdj(i,k,ii) = ygt
c
            y1 = wdj(1,k,ii)
            y2 = wdj(2,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdj(i,k,ii) = ygt
  330     continue
c
          x1 = dfloat(il)
          x2 = dfloat(il-1)
          do 350 i = ilp, iupper
            udj(i,k,ii) = 0.d0
c
            xgv = dfloat(i)
            y1 = vdj(il,k,ii)
            y2 = vdj(il-1,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdj(i,k,ii) = ygt
c
            y1 = wdj(il,k,ii)
            y2 = wdj(il-1,k,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdj(i,k,ii) = ygt
  350     continue
  370   continue
c
        do 430 i = 1, il
          x1 = 1.d0
          x2 = 2.d0
          do 390 k = 0, klower, -1
            udj(i,k,ii) = 0.d0
c
            xgv = dfloat(k)
            y1 = vdj(i,1,ii)
            y2 = vdj(i,2,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdj(i,k,ii) = ygt
c
            y1 = wdj(i,1,ii)
            y2 = wdj(i,2,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdj(i,k,ii) = ygt
  390     continue
c
          x1 = dfloat(kl)
          x2 = dfloat(kl-1)
          do 410 k = klp, kupper
            udj(i,k,ii) = 0.d0
c
            xgv = dfloat(k)
            y1 = vdj(i,kl,ii)
            y2 = vdj(i,kl-1,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdj(i,k,ii) = ygt
c
            y1 = wdj(i,kl,ii)
            y2 = wdj(i,kl-1,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdj(i,k,ii) = ygt
  410     continue
  430   continue
c
c ... four corners ...
c
        do 450 k = klower, 0
        do 450 i = ilower, 0
          udj(i,k,ii) = 0.d0
          vdj(i,k,ii) = 0.5d0*( vdj(1,k,ii)+vdj(i,1,ii) )
          wdj(i,k,ii) = 0.5d0*( wdj(1,k,ii)+wdj(i,1,ii) )
  450   continue
c
        do 470 k = klp, kupper
        do 470 i = ilower, 0
          udj(i,k,ii) = 0.d0
          vdj(i,k,ii) = 0.5d0*( vdj(1,k,ii)+vdj(i,kl,ii) )
          wdj(i,k,ii) = 0.5d0*( wdj(1,k,ii)+wdj(i,kl,ii) )
  470   continue
c
        do 490 k = klower, 0
        do 490 i = ilp, iupper
          udj(i,k,ii) = 0.d0
          vdj(i,k,ii) = 0.5d0*( vdj(il,k,ii)+vdj(i,1,ii) )
          wdj(i,k,ii) = 0.5d0*( wdj(il,k,ii)+wdj(i,1,ii) )
  490   continue
c
        do 510 k = klp, kupper
        do 510 i = ilp, iupper
          udj(i,k,ii) = 0.d0
          vdj(i,k,ii) = 0.5d0*( vdj(il,k,ii)+vdj(i,kl,ii) )
          wdj(i,k,ii) = 0.5d0*( wdj(il,k,ii)+wdj(i,kl,ii) )
  510   continue
c
c *** k surfaces ***
c
c ... four sides ...
c
        do 570 j = 1, jl
          x1 = 1.d0
          x2 = 2.d0
          do 530 i = 0, ilower, -1
            udk(i,j,ii) = 0.d0
c
            xgv = dfloat(i)
            y1 = vdk(1,j,ii)
            y2 = vdk(2,j,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdk(i,j,ii) = ygt
c
            y1 = wdk(1,j,ii)
            y2 = wdk(2,j,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdk(i,j,ii) = ygt
  530     continue
c
          x1 = dfloat(il)
          x2 = dfloat(il-1)
          do 550 i = ilp, iupper
            udk(i,j,ii) = 0.d0
c
            xgv = dfloat(i)
            y1 = vdk(il,j,ii)
            y2 = vdk(il-1,j,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdk(i,j,ii) = ygt
c
            y1 = wdk(il,j,ii)
            y2 = wdk(il-1,j,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdk(i,j,ii) = ygt
  550     continue
  570   continue
c
        do 630 i = 1, il
          x1 = 1.d0
          x2 = 2.d0
          do 590 j = 0, jlower, -1
            udk(i,j,ii) = 0.d0
c
            xgv = dfloat(j)
            y1 = vdk(i,1,ii)
            y2 = vdk(i,2,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdk(i,j,ii) = ygt
c
            y1 = wdk(i,1,ii)
            y2 = wdk(i,2,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdk(i,j,ii) = ygt
  590     continue
c
          x1 = dfloat(jl)
          x2 = dfloat(jl-1)
          do 610 j = jlp, jupper
            udk(i,j,ii) = 0.d0
c
            xgv = dfloat(j)
            y1 = vdk(i,jl,ii)
            y2 = vdk(i,jl-1,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            vdk(i,j,ii) = ygt
c
            y1 = wdk(i,jl,ii)
            y2 = wdk(i,jl-1,ii)
            call lintp(x1, x2, y1, y2, xgv, ygt)
            wdk(i,j,ii) = ygt
  610     continue
  630   continue
c
c ... four corners ...
c
        do 650 j = jlower, 0
        do 650 i = ilower, 0
          udk(i,j,ii) = 0.d0
          vdk(i,j,ii) = 0.5d0*( vdk(1,j,ii)+vdk(i,1,ii) )
          wdk(i,j,ii) = 0.5d0*( wdk(1,j,ii)+wdk(i,1,ii) )
  650   continue
c
        do 670 j = jlp, jupper
        do 670 i = ilower, 0
          udk(i,j,ii) = 0.d0
          vdk(i,j,ii) = 0.5d0*( vdk(1,j,ii)+vdk(i,jl,ii) )
          wdk(i,j,ii) = 0.5d0*( wdk(1,j,ii)+wdk(i,jl,ii) )
  670   continue
c
        do 690 j = jlower, 0
        do 690 i = ilp, iupper
          udk(i,j,ii) = 0.d0
          vdk(i,j,ii) = 0.5d0*( vdk(il,j,ii)+vdk(i,1,ii) )
          wdk(i,j,ii) = 0.5d0*( wdk(il,j,ii)+wdk(i,1,ii) )
  690   continue
c
        do 710 j = jlp, jupper
        do 710 i = ilp, iupper
          udk(i,j,ii) = 0.d0
          vdk(i,j,ii) = 0.5d0*( vdk(il,j,ii)+vdk(i,jl,ii) )
          wdk(i,j,ii) = 0.5d0*( wdk(il,j,ii)+wdk(i,jl,ii) )
  710   continue
c
  730 continue
c
      return
      end
c
c======================================================================c
