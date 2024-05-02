C======================================================================C
clist()
      subroutine mass_flow(il, jl, kl, nl, x, y, z, time, dstep,
     $           ilower, iupper, jlower, jupper, klower, kupper, q)
c
c     This subroutine is used to calculate the mass flow at the inlet
c     of rotor.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! number of equations to solve
     $     dstep               ! number of time step
c
      double precision, intent(in)::
     $     time,
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
c     LOCAL VARIABLES
c
      double precision::
     $     mflow, mflowin, mflowout, mflowl, mflowr,
     $     da1, db1, dc1, da2, db2, dc2,
     $     axi, ayi, azi, rref, uref, lref, coef, ru, rv, rw
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     ist, ied,            ! i index for inlet
     $     ilp, jlp, klp        ! il+1, jl+1, kl+1
c
c --- SUBROUTINE START ---
c
c *** some parameters ***
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c new
c     rref = 1.27924681794566d0
c     uref = 156.732850287718d0
c original
      rref = 1.01720351410267d0
      uref = 201.055912842603d0
      lref = 0.1d0
c
      coef = -0.5d0*36.d0*rref*uref*lref*lref
c
c *** mass flow at inlet ***
c
csp
      ist =  92
      ied = 111
csp
      mflow = 0.d0
      j = jlp
      do 10 k =   1, kl
      do 10 i = ist, ied
c ... area vecter
        da1 = x(i+1,j,k+1) - x(i,j,k)
        db1 = y(i+1,j,k+1) - y(i,j,k)
        dc1 = z(i+1,j,k+1) - z(i,j,k)
        da2 = x(i+1,j,k) - x(i,j,k+1)
        db2 = y(i+1,j,k) - y(i,j,k+1)
        dc2 = z(i+1,j,k) - z(i,j,k+1)
c       
        axi = db1 * dc2 - dc1 * db2
        ayi = dc1 * da2 - da1 * dc2
        azi = da1 * db2 - db1 * da2
c
        ru = q(i,j-1,k,2)+q(i,j,k,2)
        rv = q(i,j-1,k,3)+q(i,j,k,3)
        rw = q(i,j-1,k,4)+q(i,j,k,4)
c
        mflow = mflow + ru*axi+rv*ayi+rw*azi
   10 continue
      mflowin  = 0.5d0*coef*mflow
c
c *** mass flow at outlet ***
c
csp
      ist =   1
      ied =  10
csp
      mflow = 0.d0
      j = jlp
      do 70 k =   1, kl
        do 30 i = ist, ied
c ... area vecter
          da1 = x(i+1,j,k+1) - x(i,j,k)
          db1 = y(i+1,j,k+1) - y(i,j,k)
          dc1 = z(i+1,j,k+1) - z(i,j,k)
          da2 = x(i+1,j,k) - x(i,j,k+1)
          db2 = y(i+1,j,k) - y(i,j,k+1)
          dc2 = z(i+1,j,k) - z(i,j,k+1)
c
          axi = db1 * dc2 - dc1 * db2
          ayi = dc1 * da2 - da1 * dc2
          azi = da1 * db2 - db1 * da2
c
          ru = q(i,j-1,k,2)+q(i,j,k,2)
          rv = q(i,j-1,k,3)+q(i,j,k,3)
          rw = q(i,j-1,k,4)+q(i,j,k,4)
c
          mflow = mflow + ru*axi+rv*ayi+rw*azi
   30   continue
c
        do 50 i = 193, il
c ... area vecter
          da1 = x(i+1,j,k+1) - x(i,j,k)
          db1 = y(i+1,j,k+1) - y(i,j,k)
          dc1 = z(i+1,j,k+1) - z(i,j,k)
          da2 = x(i+1,j,k) - x(i,j,k+1)
          db2 = y(i+1,j,k) - y(i,j,k+1)
          dc2 = z(i+1,j,k) - z(i,j,k+1)
c
          axi = db1 * dc2 - dc1 * db2
          ayi = dc1 * da2 - da1 * dc2
          azi = da1 * db2 - db1 * da2
c
          ru = q(i,j-1,k,2)+q(i,j,k,2)
          rv = q(i,j-1,k,3)+q(i,j,k,3)
          rw = q(i,j-1,k,4)+q(i,j,k,4)
c
          mflow = mflow + ru*axi+rv*ayi+rw*azi
   50   continue
   70 continue
      mflowout = -0.5d0*coef*mflow
c
c *** mass flow at left side ***
c
csp
      ist = 112
      ied = 192
csp
      mflow = 0.d0
      j = jlp
      do 90 k =   1, kl
      do 90 i = ist, ied
c ... area vecter
        da1 = x(i+1,j,k+1) - x(i,j,k)
        db1 = y(i+1,j,k+1) - y(i,j,k)
        dc1 = z(i+1,j,k+1) - z(i,j,k)
        da2 = x(i+1,j,k) - x(i,j,k+1)
        db2 = y(i+1,j,k) - y(i,j,k+1)
        dc2 = z(i+1,j,k) - z(i,j,k+1)
c
        axi = db1 * dc2 - dc1 * db2
        ayi = dc1 * da2 - da1 * dc2
        azi = da1 * db2 - db1 * da2
c
        ru = q(i,j-1,k,2)+q(i,j,k,2)
        rv = q(i,j-1,k,3)+q(i,j,k,3)
        rw = q(i,j-1,k,4)+q(i,j,k,4)
c
        mflow = mflow + ru*axi+rv*ayi+rw*azi
   90 continue
      mflowl = 0.5d0*coef*mflow
c
c *** mass flow at left side ***
c
csp
      ist = 11
      ied = 91
csp
      mflow = 0.d0
      j = jlp
      do 110 k =   1, kl
      do 110 i = ist, ied
c ... area vecter
        da1 = x(i+1,j,k+1) - x(i,j,k)
        db1 = y(i+1,j,k+1) - y(i,j,k)
        dc1 = z(i+1,j,k+1) - z(i,j,k)
        da2 = x(i+1,j,k) - x(i,j,k+1)
        db2 = y(i+1,j,k) - y(i,j,k+1)
        dc2 = z(i+1,j,k) - z(i,j,k+1)
c
        axi = db1 * dc2 - dc1 * db2
        ayi = dc1 * da2 - da1 * dc2
        azi = da1 * db2 - db1 * da2
c
        ru = q(i,j-1,k,2)+q(i,j,k,2)
        rv = q(i,j-1,k,3)+q(i,j,k,3)
        rw = q(i,j-1,k,4)+q(i,j,k,4)
c
        mflow = mflow + ru*axi+rv*ayi+rw*azi
  110 continue
      mflowr = 0.5d0*coef*mflow
c
      write(24,1010) dstep, time, mflowin, mflowout, mflowl,
     >               mflowr
      write( *,1010) dstep, time, mflowin, mflowout, mflowl,
     >               mflowr
 1010 format( i6,1x,f7.2,1x,4(f12.6,1x) )
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine rotrslt(il, jl, kl, nl, x, y, z, time, dstep,
     $           ilower, iupper, jlower, jupper, klower, kupper, q,
     $           gamma, machinf, udj, vdj, wdj, vol)
c
c     This subroutine is used to calculate the computed results for
c     rotor 37.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! number of equations to solve
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     dstep               ! number of time step
c
      double precision, intent(in),
     $     dimension(ilower:iupper,klower:kupper,4)::
     $     udj, vdj, wdj
c
      double precision, intent(in)::
     $     time, gamma, machinf,
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(in):: ! mesh cell volume
     $     vol(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)
c
c     LOCAL VARIABLES
c
      double precision::
     $     gamma1, summ, sump, sumt, sumpt, sumtt, sumru, sumv,
     $     r, ri, u, v, w, et, p, t, pt, qq, psi, theta, as, mach2,
     $     pw(kl,2), tw(kl,2), ptw(kl,2), ttw(kl,2),
     $     ruw(kl), vw(kl), pp, tt, ppt, ttt, pep, eff, at1, at2,
     $     sum1, sum2, sum3, sum4, sum5, sum6, coef, cf1, cf2
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     ist, ied,            ! i index for inlet
     $     ilp, jlp, klp        ! il+1, jl+1, kl+1
c
c --- SUBROUTINE START ---
c
c *** some parameters ***
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
      j = jl
c
      gamma1 = gamma-1.d0
      cf1 = gamma*machinf*machinf
c
c *** inlet ***
c
csp
      ist =  92
      ied = 111
csp
      do 30 k = 1, kl
        summ  = 0.d0
        sump  = 0.d0
        sumt  = 0.d0
        sumpt = 0.d0
        sumtt = 0.d0
        do 10 i = ist, ied
          r  =    q(i,j,k,1)
          ri = 1.d0/r
          u  = ri*q(i,j,k,2)
          v  = ri*q(i,j,k,3)
          w  = ri*q(i,j,k,4)
          et =    q(i,j,k,5)
          qq = 0.5d0*( u*u+v*v+w*w )
          psi = vdj(i,k,2)*vdj(i,k,2)+wdj(i,k,2)*wdj(i,k,2)
          theta = qq-0.5d0*psi
          p = gamma1*( et-r*theta )
          t = cf1*p*ri
          as = gamma*p*ri
          v = v+vdj(i,k,2)
          w = w+wdj(i,k,2)
          qq = 0.5d0*( u*u+v*v+w*w )
          mach2 = qq/as
          coef = 1.d0+gamma1*mach2
          pt = p*coef**3.5
          tt = t*coef
c
          cf2 = r*vol(i,j,k)*u
          summ  = summ  + cf2
          sump  = sump  + cf2*p
          sumt  = sumt  + cf2*t
          sumpt = sumpt + cf2*pt
          sumtt = sumtt + cf2*tt
   10   continue
         pw(k,1) =  sump/summ
         tw(k,1) =  sumt/summ
        ptw(k,1) = sumpt/summ
        ttw(k,1) = sumtt/summ
   30 continue
c
c *** outlet ***
c
      do 90 k =   1, kl
        summ  = 0.d0
        sump  = 0.d0
        sumt  = 0.d0
        sumpt = 0.d0
        sumtt = 0.d0
        sumru = 0.d0
        sumv  = 0.d0
        do 50 i = 1, 10
          r  =    q(i,j,k,1)
          ri = 1.d0/r
          u  = ri*q(i,j,k,2)
          v  = ri*q(i,j,k,3)
          w  = ri*q(i,j,k,4)
          et =    q(i,j,k,5)
          qq = 0.5d0*( u*u+v*v+w*w )
          psi = vdj(i,k,2)*vdj(i,k,2)+wdj(i,k,2)*wdj(i,k,2)
          theta = qq-0.5d0*psi
          p = gamma1*( et-r*theta )
          t = cf1*p*ri
          as = gamma*p*ri
          v = v+vdj(i,k,2)
          w = w+wdj(i,k,2)
          qq = 0.5d0*( u*u+v*v+w*w )
          mach2 = qq/as
          coef = 1.d0+gamma1*mach2
          pt = p*coef**3.5
          tt = t*coef
c
          cf2 = r*vol(i,j,k)*u
          summ  = summ  + cf2
          sump  = sump  + cf2*p
          sumt  = sumt  + cf2*t
          sumpt = sumpt + cf2*pt
          sumtt = sumtt + cf2*tt
          sumru = sumru + cf2*r*u
          sumv  = sumv  + vol(i,j,k)
   50   continue
c
        do 70 i = 193, il
          r  =    q(i,j,k,1)
          ri = 1.d0/r
          u  = ri*q(i,j,k,2)
          v  = ri*q(i,j,k,3)
          w  = ri*q(i,j,k,4)
          et =    q(i,j,k,5)
          qq = 0.5d0*( u*u+v*v+w*w )
          psi = vdj(i,k,2)*vdj(i,k,2)+wdj(i,k,2)*wdj(i,k,2)
          theta = qq-0.5d0*psi
          p = gamma1*( et-r*theta )
          t = cf1*p*ri
          as = gamma*p*ri
          v = v+vdj(i,k,2)
          w = w+wdj(i,k,2)
          qq = 0.5d0*( u*u+v*v+w*w )
          mach2 = qq/as
          coef = 1.d0+gamma1*mach2
          pt = p*coef**3.5
          tt = t*coef
c
          cf2 = r*vol(i,j,k)*u
          summ  = summ  + cf2
          sump  = sump  + cf2*p
          sumt  = sumt  + cf2*t
          sumpt = sumpt + cf2*pt
          sumtt = sumtt + cf2*tt
          sumru = sumru + cf2*r*u
          sumv  = sumv  + vol(i,j,k)
   70   continue
         pw(k,2) =  sump/summ
         tw(k,2) =  sumt/summ
        ptw(k,2) = sumpt/summ
        ttw(k,2) = sumtt/summ
        ruw(k)   = sumru/summ
         vw(k)   = sumv
   90 continue
c
c *** averaging in radial direction ***
c
      at1 = gamma1/gamma
      at2 = gamma/gamma1
c
      sum1 = 0.d0
      sum2 = 0.d0
      sum3 = 0.d0
      sum4 = 0.d0
      sum5 = 0.d0
      sum6 = 0.d0
      do 110 k = 1, kl
        cf2 = ruw(k)*vw(k)
        sum1 = sum1+cf2
c
        coef = pw(k,2)/pw(k,1)
        sum2 = sum2+cf2*coef
c
        coef = tw(k,2)/tw(k,1)
        sum3 = sum3+cf2*coef
c
        coef = ptw(k,2)/ptw(k,1)
        sum4 = sum4+cf2*coef
c
        coef = ttw(k,2)/ttw(k,1)
        sum5 = sum5+cf2*coef
c
        coef = pw(k,2)/ptw(k,1)
        sum6 = sum6+cf2*coef
  110 continue
      pp  = sum2/sum1
      tt  = sum3/sum1
      ppt = sum4/sum1
      ttt = sum5/sum1
      pep = sum6/sum1
      eff = ( ppt**at1-1.d0 )/( ttt-1.d0 )
c
      write(25,1030) dstep, time, pp, tt, ppt, ttt, pep, eff
      write( *,1030) dstep, time, pp, tt, ppt, ttt, pep, eff
 1030 format( i6,1x,f7.2,1x,6(f10.6,1x) )
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      integer function itplk(igv,kgv)
c
c     finds the i index for tip clearance for linking.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     igv,                ! given i index
     $     kgv                 ! given k index
c
c     LOCAL VARIABLES
c
      integer::
     $     ilk(202), i, k
c
      k = kgv
      if ( kgv.gt.44) k = 44
      select case (k)
      case(40,41)
        data ilk /
     >    178, 177, 177, 177, 176, 175, 174, 173, 172, 172,
     >    171, 170, 169, 168, 167, 166, 165, 164, 164, 163,
     >    162, 161, 160, 159, 158, 157, 156, 155, 154, 153,
     >    152, 152, 151, 150, 149, 148, 147, 146, 145, 144,
     >    143, 142, 141, 140, 139, 138, 137, 136, 135, 134,
     >    133, 132, 131, 129, 128, 127, 126, 125, 124, 123,
     >    122, 121, 119, 118, 117, 116, 115, 114, 113, 111,
     >    110, 109, 107, 105, 103, 101,  98,  94,  92,  91,
     >     89,  88,  86,  85,  84,  83,  82,  82,  81,  80,
     >     80,  79,  79,  78,  78,  77,  77,  77,  77,  76,
     >     76,  76,  75,  75,  74,  74,  73,  72,  72,  71,
     >     70,  69,  69,  68,  67,  66,  65,  64,  63,  62,
     >     62,  61,  60,  59,  58,  57,  56,  55,  54,  53,
     >     53,  52,  51,  50,  49,  48,  47,  46,  45,  44,
     >     43,  42,  41,  40,  39,  38,  37,  36,  35,  34,
     >     33,  32,  31,  29,  28,  27,  26,  25,  24,  23,
     >     22,  21,  20,  19,  17,  16,  15,  14,  13,  12,
     >     11,  10,   8,   7,   6,   5,   3,   1, 200, 197,
     >    195, 193, 192, 190, 189, 187, 186, 185, 185, 184,
     >    183, 183, 182, 181, 181, 180, 180, 180, 179, 179,
     >    179, 178/
      case(42)
        data ilk /
     >    178, 177, 177, 177, 176, 175, 174, 173, 172, 172,
     >    171, 170, 169, 168, 167, 166, 165, 164, 164, 163,
     >    162, 161, 160, 159, 158, 157, 156, 155, 154, 153,
     >    152, 152, 151, 150, 149, 148, 147, 146, 145, 144,
     >    143, 142, 141, 140, 139, 138, 137, 136, 135, 134,
     >    133, 132, 131, 129, 128, 127, 126, 125, 124, 123,
     >    122, 121, 119, 118, 117, 116, 115, 114, 113, 111,
     >    110, 109, 107, 105, 103, 101,  98,  94,  92,  91,
     >     89,  88,  86,  85,  84,  83,  82,  82,  81,  80,
     >     80,  79,  79,  78,  78,  77,  77,  77,  77,  76,
     >     76,  76,  75,  75,  74,  74,  73,  72,  72,  71,
     >     70,  69,  69,  68,  67,  66,  65,  64,  63,  62,
     >     62,  61,  60,  59,  58,  57,  56,  55,  54,  53,
     >     53,  52,  51,  50,  49,  48,  47,  46,  45,  44,
     >     43,  42,  41,  40,  39,  38,  37,  36,  35,  34,
     >     33,  32,  30,  29,  28,  27,  26,  25,  24,  23,
     >     22,  21,  20,  19,  17,  16,  15,  14,  13,  12,
     >     11,  10,   8,   7,   6,   5,   3,   1, 200, 197,
     >    195, 193, 192, 190, 189, 187, 186, 185, 185, 184,
     >    183, 183, 182, 181, 181, 180, 180, 180, 179, 179,
     >    179, 178/
      case(43,44)
        data ilk /
     >    178, 177, 177, 177, 176, 175, 174, 173, 172, 172,
     >    171, 170, 169, 168, 167, 166, 165, 164, 164, 163,
     >    162, 161, 160, 159, 158, 157, 156, 155, 154, 153,
     >    152, 152, 151, 150, 149, 148, 147, 146, 145, 144,
     >    143, 142, 141, 140, 139, 138, 137, 136, 135, 134,
     >    133, 132, 131, 129, 128, 127, 126, 125, 124, 123,
     >    122, 121, 119, 118, 117, 116, 115, 114, 113, 111,
     >    110, 109, 107, 105, 103, 101,  98,  94,  92,  91,
     >     89,  88,  86,  85,  84,  83,  82,  82,  81,  80,
     >     80,  79,  79,  78,  78,  77,  77,  77,  77,  76,
     >     76,  76,  75,  75,  74,  74,  73,  72,  72,  71,
     >     70,  69,  69,  68,  67,  66,  65,  64,  63,  62,
     >     62,  61,  60,  59,  58,  57,  56,  55,  54,  54,
     >     53,  52,  51,  50,  49,  48,  47,  46,  45,  44,
     >     43,  42,  41,  40,  39,  38,  37,  36,  35,  34,
     >     33,  32,  30,  29,  28,  27,  26,  25,  24,  23,
     >     22,  21,  20,  19,  17,  16,  15,  14,  13,  12,
     >     11,  10,   8,   7,   6,   5,   3,   1, 200, 197,
     >    195, 193, 192, 190, 189, 187, 186, 185, 185, 184,
     >    183, 183, 182, 181, 181, 180, 180, 180, 179, 179,
     >    179, 178/
      case default
        write(*,*) 'k = ', k
        stop
      end select
c
      i = igv
      if ( igv.lt.  1 ) i = 202+igv
      if ( igv.gt.202 ) i = igv-202
      itplk = ilk(i)
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      integer function itplksm(igv)
c
c     finds the i index for tip clearance for linking.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     igv                 ! given i index
c
c     LOCAL VARIABLES
c
      integer::
     $     ilk(102), i
c
      data ilk /
     >  90,  90,  89,  88,  87,  86,  85,  85,  84,  83,
     >  82,  81,  80,  79,  78,  77,  76,  75,  74,  73,
     >  73,  72,  70,  69,  68,  67,  66,  65,  64,  63,
     >  62,  61,  60,  58,  57,  55,  54,  52,  49,  46,
     >  44,  43,  42,  41,  41,  40,  40,  39,  39,  39,
     >  38,  38,  37,  37,  36,  36,  35,  34,  33,  33,
     >  32,  31,  30,  29,  28,  27,  26,  25,  24,  23,
     >  22,  21,  20,  19,  18,  17,  16,  15,  14,  13,
     >  12,  11,  10,   9,   8,   6,   5,   4,   3,   1,
     > 101,  98,  97,  95,  94,  93,  93,  92,  92,  91,
     >  91,  90/
c
      i = igv
      if ( igv.lt.  1 ) i = 102+igv
      if ( igv.gt.102 ) i = igv-102
      itplksm = ilk(i)
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine prss_exit1(il, jl, kl, nl, ii, kk, x, y, z, q,
     $           ilower, iupper, jlower, jupper, klower, kupper,
     $           ronum, poutlet, p_exit)
c
c     This subroutine is used to set the exit hub static pressure
c     and solve radial equilibrium at the exit to get the pressure 
c     values at all radial locations. Specially for rotor 37 only.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! number of equations to solve
     $     ii, kk              ! index for which layer in radial direction
c
      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(in)::
     $     ronum, poutlet
c
      double precision, intent(out)::
     $     p_exit
c
c     LOCAL VARIABLES
c
      double precision::
     $     da1, db1, dc1, da2, db2, dc2, axj, ayj, azj, aaj,
     $     suma, sumr, sumv, sumi, yctr, zctr, vth, r, ri,
     $     v, w, rad, radl, dr, cfvw, cf1
c
      integer::
     $     i, j, k,             ! cell iteration index
     $     ig, kg,
     $     ilp, jlp, klp        ! il+1, jl+1, kl+1
c
c --- subroutine starts ---
c
c *** some parameters ***
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c
      ig = ii
      if ( ig.lt.1  ) ig =  1
      if ( ig.gt.il ) ig = il
c
      kg = kk
      if ( kg.lt.1  ) kg =  1
      if ( kg.gt.kl ) kg = kl
c
c *** integration ***
c
      j = jlp
      cfvw = 0.25d0*ronum
c
      sumi = poutlet
      do 50 k = 1, kg
        suma  = 0.d0
        sumr  = 0.d0
        sumv  = 0.d0
        do 10 i = 1, 10
          da1 = x(i+1,j,k+1) - x(i,j,k)
          db1 = y(i+1,j,k+1) - y(i,j,k)
          dc1 = z(i+1,j,k+1) - z(i,j,k)
          da2 = x(i+1,j,k) - x(i,j,k+1)
          db2 = y(i+1,j,k) - y(i,j,k+1)
          dc2 = z(i+1,j,k) - z(i,j,k+1)
c
          axj = db1*dc2-dc1*db2
          ayj = dc1*da2-da1*dc2
          azj = da1*db2-db1*da2
          aaj = dsqrt( axj*axj+ayj*ayj+azj*azj )
c
          r  = 0.5d0*( q(i,jl,k,1)+q(i,jlp,k,1) )
          ri = 0.5d0/r
          v  = ri*( q(i,jl,k,3)+q(i,jlp,k,3) )
          w  = ri*( q(i,jl,k,4)+q(i,jlp,k,4) )
          yctr = y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1)
          zctr = z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1)
          cf1 = dsqrt( yctr*yctr+zctr*zctr )
          vth = ( w*yctr-v*zctr )/cf1
          v  = v-cfvw*zctr
          w  = w+cfvw*yctr
c
          suma = suma + aaj
          sumr = sumr + aaj*r
          sumv = sumv + aaj*vth
   10   continue
        do 30 i = 193, il
          da1 = x(i+1,j,k+1) - x(i,j,k)
          db1 = y(i+1,j,k+1) - y(i,j,k)
          dc1 = z(i+1,j,k+1) - z(i,j,k)
          da2 = x(i+1,j,k) - x(i,j,k+1)
          db2 = y(i+1,j,k) - y(i,j,k+1)
          dc2 = z(i+1,j,k) - z(i,j,k+1)
c
          axj = db1*dc2-dc1*db2
          ayj = dc1*da2-da1*dc2
          azj = da1*db2-db1*da2
          aaj = dsqrt( axj*axj+ayj*ayj+azj*azj )
c
          r  = 0.5d0*( q(i,jl,k,1)+q(i,jlp,k,1) )
          ri = 0.5d0/r
          v  = ri*( q(i,jl,k,3)+q(i,jlp,k,3) )
          w  = ri*( q(i,jl,k,4)+q(i,jlp,k,4) )
          yctr = y(i,j,k)+y(i+1,j,k)+y(i,j,k+1)+y(i+1,j,k+1)
          zctr = z(i,j,k)+z(i+1,j,k)+z(i,j,k+1)+z(i+1,j,k+1)
          cf1 = dsqrt( yctr*yctr+zctr*zctr )
          vth = ( w*yctr-v*zctr )/cf1
          v  = v-cfvw*zctr
          w  = w+cfvw*yctr
c
          suma = suma + aaj
          sumr = sumr + aaj*r
          sumv = sumv + aaj*vth
   30   continue
        sumr = sumr/suma
        sumv = sumv/suma
c
        yctr = y(ig,j,k)+y(ig+1,j,k)+y(ig,j,k+1)+y(ig+1,j,k+1)
        zctr = z(ig,j,k)+z(ig+1,j,k)+z(ig,j,k+1)+z(ig+1,j,k+1)
        rad  = 0.25d0*dsqrt( yctr*yctr+zctr*zctr )
        if ( k.lt.2 ) then
          dr = 0.d0
        else
          yctr = y(ig,j,k-1)+y(ig+1,j,k-1)+y(ig,j,k)+y(ig+1,j,k)
          zctr = z(ig,j,k-1)+z(ig+1,j,k-1)+z(ig,j,k)+z(ig+1,j,k)
          radl = 0.25d0*dsqrt( yctr*yctr+zctr*zctr )
          dr = rad-radl
        end if
        cf1 = sumr*sumv*sumv/rad
        sumi = sumi + cf1*dr
   50 continue
      p_exit = sumi
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine prss_exit(il, jl, kl, nl, ii, kk, x, y, z, q,
     $           ilower, iupper, jlower, jupper, klower, kupper,
     $           ronum, poutlet, p_exit)
c
c     This subroutine is used to set the exit hub static pressure
c     and solve radial equilibrium at the exit to get the pressure 
c     values at all radial locations. Specially for rotor 37 only.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! number of equations to solve
     $     ii, kk              ! index for which layer in radial direction
c
      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(in)::
     $     ronum, poutlet
c
      double precision, intent(out)::
     $     p_exit
c
c     LOCAL VARIABLES
c
      double precision::
     $     r, ri, v, w, yctr, zctr, rad, radl, radm,
     $     vth, dr, sum, cfvw, cf1
c
      integer::
     $     j, k,               ! cell iteration index
     $     ig, kg,
     $     ilp, jlp, klp
c
c --- subroutine starts ---
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c
c *** some parameters ***
c
csp
      ig = ii
      if ( ig.lt.1  ) ig = il + ii
      if ( ig.gt.il ) ig = ii - il
csp
      kg = kk
      if ( kg.lt.1  ) kg =  1
      if ( kg.gt.kl ) kg = kl
c
c *** integration ***
c
      j = jlp
      cfvw = 0.25d0*ronum
c
      sum = poutlet
      do 10 k = 1, kg
        r  = 0.5d0*( q(ig,j,k,1)+q(ig,j-1,k,1) )
        ri = 0.5d0/r
        v  = ri*( q(ig,j,k,3)+q(ig,j-1,k,3) )
        w  = ri*( q(ig,j,k,4)+q(ig,j-1,k,4) )
        yctr = y(ig,j,k)+y(ig+1,j,k)+y(ig,j,k+1)+y(ig+1,j,k+1)
        zctr = z(ig,j,k)+z(ig+1,j,k)+z(ig,j,k+1)+z(ig+1,j,k+1)
        rad = dsqrt( yctr*yctr+zctr*zctr )
        v  = v-cfvw*zctr
        w  = w+cfvw*yctr
        vth = ( w*yctr-v*zctr )/rad
c
        if ( k.lt.2 ) then
          dr = 0.d0
        else
          yctr = y(ig,j,k-1)+y(ig+1,j,k-1)+y(ig,j,k)+y(ig+1,j,k)
          zctr = z(ig,j,k-1)+z(ig+1,j,k-1)+z(ig,j,k)+z(ig+1,j,k)
          radl = dsqrt( yctr*yctr+zctr*zctr )
          dr = 0.25d0*( rad-radl )
        end if
        radm = 0.125d0*( rad+radl )
        cf1 = r*vth*vth/rad
        sum = sum + cf1*dr
   10 continue
      p_exit = sum
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine uprf(il, jl, kl, nl, ii, kk, x, y, z, uk,
     $           ilower, iupper, jlower, jupper, klower, kupper)
c
c     This subroutine is used to calculate the u profile at inlet.
c     Specially for rotor 37 only.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl,                 ! number of equations to solve
     $     ii, kk              ! index for which layer in radial direction
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(out)::
     $     uk
c
c     LOCAL VARIABLES
c
      double precision::
     $     rfrm(13), ufrm(13),
     $     rad, yctr, zctr, x1, x2, y1, y2, xgv, ygt
c
      integer::
     $     ilp, jlp, klp,       ! il+1, jl+1, kl+1
     $     ig, kg,
     $     j, k,                ! cell iteration index
     $     kw
c
      data rfrm /
     >  1.75253000000000d0,
     >  1.77301314965986d0,
     >  1.86878571428571d0,
     >  1.96012948979592d0,
     >  2.04150848979592d0,
     >  2.12233389115646d0,
     >  2.19651610884354d0,
     >  2.26959112925170d0,
     >  2.34045175510204d0,
     >  2.40965158503401d0,
     >  2.47885141496599d0,
     >  2.54971204081633d0,
     >  2.56632000000000d0 /
c
      data ufrm /
     >  0.903671267984013d0,
     >  0.925141810253961d0,
     >  0.980010973832715d0,
     >  1.01293247197997d0,
     >  1.03297164476525d0,
     >  1.04322979273867d0,
     >  1.04346835431945d0,
     >  1.03440301424992d0,
     >  1.01841938833784d0,
     >  0.996471722906342d0,
     >  0.964027347920644d0,
     >  0.923710440769298d0,
     >  0.903671267984013d0 /
c
c --- subroutine starts ---
c
c *** some parameters ***
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c
      ig = ii
      if ( ig.lt.1  ) ig =  1
      if ( ig.gt.il ) ig = il
c
      kg = kk
      if ( kg.lt.1  ) kg =  1
      if ( kg.gt.kl ) kg = kl
c
c *** calculating ***
c
      j = jlp
c
      yctr = y(ig,j,kg)+y(ig+1,j,kg)+y(ig,j,kg+1)+y(ig+1,j,kg+1)
      zctr = z(ig,j,kg)+z(ig+1,j,kg)+z(ig,j,kg+1)+z(ig+1,j,kg+1)
      rad  = 0.25d0*dsqrt( yctr*yctr+zctr*zctr )
c
      do 10 k = 1, 12
        if ( rad.ge.rfrm(k).and.rad.le.rfrm(k+1) ) then
          kw = k
          go to 30
        end if
   10 continue
      uk = 0.d0
      return
   30 continue
c
      x1 = rfrm(kw)
      x2 = rfrm(kw+1)
      y1 = ufrm(kw)
      y2 = ufrm(kw+1)
      xgv = rad
      call lintp( x1,x2,y1,y2,xgv,ygt )
      uk = ygt
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      double precision function ptb(ptotal,rgv)
c
c     find total pressure per reference Pt when radius is given.
c     Bith ptb & rgv are dimensionless quantities.
c
      implicit none
c
c     interface variables
c
      double precision, intent(in)::
     $     ptotal,             ! dimensionless total pressure
     $     rgv                 ! given temperature
c
c     local variables
c
      integer::
     $     np,
     $     i, ii
c
      double precision::
     $     radius(20), ptpr(20),
     $     x1, x2, y1, y2, xgv, ygt
c
      data radius /
     > 1.70000d0, 1.79222d0, 1.83490d0, 1.87452d0, 1.91414d0, 1.95682d0,
     > 1.99644d0, 2.05435d0, 2.10922d0, 2.16713d0, 2.22504d0, 2.27990d0,
     > 2.32258d0, 2.36220d0, 2.40182d0, 2.44450d0, 2.48412d0, 2.51765d0,
     > 2.54203d0, 2.60000d0/
c
      data ptpr /
     > 0.9800d0, 0.9864d0, 1.0040d0, 1.0054d0, 1.0054d0, 1.0054d0,
     > 1.0054d0, 1.0054d0, 1.0054d0, 1.0048d0, 1.0041d0, 1.0034d0,
     > 1.0041d0, 1.0048d0, 1.0041d0, 1.0020d0, 0.9959d0, 0.9762d0,
     > 0.9435d0, 0.9400d0/
c
c --- function starts ---
c
      np = 20
c
      if ( rgv.lt.radius(1) ) then
        write(*,*) 'rgv is less than ', radius(1)
        stop
      end if
      if ( rgv.gt.radius(np) ) then
        write(*,*) 'rgv is greater than ', radius(np)
        stop
      end if
c
      do 10 i = 2, np
        if ( rgv.ge.radius(i-1).and.rgv.le.radius(i) ) then
          ii = i
          go to 30
        end if
   10 continue
   30 continue
c
      x1 = radius(ii-1)
      x2 = radius(ii)
      y1 = ptpr(ii-1)
      y2 = ptpr(ii)
      xgv = rgv
      call lintp( x1, x2, y1, y2, xgv, ygt )
      ptb = ptotal*ygt
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      double precision function ttb(ttotal,rgv)
c
c     find total temperature per reference Tt when radius is given.
c     Bith ttb & rgv are dimensionless quantities.
c
      implicit none
c
c     interface variables
c
      double precision, intent(in)::
     $     ttotal,             ! dimensionless total temperature
     $     rgv                 ! given temperature
c
c     local variables
c
      integer::
     $     np,
     $     i, ii
c
      double precision::
     $     radius(20), ttpr(20),
     $     x1, x2, y1, y2, xgv, ygt
c
      data radius /
     > 1.70000d0, 1.79222d0, 1.83490d0, 1.87452d0, 1.91414d0, 1.95682d0,
     > 1.99644d0, 2.05435d0, 2.10922d0, 2.16713d0, 2.22504d0, 2.27990d0,
     > 2.32258d0, 2.36220d0, 2.40182d0, 2.44450d0, 2.48412d0, 2.51765d0,
     > 2.54203d0, 2.60000d0/
c
      data ttpr /
     > 1.0005d0, 1.0004d0, 0.9990d0, 0.9987d0, 0.9988d0, 0.9988d0,
     > 0.9988d0, 0.9990d0, 0.9994d0, 0.9996d0, 1.0004d0, 1.0008d0,
     > 1.0006d0, 1.0002d0, 1.0000d0, 1.0002d0, 1.0004d0, 1.0004d0,
     > 1.0008d0, 1.0009d0/
c
c --- function starts ---
c
      np = 20
c
      if ( rgv.lt.radius(1) ) then
        write(*,*) 'rgv is less than ', radius(1)
        stop
      end if
      if ( rgv.gt.radius(np) ) then
        write(*,*) 'rgv is greater than ', radius(np)
        stop
      end if
c
      do 10 i = 2, np
        if ( rgv.ge.radius(i-1).and.rgv.le.radius(i) ) then
          ii = i
          go to 30
        end if
   10 continue
   30 continue
c
      x1 = radius(ii-1)
      x2 = radius(ii)
      y1 = ttpr(ii-1)
      y2 = ttpr(ii)
      xgv = rgv
      call lintp( x1, x2, y1, y2, xgv, ygt )
      ttb = ttotal*ygt
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine pbar(il, jl, kl, nl, x, y, z, q,
     $           ilower, iupper, jlower, jupper, klower, kupper,
     $           gamma, ronum, p_bar)
c
c     This subroutine is used to calculate the averaged pressure at
c     outlet cross-section. Specially for rotor 37 only.                         
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     nl                  ! number of equations to solve
c
      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
      double precision, intent(in)::
     $     gamma, ronum
c
      double precision, intent(out)::
     $     p_bar
c
c     LOCAL VARIABLES
c
      double precision::
     $     da1, db1, dc1, da2, db2, dc2, axj, ayj, azj, aaj,
     $     r, ri, u, v, w, et, qq, p, psi, theta, suma, sump,
     $     gamma1, vw, ww, yct, zct
c
      integer::
     $     js,
     $     i, j, k             ! cell iteration index
c
c --- function starts ---
c
      gamma1 = gamma-1.d0
c
      j  = jl
      js = jl+1
c
      suma = 0.d0
      sump = 0.d0
      do 50 k = 1, kl
        do 10 i = 1, 10
          da1 = x(i+1,js,k+1) - x(i,js,k)
          db1 = y(i+1,js,k+1) - y(i,js,k)
          dc1 = z(i+1,js,k+1) - z(i,js,k)
          da2 = x(i+1,js,k) - x(i,js,k+1)
          db2 = y(i+1,js,k) - y(i,js,k+1)
          dc2 = z(i+1,js,k) - z(i,js,k+1)
c
          axj = db1*dc2-dc1*db2
          ayj = dc1*da2-da1*dc2
          azj = da1*db2-db1*da2
          aaj = dsqrt( axj*axj+ayj*ayj+azj*azj )
c
          call yzsfj(il, jl, kl, x, y, z, ilower, iupper,
     >               jlower, jupper, klower, kupper,
     >               i, js, k, yct, zct)
          vw = -ronum*zct
          ww =  ronum*yct
          psi = vw*vw+ww*ww
c
          r  =    q(i,j,k,1)
          ri = 1.d0/r
          u  = ri*q(i,j,k,2)
          v  = ri*q(i,j,k,3)
          w  = ri*q(i,j,k,4)
          et =    q(i,j,k,5)
          qq = 0.5d0*( u*u+v*v+w*w )
          theta = qq-0.5d0*psi
          p = gamma1*( et-r*theta )
c
          suma = suma + aaj
          sump = sump + aaj*p
   10   continue
        do 30 i = 193, il
          da1 = x(i+1,js,k+1) - x(i,js,k)
          db1 = y(i+1,js,k+1) - y(i,js,k)
          dc1 = z(i+1,js,k+1) - z(i,js,k)
          da2 = x(i+1,js,k) - x(i,js,k+1)
          db2 = y(i+1,js,k) - y(i,js,k+1)
          dc2 = z(i+1,js,k) - z(i,js,k+1)
c
          axj = db1*dc2-dc1*db2
          ayj = dc1*da2-da1*dc2
          azj = da1*db2-db1*da2
          aaj = dsqrt( axj*axj+ayj*ayj+azj*azj )
c
          call yzsfj(il, jl, kl, x, y, z, ilower, iupper,
     >               jlower, jupper, klower, kupper,
     >               i, js, k, yct, zct)
          vw = -ronum*zct
          ww =  ronum*yct
          psi = vw*vw+ww*ww
c
          r  =    q(i,j,k,1)
          ri = 1.d0/r
          u  = ri*q(i,j,k,2)
          v  = ri*q(i,j,k,3)
          w  = ri*q(i,j,k,4)
          et =    q(i,j,k,5)
          qq = 0.5d0*( u*u+v*v+w*w )
          theta = qq-0.5d0*psi
          p = gamma1*( et-r*theta )
c
          suma = suma + aaj
          sump = sump + aaj*p
   30   continue
   50 continue
      p_bar = sump/suma
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine yzctr(il, jl, kl, x, y, z, ilower, iupper,
     $                 jlower, jupper, klower, kupper, 
     $                 igv, jgv, kgv, yct, zct)
c
c     This subroutine is used to calculate the y & z coordinates at
c     the center of each volume.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     igv, jgv, kgv       ! given index for each volume
c
      double precision, intent(out):: 
     $     yct, zct
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
c     LOCAL VARIABLES
c
      double precision::
     $     ymd(4), zmd(4),
     $     rad, sum, yctr, zctr, ang
c
      integer::
     $     i
c
c --- SUBROUTINE START ---
c
      ymd(1) = y(igv  ,jgv  ,kgv)+y(igv  ,jgv  ,kgv+1)
      ymd(2) = y(igv+1,jgv  ,kgv)+y(igv+1,jgv  ,kgv+1)
      ymd(3) = y(igv+1,jgv+1,kgv)+y(igv+1,jgv+1,kgv+1)
      ymd(4) = y(igv  ,jgv+1,kgv)+y(igv  ,jgv+1,kgv+1)
c
      zmd(1) = z(igv  ,jgv  ,kgv)+z(igv  ,jgv  ,kgv+1)
      zmd(2) = z(igv+1,jgv  ,kgv)+z(igv+1,jgv  ,kgv+1)
      zmd(3) = z(igv+1,jgv+1,kgv)+z(igv+1,jgv+1,kgv+1)
      zmd(4) = z(igv  ,jgv+1,kgv)+z(igv  ,jgv+1,kgv+1)
c
      sum  = 0.d0
      yctr = 0.d0
      zctr = 0.d0
      do 10 i = 1, 4
        rad = dsqrt( ymd(i)*ymd(i)+zmd(i)*zmd(i) )
        sum = sum + rad
        yctr = yctr + ymd(i)
        zctr = zctr + zmd(i)
   10 continue
      rad = sum/8.d0
c
      call angle(yctr,zctr,ang)
c
      yct = rad*dcos(ang)
      zct = rad*dsin(ang)
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine yzsfj(il, jl, kl, x, y, z, ilower, iupper,
     $                 jlower, jupper, klower, kupper, 
     $                 igv, jsf, kgv, yct, zct)
c
c     This subroutine is used to calculate the y & z coordinates at
c     j surface.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     igv, jsf, kgv       ! given index for each volume
c
      double precision, intent(out):: 
     $     yct, zct
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
c     LOCAL VARIABLES
c
      double precision::
     $     ymd(4), zmd(4),
     $     rad, sum, yctr, zctr, ang
c
      integer::
     $     i
c
c --- SUBROUTINE START ---
c
      ymd(1) = y(igv  ,jsf,kgv  )
      ymd(2) = y(igv+1,jsf,kgv  )
      ymd(3) = y(igv  ,jsf,kgv+1)
      ymd(4) = y(igv+1,jsf,kgv+1)
c
      zmd(1) = z(igv  ,jsf,kgv  )
      zmd(2) = z(igv+1,jsf,kgv  )
      zmd(3) = z(igv  ,jsf,kgv+1)
      zmd(4) = z(igv+1,jsf,kgv+1)
c
      sum  = 0.d0
      yctr = 0.d0
      zctr = 0.d0
      do 10 i = 1, 4
        rad = dsqrt( ymd(i)*ymd(i)+zmd(i)*zmd(i) )
        sum = sum + rad
        yctr = yctr + ymd(i)
        zctr = zctr + zmd(i)
   10 continue
      rad = sum/4.d0
c
      call angle(yctr,zctr,ang)
c
      yct = rad*dcos(ang)
      zct = rad*dsin(ang)
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      subroutine yzsfk(il, jl, kl, x, y, z, ilower, iupper,
     $                 jlower, jupper, klower, kupper, 
     $                 igv, jgv, ksf, yct, zct)
c
c     This subroutine is used to calculate the y & z coordinates at
c     k surface.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     il, jl, kl,         ! cell number in 3 directions
     $     igv, jgv, ksf       ! given index for each volume
c
      double precision, intent(out):: 
     $     yct, zct
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
c     LOCAL VARIABLES
c
      double precision::
     $     ymd(4), zmd(4),
     $     rad, sum, yctr, zctr, ang
c
      integer::
     $     i
c
c --- SUBROUTINE START ---
c
      ymd(1) = y(igv  ,jgv  ,ksf )
      ymd(2) = y(igv+1,jgv  ,ksf )
      ymd(3) = y(igv  ,jgv+1,ksf )
      ymd(4) = y(igv+1,jgv+1,ksf )
c
      zmd(1) = z(igv  ,jgv  ,ksf )
      zmd(2) = z(igv+1,jgv  ,ksf )
      zmd(3) = z(igv  ,jgv+1,ksf )
      zmd(4) = z(igv+1,jgv+1,ksf )
c
      sum  = 0.d0
      yctr = 0.d0
      zctr = 0.d0
      do 10 i = 1, 4
        rad = dsqrt( ymd(i)*ymd(i)+zmd(i)*zmd(i) )
        sum = sum + rad
        yctr = yctr + ymd(i)
        zctr = zctr + zmd(i)
   10 continue
      rad = sum/4.d0
c
      call angle(yctr,zctr,ang)
c
      yct = rad*dcos(ang)
      zct = rad*dsin(ang)
c
      return
      end
c
c======================================================================c
c======================================================================c
clist()
      double precision function angr(il, jl, kl, i, k, x, y, z,
     >          ilower, iupper, jlower, jupper, klower, kupper)
c
c     find tan( ang ) in radial direction.
c
      implicit none
c
c     interface variables
c
      integer, intent(in)::
     $     il, jl, kl,   
     $     i, k,
     $     ilower, iupper, jlower, jupper, klower, kupper
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c
c     local variables
c
      integer::
     $     jlp, klp
c
      double precision::
     $     ahub, acrn, rhub, rcrn, rad, dx, dr,
     $     ytmp, ztmp, r1, r2, x1, x2, y1, y2, xgv, ygt
c
c --- function starts ---
c
      jlp = jl + 1
      klp = kl + 1
c
c *** tan( ang ) at hub ***
c
      ytmp = y(i,jl,1)+y(i+1,jl,1)
      ztmp = z(i,jl,1)+z(i+1,jl,1)
      r2 = dsqrt( ytmp*ytmp+ztmp*ztmp )
      ytmp = y(i,jlp,1)+y(i+1,jlp,1)
      ztmp = z(i,jlp,1)+z(i+1,jlp,1)
      r1 = dsqrt( ytmp*ytmp+ztmp*ztmp )
      dr = 0.5d0*( r2-r1 )
c
      x2 = x(i,jl ,1)+x(i+1,jl ,1)
      x1 = x(i,jlp,1)+x(i+1,jlp,1)
      dx = 0.5d0*( x2-x1 )
c
      ahub = dr/dx
      y1 = y(i,jlp,1)+y(i+1,jlp,1)+y(i,jlp,2)+
     >     y(i+1,jlp,2)
      x1 = z(i,jlp,1)+z(i+1,jlp,1)+z(i,jlp,2)+
     >     z(i+1,jlp,2)
      rhub = 0.25d0*dsqrt( x1*x1+y1*y1 )
c
c *** tan( ang ) at crown ***
c
      ytmp = y(i,jl,klp)+y(i+1,jl,klp)
      ztmp = z(i,jl,klp)+z(i+1,jl,klp)
      r2 = dsqrt( ytmp*ytmp+ztmp*ztmp )
      ytmp = y(i,jlp,klp)+y(i+1,jlp,klp)
      ztmp = z(i,jlp,klp)+z(i+1,jlp,klp)
      r1 = dsqrt( ytmp*ytmp+ztmp*ztmp )
      dr = 0.5d0*( r2-r1 )
c
      x2 = x(i,jl ,klp)+x(i+1,jl ,klp)
      x1 = x(i,jlp,klp)+x(i+1,jlp,klp)
      dx = 0.5d0*( x2-x1 )
c
      acrn = dr/dx
      y1 = y(i,jlp,kl)+y(i+1,jlp,kl)+y(i,jlp,klp)+
     >     y(i+1,jlp,klp)
      x1 = z(i,jlp,kl)+z(i+1,jlp,kl)+z(i,jlp,klp)+
     >     z(i+1,jlp,klp)
      rcrn = 0.25d0*dsqrt( x1*x1+y1*y1 )
c
      y1 = y(i,jlp,k)+y(i+1,jlp,k)+y(i,jlp,k+1)+
     >     y(i+1,jlp,k+1)
      x1 = z(i,jlp,k)+z(i+1,jlp,k)+z(i,jlp,k+1)+
     >     z(i+1,jlp,k+1)
      rad = 0.25d0*dsqrt( x1*x1+y1*y1 )
c
      x1 = rhub
      x2 = rcrn
      y1 = ahub
      y2 = acrn
      xgv = rad
      call lintp( x1, x2, y1, y2, xgv, ygt )
      angr = ygt
c
      return
      end
c
c======================================================================c
