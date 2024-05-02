c======================================================================c
c
      subroutine turb_line(dist, num, uw, q, 
     $     omega, visturb, wake, pcheck, qtt,
     $     control)
c
c     compute B_L turbulent viscosity along a line normal to wall
c
      use datain 

      implicit none

      type (datain_type)::control
c
c     Interface Variables
c
      integer, intent(in)::
     $     num
c     
      double precision, dimension(num), intent(in)::
     $     dist, omega, q(num,5), qtt(num,3)
c     
      double precision, intent(in)::
     $     uw
      double precision::
     $     gamma, machinf, tref, reynolds, ronum
C     
      logical, intent(in)::
     $     wake, pcheck
      logical::
     $     suther
c
      double precision, dimension(num), intent(out)::
     $     visturb
C
c     Local Varibles
c      
      double precision, parameter::
     $     aplus = 26.d0,
     $     ccp = 1.6d0,
     $     ckleb = .3d0,
     $     cwk = .25d0,
     $     kcoef = .4d0,
     $     capk = .0168d0,
     $     cmutm = 14.d0
c
      integer::
     $     i
c
      double precision::
     $     qq, tw, vislamw, fmax, distmax, l, yplus, 
     $     stress, fwake, fkleb, udifmin, udifmax,
     $     psi, theta
c
      double precision, dimension(num)::
     $     turbin, turbout, fy, u, v, w, rho, rhoe, udif
c
      logical::
     $     max1st = .False.
c--------------------------parameter transfer
      gamma=control%gamma
      machinf=control%machinf
      reynolds=control%reynolds
      ronum=control%ronum
      tref=control%tref
 
      suther=control%suther

c--------------------------------------------------------      
c
c --- subroutine start ---
c
      rho = q(:,1)
      u = q(:,2)/rho
      v = q(:,3)/rho
      w = q(:,4)/rho
      rhoe = q(:,5)
c
      udif = 0.d0
      qq = 0.5d0*( u(1)*u(1)+v(1)*v(1)+w(1)*w(1) )
      if ( dabs(ronum).gt.1e-9 ) then
        psi = qtt(1,2)*qtt(1,2)+qtt(1,3)*qtt(1,3)
        theta = qq-0.5d0*psi
        tw = gamma*(gamma-1.d0)*machinf*machinf *
     >     ( rhoe(1)-rho(1)*theta )/rho(1)
      else
        tw = gamma*(gamma-1.d0)*machinf*machinf *
     >     ( rhoe(1)-rho(1)*qq )/rho(1)
      end if
c
      if ( suther ) then
         vislamw = dsqrt(tw)**3*(1.d0+tref)/(tw + tref)
      else
         vislamw = tw
      end if
c
      fmax = 0.d0
      stress = vislamw * uw / dist(1)
      do i = 1, num
         udif(i) = dsqrt(u(i)*u(i) +v(i)*v(i)+w(i)*w(i))
         if(wake) then
            fy(i) = dist(i)*omega(i)
         else
            yplus = dsqrt(rho(1)*stress)*dist(i)
     $           /vislamw*dsqrt(reynolds)
            l = kcoef*dist(i)*(1.d0-dexp(-yplus/aplus))
            fy(i) = dist(i)*omega(i)*(1.d0-dexp(-yplus/aplus))
            turbin(i) = rho(i)*l*l*omega(i)
            if(pcheck) print *, 1.d0-dexp(-yplus/aplus), yplus/aplus
         end if
      end do
c
      udifmin = minval(udif)
c
c$$$      if(.not.wake) udifmin = 0.d0
c
      if(max1st) then
         fmax = fy(1)
         distmax = dist(1)
         do i = 2, num
            if(fy(i).lt.fmax) exit
            fmax = fy(i)
            distmax = dist(i)
         end do
      else
         udifmax = maxval(udif)
         fmax = maxval(fy)
         distmax = sum(dist(maxloc(fy)))
      end if
c
      if (fmax.lt.1.d-12) fmax = 1.d-12
c
      do i = 1, num
         fwake = min(fmax*distmax,
     $        cwk*distmax*(udifmax-udifmin)**2/fmax)
         fkleb = 1.d0/(1.d0+5.5d0*(ckleb*dist(i)/distmax)**6)
         turbout(i) = capk*ccp*rho(i)*fwake*fkleb
      end do
c
      if(pcheck) then
         print *, 'omega  : '
         print '(10e12.4)', omega
         print *, 'dist   : '
         print '(10e12.4)', dist
         print *, 'distomg: '
         print '(10e12.4)', omega*dist
         print *, 'fy     : '
         print '(10e12.4)', fy
         print *, 'turbout: '
         print '(10e12.4)', turbout
         print *, 'turbin : '
         print '(10e12.4)', turbin
      end if
c
      visturb = turbout
c
      if (.not.wake) then
        do i = 1, num
c$$$            if(i.ge.2.and.(turbin(i)-turbout(i))*
c$$$     $           (turbin(i-1)-turbout(i-1)).lt.0.d0) then
c$$$               exit
c$$$            end if
c$$$            visturb(i) = turbin(i)
          if ( turbin(i).lt.visturb(i) ) then
            visturb(i) = turbin(i)
          else
            exit
          end if
        end do
      end if
c
      if(pcheck) then
         print *, 'visturb : '
         print '(10e12.4)', visturb
         print *, 'visturb : '
         print '(10e12.4)', visturb*reynolds
      end if
c
c$$$      if(maxval(visturb)*reynolds.lt.cmutm) visturb = 0.d0
c
      qq = 1000.d0/reynolds
      do i = 1, num
        if ( visturb(i).gt.qq ) visturb(i) = qq
      end do
c
      return
      end
c
c======================================================================c
