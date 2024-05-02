      subroutine invsid_flux(index,fl, il, jl, kl, nl, dim2,
     $     delta, capul, capur, capvl, capvr, capwl, capwr, el,
     $     er, pl, pr, ql, qr, rhoinvl, rhoinvr, rhol, rhor, rhoul,
     $     rhour, rhovl, rhovr, rhowl, rhowr, rhoel, rhoer,
     $     rhotkl, rhotkr, rhotwl, rhotwr, sqrtinv,
     $     sqrtrhol, sqrtrhor, tenthalpyl, tenthalpyr, ul, ur, vl, vr,
     $     wl, wr, lx, ly, lz, lxhat, lyhat, lzhat, mx, my, mz, mxhat,
     $     myhat, mzhat, nx, ny, nz, nxhat, nyhat, nzhat, bclower,
     $     bcupper, tkl, tkr, twl, twr, check, ke, qtl, qtr,
     $     ilower, vist,volumeinv,tl, tr, diver,qt1d,  control)

c..   This is to calculate the inviscid flux using 
c..   different upwind schemes. Zha 09/08/2002

C..   Scheme Index (rhs_scheme), Zha 10/15/2006

c     1: Roe Scheme   (Good for everything)
c     2: Zha2 Scheme (smooth, no contact discontinuity, AIAA J. No.5, 2005)
c     3: Zha ECUSPLD (Good for everything, 10/15/2006, call scheme 12 in 2004 version)
c     4: Van Leer Scheme
c     5  Edwards' LDFSS2 (Good for everything)
c     6  Zha6 Scheme (Contact discontinuity, not smooth, AIAA J., Jan. 2004)
c     7  Liou's AUSM+ Scheme
c     8  Zha's scheme-1999 (AIAA J, No. 8,1999)
c     9  Liou AUSMV
c     10 Wada-Liou, AUSMD
c     11 Van Leer-Hanel Scheme (1987)

c
c     IMPLICIT STATEMENT
c     
      use datain

      implicit none

      type (datain_type), intent(in)::control
c     
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     il, jl, kl, nl,
     $     dim2,
     $     index,
     $     bclower, bcupper,
     $     ilower
      integer::
     $     rhs_scheme, precondition,
     $     moving

      logical::
     $     inviscid,suther

      double precision, intent(in)::
     $     delta(5)
      double precision::
     $     machinf,
     $     gamma, tref,prandtl,prt,reynolds, k_prec,
     $     ronum
      
      double precision, dimension(dim2), intent(in)::
     $     capul, capur,
     $     capvl, capvr,
     $     capwl, capwr,
     $     el, er,
     $     pl, pr,
     $     ql, qr,
     $     rhoinvl, rhoinvr,
     $     rhol, rhor,
     $     rhoul, rhour,
     $     rhovl, rhovr,
     $     rhowl, rhowr,
     $     rhoel, rhoer,
     $     rhotkl, rhotkr,
     $     rhotwl, rhotwr,
     $     sqrtinv,
     $     sqrtrhol, sqrtrhor,
     $     tenthalpyl, tenthalpyr,
     $     ul, ur,
     $     vl, vr,
     $     wl, wr,
     $     tkl, tkr,
     $     twl, twr,
     $     qtl, qtr,
     $     tl,tr

      double precision, dimension(ilower:dim2)::
     $     diver
      double precision, dimension(ilower:dim2), intent(in)::
     $     vist, volumeinv

      double precision, dimension(ilower+1:dim2), intent(in)::
     $     lx, ly, lz,
     $     lxhat, lyhat, lzhat, 
     >     mx, my, mz,
     $     mxhat, myhat, mzhat,
     >     nx, ny, nz,
     $     nxhat, nyhat, nzhat

      double precision, intent(out)::
     $     fl(dim2, nl)

      logical, intent(in)::
     $     check,
     $     ke
c
      double precision, intent(in)::
     $     qt1d(ilower:dim2,3)
c     
c     LOCAL VARIABLES
c     
      double precision, dimension(dim2)::
     $     dq(dim2,nl),
     $     ax, ay, az,
     $     capu_l, capu_r

      double precision, dimension(nl,nl,dim2)::
     $     a

      integer::
     $     i,
     $     i_start,i_end

c..   Variables used for Zha and AUSM+ schemes

      double precision::
     $     area_,               ! interface area
     $     c_l, c_r,            ! critical speed of sound on left and right,
     $                          ! ausm+
     $     c_l_tilde,
     $     c_r_tilde,           ! interim speed of sound on left and right,
     $                          ! ausm+
     $     c_inter_face,        ! interface speed of sound, ausm+
     $     qml1, qmr1,          ! mach number on left and right to check if
     $                          ! the flow is supersonic or subsonic, this is
     $                          ! only used in zha scheme.
     $     c_mid,               ! intermediate speed of sound for pressure,
     $                          ! zha scheme
     $     qml0, qmr0,          ! mach number on left and right for pressue,
     $                          ! zha scheme
     $     plp, prn,            ! pressure spliting on left and right,
     $                          ! zha scheme
     $     htl, htr,            ! total enthalpy on left and right
     $     m_interface,         ! interface mach number, ausm+
     $     alfa_l, alfa_r,      ! used for to correct pressure oscillation, zha2 scheme
     $     omega_l, omega_r,    ! used for to correct pressure oscillation, zha2 scheme
     $     uu_l,uu_r,
     $     mid_rhou,
     $     mid_qt


C..   Variables for ECUSPLD, Zha3 Scheme

      double precision::
     >     QML,QMR,sign_l,sign_r,
     >     AFA_LP,AFA_RN,BTA_L,BTA_R,
     >     DETAP,DETAN,QMLP,QMLN,QMRP,QMRN,
     >     QM_MID,QM_MID_P,QM_MID_N,C_P,C_N,
     >     DLP,DRN,PHI,SL_P,SR_N,QML_BAR,QMR_BAR


C..   VARIABLES FOR VAN LEER SCHEME

      integer::
     $     iflux

      double precision::
     $     vx1,vy1,vz1,         ! interface unit vector components
     $     u_l,v_l,w_l,         ! left velocity components in x,y,z dir
     $     d_l,                 ! left density
     $     u_r,v_r,w_r,         ! right velocity components in x,y,z dir
     $     d_r,                 ! right density
     $     upo(nl),             ! calculated van leer flux
     $     qt_l, qt_r           ! grid velocities at left and right sides

      double precision, parameter::
     $     bta=0.d0,
     $     afr=0.d0             ! beta and alpha in ausm+
c     
c     *** SUBROUTINE START ***
c     
c-----------------------------------
      rhs_scheme=control%rhs_scheme
      moving=control%moving
      precondition=control%precondition

      inviscid=control%inviscid
      suther=control%suther

      gamma=control%gamma
      k_prec=control%k_prec
      machinf=control%machinf
      tref=control%tref
      prandtl=control%prandtl
      prt=control%prt
      reynolds=control%reynolds
      ronum=control%ronum

c----------------------------------------------


      select case (index)
      case (1)
         i_start =1
         i_end = il+1
         do i = i_start,i_end
            capu_l(i) = capul(i)
            capu_r(i) = capur(i)
            ax(i)    = lx(i)
            ay(i)    = ly(i)
            az(i)    = lz(i)
         end do
      case (2)
         i_start =1
         i_end = jl+1
         do i = i_start,i_end
            capu_l(i) = capvl(i)
            capu_r(i) = capvr(i)
            ax(i)    = mx(i)
            ay(i)    = my(i)
            az(i)    = mz(i)
         end do
      case (3)
         i_start =1
         i_end = kl+1
         do i = i_start,i_end
            capu_l(i) = capwl(i)
            capu_r(i) = capwr(i)
            ax(i)    = nx(i)
            ay(i)    = ny(i)
            az(i)    = nz(i)
         end do
      end select

      if ( moving.ge.1 ) then
         do i = i_start,i_end 
            capu_l(i) = capu_l(i) + qtl(i)
            capu_r(i) = capu_r(i) + qtr(i)
         end do
      end if


      SELECT CASE (rhs_scheme)
c     
c..   ROE SCHEME
c     
      CASE(1)                   ! Roe Scheme

c..   contribution to E from left and right faces
         do   i = i_start,i_end
            fl(i,1)  = rhol(i)*capu_l(i) + rhor(i)*capu_r(i)
            fl(i,2)  = rhoul(i)*capu_l(i) + ax(i)*pl(i) +
     >           rhour(i)*capu_r(i) + ax(i)*pr(i)
            fl(i,3)  = rhovl(i)*capu_l(i) + ay(i)*pl(i) +
     >           rhovr(i)*capu_r(i) + ay(i)*pr(i)
            fl(i,4)  = rhowl(i)*capu_l(i) + az(i)*pl(i) +
     >           rhowr(i)*capu_r(i) + az(i)*pr(i)
            fl(i,5)  = (rhoel(i) + pl(i))*capu_l(i) +
     >           (rhoer(i) + pr(i))*capu_r(i)
            if(nl.eq.6) then
               fl(i,6)  = rhol(i)*tkl(i)*capu_l(i) +
     >              rhor(i)*tkr(i)*capu_r(i)
            else if(nl.eq.7) then
               fl(i,6)  = rhol(i)*tkl(i)*capu_l(i) +
     >              rhor(i)*tkr(i)*capu_r(i)
               fl(i,7)  = rhol(i)*twl(i)*capu_l(i) +
     >              rhor(i)*twr(i)*capu_r(i)
            end if
           if ( moving.ge.1 ) then
             fl(i,5)  = fl(i,5) - qtl(i)*pl(i) - qtr(i)*pr(i)
           end if
         end do

c..   contribution from Roe matrix

         if (precondition.ge.1) then
            do   i = i_start,i_end
               dq(i,1)  = pl(i)    - pr(i)
               dq(i,2)  = ul(i)    - ur(i)
               dq(i,3)  = vl(i)    - vr(i)
               dq(i,4)  = wl(i)    - wr(i)
               dq(i,5)  = tl(i)    - tr(i)
               if(nl.eq.6) then
                  dq(i,6)  = rhol(i)*tkl(i) - rhor(i)*tkr(i)
               else if(nl.eq.7) then
                  dq(i,6)  = rhol(i)*tkl(i) - rhor(i)*tkr(i)
                  dq(i,7)  = rhol(i)*twl(i) - rhor(i)*twr(i)
               end if
            end do
            call roe_matrix_p(index,i_end, nl,
     $        delta, sqrtinv, sqrtrhol, sqrtrhor,
     $        ul, ur, vl, vr, wl, wr, tenthalpyl, tenthalpyr,
     $        lx, ly, lz, lxhat, lyhat, lzhat, mx, my, mz, nx, ny, nz,
     $        mxhat, myhat, mzhat, nxhat, nyhat, nzhat, tkl, tkr,
     $        twl, twr, a, ke, qtl, qtr,ilower,
     $        rhol,rhor,tl,tr,pl,pr,vist,volumeinv,
     $        diver, control)
         else
            do   i = i_start,i_end
               dq(i,1)  = rhol(i)     - rhor(i)
               dq(i,2)  = rhoul(i)    - rhour(i)
               dq(i,3)  = rhovl(i)    - rhovr(i)
               dq(i,4)  = rhowl(i)    - rhowr(i)
               dq(i,5)  = rhoel(i)    - rhoer(i)
               if (index.eq.2) then
               end if
               if(nl.eq.6) then
                  dq(i,6)  = rhol(i)*tkl(i) - rhor(i)*tkr(i)
               else if(nl.eq.7) then
                  dq(i,6)  = rhol(i)*tkl(i) - rhor(i)*tkr(i)
                  dq(i,7)  = rhol(i)*twl(i) - rhor(i)*twr(i)
               end if
            end do
            call roe_matrix(index,i_end, nl, dim2,
     $        delta, sqrtinv, sqrtrhol, sqrtrhor,
     $        ul, ur, vl, vr, wl, wr, tenthalpyl, tenthalpyr,
     $        lx, ly, lz, lxhat, lyhat, lzhat, mx, my, mz, nx, ny, nz,
     $        mxhat, myhat, mzhat, nxhat, nyhat, nzhat, tkl, tkr,
     $        twl, twr, a, ke, qtl, qtr, ilower, qt1d, control)
         end if
c         write(*,*)'invsid_flux',a(5,5,i_end-1),a(5,5,i_end)
         do i = i_start,i_end
            fl(i,:) = fl(i,:)+matmul(a(:,:,i),dq(i,:))
         end do

         fl = .5d0 * fl

c 
c     
c..   ZHA2 SCHEME
c     
      CASE (2)


c     modify the dissipation for energy equation,
C     very accurate and smooth for both temperature and velocity

         do i = i_start,i_end

c..   compute Left and Right sound speed

            area_ = dsqrt(ax(i)**2 + ay(i)**2 + az(i)**2)

            c_l = dsqrt(gamma*pl(i)/rhol(i))*area_
            c_r = dsqrt(gamma*pr(i)/rhor(i))*area_

            c_mid = 0.5*(c_l + c_r)

            qml0 = capu_l(i)/c_mid
            qmr0 = capu_r(i)/c_mid

            qml1 = qml0           
            qmr1 = qmr0

            if (qml1.le.1.0.and.qml1.ge.-1.0) then

c..   note: the coefficient 0.5 is ommited for here and muptiplied later

c     plp=pl(i)*(1.d0+qml0)
c     prn=pr(i)*(1.d0-qmr0)
               htl = (rhoel(i) + pl(i))/rhol(i)
               htr = (rhoer(i) + pr(i))/rhor(i)
               
               plp=2.0*pl(i)*(.25*(qml0+1.0)**2*(2.0-qml0)          
     >              + afr*qml0*(qml0**2 - 1.0)**2)
               prn=2.0*pr(i)*(.25*(qmr0-1.0)**2*(2.0+qmr0)
     >              - afr*qmr0*(qmr0**2 - 1.0)**2)

               alfa_l = 2.*(pl(i)/rhol(i))/(pl(i)/rhol(i)+pr(i)/rhor(i))
               alfa_r = 2.*(pr(i)/rhor(i))/(pl(i)/rhol(i)+pr(i)/rhor(i))

               uu_l=c_mid*(0.5*(qml0+abs(qml0)) +alfa_l*(0.25*(qml0+1.)
     >              **2-0.5*(qml0+abs(qml0))))
               uu_r=c_mid*(0.5*(qmr0-abs(qmr0)) +alfa_r*(-0.25*(qmr0-1.)
     >              **2-0.5*(qmr0-abs(qmr0))))

               mid_rhou = rhol(i)*uu_l + rhor(i)*uu_r

               fl(i,1)  = 2.*mid_rhou
               fl(i,2)  = mid_rhou*(ul(i) + ur(i)) 
     >              + ax(i)*plp + ax(i)*prn
               fl(i,3)  = mid_rhou*(vl(i) + vr(i)) 
     >              + ay(i)*plp + ay(i)*prn
               fl(i,4)  = mid_rhou*(wl(i) + wr(i)) 
     >              + az(i)*plp + az(i)*prn

               alfa_l =
     $              2.d0*(tenthalpyl(i)/rhol(i))/(tenthalpyl(i)/rhol(i) !very accurate and smooth
     >              +tenthalpyr(i)/rhor(i))
               alfa_r =
     $              2.d0*(tenthalpyr(i)/rhor(i))/(tenthalpyl(i)/rhol(i)
     >              +tenthalpyr(i)/rhor(i))

               uu_l=c_mid*(0.5*(qml0+abs(qml0)) +alfa_l*(0.25*(qml0+1.)
     >              **2-0.5*(qml0+abs(qml0))))
               uu_r=c_mid*(0.5*(qmr0-abs(qmr0)) +alfa_r*(-0.25*(qmr0-1.)
     >              **2-0.5*(qmr0-abs(qmr0))))

               mid_rhou = rhol(i)*uu_l + rhor(i)*uu_r

               fl(i,5)  = mid_rhou*(el(i) + er(i))
     >              + pl(i)*capu_l(i) + pr(i)*capu_r(i)

               if ( moving.ge.1 ) then
                  fl(i,5) = fl(i,5) - pl(i)*qtl(i) - pr(i)*qtr(i)
               end if
cc     
c     contribution from numerical diffusion part for convective
c     flux(pressure dissipation is already added above in plp and prn)
c     
               fl(i,2) = fl(i,2) -abs(mid_rhou)*(ur(i) - ul(i))
               fl(i,3) = fl(i,3) -abs(mid_rhou)*(vr(i) - vl(i))
               fl(i,4) = fl(i,4) -abs(mid_rhou)*(wr(i) - wl(i))
               fl(i,5) = fl(i,5) -abs(mid_rhou)*(er(i) - el(i))
     >              -(pr(i)- pl(i))*c_mid 

            else if(qml1.gt.1.0) then

               fl(i,1)  = 2.0*rhol(i)*capu_l(i) 
               fl(i,2)  = 2.0*(rhoul(i)*capu_l(i) + ax(i)*pl(i) )
               fl(i,3)  = 2.0*(rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  = 2.0*(rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  = 2.0* (rhoel(i) + pl(i))*capu_l(i) 
               if ( moving.ge.1 ) then
                  fl(i,5) = fl(i,5) - 2.d0*qtl(i)*pl(i)
               end if

            else if(qml1.lt.-1.0) then

               fl(i,1)  = 2.0*rhor(i)*capu_r(i) 
               fl(i,2)  = 2.0*(rhour(i)*capu_r(i) + ax(i)*pr(i) )
               fl(i,3)  = 2.0*(rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  = 2.0*(rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  = 2.0* (rhoer(i) + pr(i))*capu_r(i) 
               if ( moving.ge.1 ) then
                  fl(i,5) = fl(i,5) - 2.d0*qtr(i)*pr(i)
               end if

            end if

         end do

         fl = .5d0*fl

c     
c..   Zha ECUSPLD Scheme, best, 10/15/06 (called Zha12 in the 2004 version code)
c     
      case (3)

         do i = i_start,i_end

            AREA_ = SQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)

            C_L = SQRT(GAMMA*PL(I)/RHOL(I))*AREA_
            C_R = SQRT(GAMMA*PR(I)/RHOR(I))*AREA_

            C_MID = 0.5*(C_L + C_R)

            QML = capu_l(i)/C_MID
            QMR = capu_r(i)/C_MID

                 if(QML.ge.0.)then
                    sign_l=1.0
                 else
                    sign_l=-1.0
                 end if

                 if(QMR.ge.0.)then 
                    sign_r=1.0
                 else
                    sign_r=-1.0
                 end if

                 AFA_LP = 0.5*( 1.0 + sign_l ) 
                 AFA_RN = 0.5*( 1.0 - sign_r ) 

                 BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
                 BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

                 if(0.5*(QML + QMR).ge.0.)then 
                    sign_l=1.0
                 else
                    sign_l=-1.0
                 end if

                 DETAP = 0.5*(1.0 + sign_l )
                 DETAN = 0.5*(1.0 - sign_l )


                 QMLP = 1.0/4.0*(QML+1.0)**2
                 QMLN =-1.0/4.0*(QML-1.0)**2

                 QMRP = 1.0/4.0*(QMR+1.0)**2
                 QMRN =-1.0/4.0*(QMR-1.0)**2

                 QM_MID = BTA_L*DETAP*QMLN - BTA_R*DETAN*QMRP

                 PHI = (rhor(i)*C_R**2)/(rhol(i)*C_L**2)

                 QM_MID_P = QM_MID*(C_R + C_L*PHI)/(C_R +C_L)
                 QM_MID_N = QM_MID*(C_L + C_R/PHI)/(C_R +C_L)

                 C_P = AFA_LP*(1.0 +BTA_L)*QML-BTA_L*QMLP-QM_MID_P
                 C_N = AFA_RN*(1.0 +BTA_R)*QMR-BTA_R*QMRN+QM_MID_N

                 PLP=.250*(QML+1.0)**2*(2.0-QML)
                 PRN=.250*(QMR-1.0)**2*(2.0+QMR)

                 DLP = AFA_LP*(1. + BTA_L) - BTA_L*PLP
                 DRN = AFA_RN*(1. + BTA_R) - BTA_R*PRN



               fl(i,1)  = C_MID*(rhol(i)*C_P + rhor(i)*C_N)

               fl(i,2)  =C_MID*(rhoul(i)*C_P + rhour(i)*C_N)
     >              + ax(i)*(DLP*pl(i) + DRN*pr(i))
               fl(i,3)  =C_MID*(rhovl(i)*C_P + rhovr(i)*C_N)
     >              + ay(i)*(DLP*pl(i) + DRN*pr(i))
               fl(i,4)  =C_MID*(rhowl(i)*C_P + rhowr(i)*C_N)
     >              + az(i)*(DLP*pl(i) + DRN*pr(i))
               fl(i,5)  = C_MID*(rhol(i)*el(i)*C_P + rhor(i)*er(i)*C_N)
               if (nl.eq.6) then
                 fl(i,6) = C_MID*(rhotkl(i)*C_P + rhotkr(i)*C_N)
               end if

c..       Consider moving grid for the pressure term in energy eq. 
c..       NOTE: For stationary grid, qtl(i)=qtr(i)=0, the following naturally satisfy. Zha 10/19/06 
 
                 C_MID = 0.5*(C_L + C_R - qtl(i) - qtr(i))

                 QML_BAR = (capu_l(i)-qtl(i))/C_MID
                 QMR_BAR = (capu_r(i)-qtr(i))/C_MID

                 QMLP = 1.0/4.0*(QML_BAR+1.0)**2
                 QMRN =-1.0/4.0*(QMR_BAR-1.0)**2

                 SL_P = AFA_LP*(1.0 +BTA_L)*QML-BTA_L*QMLP
                 SR_N = AFA_RN*(1.0 +BTA_R)*QMR-BTA_R*QMRN               

                 fl(i,5)  = fl(i,5)  + C_MID*(SL_P*pl(i) + SR_N*pr(i)) 

         end do


c     
c..   Van Leer scheme
c     
      CASE (4)

         do i = i_start,i_end

c..   compute Left and Right sound speed

            AREA_ = DSQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)
            VX1  = AX(I)/AREA_
            VY1  = AY(I)/AREA_
            VZ1  = AZ(I)/AREA_

            C_L = DSQRT(GAMMA*PL(I)/RHOL(I))
            D_L  = rhol(i)
            U_L  = rhoul(i)/rhol(i)
            V_L  = rhovl(i)/rhol(i)
            W_L  = rhowl(i)/rhol(i)
            QML1 = (U_L*VX1+V_L*VY1+W_L*VZ1)/C_L

            C_R = DSQRT(GAMMA*PR(I)/RHOR(I))
            D_R  = rhor(i)
            U_R  = rhour(i)/rhor(i)
            V_R  = rhovr(i)/rhor(i)
            W_R  = rhowr(i)/rhor(i)
            QMR1 = (U_R*VX1+V_R*VY1+W_R*VZ1)/C_R

            IF (QML1.LE.1.0.AND.QML1.GE.-1.0) then

               if ( moving.ge.1 ) then
                  QT_L = qtl(i)/AREA_
                  QT_R = qtr(i)/AREA_
               else
                  QT_L = 0.d0
                  QT_R = 0.d0
               end if

               call van_leer_flux(GAMMA,NL,VX1,VY1,VZ1,C_L,D_L,
     >              U_L,V_L,W_L,C_R,D_R,U_R,V_R,W_R,qt_l, qt_r, UPO)
               do iflux = 1,nl
                  fl(i,iflux)  = AREA_*UPO(iflux)
               end do

            ELSE IF(QML1.GT.1.0) then

               fl(i,1)  =  rhol(i)*capu_l(i) 
               fl(i,2)  =  (rhoul(i)*capu_l(i) + ax(i)*pl(I) )
               fl(i,3)  =  (rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  =  (rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  =   (rhoel(i) + pl(i))*capu_l(i) 

               if ( moving.ge.1 ) then
                  fl(i,5) = fl(i,5) - qtl(i)*pl(i)
               end if

            ELSE IF(QML1.LT.-1.0) then
               
               fl(i,1)  =  rhor(i)*capu_r(i) 
               fl(i,2)  =  (rhour(i)*capu_r(i) + ax(i)*pr(I) )
               fl(i,3)  =  (rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  =  (rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  =  (rhoer(i) + pr(i))*capu_r(i) 
               
               if ( moving.ge.1 ) then
                  fl(i,5) = fl(i,5) - qtr(i)*pr(i)
               end if

            END IF

         end do

c     
c..   EDWARS LDFSS2 Scheme
c     
      case (5)

         do i = i_start,i_end

            AREA_ = SQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)

            htl = (rhoel(i) + pl(i))/rhol(i)
            htr = (rhoer(i) + pr(i))/rhor(i)

            C_L = SQRT(GAMMA*PL(I)/RHOL(I))*AREA_
            C_R = SQRT(GAMMA*PR(I)/RHOR(I))*AREA_

            C_MID = 0.5*(C_L + C_R)

            QML = capu_l(i)/C_MID
            QMR = capu_r(i)/C_MID


                 if(QML.ge.0.)then
                    sign_l=1.0
                 else
                    sign_l=-1.0
                 end if


                 if(QMR.ge.0.)then 
                    sign_r=1.0
                 else
                    sign_r=-1.0
                 end if

                 AFA_LP = 0.5*( 1.0 + sign_l ) 
                 AFA_RN = 0.5*( 1.0 - sign_r ) 

                 BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
                 BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

                 if(0.5*(QML + QMR).ge.0.)then 
                    sign_l=1.0
                 else
                    sign_l=-1.0
                 end if


                 DETAP = 0.5*(1.0 + sign_l )
                 DETAN = 0.5*(1.0 - sign_l )


                 QMLP = 1.0/4.0*(QML+1.0)**2
                 QMLN =-1.0/4.0*(QML-1.0)**2

                 QMRP = 1.0/4.0*(QMR+1.0)**2
                 QMRN =-1.0/4.0*(QMR-1.0)**2

                QM_MID = BTA_L*DETAP*QMLN - BTA_R*DETAN*QMRP

                 PHI = (rhor(i)*C_R**2)/(rhol(i)*C_L**2)

                 QM_MID_P = QM_MID*(C_R + C_L*PHI)/(C_R +C_L)
                 QM_MID_N = QM_MID*(C_L + C_R/PHI)/(C_R +C_L)

                 C_P = AFA_LP*(1.0 +BTA_L)*QML-BTA_L*QMLP-QM_MID_P
                 C_N = AFA_RN*(1.0 +BTA_R)*QMR-BTA_R*QMRN+QM_MID_N

                 PLP=.250*(QML+1.0)**2*(2.0-QML)
                 PRN=.250*(QMR-1.0)**2*(2.0+QMR)

                 DLP = AFA_LP*(1. + BTA_L) - BTA_L*PLP
                 DRN = AFA_RN*(1. + BTA_R) - BTA_R*PRN


               fl(i,1)  = C_MID*(rhol(i)*C_P + rhor(i)*C_N)

               fl(i,2)  =C_MID*(rhoul(i)*C_P + rhour(i)*C_N)
     >              + ax(i)*(DLP*pl(i) + DRN*pr(i))
               fl(i,3)  =C_MID*(rhovl(i)*C_P + rhovr(i)*C_N)
     >              + ay(i)*(DLP*pl(i) + DRN*pr(i))
               fl(i,4)  =C_MID*(rhowl(i)*C_P + rhowr(i)*C_N)
     >              + az(i)*(DLP*pl(i) + DRN*pr(i))
               fl(i,5)  = C_MID*(rhol(i)*htl*C_P + rhor(i)*htr*C_N)

         end do

c     
c..   ZHA6 SCHEME
c     
      CASE (6)

c     add dissipation to remove pressure oscillation, a successful scheme

         do i = i_start,i_end

c..   compute Left and Right sound speed

            AREA_ = DSQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)

            C_L = DSQRT(GAMMA*PL(I)/RHOL(I))*AREA_
            C_R = DSQRT(GAMMA*PR(I)/RHOR(I))*AREA_


            C_MID = 0.5*(C_L + C_R)

            QML0 = capu_l(i)/C_MID
            QMR0 = capu_r(i)/C_MID


            QML1 = QML0           
            QMR1 = QMR0

            IF (QML1.LE.1.0.AND.QML1.GE.-1.0) then

C..   Note: the coefficient 0.5 is ommited for here and muptiplied later

c     PLP=pl(i)*(1.D0+QML0)
c     PRN=pr(i)*(1.D0-QMR0)
               htl = (rhoel(i) + pl(i))/rhol(i)
               htr = (rhoer(i) + pr(i))/rhor(i)
               
               PLP=2.0*pl(i)*(.25*(QML0+1.0)**2*(2.0-QML0)          
     >              + AFR*QML0*(QML0**2 - 1.0)**2)
               PRN=2.0*pr(i)*(.25*(QMR0-1.0)**2*(2.0+QMR0)
     >              - AFR*QMR0*(QMR0**2 - 1.0)**2)

               ALFA_L = 2.*(PL(I)/RHOL(I))/(PL(I)/RHOL(I)+PR(I)/RHOR(I))
               ALFA_R = 2.*(PR(I)/RHOR(I))/(PL(I)/RHOL(I)+PR(I)/RHOR(I))

               uu_L=C_MID*(0.5*(QML0+abs(QML0)) +ALFA_L*(0.25*(QML0+1.)
     >              **2-0.5*(QML0+ABS(QML0))))
               uu_R=C_MID*(0.5*(QMR0-abs(QMR0)) +ALFA_R*(-0.25*(QMR0-1.)
     >              **2-0.5*(QMR0-ABS(QMR0))))

               mid_rhou = rhol(i)*uu_L + rhor(i)*uu_R

               fl(i,1)  = 2.*mid_rhou
               fl(i,2)  = mid_rhou*(ul(i) + ur(i)) 
     >              + ax(i)*PLP + ax(i)*PRN
               fl(i,3)  = mid_rhou*(vl(i) + vr(i)) 
     >              + ay(i)*PLP + ay(i)*PRN
               fl(i,4)  = mid_rhou*(wl(i) + wr(i)) 
     >              + az(i)*PLP + az(i)*PRN
               fl(i,5)  = mid_rhou*(el(i) + er(i))
     >              + pl(i)*capu_l(i) + pr(i)*capu_r(i)
               if ( moving.ge.1 ) then
                  mid_qt = 0.5d0*( qtl(i) + qtr(i) )
                  fl(i,5) = fl(i,5) - pl(i)*qtl(i) - pr(i)*qtr(i) +
     >                 ( pr(i) - pl(i) )*mid_qt
               end if
c     
c     contribution from NUMERICAL DIFFUSION PART FOR CONVECTIVE
C     FLUX(PRESSURE DISSIPATION IS ALREADY ADDED ABOVE IN PLP AND PRN)
c     

               fl(i,2) = fl(i,2) -abs(mid_rhou)*(ur(i) - ul(i))
               fl(i,3) = fl(i,3) -abs(mid_rhou)*(vr(i) - vl(i))
               fl(i,4) = fl(i,4) -abs(mid_rhou)*(wr(i) - wl(i))
               fl(i,5) = fl(i,5) -abs(mid_rhou)*(er(i) - el(i))
     >              -(PR(I)- PL(I))*C_MID 


            ELSE IF(QML1.GT.1.0) then

               fl(i,1)  = 2.0*rhol(i)*capu_l(i) 
               fl(i,2)  = 2.0*(rhoul(i)*capu_l(i) + ax(i)*pl(I) )
               fl(i,3)  = 2.0*(rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  = 2.0*(rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  = 2.0* (rhoel(i) + pl(i))*capu_l(i) 
               if ( moving.ge.1 ) then
                  fl(i,5) = fl(i,5) - 2.d0*qtl(i)*pl(i)
               end if

            ELSE IF(QML1.LT.-1.0) then

               fl(i,1)  = 2.0*rhor(i)*capu_r(i) 
               fl(i,2)  = 2.0*(rhour(i)*capu_r(i) + ax(i)*pr(I) )
               fl(i,3)  = 2.0*(rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  = 2.0*(rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  = 2.0* (rhoer(i) + pr(i))*capu_r(i) 
               if ( moving.ge.1 ) then
                  fl(i,5) = fl(i,5) - 2.d0*qtr(i)*pr(i)
               end if

            END IF

         end do

         fl = .5d0 * fl
  
c     
c..   AUSM+ SCHEME
c     
      case (7)

         do i = i_start,i_end

c..   compute Left and Right sound speed

            AREA_ = DSQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)

            htl = (rhoel(i) + pl(i))/rhol(i)
            htr = (rhoer(i) + pr(i))/rhor(i)

            C_L = DSQRT(2.*(GAMMA-1)/(GAMMA+1)*htl)*AREA_
            C_R = DSQRT(2.*(GAMMA-1)/(GAMMA+1)*htr)*AREA_

            C_L_TILDE = C_L**2/MAX(C_L,ABS(capu_l(i)))
            C_R_TILDE = C_R**2/MAX(C_R,ABS(capu_r(i)))

            C_INTER_FACE = MIN(C_L_TILDE,C_R_TILDE)

            QML1 = capu_l(i)/C_INTER_FACE
            QMR1 = capu_r(i)/C_INTER_FACE

            IF (QML1.LE.1.0.AND.QML1.GE.-1.0) then

               M_INTERFACE = 0.25*(QML1+1.0)**2 + 
     >              BTA*(QML1**2 - 1.0)**2
     >              - 0.25*(QMR1-1.0)**2 - 
     >              BTA*(QMR1**2 - 1.0)**2

               PLP=pl(i)*(.25*(QML1+1.0)**2*(2.0-QML1)
     >              + AFR*QML1*(QML1**2 - 1.0)**2)
               PRN=pr(i)*(.25*(QMR1-1.0)**2*(2.0+QMR1)
     >              - AFR*QMR1*(QMR1**2 - 1.0)**2)

               fl(i,1)  =
     $              C_INTER_FACE*(0.5*M_INTERFACE*(rhol(i)+rhor(i))
     >              -0.5*ABS(M_INTERFACE)*(rhor(i)-rhol(i)))   
               fl(i,2)  =
     $              C_INTER_FACE*(0.5*M_INTERFACE*(rhoul(i)+rhour(i))
     >              -0.5*ABS(M_INTERFACE)*(rhour(i)-rhoul(i)))   
     >              + ax(i)*(PLP + PRN)
               fl(i,3)  =
     $              C_INTER_FACE*(0.5*M_INTERFACE*(rhovl(i)+rhovr(i))
     >              -0.5*ABS(M_INTERFACE)*(rhovr(i)-rhovl(i)))   
     >              + ay(i)*(PLP + PRN)
               fl(i,4)  =
     $              C_INTER_FACE*(0.5*M_INTERFACE*(rhowl(i)+rhowr(i))
     >              -0.5*ABS(M_INTERFACE)*(rhowr(i)-rhowl(i)))   
     >              + az(i)*(PLP + PRN)
               fl(i,5)  = C_INTER_FACE*(0.5*M_INTERFACE*
     >              (rhol(i)*htl+rhor(i)*htr)
     >              -0.5*ABS(M_INTERFACE)*(rhor(i)*htr-rhol(i)*htl))   

            ELSE IF(QML1.GT.1.0) then

               fl(i,1)  =  rhol(i)*capu_l(i) 
               fl(i,2)  =  (rhoul(i)*capu_l(i) + ax(i)*pl(I) )
               fl(i,3)  =  (rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  =  (rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  =   (rhoel(i) + pl(i))*capu_l(i) 

            ELSE IF(QML1.LT.-1.0) then

               fl(i,1)  =  rhor(i)*capu_r(i) 
               fl(i,2)  =  (rhour(i)*capu_r(i) + ax(i)*pr(I) )
               fl(i,3)  =  (rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  =  (rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  =  (rhoer(i) + pr(i))*capu_r(i) 

            END IF

         end do



c     
c..   Zha scheme
c     
      CASE (8)

         do i = i_start,i_end

c..   compute Left and Right sound speed

            AREA_ = DSQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)

            C_L = DSQRT(GAMMA*PL(I)/RHOL(I))*AREA_
            C_R = DSQRT(GAMMA*PR(I)/RHOR(I))*AREA_


            C_MID = 0.5*(C_L + C_R)

            QML0 = capu_l(i)/C_MID
            QMR0 = capu_r(i)/C_MID


            QML1 = QML0           
            QMR1 = QMR0

            IF (QML1.LE.1.0.AND.QML1.GE.-1.0) then

C..   Note: the coefficient 0.5 is ommited for here and muptiplied later

c     PLP=pl(i)*(1.D0+QML0)
c     PRN=pr(i)*(1.D0-QMR0)

               
               PLP=2.0*pl(i)*(.25*(QML0+1.0)**2*(2.0-QML0)          
     >              + AFR*QML0*(QML0**2 - 1.0)**2)
               PRN=2.0*pr(i)*(.25*(QMR0-1.0)**2*(2.0+QMR0)
     >              - AFR*QMR0*(QMR0**2 - 1.0)**2)


               fl(i,1)  = rhol(i)*capu_l(i) + rhor(i)*capu_r(i)
               fl(i,2)  = rhoul(i)*capu_l(i) + ax(i)*PLP +
     >              rhour(i)*capu_r(i) + ax(i)*PRN
               fl(i,3)  = rhovl(i)*capu_l(i) + ay(i)*PLP +
     >              rhovr(i)*capu_r(i) + ay(i)*PRN
               fl(i,4)  = rhowl(i)*capu_l(i) + az(i)*PLP +
     >              rhowr(i)*capu_r(i) + az(i)*PRN
               fl(i,5)  = (rhoel(i) + pl(i))*capu_l(i) +
     >              (rhoer(i) + pr(i))*capu_r(i)

c     
c     contribution from NUMERICAL DIFFUSION PART FOR CONVECTIVE
C     FLUX(PRESSURE DISSIPATION IS ALREADY ADDED ABOVE IN PLP AND PRN)
c     



               fl(i,1) = fl(i,1) -( DABS(CAPU_R(I)) * RHOR(I)
     >              -DABS(CAPU_L(I)) * RHOL(I) ) 
               fl(i,2) = fl(i,2) -( DABS(CAPU_R(I)) * RHOUR(I)
     >              -DABS(CAPU_L(I)) * RHOUL(I) )
               fl(i,3) = fl(i,3) -( DABS(CAPU_R(I)) * RHOVR(I)
     >              -DABS(CAPU_L(I)) * RHOVL(I) )
               fl(i,4) = fl(i,4)-(  DABS(CAPU_R(I)) * RHOWR(I)
     >              -DABS(CAPU_L(I)) * RHOWL(I) )
               fl(i,5) = fl(i,5)-(  DABS(CAPU_R(I)) * RHOER(I)
     >              -DABS(CAPU_L(I)) * RHOEL(I) )
     >              -( PR(I) - PL(I))*C_MID 


            ELSE IF(QML1.GT.1.0) then

               fl(i,1)  = 2.0*rhol(i)*capu_l(i) 
               fl(i,2)  = 2.0*(rhoul(i)*capu_l(i) + ax(i)*pl(I) )
               fl(i,3)  = 2.0*(rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  = 2.0*(rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  = 2.0* (rhoel(i) + pl(i))*capu_l(i) 


            ELSE IF(QML1.LT.-1.0) then

               fl(i,1)  = 2.0*rhor(i)*capu_r(i) 
               fl(i,2)  = 2.0*(rhour(i)*capu_r(i) + ax(i)*pr(I) )
               fl(i,3)  = 2.0*(rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  = 2.0*(rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  = 2.0* (rhoer(i) + pr(i))*capu_r(i) 


            END IF

         end do

         fl = .5d0*fl
c     
c..   AUSMV
c     
      CASE (9)

         do i = i_start,i_end

c..   compute Left and Right sound speed

            AREA_ = DSQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)

            C_L = DSQRT(GAMMA*PL(I)/RHOL(I))*AREA_
            C_R = DSQRT(GAMMA*PR(I)/RHOR(I))*AREA_


            C_MID = 0.5*(C_L + C_R)

            QML0 = capu_l(i)/C_MID
            QMR0 = capu_r(i)/C_MID


            QML1 = QML0           
            QMR1 = QMR0
            htl = (rhoel(i) + pl(i))/rhol(i)
            htr = (rhoer(i) + pr(i))/rhor(i)

            IF (QML1.LE.1.0.AND.QML1.GE.-1.0) then

C..   Note: the coefficient 0.5 is ommited for here and muptiplied later

c     PLP=pl(i)*(1.D0+QML0)
c     PRN=pr(i)*(1.D0-QMR0)

               
               PLP=2.0*pl(i)*(.25*(QML0+1.0)**2*(2.0-QML0)          
     >              + AFR*QML0*(QML0**2 - 1.0)**2)
               PRN=2.0*pr(i)*(.25*(QMR0-1.0)**2*(2.0+QMR0)
     >              - AFR*QMR0*(QMR0**2 - 1.0)**2)


               fl(i,1)  = rhol(i)*capu_l(i) + rhor(i)*capu_r(i)
               fl(i,2)  = rhoul(i)*capu_l(i) + ax(i)*PLP +
     >              rhour(i)*capu_r(i) + ax(i)*PRN
               fl(i,3)  = rhovl(i)*capu_l(i) + ay(i)*PLP +
     >              rhovr(i)*capu_r(i) + ay(i)*PRN
               fl(i,4)  = rhowl(i)*capu_l(i) + az(i)*PLP +
     >              rhowr(i)*capu_r(i) + az(i)*PRN
               fl(i,5)  = (rhoel(i) + pl(i))*capu_l(i) +
     >              (rhoer(i) + pr(i))*capu_r(i)

c     
c     contribution from NUMERICAL DIFFUSION PART FOR CONVECTIVE
C     FLUX(PRESSURE DISSIPATION IS ALREADY ADDED ABOVE IN PLP AND PRN)
c     

               fl(i,1) = fl(i,1) -( DABS(CAPU_R(I)) * RHOR(I)
     >              -DABS(CAPU_L(I)) * RHOL(I) ) 
               fl(i,2) = fl(i,2) -( DABS(CAPU_R(I)) * RHOUR(I)
     >              -DABS(CAPU_L(I)) * RHOUL(I) )
               fl(i,3) = fl(i,3) -( DABS(CAPU_R(I)) * RHOVR(I)
     >              -DABS(CAPU_L(I)) * RHOVL(I) )
               fl(i,4) = fl(i,4)-(  DABS(CAPU_R(I)) * RHOWR(I)
     >              -DABS(CAPU_L(I)) * RHOWL(I) )
               fl(i,5) = fl(i,5)-(  DABS(CAPU_R(I)) * RHOR(I)*htr
     >              -DABS(CAPU_L(I)) * RHOL(I)*htl )


c..   add dissipation to correct pressure oscillation 


               ALFA_L = 2.*(PL(I)/RHOL(I))/(PL(I)/RHOL(I)+PR(I)/RHOR(I))
               ALFA_R = 2.*(PR(I)/RHOR(I))/(PL(I)/RHOL(I)+PR(I)/RHOR(I))

               OMEGA_L=C_MID*ALFA_L*
     $              (0.25*(QML0+1.)**2-0.5*(QML0+ABS(QML0)))
               OMEGA_R=C_MID*ALFA_R*
     $              (0.25*(QMR0-1.)**2+0.5*(QMR0-ABS(QMR0)))

               fl(i,1) =
     $              fl(i,1) -2.0*(OMEGA_R*RHOR(I) - OMEGA_L*RHOL(I))
               fl(i,2) =
     $              fl(i,2) -2.0*(OMEGA_R*RHOUR(I) - OMEGA_L*RHOUL(I))
               fl(i,3) =
     $              fl(i,3) -2.0*(OMEGA_R*RHOVR(I) - OMEGA_L*RHOVL(I))
               fl(i,4) =
     $              fl(i,4) -2.0*(OMEGA_R*RHOWR(I) - OMEGA_L*RHOWL(I))
               fl(i,5) =
     $              fl(i,5) -2.0*(OMEGA_R*RHOR(I)*htr - 
     >              OMEGA_L*RHOL(I)*htl)

            ELSE IF(QML1.GT.1.0) then

               fl(i,1)  = 2.0*rhol(i)*capu_l(i) 
               fl(i,2)  = 2.0*(rhoul(i)*capu_l(i) + ax(i)*pl(I) )
               fl(i,3)  = 2.0*(rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  = 2.0*(rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  = 2.0* (rhoel(i) + pl(i))*capu_l(i) 


            ELSE IF(QML1.LT.-1.0) then

               fl(i,1)  = 2.0*rhor(i)*capu_r(i) 
               fl(i,2)  = 2.0*(rhour(i)*capu_r(i) + ax(i)*pr(I) )
               fl(i,3)  = 2.0*(rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  = 2.0*(rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  = 2.0* (rhoer(i) + pr(i))*capu_r(i) 


            END IF

         end do

         fl = .5d0*fl
c     
c..   AUSMD
c     
      CASE (10)

         do i = i_start,i_end

c..   compute Left and Right sound speed

            AREA_ = DSQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)

            C_L = DSQRT(GAMMA*PL(I)/RHOL(I))*AREA_
            C_R = DSQRT(GAMMA*PR(I)/RHOR(I))*AREA_


            C_MID = 0.5*(C_L + C_R)

            QML0 = capu_l(i)/C_MID
            QMR0 = capu_r(i)/C_MID


            QML1 = QML0           
            QMR1 = QMR0
            htl = (rhoel(i) + pl(i))/rhol(i)
            htr = (rhoer(i) + pr(i))/rhor(i)

            IF (QML1.LE.1.0.AND.QML1.GE.-1.0) then

C..   Note: the coefficient 0.5 is ommited for here and muptiplied later

c     PLP=pl(i)*(1.D0+QML0)
c     PRN=pr(i)*(1.D0-QMR0)

               
               PLP=2.0*pl(i)*(.25*(QML0+1.0)**2*(2.0-QML0)          
     >              + AFR*QML0*(QML0**2 - 1.0)**2)
               PRN=2.0*pr(i)*(.25*(QMR0-1.0)**2*(2.0+QMR0)
     >              - AFR*QMR0*(QMR0**2 - 1.0)**2)
               htl = (rhoel(i) + pl(i))/rhol(i)
               htr = (rhoer(i) + pr(i))/rhor(i)

               ALFA_L = 2.*(PL(I)/RHOL(I))/(PL(I)/RHOL(I)+PR(I)/RHOR(I))
               ALFA_R = 2.*(PR(I)/RHOR(I))/(PL(I)/RHOL(I)+PR(I)/RHOR(I))

               uu_L=C_MID*(0.5*(QML0+abs(QML0)) +ALFA_L*(0.25*(QML0+1.)
     >              **2-0.5*(QML0+ABS(QML0))))
               uu_R=C_MID*(0.5*(QMR0-abs(QMR0)) +ALFA_R*(-0.25*(QMR0-1.)
     >              **2-0.5*(QMR0-ABS(QMR0))))

               mid_rhou = rhol(i)*uu_L + rhor(i)*uu_R

               fl(i,1)  = 2.*mid_rhou
               fl(i,2)  = mid_rhou*(ul(i) + ur(i)) 
     >              + ax(i)*PLP + ax(i)*PRN
               fl(i,3)  = mid_rhou*(vl(i) + vr(i)) 
     >              + ay(i)*PLP + ay(i)*PRN
               fl(i,4)  = mid_rhou*(wl(i) + wr(i)) 
     >              + az(i)*PLP + az(i)*PRN
               fl(i,5)  = mid_rhou*(htl + htr)

c     
c     contribution from NUMERICAL DIFFUSION PART FOR CONVECTIVE
C     FLUX(PRESSURE DISSIPATION IS ALREADY ADDED ABOVE IN PLP AND PRN)
c     


               fl(i,2) = fl(i,2) -abs(mid_rhou)*(ur(i) - ul(i))
               fl(i,3) = fl(i,3) -abs(mid_rhou)*(vr(i) - vl(i))
               fl(i,4) = fl(i,4) -abs(mid_rhou)*(wr(i) - wl(i))
               fl(i,5) = fl(i,5) -abs(mid_rhou)*(htr -htl)

            ELSE IF(QML1.GT.1.0) then

               fl(i,1)  = 2.0*rhol(i)*capu_l(i) 
               fl(i,2)  = 2.0*(rhoul(i)*capu_l(i) + ax(i)*pl(I) )
               fl(i,3)  = 2.0*(rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  = 2.0*(rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  = 2.0* (rhoel(i) + pl(i))*capu_l(i) 


            ELSE IF(QML1.LT.-1.0) then

               fl(i,1)  = 2.0*rhor(i)*capu_r(i) 
               fl(i,2)  = 2.0*(rhour(i)*capu_r(i) + ax(i)*pr(I) )
               fl(i,3)  = 2.0*(rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  = 2.0*(rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  = 2.0* (rhoer(i) + pr(i))*capu_r(i) 


            END IF

         end do

         fl = .5d0*fl

c     
c..   VAN LEER-HANEL SCHEME
c     
      CASE (11)

         do i = i_start,i_end

c..   compute Left and Right sound speed

            AREA_ = DSQRT(AX(I)**2 + AY(I)**2 + AZ(I)**2)
            VX1  = AX(I)/AREA_
            VY1  = AY(I)/AREA_
            VZ1  = AZ(I)/AREA_

            C_L = DSQRT(GAMMA*PL(I)/RHOL(I))
            D_L  = rhol(i)
            U_L  = rhoul(i)/rhol(i)
            V_L  = rhovl(i)/rhol(i)
            W_L  = rhowl(i)/rhol(i)
            QML1 = (U_L*VX1+V_L*VY1+W_L*VZ1)/C_L
            htl = (rhoel(i) + pl(i))/rhol(i)

            C_R = DSQRT(GAMMA*PR(I)/RHOR(I))
            D_R  = rhor(i)
            U_R  = rhour(i)/rhor(i)
            V_R  = rhovr(i)/rhor(i)
            W_R  = rhowr(i)/rhor(i)
            QMR1 = (U_R*VX1+V_R*VY1+W_R*VZ1)/C_R
            htr = (rhoer(i) + pr(i))/rhor(i)

            IF (QML1.LE.1.0.AND.QML1.GE.-1.0) then

               call vanleer_hanel_flux(GAMMA,NL,VX1,VY1,VZ1,C_L,D_L,
     >              U_L,V_L,W_L,htl,C_R,D_R,U_R,V_R,W_R,htr,UPO)

               do iflux = 1,nl
                  fl(i,iflux)  = AREA_*UPO(iflux)
               end do

            ELSE IF(QML1.GT.1.0) then

               fl(i,1)  =  rhol(i)*capu_l(i) 
               fl(i,2)  =  (rhoul(i)*capu_l(i) + ax(i)*pl(I) )
               fl(i,3)  =  (rhovl(i)*capu_l(i) + ay(i)*pl(i) )  
               fl(i,4)  =  (rhowl(i)*capu_l(i) + az(i)*pl(i) )
               fl(i,5)  =   (rhoel(i) + pl(i))*capu_l(i) 

            ELSE IF(QML1.LT.-1.0) then

               fl(i,1)  =  rhor(i)*capu_r(i) 
               fl(i,2)  =  (rhour(i)*capu_r(i) + ax(i)*pr(I) )
               fl(i,3)  =  (rhovr(i)*capu_r(i) + ay(i)*pr(i) )  
               fl(i,4)  =  (rhowr(i)*capu_r(i) + az(i)*pr(i) )
               fl(i,5)  =  (rhoer(i) + pr(i))*capu_r(i) 

            END IF
         end do
      END SELECT                ! Scheme OPtion Ends
c     
      if (moving.eq.0) then
        select case (index)
          case(1)
            if(bclower.eq.3.or.bclower.eq.19 .or.
     $              (bclower.ge.101 .and. bclower.le.110)) then
               i = 1
               fl(i,:) = 0.d0
               fl(i,2) = lx(i)*pr(i)
               fl(i,3) = ly(i)*pr(i)
               fl(i,4) = lz(i)*pr(i)
            end if
            if(bcupper.eq.3.or.bcupper.eq.19 .or.
     $              (bcupper.ge.101 .and. bcupper.le.110)) then
               i = il+1
               fl(i,:) = 0.d0
               fl(i,2) = lx(i)*pl(i)
               fl(i,3) = ly(i)*pl(i)
               fl(i,4) = lz(i)*pl(i)
            end if
          case(2)
            if(bclower.eq.3.or.bclower.eq.19 .or.
     $              (bclower.ge.101 .and. bclower.le.110)) then
               i = 1
               fl(i,:) = 0.d0
               fl(i,2) = mx(i)*pr(i)
               fl(i,3) = my(i)*pr(i)
               fl(i,4) = mz(i)*pr(i)
            end if
            if(bcupper.eq.3.or.bcupper.eq.19.or.
     $              (bcupper.ge.101 .and. bcupper.le.110)) then
               i = jl+1
               fl(i,:) = 0.d0
               fl(i,2) = mx(i)*pl(i)
               fl(i,3) = my(i)*pl(i)
               fl(i,4) = mz(i)*pl(i)
            end if
          case(3)
            if(bclower.eq.3.or.bclower.eq.19 .or.
     $              (bclower.ge.101 .and. bclower.le.110)) then
               i = 1
               fl(i,:) = 0.d0
               fl(i,2) = nx(i)*pr(i)
               fl(i,3) = ny(i)*pr(i)
               fl(i,4) = nz(i)*pr(i)
            end if
            if(bcupper.eq.3.or.bcupper.eq.19 .or.
     $              (bcupper.ge.101 .and. bcupper.le.110)) then
               i = kl+1
               fl(i,:) = 0.d0
               fl(i,2) = nx(i)*pl(i)
               fl(i,3) = ny(i)*pl(i)
               fl(i,4) = nz(i)*pl(i)
            end if
          end select
      end if
c
      return
      end
