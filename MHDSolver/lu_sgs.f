C-------------------------------------------------------------
        SUBROUTINE lu_sgs(dt, x, y, z, vol, q, q_n, rflux, sflux, tflux,
     $     bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta, 
     $     rhs, il, jl,
     $     kl, nl, dim1, dim2, ilower, iupper, jlower, jupper, klower,
     $     kupper, visturb,  dstep, checktime, wallfunc, ke,
     $     qt,  sr,  udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $     wdk, cst_nri, ced_nri, cst_nrj, ced_nrj, control)

c     use LU-SGS
c     to invert the implicit LHS
C     parameters are same as the subroutine ursn.f
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
     $     il, jl, kl, nl, dim2, dim1,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     dstep
      integer::
     $     dual_t,
     $     moving, gcl

      double precision, intent(in)::
     $     rhs(il,jl,kl,nl),
     $     sr(il,jl,kl)

      double precision, dimension(jl,kl,4), intent(in)::
     $     udi, vdi, wdi

      double precision, dimension(il,kl,4), intent(in)::
     $     udj, vdj, wdj

      double precision, dimension(il,jl,4), intent(in)::
     $     udk, vdk, wdk

      double precision, intent(IN)::
     $     cst_nri(nl,nl,jl,kl,5), ced_nri(nl,nl,jl,kl,5),
     $     cst_nrj(nl,nl,il,kl,5), ced_nrj(nl,nl,il,kl,5)

      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z

      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol

      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q

      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, 3)::
     $     qt

      double precision, dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1), intent(in)::
     $     visturb

      double precision, intent(in)::
     $     dt(il,jl,kl)
      double precision::
     $     gamma, prandtl, prt, tref, ptotal, ttotal,
     $     angl1, angl2, poutlet

      integer, intent(in)::
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl)
      
      double precision, intent(in)::
     $     delta(5)

      integer::
     $     rhs_order, limiter, rhs_scheme, idimen,
     $     lhs_order, iter_gs, lhs_scheme
      
      double precision::
     $     machinf, epsfactor, kfactor, reynolds,
     $     tintvl,
     $     ronum

      logical, intent(in)::
     $     checktime,           ! if systime check necessary
     $     wallfunc,
     $     ke
      logical::
     $     inviscid, suther

      double precision, dimension(il+1, jl+1, kl+1, nl), intent(in)::
     $     rflux, sflux, tflux

      double precision, intent(out), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q_n

c     
c     LOCAL VARIABLES
c     
      double precision, dimension(nl,nl,dim1)::
     $     capa, capb, capc,
     $     f(nl,dim1),
     $     soln(nl,dim1)

      double precision::
     $     qtline(-1:dim2,3)

      double precision, dimension(nl,nl,3,2)::
     $     ci_nr

c..   coefficient matrices

      double precision, dimension(nl,nl,il,jl,kl)::
     $     a_minus, a_plus,
     $     b_minus, b_plus,
     $     c_minus, c_plus,
     $     b_bar

      integer::
     $     n,
     $     ii, jj,
     $     i, j, k,
     $     nblocks,
     $     is,js,
     $     m,
     $     nj,
     $     bclower, bcupper

      double precision::
     $     cx1,cx2,
     $     cy1,cy2,
     $     cz1,cz2

      integer, parameter::
     $     bcperi = 10,
     $     bcwall = 3
c-----------next is for metric
      double precision, dimension(dim2)::
     $     lxhat, lyhat, lzhat,
     $     mxhat, myhat, mzhat,
     $     nxhat, nyhat, nzhat,
     $     lx, ly, lz,
     $     mx, my, mz,
     $     nx, ny, nz
c     
c-----------end of metric



	integer:: index,imax,ip,jp,kp,l,ir,jr,kr
	double precision:: xx,xy,xz,yx,yy,yz,zx,zy,zz
	double precision:: r0,u0,v0,w0,e0,p0,t,as,qq,qq1,qq2,qq3,
     $                   qqx,qqy,qqz,rr,rr1,rr2,rr3,rho,spec,
     $                   s1,s2,s3,s4,s5,a2,a5,uvw,cs,phi,phi1,vnu
	double precision, dimension(il,jl,kl)::
     $                  rh

	double precision, dimension(il,kl,nl)::
     $                  a,b,c

	double precision, 
     $ dimension(3,4,ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $                  xk

      double precision, dimension(0:il+1,0:jl+1,0:kl+1)::
     $     ZMUL
c--------------------------parameter transfer
     
      dual_t=control%dual_t
      gcl=control%gcl
      idimen=control%idimen
      iter_gs=control%iter_gs
      lhs_order=control%lhs_order
      lhs_scheme=control%lhs_scheme
      limiter=control%limiter
      moving=control%moving
      rhs_order=control%rhs_order
      rhs_scheme=control%rhs_scheme


      inviscid=control%inviscid
      suther=control%suther

      angl1=control%angl1
      angl2=control%angl2
      epsfactor=control%epsfactor
      gamma=control%gamma
      kfactor=control%kfactor
      machinf=control%machinf
      poutlet=control%poutlet
      ptotal=control%ptotal
      prandtl=control%prandtl
      prt=control%prt
      reynolds=control%reynolds
      ronum=control%ronum
      tref=control%tref
      tintvl=control%tintvl
      ttotal=control%ttotal
c--------------------------------------------------------
c
c ... the array q_n is used for dq in this subroutine
c
c *** SUBROUTINE START ***
c     
      q_n = 0.d0
c$$$      a_minus = 0.d0
c$$$      a_plus = 0.d0
c$$$      b_minus = 0.d0
c$$$      b_plus = 0.d0
c$$$      c_minus = 0.d0
c$$$      c_plus = 0.d0

c ... calculate the matrices in three direction

	index=1
	imax = il+1
	do j=1,jl
	do k=1,kl
      call metric(index,j, k, imax, il, jl, kl, x, y, z, lx, ly, lz,
     $     mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat, mxhat, myhat,
     $     mzhat, nxhat, nyhat, nzhat, ilower, iupper, jlower, jupper,
     $     klower, kupper, bclower, bcupper, dim2)
      do i=1,imax
	xk(1,1,i,j,k) = lx(i)	
	xk(1,2,i,j,k) = ly(i)
	xk(1,3,i,j,k) = lz(i)
      end do
	end do
	end do
	
	index=2
	imax = jl+1
	do i=1,il
	do k=1,kl
      call metric(index,i, k, imax, il, jl, kl, x, y, z, lx, ly, lz,
     $     mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat, mxhat, myhat,
     $     mzhat, nxhat, nyhat, nzhat, ilower, iupper, jlower, jupper,
     $     klower, kupper, bclower, bcupper, dim2)
      do j=1,imax
	xk(2,1,i,j,k) = mx(j)	
	xk(2,2,i,j,k) = my(j)
	xk(2,3,i,j,k) = mz(j)
      end do
	end do
	end do

	index= 3
	imax = kl+1
	do i=1,il
	do j=1,jl
      call metric(index,i, j, imax, il, jl, kl, x, y, z, lx, ly, lz,
     $     mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat, mxhat, myhat,
     $     mzhat, nxhat, nyhat, nzhat, ilower, iupper, jlower, jupper,
     $     klower, kupper, bclower, bcupper, dim2)
      do k=1,imax
	xk(3,1,i,j,k) = nx(k)	
	xk(3,2,i,j,k) = ny(k)
	xk(3,3,i,j,k) = nz(k)
      end do
	end do
	end do	
c         end of the matrices
c---------------------------------------------			

        CS=110.6/TREF
        PHI=0.
	  PHI1=1./(1.+PHI)

        DO J=1,JL
        DO I=1,IL
	  DO K=1,KL
        IP=I+1
	  JP=J+1
	  KP=K+1
        R0=Q(I,J,K,1)
        U0=Q(I,J,K,2)/R0
        V0=Q(I,J,K,3)/R0
        W0=Q(I,J,K,4)/R0
        E0=Q(I,J,K,5)/R0
        P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2))
        AS=SQRT(GAMMA*P0/R0)
        XX=0.5*(XK(1,1,I,J,K)+XK(1,1,IP,J,K))
        XY=0.5*(XK(1,2,I,J,K)+XK(1,2,IP,J,K))
        XZ=0.5*(XK(1,3,I,J,K)+XK(1,3,IP,J,K))
C        XT=0.5*(XK(1,4,I,J,K)+XK(1,4,IP,J,K))
        QQ1=XX*U0+XY*V0+XZ*W0
c	  QQX=ABS(XT+QQ1)
	  QQX=ABS(QQ1)
	  RR1=SQRT(XX**2+XY**2+XZ**2)
        YX=0.5*(XK(2,1,I,J,K)+XK(2,1,I,JP,K))
        YY=0.5*(XK(2,2,I,J,K)+XK(2,2,I,JP,K))
        YZ=0.5*(XK(2,3,I,J,K)+XK(2,3,I,JP,K))
C        YT=0.5*(XK(2,4,I,J,K)+XK(2,4,I,JP,K))
        QQ2=YX*U0+YY*V0+YZ*W0
c	  QQY=ABS(YT+QQ2)
	  QQY=ABS(QQ2)
	  RR2=SQRT(YX**2+YY**2+YZ**2)
        ZX=0.5*(XK(3,1,I,J,K)+XK(3,1,I,J,KP))
        ZY=0.5*(XK(3,2,I,J,K)+XK(3,2,I,J,KP))
        ZZ=0.5*(XK(3,3,I,J,K)+XK(3,3,I,J,KP))
C        ZT=0.5*(XK(3,4,I,J,K)+XK(3,4,I,J,KP))
        QQ3=ZX*U0+ZY*V0+ZZ*W0
c	  QQZ=ABS(ZT+QQ3)
	  QQZ=ABS(QQ3)
	  RR3=SQRT(ZX**2+ZY**2+ZZ**2)
        T=GAMMA*MACHINF*MACHINF*P0/R0
        ZMUL(I,J,K)=(1.+TREF)/(T+TREF)*T**1.5
        VNU=2.*(ZMUL(I,J,K)+VISTURB(I,J,K)*REYNOLDS)*
     &      (RR1**2+RR2**2+RR3**2)/(R0*REYNOLDS)
C        VNU=0.
	  RH(I,J,K)=1.0/(1.0+PHI1*DT(I,J,K)/VOL(I,J,K)
     &	  *(QQX+QQY+QQZ+AS*(RR1+RR2+RR3))+VNU)

c	  RH(I,J,K)=1.0/(1.0+PHI1*VOL(I,J,K)/DT(I,J,K)
c     &	  *(QQX+QQY+QQZ+AS*(RR1+RR2+RR3))+VNU)
        END DO
	  END DO
	  END DO
	          
        DO 100 J=1,JL
	  DO 110 I=1,IL
	  DO 120 K=1,KL
C---------------------------------J------------------------	  
        JR=J-1
	  IF(J.EQ.1) THEN
        DO L=1,NL
	  B(I,K,L)=0.0
	  END DO
	  ELSE
        R0=Q(I,JR,K,1)
        U0=Q(I,JR,K,2)/R0
        V0=Q(I,JR,K,3)/R0
        W0=Q(I,JR,K,4)/R0
        E0=Q(I,JR,K,5)/R0
        P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2))
        AS=SQRT(GAMMA*P0/R0)
	  UVW=0.5*(U0*U0+V0*V0+W0*W0)
        XX=0.5*(XK(2,1,I,JR,K)+XK(2,1,I,J,K))
        YY=0.5*(XK(2,2,I,JR,K)+XK(2,2,I,J,K))
        ZZ=0.5*(XK(2,3,I,JR,K)+XK(2,3,I,J,K))
C        TT=0.5*(XK(2,4,I,JR,K)+XK(2,4,I,J,K))
        QQ=XX*U0+YY*V0+ZZ*W0
	  RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
C        RHO=ABS(TT+QQ)+AS*RR
        RHO=ABS(QQ)+AS*RR
        VNU=2.*(ZMUL(I,JR,K)+VISTURB(I,JR,K)*REYNOLDS)*
     &      (RR1**2+RR2**2+RR3**2)/(R0*REYNOLDS)
C        VNU=0.
C	  SPEC=0.5*(TT+QQ+RHO)
	  SPEC=0.5*(QQ+RHO)+0.5*VNU
c	  SPEC=0.5*RHO+0.5*VNU
	  S1=Q_N(I,JR,K,1)
	  S2=Q_N(I,JR,K,2)
	  S3=Q_N(I,JR,K,3)
	  S4=Q_N(I,JR,K,4)
	  S5=Q_N(I,JR,K,5)
	  A5=0.5*(XX*S2+YY*S3+ZZ*S4-QQ*S1)
	  A2=(GAMMA-1.)/2.*(UVW*S1-(U0*S2+V0*S3+W0*S4)+S5)
        B(I,K,1)=A5+SPEC*S1
	  B(I,K,2)=XX*A2+U0*A5+SPEC*S2
	  B(I,K,3)=YY*A2+V0*A5+SPEC*S3
	  B(I,K,4)=ZZ*A2+W0*A5+SPEC*S4
	  B(I,K,5)=QQ*A2+(GAMMA*E0-(GAMMA-1.)*UVW)*A5+SPEC*S5
        END IF
C-------------------------I------------------------------
	  IR=I-1
	  IF(I.EQ.1) THEN
        DO L=1,NL
	  A(I,K,L)=0.
	  END DO
	  ELSE
        R0=Q(IR,J,K,1)
        U0=Q(IR,J,K,2)/R0
        V0=Q(IR,J,K,3)/R0
        W0=Q(IR,J,K,4)/R0
        E0=Q(IR,J,K,5)/R0
        P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2))
        AS=SQRT(GAMMA*P0/R0)
	  UVW=0.5*(U0*U0+V0*V0+W0*W0)
        XX=0.5*(XK(1,1,IR,J,K)+XK(1,1,I,J,K))
        YY=0.5*(XK(1,2,IR,J,K)+XK(1,2,I,J,K))
        ZZ=0.5*(XK(1,3,IR,J,K)+XK(1,3,I,J,K))
C        TT=0.5*(XK(1,4,IR,J,K)+XK(1,4,I,J,K))
        QQ=XX*U0+YY*V0+ZZ*W0
	  RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
C        RHO=ABS(TT+QQ)+AS*RR
        RHO=ABS(QQ)+AS*RR
        VNU=2.*(ZMUL(IR,J,K)+VISTURB(IR,J,K)*REYNOLDS)*
     &      (RR1**2+RR2**2+RR3**2)/(R0*REYNOLDS)
C        VNU=0.
C	  SPEC=0.5*(TT+QQ+RHO)
	  SPEC=0.5*(QQ+RHO)+0.5*VNU
c	  SPEC=0.5*RHO+0.5*VNU
	  S1=Q_N(IR,J,K,1)
	  S2=Q_N(IR,J,K,2)
	  S3=Q_N(IR,J,K,3)
	  S4=Q_N(IR,J,K,4)
	  S5=Q_N(IR,J,K,5)
	  A5=0.5*(XX*S2+YY*S3+ZZ*S4-QQ*S1)
	  A2=(GAMMA-1.)/2.*(UVW*S1-(U0*S2+V0*S3+W0*S4)+S5)
        A(I,K,1)=A5+SPEC*S1
	  A(I,K,2)=XX*A2+U0*A5+SPEC*S2
	  A(I,K,3)=YY*A2+V0*A5+SPEC*S3
	  A(I,K,4)=ZZ*A2+W0*A5+SPEC*S4
	  A(I,K,5)=QQ*A2+(GAMMA*E0-(GAMMA-1.)*UVW)*A5+SPEC*S5
        END IF
C------------------------K--------------------------
	  KR=K-1
        IF(K.EQ.1) THEN
	  DO L=1,NL
        C(I,K,L)=0.
	  END DO
	  ELSE
        R0=Q(I,J,KR,1)
        U0=Q(I,J,KR,2)/R0
        V0=Q(I,J,KR,3)/R0
        W0=Q(I,J,KR,4)/R0
        E0=Q(I,J,KR,5)/R0
        P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2)) 
        AS=SQRT(GAMMA*P0/R0)
	  UVW=0.5*(U0*U0+V0*V0+W0*W0)
        XX=0.5*(XK(3,1,I,J,KR)+XK(3,1,I,J,K))
        YY=0.5*(XK(3,2,I,J,KR)+XK(3,2,I,J,K))
        ZZ=0.5*(XK(3,3,I,J,KR)+XK(3,3,I,J,K))
C        TT=0.5*(XK(3,4,I,J,KR)+XK(3,4,I,J,K))
        QQ=XX*U0+YY*V0+ZZ*W0
	  RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
C        RHO=ABS(TT+QQ)+AS*RR
        RHO=ABS(QQ)+AS*RR
        VNU=2.*(ZMUL(I,J,KR)+VISTURB(I,J,KR)*REYNOLDS)*
     &      (RR1**2+RR2**2+RR3**2)/(R0*REYNOLDS)
C        VNU=0.
C	  SPEC=0.5*(TT+QQ+RHO)
	  SPEC=0.5*(QQ+RHO)+0.5*VNU
c	  SPEC=0.5*RHO+0.5*VNU
	  S1=Q_N(I,J,KR,1)
	  S2=Q_N(I,J,KR,2)
	  S3=Q_N(I,J,KR,3)
	  S4=Q_N(I,J,KR,4)
	  S5=Q_N(I,J,KR,5)
	  A5=0.5*(XX*S2+YY*S3+ZZ*S4-QQ*S1)
	  A2=(GAMMA-1.)/2.*(UVW*S1-(U0*S2+V0*S3+W0*S4)+S5)
        C(I,K,1)=A5+SPEC*S1
	  C(I,K,2)=XX*A2+U0*A5+SPEC*S2
	  C(I,K,3)=YY*A2+V0*A5+SPEC*S3
	  C(I,K,4)=ZZ*A2+W0*A5+SPEC*S4
	  C(I,K,5)=QQ*A2+(GAMMA*E0-(GAMMA-1.)*UVW)*A5+SPEC*S5
        END IF
	  PHI=PHI1*DT(I,J,K)/VOL(I,J,K)

c	  PHI=PHI1*VOL(I,J,K)/DT(I,J,K)

        DO L=1,NL
	  Q_N(I,J,K,L)=(RHS(I,J,K,L)
     $     +PHI*(A(I,K,L)+B(I,K,L)+C(I,K,L)))*RH(I,J,K)
        END DO
120     CONTINUE
110     CONTINUE
100     CONTINUE

        DO 200 J=JL,1,-1
	  DO 210 I=IL,1,-1
	  DO 220 K=KL,1,-1
C---------------------------------J------------------------	  
        JR=J+1
	  IF(J.EQ.JL) THEN
        DO L=1,NL
	  B(I,K,l)=0.0
	  END DO
	  ELSE
        R0=Q(I,JR,K,1)
        U0=Q(I,JR,K,2)/R0
        V0=Q(I,JR,K,3)/R0
        W0=Q(I,JR,K,4)/R0
        E0=Q(I,JR,K,5)/R0
        P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2))
        AS=SQRT(GAMMA*P0/R0)
	  UVW=0.5*(U0*U0+V0*V0+W0*W0)
        XX=0.5*(XK(2,1,I,JR,K)+XK(2,1,I,J,K))
        YY=0.5*(XK(2,2,I,JR,K)+XK(2,2,I,J,K))
        ZZ=0.5*(XK(2,3,I,JR,K)+XK(2,3,I,J,K))
C        TT=0.5*(XK(2,4,I,JR,K)+XK(2,4,I,J,K))
        QQ=XX*U0+YY*V0+ZZ*W0
	  RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
C        RHO=ABS(TT+QQ)+AS*RR
        RHO=ABS(QQ)+AS*RR
        VNU=2.*(ZMUL(I,JR,K)+VISTURB(I,JR,K)*REYNOLDS)*
     &      (RR1**2+RR2**2+RR3**2)/(R0*REYNOLDS)
C        VNU=0.
C	  SPEC=0.5*(TT+QQ-RHO)
	  SPEC=0.5*(QQ-RHO)-0.5*VNU
c	  SPEC=-0.5*RHO-0.5*VNU
	  S1=Q_N(I,JR,K,1)
	  S2=Q_N(I,JR,K,2)
	  S3=Q_N(I,JR,K,3)
	  S4=Q_N(I,JR,K,4)
	  S5=Q_N(I,JR,K,5)
	  A5=0.5*(XX*S2+YY*S3+ZZ*S4-QQ*S1)
	  A2=(GAMMA-1.)/2.*(UVW*S1-(U0*S2+V0*S3+W0*S4)+S5)
        B(I,K,1)=A5+SPEC*S1
	  B(I,K,2)=XX*A2+U0*A5+SPEC*S2
	  B(I,K,3)=YY*A2+V0*A5+SPEC*S3
	  B(I,K,4)=ZZ*A2+W0*A5+SPEC*S4
	  B(I,K,5)=QQ*A2+(GAMMA*E0-(GAMMA-1.)*UVW)*A5+SPEC*S5
        END IF
C-------------------------I------------------------------
	  IR=I+1
	  IF(I.EQ.IL) THEN
        DO L=1,NL
	  A(I,K,L)=0.
	  END DO
	  ELSE
        R0=Q(IR,J,K,1)
        U0=Q(IR,J,K,2)/R0
        V0=Q(IR,J,K,3)/R0
        W0=Q(IR,J,K,4)/R0
        E0=Q(IR,J,K,5)/R0
        P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2))
        AS=SQRT(GAMMA*P0/R0)
	  UVW=0.5*(U0*U0+V0*V0+W0*W0)
        XX=0.5*(XK(1,1,IR,J,K)+XK(1,1,I,J,K))
        YY=0.5*(XK(1,2,IR,J,K)+XK(1,2,I,J,K))
        ZZ=0.5*(XK(1,3,IR,J,K)+XK(1,3,I,J,K))
C        TT=0.5*(XK(1,4,IR,J,K)+XK(1,4,I,J,K))
        QQ=XX*U0+YY*V0+ZZ*W0
	  RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
C        RHO=ABS(TT+QQ)+AS*RR
        RHO=ABS(QQ)+AS*RR
        VNU=2.*(ZMUL(IR,J,K)+VISTURB(IR,J,K)*REYNOLDS)*
     &      (RR1**2+RR2**2+RR3**2)/(R0*REYNOLDS)
C        VNU=0.
C	  SPEC=0.5*(TT+QQ-RHO)
	  SPEC=0.5*(QQ-RHO)-0.5*VNU
c	  SPEC=-0.5*RHO-0.5*VNU
	  S1=Q_N(IR,J,K,1)
	  S2=Q_N(IR,J,K,2)
	  S3=Q_N(IR,J,K,3)
	  S4=Q_N(IR,J,K,4)
	  S5=Q_N(IR,J,K,5)
	  A5=0.5*(XX*S2+YY*S3+ZZ*S4-QQ*S1)
	  A2=(GAMMA-1.)/2.*(UVW*S1-(U0*S2+V0*S3+W0*S4)+S5)
        A(I,K,1)=A5+SPEC*S1
	  A(I,K,2)=XX*A2+U0*A5+SPEC*S2
	  A(I,K,3)=YY*A2+V0*A5+SPEC*S3
	  A(I,K,4)=ZZ*A2+W0*A5+SPEC*S4
	  A(I,K,5)=QQ*A2+(GAMMA*E0-(GAMMA-1.)*UVW)*A5+SPEC*S5
        END IF
C------------------------K--------------------------
	  KR=K+1
        IF(K.EQ.KL) THEN
	  DO L=1,NL
        C(I,K,L)=0.
	  END DO
	  ELSE
        R0=Q(I,J,KR,1)
        U0=Q(I,J,KR,2)/R0
        V0=Q(I,J,KR,3)/R0
        W0=Q(I,J,KR,4)/R0
        E0=Q(I,J,KR,5)/R0
        P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2)) 
        AS=SQRT(GAMMA*P0/R0)
	  UVW=0.5*(U0*U0+V0*V0+W0*W0)
        XX=0.5*(XK(3,1,I,J,KR)+XK(3,1,I,J,K))
        YY=0.5*(XK(3,2,I,J,KR)+XK(3,2,I,J,K))
        ZZ=0.5*(XK(3,3,I,J,KR)+XK(3,3,I,J,K))
C        TT=0.5*(XK(3,4,I,J,KR)+XK(3,4,I,J,K))
        QQ=XX*U0+YY*V0+ZZ*W0
	  RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
C        RHO=ABS(TT+QQ)+AS*RR
        RHO=ABS(QQ)+AS*RR
        VNU=2.*(ZMUL(I,J,KR)+VISTURB(I,J,KR)*REYNOLDS)*
     &      (RR1**2+RR2**2+RR3**2)/(R0*REYNOLDS)
C        VNU=0.
C	  SPEC=0.5*(TT+QQ-RHO)
	  SPEC=0.5*(QQ-RHO)-0.5*VNU
c	  SPEC=-0.5*RHO-0.5*VNU
	  S1=Q_N(I,J,KR,1)
	  S2=Q_N(I,J,KR,2)
	  S3=Q_N(I,J,KR,3)
	  S4=Q_N(I,J,KR,4)
	  S5=Q_N(I,J,KR,5)
	  A5=0.5*(XX*S2+YY*S3+ZZ*S4-QQ*S1)
	  A2=(GAMMA-1.)/2.*(UVW*S1-(U0*S2+V0*S3+W0*S4)+S5)
        C(I,K,1)=A5+SPEC*S1
	  C(I,K,2)=XX*A2+U0*A5+SPEC*S2
	  C(I,K,3)=YY*A2+V0*A5+SPEC*S3
	  C(I,K,4)=ZZ*A2+W0*A5+SPEC*S4
	  C(I,K,5)=QQ*A2+(GAMMA*E0-(GAMMA-1.)*UVW)*A5+SPEC*S5
        END IF
        PHI=PHI1*DT(I,J,K)/VOL(I,J,K)

c        PHI=PHI1*VOL(I,J,K)/DT(I,J,K)
	  DO L=1,NL
	  Q_N(I,J,K,L)=Q_N(I,J,K,L)-PHI*
     $	          (A(I,K,L)+B(I,K,L)+C(I,K,L))*RH(I,J,K)
        END DO
220     CONTINUE
210     CONTINUE
200     CONTINUE
        
        RETURN
        END
C-----------------------------------------------------------------------
