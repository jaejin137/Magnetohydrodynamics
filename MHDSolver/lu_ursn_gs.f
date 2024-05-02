      subroutine lu_ursn_gs(dt, x, y, z, vol, q, q_n, 
     $     rflux, sflux, tflux,
     $     bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, delta,
     $     rhs, il, jl,
     $     kl, nl, dim1, dim2, ilower, iupper, jlower, jupper, klower,
     $     kupper, visturb, dstep, checktime, wallfunc, ke,
     $     qt, sr, udi, vdi, wdi, udj, vdj, wdj, udk, vdk,
     $     wdk, cst_nri, ced_nri, cst_nrj, ced_nrj, tko, qiii, diver, 
     $     control)
c     
c     use Upwind Relaxation Sweeping(urs)  with Gauss-Seidel
c     to invert the implicit LHS, renamed from original urs.f
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

      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z

      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol

      double precision, intent(in), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q,qiii

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
      
      double precision,intent(in)::
     $     diver(ilower:iupper,jlower:jupper,klower:kupper)

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

c
c     SA 1eq turbulent model contants
c
      double precision, intent(in):: tko

      double precision, intent(out), dimension(ilower:iupper,
     $     jlower:jupper, klower:kupper, nl)::
     $     q_n

      double precision, intent(out)::
     $     cst_nri(nl,nl,jl,kl,5), ced_nri(nl,nl,jl,kl,5),
     $     cst_nrj(nl,nl,il,kl,5), ced_nrj(nl,nl,il,kl,5)
c     
c     LOCAL VARIABLES
c     
      double precision, dimension(nl,nl,dim1)::
     $     capa, capb, capc,
     $     f(nl,dim1),
     $     soln(nl,dim1)

c..   coefficient matrices

      double precision, dimension(nl,nl,il,jl,kl)::
     $     a_minus, a_plus,
     $     b_minus, b_plus,
     $     c_minus, c_plus,
     $     b_bar

      integer::
     $     n,
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


	integer:: index,imax,ip,jp,kp,l,l2
	double precision:: xx,xy,xz
	double precision:: r0,u0,v0,w0,e0,p0,t0,m0,nu,as,qq,
     $                   rr,specx,specy,specz,dtvol,
     $                   cu,gam1,gam2

	double precision, dimension(0:il+1,0:jl+1,0:kl+1)::
     $                  d_bar

	double precision, 
     $ dimension(3,4,ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $                  xk

        DOUBLE PRECISION FIRST(NL)
      
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
      a_minus = 0.d0
      a_plus = 0.d0
      b_minus = 0.d0
      b_plus = 0.d0
      c_minus = 0.d0
      c_plus = 0.d0
      b_bar = 0.d0

c ... calculate the matrices in xi-direction




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
        
c-----------------------evaluate A,B,C---------------------

        gam1=gamma-1.
        gam2=2.-gamma
        do i=1,il
           ip=i+1
           do j=1,jl
              jp=j+1
              do k=1,kl
                 kp=k+1
                 R0=Q(I,J,K,1)
                 U0=Q(I,J,K,2)/R0
                 V0=Q(I,J,K,3)/R0
                 W0=Q(I,J,K,4)/R0
                 E0=Q(I,J,K,5)/R0
                 XX=0.5*(XK(1,1,I,J,K)+XK(1,1,IP,J,K))
                 XY=0.5*(XK(1,2,I,J,K)+XK(1,2,IP,J,K))
                 XZ=0.5*(XK(1,3,I,J,K)+XK(1,3,IP,J,K))
                 P0=(GAMMA-1.)*(R0*E0-0.5*R0*(U0**2+V0**2+W0**2))
                 AS=SQRT(GAMMA*P0/R0)
                 RR=SQRT(XX*XX+XY*XY+XZ*XZ)
                 QQ=0.5*(U0**2+V0**2+W0**2)
                 CU=XX*U0+XY*V0+XZ*W0
                 b_bar(1,1,i,j,k)=0.
                 b_bar(1,2,i,j,k)=xx
                 b_bar(1,3,i,j,k)=xy
                 b_bar(1,4,i,j,k)=xz
                 b_bar(1,5,i,j,k)=0.

                  
                 b_bar(2,1,i,j,k)=-U0*CU+gam1*XX*QQ
                 b_bar(2,2,i,j,k)=CU+XX*gam2*U0
                 b_bar(2,3,i,j,k)=U0*XY-GAM1*XX*V0
                 b_bar(2,4,i,j,k)=U0*XZ-GAM1*XX*W0
                 b_bar(2,5,i,j,k)=GAM1*XX

                 b_bar(3,1,i,j,k)=-V0*CU+gam1*XY*QQ
                 b_bar(3,2,i,j,k)=V0*XX-GAM1*XY*U0
                 b_bar(3,3,i,j,k)=CU+GAM2*XY*V0
                 b_bar(3,4,i,j,k)=V0*XZ-GAM1*XY*W0
                 b_bar(3,5,i,j,k)=GAM1*XY

                 b_bar(4,1,i,j,k)=-W0*CU+gam1*XZ*QQ
                 b_bar(4,2,i,j,k)=W0*XX-GAM1*XZ*U0
                 b_bar(4,3,i,j,k)=W0*XY-GAM1*XZ*V0
                 b_bar(4,4,i,j,k)=CU+GAM2*XZ*W0
                 b_bar(4,5,i,j,k)=GAM1*XZ

                 b_bar(5,1,i,j,k)=-GAMMA*CU*E0+2*GAM1*CU*QQ
                 b_bar(5,2,i,j,k)=-GAM1*CU*U0+XX*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,3,i,j,k)=-GAM1*CU*V0+XY*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,4,i,j,k)=-GAM1*CU*W0+XZ*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,5,i,j,k)=GAMMA*CU

                 do l=1,nl
                    do l2=1,nl
                       a_minus(l,l2,i,j,k)=0.5*b_bar(l,l2,i,j,k) 
                       a_plus(l,l2,i,j,k)=0.5*b_bar(l,l2,i,j,k)
                    end do
                 end do


                 t0 = p0 * gamma * machinf * machinf / r0
                 if ( suther )  then
                    m0 = dsqrt(t0)**3 * (1.d0 + tref) / (t0 + tref)
                 else
                    m0 = t0
                 end if
                 nu=m0+visturb(i,j,k)*reynolds


c                 nu=nu*rr*2/r0/reynolds/vol(i,j,k)
                 nu=nu*rr*2/r0/reynolds
                 specx=0.5*(abs(cu)+rr*as)+nu
                 do l=1,nl
                    a_minus(l,l,i,j,k)=a_minus(l,l,i,j,k)+specx
                    a_plus(l,l,i,j,k)=a_plus(l,l,i,j,k)-specx
                 end do

c-----------------y
                 XX=0.5*(XK(2,1,I,J,K)+XK(2,1,I,JP,K))
                 XY=0.5*(XK(2,2,I,J,K)+XK(2,2,I,JP,K))
                 XZ=0.5*(XK(2,3,I,J,K)+XK(2,3,I,JP,K))
                 RR=SQRT(XX*XX+XY*XY+XZ*XZ)
                 CU=XX*U0+XY*V0+XZ*W0
                 b_bar(1,1,i,j,k)=0.
                 b_bar(1,2,i,j,k)=xx
                 b_bar(1,3,i,j,k)=xy
                 b_bar(1,4,i,j,k)=xz
                 b_bar(1,5,i,j,k)=0.

                  
                 b_bar(2,1,i,j,k)=-U0*CU+gam1*XX*QQ
                 b_bar(2,2,i,j,k)=CU+XX*gam2*U0
                 b_bar(2,3,i,j,k)=U0*XY-GAM1*XX*V0
                 b_bar(2,4,i,j,k)=U0*XZ-GAM1*XX*W0
                 b_bar(2,5,i,j,k)=GAM1*XX

                 b_bar(3,1,i,j,k)=-V0*CU+gam1*XY*QQ
                 b_bar(3,2,i,j,k)=V0*XX-GAM1*XY*U0
                 b_bar(3,3,i,j,k)=CU+GAM2*XY*V0
                 b_bar(3,4,i,j,k)=V0*XZ-GAM1*XY*W0
                 b_bar(3,5,i,j,k)=GAM1*XY

                 b_bar(4,1,i,j,k)=-W0*CU+gam1*XZ*QQ
                 b_bar(4,2,i,j,k)=W0*XX-GAM1*XZ*U0
                 b_bar(4,3,i,j,k)=W0*XY-GAM1*XZ*V0
                 b_bar(4,4,i,j,k)=CU+GAM2*XZ*W0
                 b_bar(4,5,i,j,k)=GAM1*XZ

                 b_bar(5,1,i,j,k)=-GAMMA*CU*E0+2*GAM1*CU*QQ
                 b_bar(5,2,i,j,k)=-GAM1*CU*U0+XX*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,3,i,j,k)=-GAM1*CU*V0+XY*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,4,i,j,k)=-GAM1*CU*W0+XZ*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,5,i,j,k)=GAMMA*CU

                 do l=1,nl
                    do l2=1,nl
                       b_minus(l,l2,i,j,k)=0.5*b_bar(l,l2,i,j,k) 
                       b_plus(l,l2,i,j,k)=0.5*b_bar(l,l2,i,j,k)
                    end do
                 end do
                 
c                 nu=nu*rr*2/r0/reynolds/vol(i,j,k)
                 nu=nu*rr*2/r0/reynolds
                 specy=0.5*(abs(cu)+rr*as)+nu
                 do l=1,nl
                    b_minus(l,l,i,j,k)=b_minus(l,l,i,j,k)+specy
                    b_plus(l,l,i,j,k)=b_plus(l,l,i,j,k)-specy
                 end do
c--------------------z
                 XX=0.5*(XK(3,1,I,J,K)+XK(3,1,I,J,KP))
                 XY=0.5*(XK(3,2,I,J,K)+XK(3,2,I,J,KP))
                 XZ=0.5*(XK(3,3,I,J,K)+XK(3,3,I,J,KP))
                 RR=SQRT(XX*XX+XY*XY+XZ*XZ)
                 CU=XX*U0+XY*V0+XZ*W0
                 b_bar(1,1,i,j,k)=0.
                 b_bar(1,2,i,j,k)=xx
                 b_bar(1,3,i,j,k)=xy
                 b_bar(1,4,i,j,k)=xz
                 b_bar(1,5,i,j,k)=0.

                  
                 b_bar(2,1,i,j,k)=-U0*CU+gam1*XX*QQ
                 b_bar(2,2,i,j,k)=CU+XX*gam2*U0
                 b_bar(2,3,i,j,k)=U0*XY-GAM1*XX*V0
                 b_bar(2,4,i,j,k)=U0*XZ-GAM1*XX*W0
                 b_bar(2,5,i,j,k)=GAM1*XX

                 b_bar(3,1,i,j,k)=-V0*CU+gam1*XY*QQ
                 b_bar(3,2,i,j,k)=V0*XX-GAM1*XY*U0
                 b_bar(3,3,i,j,k)=CU+GAM2*XY*V0
                 b_bar(3,4,i,j,k)=V0*XZ-GAM1*XY*W0
                 b_bar(3,5,i,j,k)=GAM1*XY

                 b_bar(4,1,i,j,k)=-W0*CU+gam1*XZ*QQ
                 b_bar(4,2,i,j,k)=W0*XX-GAM1*XZ*U0
                 b_bar(4,3,i,j,k)=W0*XY-GAM1*XZ*V0
                 b_bar(4,4,i,j,k)=CU+GAM2*XZ*W0
                 b_bar(4,5,i,j,k)=GAM1*XZ

                 b_bar(5,1,i,j,k)=-GAMMA*CU*E0+2*GAM1*CU*QQ
                 b_bar(5,2,i,j,k)=-GAM1*CU*U0+XX*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,3,i,j,k)=-GAM1*CU*V0+XY*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,4,i,j,k)=-GAM1*CU*W0+XZ*(GAMMA*E0-GAM1*QQ)
                 b_bar(5,5,i,j,k)=GAMMA*CU

                 do l=1,nl
                    do l2=1,nl
                       c_minus(l,l2,i,j,k)=0.5*b_bar(l,l2,i,j,k) 
                       c_plus(l,l2,i,j,k)=0.5*b_bar(l,l2,i,j,k)
                    end do
                 end do

c                 nu=nu*rr*2/r0/reynolds/vol(i,j,k)
                 nu=nu*rr*2/r0/reynolds
                 specz=0.5*(abs(cu)+rr*as)+nu
                 do l=1,nl
                    c_minus(l,l,i,j,k)=c_minus(l,l,i,j,k)+specz
                    c_plus(l,l,i,j,k)=c_plus(l,l,i,j,k)-specz
                 end do
c-------------------------------
                 d_bar(i,j,k)=2.*(specx+specy+specz)

                 dtvol=dt(i,j,k)/vol(i,j,k)
                 a_minus(:,:,i,j,k)=-a_minus(:,:,i,j,k)*dtvol
                 a_plus(:,:,i,j,k)=a_plus(:,:,i,j,k)*dtvol
                 b_minus(:,:,i,j,k)=-b_minus(:,:,i,j,k)*dtvol
                 b_plus(:,:,i,j,k)=b_plus(:,:,i,j,k)*dtvol
                 c_minus(:,:,i,j,k)=-c_minus(:,:,i,j,k)*dtvol
                 c_plus(:,:,i,j,k)=c_plus(:,:,i,j,k)*dtvol
                 d_bar(i,j,k)=d_bar(i,j,k)*dtvol
c------------------------------
              end do
           end do
        end do

        b_bar=0.

        do l=1,nl
           b_bar(l,l,:,:,:)=1.0+d_bar(:,:,:)
        end do



c..   start Gauss-Seidel Relaxation Sweeping
      
c..   eta-constant line

      if(jl.gt.1) then

         do m=1,iter_gs
            do k = 1, kl
               do i = 1, il
                  
                  if(mod(m,2).eq.1) then
                     is =i
                  else
                     is =il+1-i
                  end if
                  
                  do j=1,jl
                  capb(:,:,j) = b_minus(:,:,is,j-1,k) 
                  capc(:,:,j) = b_plus(:,:,is,j+1,k) 
                  capa(:,:,j) = b_bar(:,:,is,j,k)
                  end do
c     
c     form right-hand-side vector f(n,j)
c     
c..   initilize zero
                  f = 0.d0

                  cz1=1.0
                  cz2=1.0
                  cx1=1.0
                  cx2=1.0

                  if(k.eq.1)  cz1=0.0
                  if(k.eq.kl) cz2=0.0
                  if(i.eq.1)  cx1=0.0
                  if(i.eq.il) cx2=0.0

                  do j = 1,jl
                     do n = 1,nl
                        do nj = 1,nl
                           f(n,j) = f(n,j) 
     >                          -cx1*
     $                          a_minus(n,nj,is-1,j,k)*q_n(is-1,j,k,nj)
     >                          -cx2*
     $                          a_plus(n,nj,is+1,j,k)*q_n(is+1,j,k,nj)
     >                          -cz1*
     $                          c_minus(n,nj,is,j,k-1)*q_n(is,j,k-1,nj)
     >                          -cz2*
     $                          c_plus(n,nj,is,j,k+1)*q_n(is,j,k+1,nj)
                        end do
                     end do
                  end do
                  
                  do n = 1, nl
                     f(n,:jl) = f(n,:jl) + rhs(is,:,k,n) 
                  end do

c     solve linear system
c     
                  nblocks = jl

                  if(bc_eta_lower(i,k).ne.bcperi
     $               .or.bc_eta_lower(i,k).ne.20) then

C------------------------->				 
		 DO J=1,JL
	                FIRST(:)=0.
	                IF (J.EQ.1) THEN
	                   FIRST(:)=0.
		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $				  CAPB(N,NJ,J)*Q_N(IS,J-1,K,NJ)
	                      END DO
	                   END DO
	                END IF
	                FIRST(:)=F(:,J)-FIRST(:)
			DO N=1,NL
	                   Q_N(IS,J,K,N)=FIRST(N)/B_BAR(N,N,IS,J,K)
	                END DO
	             ENDDO  !END OF J

C-----<-------------------------
			 DO J=JL,1,-1
	                FIRST(:)=0.
	                IF (J.EQ.JL) THEN
	                   FIRST(:)=0.
		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $				  CAPC(N,NJ,J)*Q_N(IS,J+1,K,NJ)
	                      END DO
	                   END DO
	                END IF
			DO N=1,NL
	                   Q_N(IS,J,K,N)=Q_N(IS,J,K,N)-
     $                             FIRST(N)/B_BAR(N,N,IS,J,K)
	                END DO
	             ENDDO  !END OF J

C                     call block(nblocks, nl, dim1, capa, capb, capc,
C     $                    soln, f)
                  else
                     soln(:,:jl) = f
                     call pbtrip(capb(:,:,:jl), capa(:,:,:jl),
     $                    capc(:,:,:jl), soln(:,:jl), 1, jl, nl)

                  do n = 1,nl
                     q_n(is,1:jl,k,n) = soln(n,1:jl)
                  end do

                  end if
c     
               end do
            end do
         end do

      end if


c      return



c..   xi-constant line

      if(il.gt.1) then

         do m=1,iter_gs
            do k = 1, kl
               do j = 1, jl
                  
                  if(mod(m,2).eq.1) js =j
                  if(mod(m,2).eq.0) js =jl+1-j

                  do i = 1,il
                     capb(:,:,i) = a_minus(:,:,i-1,js,k) 
                     capc(:,:,i) = a_plus(:,:,i+1,js,k) 
                     capa(:,:,i) = b_bar(:,:,i,js,k)
                  end do
c     
c     form right-hand-side vector f(n,j)
c     
                  do n = 1,nl
                     f(n,:il) = rhs(:il,js,k,n) 
                  end do

                  cy1=1.0
                  cy2=1.0
                  cz1=1.0
                  cz2=1.0
                  
                  if(j.eq.1)  cy1=0.0
                  if(j.eq.jl) cy2=0.0
                  if(k.eq.1)  cz1=0.0
                  if(k.eq.kl) cz2=0.0

                  do i = 1,il
                     do n = 1,nl
                        do nj = 1,nl
                           f(n,i) = f(n,i) 
     >                          -cy1*
     $                          b_minus(n,nj,i,js-1,k)*q_n(i,js-1,k,nj)
     >                          -cy2*
     $                          b_plus(n,nj,i,js+1,k)*q_n(i,js+1,k,nj)
     >                          -cz1*
     $                          c_minus(n,nj,i,js,k-1) *q_n(i,js,k-1,nj)
     >                          -cz2*
     $                          c_plus(n,nj,i,js,k+1)*q_n(i,js,k+1,nj)
                        end do
                     end do
                  end do

c     
c     solve linear system
c     
                  nblocks = il

                  if(bc_xi_lower(j,k).ne.bcperi
     $               .or.bc_xi_lower(j,k).ne.20) then

C------------------------->				 
		 DO I=1,IL
	                FIRST(:)=0.
	                IF (I.EQ.1) THEN
	                   FIRST(:)=0.
		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $				   CAPB(N,NJ,I)*Q_N(I-1,JS,K,NJ)
	                      END DO
	                   END DO
	                END IF
	                FIRST(:)=F(:,I)-FIRST(:)
			DO N=1,NL
	                   Q_N(I,JS,K,N)=FIRST(N)/B_BAR(N,N,I,JS,K)
	                END DO
	             ENDDO  !END OF I

C-----<-------------------------
		 DO I=IL,1,-1
	                FIRST(:)=0.
	                IF (I.EQ.IL) THEN
	                   FIRST(:)=0.
		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $				   CAPC(N,NJ,I)*Q_N(I+1,JS,K,NJ)
	                      END DO
	                   END DO
	                END IF
			DO N=1,NL
	                   Q_N(I,JS,K,N)=Q_N(I,JS,K,N)-
     $                            FIRST(N)/B_BAR(N,N,I,JS,K)
	                END DO
	             ENDDO  !END OF I


C                     call block(nblocks, nl, dim1, capa, capb, capc,
C     $                    soln, f)
                  else
                     soln(:,:il) = f
                     call  pbtrip(capb(:,:,:il),capa(:,:,:il),
     $                    capc(:,:,:il), soln(:,:il),1, il, nl)


                  do   i = 1,il
                     do   n = 1,nl
                        q_n(i,js,k,n) = soln(n,i)
                     end do
                  end do

                  end if

C                  do   i = 1,il
C                     do   n = 1,nl
C                        q_n(i,js,k,n) = soln(n,i)
C                     end do
C                  end do

c     
               end do
            end do
         end do

      end if    

c..   zeta-constant line

      if(kl.gt.1) then

         do m=1,iter_gs
            do j = 1, jl
               do i = 1, il
                  
                  if(mod(m+j,2).eq.1) is =i
                  if(mod(m+j,2).eq.0) is =il+1-i

                  do k = 1, kl
                     capb(:,:,k) = c_minus(:,:,is,j,k-1) 
                     capc(:,:,k) = c_plus(:,:,is,j,k+1) 
                     capa(:,:,k) = b_bar(:,:,is,j,k)
                  end do
c     
c     form right-hand-side vector f(n,j)
c     
c..   initilize zero

                  f = 0.d0
                  
                  cx1=1.0
                  cx2=1.0
                  cy1=1.0
                  cy2=1.0
                  
                  if(i.eq.1)  cx1=0.0
                  if(i.eq.il) cx2=0.0
                  if(j.eq.1)  cy1=0.0
                  if(j.eq.jl) cy2=0.0

                  do  k = 1, kl
                     do   n = 1,nl
                        do   nj = 1,nl
                           f(n,k) = f(n,k) 
     >                          -cx1*
     $                          a_minus(n,nj,is-1,j,k)*q_n(is-1,j,k,nj)
     >                          -cx2*
     $                          a_plus(n,nj,is+1,j,k)*q_n(is+1,j,k,nj)
     >                          -cy1*
     $                          b_minus(n,nj,is,j-1,k) *q_n(is,j-1,k,nj)
     >                          -cy2*
     $                          b_plus(n,nj,is,j+1,k)*q_n(is,j+1,k,nj)
                        end do
                     end do
                  end do

                  do   k = 1, kl
                     do   n = 1, nl
                        f(n,k) = f(n,k) + rhs(is,j,k,n) 
                     end do
                  end do
c     
c     solve linear system
c     
                  nblocks = kl

                  if(bc_zeta_lower(i,j).ne.bcperi
     $               .or.bc_zeta_lower(i,j).ne.20) then


C------------------------->				 
		 DO K=1,KL
	                FIRST(:)=0.
	                IF (K.EQ.1) THEN
	                   FIRST(:)=0.
		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $				   CAPB(N,NJ,K)*Q_N(IS,J,K-1,NJ)
	                      END DO
	                   END DO
	                END IF
	                FIRST(:)=F(:,K)-FIRST(:)
			DO N=1,NL
	                   Q_N(IS,J,K,N)=FIRST(N)/B_BAR(N,N,IS,J,K)
	                END DO
	             ENDDO  !END OF K

C-----<-------------------------
			 DO K=KL,1,-1
	                FIRST(:)=0.
	                IF (K.EQ.KL) THEN
	                   FIRST(:)=0.
		            ELSE
	                   DO N=1,NL
	                      FIRST(N)=0.
	                      DO NJ=1,NL
	                         FIRST(N)=FIRST(N)+
     $				  CAPC(N,NJ,K)*Q_N(IS,J,K+1,NJ)
	                      END DO
	                   END DO
	                END IF
                     	DO N=1,NL
	                   Q_N(IS,J,K,N)=Q_N(IS,J,K,N)-
     $                            FIRST(N)/B_BAR(N,N,IS,J,K)
	                END DO
	             ENDDO  !END OF K





C                     call block(nblocks, nl, dim1, capa, capb, capc,
C     $                    soln, f)
                  else
                     soln(:,:kl) = f
                     call pbtrip(capb(:,:,:kl), capa(:,:,:kl),
     $                    capc(:,:,:kl), soln(:,:kl), 1, kl, nl)
                  do   k = 1,kl
                     do   n = 1,nl
                        q_n(is,j,k,n) = soln(n,k)
                     end do
                  end do

                  end if
c     
C                  do   k = 1,kl
C                     do   n = 1,nl
C                        q_n(is,j,k,n) = soln(n,k)
C                     end do
C                  end do
c     
               end do
            end do
         end do

      end if

      return
      end
