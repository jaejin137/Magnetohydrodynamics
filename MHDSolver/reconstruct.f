      subroutine reconstruct(index,imax, iside, dim2, il, jl, kl, nl,
     $     lx, ly, lz,
     $     mx, my, mz, nx, ny, nz, rho, rhou, rhov, rhow, rhoe,
     $     capul, capur, capvl, capvr, capwl, capwr, el, er, pl, pr,
     $     ql, qr, ul, ur, vl, vr, wl, wr, rhoinvl, rhoinvr, rhol, rhor,
     $     rhoul, rhour, rhovl, rhovr, rhowl, rhowr, rhoel, rhoer,
     $     rhotkl, rhotkr, rhotwl, rhotwr,
     $     sqrtinv, sqrtrhol, sqrtrhor, tenthalpyl, tenthalpyr,
     $     bc_lower, bc_upper, rhotk, rhotw, tkl, tkr, twl, twr, check,
     $     ke, qt, qtl, qtr, ilower, tl, tr,
     $     diver, control)
c
c     linear reconstruction to the cell faces 
c     using conservative variables
c
      use datain

      implicit none

      type (datain_type), intent(in)::control
      
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     index,
     $     imax,
     $     iside,               ! MUSCL order control (due to lhs and rhs)
     $     dim2,
     $     il, jl, kl,
     $     nl, ilower,
     $     bc_lower, bc_upper
            
      integer::
     $     limiter,
     $     moving,               ! 0, rest, 1, forced, 2, self
     $     precondition              

      double precision::
     $     epsfactor,
     $     kfactor,             ! control order when fi ne 0
     $     gamma, machinf,
     $     ronum

      double precision, dimension(ilower+1:dim2), intent(in)::
     $     lx, ly, lz,
     $     mx, my, mz,
     $     nx, ny, nz
 
      double precision, dimension(ilower:dim2), intent(in)::
     $     rho,
     $     rhou,
     $     rhov,
     $     rhow,
     $     rhoe,
     $     rhotk,
     $     rhotw,
     $     qt(ilower:dim2,3),diver

      double precision, dimension(dim2), intent(out)::
     $     capul, capur,
     $     capvl, capvr,
     $     capwl, capwr,        
     >     el, er,
     $     pl, pr,
     $     ql, qr,
     >     rhoinvl, rhoinvr,
     $     rhol, rhor,
     $     rhoul, rhour,
     $     rhovl, rhovr,
     $     rhowl, rhowr,
     $     rhoel, rhoer,
     $     rhotkl, rhotkr,
     $     rhotwl, rhotwr,
     $     sqrtinv,
     $     sqrtrhol, sqrtrhor,
     >     tenthalpyl, tenthalpyr,
     $     ul, ur,
     $     vl, vr,
     $     wl, wr,
     $     tkl, tkr,
     $     twl, twr,
     $     qtl, qtr,
     $     tl,tr  

      logical, intent(in)::
     $     check

      logical, intent(in)::
     $     ke
c
c     LOCAL VARIABLES
c
      double precision, dimension(ilower:dim2,nl+3)::
     $     dqm, dqp, slimit, t1, t2, t3, rp, rm, phip, phim

      double precision, dimension(dim2,3)::
     $     qll, qrr
c
      integer 
     $     i, ii,               ! iteration index for faces
     $     n,                   ! iteration index for variables
     $     first                ! first order index region on boundaries

      double precision::
     $     fi,                  ! fi=1, 1-order; fi>=2, 2- or higher o$|+ c
     $     omega                ! coefficient of minmod limiter        |        integer

      logical::
     $     bndl, bndu
c
      double precision::
     $     psil, psir, thetal, thetar
c
c *** SUBROUTINE START ***
c
      limiter=control%limiter
      moving=control%moving
      precondition=control%precondition
      epsfactor=control%epsfactor
      gamma=control%gamma
      kfactor=control%kfactor
      machinf=control%machinf
      ronum=control%ronum
c--------------------------------------------------------
      qtl = 0.d0
      qtr = 0.d0
      if (precondition.ge.1.and.nl.ge.6) then
         write(*,*)'from reconstruct'
         write(*,*)'for this case, the code is not finished'
         stop
      end if
c
      if (iside.ge.3) then
        do i = ilower+1,imax-ilower
          dqm(i,1)  = rho(i) - rho(i-1)
          dqm(i,2)  = rhou(i) - rhou(i-1)
          dqm(i,3)  = rhov(i) - rhov(i-1)
          dqm(i,4)  = rhow(i) - rhow(i-1)
          dqm(i,5)  = rhoe(i) - rhoe(i-1)
          if(nl.eq.6) then
            dqm(i,6)  = rhotk(i) - rhotk(i-1)
          else if(nl.eq.7) then
            dqm(i,6)  = rhotk(i) - rhotk(i-1)
            dqm(i,7)  = rhotw(i) - rhotw(i-1)
          end if
          if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            dqm(i,nl+1)  = qt(i,1) - qt(i-1,1)
            dqm(i,nl+2)  = qt(i,2) - qt(i-1,2)
            dqm(i,nl+3)  = qt(i,3) - qt(i-1,3)
          end if
        end do
        call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                rho,dqm(:,1),rhol,rhor,control)
        call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                rhou,dqm(:,2),rhoul,rhour,control)
        call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                rhov,dqm(:,3),rhovl,rhovr,control)
        call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                rhow,dqm(:,4),rhowl,rhowr,control)
        call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                rhoe,dqm(:,5),rhoel,rhoer,control)
        if(nl.eq.6) then
          call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                  rhotk,dqm(:,6),rhotkl,rhotkr,control)
        else if(nl.eq.7) then
          call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                  rhotk,dqm(:,6),rhotkl,rhotkr,control)
          call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                  rhotw,dqm(:,7),rhotwl,rhotwr,control)
        end if
        if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
          call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                 qt(:,1),dqm(:,nl+1),qll(:,1),qrr(:,1),control)
          call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                 qt(:,2),dqm(:,nl+2),qll(:,2),qrr(:,2),control)
          call weno_5(iside,imax,ilower,dim2,bc_lower,bc_upper,
     &                 qt(:,3),dqm(:,nl+3),qll(:,3),qrr(:,3),control)
        end if
      else
        do i = 1,imax
          dqm(i,1)  = rho(i) - rho(i-1)
          dqm(i,2)  = rhou(i) - rhou(i-1)
          dqm(i,3)  = rhov(i) - rhov(i-1)
          dqm(i,4)  = rhow(i) - rhow(i-1)
          dqm(i,5)  = rhoe(i) - rhoe(i-1)
          if(nl.eq.6) then
            dqm(i,6)  = rhotk(i) - rhotk(i-1)
          else if(nl.eq.7) then
            dqm(i,6)  = rhotk(i) - rhotk(i-1)
            dqm(i,7)  = rhotw(i) - rhotw(i-1)
          end if
          if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            dqm(i,nl+1)  = qt(i,1) - qt(i-1,1)
            dqm(i,nl+2)  = qt(i,2) - qt(i-1,2)
            dqm(i,nl+3)  = qt(i,3) - qt(i-1,3)
          end if
        end do
c
        if (bc_lower.ne.7.and.bc_lower.ne.8
     $      .and.bc_lower.ne.10.and.bc_lower.ne.20) then
          do n = 1,nl
            dqm(0,n) = 0.d0
          end do
        else
          i = 0
          dqm(i,1)  = rho(i) - rho(i-1)
          dqm(i,2)  = rhou(i) - rhou(i-1)
          dqm(i,3)  = rhov(i) - rhov(i-1)
          dqm(i,4)  = rhow(i) - rhow(i-1)
          dqm(i,5)  = rhoe(i) - rhoe(i-1)
          if(nl.eq.6) then
            dqm(i,6)  = rhotk(i) - rhotk(i-1)
          else if(nl.eq.7) then
            dqm(i,6)  = rhotk(i) - rhotk(i-1)
            dqm(i,7)  = rhotw(i) - rhotw(i-1)
          end if
          if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            dqm(i,nl+1)  = qt(i,1) - qt(i-1,1)
            dqm(i,nl+2)  = qt(i,2) - qt(i-1,2)
            dqm(i,nl+3)  = qt(i,3) - qt(i-1,3)
          end if
        end if
c
        do n = 1,nl+3
          do i = 0,imax-1
            dqp(i,n) = dqm(i+1,n)
          end do
        end do
c
        if (bc_upper.ne.7.and.bc_upper.ne.8
     $      .and.bc_upper.ne.10.and.bc_upper.ne.20) then
          do n = 1,nl
            dqp(imax,n) = 0.d0
          end do
        else
          i = imax
          dqp(i,1)  = rho(i+1)- rho(i)
          dqp(i,2)  = rhou(i+1) - rhou(i)
          dqp(i,3)  = rhov(i+1) - rhov(i)
          dqp(i,4)  = rhow(i+1) - rhow(i)
          dqp(i,5)  = rhoe(i+1) - rhoe(i)
          if(nl.eq.6) then
            dqp(i,6)  = rhotk(i+1) - rhotk(i)
          else if(nl.eq.7) then
            dqp(i,6)  = rhotk(i+1) - rhotk(i)
            dqp(i,7)  = rhotw(i+1) - rhotw(i)
          end if
          if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            dqp(i,nl+1)  = qt(i+1,1) - qt(i,1)
            dqp(i,nl+2)  = qt(i+1,2) - qt(i,2)
            dqp(i,nl+3)  = qt(i+1,3) - qt(i,3)
          end if
        end if
c
        omega = 1.d0
        phip = 1.d0   
        phim = 1.d0   
        slimit = 1.d0
c
c *** Use a TVD limiter ***
c
        if (limiter.eq.1.or.limiter.eq.2) then
c
          do ii = 1, nl+3
            if(moving.eq.0.and.ii.gt.nl) exit
            do i = 1, imax-1
               rm(i,ii) = dqm(i,ii)/(dqp(i,ii) + epsfactor)    
               rp(i,ii) = dqp(i,ii)/(dqm(i,ii) + epsfactor)    
            end do
          end do
c
          if ( dabs(ronum).gt.1e-9 ) then
            do ii = nl+1, nl+3
              do i  = 1, imax-1
                rm(i,ii) = dqm(i,ii)/(dqp(i,ii)+epsfactor)    
                rp(i,ii) = dqp(i,ii)/(dqm(i,ii)+epsfactor)    
              end do
            end do
          end if
c
        end if
c
        if (limiter.eq.1) then
c
c ... MINMOD Limiter ...
c
          do ii = 1, nl+3
            if(moving.eq.0.and.ii.gt.nl) exit
            do i = 1, imax-1
               phip(i,ii) = max(0.d0, min(1.d0,omega*rp(i,ii) ) )    
               phim(i,ii) = max(0.d0, min(1.d0,omega*rm(i,ii) ) )    
            end do
          end do
c
          if ( dabs(ronum).gt.1e-9 ) then
            do ii = nl+1, nl+3
              do i  = 1, imax-1
                phip(i,ii) = max(0.d0, min(1.d0,omega*rp(i,ii) ) )    
                phim(i,ii) = max(0.d0, min(1.d0,omega*rm(i,ii) ) )    
              end do
            end do
          end if
c
        else if(limiter.eq.2) then
c
c ... SUPERBEE Limiter ...
c
          do ii = 1, nl+3
            if(moving.eq.0.and.ii.gt.nl) exit
            do i = 1, imax-1
               phip(i,ii) = max(0.d0, min( 2.d0*rp(i,ii),1.d0 ),
     $              min( rp(i,ii),2.d0 ) )    
               phim(i,ii) = max(0.d0, min( 2.d0*rm(i,ii),1.d0 ),    
     $              min( rm(i,ii),2.d0 ) )    
            end do
          end do
c
          if ( dabs(ronum).gt.1e-9 ) then
            do ii = nl+1, nl+3
              do i  = 1, imax-1
                phip(i,ii) = max(0.d0, min( 2.d0*rp(i,ii),1.d0 ),
     $                       min( rp(i,ii),2.d0 ) )    
                phim(i,ii) = max(0.d0, min( 2.d0*rm(i,ii),1.d0 ),    
     $                       min( rm(i,ii),2.d0 ) )    
              end do
            end do
          end if

        else if(limiter.eq.3) then
c
c ... Anderson, Thomas & Van Leer ...
c
          do n = 1,nl+3
            if(moving.eq.0.and.n.gt.nl) exit
            do i = 0,imax
               slimit(i,n) = dqp(i,n)**2 + dqm(i,n)**2 + epsfactor
               slimit(i,n) = (2.d0*dqp(i,n)*dqm(i,n) + epsfactor)/
     >              slimit(i,n)
            end do
          end do
c
          if ( dabs(ronum).gt.1e-9 ) then
            do n = nl+1,nl+3
              do i = 0,imax
                slimit(i,n) = dqp(i,n)**2 + dqm(i,n)**2 + epsfactor
                slimit(i,n) = (2.d0*dqp(i,n)*dqm(i,n) + epsfactor)/
     >                        slimit(i,n)
              end do
            end do
          end if

        end if
c
        do n = 1,nl+3
          if(moving.eq.0.and.n.gt.nl) exit
          do i = 0,imax
            t1(i,n) = 0.25d0*slimit(i,n)
            t2(i,n) = 1.d0 - kfactor * slimit(i,n)
            t3(i,n) = 1.d0 + kfactor * slimit(i,n)
          end do
        end do
c
        if ( dabs(ronum).gt.1e-9 ) then
          do n = nl+1,nl+3
            do i = 0,imax
              t1(i,n) = 0.25d0*slimit(i,n)
              t2(i,n) = 1.d0 - kfactor * slimit(i,n)
              t3(i,n) = 1.d0 + kfactor * slimit(i,n)
            end do
          end do
        end if
c
c ... left side ...
c
        first = 3
        do i = 0,imax-1
c        
          fi=1.d0
          if(iside.eq.0) fi=0.d0
          bndl = i.le.first-1.and.(bc_lower.ne.7.and.bc_lower.ne.8
     $           .and.bc_lower.ne.10.and.bc_lower.ne.20)
          bndu = i.ge.(imax-first).and.(bc_upper.ne.7.and.
     $           bc_lower.ne.8.and.bc_upper.ne.10.and.bc_upper.ne.20)
          if(bndl.or.bndu) fi = 0.d0
                    
          rhol(i+1)    = rho(i)   +  fi*t1(i,1)*(t2(i,1)*phip(i,1)*
     >        dqm(i,1) +  t3(i,1)*phim(i,1)*dqp(i,1))
          rhoul(i+1)   = rhou(i)  +  fi*t1(i,2)*(t2(i,2)*phip(i,2)*
     >        dqm(i,2) +  t3(i,2)*phim(i,2)*dqp(i,2))
          rhovl(i+1)   = rhov(i)  +  fi*t1(i,3)*(t2(i,3)*phip(i,3)*
     >        dqm(i,3) +  t3(i,3)*phim(i,3)*dqp(i,3))
          rhowl(i+1)   = rhow(i)  +  fi*t1(i,4)*(t2(i,4)*phip(i,4)*
     >        dqm(i,4) +  t3(i,4)*phim(i,4)*dqp(i,4))
          rhoel(i+1)   = rhoe(i)  +  fi*t1(i,5)*(t2(i,5)*phip(i,5)*
     >        dqm(i,5) +  t3(i,5)*phim(i,5)*dqp(i,5))
          if(nl.eq.6) then
            rhotkl(i+1)   = rhotk(i)  +  fi*t1(i,6)*(t2(i,6)*phip(i,6)*
     >           dqm(i,6) +  t3(i,6)*phim(i,6)*dqp(i,6))
          else if(nl.eq.7) then
            rhotkl(i+1)   = rhotk(i)  +  fi*t1(i,6)*(t2(i,6)*phip(i,6)*
     >           dqm(i,6) +  t3(i,6)*phim(i,6)*dqp(i,6))
            rhotwl(i+1)   = rhotw(i)  +  fi*t1(i,7)*(t2(i,7)*phip(i,7)*
     >           dqm(i,7) +  t3(i,7)*phim(i,7)*dqp(i,7))
          end if
          if( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            qll(i+1,1) = qt(i,1)+fi*t1(i,nl+1)*(t2(i,nl+1)*phip(i,nl+1)*
     >           dqm(i,nl+1) +  t3(i,nl+1)*phim(i,nl+1)*dqp(i,nl+1))
            qll(i+1,2) = qt(i,2)+fi*t1(i,nl+2)*(t2(i,nl+2)*phip(i,nl+2)*
     >           dqm(i,nl+2) +  t3(i,nl+2)*phim(i,nl+2)*dqp(i,nl+2))
            qll(i+1,3) = qt(i,3)+fi*t1(i,nl+3)*(t2(i,nl+3)*phip(i,nl+3)*
     >           dqm(i,nl+3) +  t3(i,nl+3)*phim(i,nl+3)*dqp(i,nl+3))
          end if
        end do
c
c ... right side
c
        do i = 1,imax
c
          fi=1.d0
          if(iside.eq.0) fi=0.d0
          bndl = i.le.first.and.(bc_lower.ne.7.and.bc_lower.ne.8
     >          .and.bc_lower.ne.10.and.bc_lower.ne.20)
          bndu = i.ge.(imax-first+1).and.(bc_upper.ne.7
     $        .and.bc_upper.ne.8.and.bc_upper.ne.10.and.bc_upper.ne.20)
          if(bndl.or.bndu) fi = 0.d0
c
          rhor(i) = rho(i) - fi*t1(i,1)*( t2(i,1)*phim(i,1)*
     >        dqp(i,1) + t3(i,1)*phip(i,1)*dqm(i,1) )
          rhour(i) = rhou(i) - fi*t1(i,2)*(t2(i,2)*phim(i,2)*
     >        dqp(i,2) + t3(i,2)*phip(i,2)*dqm(i,2))
          rhovr(i) = rhov(i) - fi*t1(i,3)*(t2(i,3)*phim(i,3)*
     >        dqp(i,3) + t3(i,3)*phip(i,3)*dqm(i,3))
          rhowr(i) = rhow(i) - fi*t1(i,4)*(t2(i,4)*phim(i,4)*
     >        dqp(i,4) + t3(i,4)*phip(i,4)*dqm(i,4))
          rhoer(i) = rhoe(i) - fi*t1(i,5)*(t2(i,5)*phim(i,5)*
     >        dqp(i,5) + t3(i,5)*phip(i,5)*dqm(i,5))
          if(nl.eq.6) then
            rhotkr(i)   = rhotk(i)  -  fi*t1(i,6)*(t2(i,6)*phim(i,6)*
     >           dqp(i,6) +  t3(i,6)*phip(i,6)*dqm(i,6))
          else if(nl.eq.7) then
            rhotkr(i)   = rhotk(i)  -  fi*t1(i,6)*(t2(i,6)*phim(i,6)*
     >           dqp(i,6) +  t3(i,6)*phip(i,6)*dqm(i,6))
            rhotwr(i)   = rhotw(i)  -  fi*t1(i,7)*(t2(i,7)*phim(i,7)*
     >           dqp(i,7) +  t3(i,7)*phip(i,7)*dqm(i,7))
          end if
          if ( moving.ge.1.or.dabs(ronum).gt.1e-9 ) then
            qrr(i,1) = qt(i,1) - fi*t1(i,nl+1)*(t2(i,nl+1)*phim(i,nl+1)*
     >        dqp(i,nl+1) + t3(i,nl+1)*phip(i,nl+1)*dqm(i,nl+1))
            qrr(i,2) = qt(i,2) - fi*t1(i,nl+2)*(t2(i,nl+2)*phim(i,nl+2)*
     >        dqp(i,nl+2) + t3(i,nl+2)*phip(i,nl+2)*dqm(i,nl+2))
            qrr(i,3) = qt(i,3) - fi*t1(i,nl+3)*(t2(i,nl+3)*phim(i,nl+3)*
     >        dqp(i,nl+3) + t3(i,nl+3)*phip(i,nl+3)*dqm(i,nl+3))
          end if
        end do
      end if
c
c----------------------------------------------------------
c
      if (precondition.ge.2) then
      do i = 1,imax
c
        pl(i) = rhol(i)
        ul(i) = rhoul(i)
        vl(i) = rhovl(i)
        wl(i) = rhowl(i)
        tl(i) = rhoel(i)
c
        pr(i) = rhor(i)
        ur(i) = rhour(i)
        vr(i) = rhovr(i)
        wr(i) = rhowr(i)
        tr(i) = rhoer(i)
c
        ql(i) = 0.5d0 * (ul(i)**2 + vl(i)**2 + wl(i)**2)
        qr(i) = 0.5d0 * (ur(i)**2 + vr(i)**2 + wr(i)**2)
c
c----------------------------------
c
	call qstat(pl(i),tl(i),diver(i),rhoinvl(i))
	call qstat(pr(i),tr(i),diver(i),rhoinvr(i))
        rhol(i)=1.0d0/rhoinvl(i)
        rhor(i)=1.0d0/rhoinvr(i)
c
        sqrtrhol(i) = dsqrt(rhol(i))
        sqrtrhor(i) = dsqrt(rhor(i))
        sqrtinv(i) = 1.d0 / (sqrtrhol(i) + sqrtrhor(i))
c
	call enthalpy(tenthalpyl(i),tl(i),rhoinvl(i),pl(i))
	call enthalpy(tenthalpyr(i),tr(i),rhoinvr(i),pr(i))
	tenthalpyl(i)=tenthalpyl(i)+ql(i)
	tenthalpyr(i)=tenthalpyr(i)+qr(i)
	er(i)=tenthalpyr(i)-pr(i)*rhoinvr(i)
	el(i)=tenthalpyl(i)-pl(i)*rhoinvl(i)
c
        rhoul(i)=ul(i) * rhol(i)
        rhovl(i)=vl(i) * rhol(i) 
        rhowl(i)=wl(i) * rhol(i) 
        rhoel(i)=el(i) * rhol(i) 
        rhour(i)=ur(i) * rhor(i) 
        rhovr(i)=vr(i) * rhor(i) 
        rhowr(i)=wr(i) * rhor(i) 
        rhoer(i)=er(i) * rhor(i) 
c
c----------------------------------
c
        if(nl.eq.6) then
           tkl(i) = rhoinvl(i) * rhotkl(i)
           tkr(i) = rhoinvr(i) * rhotkr(i)
        else if(nl.eq.7) then
           tkl(i) = rhoinvl(i) * rhotkl(i)
           twl(i) = rhoinvl(i) * rhotwl(i)
c
           tkr(i) = rhoinvr(i) * rhotkr(i)
           twr(i) = rhoinvr(i) * rhotwr(i)
c
           if(ke) then
              pl(i) = (gamma-1.d0)*(rhoel(i)-rhotkl(i)-rhol(i)*ql(i))
              pr(i) = (gamma-1.d0)*(rhoer(i)-rhotkr(i)-rhor(i)*qr(i))
           end if
        end if
      end do
c
      else
c
c----------------------------------------------------------
c
      do i = 1,imax
        rhoinvl(i)  = 1.d0/rhol(i)
        rhoinvr(i)  = 1.d0/rhor(i)
        sqrtrhol(i) = dsqrt(rhol(i))
        sqrtrhor(i) = dsqrt(rhor(i))
        sqrtinv(i) = 1.d0 / (sqrtrhol(i) + sqrtrhor(i))
c
        ul(i) = rhoinvl(i) * rhoul(i)
        vl(i) = rhoinvl(i) * rhovl(i)
        wl(i) = rhoinvl(i) * rhowl(i)
        el(i) = rhoinvl(i) * rhoel(i)
c
        ur(i) = rhoinvr(i) * rhour(i)
        vr(i) = rhoinvr(i) * rhovr(i)
        wr(i) = rhoinvr(i) * rhowr(i)
        er(i) = rhoinvr(i) * rhoer(i)
c
        ql(i) = 0.5d0*( ul(i)**2+vl(i)**2+wl(i)**2 )
        qr(i) = 0.5d0*( ur(i)**2+vr(i)**2+wr(i)**2 )
c
        if ( dabs(ronum).gt.1e-9 ) then
          psil = qll(i,2)*qll(i,2)+qll(i,3)*qll(i,3)
          psir = qrr(i,2)*qrr(i,2)+qrr(i,3)*qrr(i,3)
          thetal = ql(i)-0.5d0*psil
          thetar = qr(i)-0.5d0*psir
          pl(i) = (gamma-1.d0)*( rhoel(i) - rhol(i)*thetal )
          pr(i) = (gamma-1.d0)*( rhoer(i) - rhor(i)*thetar )
        else
          pl(i) = (gamma-1.d0)*(rhoel(i) - rhol(i)*ql(i))
          pr(i) = (gamma-1.d0)*(rhoer(i) - rhor(i)*qr(i))
        end if
c
        if(nl.eq.6) then
           tkl(i) = rhoinvl(i) * rhotkl(i)
           tkr(i) = rhoinvr(i) * rhotkr(i)
        else if(nl.eq.7) then
           tkl(i) = rhoinvl(i) * rhotkl(i)
           twl(i) = rhoinvl(i) * rhotwl(i)

           tkr(i) = rhoinvr(i) * rhotkr(i)
           twr(i) = rhoinvr(i) * rhotwr(i)

           if(ke) then
              pl(i) = (gamma-1.d0)*(rhoel(i)-rhotkl(i)-rhol(i)*ql(i))
              pr(i) = (gamma-1.d0)*(rhoer(i)-rhotkr(i)-rhor(i)*qr(i))
           end if
        end if
c
        tl(i)=gamma*machinf**2*pl(i)/rhol(i)
        tr(i)=gamma*machinf**2*pr(i)/rhor(i)
c
        tenthalpyl(i) = rhoinvl(i)*( rhoel(i)+pl(i) )
        tenthalpyr(i) = rhoinvr(i)*( rhoer(i)+pr(i) )
      end do

      end if
c
c-----------end of precondtion--------------------------------------
c
      select case (index)
      case (1)
        do i = 1,imax
          capul(i) = ul(i)*lx(i) + vl(i)*ly(i) + wl(i)*lz(i)
          capur(i) = ur(i)*lx(i) + vr(i)*ly(i) + wr(i)*lz(i)
          if(moving.ge.1) then
            qtl(i) = qll(i,1)*lx(i) + qll(i,2)*ly(i) + qll(i,3)*lz(i)
            qtr(i) = qrr(i,1)*lx(i) + qrr(i,2)*ly(i) + qrr(i,3)*lz(i)
          end if
        end do
      case (2)
        do i = 1,imax
          capvl(i) = ul(i)*mx(i) + vl(i)*my(i) + wl(i)*mz(i)
          capvr(i) = ur(i)*mx(i) + vr(i)*my(i) + wr(i)*mz(i)
          if(moving.ge.1) then
            qtl(i) = qll(i,1)*mx(i) + qll(i,2)*my(i) + qll(i,3)*mz(i)
            qtr(i) = qrr(i,1)*mx(i) + qrr(i,2)*my(i) + qrr(i,3)*mz(i)
          end if
        end do
      case (3)
        do i = 1,imax
          capwl(i) = ul(i)*nx(i) + vl(i)*ny(i) + wl(i)*nz(i)
          capwr(i) = ur(i)*nx(i) + vr(i)*ny(i) + wr(i)*nz(i)
          if(moving.ge.1) then
            qtl(i) = qll(i,1)*nx(i) + qll(i,2)*ny(i) + qll(i,3)*nz(i)
            qtr(i) = qrr(i,1)*nx(i) + qrr(i,2)*ny(i) + qrr(i,3)*nz(i)
          end if
        end do
      end select
c
      return
      end
