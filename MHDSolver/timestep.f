c======================================================================c
c
      subroutine timestep(dt, il, jl, kl, nl, x, y, z, q, dim2,
     $     vol, ilower, iupper, jlower, jupper, klower,
     $     kupper, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper,
     $     qiii, diver, qt, control)
c
c     compute local time step on each cell
c     
c     IMPLICIT STATEMENT
c     
      use datain
      implicit none
      type (datain_type)::control
c     
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     il, jl, kl,          ! cell number in 3 directions
     $     nl,                  ! equation number
     $     dim2,                ! maximum face number in 3D
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper
      integer::
     $     idimen,              ! dimension number
     $     unidt, precondition
cc
      integer, intent(in)::
     $     bc_xi_lower(jl,kl), bc_xi_upper(jl,kl),
     $     bc_eta_lower(il,kl), bc_eta_upper(il,kl),
     $     bc_zeta_lower(il,jl), bc_zeta_upper(il,jl)
c
      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z
c     
      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol
c
      double precision, intent(in),
     $     dimension(ilower:iupper,jlower:jupper,klower:kupper,3)::
     $     qt
c
      double precision, intent(in)::
     $     q(ilower:iupper, jlower:jupper, klower:kupper, nl),
     $     qiii(ilower:iupper, jlower:jupper, klower:kupper, nl),
     $     diver(ilower:iupper, jlower:jupper, klower:kupper)
c
      double precision::
     $     cfl,                 ! cfl number
     $     gamma, machinf ,     ! specific heat ratio
     $     ronum,
     $     k_prec
c     
      double precision, intent(out)::
     $     dt(il,jl,kl)         ! time marching interval
c     
C     LOCAL VARIABLES
c     
      integer::
     $     i, j, k,             ! cell iteration index
     $     index,               ! direction index
     $     bclower, bcupper
c
      double precision, dimension(ilower+1:dim2)::
     $     lxhat, lyhat, lzhat,
     $     mxhat, myhat, mzhat,
     $     nxhat, nyhat, nzhat,
     $     lx, ly, lz,
     $     mx, my, mz,
     $     nx, ny, nz
c
      double precision::
     $     capvel, sound, rho, rinv,
     $     u, v, w, p, t
c     
      double precision
     $     eigenmax(il,jl,kl)
c
      double precision::
     $     vel, urr, temp, capup, drhodp, drhodt, cp,
     $     betap, alpha, capcpp, gamma1, psi, theta, qq
c--------------------------parameter transfer

      idimen=control%idimen
      precondition=control%precondition
      unidt=control%unidt

      cfl=control%cfl
      gamma=control%gamma
      k_prec=control%k_prec
      machinf=control%machinf
      ronum=control%ronum

c     
c     --- SUBROUTINE START ---
c     
      eigenmax = 0.d0
      dt = 0.d0
      gamma1 = gamma - 1.d0
c
c *** xi direction ***
c
      index    = 1
c
      do k = 1,kl
        do j = 1,jl
          bclower = bc_xi_lower(j,k)
          bcupper = bc_xi_upper(j,k)
c
          call metric(index,j,k, il+1, il, jl, kl, x, y, z, lx, ly,
     $         lz, mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat, mxhat,
     $         myhat, mzhat, nxhat, nyhat, nzhat, ilower, iupper,
     $         jlower, jupper, klower, kupper, bclower, bcupper)
c           
          if (precondition.ge.1) then
            do i = 1,il
              if (precondition.ge.2) then
                u=qiii(i,j,k,2)
                v=qiii(i,j,k,3)
                w=qiii(i,j,k,4)
                p=qiii(i,j,k,1)
                t=qiii(i,j,k,5)
                call qstat(p,t,diver(i,j,k),temp)
                rho=1.0d0/temp                 
                call soundd(temp,t,sound)
                call cpt(t,cp)
                call drhdtpp(temp,t,drhodt)
                call drhdptt(temp,t,drhodp)
              else
                rho = q(i,j,k,1)
                rinv = 1.d0 / rho
                u = rinv * q(i,j,k,2)
                v = rinv * q(i,j,k,3)
                w = rinv * q(i,j,k,4)
                if (nl.eq.7) then
                  p = (gamma-1.d0)*(q(i,j,k,5)-q(i,j,k,6) -
     >                0.5d0*rho*(u**2+v**2+w**2))
                else
                  p = (gamma-1.d0)*(q(i,j,k,5) - 0.5d0*rho*
     >                (u**2+v**2+w**2))
                end if
                sound=dsqrt(gamma *p * rinv)
                temp=gamma*machinf**2*p/rho
                drhodp=gamma/sound**2
                drhodt=-rho/temp
                cp=1.0/((gamma-1.0)*machinf**2)
              end if
                vel=dsqrt(u**2+v**2+w**2)
                if (vel.gt.sound) then
                  urr=sound
                elseif (vel.lt.sound/100000.0) then
                  urr=sound/100000.0
                else
                  urr=vel
                end if
c               urr=sound
                urr=dmin1(sound,dmax1(vel,k_prec))
                capvel = dabs(lx(i)*u + ly(i)*v + lz(i)*w)
                betap=(drhodp+drhodt/(rho*cp))
                alpha=0.5*(1.0-betap*urr**2)
                capcpp=dsqrt((alpha*capvel)**2+
     $                 urr*urr*(lx(i)**2 + ly(i)**2 + lz(i)**2))
                capup=capvel*(1.0-alpha)
                eigenmax(i,j,k)=dmax1(capvel,capup+capcpp)/vol(i,j,k)
            end do
          else
            do i = 1,il
              rho = q(i,j,k,1)
              rinv = 1.d0 / rho
              u = rinv * q(i,j,k,2)
              v = rinv * q(i,j,k,3)
              w = rinv * q(i,j,k,4)
              qq = 0.5d0*( u*u+v*v+w*w )
              if (nl.eq.7) then
                p = gamma1*( q(i,j,k,5)-q(i,j,k,6) - rho*qq )
              else
                if ( dabs(ronum).gt.1e-9 ) then
                  psi = qt(i,j,k,2)*qt(i,j,k,2)+qt(i,j,k,3)*qt(i,j,k,3)
                  theta = qq-0.5d0*psi
                  p = gamma1*( q(i,j,k,5)-rho*theta )
                else
                  p = gamma1*( q(i,j,k,5)-rho*qq )
                end if
              end if
              sound = dsqrt( (lx(i)**2+ly(i)**2+lz(i)**2) * gamma *
     >                p * rinv )
              capvel = dabs( lx(i)*u + ly(i)*v + lz(i)*w )
              eigenmax(i,j,k) = dmax1(eigenmax(i,j,k),
     $                          ((capvel+sound)/vol(i,j,k)))
            end do
          end if
        end do
      end do
c
c *** eta direction ***
c
      if (idimen.ge.2) then
c
        index    = 2
c
        do k = 1,kl
          do i = 1,il
            bclower = bc_eta_lower(i,k)
            bcupper = bc_eta_upper(i,k)
c
            call metric(index,i,k, jl+1, il, jl, kl, x, y, z, lx, ly,
     $              lz, mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat,
     $              mxhat, myhat, mzhat, nxhat, nyhat, nzhat, ilower,
     $              iupper, jlower, jupper, klower, kupper, bclower,
     $              bcupper)
c     
            if (precondition.ge.1) then
              do j = 1,jl
                if (precondition.ge.2) then
                     u=qiii(i,j,k,2)
                     v=qiii(i,j,k,3)
                     w=qiii(i,j,k,4)
                     p=qiii(i,j,k,1)
                     t=qiii(i,j,k,5)
                     call qstat(p,t,diver(i,j,k),temp)
                     rho=1.0d0/temp                 
                     call soundd(temp,t,sound)
                     call cpt(t,cp)
                     call drhdtpp(temp,t,drhodt)
                     call drhdptt(temp,t,drhodp)
                else
                     rho = q(i,j,k,1)
                     rinv = 1.d0 / rho
                     u = rinv * q(i,j,k,2)
                     v = rinv * q(i,j,k,3)
                     w = rinv * q(i,j,k,4)
                     if(nl.eq.7) then
                       p = (gamma-1.d0)*(q(i,j,k,5)-q(i,j,k,6) -
     >                     0.5d0*rho*(u**2+v**2+w**2))
                     else
                       p = (gamma-1.d0)*(q(i,j,k,5) - 0.5d0*rho*
     >                    (u**2+v**2+w**2))
                     end if
                     sound=dsqrt(gamma *p * rinv)
                     temp=gamma*machinf**2*p/rho
                     drhodp=gamma/sound**2
                     drhodt=-rho/temp
                     cp=1.0/((gamma-1.0)*machinf**2)
                end if
                  vel=dsqrt(u**2+v**2+w**2)
                  if (vel.gt.sound) then
                     urr=sound
                  elseif (vel.lt.sound/100000.0) then
                     urr=sound/100000.0
                  else
                     urr=vel
                  end if
c                 urr=sound
                  urr=dmin1(sound,dmax1(vel,k_prec))
                  capvel = dabs(mx(j)*u + my(j)*v + mz(j)*w)
                  betap=(drhodp+drhodt/(rho*cp))
                  alpha=0.5*(1.0-betap*urr**2)
                  capcpp=dsqrt((alpha*capvel)**2+
     $              urr*urr*(mx(j)**2 + my(j)**2 + mz(j)**2))
                  capup=capvel*(1.0-alpha)
                  eigenmax(i,j,k)=dmax1(capvel,capup+capcpp)/vol(i,j,k)
              end do
            else
              do j = 1,jl
                rho = q(i,j,k,1)
                rinv = 1.d0 / rho
                u = rinv * q(i,j,k,2)
                v = rinv * q(i,j,k,3)
                w = rinv * q(i,j,k,4)
                qq = 0.5d0*( u*u+v*v+w*w )
                if (nl.eq.7) then
                  p = gamma1*( q(i,j,k,5)-q(i,j,k,6)-rho*qq )
                else
                  if ( dabs(ronum).gt.1e-9 ) then
                    psi = qt(i,j,k,2)*qt(i,j,k,2)+
     >                    qt(i,j,k,3)*qt(i,j,k,3)
                    theta = qq-0.5d0*psi
                    p = gamma1*( q(i,j,k,5)-rho*theta )
                  else
                    p = gamma1*( q(i,j,k,5)-rho*qq )
                  end if
                end if
                sound = dsqrt((mx(j)**2+my(j)**2+mz(j)**2) * gamma *
     >                  p * rinv)
                capvel = dabs(mx(j)*u + my(j)*v + mz(j)*w)
                eigenmax(i,j,k) = dmax1(eigenmax(i,j,k),
     $                            ((capvel+sound)/vol(i,j,k)))
              end do
            end if
          end do
        end do
      end if
c
c *** zeta direction ***
c
      if (idimen.eq.3) then
c
         index    = 3
c
         do j = 1,jl
           do i = 1,il
             bclower = bc_zeta_lower(i,j)
             bcupper = bc_zeta_upper(i,j)
c
             call metric(index,i,j, kl+1, il, jl, kl, x, y, z, lx, ly,
     $              lz, mx, my, mz, nx, ny, nz, lxhat, lyhat, lzhat,
     $              mxhat, myhat, mzhat, nxhat, nyhat, nzhat, ilower,
     $              iupper, jlower, jupper, klower, kupper, bclower,
     $              bcupper)
c     

             if (precondition.ge.1) then
               do k = 1,kl
                 if (precondition.ge.2) then
                     u=qiii(i,j,k,2)
                     v=qiii(i,j,k,3)
                     w=qiii(i,j,k,4)
                     p=qiii(i,j,k,1)
                     t=qiii(i,j,k,5)
                     call qstat(p,t,diver(i,j,k),temp)
                     rho=1.0d0/temp                 
                     call soundd(temp,t,sound)
                     call cpt(t,cp)
                     call drhdtpp(temp,t,drhodt)
                     call drhdptt(temp,t,drhodp)
                 else
                     rho = q(i,j,k,1)
                     rinv = 1.d0 / rho
                     u = rinv * q(i,j,k,2)
                     v = rinv * q(i,j,k,3)
                     w = rinv * q(i,j,k,4)
                     if(nl.eq.7) then
                       p = (gamma-1.d0)*(q(i,j,k,5)-q(i,j,k,6) -
     >                     0.5d0*rho*(u**2+v**2+w**2))
                     else
                       p = (gamma-1.d0)*(q(i,j,k,5) - 0.5d0*rho*
     >                    (u**2+v**2+w**2))
                     end if
                     sound=dsqrt(gamma *p * rinv)
                     temp=gamma*machinf**2*p/rho
                     drhodp=gamma/sound**2
                     drhodt=-rho/temp
                     cp=1.0/((gamma-1.0)*machinf**2)
                 end if
                 vel=dsqrt(u**2+v**2+w**2)
                 if (vel.gt.sound) then
                     urr=sound
                 elseif (vel.lt.sound/100000.0) then
                     urr=sound/100000.0
                 else
                     urr=vel
                 end if
c                urr=sound
                 urr=dmin1(sound,dmax1(vel,k_prec))
                 capvel = dabs(nx(k)*u + ny(k)*v + nz(k)*w)
                 betap=(drhodp+drhodt/(rho*cp))
                 alpha=0.5*(1.0-betap*urr**2)
                 capcpp=dsqrt((alpha*capvel)**2+
     $                  urr*urr*(nx(k)**2 + ny(k)**2 + nz(k)**2))
                 capup=capvel*(1.0-alpha)
                 eigenmax(i,j,k)=dmax1(capvel,capup+capcpp)/vol(i,j,k)
               end do
             else
               do k = 1,kl
                 rho = q(i,j,k,1)
                 rinv = 1.d0 / rho
                 u = rinv * q(i,j,k,2)
                 v = rinv * q(i,j,k,3)
                 w = rinv * q(i,j,k,4)
                 qq = 0.5d0*( u*u+v*v+w*w )
                 if (nl.eq.7) then
                   p = (gamma-1.d0)*(q(i,j,k,5)-q(i,j,k,6)-0.5d0*rho*
     >                 (u**2+v**2+w**2))
                 else
                   if ( dabs(ronum).gt.1e-9 ) then
                     psi = qt(i,j,k,2)*qt(i,j,k,2)+
     >                     qt(i,j,k,3)*qt(i,j,k,3)
                     theta = qq-0.5d0*psi
                     p = gamma1*( q(i,j,k,5)-rho*theta )
                   else
                     p = gamma1*( q(i,j,k,5)-rho*qq )
                   end if
                 end if
                 sound = dsqrt( (nx(k)**2+ny(k)**2+nz(k)**2) * gamma *
     >                   p * rinv )
                 capvel = dabs( nx(k)*u + ny(k)*v + nz(k)*w )
                 eigenmax(i,j,k) = dmax1(eigenmax(i,j,k),
     $                             ((capvel+sound)/vol(i,j,k)))
               end do
             end if
           end do
        end do
      end if
c
      if(unidt.gt.0) eigenmax = maxval(eigenmax) ! local time step switch
c
      dt = cfl / eigenmax
c     
      return
      end
c
c======================================================================c
