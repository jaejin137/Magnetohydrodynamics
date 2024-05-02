      subroutine zha_matrix(index,imax, dim2, nl,
     $     delta, gamma, ul, ur, vl, vr, wl, wr,rhol, rhor, 
     $     pl, pr,rhoel, rhoer,ql,qr,capul, capur, capvl, capvr, capwl, 
     $     capwr, lx, ly, lz,  mx, my, mz, nx, ny, nz,
     >     al,ar,dhat_l,dhat_r, ilower)
c
c     compute Zha matrix, dhat_l and dhat_r
c
c       index    imax        function
c      
c         1       ilp     matrix for E
c         2       jlp     matrix for F
c         3       klp     matrix for G
c
c     The matrix are stored in dhat_l,dhat_r

c     !! NOTE: When call this subroutine, there are 3 dummy variables need to
c              be set in the upper level subroutine according to the different direction, 
c              imax, al, ar
c              in xi-direction,   imax, al, ar = imax, al, ar
c              in eta-direction,  imax, al, ar = jmax, bl, br
c              in zeta-direction, imax, al, ar = kmax, cl, cr


c
c     IMPLICIT STATEMENT
c
      implicit none
c
c     DEFINITION OF LOCAL VARIABLES
c
c     When index=1, ax, ay, az =lx, ly, lz; capu_l=capul,capu_r=capur
c     When index=2, ax, ay, az =mx, my, mz; capu_l=capvl,capu_r=capvr
c     When index=3, ax, ay, az =nx, ny, nz; capu_l=capwl,capu_r=capwr
c
c       q(i)           0.5 * (u(i)**2 + v(i)**2 + w(i)**2)
c
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     index, imax, dim2, nl, ilower
c
      double precision, intent(in)::
     $     delta(5), gamma
c
      double precision, dimension(dim2), intent(in)::
     $     ul, ur, vl, vr, wl, wr,rhol, rhor,
     $     pl, pr,rhoel, rhoer, ql, qr,
     >     capul, capur, capvl, capvr, 
     >     capwl,capwr
c
      double precision, dimension(ilower+1:dim2), intent(in)::
     $     lx, ly, lz,
     $     mx, my, mz,
     $     nx, ny, nz

      double precision, dimension(imax, nl, nl),intent(in)::
     $     al, ar
c
      double precision, dimension(imax, nl, nl),intent(out)::
     $     dhat_l, dhat_r



c     LOCAL Variables

      double precision        ! temporary variables
     > temp0_l,temp1_l,temp3_l,temp4_l,temp5_l,temp6_l,temp7_l,
     > temp0_r,temp1_r,temp3_r,temp4_r,temp5_r,temp6_r,temp7_r,
     > beta, area,f

      integer i, n1,n2

      double precision 
     >     c_l,      ! left speed of sound
     >     c_r,      ! right speed of sound
     >     c_interface, ! interface speed of sound
     >     mach_l,   ! left normal Mach;
     >     mach_r    ! right normal Mach;

      double precision ax(imax),ay(imax),az(imax),
     >                 capu_l(imax),capu_r(imax)

c
c
c     the metrics were defined in subroutine metric
c     
c
c *** SUBROUTINE START ***
c

      beta = gamma*(gamma-1.) 

c..   Initialize the matrices
      do i = 1,imax
         do n1=1,nl
            do n2=1,nl
               dhat_l(i,n1,n2)=0.0
               dhat_r(i,n1,n2)=0.0
            end do
         end do
      end do



      if(index.eq.1) then          ! xi-direction

         do i = 1,imax
            ax(i)=lx(i)
            ay(i)=ly(i)
            az(i)=lz(i)
            capu_l(i)=capul(i)
            capu_r(i)=capur(i)
         end do

      else  if(index.eq.2) then    ! eta-direction 

         do i = 1,imax
            ax(i)=mx(i)
            ay(i)=my(i)
            az(i)=mz(i)
            capu_l(i)=capvl(i)
            capu_r(i)=capvr(i)
         end do

      else  if(index.eq.3) then    ! zeta-direction 

         do i = 1,imax
            ax(i)=nx(i)
            ay(i)=ny(i)
            az(i)=nz(i)
            capu_l(i)=capwl(i)
            capu_r(i)=capwr(i)
         end do
      
      end if


          DO   i = 1,imax

c..   calculate the interface  variables

            c_l = sqrt(beta*( rhoel(i)/rhol(i) -ql(i) ))
            c_r = sqrt(beta*( rhoer(i)/rhor(i) -qr(i) ))

            area=sqrt(ax(i)**2+ay(i)**2+az(i)**2)
            c_interface = 0.5*(c_l + c_r)*area

            mach_l = capu_l(i)/c_interface
            mach_r = capu_r(i)/c_interface

            IF(abs(mach_l).lt.1.0) Then ! subsonic

            IF(mach_l.ge.0.0) f=1.0
            IF(mach_l.lt.0.0) f=-1.0

c..         calculate the temporary variables 

            temp0_l=-rhoel(i)/rhol(i)**2 + 2.*ql(i)/rhol(i)
            temp0_r=-rhoer(i)/rhor(i)**2 + 2.*qr(i)/rhor(i)

            temp1_l =mach_l*((gamma-1.)*ql(i) -pl(i)/rhol(i)
     >           - area*pl(i)*beta/(4.*c_l*c_interface)*temp0_l)
            temp1_r =mach_r*((gamma-1.)*qr(i) -pr(i)/rhor(i)
     >           - area*pr(i)*beta/(4.*c_r*c_interface)*temp0_r)

            temp3_l = area*beta*capu_l(i)*pl(i)/(4.*rhol(i)*c_l
     >           *c_interface**2) 
            temp3_r = area*beta*capu_r(i)*pr(i)/(4.*rhor(i)*c_r
     >           *c_interface**2) 

            temp4_l = ((ax(i)*pl(i))/(rhol(i)*c_interface)
     >         + ul(i)*temp3_l - (gamma-1.)*ul(i)*mach_l )
            temp4_r = ((ax(i)*pr(i))/(rhor(i)*c_interface)
     >         + ur(i)*temp3_r - (gamma-1.)*ur(i)*mach_r) 

            temp5_l = ((ay(i)*pl(i))/(rhol(i)*c_interface)
     >         + vl(i)*temp3_l - (gamma-1.)*vl(i)*mach_l)
            temp5_r = ((ay(i)*pr(i))/(rhor(i)*c_interface)
     >         + vr(i)*temp3_r - (gamma-1.)*vr(i)*mach_r)

            temp6_l = ((az(i)*pl(i))/(rhol(i)*c_interface)
     >         + wl(i)*temp3_l - (gamma-1.)*wl(i)*mach_l)
            temp6_r =((az(i)*pr(i))/(rhor(i)*c_interface)
     >         + wr(i)*temp3_r - (gamma-1.)*wr(i)*mach_r) 

            temp7_l = (gamma-1.)*mach_l*(1.-
     >          area*gamma*pl(i)/(4.*rhol(i)*c_l*c_interface))
            temp7_r = (gamma-1.)*mach_r*(1.-
     >          area*gamma*pr(i)/(4.*rhor(i)*c_r*c_interface)) 

c..   calculate left and right matrices

c..         row 1

            dhat_l(i,1,2) = f*ax(i)
            dhat_r(i,1,2) = f*ax(i)

            dhat_l(i,1,3) = f*ay(i)
            dhat_r(i,1,3) = f*ay(i)

            dhat_l(i,1,4) = f*az(i)
            dhat_r(i,1,4) = f*az(i)

c..         row 2

            dhat_l(i,2,1) = ax(i)*temp1_l - f*(ul(i)*capu_l(i))
            dhat_r(i,2,1) = ax(i)*temp1_r - f*(ur(i)*capu_r(i))

            dhat_l(i,2,2) = ax(i)*temp4_l + f*(capu_l(i)+ul(i)*ax(i))
            dhat_r(i,2,2) = ax(i)*temp4_r + f*(capu_r(i)+ur(i)*ax(i))

            dhat_l(i,2,3) = ax(i)*temp5_l + f*(ul(i)*ay(i))
            dhat_r(i,2,3) = ax(i)*temp5_r + f*(ur(i)*ay(i))

            dhat_l(i,2,4) = ax(i)*temp6_l + f*(ul(i)*az(i))
            dhat_r(i,2,4) = ax(i)*temp6_r + f*(ur(i)*az(i))

            dhat_l(i,2,5) = ax(i)*temp7_l
            dhat_r(i,2,5) = ax(i)*temp7_r

c..         row 3

            dhat_l(i,3,1) = ay(i)*temp1_l - f*(vl(i)*capu_l(i))
            dhat_r(i,3,1) = ay(i)*temp1_r - f*(vr(i)*capu_r(i))

            dhat_l(i,3,2) = ay(i)*temp4_l + f*(vl(i)*ax(i))
            dhat_r(i,3,2) = ay(i)*temp4_r + f*(vr(i)*ax(i))

            dhat_l(i,3,3) = ay(i)*temp5_l + f*(capu_l(i)+vl(i)*ay(i))
            dhat_r(i,3,3) = ay(i)*temp5_r + f*(capu_r(i)+vr(i)*ay(i))

            dhat_l(i,3,4) = ay(i)*temp6_l + f*(vl(i)*az(i))
            dhat_r(i,3,4) = ay(i)*temp6_r + f*(vr(i)*az(i))

            dhat_l(i,3,5) = ay(i)*temp7_l
            dhat_r(i,3,5) = ay(i)*temp7_r

c..         row 4

            dhat_l(i,4,1) = az(i)*temp1_l - f*(wl(i)*capu_l(i))
            dhat_r(i,4,1) = az(i)*temp1_r - f*(wr(i)*capu_r(i))

            dhat_l(i,4,2) = az(i)*temp4_l + f*(wl(i)*ax(i))
            dhat_r(i,4,2) = az(i)*temp4_r + f*(wr(i)*ax(i))

            dhat_l(i,4,3) = az(i)*temp5_l + f*(wl(i)*ay(i))
            dhat_r(i,4,3) = az(i)*temp5_r + f*(wr(i)*ay(i))

            dhat_l(i,4,4) = az(i)*temp6_l + f*(capu_l(i)+wl(i)*az(i))
            dhat_r(i,4,4) = az(i)*temp6_r + f*(capu_r(i)+wr(i)*az(i))

            dhat_l(i,4,5) = az(i)*temp7_l
            dhat_r(i,4,5) = az(i)*temp7_r

c..         row 5

            temp7_l=c_interface + area*gamma*pl(i)/(4.*rhol(i)*c_l)
            temp7_r=c_interface + area*gamma*pr(i)/(4.*rhor(i)*c_r)

            dhat_l(i,5,1) = (gamma-1.)*ql(i)*c_interface
     >         + area*beta*pl(i)/(4.*c_l)*temp0_l 
     >           - f*(rhoel(i)*capu_l(i)/rhol(i))
            dhat_r(i,5,1) = (gamma-1.)*qr(i)*c_interface
     >         + area*beta*pr(i)/(4.*c_r)*temp0_r
     >           - f*(rhoer(i)*capu_r(i)/rhor(i))

            dhat_l(i,5,2) = -(gamma-1.)*ul(i)*temp7_l 
     >           + f*(rhoel(i)*ax(i)/rhol(i))
            dhat_r(i,5,2) = -(gamma-1.)*ur(i)*temp7_r
     >           + f*(rhoer(i)*ax(i)/rhor(i))

            dhat_l(i,5,3) = -(gamma-1.)*vl(i)*temp7_l
     >           + f*(rhoel(i)*ay(i)/rhol(i))
            dhat_r(i,5,3) = -(gamma-1.)*vr(i)*temp7_r
     >           + f*(rhoer(i)*ay(i)/rhor(i))

            dhat_l(i,5,4) = -(gamma-1.)*wl(i)*temp7_l
     >           + f*(rhoel(i)*az(i)/rhol(i))
            dhat_r(i,5,4) = -(gamma-1.)*wr(i)*temp7_r
     >           + f*(rhoer(i)*az(i)/rhor(i))

            dhat_l(i,5,5) = (gamma-1.)*temp7_l + f*capu_l(i)
            dhat_r(i,5,5) = (gamma-1.)*temp7_r + f*capu_r(i)


            ELSE IF(mach_l.ge.1.0) Then   ! supersonic

               do n1=1,nl
               do n2=1,nl
                  dhat_l(i,n1,n2) = al(i,n1,n2)
                  dhat_r(i,n1,n2) = ar(i,n1,n2)
               end do
               end do

            ELSE IF(mach_l.le.-1.0) Then    ! supersonic

               do n1=1,nl
               do n2=1,nl
                  dhat_l(i,n1,n2) = -al(i,n1,n2)
                  dhat_r(i,n1,n2) = -ar(i,n1,n2)
               end do
               end do

            END IF

           if(i.eq.7000)then
              open(21,file='zha_matrices')
              write(21,*)'dhat_l,mach_l=', mach_l
           do n1= 1,nl
              write(21,2156)(dhat_l(i,n1,n2),n2=1,nl)
           end do 
              write(21,*)'dhat_r,index=',index
           do n1= 1,nl
              write(21,2156)(dhat_r(i,n1,n2),n2=1,nl)
           end do 
           end if
 2156      format(5e20.12)


          END DO
c
c
      return
      end
