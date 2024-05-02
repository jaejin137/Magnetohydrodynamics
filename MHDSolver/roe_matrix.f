      subroutine roe_matrix(index,imax,nl,dim2,delta,sqrtinv,
     $     sqrtrhol, sqrtrhor, ul, ur, vl, vr, wl, wr, tenthalpyl,
     $     tenthalpyr, lx, ly, lz, lxhat, lyhat, lzhat, mx, my, mz, nx,
     $     ny, nz, mxhat, myhat, mzhat, nxhat, nyhat, nzhat, tkl, tkr,
     $     twl, twr, a, ke, qtl, qtr, ilower, qt1d, control)
c     
c     compute Roe matrix
c
c       index    imax        function
c      
c         1       ilp     Roe matrix for E
c         2       jlp     Roe matrix for F
c         3       klp     Roe matrix for G
c
c     Roe matrix is stored in a11( ), etc
c
c     IMPLICIT STATEMENT
c     
      use datain

      implicit none

      type (datain_type), intent(in)::control
c
c     DEFINITION OF LOCAL VARIABLES
c
c       NOTE:  The underscore "_" indicates a Roe-averaged or
c              Roe-defined variable
c
c       beta(i)         (gamma-1)/(2*c_(i)**2)
c
c       c_(i)           speed of sound
c
c       capc_(i)        c_(i) * sqrt(lx**2 + ly**2 + lz**2) for index = 1
c                       c_(i) * sqrt(mx**2 + my**2 + mz**2) for index = 2
c                       c_(i) * sqrt(nx**2 + ny**2 + nz**2) for index = 3
c
c       caph_(i)        total enthalpy
c
c       capu_(i)        lx(i)*u_(i) + ly(i)*v_(i) + lz(i)*w_(i)
c       capv_(i)        mx(i)*u_(i) + my(i)*v_(i) + mz(i)*w_(i)
c       capw_(i)        nx(i)*u_(i) + ny(i)*v_(i) + nz(i)*w_(i)
c
c       capuhat_(i)     lxhat(i)*u_(i) + lyhat(i)*v_(i) + lzhat(i)*w_(i)
c       capvhat_(i)     mxhat(i)*u_(i) + myhat(i)*v_(i) + mzhat(i)*w_(i)
c       capwhat_(i)     nxhat(i)*u_(i) + nyhat(i)*v_(i) + nzhat(i)*w_(i)
c
c       cinv_(i)        1./c_(i)
c
c
c       g11(i), ...     elements of (Lambda_tilda * L_tilda) where
c                       Lambda_tilda is diagonal matrix containing the
c                       absolute value of the Roe eigenvalues, and
c                       L_tilda is the matrix of left eigenvectors
c
c       h_(i)           static enthalpy
c
c       lambda_(i,n)    absolute value of n-th eigenvalue
c
c       llxhat(i)       components of unit vector 
c       llyhat(i)       whose actual definition 
c       llzhat(i)       depends on the value of "index"
c
c       mmxhat(i)       components of unit vector 
c       mmyhat(i)       whose actual definition 
c       mmzhat(i)       depends on the value of "index"
c
c       nnxhat(i)       components of unit vector 
c       nnyhat(i)       whose actual definition 
c       nnzhat(i)       depends on the value of "index"
c
c       q_(i)           0.5 * (u_(i)**2 + v_(i)**2 + w_(i)**2)
c
c       r11(i), ...     elements of R_tilda matrix (right eigenvectors)
c
c
c       sqrtrholinv(i)  1./sqrt(rhol)
c       sqrtrhorinv(i)  1./sqrt(rhor)
c
c       uuhat_(i)       inner product of unit normal and
c       vvhat_(i)       Roe-averaged velocity vector where
c       wwhat_(i)       the definition depends on the value of "index"
c
c       u_(i)           velocity in x-direction
c
c       v_(i)           velocity in y-direction
c
c       w_(i)           velocity in z-direction
c
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     index,
     $     imax,
     $     nl, dim2,
     $     ilower
      integer::
     $     moving
c     
      double precision, intent(in)::
     $     delta(5)
      double precision::
     $     gamma, ronum
c     
      double precision, dimension(imax), intent(in)::
     $     sqrtinv, sqrtrhol, sqrtrhor,
     $     ul, ur, vl, vr, wl, wr, tenthalpyl, tenthalpyr,
     $     qtl, qtr
c     
      double precision, dimension(ilower+1:imax-ilower), intent(in)::
     >     lx,    ly,    lz,
     >     lxhat, lyhat, lzhat, 
     >     mx,    my,    mz,
     >     mxhat, myhat, mzhat, 
     >     nx,    ny,    nz,
     >     nxhat, nyhat, nzhat
c     
      double precision, dimension(nl,nl,imax), intent(out)::
     $     a
c
c ... komega model
c
      double precision, dimension(imax), intent(in)::
     $     tkl, tkr,
     $     twl, twr
c
      logical, intent(in)::
     $     ke
c
      double precision, intent(in)::
     $     qt1d(ilower:dim2,3)
c     
c     LOCAL VARIABLES
c     
      integer::
     $     i, n
c     
      double precision::
     $     deltainv2, deltasq, c1
c     
      double precision, dimension(imax)::
     $     beta,     c_,        capc_,
     >     caph_,    capu_,     capv_,  
     >     capw_,    capuhat_,  capvhat_, 
     >     capwhat_, cinv_, 
     >     h_,       lambda_(imax,nl), 
     >     llxhat,   llyhat,    llzhat,
     >     mmxhat,   mmyhat,    mmzhat,   
     >     nnxhat,   nnyhat,    nnzhat,
     >     q_,       sqrtrholinv,  
     >     sqrtrhorinv,               u_, 
     >     uuhat_,   v_,        vvhat_,
     >     w_,       wwhat_,
     $     qt_
      double precision, dimension(nl,nl,imax)::
     $     r, g
c
      double precision, dimension(imax)::
     $     tk_, tw_
c
      double precision::
     $     ve, we, psi, phi, theta, gamma1
c     
c     --- SUBROUTINE START ---
c     
c---------------------------------------------
      moving=control%moving
      gamma=control%gamma
      ronum=control%ronum
c---------------------------------------------
c
      gamma1 = gamma - 1.d0
c
      do i = 1,imax
         u_(i)    = sqrtinv(i) * (sqrtrhol(i)*ul(i) + sqrtrhor(i)*ur(i))
         v_(i)    = sqrtinv(i) * (sqrtrhol(i)*vl(i) + sqrtrhor(i)*vr(i))
         w_(i)    = sqrtinv(i) * (sqrtrhol(i)*wl(i) + sqrtrhor(i)*wr(i))
         if(moving.ge.1) qt_(i)
     $        = sqrtinv(i) * (sqrtrhol(i)*qtl(i)+ sqrtrhor(i)*qtr(i))
         caph_(i) = sqrtinv(i) * (sqrtrhol(i)*tenthalpyl(i) + 
     >        sqrtrhor(i)*tenthalpyr(i))
         if(nl.eq.6) then
            tk_(i) = sqrtinv(i)*(sqrtrhol(i)*tkl(i)+sqrtrhor(i)*tkr(i))
         else if(nl.eq.7) then
            tk_(i) = sqrtinv(i)*(sqrtrhol(i)*tkl(i)+sqrtrhor(i)*tkr(i))
            tw_(i) = sqrtinv(i)*(sqrtrhol(i)*twl(i)+sqrtrhor(i)*twr(i))
         end if
      end do

      do i = 1,imax
         sqrtrholinv(i) = 1.d0 / sqrtrhol(i)
         sqrtrhorinv(i) = 1.d0 / sqrtrhor(i)
      end do
c     
      do i = 1,imax
         capu_(i)    = lx(i)*u_(i)    + ly(i)*v_(i)    + lz(i)*w_(i)
         capv_(i)    = mx(i)*u_(i)    + my(i)*v_(i)    + mz(i)*w_(i)
         capw_(i)    = nx(i)*u_(i)    + ny(i)*v_(i)    + nz(i)*w_(i)
         capuhat_(i) = lxhat(i)*u_(i) + lyhat(i)*v_(i) + lzhat(i)*w_(i)
         capvhat_(i) = mxhat(i)*u_(i) + myhat(i)*v_(i) + mzhat(i)*w_(i)
         capwhat_(i) = nxhat(i)*u_(i) + nyhat(i)*v_(i) + nzhat(i)*w_(i)
      end do
c     
      do i = 1,imax
         q_(i) = 0.5d0 * (u_(i)**2 + v_(i)**2 + w_(i)**2)
         if(nl.eq.7.and.ke) then
            c_(i) = dsqrt( gamma1*(caph_(i)-q_(i)-tk_(i)) )
         else
           if ( dabs(ronum).gt.1e-9 ) then
             ve = qt1d(i,2)+qt1d(i-1,2)
             we = qt1d(i,3)+qt1d(i-1,3)
             psi = ve*ve+we*we
             theta = q_(i)-0.125d0*psi
             c_(i) = dsqrt( gamma1*(caph_(i)-theta) )
           else
             c_(i) = dsqrt( gamma1*(caph_(i) - q_(i)) )
           end if
         end if
         cinv_(i) = 1.d0 / c_(i)
         beta(i) = 0.5d0 * gamma1 / (c_(i)**2)
         h_(i)   = c_(i)*c_(i)/gamma1
      end do
c     
      SELECT CASE (index)

      CASE (1)

         do i = 1,imax
            capc_(i) = c_(i) * dsqrt(lx(i)**2 + ly(i)**2 + lz(i)**2)
         end do
         do i = 1, imax
            if(moving.eq.0) then
               lambda_(i,1) = dabs(capu_(i) + capc_(i))
               lambda_(i,2) = dabs(capu_(i) - capc_(i))
               lambda_(i,3:nl) = dabs(capu_(i))
            else
               lambda_(i,1) = dabs(capu_(i) + capc_(i)+qt_(i))
               lambda_(i,2) = dabs(capu_(i) - capc_(i)+qt_(i))
               lambda_(i,3:nl) = dabs(capu_(i)+qt_(i))
            end if
         end do
c     
c     limit eigenvalues
c     
         do n = 1, 5
            if  (delta(n).le.0.d0) cycle
            deltainv2 = 0.5d0/delta(n)
            deltasq   = delta(n)**2
            do i = 1,imax
               if  ( lambda_(i,n).lt.delta(n) )  then
                  lambda_(i,n) = deltainv2*(lambda_(i,n)**2+deltasq)
               end if
            end do
         end do
c     
         do i = 1,imax
            uuhat_(i) = capuhat_(i)
            vvhat_(i) = capvhat_(i)
            wwhat_(i) = capwhat_(i)
            llxhat(i) = lxhat(i)
            llyhat(i) = lyhat(i)
            llzhat(i) = lzhat(i)
            mmxhat(i) = mxhat(i)
            mmyhat(i) = myhat(i)
            mmzhat(i) = mzhat(i)
            nnxhat(i) = nxhat(i)
            nnyhat(i) = nyhat(i)
            nnzhat(i) = nzhat(i)
         end do
c     
      CASE (2)

         do i = 1,imax
            capc_(i) = c_(i) * dsqrt(mx(i)**2 + my(i)**2 + mz(i)**2)
         end do
         do i = 1, imax
            if(moving.eq.0) then
               lambda_(i,1) = dabs(capv_(i) + capc_(i))
               lambda_(i,2) = dabs(capv_(i) - capc_(i))
               lambda_(i,3:nl) = dabs(capv_(i))
            else
               lambda_(i,1) = dabs(capv_(i) + capc_(i)+qt_(i))
               lambda_(i,2) = dabs(capv_(i) - capc_(i)+qt_(i))
               lambda_(i,3:nl) = dabs(capv_(i)+qt_(i))
            end if
         end do
c     
c     limit eigenvalues
c     
         do n = 1, 5
            if  (delta(n).le.0.d0) cycle
            deltainv2 = 0.5d0/delta(n)
            deltasq   = delta(n)**2
            do i = 1,imax
               if  ( lambda_(i,n).lt.delta(n) )  then
                  lambda_(i,n) = deltainv2*(lambda_(i,n)**2+deltasq)
               end if
            end do
         end do
c     
         do i = 1,imax
            uuhat_(i) = capvhat_(i)
            vvhat_(i) = capwhat_(i)
            wwhat_(i) = capuhat_(i)
            llxhat(i) = mxhat(i)
            llyhat(i) = myhat(i)
            llzhat(i) = mzhat(i)
            mmxhat(i) = nxhat(i)
            mmyhat(i) = nyhat(i)
            mmzhat(i) = nzhat(i)
            nnxhat(i) = lxhat(i)
            nnyhat(i) = lyhat(i)
            nnzhat(i) = lzhat(i)
         end do
c     
         CASE (3)

         do i = 1,imax
            capc_(i) = c_(i) * dsqrt(nx(i)**2 + ny(i)**2 + nz(i)**2)
         end do

         do i = 1, imax
            if(moving.eq.0) then
               lambda_(i,1) = dabs(capw_(i) + capc_(i))
               lambda_(i,2) = dabs(capw_(i) - capc_(i))
               lambda_(i,3:nl) = dabs(capw_(i))
            else
               lambda_(i,1) = dabs(capw_(i) + capc_(i)+qt_(i))
               lambda_(i,2) = dabs(capw_(i) - capc_(i)+qt_(i))
               lambda_(i,3:nl) = dabs(capw_(i)+qt_(i))
            end if
         end do
c     
c     limit eigenvalues
c     
         do n = 1, 5
            if  (delta(n).le.0.d0) cycle
            deltainv2 = 0.5d0/delta(n)
            deltasq   = delta(n)**2
            do i = 1,imax
               if  ( lambda_(i,n).lt.delta(n) )  then
                  lambda_(i,n) = deltainv2*(lambda_(i,n)**2+deltasq)
               end if
            end do
         end do
c     
         do i = 1,imax
            uuhat_(i) = capwhat_(i)
            vvhat_(i) = capuhat_(i)
            wwhat_(i) = capvhat_(i)
            llxhat(i) = nxhat(i)
            llyhat(i) = nyhat(i)
            llzhat(i) = nzhat(i)
            mmxhat(i) = lxhat(i)
            mmyhat(i) = lyhat(i)
            mmzhat(i) = lzhat(i)
            nnxhat(i) = mxhat(i)
            nnyhat(i) = myhat(i)
            nnzhat(i) = mzhat(i)
         end do
c     
      CASE DEFAULT

         write(6,"('*** index = ',i4,' in roe_matrix.  Stop.')") index
         stop

      END SELECT
c     
c     
      c1 = 1.d0 / gamma1
c     
c     MATRIX R
c     

      r = 0.d0

      do i = 1,imax
         r(1,1,i) = beta(i)
         r(1,2,i) = beta(i)
         r(1,5,i) = -2.d0*beta(i)
         if(nl.eq.7.and.ke) r(1,6,i) = 2.d0*beta(i)

         r(2,1,i) = beta(i)*(u_(i) + c_(i)*llxhat(i))
         r(2,2,i) = beta(i)*(u_(i) - c_(i)*llxhat(i))
         r(2,3,i) = mmxhat(i)
         r(2,4,i) = nnxhat(i)
         r(2,5,i) = -2.d0*beta(i)*u_(i)
         if(nl.eq.7.and.ke) r(2,6,i) = 2.d0*beta(i)*u_(i)

         r(3,1,i) = beta(i)*(v_(i) + c_(i)*llyhat(i))
         r(3,2,i) = beta(i)*(v_(i) - c_(i)*llyhat(i))
         r(3,3,i) = mmyhat(i)
         r(3,4,i) = nnyhat(i)
         r(3,5,i) = -2.d0*beta(i)*v_(i)
         if(nl.eq.7.and.ke) r(3,6,i) = 2.d0*beta(i)*v_(i)

         r(4,1,i) = beta(i)*(w_(i) + c_(i)*llzhat(i))
         r(4,2,i) = beta(i)*(w_(i) - c_(i)*llzhat(i))
         r(4,3,i) = mmzhat(i)
         r(4,4,i) = nnzhat(i)
         r(4,5,i) = -2.d0*beta(i)*w_(i)
         if(nl.eq.7.and.ke) r(4,6,i) = 2.d0*beta(i)*w_(i)

         r(5,1,i) = beta(i)*(caph_(i) + c_(i)*uuhat_(i))
         r(5,2,i) = beta(i)*(caph_(i) - c_(i)*uuhat_(i))
         r(5,3,i) = vvhat_(i)
         r(5,4,i) = wwhat_(i)
         r(5,5,i) = 1.d0 - 2.d0*beta(i)*caph_(i)
         if(nl.eq.7.and.ke) r(5,6,i) = 2.d0*beta(i)*caph_(i)

         if(nl.eq.6) then
            r(6,1,i) = beta(i)*tk_(i)
            r(6,2,i) = beta(i)*tk_(i)
            r(6,5,i) =-tk_(i)/h_(i)
            r(6,6,i) = 1.d0
         else if(nl.eq.7) then
            r(6,1,i) = beta(i)*tk_(i)
            r(6,2,i) = beta(i)*tk_(i)
            r(6,5,i) =-2.d0*beta(i)*tk_(i)
            if(ke) then
               r(6,6,i) = 2.d0*beta(i)*(tk_(i)+h_(i))
            else
               r(6,6,i) = 1.d0
            end if
            r(7,1,i) = beta(i)*tw_(i)
            r(7,2,i) = beta(i)*tw_(i)
            r(7,5,i) =-2.d0*beta(i)*tw_(i)
            if(ke) r(7,6,i) = 2.d0*beta(i)*tw_(i)
            r(7,7,i) = 1.d0
         end if
      end do

c     
c     MATRIX LAMBDA * L
c     
      g = 0.d0
c
      do i = 1,imax
         if ( dabs(ronum).gt.1e-9 ) then
           ve = qt1d(i,2)+qt1d(i-1,2)
           we = qt1d(i,3)+qt1d(i-1,3)
           psi = ve*ve+we*we
           phi = q_(i)+0.125d0*psi
         else
           phi = q_(i)
         end if
c
         g(1,1,i) = lambda_(i,1)*(phi - c1*c_(i)*uuhat_(i))
         g(1,2,i) = lambda_(i,1)*(c1*c_(i)*llxhat(i) - u_(i))
         g(1,3,i) = lambda_(i,1)*(c1*c_(i)*llyhat(i) - v_(i))
         g(1,4,i) = lambda_(i,1)*(c1*c_(i)*llzhat(i) - w_(i))
         g(1,5,i) = lambda_(i,1)
         if(nl.eq.7.and.ke) g(1,6,i) = -lambda_(i,1)
c
         g(2,1,i) =  lambda_(i,2)*(phi + c1*c_(i)*uuhat_(i))
         g(2,2,i) = -lambda_(i,2)*(u_(i) + c1*c_(i)*llxhat(i))
         g(2,3,i) = -lambda_(i,2)*(v_(i) + c1*c_(i)*llyhat(i))
         g(2,4,i) = -lambda_(i,2)*(w_(i) + c1*c_(i)*llzhat(i))
         g(2,5,i) =  lambda_(i,2)
         if(nl.eq.7.and.ke) g(2,6,i) = -lambda_(i,2)
c
         g(3,1,i) = -lambda_(i,3)*vvhat_(i)
         g(3,2,i) =  lambda_(i,3)*mmxhat(i)
         g(3,3,i) =  lambda_(i,3)*mmyhat(i)
         g(3,4,i) =  lambda_(i,3)*mmzhat(i)
c
         g(4,1,i) = -lambda_(i,4)*wwhat_(i)
         g(4,2,i) =  lambda_(i,4)*nnxhat(i)
         g(4,3,i) =  lambda_(i,4)*nnyhat(i)
         g(4,4,i) =  lambda_(i,4)*nnzhat(i)
c
         if(nl.eq.7.and.ke) then
            g(5,1,i) =  lambda_(i,5)*(q_(i) - h_(i) - tk_(i))
         else
            g(5,1,i) =  lambda_(i,5)*(phi - h_(i))
         end if
         g(5,2,i) = -lambda_(i,5)*u_(i)
         g(5,3,i) = -lambda_(i,5)*v_(i)
         g(5,4,i) = -lambda_(i,5)*w_(i)
         g(5,5,i) =  lambda_(i,5)
c
         if(nl.eq.6) then
            g(6,1,i) = -lambda_(i,6) * tk_(i)
            g(6,6,i) =  lambda_(i,6)
         else if(nl.eq.7) then
            g(6,1,i) = -lambda_(i,6) * tk_(i)
            g(6,6,i) =  lambda_(i,6)
            g(7,1,i) = -lambda_(i,7) * tw_(i)
            g(7,7,i) =  lambda_(i,7)
         end if
      end do
c     
c     ROE MATRIX
c     
      do i = 1, imax
         a(:,:,i) = matmul(r(:,:,i), g(:,:,i))
      end do
c     
      return
      end
