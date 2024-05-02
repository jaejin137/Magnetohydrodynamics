      subroutine initflow(rank, il, jl, kl, nl, ilower, iupper,
     $     jlower, jupper, klower, kupper, snum, 
     $     ntotal, time, 
     $     bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $     bc_zeta_lower, bc_zeta_upper, xsl, xs, ys, visturb,
     $     x, y, z, xo, yo, zo, q, qi, qii, qiii, control)
c
c     read rstart to initialize flowfield
c
c     IMPLICIT STATEMENT
c
      use datain

      implicit none

      type (datain_type), intent(in)::control
c
c     INTERFACE VARIABLES       ! specified in aamain.f 
c
      integer, intent(in)::
     $     rank,
     $     il, jl, kl,
     $     nl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     snum
      integer::
     $     strtp,
     $     moving,
     $     dual_t,
     $     integrate, precondition,
     $     turb, main_dir

      double precision::
     $     tintvl, gamma, machinf

      integer, intent(out)::
     $     ntotal

      double precision, intent(out)::
     $     time

      integer, intent(out)::
     $     bc_xi_lower(jl,kl), bc_xi_upper(jl,kl),
     $     bc_eta_lower(il,kl), bc_eta_upper(il,kl),
     $     bc_zeta_lower(il,jl), bc_zeta_upper(il,jl)

      double precision, intent(out)::
     $  xsl(4,3),
     $  xs(il+1,jl+1,snum),
     $  ys(il+1,jl+1,snum),
     $  visturb(ilower+1:iupper-1,jlower+1:jupper-1,klower+1:kupper-1)

      double precision, intent(out),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z,
     $     xo, yo, zo

      double precision, intent(out),
     $     dimension(ilower:iupper,jlower:jupper,klower:kupper,nl)::
     $     q, qi, qii, qiii
c
c     Local Variables
c
      character(len=20)::      !original
     $     file_rstart,        !original
     $     filename            !original

      logical::
     $     ex

      integer::
     $     iltp, jltp, kltp,
     $     i, j, k, n, blen,
     $     sts
      
      integer, parameter::
     $     file1 = 1,
     $     fileosc = 11

      double precision::
     $     p,t,r,th,v2

c--------------------------parameter transfer

      integrate=control%integrate
      main_dir=control%main_dir 
      moving=control%moving
      dual_t=control%dual_t
      turb=control%turb
      precondition=control%precondition
      strtp=control%strtp

      gamma=control%gamma
      machinf=control%machinf
      tintvl=control%tintvl
c--------------------------------------------------------

c
c *** begin
c      
c *** file name and existance check
c
      write(file_rstart, '("rstart",i3.3)') rank
      print *, "rstart file name is ",file_rstart
      inquire(file=file_rstart, exist=ex)
      if(ex) then
         open(file1, file=file_rstart, status='unknown', 
     >        access='sequential', form= 'unformatted')
      else
         write(*,*) '"',trim(file_rstart),'" not exist'

         stop
      end if

c..   read rstart

      if (rank.eq.1) write(*,*) 'Start reading rstart ...'

      read(file1)  ntotal, time, blen

      read(file1)  iltp, jltp, kltp

      do k = klower+1,kupper
         do j = jlower+1,jupper
            do i = ilower+1,iupper
               read(file1)  x(i,j,k), y(i,j,k), z(i,j,k)
            end do
         end do
      end do

      do n = 1, nl
         do k = klower, kupper
            do j = jlower, jupper
               do i = ilower, iupper
                  read(file1)  q(i,j,k,n)
               end do
            end do
         end do
      end do

      do k = 1,kl
         do j = 1,jl
            read(file1)  bc_xi_lower(j,k), bc_xi_upper(j,k)
         end do
      end do

      do k = 1,kl
         do i = 1,il
            read(file1)  bc_eta_lower(i,k), bc_eta_upper(i,k)
         end do
      end do

      do j = 1,jl
         do i = 1,il
            read(file1)  bc_zeta_lower(i,j), bc_zeta_upper(i,j)
         end do
      end do

      read(file1) ((xsl(i,j), i=1,4), j=1,3)

      close(file1)

c..   read oscillating mesh info

      if(strtp.eq.3.and.moving.eq.1) then

         write(*,*) 'reading oscillating mesh ..'

         write(filename, '("bin_2d_r",i2.2,".dat")') rank

         inquire(file=filename, exist=ex)

         if(ex) then
            write(*,*) 'oscillating grd on rank ', rank
            write(*,*) 'read ', filename
            open(fileosc, file=filename, form='unformatted')
            read(fileosc) iltp, jltp
            if(iltp.ne.il.or.jltp.ne.jl) then
               write(*,*) 'il, jl differ in datain and osci mesh ',
     $              filename
               write(*,*) 'datain : il=', il, ' jl=', jl
               write(*,*) "oscmesh: il=", iltp, ' jl=', jltp
               stop
            end if
            do n = 1, snum
               do j = 1, jl+1
                  do i = 1, il+1
                     read(fileosc, iostat=sts) xs(i,j,n), ys(i,j,n)
                  end do
               end do
            end do
            close(fileosc)
         else                   ! mesh unchanged for fixed passages
            write(*,*) 'steady grd      on rank ', rank
            do n = 1, snum
               xs(:,:,n) = x(1:il+1,1:jl+1,1)
               ys(:,:,n) = y(1:il+1,1:jl+1,1)
            end do
         end if
      end if
c
c *** read previous step for dual time stepping
c
      if(dual_t.eq.1) then      ! qi = q, qii = previous step
         qi = q
         if(moving.ge.1) then   ! previous mesh
            xo = x
            yo = y
            zo = z
         end if
         if(time.gt.tintvl) then ! resume computation
            write(*,*) 'read rstart at previous step ..'
            write(file_rstart, '("rstart", i3.3, "p")') rank
            write(*,*) file_rstart
            inquire(file=file_rstart, exist=ex)
            if(.not.ex) then
               write(*,*) 'previous step ', file_rstart, ' missing'
               write(*,*) 'initial time = 0, initial step = 1'
               time = 0.
            else
               open(file1,file=file_rstart,form='unformatted')
               read(file1)      ! ntotal, time
               read(file1)      ! il, jl, kl
               if(moving.eq.0) then
                 do k = klower+1,kupper
                   do j = jlower+1,jupper
                     do i = ilower+1,iupper
                       read(file1)
                     end do
                   end do
                 end do
               else
                 do k = klower+1,kupper
                   do j = jlower+1,jupper
                     do i = ilower+1,iupper
                       read(file1)  xo(i,j,k), yo(i,j,k), zo(i,j,k)
                     end do
                   end do
                 end do
               end if
               do n = 1, nl
                  do k = klower, kupper
                     do j = jlower, jupper
                        do i = ilower, iupper
                           read(file1)  qii(i,j,k,n)
                        end do
                     end do
                  end do
               end do
               close(file1)
            end if
         end if
      end if
c
      do k = klower+1,kupper-1
        do j = jlower+1,jupper-1
          do i = ilower+1,iupper-1
            visturb(i,j,k) = 0.d0
          end do
        end do
      end do
c
      if (precondition.ge.1) then
         if (precondition.eq.1) then
            do k = klower, kupper
               do j = jlower, jupper
                  do i = ilower, iupper
                     qiii(i,j,k,2)=q(i,j,k,2)/q(i,j,k,1)
                     qiii(i,j,k,3)=q(i,j,k,3)/q(i,j,k,1)
                     qiii(i,j,k,4)=q(i,j,k,4)/q(i,j,k,1)
                     v2=qiii(i,j,k,2)**2+qiii(i,j,k,3)**2+
     $                  qiii(i,j,k,4)**2
                     qiii(i,j,k,1)=(gamma-1.0)*
     $                  (q(i,j,k,5)-0.5*q(i,j,k,1)*v2)
                     qiii(i,j,k,5)=qiii(i,j,k,1)*gamma*
     $                          machinf*machinf/q(i,j,k,1)
                  end do
               end do
            end do
         else
            write(file_rstart, '("pre_rstart",i3.3)') rank
            inquire(file=file_rstart, exist=ex)
            if(ex) then
              open(file1, file=file_rstart, status='unknown', 
     >           access='sequential', form= 'unformatted')
            else
               write(*,*) '"',trim(file_rstart),'" not exist'
               stop
            end if
            do n = 1, nl
               do k = klower, kupper
                  do j = jlower, jupper
                     do i = ilower, iupper
                        read(file1)  qiii(i,j,k,n)
                     end do
                  end do
               end do
            end do
            close(file1)
c---------------------
               do k=1,kl
                 do j=1,jl
                  do i=1,il
                     p=qiii(i,j,k,1)
                     t=qiii(i,j,k,5)
                     call qstat(p,t,0.0,r)
                     q(i,j,k,1) = 1.0d0/r
                     q(i,j,k,2) = qiii(i,j,k,2)*q(i,j,k,1)
                     q(i,j,k,3) = qiii(i,j,k,3)*q(i,j,k,1)
                     q(i,j,k,4) = qiii(i,j,k,4)*q(i,j,k,1)
                     v2 = 0.5*(qiii(i,j,k,2)**2+qiii(i,j,k,3)**2+
     $                    qiii(i,j,k,4)**2)
                     call enthalpy(th,t,r,p)
                     q(i,j,k,5) =  q(i,j,k,1)*(th+v2)-p                   
                   end do
                 end do
               end do             
         end if   
      end if
c
      return
      end
