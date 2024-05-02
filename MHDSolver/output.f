      subroutine output(rank, dual_t, dstep, ntotal, il, jl, kl, nl,
     $     ilower, iupper, jlower, jupper, klower, kupper, x, y, z, xo,
     $     yo, zo, xsl, bc_xi_lower, bc_xi_upper, bc_eta_lower,
     $     bc_eta_upper, bc_zeta_lower, bc_zeta_upper, q, qii, time,
     $     tintvl, turb, visturb, blen,precondition,qiii)
c     
c     write flowfield file and tecplot file ( for 2d airfoil )
c     
c     IMPLICIT STATEMENT
c     
      implicit none
c     
c     INTERFACE VARIABLES
c     
      integer, intent(in)::
     $     rank,
     $     dual_t,
     $     dstep,
     $     ntotal,
     $     il, jl, kl,
     $     nl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     turb,precondition
c    
      double precision, intent(in)::
     $     time,
     $     tintvl
c
      double precision, intent(in), dimension(ilower+1:iupper,
     $     jlower+1:jupper, klower+1:kupper)::
     $     x, y, z, xo, yo, zo
c
      double precision, intent(in)::
     $ visturb(ilower+1:iupper-1,jlower+1:jupper-1, klower+1:kupper-1)
c     
      double precision, intent(in)::
     $     xsl(4,3)
c
      double precision, intent(in),
     $     dimension(ilower:iupper,jlower:jupper, klower:kupper, nl)::
     $     q, qii,qiii
c     
      integer, dimension(jl,kl), intent(in)::
     $     bc_xi_lower, bc_xi_upper
c     
      integer, dimension(il, kl), intent(in)::
     $     bc_eta_lower, bc_eta_upper
c     
      integer, dimension(il, jl), intent(in)::
     $     bc_zeta_lower, bc_zeta_upper
c     
c     LOCAL VARIABLES
c     
      integer::
     $     i, j, k, blen,       ! iteration index for cells
     $     n                    ! iteration index for variables
c     
      integer, parameter::
     $     file2 = 2,           ! raw data file
     $     file3 = 3           ! results at previous step
c
      character(len=20)::
     $     file_rstore,
     $     file_rstorep,
     $     pre_file_rstore
c
c     *** SUBROUTINE START ***
c     
c-----------------------
      if (precondition.ge.2) then
        write(pre_file_rstore, '("pre_rstore",i3.3)') rank
        open(file2, file=pre_file_rstore, form= 'unformatted')
        do n = 1, nl
          do k = klower,kupper
            do j = jlower, jupper
              do i = ilower,iupper
                write(file2) qiii(i,j,k,n)
              end do
            end do
          end do
        end do
        close(file2)
      end if
c-------------------------

      if(dual_t.eq.0) then
        write(file_rstore, '("rstore",i3.3)') rank
      else
        write(file_rstore, '("rstore",i3.3,"_",i7.7)')
     $    rank, dstep
        write(file_rstorep, '("rstore",i3.3,"_",i7.7)')
     $    rank, dstep-1
      end if
c
      open(file2, file=file_rstore, form= 'unformatted')
      write(file2) ntotal, time, blen
      write(file2) il, jl, kl
c
      do k = klower+1,kupper
        do j = jlower+1,jupper
          do i = ilower+1,iupper
            write(file2) x(i,j,k), y(i,j,k), z(i,j,k)
          end do
        end do
      end do
c
      do n = 1, nl
        do k = klower,kupper
          do j = jlower, jupper
            do i = ilower,iupper
              write(file2) q(i,j,k,n)
            end do
          end do
        end do
      end do
c
      do k = 1,kl
        do j = 1,jl
          write(file2) bc_xi_lower(j,k), bc_xi_upper(j,k)
        end do
      end do
c
      do k = 1,kl
        do i = 1,il
          write(file2) bc_eta_lower(i,k), bc_eta_upper(i,k)
        end do
      end do
c
      do j = 1,jl
        do i = 1,il
          write(file2) bc_zeta_lower(i,j), bc_zeta_upper(i,j)
        end do
      end do
c
      write(file2) ((xsl(i,j), i=1,4), j=1,3)
c
      if (turb.gt.0) then
        do k = klower+1, kupper-1
          do j = jlower+1, jupper-1
            do i = ilower+1, iupper-1
              write(file2) visturb(i,j,k)
            end do
          end do
        end do
      end if
c
      close(file2)
c
c *** unsteady
c
      if (dual_t.eq.1) then
        open(file3, file=file_rstorep, form='unformatted')
        write(file3) ntotal, time-tintvl, blen ! ntotal is not correct
        write(file3) il, jl, kl
c
        do k = klower+1,kupper
          do j = jlower+1,jupper
            do i = ilower+1,iupper
              write(file3) xo(i,j,k), yo(i,j,k), zo(i,j,k)
            end do
          end do
        end do
c
        do n = 1, nl
          do k = klower,kupper
            do j = jlower, jupper
              do i = ilower,iupper
                write(file3) qii(i,j,k,n)
              end do
            end do
          end do
        end do
        close(file3)
      end if
c     
      if (rank.eq.0) write(*,*) 'wrote ', file_rstore
      if (dual_t.eq.1.and.rank.eq.0) write(*,*) 'wrote ', file_rstorep
c
      return
      end
