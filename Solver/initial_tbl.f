c======================================================================c
c
      subroutine initial_tbl(rank, il, jl, kl, nl, ilower, iupper,
     $     jlower, jupper, klower, kupper, nstatis,time,
     $       rave,uave,vave,wave,eave)
c
c     read rstart to initialize flowfield
c
c     IMPLICIT STATEMENT
c
      implicit none
c
c     INTERFACE VARIABLES       ! specified in aamain.f 
c
      integer, intent(in)::
     $     rank,
     $     il, jl, kl,
     $     nl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper

      integer, intent(inout)::
     $     nstatis

c
      double precision, intent(out)::
     $     time
c
c
      double precision, intent(out),
     $     dimension(ilower:iupper,jlower:jupper,klower:kupper)::
     $     rave, uave, vave, wave,eave
c
c     Local Variables
c
      character(len=20)::
     $     filename                
c
      logical::
     $     ex
c
      integer::
     $     i, j, k
c     
      integer, parameter::
     $     file1 = 1
c
c --- Begin ---
c      
c ... file name and existance check
c
      write(filename, '("statistical",i3.3)') rank
      inquire(file=filename, exist=ex)
      if (ex) then
        open(file1, file=filename, status='unknown', 
     >       access='sequential', form= 'unformatted')
c
c ... read rstart
c
      if (rank.eq.1) write(6,*) 'Start reading statistical average ...'
c
      read(file1)  nstatis, time
c
      do k = klower+1,kupper
        do j = jlower+1,jupper
          do i = ilower+1,iupper
            read(file1) rave(i,j,k), uave(i,j,k), 
     $                  vave(i,j,k),wave(i,j,k),eave(i,j,k)
          end do
        end do
      end do
      close(file1)

      else
        write(*,*) 'not statistical values'
        nstatis = 0
        rave = 0.d0
        uave = 0.d0
        vave = 0.d0
        wave = 0.d0
        eave = 0.d0
      end if
      return
      end
c
c======================================================================c
