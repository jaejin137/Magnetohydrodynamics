      subroutine exchange_res(rank, np, it0, resdulave, resdulmax, loc)
c
c     exchange time step in mpi
c
c     IMPLICIT STATEMENT
c
      implicit none
c
c     include header file
c
      include "/usr/local/include/mpif.h"
c
c     INTERFACE VARIABLES
c
      double precision, intent(inout):: resdulave, resdulmax
      integer, intent(inout):: loc(4)
      integer, intent(in):: rank, np, it0
c
c     LOCAL VARIABLES
c     
      double precision:: res0, res1
      integer:: status(mpi_status_size), ierr
      integer:: tag0 = 1, tag1 = 2, tag2 = 3
      integer:: i, loc0(4)
c
c *** START SUBROUTINE ***      
c
      if(rank.ne.0) then
         call mpi_send(resdulave, 1, mpi_double_precision, 0, tag0,
     $        mpi_comm_world, ierr)
         call mpi_send(resdulmax, 1, mpi_double_precision, 0, tag1,
     $        mpi_comm_world, ierr)
         call mpi_send(loc, 4, mpi_integer, 0, tag2,
     $        mpi_comm_world, ierr)
         call mpi_recv(resdulave, 1, mpi_double_precision, 0,
     $        tag0, mpi_comm_world, status, ierr)
         call mpi_recv(resdulmax, 1, mpi_double_precision, 0,
     $        tag1, mpi_comm_world, status, ierr)
         call mpi_recv(loc, 4, mpi_integer, 0,
     $        tag2, mpi_comm_world, status, ierr)
      else
         do i = 1, np-1
           call mpi_recv(res0, 1, mpi_double_precision, i,
     $          tag0, mpi_comm_world, status, ierr)
           resdulave = resdulave + res0
           call mpi_recv(res1, 1, mpi_double_precision, i,
     $          tag1, mpi_comm_world, status, ierr)
           call mpi_recv(loc0, 4, mpi_integer, i,
     $          tag2, mpi_comm_world, status, ierr)
           if (resdulmax.lt.res1) then
             resdulmax = res1
             loc = loc0
           end if
         end do
         resdulave = resdulave/dfloat(it0)
         do i = 1, np-1
           call mpi_send(resdulave, 1, mpi_double_precision, i, tag0,
     $           mpi_comm_world, ierr)
           call mpi_send(resdulmax, 1, mpi_double_precision, i, tag1,
     $           mpi_comm_world, ierr)
           call mpi_send(loc, 4, mpi_integer, i, tag2,
     $           mpi_comm_world, ierr)
         end do

      end if
      return
      end
