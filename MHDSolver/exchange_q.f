      subroutine exchange_q(fn,nb,blen,nl,
     $     bc_max,bnum,btb,bnb,bdir,bstart,bend,order,idimen,
     $     ilow,iupp,jlow,jupp,klow,kupp,x)
c
c     update q0 by exchange the values of inner BC
c
      implicit none
c
      include "/usr/local/include/mpif.h"
      integer, intent(in)::
     $     nb, fn, idimen,
     $     nl, bc_max 
c
c *** the variables of block connection BC
c
      integer, intent(in):: ilow,iupp,jlow,jupp,
     $                      klow,kupp
      integer, intent(in):: bnum(nb),btb(nb,bc_max),bnb(nb,bc_max),
     $   bdir(nb,bc_max),bstart(3,nb,bc_max),bend(3,nb,bc_max),
     $   order(nb,bc_max)
c
c *** the exchange variables
c
      double precision,intent(inout)::
     $  x(ilow:iupp,jlow:jupp,klow:kupp,nl)
c     
c     Local
c     
c *** intermidiate array
c
      double precision, allocatable, dimension(:):: bcb
c
      integer::
     $     start(3),end(3),
     $     blen,ic(3),sb(3),eb(3)
c
      integer:: i, j, k, n, l, l1, ib, inum, ii, i1, ijk_max
      integer:: status(mpi_status_size), ierr
c     
      if (bnum(fn).gt.0) then
c
c *** begin exchange
c     
        l=blen
        l1=blen-1
        if(idimen.eq.2) then
c     
c *** 2-D problem
c
          do ii=1,bnum(fn)
            ib=btb(fn,ii)
            inum=bnb(fn,ii)
            do i=1,2
              start(i)=bstart(i,fn,ii)
              end(i)=bend(i,fn,ii)
              ic(i)=abs(end(i)-start(i))+1+2*blen
            end do 
c
c ***allocate the memory for exchange array bcb
c
            ijk_max=blen*nl*ic(1)*ic(2)
            allocate(bcb(ijk_max))
            do i=1,2
              if (end(i).lt.start(i)) then
                ic(i)=-1
                eb(i)=-blen
                sb(i)=blen
              else if (end(i).gt.start(i)) then
                ic(i)=1
                eb(i)=blen
                sb(i)=-blen
              end if
            end do
            i1=0
            if (fn.lt.ib) then
c
c ***send the variables to the connected block ib
c
c ***xi direction
              if(bdir(fn,ii).eq.1) then
                do i=start(1),start(1)+l1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.2) then
                do i=start(1),start(1)-l1,-1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
c ***eta direction
              else if (bdir(fn,ii).eq.3) then
                do j=start(2),start(2)+l1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.4) then
                do j=start(2),start(2)-l1,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
              else
                write(*,*)'inner BC define error in block',ib
                write(*,*)'bdir should be 4 or less in 2-D problem'
                stop
              end if
              call mpi_send(bcb,ijk_max,mpi_double_precision,ib-1,inum,
     $             mpi_comm_world,ierr)
c
c ***send complete 
c
            else if (fn.gt.ib) then
c
c ***receive messages from connected blocks
c
              call mpi_recv(bcb,ijk_max,mpi_double_precision,ib-1,ii,
     $             mpi_comm_world,status,ierr)
              if(bdir(fn,ii).eq.1) then
                do i=start(1)-1,start(1)-l,-1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else if(bdir(fn,ii).eq.2) then
                do i=start(1)+1,start(1)+l
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else if(bdir(fn,ii).eq.3) then
                do j=start(2)-1,start(2)-l,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else if(bdir(fn,ii).eq.4) then
                do j=start(2)+1,start(2)+l
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else
                write(*,*)'inner BC define error in block',fn
                write(*,*)'bdir should be 4 or less in 2-D problem'
                stop
              end if
            end if
            deallocate(bcb)
          end do

          do ii=1,bnum(fn)
            ib=btb(fn,ii)
            inum=bnb(fn,ii)
            do i=1,2
              start(i)=bstart(i,fn,ii)
              end(i)=bend(i,fn,ii)
              ic(i)=abs(end(i)-start(i))+1+2*blen
            end do 
c
c ***allocate the memory for exchange array bcb
c
            ijk_max=blen*nl*ic(1)*ic(2)
            allocate(bcb(ijk_max))
            do i=1,2
              if (end(i).lt.start(i)) then
                ic(i)=-1
                eb(i)=-blen
                sb(i)=blen
              else if (end(i).gt.start(i)) then
                ic(i)=1
                eb(i)=blen
                sb(i)=-blen
              end if
            end do
            i1=0
c
            if (fn.gt.ib) then
c
c ***send the variables to the connected block ib
c
c ***xi direction
              if(bdir(fn,ii).eq.1) then
                do i=start(1),start(1)+l1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.2) then
                do i=start(1),start(1)-l1,-1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
c ***eta direction
              else if (bdir(fn,ii).eq.3) then
                do j=start(2),start(2)+l1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.4) then
                do j=start(2),start(2)-l1,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      bcb(i1)=x(i,j,1,n)
                    end do
                  end do
                end do
              else
                write(*,*)'inner BC define error in block',ib
                write(*,*)'bdir should be 4 or less in 2-D problem'
                stop
              end if
              call mpi_send(bcb,ijk_max,mpi_double_precision,ib-1,inum,
     $             mpi_comm_world,ierr)
c
c ***send complete 
c
            else if (fn.lt.ib) then
c
c ***receive messages from connected blocks
c
              call mpi_recv(bcb,ijk_max,mpi_double_precision,ib-1,ii,
     $             mpi_comm_world,status,ierr)
              if(bdir(fn,ii).eq.1) then
                do i=start(1)-1,start(1)-l,-1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else if(bdir(fn,ii).eq.2) then
                do i=start(1)+1,start(1)+l
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else if(bdir(fn,ii).eq.3) then
                do j=start(2)-1,start(2)-l,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else if(bdir(fn,ii).eq.4) then
                do j=start(2)+1,start(2)+l
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do n=1,nl
                      i1=i1+1
                      x(i,j,1,n)=bcb(i1)
                    end do
                  end do
                end do
              else
                write(*,*)'inner BC define error in block',fn
                write(*,*)'bdir should be 4 or less in 2-D problem'
                stop
              end if
            end if
            deallocate(bcb)
          end do
c
c ***exchange complete 
c
          do j=jlow,jupp
            do i=ilow,iupp
              do n=1,nl
                x(i,j,2,n)=x(i,j,1,n)
              end do
            end do
          end do
c
        else if(idimen.eq.3) then
c
c ***exchange information of the connected blocks
c
          do ii=1,bnum(fn)
            ib=btb(fn,ii)
            inum=bnb(fn,ii)
            do i=1,3
              start(i)=bstart(i,fn,ii)
              end(i)=bend(i,fn,ii)
              ic(i)=abs(end(i)-start(i))+1+2*blen
            end do 
c
c ***allocate the memory for exchange array bcb
c
            ijk_max=blen*nl*ic(1)*ic(2)*ic(3)
            allocate(bcb(ijk_max))
            do i=1,3
              if (end(i).lt.start(i)) then
                ic(i)=-1
                eb(i)=-blen
                sb(i)=blen
              else if (end(i).gt.start(i)) then
                ic(i)=1
                eb(i)=blen
                sb(i)=-blen
              end if
            end do
            i1=0
c
            if (fn.lt.ib) then
c
c ***send connective BC values to ib 
c
              if (bdir(fn,ii).eq.1) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.3.or.
     $              order(fn,ii).eq.5) then
                  do i=start(1),start(1)+l1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.2.or.
     $                   order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.6) then
                  do i=start(1),start(1)+l1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.2) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.3.or.
     $              order(fn,ii).eq.5) then
                  do i=start(1),start(1)-l1,-1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.2.or.
     $                   order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.6) then
                  do i=start(1),start(1)-l1,-1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.3) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.3) then
                  do j=start(2),start(2)+l1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do j=start(2),start(2)+l1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.4) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.3) then
                  do j=start(2),start(2)-l1,-1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do j=start(2),start(2)-l1,-1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.5) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.4) then
                  do k=start(3),start(3)+l1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.3.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do k=start(3),start(3)+l1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.6) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.4) then
                  do k=start(3),start(3)-l1,-1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.3.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do k=start(3),start(3)-l1,-1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              end if
              call mpi_send(bcb,ijk_max,mpi_double_precision,ib-1,inum,
     $             mpi_comm_world,ierr)

            else if (fn.gt.ib) then

              call mpi_recv(bcb,ijk_max,mpi_double_precision,ib-1,ii,
     $             mpi_comm_world,status,ierr)
              if (bdir(fn,ii).eq.1) then
                do i=start(1)-1,start(1)-l,-1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.2) then
                do i=start(1)+1,start(1)+l
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.3) then
                do j=start(2)-1,start(2)-l,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.4) then
                do j=start(2)+1,start(2)+l
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.5) then
                do k=start(3)-1,start(3)-l,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.6) then
                do k=start(3)+1,start(3)+l
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else
                write(*,*)'error in bdir,block,number=',fn,ii
              end if
            end if
            deallocate(bcb)
          end do

          do ii=1,bnum(fn)
            ib=btb(fn,ii)
            inum=bnb(fn,ii)
            do i=1,3
              start(i)=bstart(i,fn,ii)
              end(i)=bend(i,fn,ii)
              ic(i)=abs(end(i)-start(i))+1+2*blen
            end do 
c
c ***allocate the memory for exchange array bcb
c
            ijk_max=blen*nl*ic(1)*ic(2)*ic(3)
            allocate(bcb(ijk_max))
            do i=1,3
              if (end(i).lt.start(i)) then
                ic(i)=-1
                eb(i)=-blen
                sb(i)=blen
              else if (end(i).gt.start(i)) then
                ic(i)=1
                eb(i)=blen
                sb(i)=-blen
              end if
            end do
            i1=0
c
            if (fn.gt.ib) then
c
c ***send connective BC values to ib 
c
              if (bdir(fn,ii).eq.1) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.3.or.
     $              order(fn,ii).eq.5) then
                  do i=start(1),start(1)+l1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.2.or.
     $                   order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.6) then
                  do i=start(1),start(1)+l1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.2) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.3.or.
     $              order(fn,ii).eq.5) then
                  do i=start(1),start(1)-l1,-1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.2.or.
     $                   order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.6) then
                  do i=start(1),start(1)-l1,-1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.3) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.3) then
                  do j=start(2),start(2)+l1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do j=start(2),start(2)+l1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.4) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.3) then
                  do j=start(2),start(2)-l1,-1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.4.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do j=start(2),start(2)-l1,-1
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.5) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.4) then
                  do k=start(3),start(3)+l1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.3.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do k=start(3),start(3)+l1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              else if (bdir(fn,ii).eq.6) then
                if (order(fn,ii).eq.1.or.
     $              order(fn,ii).eq.2.or.
     $              order(fn,ii).eq.4) then
                  do k=start(3),start(3)-l1,-1
                    do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                      do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else if (order(fn,ii).eq.3.or.
     $                   order(fn,ii).eq.5.or.
     $                   order(fn,ii).eq.6) then
                  do k=start(3),start(3)-l1,-1
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                        do n=1,nl
                          i1=i1+1
                          bcb(i1)=x(i,j,k,n)
                        end do
                      end do
                    end do
                  end do
                else
                  write(*,*)'error in iorder(ib,inum).ib,inum=',ib,inum
                end if
              end if
              call mpi_send(bcb,ijk_max,mpi_double_precision,ib-1,inum,
     $             mpi_comm_world,ierr)

            else if (fn.lt.ib) then

              call mpi_recv(bcb,ijk_max,mpi_double_precision,ib-1,ii,
     $             mpi_comm_world,status,ierr)
              if (bdir(fn,ii).eq.1) then
                do i=start(1)-1,start(1)-l,-1
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.2) then
                do i=start(1)+1,start(1)+l
                  do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.3) then
                do j=start(2)-1,start(2)-l,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.4) then
                do j=start(2)+1,start(2)+l
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do k=start(3)+sb(3),end(3)+eb(3),ic(3)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.5) then
                do k=start(3)-1,start(3)-l,-1
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else if (bdir(fn,ii).eq.6) then
                do k=start(3)+1,start(3)+l
                  do i=start(1)+sb(1),end(1)+eb(1),ic(1)
                    do j=start(2)+sb(2),end(2)+eb(2),ic(2)
                      do n=1,nl
                        i1=i1+1
                        x(i,j,k,n)=bcb(i1)
                      end do
                    end do
                  end do
                end do
              else
                write(*,*)'error in bdir,block,number=',fn,ii
              end if
            end if
            deallocate(bcb)
          end do
        end if
      end if

      return
      end
