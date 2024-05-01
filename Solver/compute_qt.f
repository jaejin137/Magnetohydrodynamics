      subroutine compute_qt(il, jl, kl, ilower, iupper, jlower, jupper,
     $     klower, kupper, x, y, z, xo, yo, zo, tintvl, dim2,
     $     qt, bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper,
     $     bc_zeta_lower, bc_zeta_upper)
c
c     compute moving grid velocity at each computation cell center
c
      implicit none
c
      integer, intent(in)::
     $     il, jl, kl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     dim2

      double precision, intent(in)::
     $     tintvl

      double precision, intent(inout),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z,
     $     xo, yo, zo

      integer, intent(in)::
     $     bc_xi_lower(jl, kl), bc_xi_upper(jl, kl),
     $     bc_eta_lower(il, kl), bc_eta_upper(il, kl),
     $     bc_zeta_lower(il, jl), bc_zeta_upper(il, jl)

      double precision, intent(out)::
     $     qt(ilower:iupper,jlower:jupper,klower:kupper,3)
c
      integer::
     $     index,
     $     i, j, k

      double precision::
     $     qtline(ilower:dim2,3),
     $     bc
c
      integer, parameter::
     $     peri = 10,           ! periodical boundary
     $     innr = 7             ! inner boundary for mpi
c
c     subroutine start
c
      index = 1
      do j = 1, jl
         do k = 1, kl
            call spvel(index, j, k, dim2, il, jl, kl, x, y, z, xo,
     $           yo, zo, ilower, iupper, jlower, jupper, klower, kupper,
     $           tintvl, qtline)
            qt(1:il,j,k,:) = qtline(1:il,:)
         end do
      end do
      
c..   boundaries

c..   xi_lower

      do j = 1, jl
        do k = 1, kl
          bc = bc_xi_lower(j,k)
          if(bc.eq.peri.or.bc.eq.20) then
            qt(ilower:0,j,k,:) = qt(il+ilower:il,j,k,:)
          else if(bc.eq.innr) then
     $              ! assigned later by exchange_qt
          else
            do i=ilower,0
              qt(i,j,k,:) = qt(1,j,k,:)
            end do
          end if
        end do
      end do

c..   xi_upper

      do j = 1, jl
        do k = 1, kl
          bc = bc_xi_upper(j,k)
          if(bc.eq.peri.or.bc.eq.20) then
            qt(il+1:iupper,j,k,:) = qt(1:iupper-il,j,k,:)
          else if(bc.eq.innr) then
          else
            do i=il+1,iupper
              qt(i,j,k,:) = qt(il,j,k,:)
            end do
          end if
        end do
      end do

c..   eta_lower

      do i = 1, il
        do k = 1, kl
          bc = bc_eta_lower(i,k)
          if(bc.eq.peri.or.bc.eq.20) then
            qt(i,jlower:0,k,:) = qt(i,jl+jlower:jl,k,:)
          else if(bc.eq.innr) then
          else
            do j=jlower,0
              qt(i,j,k,:) = qt(i,1,k,:)
            end do
          end if
        end do
      end do

c..   eta_upper

      do i = 1, il
        do k = 1, kl
          bc = bc_eta_upper(i,k)
          if(bc.eq.peri.or.bc.eq.20) then
            qt(i,jl+1:jupper,k,:) = qt(i,1:jupper-jl:jl,k,:)
          else if(bc.eq.innr) then
          else
            do j=jl+1,jupper
              qt(i,j,k,:) = qt(i,jl,k,:)
            end do
          end if
        end do
      end do

c..   zeta_lower

      do i = 1, il
        do j = 1, jl
          bc = bc_zeta_lower(i,j)
          if(bc.eq.peri.or.bc.eq.20) then
            qt(i,j,klower:0,:) = qt(i,j,kl+klower:kl,:)
          else if(bc.eq.innr) then
          else
            do k=klower,0
              qt(i,j,k,:) = qt(i,j,1,:)
            end do
          end if
        end do
      end do

c..   zeta_upper

      do i = 1, il
        do j = 1, jl
          bc = bc_zeta_upper(i,j)
          if(bc.eq.peri.or.bc.eq.20) then
            qt(i,j,kl+1:kupper,:) = qt(i,j,1:kupper-kl,:)
          else if(bc.eq.innr) then
          else
            do k=kl+1,kupper
              qt(i,j,k,:) = qt(i,j,kl,:)
            end do
          end if
        end do
      end do

      end
