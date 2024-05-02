      subroutine source_move(dim2, il, jl, kl, x, y, z, vol,
     $     ilower, iupper, jlower, jupper, klower, kupper,
     $     idimen, dt, tintvl, qt, sr)
c
c     This subroutine is used to calculate the source term due to
c     the geometric conservation law.
c
      implicit none
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     dim2,                ! maxium face number in 3 directions
     $     il, jl, kl,          ! cell number in 3 directions
     $     idimen,              ! number of dimensions
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper

      double precision, intent(in):: 
     $     dt(il,jl,kl),
     $     tintvl,
     $     qt(ilower:iupper,jlower:jupper,klower:kupper,3)

      double precision, dimension(il,jl,kl), intent(out):: sr

      double precision, intent(in),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z

      double precision, intent(in), dimension(ilower+1:iupper-1,
     $     jlower+1:jupper-1, klower+1:kupper-1)::
     $     vol
c
c     LOCAL VARIABLES
c
      double precision::
     $     tinv, xdt(0:dim2), ydt(0:dim2), zdt(0:dim2), st(0:dim2),
     $     dx1, dy1, dz1, dx2, dy2, dz2, ax, ay, az,
     $     cv2x                 ! coef to convert negative vel to distance

      integer::
     $     i, j, k,             ! cell iteration index
     $     ip, jp, kp,          ! i+1, j+1, k+1
     $     ilp, jlp, klp        ! il+1, jl+1, kl+1
c
c *** SUBROUTINE START ***
c
c *** set up some parameters
c
      ilp = il + 1
      jlp = jl + 1
      klp = kl + 1
c
      tinv = 0.03125d0/tintvl
      cv2x = -8.d0*tintvl
c
c *** xi direction
c     
      do k = 1,kl
        do j = 1,jl
c     
c     ... grid moving velocity
c     
          do i = 0,il+1
            xdt(i) = qt(i,j,k,1)*cv2x
            ydt(i) = qt(i,j,k,2)*cv2x
            zdt(i) = qt(i,j,k,3)*cv2x
          end do
c     
c     ... area vector at xi surfece
c     
          jp = j + 1
          kp = k + 1
          do i = 1,ilp
            dx1 = x(i,jp,kp) - x(i,j,k)
            dy1 = y(i,jp,kp) - y(i,j,k)
            dz1 = z(i,jp,kp) - z(i,j,k)
            dx2 = x(i,j,kp)  - x(i,jp,k)
            dy2 = y(i,j,kp)  - y(i,jp,k)
            dz2 = z(i,j,kp)  - z(i,jp,k)
c     
            ax = dy1 * dz2 - dz1 * dy2
            ay = dz1 * dx2 - dx1 * dz2 
            az = dx1 * dy2 - dy1 * dx2 
c     
            st(i) = ((xdt(i)+xdt(i-1))*ax
     $              +(ydt(i)+ydt(i-1))*ay+(zdt(i)+zdt(i-1))*az)*tinv
          end do
c     
c     ... integrated on the xi surface for each volume
c     
          do i = 1,il
            sr(i,j,k) = st(i) - st(i+1)
          end do
        end do
      end do
c
c *** eta direction
c     
      do k = 1,kl
        do i = 1,il
c     
c     ... grid moving velocity
c     
          do j = 0,jl+1
            xdt(j) = qt(i,j,k,1)*cv2x
            ydt(j) = qt(i,j,k,2)*cv2x
            zdt(j) = qt(i,j,k,3)*cv2x
          end do
c     
c     ... area vector at eta surfece
c     
          ip = i + 1
          kp = k + 1
          do j = 1,jlp
            dx1 = x(ip,j,kp) - x(i,j,k)
            dy1 = y(ip,j,kp) - y(i,j,k)
            dz1 = z(ip,j,kp) - z(i,j,k)
            dx2 = x(ip,j,k)  - x(i,j,kp)
            dy2 = y(ip,j,k)  - y(i,j,kp)
            dz2 = z(ip,j,k)  - z(i,j,kp)
c     
            ax = dy1 * dz2 - dz1 * dy2
            ay = dz1 * dx2 - dx1 * dz2
            az = dx1 * dy2 - dy1 * dx2
c     
            st(j) = ((xdt(j)+xdt(j-1))*ax
     $              +(ydt(j)+ydt(j-1))*ay+(zdt(j)+zdt(j-1))*az)*tinv
          end do
c     
c     ... integrated on the eta surface for each volume
c     
          do j = 1,jl
            sr(i,j,k) = sr(i,j,k) + st(j)-st(j+1)
          end do
        end do
      end do
c
c *** zeta direction
c     
      if (idimen.eq.3) then
        do i = 1, il
          do j = 1, jl
c     
c     ... grid moving velocity
c     
            do k = 0,kl+1
              xdt(k) = qt(i,j,k,1)*cv2x
              ydt(k) = qt(i,j,k,2)*cv2x
              zdt(k) = qt(i,j,k,3)*cv2x
            end do
c     
c     ... area vector at zeta surfece
c     
            ip = i + 1
            jp = j + 1
            do k = k, klp
              dx1 = x(ip,jp,k) - x(i,j,k)
              dy1 = y(ip,jp,k) - y(i,j,k)
              dz1 = z(ip,jp,k) - z(i,j,k)
              dx2 = x(i,jp,k)  - x(ip,j,k)
              dy2 = y(i,jp,k)  - y(ip,j,k)
              dz2 = z(i,jp,k)  - z(ip,j,k)
c     
              ax = dy1 * dz2 - dz1 * dy2
              ay = dz1 * dx2 - dx1 * dz2
              az = dx1 * dy2 - dy1 * dx2
c     
              st(k) = ((xdt(k)+xdt(k-1))*ax
     $                 +(ydt(k)+ydt(k-1))*ay+(zdt(k)+zdt(k-1))*az)*tinv
            end do
c     
c     ... integrated on the eta surface for each volume
c     
            do k = 1, kl
              sr(i,j,k) = sr(i,j,k) + st(k)-st(k+1)
            end do
          end do
        end do
      end if
c
      do k = 1,kl
        do j = 1,jl
          do i = 1,il
            sr(i,j,k) = sr(i,j,k)/vol(i,j,k)
          end do
        end do
      end do
c
      return
      end
