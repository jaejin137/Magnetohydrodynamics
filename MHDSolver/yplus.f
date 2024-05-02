      function distoplane(x0, y0, z0, x, y, z)
c
c     compute distance from point (x0,y0,z0) to surface (x,y,z)
c
      implicit none
c
      double precision, intent(in)::
     $     x0, y0, z0           ! difining point
      
      double precision, dimension(1), intent(in)::
     $     x, y, z              ! defining surface

      double precision::
     $     distoplane

      double precision, dimension(3)::
     $     line, x1, x2, n
      
      double precision, dimension(3,3)::
     $     p
c
      p(:,1) = x(:3)
      p(:,2) = y(:3)
      p(:,3) = z(:3)

      x1 = p(2,:)-p(1,:)
      x2 = p(3,:)-p(1,:)
      
      call cross_product(x1,x2,n)
      n = n/dsqrt(dot_product(n,n))

      line(1) = x0
      line(2) = y0
      line(3) = z0
      line = line-p(1,:)

      distoplane = abs(dot_product(n,line))
      
      end

c-------------------------

      subroutine cross_product(x1,x2, cross)
c
      implicit none
c
      double precision, dimension(3), intent(in)::
     $     x1, x2

      double precision, dimension(3)::
     $     cross
c
      cross(1) = x1(2)*x2(3)-x1(3)*x2(2)
      cross(2) = x1(3)*x2(1)-x1(1)*x2(3)
      cross(3) = x1(1)*x2(2)-x1(2)*x2(1)
      
      end
