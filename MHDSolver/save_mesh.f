      subroutine save_mesh(il, jl, kl, ilower, iupper, jlower, jupper,
     $     klower, kupper, x, y, z, xo, yo, zo, udi, vdi, wdi, udj, vdj,
     $     wdj, udk, vdk, wdk, xsl, xzl, time, nmode)
c
c     save mesh info of the previous time step
c
      implicit none
c
c     interface variables
c
      integer, intent(in)::
     $     il, jl, kl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     nmode
c
      double precision, intent(in)::
     $     time
c
      double precision, intent(inout),
     $     dimension(ilower+1:iupper,jlower+1:jupper,klower+1:kupper)::
     $     x, y, z, xo, yo, zo
c
      double precision, intent(inout)::
     $     xsl(4,3),
     $     xzl(2,3,6)
c
      double precision, dimension(jlower:jupper,
     $     klower:kupper,4), intent(inout):: udi, vdi, wdi
c
      double precision, dimension(ilower:iupper,
     $     klower:kupper,4), intent(inout):: udj, vdj, wdj
c
      double precision, dimension(ilower:iupper,
     $     jlower:jupper,4), intent(inout):: udk, vdk, wdk
c
c     local variables
c
      integer::
     $     j
c
c     subroutine begin
c
      xsl(1,3) = xsl(1,2)
      xsl(2,3) = xsl(2,2)
      xsl(3,3) = xsl(3,2)
      xsl(4,3) = xsl(4,2)
c
      xsl(1,2) = xsl(1,1)
      xsl(2,2) = xsl(2,1)
      xsl(3,2) = xsl(3,1)
      xsl(4,2) = xsl(4,1)
c
      do 10 j = 1, nmode
        xzl(1,3,j) = xzl(1,2,j)
        xzl(2,3,j) = xzl(2,2,j)
c
        xzl(1,2,j) = xzl(1,1,j)
        xzl(2,2,j) = xzl(2,1,j)
   10 continue
c
      xo = x
      yo = y
      zo = z
c
      if (time.lt.0.001) then
        udi = 0.d0
        vdi = 0.d0
        wdi = 0.d0
c
        udj = 0.d0
        vdj = 0.d0
        wdj = 0.d0
c
        udk = 0.d0
        vdk = 0.d0
        wdk = 0.d0
      else
        udi(:,:,3:4) = udi(:,:,1:2)
        vdi(:,:,3:4) = vdi(:,:,1:2)
        wdi(:,:,3:4) = wdi(:,:,1:2)
c
        udj(:,:,3:4) = udj(:,:,1:2)
        vdj(:,:,3:4) = vdj(:,:,1:2)
        wdj(:,:,3:4) = wdj(:,:,1:2)
c
        udk(:,:,3:4) = udk(:,:,1:2)
        vdk(:,:,3:4) = vdk(:,:,1:2)
        wdk(:,:,3:4) = wdk(:,:,1:2)
      end if
c
      end subroutine save_mesh
c
c======================================================================c
