      subroutine record_p(dstep, rank, time, il, jl, kl, nl, ilower,
     $     iupper, jlower, jupper, klower, kupper, x, y, q, gamma,
     $     dual_t)
c
c     record pressure history for oscillating cascade
c
      implicit none
c
c     Interface Variables
c
      integer, intent(in)::
     $     dstep, rank,
     $     il, jl, kl, nl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper,
     $     dual_t

      double precision, intent(in)::
     $     x(ilower+1:iupper,jlower+1:jupper,klower+1:kupper),
     $     y(ilower+1:iupper,jlower+1:jupper,klower+1:kupper),
     $     q(ilower:iupper,jlower:jupper,klower:kupper,nl),
     $     gamma,
     $     time
c
c     Local Variables
c      
      double precision, dimension(jl)::
     $     r, u, v, vel, pin, area, m

      double precision, dimension(il)::
     $     psurl, psuru
      
      character(len=20)::
     $     filename
c
c     Begin
c
c..   compute inlet r, vel, p, area and m
      
      r = q(1,1:jl,1,1)
      u = q(1,1:jl,1,2)/r
      v = q(1,1:jl,1,3)/r
      vel = dsqrt(u*u+v*v)
      pin = (gamma-1.d0)*(q(1,1:jl,1,5)-0.5d0*r*(u*u+v*v))
      area = y(1,2:jl+1,1)-y(1,1:jl,1)
      m = r*u*area
      
c..   compute surface pressure

      psurl = (gamma-1.d0)*
     $     (q(1:il,1,1,5)-(q(1:il,1,1,2)**2+q(1:il,1,1,3**2))
     $     /(2.*q(1:il,1,1,1)))
      psuru = (gamma-1.d0)*
     $     (q(1:il,jl,1,5)-(q(1:il,jl,1,2)**2+q(1:il,jl,1,3**2))
     $     /(2.*q(1:il,jl,1,1)))

c..   record unsteady parameter

      write(filename, '("pre",i2.2,".dat")') rank

      open(1,file=filename, position='append', form='unformatted')
      
      if(dual_t.eq.1) then
         write(1) dstep,r, vel, pin, area, m, psurl, psuru,
     $        x(1:il+1,1,1), y(1:il+1,1,1), x(1:il+1,jl+1,1),
     $        y(1:il+1,jl+1,1)
      else
         write(1) dstep,r, vel, pin, area, m, psurl, psuru,
     $        x(1:il+1,1,1), y(1:il+1,1,1), x(1:il+1,jl+1,1),
     $        y(1:il+1,jl+1,1), time
      end if

      close(1)

      end
