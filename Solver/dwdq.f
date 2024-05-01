        subroutine dwdq(nl,cp,rhop,rhot,th,rho,u,v,w,k)

c         calculate (dw/dq)  w is conservational variables, 
c                       q is original variables.

c    rhop d(rho)/dp
c    th is total enthalpy
c    k  is  (dw/dq)
c    rhot is d(rho)/dT 

      implicit none
      integer,intent(in):: 
     $    nl

      double precision,intent(in)::
     $    cp,rhop,rhot,th,rho,u,v,w

      double precision,intent(out)::  
     $    k(nl,nl)


	  k(1,1)=rhop
	  k(1,2)=0.0d0
	  k(1,3)=0.0d0
	  k(1,4)=0.0d0
	  k(1,5)=rhot
	  
	  k(2,1)=rhop*u
	  k(2,2)=rho
	  k(2,3)=0.0d0
	  k(2,4)=0.0d0
	  k(2,5)=rhot*u

	  k(3,1)=rhop*v
	  k(3,2)=0.0d0
	  k(3,3)=rho
	  k(3,4)=0.0d0
	  k(3,5)=rhot*v

	  k(4,1)=rhop*w
	  k(4,2)=0.0d0
	  k(4,3)=0.0d0
	  k(4,4)=rho
	  k(4,5)=rhot*w

	  k(5,1)=rhop*th-1.0d0
	  k(5,2)=rho*u
	  k(5,3)=rho*v
	  k(5,4)=rho*w
	  k(5,5)=rhot*th+cp*rho

	  return
	end

