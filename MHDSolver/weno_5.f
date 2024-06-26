      subroutine weno_5(iside,imax,dim1,dim2,
     $     bc_lower,bc_upper,u,du,ul,ur, control)

c------------iside=0   1-order scheme-------------------------
c------------iside=2   2-order scheme-------------------------
c------------iside=3   3-order scheme-------------------------
c------------iside=4   weno3-order scheme-------------------------
c------------iside=5   5-order scheme-------------------------
c------------iside=6   weno5-order scheme-------------------------
c------------iside=106 Borges weno5-order scheme-------------------------

c
c     reconstruction to the cell faces 
c     using conservative variables
c
      use datain

      implicit none

      type (datain_type), intent(in)::control
      
c
c     INTERFACE VARIABLES
c
      integer, intent(in)::
     $     iside, imax,
     $     dim1, dim2,bc_lower,bc_upper

      double precision, dimension(dim1:dim2), intent(in)::
     $     u, du

      double precision, dimension(dim2), intent(out)::
     $     ul, ur
 
      integer 
     $     i,first,ii,j, wall_order
      logical::
     $     bndl, bndu


c-------------for 5-weno----------------
      double precision::
     $     eps,c20,c21,c22,is0,is1,is2,
     $     a0,a1,a2,w0,w1,w2,w3,w4,q0,q1,q2,
     $     b0

c------------for 7-weno------------------
      double precision::
     $   c40,c41,c42,c43,sss,
     $   q40,q41,q42,q43,
     $   c0,c1,c2,c3,c4,
     $   r0,r1,r2,r3,r4

c-----------for 9-weno-------------------
      double precision::
     $   c50,c51,c52,c53,c54,
     $   q50,q51,q52,q53,q54

	double precision,dimension(:)::
     $   A9(5,4),B9(5,4),C9(5,4),D9(5,4),E9(5,4),
     $   CC(4),CCC(5)
c------------for Borges WENO5
      double precision::
     $   tao5
c-----------------------------------------    
        wall_order=control%wall_order

	eps=control%varepsilon
c-----------------------------------------
        sss=eps
	c20=0.1
	c21=0.6
	c22=0.3
	b0=13./12.

c------------------note: du(i)=u(i)-u(i-1)
        select case(iside)

c-----------
        case (0)
c----------------------first-order------------------        
        if (iside.eq.0) then
           do i=1,imax
              ul(i)=u(i-1)
              ur(i)=u(i)
           end do
           return
        end if

c-----------
        case (2)
c---------------------------TVD----------------------
      first = 0
      do i = first+1,imax-first
         if (abs(du(i-1)).le.abs(du(i))) then 
          ul(i)=u(i-1)+0.5*du(i-1)
         else
          ul(i)=u(i-1)+0.5*du(i)
         end if
         
         if (abs(du(i)).le.abs(du(i+1))) then
          ur(i)=u(i)-0.5*du(i)
         else
          ur(i)=u(i)-0.5*du(i+1)
         end if
       end do

c-----------
        case (3)       
c-------------------------third-order---------------------
      first = 0
      do i = first+1,imax-first
           q0=u(i-1)+0.5*du(i-1)
           q1=u(i-1)+0.5*du(i)
           ul(i)=0.33333333d0*q0+0.66666667d0*q1
c           ul(i)=1.0d0/3.0d0*q0+2.0d0/3.0d0*q1

           q0=u(i)-0.5*du(i)
           q1=u(i)-0.5*du(i+1)
           ur(i)=0.66666667d0*q0+0.33333333d0*q1
c           ur(i)=2.0d0/3.0d0*q0+1.0d0/3.0d0*q1
      end do
c-----------
        case (4)  
c-------------------------WENO-3---------------------
      first = 0
      do i = first+1,imax-first

           is0=du(i-1)**2
           is1=du(i)**2
           a0=0.333333/(eps+is0)**2
           a1=0.666667/(eps+is1)**2
           w0=a0/(a0+a1)
           w1=a1/(a0+a1)
           q0=u(i-1)+0.5*du(i-1)
           q1=u(i-1)+0.5*du(i)
           ul(i)=w0*q0+w1*q1

           is0=du(i)**2
           is1=du(i+1)**2
           a0=0.666667/(eps+is0)**2
           a1=0.333333/(eps+is1)**2
           w0=a0/(a0+a1)
           w1=a1/(a0+a1)
           q0=u(i)-0.5*du(i)
           q1=u(i)-0.5*du(i+1)
           ur(i)=w0*q0+w1*q1
      end do

c-----------
        case (5)  
c-------------------------fifth-order---------------------

       first=0
       do i = first+1,imax-first
	   q0=u(i-1)+0.5*du(i-1)+0.333333*(du(i-1)-du(i-2))
	   q1=u(i-1)+0.5*du(i-1)+0.333333*(du(i)-du(i-1))
	   q2=u(i-1)+0.5*du(i)-0.166667*(du(i+1)-du(i))
           ul(i)=c20*q0+c21*q1+c22*q2

	   q0=u(i)-0.5*du(i)-0.166667*(du(i)-du(i-1))
	   q1=u(i)-0.5*du(i+1)+0.333333*(du(i+1)-du(i))
	   q2=u(i)-0.5*du(i+1)+0.333333*(du(i+2)-du(i+1))
           ur(i)=c22*q0+c21*q1+c20*q2
	end do


c-----------
        case (6)  
c-------------------------weno5---------------------
       first=0
       do i = first+1,imax-first
	   is0=b0*(du(i-1)-du(i-2))**2+0.25*(3.*du(i-1)-du(i-2))**2
           is1=b0*(du(i)-du(i-1))**2+0.25*(du(i)+du(i-1))**2
	   is2=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)-3.*du(i))**2
	   a0=c20/(eps+is0)**2
	   a1=c21/(eps+is1)**2
	   a2=c22/(eps+is2)**2
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i-1)+0.5*du(i-1)+0.333333*(du(i-1)-du(i-2))
	   q1=u(i-1)+0.5*du(i-1)+0.333333*(du(i)-du(i-1))
	   q2=u(i-1)+0.5*du(i)-0.166667*(du(i+1)-du(i))
	   ul(i)=w0*q0+w1*q1+w2*q2

	   is0=b0*(du(i)-du(i-1))**2+0.25*(3.*du(i)-du(i-1))**2
	   is1=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)+du(i))**2
	   is2=b0*(du(i+2)-du(i+1))**2+0.25*(du(i+2)-3.*du(i+1))**2
	   a0=c22/(eps+is0)**2
	   a1=c21/(eps+is1)**2
	   a2=c20/(eps+is2)**2
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i)-0.5*du(i)-0.166667*(du(i)-du(i-1))
	   q1=u(i)-0.5*du(i+1)+0.333333*(du(i+1)-du(i))
	   q2=u(i)-0.5*du(i+1)+0.333333*(du(i+2)-du(i+1))
	   ur(i)=w0*q0+w1*q1+w2*q2
	end do

c-----------
        case (7)  
c-------------------------weno7(shu's)---------------------
       first=4
	C40=1./35.
	C41=12./35.
	C42=18./35.
        C43=4./35.
       do i = first+1-1,imax-first-1
          ur(i+1)=-1./140.*u(i+4)+5./84.*u(i+3)-101./420.*u(i+2)
     $       +319./420.*u(i+1)+107./210.*u(i)-19./210.*u(i-1)
     $       +1./105.*u(i-2)


        ul(i+1)=-1./140.*u(i-3)+5./84.*u(i-2)-101./420.*u(i-1)
     $       +319./420.*u(i)+107./210.*u(i+1)-19./210.*u(i+2)
     $       +1./105.*u(i+3)
	end do

c-----------
        case (8)  
c-------------------------weno7(shu's)---------------------
       first=4
	C40=1./35.
	C41=12./35.
	C42=18./35.
        C43=4./35.
       do i = first+1-1,imax-first-1
               q40=(-1./4.*u(i+4)+13./12.*u(i+3)
     $              -23./12.*u(i+2)+25./12.*u(i+1))
	       q41=(1./12.*u(i+3)-5./12.*u(i+2)
     $              +13./12.*u(i+1)+1./4.*u(i))
	       q42=(-1./12.*u(i+2)+7./12.*u(i+1)
     $               +7./12.*u(i)-1./12.*u(i-1))
	       q43=(1./4.*u(i+1)+13./12.*u(i)
     $              -5./12.*u(i-1)+1./12.*u(i-2))

	       C0=u(i+4)*(547*u(i+4)-3882*u(i+3)
     $                         +4642*u(i+2)-1854*u(i+1))
     $           +u(i+3)*(7043*u(i+3)-17246*u(i+2)
     $                         +7042*u(i+1))
     $           +u(i+2)*(11003*u(i+2)-9402*u(i+1))
     $                         +2107*u(i+1)**2
	       C1=u(i+3)*(267*u(i+3)-1642*u(i+2)
     $                         +1602*u(i+1)-494*u(i))
     $           +u(i+2)*(2843*u(i+2)-5966*u(i+1)
     $                         +1922*u(i))
     $           +u(i+1)*(3443*u(i+1)-2522*u(i))
     $                         +547*u(i)**2
	       C2=u(i+2)*(547*u(i+2)-2522*u(i+1)
     $                         +1922*u(i)-494*u(i-1))
     $           +u(i+1)*(3443*u(i+1)-5966*u(i)
     $                         +1602*u(i-1))
     $           +u(i)*(2843*u(i)-1642*u(i-1))
     $                         +267*u(i-1)**2
	       C3=u(i+1)*(2107*u(i+1)-9402*u(i)
     $                         +7042*u(i-1)-1854*u(i-2))
     $           +u(i)*(11003*u(i)-17246*u(i-1)
     $                         +4642*u(i-2))
     $           +u(i-1)*(7043*u(i-1)-3882*u(i-2))
     $                         +547*u(i-2)**2

	   R0=C40/(SSS+C0)**2
           R1=C41/(SSS+C1)**2
           R2=C42/(SSS+C2)**2
           R3=C43/(SSS+C3)**2
        W0=R0/(R0+R1+R2+R3)
	W1=R1/(R0+R1+R2+R3)
	W2=R2/(R0+R1+R2+R3)
	W3=R3/(R0+R1+R2+R3)
	ur(i+1)=(W0*q40+W1*q41+W2*q42+W3*q43)


        q40=(-1./4.*u(i-3)+13./12.*u(i-2)
     $       -23./12.*u(i-1)+25./12.*u(i))
	q41=(1./12.*u(i-2)-5./12.*u(i-1)
     $       +13./12.*u(i)+1./4.*u(i+1))
	q42=(-1./12.*u(i-1)+7./12.*u(i)
     $       +7./12.*u(i+1)-1./12.*u(i+2))
	q43=(1./4.*u(i)+13./12.*u(i+1)
     $       -5./12.*u(i+2)+1./12.*u(i+3))

	  C0=u(i-3)*(547*u(i-3)-3882*u(i-2)
     $                  +4642*u(i-1)-1854*u(i))
     $      +u(i-2)*(7043*u(i-2)-17246*u(i-1)
     $                  +7042*u(i))
     $      +u(i-1)*(11003*u(i-1)-9402*u(i))
     $      +2107*u(i)**2
	  C1=u(i-2)*(267*u(i-2)-1642*u(i-1)
     $                  +1602*u(i)-494*u(i+1))
     $      +u(i-1)*(2843*u(i-1)-5966*u(i)
     $                  +1922*u(i+1))
     $      +u(i)*(3443*u(i)-2522*u(i+1))
     $      +547*u(i+1)**2
	  C2=u(i-1)*(547*u(i-1)-2522*u(i)
     $                  +1922*u(i+1)-494*u(i+2))
     $      +u(i)*(3443*u(i)-5966*u(i+1)
     $                  +1602*u(i+2))
     $      +u(i+1)*(2843*u(i+1)-1642*u(i+2))
     $      +267*u(i+2)**2
	  C3=u(i)*(2107*u(i)-9402*u(i+1)
     $                  +7042*u(i+2)-1854*u(i+3))
     $      +u(i+1)*(11003*u(i+1)-17246*u(i+2)
     $                  +4642*u(i+3))
     $      +u(i+2)*(7043*u(i+2)-3882*u(i+3))
     $      +547*u(i+3)**2

           R0=C40/(SSS+C0)**2
           R1=C41/(SSS+C1)**2
           R2=C42/(SSS+C2)**2
           R3=C43/(SSS+C3)**2

        W0=R0/(R0+R1+R2+R3)
	W1=R1/(R0+R1+R2+R3)
	W2=R2/(R0+R1+R2+R3)
	W3=R3/(R0+R1+R2+R3)

	ul(i+1)=(W0*q40+W1*q41+W2*q42+W3*q43)
	end do

c-----------
        case (81)  
c-------------------------weno7(new1)---------------------
       first=4
	C40=1./35.
	C41=12./35.
	C42=18./35.
        C43=4./35.
       do i = first+1-1,imax-first-1

               q40=(-1./4.*u(i+4)+13./12.*u(i+3)
     $              -23./12.*u(i+2)+25./12.*u(i+1))
	       q41=(1./12.*u(i+3)-5./12.*u(i+2)
     $              +13./12.*u(i+1)+1./4.*u(i))
	       q42=(-1./12.*u(i+2)+7./12.*u(i+1)
     $               +7./12.*u(i)-1./12.*u(i-1))
	       q43=(1./4.*u(i+1)+13./12.*u(i)
     $              -5./12.*u(i-1)+1./12.*u(i-2))
	  C0=(-u(i+4)+3.*u(i+3)-3.*u(i+2)
     $         +u(i+1))**2/36.
     $       +(-u(i+4)+4.*u(i+3)-5.*u(i+2)
     $         +2.*u(i+1))**2/4.
     $       +(-2.*u(i+4)+9.*u(i+3)-18*u(i+2)
     $          +11.*u(i+1))**2/36.
	  C1=(-u(i+3)+3.*u(i+2)-3.*u(i+1)
     $         +u(i))**2/36.
     $       +(u(i+2)-2.*u(i+1)
     $         +u(i))**2/4.
     $       +(u(i+3)-6.*u(i+2)+3.*u(i+1)
     $          +2.*u(i))**2/36.
	  C2=(-u(i+2)+3.*u(i+1)-3.*u(i)
     $         +u(i-1))**2/36.
     $       +(u(i+2)-2.*u(i+1)+u(i))**2/4.
     $       +(-2.*u(i+2)-3.*u(i+1)+6.*u(i)
     $          -u(i-1))**2/36.
	  C3=(-u(i+1)+3.*u(i)-3.*u(i-1)
     $         +u(i-2))**2/36.
     $       +(2.*u(i+1)-5.*u(i)+4.*u(i-1)
     $          -3.*u(i-2))**2/4.
     $       +(-11.*u(i+1)+18.*u(i)-9.*u(i-1)
     $          +2.*u(i-2))**2/36.

	   R0=C40/(SSS+C0)**2
           R1=C41/(SSS+C1)**2
           R2=C42/(SSS+C2)**2
           R3=C43/(SSS+C3)**2
        W0=R0/(R0+R1+R2+R3)
	W1=R1/(R0+R1+R2+R3)
	W2=R2/(R0+R1+R2+R3)
	W3=R3/(R0+R1+R2+R3)
	ur(i+1)=(W0*q40+W1*q41+W2*q42+W3*q43)


        q40=(-1./4.*u(i-3)+13./12.*u(i-2)
     $       -23./12.*u(i-1)+25./12.*u(i))
	q41=(1./12.*u(i-2)-5./12.*u(i-1)
     $       +13./12.*u(i)+1./4.*u(i+1))
	q42=(-1./12.*u(i-1)+7./12.*u(i)
     $       +7./12.*u(i+1)-1./12.*u(i+2))
	q43=(1./4.*u(i)+13./12.*u(i+1)
     $       -5./12.*u(i+2)+1./12.*u(i+3))
	  C0=(-u(i-3)+3.*u(i-2)-3.*u(i-1)
     $         +u(i))**2/36.
     $       +(-u(i-3)+4.*u(i-2)-5.*u(i-1)
     $         +2.*u(i))**2/4.
     $       +(-2.*u(i-3)+9.*u(i-2)-18*u(i-1)
     $          +11.*u(i))**2/36.
	  C1=(-u(i-2)+3.*u(i-1)-3.*u(i)
     $         +u(i+1))**2/36.
     $       +(u(i-1)-2.*u(i)+u(i+1))**2/4.
     $       +(u(i-2)-6.*u(i-1)+3.*u(i)
     $          +2.*u(i+1))**2/36.
	  C2=(-u(i-1)+3.*u(i)-3.*u(i+1)
     $         +u(i+2))**2/36.
     $       +(u(i-1)-2.*u(i)+u(i+1))**2/4.
     $       +(-2.*u(i-1)-3.*u(i)+6.*u(i+1)
     $          -u(i+2))**2/36.
	  C3=(-u(i)+3.*u(i+1)-3.*u(i+2)
     $         +u(i+3))**2/36.
     $       +(2.*u(i)-5.*u(i+1)
     $         +4.*u(i+2)-u(i+3))**2/4.
     $       +(-11.*u(i)+18.*u(i+1)-9.*u(i+2)
     $          +2.*u(i+3))**2/36.

           R0=C40/(SSS+C0)**2
           R1=C41/(SSS+C1)**2
           R2=C42/(SSS+C2)**2
           R3=C43/(SSS+C3)**2

        W0=R0/(R0+R1+R2+R3)
	W1=R1/(R0+R1+R2+R3)
	W2=R2/(R0+R1+R2+R3)
	W3=R3/(R0+R1+R2+R3)

	ul(i+1)=(W0*q40+W1*q41+W2*q42+W3*q43)

	end do

c-----------
        case (9)  
c-------------------------weno7(new2)---------------------
       first=4
	C40=1./35.
	C41=12./35.
	C42=18./35.
        C43=4./35.
       do i = first+1-1,imax-first-1
              q40=(-1./4.*u(i+4)+13./12.*u(i+3)
     $              -23./12.*u(i+2)+25./12.*u(i+1))
	       q41=(1./12.*u(i+3)-5./12.*u(i+2)
     $              +13./12.*u(i+1)+1./4.*u(i))
	       q42=(-1./12.*u(i+2)+7./12.*u(i+1)
     $               +7./12.*u(i)-1./12.*u(i-1))
	       q43=(1./4.*u(i+1)+13./12.*u(i)
     $              -5./12.*u(i-1)+1./12.*u(i-2))

	  C0=(-u(i+4)+4.*u(i+3)-5.*u(i+2)
     $         +2.*u(i+1))**2/4.
     $       +(-2.*u(i+4)+9.*u(i+3)-18*u(i+2)
     $          +11.*u(i+1))**2/36.
	  C1=(u(i+2)-2.*u(i+1)
     $         +u(i))**2/4.
     $       +(u(i+3)-6.*u(i+2)+3.*u(i+1)
     $          +2.*u(i))**2/36.
	  C2=(u(i+2)-2.*u(i+1)+u(i))**2/4.
     $       +(-2.*u(i+2)-3.*u(i+1)+6.*u(i)
     $          -u(i-1))**2/36.
	  C3=(2.*u(i+1)-5.*u(i)+4.*u(i-1)
     $          -3.*u(i-2))**2/4.
     $       +(-11.*u(i+1)+18.*u(i)-9.*u(i-1)
     $          +2.*u(i-2))**2/36.

	   R0=C40/(SSS+C0)**2
           R1=C41/(SSS+C1)**2
           R2=C42/(SSS+C2)**2
           R3=C43/(SSS+C3)**2
        W0=R0/(R0+R1+R2+R3)
	W1=R1/(R0+R1+R2+R3)
	W2=R2/(R0+R1+R2+R3)
	W3=R3/(R0+R1+R2+R3)
	ur(i+1)=(W0*q40+W1*q41+W2*q42+W3*q43)


        q40=(-1./4.*u(i-3)+13./12.*u(i-2)
     $       -23./12.*u(i-1)+25./12.*u(i))
	q41=(1./12.*u(i-2)-5./12.*u(i-1)
     $       +13./12.*u(i)+1./4.*u(i+1))
	q42=(-1./12.*u(i-1)+7./12.*u(i)
     $       +7./12.*u(i+1)-1./12.*u(i+2))
	q43=(1./4.*u(i)+13./12.*u(i+1)
     $       -5./12.*u(i+2)+1./12.*u(i+3))


	  C0=(-u(i-3)+4.*u(i-2)-5.*u(i-1)
     $         +2.*u(i))**2/4.
     $       +(-2.*u(i-3)+9.*u(i-2)-18*u(i-1)
     $          +11.*u(i))**2/36.
	  C1=(u(i-1)-2.*u(i)+u(i+1))**2/4.
     $       +(u(i-2)-6.*u(i-1)+3.*u(i)
     $          +2.*u(i+1))**2/36.
	  C2=(u(i-1)-2.*u(i)+u(i+1))**2/4.
     $       +(-2.*u(i-1)-3.*u(i)+6.*u(i+1)
     $          -u(i+2))**2/36.
	  C3=(2.*u(i)-5.*u(i+1)
     $         +4.*u(i+2)-u(i+3))**2/4.
     $       +(-11.*u(i)+18.*u(i+1)-9.*u(i+2)
     $          +2.*u(i+3))**2/36.


           R0=C40/(SSS+C0)**2
           R1=C41/(SSS+C1)**2
           R2=C42/(SSS+C2)**2
           R3=C43/(SSS+C3)**2

        W0=R0/(R0+R1+R2+R3)
	W1=R1/(R0+R1+R2+R3)
	W2=R2/(R0+R1+R2+R3)
	W3=R3/(R0+R1+R2+R3)

	ul(i+1)=(W0*q40+W1*q41+W2*q42+W3*q43)

	end do

c-----------
        case (10)  
c-------------------------weno9(shu's)---------------------
       first=5
	C50=1./126.
	C51=10./63.
	C52=10./21.
        C53=20./63.
        C54=5./126.
       do i = first+1-1,imax-first-1
        q50=(1./5.*u(i+5)-21./20.*u(i+4)
     $      +137./60.*u(i+3)
     $       -163./60.*u(i+2)+137./60.*u(i+1))
        q51=(-1./20.*u(i+4)+17./60.*u(i+3)
     $       -43./60.*u(i+2)+77./60.*u(i+1)
     $       +1./5.*u(i))
	q52=(1./30.*u(i+3)-13./60.*u(i+2)
     $       +47./60.*u(i+1)+9./20.*u(i)
     $       -1./20.*u(i-1))
	q53=(-1./20.*u(i+2)+9./20.*u(i+1)
     $       +47./60.*u(i)-13./60.*u(i-1)
     $       +1./30.*u(i-2))
	q54=(1./5.*u(i+1)+77./60.*u(i)
     $       -43./60.*u(i-1)+17./60.*u(i-2)
     $       -1./20*u(i-3))

	       C0=u(i+5)*(22658*u(i+5)-208501*u(i+4)
     $                         +364863*u(i+3)
     $                         -288007*u(i+2)+86329*u(i+1))
     $           +u(i+4)*(482963*u(i+4)-1704396*u(i+3)
     $                         +1358458*u(i+2)
     $                         -411487*u(i+1))
     $           +u(i+3)*(1521393*u(i+3)-2462076*u(i+2)
     $                         +758823*u(i+1))
     $           +u(i+2)*(1020563*u(i+2)-649501*u(i+1))
     $                         +107918*u(i+1)**2
	       C1=u(i+4)*(6908*u(i+4)-60871*u(i+3)
     $                         +99213*u(i+2)-70237*u(i+1)
     $                         +18079*u(i))
     $           +u(i+3)*(138563*u(i+3)-464976*u(i+2)
     $                         +337018*u(i+1)-88297*u(i))
     $           +u(i+2)*(406293*u(i+2)-611976*u(i+1)
     $                         +165153*u(i))
     $           +u(i+1)*(242723*u(i+1)-140251*u(i))
     $                         +22658*u(i)**2
	       C2=u(i+3)*(6908*u(i+3)-51001*u(i+2)
     $                         +67923*u(i+1)-38947*u(i)
     $                          +8209*u(i-1))
     $           +u(i+2)*(104963*u(i+2)-299076*u(i+1)
     $                         +179098*u(i)-38947*u(i-1))
     $           +u(i+1)*(231153*u(i+1)-299076*u(i)
     $                         +67923*u(i-1))
     $           +u(i)*(104963*u(i)-51001*u(i-1))
     $                         +6908*u(i-1)**2
	       C3=u(i+2)*(22658*u(i+2)-140251*u(i+1)
     $                         +165153*u(i)-88297*u(i-1)
     $                         +18079*u(i-2))
     $           +u(i+1)*(242723*u(i+1)-611976*u(i)
     $                         +337018*u(i-1)-70237*u(i-2))
     $           +u(i)*(406293*u(i)-464976*u(i-1)
     $                       +99213*u(i-2))
     $           +u(i-1)*(138563*u(i-1)-60871*u(i-2))
     $                         +6908*u(i-2)**2
	       C4=u(i+1)*(107918*u(i+1)-649501*u(i)
     $                         +758823*u(i-1)-411487*u(i-2)
     $                         +86329*u(i-3))
     $           +u(i)*(1020563*u(i)-2462076*u(i-1)
     $                         +1358458*u(i-2)-288007*u(i-3))
     $           +u(i-1)*(1521393*u(i-1)-1704396*u(i-2)
     $                         +364863*u(i-3))
     $           +u(i-2)*(482963*u(i-2)-208501*u(i-3))
     $                         +22658*u(i-3)**2

	   R0=C50/(SSS+C0)**2
           R1=C51/(SSS+C1)**2
           R2=C52/(SSS+C2)**2
           R3=C53/(SSS+C3)**2
           R4=C54/(SSS+C4)**2
        W0=R0/(R0+R1+R2+R3+R4)
	W1=R1/(R0+R1+R2+R3+R4)
	W2=R2/(R0+R1+R2+R3+R4)
	W3=R3/(R0+R1+R2+R3+R4)
	W4=R4/(R0+R1+R2+R3+R4)
	ur(i+1)=(W0*q50+W1*q51+W2*q52+W3*q53+W4*q54)


        q50=(1./5.*u(i-4)-21./20.*u(i-3)
     $      +137./60.*u(i-2)
     $       -163./60.*u(i-1)+137./60.*u(i))
        q51=(-1./20.*u(i-3)+17./60.*u(i-2)
     $       -43./60.*u(i-1)+77./60.*u(i)
     $       +1./5.*u(i+1))
	q52=(1./30.*u(i-2)-13./60.*u(i-1)
     $       +47./60.*u(i)+9./20.*u(i+1)
     $       -1./20.*u(i+2))
	q53=(-1./20.*u(i-1)+9./20.*u(i)
     $       +47./60.*u(i+1)-13./60.*u(i+2)
     $       +1./30.*u(i+3))
	q54=(1./5.*u(i)+77./60.*u(i+1)
     $       -43./60.*u(i+2)+17./60.*u(i+3)
     $       -1./20*u(i+4))


	       C0=u(i-4)*(22658*u(i-4)-208501*u(i-3)
     $                         +364863*u(i-2)
     $                         -288007*u(i-1)+86329*u(i))
     $           +u(i-3)*(482963*u(i-3)-1704396*u(i-2)
     $                         +1358458*u(i-1)
     $                         -411487*u(i))
     $           +u(i-2)*(1521393*u(i-2)-2462076*u(i-1)
     $                         +758823*u(i))
     $           +u(i-1)*(1020563*u(i-1)-649501*u(i))
     $                         +107918*u(i)**2
	       C1=u(i-3)*(6908*u(i-3)-60871*u(i-2)
     $                         +99213*u(i-1)-70237*u(i)
     $                         +18079*u(i+1))
     $           +u(i-2)*(138563*u(i-2)-464976*u(i-1)
     $                         +337018*u(i)-88297*u(i+1))
     $           +u(i-1)*(406293*u(i-1)-611976*u(i)
     $                         +165153*u(i+1))
     $           +u(i)*(242723*u(i)-140251*u(i+1))
     $                         +22658*u(i+1)**2
	       C2=u(i-2)*(6908*u(i-2)-51001*u(i-1)
     $                         +67923*u(i)-38947*u(i+1)
     $                          +8209*u(i+2))
     $           +u(i-1)*(104963*u(i-1)-299076*u(i)
     $                         +179098*u(i+1)-38947*u(i+2))
     $           +u(i)*(231153*u(i)-299076*u(i+1)
     $                         +67923*u(i+2))
     $           +u(i+1)*(104963*u(i+1)-51001*u(i+2))
     $                         +6908*u(i+2)**2
	       C3=u(i-1)*(22658*u(i-1)-140251*u(i)
     $                         +165153*u(i+1)-88297*u(i+2)
     $                         +18079*u(i+3))
     $           +u(i)*(242723*u(i)-611976*u(i+1)
     $                         +337018*u(i+2)-70237*u(i+3))
     $           +u(i+1)*(406293*u(i+1)-464976*u(i+2)
     $                       +99213*u(i+2))
     $           +u(i+2)*(138563*u(i+2)-60871*u(i+3))
     $                         +6908*u(i+3)**2
	       C4=u(i)*(107918*u(i)-649501*u(i+1)
     $                         +758823*u(i+2)-411487*u(i+3)
     $                         +86329*u(i+4))
     $           +u(i+1)*(1020563*u(i+1)-2462076*u(i+2)
     $                         +1358458*u(i+3)-288007*u(i+4))
     $           +u(i+2)*(1521393*u(i+2)-1704396*u(i+3)
     $                         +364863*u(i+4))
     $           +u(i+3)*(482963*u(i+3)-208501*u(i+4))
     $                         +22658*u(i+4)**2

           R0=C50/(SSS+C0)**2
           R1=C51/(SSS+C1)**2
           R2=C52/(SSS+C2)**2
           R3=C53/(SSS+C3)**2
           R4=C54/(SSS+C4)**2

        W0=R0/(R0+R1+R2+R3+R4)
	W1=R1/(R0+R1+R2+R3+R4)
	W2=R2/(R0+R1+R2+R3+R4)
	W3=R3/(R0+R1+R2+R3+R4)
	W4=R4/(R0+R1+R2+R3+R4)

	ul(i+1)=(W0*q50+W1*q51+W2*q52+W3*q53+W4*q54)



	end do

c-----------
        case (11)  
c-------------------------weno9(new)---------------------
       first=5
           DO II=1,4
              CC(II)=0.
           END DO

	C50=1./126.
	C51=10./63.
	C52=10./21.
        C53=20./63.
        C54=5./126.

	   A9(1,1)=1./4.
           B9(1,1)=-4./3.
           C9(1,1)=3.
           D9(1,1)=-4.
           E9(1,1)=25./12.

	   A9(1,2)=11./12.
           B9(1,2)=-14./3.
           C9(1,2)=19./2
           D9(1,2)=-26./3
           E9(1,2)=35./12.

	   A9(1,3)=3./2.
           B9(1,3)=-7.
           C9(1,3)=12.
           D9(1,3)=-9.
           E9(1,3)=5./2.

	   A9(1,4)=1.
           B9(1,4)=-4.
           C9(1,4)=6.
           D9(1,4)=-4.
           E9(1,4)=1.

	   A9(2,1)=-1./12.
           B9(2,1)=1./2.
           C9(2,1)=-3./2.
           D9(2,1)=5./6.
           E9(2,1)=1./4

	   A9(2,2)=-1./12.
           B9(2,2)=1./3.
           C9(2,2)=1./2.
           D9(2,2)=-5./3.
           E9(2,2)=11./12.

	   A9(2,3)=1./2
           B9(2,3)=-3.
           C9(2,3)=6.
           D9(2,3)=-5.
           E9(2,3)=3./2.

	   A9(2,4)=1.
           B9(2,4)=-4.
           C9(2,4)=6.
           D9(2,4)=-4.
           E9(2,4)=1.

	   A9(3,1)=1./12.
           B9(3,1)=-2./3
           C9(3,1)=0.
           D9(3,1)=2./3.
           E9(3,1)=-1./12.

	   A9(3,2)=-1./12.
           B9(3,2)=4./3.
           C9(3,2)=-5./2.
           D9(3,2)=4./3.
           E9(3,2)=-1./12.

	   A9(3,3)=-1./2.
           B9(3,3)=1.
           C9(3,3)=0.
           D9(3,3)=-1.
           E9(3,3)=1./2.

	   A9(3,4)=1.
           B9(3,4)=-4.
           C9(3,4)=6.
           D9(3,4)=-4.
           E9(3,4)=1.

	   A9(4,1)=-1./4.
           B9(4,1)=-5./6.
           C9(4,1)=3./2.
           D9(4,1)=-1./2.
           E9(4,1)=1./12.

	   A9(4,2)=11./12.
           B9(4,2)=-5./3.
           C9(4,2)=1./2.
           D9(4,2)=1./3.
           E9(4,2)=-1./12.

	   A9(4,3)=-3./2.
           B9(4,3)=5.
           C9(4,3)=-6.
           D9(4,3)=3.
           E9(4,3)=-1./2.

	   A9(4,4)=1.
           B9(4,4)=-4.
           C9(4,4)=6.
           D9(4,4)=-4.
           E9(4,4)=1.

	   A9(5,1)=-25./12.
           B9(5,1)=4.
           C9(5,1)=-3.
           D9(5,1)=4./3.
           E9(5,1)=-1./4.

	   A9(5,2)=35./12
           B9(5,2)=-26./3.
           C9(5,2)=19./2.
           D9(5,2)=-14./3.
           E9(5,2)=11./12.

	   A9(5,3)=-5./2.
           B9(5,3)=9.
           C9(5,3)=-12.
           D9(5,3)=7.
           E9(5,3)=-3./2.

	   A9(5,4)=1.
           B9(5,4)=-4.
           C9(5,4)=6.
           D9(5,4)=-4.
           E9(5,4)=1.

       do i = first+1-1,imax-first-1
        q50=(1./5.*u(i+5)-21./20.*u(i+4)
     $      +137./60.*u(i+3)
     $       -163./60.*u(i+2)+137./60.*u(i+1))
        q51=(-1./20.*u(i+4)+17./60.*u(i+3)
     $       -43./60.*u(i+2)+77./60.*u(i+1)
     $       +1./5.*u(i))
	q52=(1./30.*u(i+3)-13./60.*u(i+2)
     $       +47./60.*u(i+1)+9./20.*u(i)
     $       -1./20.*u(i-1))
	q53=(-1./20.*u(i+2)+9./20.*u(i+1)
     $       +47./60.*u(i)-13./60.*u(i-1)
     $       +1./30.*u(i-2))
	q54=(1./5.*u(i+1)+77./60.*u(i)
     $       -43./60.*u(i-1)+17./60.*u(i-2)
     $       -1./20*u(i-3))
        DO II=1,5
           DO J=1,3
              CC(J)= A9(II,J)*u(i+6-II)
     $              +B9(II,J)*u(i+5-II)
     $              +C9(II,J)*u(i+4-II)
     $              +D9(II,J)*u(i+3-II)
     $              +E9(II,J)*u(i+2-II)
           END DO
           CCC(II)=CC(1)**2+CC(2)**2/4.+CC(3)**2/36.+CC(4)**2/576.
        END DO
        C0=CCC(1)
        C1=CCC(2)
        C2=CCC(3)
        C3=CCC(4)
        C4=CCC(5)


	   R0=C50/(SSS+C0)**2
           R1=C51/(SSS+C1)**2
           R2=C52/(SSS+C2)**2
           R3=C53/(SSS+C3)**2
           R4=C54/(SSS+C4)**2
        W0=R0/(R0+R1+R2+R3+R4)
	W1=R1/(R0+R1+R2+R3+R4)
	W2=R2/(R0+R1+R2+R3+R4)
	W3=R3/(R0+R1+R2+R3+R4)
	W4=R4/(R0+R1+R2+R3+R4)
	ur(i+1)=(W0*q50+W1*q51+W2*q52+W3*q53+W4*q54)


        q50=(1./5.*u(i-4)-21./20.*u(i-3)
     $      +137./60.*u(i-2)
     $       -163./60.*u(i-1)+137./60.*u(i))
        q51=(-1./20.*u(i-3)+17./60.*u(i-2)
     $       -43./60.*u(i-1)+77./60.*u(i)
     $       +1./5.*u(i+1))
	q52=(1./30.*u(i-2)-13./60.*u(i-1)
     $       +47./60.*u(i)+9./20.*u(i+1)
     $       -1./20.*u(i+2))
	q53=(-1./20.*u(i-1)+9./20.*u(i)
     $       +47./60.*u(i+1)-13./60.*u(i+2)
     $       +1./30.*u(i+3))
	q54=(1./5.*u(i)+77./60.*u(i+1)
     $       -43./60.*u(i+2)+17./60.*u(i+3)
     $       -1./20*u(i+4))

        DO II=1,5
           DO J=1,3
              CC(J)= A9(II,J)*u(i-5+II)
     $           +B9(II,J)*u(i-4+II)
     $           +C9(II,J)*u(i-3+II)
     $           +D9(II,J)*u(i-2+II)
     $           +E9(II,J)*u(i-1+II)
           END DO
           CCC(II)=CC(1)**2+CC(2)**2/4.+CC(3)**2/36.+CC(4)**2/576.
        END DO
        C0=CCC(1)
        C1=CCC(2)
        C2=CCC(3)
        C3=CCC(4)
        C4=CCC(5)


           R0=C50/(SSS+C0)**2
           R1=C51/(SSS+C1)**2
           R2=C52/(SSS+C2)**2
           R3=C53/(SSS+C3)**2
           R4=C54/(SSS+C4)**2

        W0=R0/(R0+R1+R2+R3+R4)
	W1=R1/(R0+R1+R2+R3+R4)
	W2=R2/(R0+R1+R2+R3+R4)
	W3=R3/(R0+R1+R2+R3+R4)
	W4=R4/(R0+R1+R2+R3+R4)

	ul(i+1)=(W0*q50+W1*q51+W2*q52+W3*q53+W4*q54)



	end do


        case (106)  
c-----------------Borges weno5---------------------
       first=0
       eps=1e-20
       do i = first+1,imax-first
	   is0=b0*(du(i-1)-du(i-2))**2+0.25*(3.*du(i-1)-du(i-2))**2
           is1=b0*(du(i)-du(i-1))**2+0.25*(du(i)+du(i-1))**2
	   is2=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)-3.*du(i))**2
           tao5=abs(is0-is2)
c	   a0=c20/(eps+is0)**2
c	   a1=c21/(eps+is1)**2
c	   a2=c22/(eps+is2)**2
           a0=c20*(1+tao5/(is0+eps))
           a1=c21*(1+tao5/(is1+eps))
           a2=c22*(1+tao5/(is2+eps))
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i-1)+0.5*du(i-1)+0.333333*(du(i-1)-du(i-2))
	   q1=u(i-1)+0.5*du(i-1)+0.333333*(du(i)-du(i-1))
	   q2=u(i-1)+0.5*du(i)-0.166667*(du(i+1)-du(i))
	   ul(i)=w0*q0+w1*q1+w2*q2

	   is0=b0*(du(i)-du(i-1))**2+0.25*(3.*du(i)-du(i-1))**2
	   is1=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)+du(i))**2
	   is2=b0*(du(i+2)-du(i+1))**2+0.25*(du(i+2)-3.*du(i+1))**2
           tao5=abs(is0-is2)
c	   a0=c22/(eps+is0)**2
c	   a1=c21/(eps+is1)**2
c	   a2=c20/(eps+is2)**2
           a0=c22*(1+tao5/(is0+eps))
           a1=c21*(1+tao5/(is1+eps))
           a2=c20*(1+tao5/(is2+eps))
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i)-0.5*du(i)-0.166667*(du(i)-du(i-1))
	   q1=u(i)-0.5*du(i+1)+0.333333*(du(i+1)-du(i))
	   q2=u(i)-0.5*du(i+1)+0.333333*(du(i+2)-du(i+1))
	   ur(i)=w0*q0+w1*q1+w2*q2
	end do

        end select

c-----------------------------------------------------------
c--------------------near boundary for >7-weno--------------
        if (iside.ge.7.and.iside.le.9) then
           first=4
        else if (iside.ge.10.and.iside.le.11) then
          first=5
        else
          first=0
        end if
       if (iside.ge.7) then
       do i = 1,first
	   is0=b0*(du(i-1)-du(i-2))**2+0.25*(3.*du(i-1)-du(i-2))**2
           is1=b0*(du(i)-du(i-1))**2+0.25*(du(i)+du(i-1))**2
	   is2=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)-3.*du(i))**2
	   a0=c20/(eps+is0)**2
	   a1=c21/(eps+is1)**2
	   a2=c22/(eps+is2)**2
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i-1)+0.5*du(i-1)+0.333333*(du(i-1)-du(i-2))
	   q1=u(i-1)+0.5*du(i-1)+0.333333*(du(i)-du(i-1))
	   q2=u(i-1)+0.5*du(i)-0.166667*(du(i+1)-du(i))
	   ul(i)=w0*q0+w1*q1+w2*q2

	   is0=b0*(du(i)-du(i-1))**2+0.25*(3.*du(i)-du(i-1))**2
	   is1=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)+du(i))**2
	   is2=b0*(du(i+2)-du(i+1))**2+0.25*(du(i+2)-3.*du(i+1))**2
	   a0=c22/(eps+is0)**2
	   a1=c21/(eps+is1)**2
	   a2=c20/(eps+is2)**2
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i)-0.5*du(i)-0.166667*(du(i)-du(i-1))
	   q1=u(i)-0.5*du(i+1)+0.333333*(du(i+1)-du(i))
	   q2=u(i)-0.5*du(i+1)+0.333333*(du(i+2)-du(i+1))
	   ur(i)=w0*q0+w1*q1+w2*q2
	end do
       do i = imax-first+1,imax
	   is0=b0*(du(i-1)-du(i-2))**2+0.25*(3.*du(i-1)-du(i-2))**2
           is1=b0*(du(i)-du(i-1))**2+0.25*(du(i)+du(i-1))**2
	   is2=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)-3.*du(i))**2
	   a0=c20/(eps+is0)**2
	   a1=c21/(eps+is1)**2
	   a2=c22/(eps+is2)**2
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i-1)+0.5*du(i-1)+0.333333*(du(i-1)-du(i-2))
	   q1=u(i-1)+0.5*du(i-1)+0.333333*(du(i)-du(i-1))
	   q2=u(i-1)+0.5*du(i)-0.166667*(du(i+1)-du(i))
	   ul(i)=w0*q0+w1*q1+w2*q2

	   is0=b0*(du(i)-du(i-1))**2+0.25*(3.*du(i)-du(i-1))**2
	   is1=b0*(du(i+1)-du(i))**2+0.25*(du(i+1)+du(i))**2
	   is2=b0*(du(i+2)-du(i+1))**2+0.25*(du(i+2)-3.*du(i+1))**2
	   a0=c22/(eps+is0)**2
	   a1=c21/(eps+is1)**2
	   a2=c20/(eps+is2)**2
	   w0=a0/(a0+a1+a2)
	   w1=a1/(a0+a1+a2)
	   w2=a2/(a0+a1+a2)
	   q0=u(i)-0.5*du(i)-0.166667*(du(i)-du(i-1))
	   q1=u(i)-0.5*du(i+1)+0.333333*(du(i+1)-du(i))
	   q2=u(i)-0.5*du(i+1)+0.333333*(du(i+2)-du(i+1))
	   ur(i)=w0*q0+w1*q1+w2*q2
	end do
      end if

c------------------------------------------------------------
c--------------------boundary--------------------------------
        if (iside.ge.5 .and. iside.le.11) then
           first=3
c           first=2
        else
           first=2
           first=3
        end if
        if (bc_lower.ne.7.and.bc_lower.ne.8
     $      .and.bc_lower.ne.10.and.bc_lower.ne.20) then
           do i=1,first
              ul(i)=u(i-1)
              ur(i)=u(i)
           end do
        end if
        if (bc_upper.ne.7.and.bc_upper.ne.8
     $      .and.bc_upper.ne.10.and.bc_upper.ne.20) then
           do i=imax-first+1,imax
              ul(i)=u(i-1)
              ur(i)=u(i)
           end do
        end if
c------------------end of boundary---------------------------
c------------------------------------------------------------

        if (wall_order.eq.1) then
           return
        elseif (wall_order.eq.2) then  ! 2nd-order wall accuracy  
c-----------------wall
           i=imax-1
           if (bc_upper.eq.3.or.bc_upper.eq.19 .or.
     $              (bc_upper.ge.101 .and. bc_upper.le.110)) then
               ul(i+1) = (3.*u(i)-u(i-1))/2.
           end if
           i=1
           if (bc_lower.eq.3.or.bc_lower.eq.19 .or.
     $              (bc_lower.ge.101 .and. bc_lower.le.110)) then
             ur(i) = (3.*u(i)-u(i+1))/2.
           end if
c--------------------end of wall---------------------------
c----------------------near wall----------------------------------
           i=imax-2
           if (bc_upper.eq.3.or.bc_upper.eq.19 .or.
     $              (bc_upper.ge.101 .and. bc_upper.le.110)) then
               ul(i+1) = (3.*u(i)-u(i-1))/2.
               ur(i+1) = (3.*u(i+1)-u(i))/2.
           end if

           i=2
           if (bc_lower.eq.3.or.bc_lower.eq.19 .or.
     $              (bc_lower.ge.101 .and. bc_lower.le.110)) then
               ul(i) = (3.*u(i-1)-u(i))/2.
               ur(i) = (3.*u(i)-u(i+1))/2.
           end if
c--------------------end of near wall---------------------------
c----------------------next-near wall----------------------------------
           i=imax-3
           if (bc_upper.eq.3.or.bc_upper.eq.19 .or.
     $              (bc_upper.ge.101 .and. bc_upper.le.110)) then
c-----------------------third order------------
               ul(i+1) = (2.*u(i+1)+5.*u(i)-u(i-1))/6.
               ur(i+1) = (-u(i+2)+5.*u(i+1)+2.*u(i))/6.
           end if

           i=3
           if (bc_lower.eq.3.or.bc_lower.eq.19 .or.
     $              (bc_lower.ge.101 .and. bc_lower.le.110)) then
c-----------------------third order-----------
               ul(i) = (-u(i-2)+5.*u(i-1)+2.*u(i))/6.
               ur(i) = (2.*u(i-1)+5.*u(i)-u(i+1))/6.
           end if
c--------------------end of next-near wall---------------------------
        elseif (wall_order.eq.3) then  ! 3rd-order wall accuracy 
c  note: upper: h_(i+1/2)
c      : lower: h_(i-1/2)
c
c----------------------wall----------------------------------
c----------------------wall boundary-----------
           i=imax-1
           if (bc_upper.eq.3.or.bc_upper.eq.19 .or.
     $              (bc_upper.ge.101 .and. bc_upper.le.110)) then
               ul(i+1) = (11.*u(i)-7.*u(i-1)+2.*u(i-2))/6.
           end if

           i=1
           if (bc_lower.eq.3.or.bc_lower.eq.19 .or.
     $              (bc_lower.ge.101 .and. bc_lower.le.110)) then
               ur(i) = (2.*u(i+2)-7.*u(i+1)+11.*u(i))/6.
           end if
c--------------------end of wall---------------------------
c----------------------near wall----------------------------------
           i=imax-2
           if (bc_upper.eq.3.or.bc_upper.eq.19 .or.
     $              (bc_upper.ge.101 .and. bc_upper.le.110)) then
               ul(i+1) = (11.*u(i)-7.*u(i-1)+2.*u(i-2))/6.
               ur(i+1) = (2.*u(i+1)+5.*u(i)-u(i-1))/6.
           end if

           i=2
           if (bc_lower.eq.3.or.bc_lower.eq.19 .or.
     $              (bc_lower.ge.101 .and. bc_lower.le.110)) then
               ul(i) = (2.*u(i-1)+5.*u(i)-u(i+1))/6.
               ur(i) = (11.*u(i)-7.*u(i+1)+2.*u(i+2))/6.
           end if
c--------------------end of near wall---------------------------
           if (first.eq.3) then
c----------------------next-near wall----------------------------------
               i=imax-3
               if (bc_upper.eq.3.or.bc_upper.eq.19 .or.
     $              (bc_upper.ge.101 .and. bc_upper.le.110)) then
c-----------------------third order------------
c                  ul(i+1) = (2.*u(i+1)+5.*u(i)-u(i-1))/6.
c                  ur(i+1) = (-u(i+2)+5.*u(i+1)+2.*u(i))/6.
c-------------above do not work CFL=20 for nozzle
                  ul(i+1) = (11.*u(i)-7.*u(i-1)+2.*u(i-2))/6.
                  ur(i+1) = (2.*u(i+1)+5.*u(i)-u(i-1))/6.
c-------------above do not work CFL=20 for nozzle
c-----------------------fourth order------------
c                  ul(i+1) = (3.*u(i+1)+13.*u(i)-5.*u(i-1)+u(i-2))/12.
c                  ur(i+1) = (-u(i+2)+7.*u(i+1)+7.*u(i)-u(i-1))/12.
c----------------------------------------------------
               end if
               i=3
               if (bc_lower.eq.3.or.bc_lower.eq.19 .or.
     $              (bc_lower.ge.101 .and. bc_lower.le.110)) then
c-----------------------third order-----------
c                    ul(i) = (-u(i-2)+5.*u(i-1)+2.*u(i))/6.
c                    ur(i) = (2.*u(i-1)+5.*u(i)-u(i+1))/6.
c-----------above is also right for plate
                    ul(i) = (2.*u(i-1)+5.*u(i)-u(i+1))/6.
                    ur(i) = (11.*u(i)-7.*u(i+1)+2.*u(i+2))/6.
c-----------------------fourth order(bad)-----------
c                    ul(i) = (-u(i-2)+7.*u(i-1)+7.*u(i)-u(i+1))/12.
c                    ur(i) = (3.*u(i-1)+13.*u(i)-5.*u(i+1)+u(i+2))/12.
c---------------------------------------------
               end if
c--------------------end of next-near wall---------------------------
           end if
           return
        else       ! wall accuracy
           write(*,*)' The other accuracy order is not avaiable'
           stop
        end if




        end
