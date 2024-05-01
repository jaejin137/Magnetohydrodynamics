      subroutine generate_medium(choice, ii,nl,
     $     ilower,iupper,jlower, jupper, klower, kupper, 
     $     control)
c     generate initial condition for rse code
c
c     IMPLICIT STATEMENT
      use datain

      implicit none

      type (datain_type), intent(in)::control

      integer medium,i,j,k,n
      include "medium.f"
      double precision::
     $     poutlet, ptotal, ttotal, velinit,angl1,angl2

      character(len=1):: choice

      integer, intent(in)::
     $     ii,nl,
     $     ilower, iupper,
     $     jlower, jupper,
     $     klower, kupper

      double precision q(ilower:iupper,jlower:jupper,klower:kupper,nl)

      character*20:: 
     $         filename,pre_file_rstart

      double precision:: vtotal,temp1,vtemp
c     data statements
c--------------------------parameter transfer
      angl1=control%angl1
      angl2=control%angl2
      poutlet=control%poutlet
      ptotal=control%ptotal
      ttotal=control%ttotal
      velinit=control%velinit
cc

	write(*,*) ' Input filename:init_medium.dat'
c	READ(*,'(A)') filename
	filename='init_medium.dat'
	open(5,file=filename)

	read(5,*) medium,istate
        read(5,*)  UW
        read(5,*)  RHOWW
        read(5,*)  TW

	print *, 'medium 1-air 2-O2 3-H2 4-H2O 5-CO:',medium
c   istate  1-- ideal gas    2-real gas    3-  table
        print *,'Ref. Velocity:',UW
        print *,'Ref. density:',RHOWW
        print *,'Ref. temperature:',TW
        close(5)
	call cppp(medium)

c  read from table
	if (istate.eq.3) then
	   open(16, file='State_table.dat',form= 'unformatted')
 	   do j=1,401
              READ(16) t_t(j),t_v(j),p_c(j),t_mu(j),t_cp(j),
     $                v_c(j,1),v_c(j,2)
	      do i=1,401
                 READ(16) t_p(i,j),t_h(i,j)
	      end do
	   end do
	   close(16)
	endif

	if (istate.eq.4) then
	   sol=1200.0d0/UW
	end if
	if (istate.eq.5) then
 	   open(19,file='t_c_p_c.dat')
	   do i=1,101
	      read(19,*) t_cp(i),p_c(i)
	   end do
	   close(19)	
	   do i=1,101
	      p_c(i)=p_c(i)/(UW*UW*RHOWW)
	      t_cp(i)=t_cp(i)/TW
	   end do
	   sol=1200.0d0/UW
	   sog=350.0d0/UW
	   svl=1.0d0
	   svg=1000.d0
	end if

        if (choice.eq.'n') then
  	   call qstat(ptotal,ttotal,0.0,vtotal)
           write(*,*)'ptotal,ttotal,vtotal=',ptotal,ttotal,vtotal
           do k = klower,kupper
             do j = jlower, jupper
               do i = ilower,iupper
                 q(i,j,k,1) = poutlet
	         q(i,j,k,2) = velinit
                 q(i,j,k,3) = velinit*dtan(angl1)
                 q(i,j,k,4) = velinit*dtan(angl2)
             	 temp1=0.5*velinit*velinit	
 	         call htvp(ttotal,vtotal,ptotal,
     $		 temp1,q(i,j,k,5),q(i,j,k,1),vtemp)
               end do
             end do
           end do
           
           write(pre_file_rstart, '("pre_rstart",i3.3)') ii
           open(12, file=pre_file_rstart, form= 'unformatted')
           do n = 1, nl
             do k = klower,kupper
               do j = jlower, jupper
                 do i = ilower,iupper
                   write(12) q(i,j,k,n)
                 end do
               end do
             end do
           end do
           close(12)
           write(*,*)'precondition variables finished'
        end if

        return
	end


