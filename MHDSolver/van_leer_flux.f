      subroutine van_leer_flux(GM,NL,VX1,VY1,VZ1,C_L,D_L,
     >                         U_L,V_L,W_L,C_R,D_R,U_R,V_R,W_R,
     >                         QT_L, QT_R, UPO)
c
c *** This is to calculate Van Leer flux vector splitting at subsonic
c
c ... implicit statement
c
      implicit none
c
c ... interface variables
c
      integer, intent(in)::
     $     NL                                ! number of equations
c
      double precision, intent(in)::
     $     GM,                               ! Specific heat ratio =1.4 for ideal gas
     $     VX1,VY1,VZ1,                      ! interface unit vector components
     $     U_L,V_L,W_L,                      ! LEFT Velocity components in x,y,z dir
     $     C_L,                              ! LEFT Speed of Sound
     $     D_L,                              ! LEFT DENSITY
     $     U_R,V_R,W_R,                      ! RIGHT Velocity components in x,y,z dir
     $     C_R,                              ! RIGHT Speed of Sound
     $     D_R,                              ! RIGHT DENSITY
     $     QT_L, QT_R                        ! GRID VELOCITIES AT LEFT AND RIGHT SIDES
c
      double precision, intent(out)::
     $     UPO(nl)                           ! Split flux and total flux
c
c ... local variables
c
      double precision      U1               ! total velocity
      double precision      U05,V05,W05      ! Velocity components in x,y,z dir
      double precision      C05              ! Speed of Sound
      double precision      QM1              ! Mach Number
      double precision      RA               ! =1 Flux+(Left), =-1 Flux-(Right)
      double precision      FMS              ! Van Leer mass flux
      double precision      D05              ! Density
      double precision      GM1              ! GM - 1
c
c *** Calculate Left flux
c
      GM1 = GM-1.d0
      RA=1.d0
      U05=U_L
      V05=V_L
      W05=W_L
      C05=C_L
      D05=D_L

      U1=VX1*U05+VY1*V05+VZ1*W05+QT_L
      QM1=U1/C05
      FMS=RA*D05*C05*(0.5d0*(QM1+RA))**2
      UPO(1)=FMS
      UPO(2)=FMS*(VX1*(-U1+RA*2.D0*C05)/GM+U05)
      UPO(3)=FMS*(VY1*(-U1+RA*2.D0*C05)/GM+V05)
      UPO(4)=FMS*(VZ1*(-U1+RA*2.D0*C05)/GM+W05)
      UPO(5)=FMS*((-GM1*U1**2+RA*2.D0*GM1*U1*C05
     $ +2.D0*C05**2)/(GM**2-1.D0)+0.5D0*(U05**2+V05**2+W05**2))

c..   Add the right flux

      RA=-1.0
      U05=U_R
      V05=V_R
      W05=W_R
      C05=C_R
      D05=D_R

      U1=VX1*U05+VY1*V05+VZ1*W05+QT_R
      QM1=U1/C05
      FMS=RA*D05*C05*(0.5d0*(QM1+RA))**2
      UPO(1)=UPO(1)+ FMS
      UPO(2)=UPO(2)+ FMS*(VX1*(-U1+RA*2.D0*C05)/GM+U05)
      UPO(3)=UPO(3)+ FMS*(VY1*(-U1+RA*2.D0*C05)/GM+V05)
      UPO(4)=UPO(4)+ FMS*(VZ1*(-U1+RA*2.D0*C05)/GM+W05)
      UPO(5)=UPO(5)+ FMS*((-GM1*U1**2+RA*2.D0*GM1*U1*C05
     $ +2.D0*C05**2)/(GM**2-1.D0)+0.5D0*(U05**2+V05**2+W05**2))

      RETURN
      END
