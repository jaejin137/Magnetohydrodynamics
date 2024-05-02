      subroutine van_leer_matrix(GM,NL,VX1,VY1,VZ1,C_L,D_L,U_L,V_L,W_L,
     >                    E_L,C_R,D_R,U_R,V_R,W_R,E_R,TK_L,TK_R,
     >                    AL,AR,QT_L,QT_R)

c..   This is to calculate Van Leer flux vector splitting at subsonic

      implicit none

c..   input

      integer,         intent(in)::      NL               ! number of equations
      double precision,intent(in)::      GM               ! Specific heat ratio =1.4 for ideal gas      
      double precision,intent(in)::      VX1,VY1,VZ1      ! interface unit vector components
      double precision,intent(in)::      U_L,V_L,W_L      ! LEFT Velocity components in x,y,z dir
      double precision,intent(in)::      C_L              ! LEFT Speed of Sound
      double precision,intent(in)::      D_L              ! LEFT DENSITY
      double precision,intent(in)::      E_L              ! LEFT Total Energy
      double precision,intent(in)::      U_R,V_R,W_R      ! RIGHT Velocity components in x,y,z dir
      double precision,intent(in)::      C_R              ! RIGHT Speed of Sound
      double precision,intent(in)::      D_R              ! RIGHT DENSITY
      double precision,intent(in)::      E_R              ! Right Total Energy
      double precision,intent(in)::      TK_L, TK_R       ! turbulence viscous
      double precision,intent(in)::      QT_L, QT_R       ! GRID VELOCITIES AT LEFT AND RIGHT SIDES

c...  output

      double precision,intent(out)::      AL(nl,nl)       ! LEFT Jacobian Matrix
      double precision,intent(out)::      AR(nl,nl)       ! LEFT Jacobian Matrix

c..   local working variables

      double precision      U1               ! total velocity
      double precision      U05,V05,W05      ! Velocity components in x,y,z dir
      double precision      C05              ! Speed of Sound
      double precision      QM1              ! Mach Number
      double precision      RA               ! =1 Flux+(Left), =-1 Flux-(Right)
      double precision      FMS              ! Van Leer mass flux
      double precision      D05              ! Density
      double precision      E05              ! Total Energy
      integer               LTY,LTV
      double precision 
     >     BETA,ZM,YU,SA1,SA2,SA3,SA4,SA5,YM,RZ,AE,UVW

c...  Calculate Left Positive matrix

      BETA=GM-1.D0

C..   CALCULATE THE LEFT JACOBIAN MATRIX 

      RA=1.0
      U05=U_L
      V05=V_L
      W05=W_L
      C05=C_L
      D05=D_L
      E05=E_L

      U1=VX1*U05+VY1*V05+VZ1*W05+QT_L
      QM1=U1/C05
      FMS=RA*D05*C05*(0.5D0*(QM1+RA))**2

      IF(QM1.LE.1.D0) Then           ! subsonic

      ZM=D05**2*(U05**2+V05**2+W05**2)
      AL(1,1)=RA*1.D0/8.D0*GM*BETA*E05/D05/C05*(-QM1**2+1.D0)
      AL(1,2)=RA*1.D0/4.D0*(2.D0*VX1*(QM1+RA)
     $            +.5D0*GM*(GM-1.D0)*U05/C05*(QM1**2 -1.D0))
      AL(1,3)=RA*1.D0/4.D0*(2.D0*VY1*(QM1+RA)+.5*GM*(GM-1.D0)
     $                           *V05/C05*(QM1**2-1.D0))
      AL(1,4)=RA*1.D0/4.D0*(2.D0*VZ1*(QM1+RA)+.5D0*GM*(GM-1.D0)
     $                  *W05/C05*(QM1**2-1.D0))
      AL(1,5)=RA*1.D0/8.D0*GM*BETA/C05*(-QM1**2+1.D0)
      YU=VX1*(-U1+RA*2.D0*C05)/GM+U05
      SA1=(U1/D05+RA*GM*BETA*(-E05/D05**2+ZM/D05**3)/C05)/GM
      SA2=(-VX1/D05-RA*GM*BETA*U05/D05/C05)/GM
      SA3=(-VY1/D05-RA*GM*BETA*V05/D05/C05)/GM
      SA4=(-VZ1/D05-RA*GM*BETA*W05/D05/C05)/GM
      AL(2,1)=AL(1,1)*YU+FMS*(VX1*SA1-U05/D05)
      AL(2,2)=AL(1,2)*YU+FMS*(VX1*SA2+1.D0/D05)
      AL(2,3)=AL(1,3)*YU+FMS*VX1*SA3
      AL(2,4)=AL(1,4)*YU+FMS*VX1*SA4
      AL(2,5)=AL(1,5)*YU+RA*FMS*VX1*BETA/D05/C05
      YM=VY1*(-U1+RA*2.D0*C05)/GM+V05
      AL(3,1)=AL(1,1)*YM+FMS*(VY1*SA1-V05/D05)
      AL(3,2)=AL(1,2)*YM+FMS*VY1*SA2
      AL(3,3)=AL(1,3)*YM+FMS*(VY1*SA3+1.D0/D05)
      AL(3,4)=AL(1,4)*YM+FMS*VY1*SA4
      AL(3,5)=AL(1,5)*YM+FMS*VY1*RA*BETA/D05/C05
      RZ=VZ1*(-U1+RA*2.D0*C05)/GM+W05
      AL(4,1)=AL(1,1)*RZ+FMS*(VZ1*SA1-W05/D05)
      AL(4,2)=AL(1,2)*RZ+FMS*VZ1*SA2
      AL(4,3)=AL(1,3)*RZ+FMS*VZ1*SA3
      AL(4,4)=AL(1,4)*RZ+FMS*(VZ1*SA4+1.D0/D05)
      AL(4,5)=AL(1,5)*RZ+FMS*VZ1*RA*BETA/D05/C05
      AE=(-BETA*U1**2+RA*2.D0*BETA*U1*C05+2.D0*C05**2)/(GM**2-1.D0)
     *   +.5D0*(U05**2+V05**2+W05**2)
      AL(5,1)=AL(1,1)*AE+FMS*((2.D0*BETA*U1**2/D05+RA*GM*BETA**2*
     *  U1*(-3.D0*E05+2.D0*ZM/D05)/D05**2/C05+2.D0*GM*BETA/D05**2*(-E05
     *   +ZM/D05))/(GM**2-1.D0)-ZM/D05**3)
      AL(5,2)=AL(1,2)*AE+FMS*((-2.D0*BETA*U1*VX1/D05+RA*GM*BETA**2*(
     * 2.D0*E05*VX1/D05**2-U05*U1/D05-ZM*VX1/D05**3)/C05-2.D0*GM*BETA*
     *   U05/D05)/(GM**2-1.D0)+U05/D05)
      AL(5,3)=AL(1,3)*AE+FMS*((-2.D0*BETA*U1*VY1/D05+RA*GM*BETA**2*(
     * 2.D0*E05*VY1/D05**2-U1*V05/D05-ZM*VY1/D05**3)/C05-2.D0*GM*BETA*
     *   V05/D05)/(GM**2-1.D0)+V05/D05)
      AL(5,4)=AL(1,4)*AE+FMS*((-2.D0*BETA*U1*VZ1/D05+RA*GM*BETA**2*(
     * 2.D0*E05*VZ1/D05**2-U1*W05/D05-ZM*VZ1/D05**3)/C05-2.D0*GM*BETA*
     *   W05/D05)/(GM**2-1.D0)+W05/D05)
      AL(5,5)=AL(1,5)*AE+FMS*((RA*GM*BETA**2*U1/D05/C05+2.D0*GM*BETA/D05
     *  )/(GM**2-1.D0))
C
        IF (NL.EQ.6) THEN
          AL(1,6)=0.D0
          AL(2,6)=0.D0
          AL(3,6)=0.D0
          AL(4,6)=0.D0
          AL(5,6)=0.D0
          AL(6,1)=0.25D0*RA*TK_L*((0.5D0*GM*BETA*E05/(C05*D05)-C05)*
     $            (1.D0-QM1**2)-2.D0*C05*QM1*(QM1+RA))
          AL(6,2)=0.25D0*RA*TK_L*(2.D0*VX1*(QM1+RA)+
     $            .5D0*GM*BETA*U05/C05*(QM1**2-1.D0))
          AL(6,3)=0.25D0*RA*TK_L*(2.D0*VY1*(QM1+RA)+
     $            .5d0*GM*BETA*V05/C05*(QM1**2-1.D0))
          AL(6,4)=0.25D0*RA*TK_L*(2.D0*VZ1*(QM1+RA)+
     $            .5D0*GM*BETA*W05/C05*(QM1**2-1.D0))
          AL(6,5)=0.125D0*RA*TK_L*GM*BETA/C05*(-QM1**2+1.D0)
          AL(6,6)=0.25D0*RA*C05*(QM1+RA)**2
        END IF
c
      ELSE IF(QM1.GT.1) THEN

      UVW=U05**2+V05**2+W05**2
      AL(1,1)=0.D0
      AL(1,2)=VX1
      AL(1,3)=VY1
      AL(1,4)=VZ1
      AL(1,5)=0.D0
      AL(2,1)=-(U05**2*VX1+U05*V05*VY1+U05*W05*VZ1)+BETA*VX1/2.D0*UVW
      AL(2,2)=U05*VX1*(3.D0-GM)+V05*VY1+W05*VZ1
      AL(2,3)=U05*VY1-BETA*V05*VX1
      AL(2,4)=U05*VZ1-BETA*W05*VX1
      AL(2,5)=BETA*VX1
      AL(3,1)=-(U05*V05*VX1+V05**2*VY1+V05*W05*VZ1)+VY1*BETA/2.D0*UVW
      AL(3,2)=V05*VX1-U05*VY1*BETA
      AL(3,3)=U05*VX1+2.D0*V05*VY1+W05*VZ1-VY1*BETA*V05
      AL(3,4)=V05*VZ1-BETA*W05*VY1
      AL(3,5)=BETA*VY1
      AL(4,1)=-(U05*W05*VX1+V05*W05*VY1+W05**2*VZ1)+VZ1*BETA/2.D0*UVW
      AL(4,2)=W05*VX1-U05*VZ1*BETA
      AL(4,3)=VY1*W05-VZ1*V05*BETA
      AL(4,4)=W05*VZ1*(3.D0-GM)+U05*VX1+V05*VY1
      AL(4,5)=BETA*VZ1
      AL(5,1)=U1*(BETA*UVW-GM*E05/D05)
      SA5=GM*E05-BETA*D05*UVW/2.D0
      AL(5,2)=-BETA*U05*U1+VX1/D05*SA5
      AL(5,3)=-BETA*V05*U1+VY1/D05*SA5
      AL(5,4)=-BETA*W05*U1+VZ1/D05*SA5
      AL(5,5)=GM*U1
c
        IF (NL.EQ.6) THEN
          AL(1,6)=0.D0
          AL(2,6)=0.D0
          AL(3,6)=0.D0
          AL(4,6)=0.D0
          AL(5,6)=0.D0
          AL(6,1)=-TK_L*U1
          AL(6,2)=TK_L*VX1
          AL(6,3)=TK_L*VY1
          AL(6,4)=TK_L*VZ1
          AL(6,5)=0.d0
          AL(6,6)=U1
        END IF
      END IF

C..   CALCULATE THE RIGHT JACOBIAN MATRIX 

      RA=-1.0
      U05=U_R
      V05=V_R
      W05=W_R
      C05=C_R
      D05=D_R
      E05=E_R

      U1=VX1*U05+VY1*V05+VZ1*W05+QT_R
      QM1=U1/C05
      FMS=RA*D05*C05*(.5*(QM1+RA))**2

      IF(QM1.LE.1.D0) Then           ! subsonic

      ZM=D05**2*(U05**2+V05**2+W05**2)
      AR(1,1)=RA*1.D0/8.D0*GM*BETA*E05/D05/C05*(-QM1**2+1.D0)
      AR(1,2)=RA*1.D0/4.D0*(2.D0*VX1*(QM1+RA)
     $            +.5D0*GM*(GM-1.D0)*U05/C05*(QM1**2 -1.D0))
      AR(1,3)=RA*1.D0/4.D0*(2.D0*VY1*(QM1+RA)+.5*GM*(GM-1.D0)
     $                           *V05/C05*(QM1**2-1.D0))
      AR(1,4)=RA*1.D0/4.D0*(2.D0*VZ1*(QM1+RA)+.5D0*GM*(GM-1.D0)
     $                  *W05/C05*(QM1**2-1.D0))
      AR(1,5)=RA*1.D0/8.D0*GM*BETA/C05*(-QM1**2+1.D0)
      YU=VX1*(-U1+RA*2.D0*C05)/GM+U05
      SA1=(U1/D05+RA*GM*BETA*(-E05/D05**2+ZM/D05**3)/C05)/GM
      SA2=(-VX1/D05-RA*GM*BETA*U05/D05/C05)/GM
      SA3=(-VY1/D05-RA*GM*BETA*V05/D05/C05)/GM
      SA4=(-VZ1/D05-RA*GM*BETA*W05/D05/C05)/GM
      AR(2,1)=AR(1,1)*YU+FMS*(VX1*SA1-U05/D05)
      AR(2,2)=AR(1,2)*YU+FMS*(VX1*SA2+1.D0/D05)
      AR(2,3)=AR(1,3)*YU+FMS*VX1*SA3
      AR(2,4)=AR(1,4)*YU+FMS*VX1*SA4
      AR(2,5)=AR(1,5)*YU+RA*FMS*VX1*BETA/D05/C05
      YM=VY1*(-U1+RA*2.D0*C05)/GM+V05
      AR(3,1)=AR(1,1)*YM+FMS*(VY1*SA1-V05/D05)
      AR(3,2)=AR(1,2)*YM+FMS*VY1*SA2
      AR(3,3)=AR(1,3)*YM+FMS*(VY1*SA3+1.D0/D05)
      AR(3,4)=AR(1,4)*YM+FMS*VY1*SA4
      AR(3,5)=AR(1,5)*YM+FMS*VY1*RA*BETA/D05/C05
      RZ=VZ1*(-U1+RA*2.D0*C05)/GM+W05
      AR(4,1)=AR(1,1)*RZ+FMS*(VZ1*SA1-W05/D05)
      AR(4,2)=AR(1,2)*RZ+FMS*VZ1*SA2
      AR(4,3)=AR(1,3)*RZ+FMS*VZ1*SA3
      AR(4,4)=AR(1,4)*RZ+FMS*(VZ1*SA4+1.D0/D05)
      AR(4,5)=AR(1,5)*RZ+FMS*VZ1*RA*BETA/D05/C05
      AE=(-BETA*U1**2+RA*2.D0*BETA*U1*C05+2.D0*C05**2)/(GM**2-1.D0)
     *   +.5D0*(U05**2+V05**2+W05**2)
      AR(5,1)=AR(1,1)*AE+FMS*((2.D0*BETA*U1**2/D05+RA*GM*BETA**2*
     *  U1*(-3.D0*E05+2.D0*ZM/D05)/D05**2/C05+2.D0*GM*BETA/D05**2*(-E05
     *   +ZM/D05))/(GM**2-1.D0)-ZM/D05**3)
      AR(5,2)=AR(1,2)*AE+FMS*((-2.D0*BETA*U1*VX1/D05+RA*GM*BETA**2*(
     * 2.D0*E05*VX1/D05**2-U05*U1/D05-ZM*VX1/D05**3)/C05-2.D0*GM*BETA*
     *   U05/D05)/(GM**2-1.D0)+U05/D05)
      AR(5,3)=AR(1,3)*AE+FMS*((-2.D0*BETA*U1*VY1/D05+RA*GM*BETA**2*(
     * 2.D0*E05*VY1/D05**2-U1*V05/D05-ZM*VY1/D05**3)/C05-2.D0*GM*BETA*
     *   V05/D05)/(GM**2-1.D0)+V05/D05)
      AR(5,4)=AR(1,4)*AE+FMS*((-2.D0*BETA*U1*VZ1/D05+RA*GM*BETA**2*(
     * 2.D0*E05*VZ1/D05**2-U1*W05/D05-ZM*VZ1/D05**3)/C05-2.D0*GM*BETA*
     *   W05/D05)/(GM**2-1.D0)+W05/D05)
      AR(5,5)=AR(1,5)*AE+FMS*((RA*GM*BETA**2*U1/D05/C05+2.D0*GM*BETA/D05
     *  )/(GM**2-1.D0))
C
        IF (NL.EQ.6) THEN
          AR(1,6)=0.D0
          AR(2,6)=0.D0
          AR(3,6)=0.D0
          AR(4,6)=0.D0
          AR(5,6)=0.D0
          AR(6,1)=0.25D0*RA*TK_R*((0.5D0*GM*BETA*E05/(C05*D05)-C05)*
     $            (1.D0-QM1**2)-2.D0*C05*QM1*(QM1+RA))
          AR(6,2)=0.25D0*RA*TK_R*(2.D0*VX1*(QM1+RA)+
     $            .5D0*GM*BETA*U05/C05*(QM1**2-1.D0))
          AR(6,3)=0.25D0*RA*TK_R*(2.D0*VY1*(QM1+RA)+
     $            .5d0*GM*BETA*V05/C05*(QM1**2-1.D0))
          AR(6,4)=0.25D0*RA*TK_R*(2.D0*VZ1*(QM1+RA)+
     $            .5D0*GM*BETA*W05/C05*(QM1**2-1.D0))
          AR(6,5)=0.125D0*RA*TK_R*GM*BETA/C05*(-QM1**2+1.D0)
          AR(6,6)=0.25D0*RA*C05*(QM1+RA)**2
        END IF
c
      ELSE IF(QM1.GT.1) THEN
        DO  LTY=1,NL
          DO  LTV=1,NL
            AR(LTY,LTV)=0.D0
          END DO
        END DO
      END IF

      RETURN
      END
