       SUBROUTINE ivsnr(A,N,EP)

c..     THIS IS A SUBROUTINE TO CALCULATE THE
C..     INVERSE OF A MATRIX BY PIVOT GAUSS-JODON ELIMENATION
C..     A: INPUT:  MATRIX OF SIZE NxN,
c..        OUTPUT: INVERSE  MATRIX OF A
C..     B,CD: WORKING ARRAYS
C..     N: INPUT, THE ORDER OF THE MATRIX
C..     ME: WORKING ARRAY FOR INTEGERS
C..     DE: OUTPUT, THE VALUE OF THE DETERMINANT
C..     EP: INPUT, =1E-06, IF ABS(PIOVT)<EP,
C..     THE MATRIX IS CONSIDERED TO HAVE SINGULARITY

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION A(N,N),ME(N),B(N),CD(N)
      DE=1.d0
      DO 10 J=1,N
10    ME(J)=J
      DO 20 I=1,N
      Y=0.d0
      DO 30 J=I,N
      IF(DABS(A(I,J)).LE.DABS(Y)) GO TO 30
      K=J
      Y=A(I,J)
30    CONTINUE
      DE=DE*Y
      IF(DABS(Y).LT.EP) GO TO 32
      Y=1.D0/Y
      DO 40 J=1,N
      CD(J)=A(J,K)
      A(J,K)=A(J,I)
      A(J,I)=-CD(J)*Y
      B(J)=A(I,J)*Y
40    A(I,J)=A(I,J)*Y
      A(I,I)=Y
      J=ME(I)
      ME(I)=ME(K)
      ME(K)=J
      DO 11 K=1,N
      IF(K.EQ.I) GO TO 11
      DO 12 J=1,N
      IF(J.EQ.I) GO TO 12
      A(K,J)=A(K,J)-B(J)*CD(K)
12    CONTINUE
11    CONTINUE
20    CONTINUE
      DO 33 I=1,N
      DO 44 K=1,N
      IF(ME(K).EQ.I) GO TO 55
44    CONTINUE
55    IF(K.EQ.I) GO TO 33
      DO 66 J=1,N
      W=A(I,J)
      A(I,J)=A(K,J)
66    A(K,J)=W
      IW=ME(I)
      ME(I)=ME(K)
      ME(K)=IW
      DE=-DE
33    CONTINUE
      RETURN
32    DE=0.D0
      RETURN
      END

