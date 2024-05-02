      subroutine pbtrip(a, b, c, d, il, iu, order)
c
c     subroutine to solve periodic block tridiagonal
c     system of equations without pivoting strategy
c     each block matrix may be of dimension n with n
c     any number greater than 1
c      
c     solves block tridiagonal system 
c
c                        A x = D
c
c     where
c
c            |  B1  C1          A1  |
c            |  A2  B2  C2          |
c     A  =   |      A3  B3  C3      |
c            |          ...         |
c            |  CN          AN  BN  |
c
c
      implicit none
c
c     Interface Variables

      integer, intent(in)::
     $     order, il, iu

      double precision, intent(inout)::
     $     a(1), b(1), c(1), d(1)
c
c     Local Variables

      integer::
     $     ordsq,
     $     is, ie, iuvec, iumat, iemat, ievec,
     $     i, i0mat, i0matj, i1mat, i0vec, i1vec,
     $     iumatj, i1matj, ibac, iematj,
     $     j

      double precision ad(order*order), cd(order*order)
c
c     a(n,n,ni) = sub diagonal matrix
c     b(n,n,ni) =     diagonal matrix
c     c(n,n,ni) = sup diagonal matrix
c     d(n,ni) = right hand side vector
c     il = lower value of index for which matrices are defined
c     iu = upper value of index for which matrices are defined
c          (solution is sought for btri(a,b,c)*x=d
c          for indices of x between il and iu(inclusive).
c          solution written in d vector (orginal contents are
c          overwritten))
c     iu-il+1 = ni
c     order = n, order of a, b, c matrices and length of d vector
c             at each point renoted by index i (order can be
c             any integer greater than 1).
c             (arrays ad and cd must be at least of length
c              order**2, current length of 25 anticipates
c              maximum order of 5)

      is = il+1
      ie = iu-1
      ordsq = order**2
      iumat = 1+(iu-1)*ordsq
      iuvec = 1+(iu-1)*order
      iemat = 1+(ie-1)*ordsq
      ievec = 1+(ie-1)*order
c
c     forward elimination
c
      i = il
      i0mat = 1+(i-1)*ordsq
      i0vec = 1+(i-1)*order
      call ludeco(b(i0mat),order)
      call lusolv(b(i0mat),d(i0vec),d(i0vec),order)
      do 10 j = 1,order
         i0matj = i0mat+(j-1)*order
         call lusolv(b(i0mat),c(i0matj),c(i0matj),order)
         call lusolv(b(i0mat),a(i0matj),a(i0matj),order)
 10   continue
c     
      do 200 i = is,ie
         i0mat = 1+(i-1)*ordsq
         i0vec = 1+(i-1)*order
         i1mat = i0mat-ordsq
         i1vec = i0vec-order
         do 20 j = 1,ordsq
            i0matj = j-1+i0mat
            iumatj = j-1+iumat
            ad(j) = a(i0matj)
            cd(j) = c(iumatj)
            a(i0matj) = 0.d0
            c(iumatj) = 0.d0
 20      continue
         call mulput(ad,d(i1vec),d(i0vec),order)
         do 22 j = 1,order
            i0matj = i0mat+(j-1)*order
            i1matj = i1mat+(j-1)*order
            call mulput(ad,c(i1matj),b(i0matj),order)
            call mulput(ad,a(i1matj),a(i0matj),order)
 22      continue
         call ludeco(b(i0mat),order)
         call lusolv(b(i0mat),d(i0vec),d(i0vec),order)
         do 24 j = 1,order
            i0matj = i0mat+(j-1)*order
            call lusolv(b(i0mat),c(i0matj),c(i0matj),order)
            call lusolv(b(i0mat),a(i0matj),a(i0matj),order)
 24      continue
         call mulput(cd,d(i1vec),d(iuvec),order)
         do 26 j = 1,order
            iumatj = iumat+(j-1)*order
            i1matj = i1mat+(j-1)*order
            call mulput(cd,a(i1matj),b(iumatj),order)
            call mulput(cd,c(i1matj),c(iumatj),order)
 26      continue
 200  continue
c
      do 30 j = 1,ordsq
         iumatj = j-1+iumat
         ad(j) = a(iumatj)+c(iumatj)
 30   continue
      call mulput(ad,d(ievec),d(iuvec),order)
      do 32 j = 1,order
         iumatj = iumat+(j-1)*order
         iematj = iemat+(j-1)*order
         call mulput(ad,c(iematj),b(iumatj),order)
         call mulput(ad,a(iematj),b(iumatj),order)
 32   continue
      call ludeco(b(iumat),order)
      call lusolv(b(iumat),d(iuvec),d(iuvec),order)
c
c     back substitution
c
      do 40 ibac = il,ie
         i = ie-ibac+il
         i0mat = 1+(i-1)*ordsq
         i0vec = 1+(i-1)*order
         i1vec = i0vec+order
         call mulput(a(i0mat),d(iuvec),d(i0vec),order)
         call mulput(c(i0mat),d(i1vec),d(i0vec),order)
 40   continue
c
      return
      end

      subroutine ludeco(a, order)
c
c     subroutine to calculate L-U decomposition
c     of a given matrix a and store result in a
c     (no pivoting strategy is employed)
c
      implicit none
c
c     Interface Variables
      
      integer, intent(in)::
     $     order

      double precision, intent(inout)::
     $     a(order,1)
c
c     Local Variables

      integer::
     $     jrjc, jrjcm1, jrjcp1, jr, jm, jc
      
      double precision::
     $     sum
c
c     Begin
c
      do 8 jc = 2,order
         a(1,jc)=a(1,jc)/a(1,1)
 8    continue
      jrjc = 1
 10   continue
      jrjc = jrjc+1
      jrjcm1 = jrjc-1
      jrjcp1 = jrjc+1
      do 14 jr=jrjc,order
         sum = a(jr,jrjc)
         do 12 jm = 1,jrjcm1
            sum = sum-a(jr,jm)*a(jm,jrjc)
 12      continue
         a(jr,jrjc) = sum
 14   continue
      if(jrjc.eq.order) return
      do 18 jc = jrjcp1,order
         sum = a(jrjc,jc)
         do 16 jm = 1, jrjcm1
            sum = sum-a(jrjc,jm)*a(jm,jc)
 16      continue
         a(jrjc,jc) = sum/a(jrjc,jrjc)
 18   continue
      go to 10
      end

      subroutine mulput(a, b, c, order)
c
c     subroutine to multply a vector b by a matrix a
c     subtract result froma another vector c and store
c     result in c. thus vector c is overwritten
c
      implicit none
c
c     Interface Variables
c      
      integer, intent(in)::
     $     order

      double precision, intent(in)::
     $     a(1), b(1)

      double precision, intent(inout)::
     $     c(1)
c
c     Local Variables
c
      integer::
     $     jr, jc, ia

      double precision::
     $     sum
c      
      do 200 jr = 1,order
         sum = 0.d0
         do 100 jc = 1,order
            ia = jr+(jc-1)*order
            sum = sum+a(ia)*b(jc)
 100     continue
         c(jr) = c(jr)-sum
 200  continue
c
      return
      end

      subroutine lusolv(a, b, c, order)
c
c     subroutine to solve linear algebric system of
c     equations a*c=b and store results in vector c.
c     matrix a is input in L-U decomposition form.
c     (no pivoting strategy has been employed to
c     compute the L-U decomposition of the matrix a).
c     
      implicit none
c
c     Interface Variables
      
      integer, intent(in)::
     $     order

      double precision, intent(in)::
     $     a(order,1), b(1)

      double precision, intent(inout)::
     $     c(1)
c
c     Local Variables

      integer::
     $     jr, jrm1, jm, jrjr, jrp1, jmjm
      
      double precision::
     $     sum
c
c     first l(inv)*b
c
      c(1) = c(1)/a(1,1)
      do 14 jr = 2,order
         jrm1 = jr-1
         sum = b(jr)
         do 12 jm = 1,jrm1
            sum = sum-a(jr,jm)*c(jm)
 12      continue
         c(jr) = sum/a(jr,jr)
 14   continue
c
c     next u(inv) of l(inv)*b
c         
      do 18 jrjr=2,order
         jr = order-jrjr+1
         jrp1 = jr+1
         sum = c(jr)
         do 16 jmjm = jrp1,order
            jm = order-jmjm+jrp1
            sum = sum-a(jr,jm)*c(jm)
 16      continue
         c(jr)=sum
 18   continue
c
      return
      end
