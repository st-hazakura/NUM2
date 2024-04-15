!
! Jacobiho metoda pro vypocet v.c. a v.v.
!      
program Jacobi_Calc
    implicit none
   !
   ! deklarace promennych
   !
    integer :: MAXIT=1000
    real :: A(3,3),RVI(3,3),  &
            TOL=1.0E-6
   !
    A=RESHAPE((/ 3.0,-1.0, 1.0  &
               ,-1.0, 5.0,-1.0  &
               , 1.0,-1.0, 3.0/),(/3,3/))
   ! A je   -3 -1  1
   !        -1  5 -1
   !         1 -1  3
    RVI=0.0
   !	
    CALL JACOBI(A,RVI,3,TOL,MAXIT)
   !
    PRINT*," *** VLASTNI CISLA ***"
    PRINT*," lambda_1 =",A(1,1)
    PRINT*," lambda_2 =",A(2,2)
    PRINT*," lambda_3 =",A(3,3)
    print*," "
    PRINT*," *** VLASTNI VEKTORY ***"
    PRINT*," x1 (size):",RVI(1,1),RVI(2,1),RVI(3,1),sqrt(RVI(1,1)**2.0+RVI(2,1)**2.0+RVI(3,1)**2.0)
    PRINT*," x2 (size):",RVI(1,2),RVI(2,2),RVI(3,2),sqrt(RVI(1,2)**2.0+RVI(2,2)**2.0+RVI(3,2)**2.0)
    PRINT*," x3 (size):",RVI(1,3),RVI(2,3),RVI(3,3),sqrt(RVI(1,3)**2.0+RVI(2,3)**2.0+RVI(3,3)**2.0)
    print*," " 
   !
    stop " Jacobi_Calc: Konec vypoctu!"
   end program Jacobi_Calc
   !
   ! DIAGONALIZACE REALNE SYMETRICKE MATICE JACOBIHO METODOU
   !
   SUBROUTINE JACOBI(A,V,N,EPS,MAXIT)
    implicit none
   !
   ! deklarace promennych
   !
    INTEGER :: N,MAXIT,IT,I,J,R,S
    REAL :: EPS,A(N,N),V(N,N),QT(N,N),AMAX,tSUM,C,SGN_Ars,Z,  &
            OMG,RMU,RNU,SGN_mu,TEMP1,TEMP2,TEMP3 
   !
   ! CYCLUS PRES ITERACE
   !
    DO IT=1,MAXIT
   !
   ! VYNULUJ AKUMULATORY
   !
     AMAX=0.0
     tSUM=0.0
   !
   ! NAJDI MAX MIMODIAGONALNI PRVEK
   !
     DO I=1,N-1
      DO J=I+1,N
       tSUM=tSUM+A(I,J)**2
       IF(ABS(A(I,J)) > AMAX) THEN
        AMAX=ABS(A(I,J))
        R=I
        S=J
       ENDIF
      ENDDO
     ENDDO
   !
   ! TEST KONVERGENCE
   !
     tSUM=2.0*tSUM 
     IF(tSUM < EPS) RETURN
   !
   ! SPOCTI c=cos(Phi) A z=sin(Phi)
   !
     IF(A(R,R) == A(S,S)) THEN
      C=SQRT(2.0)/2.0
      SGN_Ars=1.0
      IF(A(R,S) < 0.0) SGN_Ars=-1.0
      Z=SGN_Ars*C
     ELSE
      OMG=(-1.0)*A(R,S)
      RMU=(A(R,R)-A(S,S))/2.0
      RNU=SQRT(OMG**2+RMU**2)
      C=SQRT((RNU+ABS(RMU))/(2.0*RNU))
      SGN_mu=1.0
      IF(RMU < 0.0) SGN_mu=-1.0
      Z=SGN_mu*OMG/(2.0*RNU*C)
      ENDIF
      PRINT*, " Z = ", Z
      PRINT*, " C = ", C 
   !
   ! SPOCTI NOVE (ZMENENE) PRVKY MATICE A
   !
     DO I=1,N
      IF((I == R).OR.(I == S)) cycle
      TEMP1=C*A(I,R)-Z*A(I,S)
      TEMP2=C*A(I,S)+Z*A(I,R)
      A(I,R)=TEMP1
      A(R,I)=TEMP1
      A(I,S)=TEMP2
      A(S,I)=TEMP2
     ENDDO
     TEMP3=C**2*A(R,R)+Z**2*A(S,S)-2.0*C*Z*A(R,S)
     A(S,S)=C**2*A(S,S)+Z**2*A(R,R)+2.0*C*Z*A(R,S)
     A(R,R)=TEMP3
     A(R,S)=0.0
     A(S,R)=0.0
   !
   !  SPOCTI V.V
   !
     QT=0.0
     I=0
     DO
      I=I+1
      IF(I > N) exit
      QT(I,I)=1.0
     ENDDO
     QT(R,R)=C
     QT(S,S)=C
     QT(R,S)=Z
     QT(S,R)=(-1.0)*Z
     IF(IT == 1) THEN
      V=QT
     ELSE
      V=MATMUL(V,QT)
     ENDIF
    ENDDO
   END SUBROUTINE JACOBI