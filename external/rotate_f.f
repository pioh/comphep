      real*8 function 
     & rotate5(v11,v12,v13,v14,v15,v22,v23,v24,v25,v33,
     &          v34,v35,v44,v45,v55,ri,rj)

	 implicit none

	 REAL*8 AR(5,5),WR(5),ZR(5,5),WORK(10000),RI,RJ
	 
	 REAL*8 v11,v12,v13,v14,v15,v22,v23,v24,v25,v33,
     &          v34,v35,v44,v45,v55
	 integer I,J,IERR,ND,K
	 REAL*8 TEMP

	 ND=5

	 AR(1,1)=v11
	 AR(1,2)=v12
	 AR(1,3)=v13
	 AR(1,4)=v14
	 AR(1,5)=v14
	 
	 AR(2,1)=v12
	 AR(2,2)=v22
	 AR(2,3)=v23
	 AR(2,4)=v24	 
	 AR(2,5)=v25
	 	 
	 AR(3,1)=v13
	 AR(3,2)=v23
	 AR(3,3)=v33
	 AR(3,4)=v34
	 AR(3,5)=v35
	 
	 AR(4,1)=v14
	 AR(4,2)=v24
	 AR(4,3)=v34
	 AR(4,4)=v44
	 AR(4,5)=v45
	 	
	 AR(5,1)=v15
	 AR(5,2)=v25
	 AR(5,3)=v35
	 AR(5,4)=v45
	 AR(5,5)=v55
	 
	 
         CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)
 100	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	    GOTO 100
	   END IF  
	  ENDDO
	 
	 I=RI
	 J=RJ
	 rotate5=ZR(I,J)
	 return
	 END
ccccccccccccccc
	 real*8 function 
     & pmass5(v11,v12,v13,v14,v15,v22,v23,v24,v25,v33,
     &          v34,v35,v44,v45,v55,ri)
	
	 implicit none
         
	 REAL*8 AR(5,5),WR(5),ZR(5,5),WORK(10000),RI,RJ
	 
	 REAL*8 v11,v12,v13,v14,v15,v22,v23,v24,v25,v33,
     &          v34,v35,v44,v45,v55
	 integer I,J,IERR,ND,K
	 
	 REAL*8 TEMP

	 ND=5

	 AR(1,1)=v11
	 AR(1,2)=v12
	 AR(1,3)=v13
	 AR(1,4)=v14
	 AR(1,5)=v14
	 
	 AR(2,1)=v12
	 AR(2,2)=v22
	 AR(2,3)=v23
	 AR(2,4)=v24
	 AR(2,5)=v25

	 AR(3,1)=v13
	 AR(3,2)=v23
	 AR(3,3)=v33
	 AR(3,4)=v34
	 AR(3,5)=v35

	 AR(4,1)=v14
	 AR(4,2)=v24
	 AR(4,3)=v34
	 AR(4,4)=v44
	 AR(4,5)=v45

	 AR(5,1)=v15
	 AR(5,2)=v25
	 AR(5,3)=v35
	 AR(5,4)=v45
	 AR(5,5)=v55

         CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)
 100	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	    GOTO 100
	   END IF  
	  ENDDO
	 
	 I=RI
	 J=RJ
	 pmass5=sqrt(WR(I))
	 
	 
	 return
	 END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
	real*8 function 
     &    rotate4(v11,v12,v13,v14,v22,v23,v24,v33,v34,v44,ri,rj)
       
	implicit none
	 
	REAL*8 AR(4,4),WR(4),ZR(4,4),WORK(10000),RI,RJ
	
	REAL*8 v11,v12,v13,v14,v22,v23,v24,v33,v34,v44
	integer I,J,IERR,IC,ND,K
        REAL*8 TEMP
	
	ND=4
	
	
	AR(1,1)=v11
	AR(1,2)=v12
	AR(1,3)=v13
	AR(1,4)=v14
	
	AR(2,1)=v12
	AR(2,2)=v22
	AR(2,3)=v23
	AR(2,4)=v24
		
	AR(3,1)=v13
	AR(3,2)=v23
	AR(3,3)=v33
	AR(3,4)=v34
	
	AR(4,1)=v14
	AR(4,2)=v24
	AR(4,3)=v34
	AR(4,4)=v44
	
	 CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)
 
 100    IC=0
	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    IC=1
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	   END IF  
	  ENDDO
	  IF(IC.ne.0) GOTO 100

	I=RI
	J=RJ
	rotate4=ZR(I,J)

	return
	END

cccccccccccccc

	 real*8 function 
     &     pmass4(v11,v12,v13,v14,v22,v23,v24,v33,
     &          v34,v44,ri)
	
	 implicit none
         
	 REAL*8 AR(4,4),WR(4),ZR(4,4),WORK(10000),RI,RJ
	 
	 REAL*8 v11,v12,v13,v14,v15,v22,v23,v24,v25,v33,
     &          v34,v35,v44,v45,v55
	 integer I,J,IERR,ND,K
	 REAL TEMP
	
	 ND=4

	 AR(1,1)=v11
	 AR(1,2)=v12
	 AR(1,3)=v13
	 AR(1,4)=v14

	 AR(2,1)=v12
	 AR(2,2)=v22
	 AR(2,3)=v23
	 AR(2,4)=v24

	 AR(3,1)=v13
	 AR(3,2)=v23
	 AR(3,3)=v33
	 AR(3,4)=v34
	 
	 AR(4,1)=v14
	 AR(4,2)=v24
	 AR(4,3)=v34
	 AR(4,4)=v44
	 
         CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)
 
 100	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	    GOTO 100
	   END IF
	  ENDDO

	 I=RI
	 J=RJ
	 pmass4=sqrt(abs(WR(I)))

	 return
	 END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 real*8 function 
     & rotate3(v11,v12,v13,v22,v23,v33,ri,rj)

	 implicit none

	 REAL*8 AR(3,3),WR(3),ZR(3,3),WORK(10000),RI,RJ

	 REAL*8 v11,v12,v13,v22,v23,v33
	 REAL*8 TEMP,pmass3
	 integer I,J,IERR,K,ND

	 ND=3

	 AR(1,1)=v11
	 AR(1,2)=v12
	 AR(1,3)=v13
	 
	 AR(2,1)=v12
	 AR(2,2)=v22
	 AR(2,3)=v23
	 
	 AR(3,1)=v13
	 AR(3,2)=v23
	 AR(3,3)=v33

         CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)

 100	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	    GOTO 100
	   END IF
	  ENDDO

	 I=RI
	 J=RJ
	 rotate3=ZR(I,J)

	 return
	 END

cccccccccccccc

	 real*8 function 
     &     pmass3(v11,v12,v13,v22,v23,v33,ri)

	 implicit none

	 REAL*8 AR(3,3),WR(3),ZR(3,3),WORK(10000),RI,RJ

	 REAL*8 v11,v12,v13,v14,v15,v22,v23,v24,v25,v33,
     &          v34,v35,v44,v45,v55
	 integer I,J,IERR,ND,K
	 REAL*8 TEMP

	 ND=3

	 AR(1,1)=v11
	 AR(1,2)=v12
	 AR(1,3)=v13

	 AR(2,1)=v12
	 AR(2,2)=v22
	 AR(2,3)=v23

	 AR(3,1)=v13
	 AR(3,2)=v23
	 AR(3,3)=v33
	 
         CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)
 
 100	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	    GOTO 100
	   END IF
	  ENDDO

	 I=RI
	 J=RJ
	 pmass3=sqrt(abs(WR(I)))

	 return
	 END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 real*8 function 
     & rotate2(v11,v12,v22,ri,rj)

	 implicit none
         
	 REAL*8 AR(2,2),WR(2),ZR(2,2),WORK(10000),RI,RJ
	 
	 REAL*8 v11,v12,v22
	 integer I,J,IERR,ND,K
	 REAL*8 TEMP

	 ND=2

	 AR(1,1)=v11
	 AR(1,2)=v12
	 AR(2,1)=v12
	 AR(2,2)=v22

         CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)
 100	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	    GOTO 100
	   END IF
	  ENDDO

	 I=RI
	 J=RJ
	 rotate2=ZR(I,J)

	 return
	 END

	 real*8 function 
     & pmass2(v11,v12,v22,ri)

	 implicit none

	 REAL*8 AR(2,2),WR(2),ZR(2,2),WORK(10000),RI,RJ
	 
	 REAL*8 v11,v12,v22
	 integer I,J,IERR,K,ND
	 REAL*8 TEMP

	 ND=2

	 AR(1,1)=v11
	 AR(1,2)=v12
	 AR(2,1)=v12
	 AR(2,2)=v22
	
	 
         CALL EISRS1(ND,ND,AR,WR,ZR,IERR,WORK)
 100	DO I=1,ND-1
 	  IF (ABS(WR(I)).GT.ABS(WR(I+1))) THEN
	    TEMP=WR(I+1)
	    WR(I+1)=WR(I)
	    WR(I)=TEMP
	    DO K=1,ND
	      TEMP=ZR(K,I+1)
	      ZR(K,I+1)=ZR(K,I)
	      ZR(K,I)=TEMP
	    ENDDO
	    GOTO 100
	   END IF  
	  ENDDO

	 I=RI
	 J=RJ
	 pmass2=sqrt(abs(WR(I)))

	 return
	 END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE EISRS1(NM,N,AR,WR,ZR,IERR,WORK)
      IMPLICIT NONE
      INTEGER NM,N,IERR
      REAL*8 AR(NM,NM),WR(N),ZR(NM,NM),WORK(10000)
      CALL TRED2(NM,N,AR,WR,WORK,ZR)
      CALL TQL2(NM,N,WR,WORK,ZR,IERR)
CAB      
      CALL EIGEN_ORDER(N,WR,ZR)
      RETURN
      END

      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C     FROM CERN PROGRAM LIBRARY

      IMPLICIT NONE
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
      REAL*8 F,G,H,HH,SCALE
      
      DO 100 I = 1, N
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
      IF (N .EQ. 1) GO TO 320
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0
         SCALE = 0.0
         IF (L .LT. 2) GO TO 130
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
         IF (SCALE .NE. 0.0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / (SCALE * H)
            G = 0.0
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
         HH = F / (H + H)
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
         DO 280 K = 1, L
  280    Z(I,K) = SCALE * Z(I,K)
  290    D(I) = H
  300 CONTINUE
  320 D(1) = 0.0
      E(1) = 0.0
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0) GO TO 380
         DO 360 J = 1, L
            G = 0.0
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0
         IF (L .LT. 1) GO TO 500
         DO 400 J = 1, L
            Z(I,J) = 0.0
            Z(J,I) = 0.0
  400    CONTINUE
  500 CONTINUE
      RETURN
      END

cccccccccccccccccccccccccccccccccccccc

      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C     FROM CERN PROGRAM LIBRARY
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      MACHEP=2.**(-23)
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      F = 0.0
      B = 0.0
      E(N) = 0.0
      DO 240 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H
         DO 110 M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
  110    CONTINUE
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
         P = (D(L+1) - D(L)) / (2.0 * E(L))
         R = SQRT(P*P+1.0)
         H = D(L) - E(L) / (P + SIGN(R,P))
         DO 140 I = L, N
  140    D(I) = D(I) - H
         F = F + H
         P = D(M)
         C = 1.0
         S = 0.0
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0)
            E(I+1) = S * E(I) * R
            S = 1.0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
  200    CONTINUE
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
      END

cccccccccccccccccccccccccccccccccccccc

      SUBROUTINE EIGEN_ORDER(N,WR,ZR)
C     SASHA
      IMPLICIT NONE
      INTEGER N
      INTEGER I,J,JMAX(N)
      REAL*8 WR(N),ZR(N,N),WRTMP(N),ZRTMP(N,N),MAXIJ
      
      DO J=1,N
        MAXIJ=ZR(J,J)**2
	JMAX(J)=J
        DO I=1,N
	  IF(MAX(ZR(I,J)**2,MAXIJ).gt.MAXIJ) THEN
	    JMAX(J)=I
	  ENDIF
	ENDDO
  	DO I=1,N
	  ZRTMP(I,JMAX(J)) = ZR(I,J)
	ENDDO
	WRTMP(JMAX(J)) = WR(J)
      ENDDO
      
      DO I=1,N
        DO J=1,N
	  ZR(I,J)=ZRTMP(I,J)
	ENDDO
	WR(I) = WRTMP(I)
      ENDDO
      RETURN
      END
