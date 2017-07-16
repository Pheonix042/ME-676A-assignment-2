       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C       DEFINITIONS
C       -----------------------------------------
C       ROMIL KADIA(16105045)
C       ANKUR MAURYA(13124)
C       ROHIT KUMAVAT(13587)
C       -----------------------------------------
C     DEFINING DIMENSION OF VARIABLES
      DIMENSION STT(6),DS(6),DSIG(6),DQQ(6),DLAM(6),SIGV(6)
C     STT=TRIAN STRAIN, DS=DEVIATORIC TRIAL STRAIN
C     DSIG=DEVIATORIC OF STRESS TRIAL, DQQ=DERIVATIVE OF QQ
C     DLAM=DERIVATIVE OF DELTALAMBDA, SIGV= VOL STRESS
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     1 SIX=6.0D0)
C     -----------------------------------------
C     DIVEN DATA
C     -----------------------------------------
C     DEFINING CONSTANTS
      SP=STATEV(1)
      H=10.0D0
      ZZZ=1.0D0
      E=1000.0D0
      MU=0.3D0
      GG=E/(TWO*(ONE+MU))
      KK=E/(THREE*(ONE-TWO*MU))
      LAMBDA=E*MU/((ONE+MU)*(ONE-TWO*MU))
C     -----------------------------------------
C     DATA MANIPULATION
      A = ZERO
      DO I=1,6
            STT(I)=STRAN(I) + DSTRAN(I)
      END DO
      TEMP=(STT(1)+STT(2)+STT(3))/THREE
      DO I=1,3
            DS(I)=STT(I)-TEMP
      END DO
      DS(4)=STT(4)
      DS(5)=STT(5)
      DS(6)=STT(6)
      DO I=1,6
            DSIG(I)=TWO*GG*DS(I)
      END DO
      DO I=1,6
            A=A+((STT(I))**2)
      END DO
      QQ=TWO*(SQRT(THREE))*GG*(SQRT(A))
      LAM=(QQ-(H*SP)*ONE)/(H+THREE*GG)
      QCON=(TWO*SQRT(THREE)*GG)/SQRT(A)
      DO I=1,6
            DQQ(I)=QCON*DS(I)
      END DO
      LAMCON=(THREE*GG)/((THREE*GG+H)*SQRT(A))
      DO I=1,6
            DLAM(I)=LAMCON*DS(I)
      END DO
      GCON=(GG**2)/(H + (THREE*GG))
      VOL=(STT(1)+STT(2)+STT(3))/3
      DO I=1,3
            SIGV(I)=KK*VOL
      END DO
      DO I=4,6
            SIGV(I)=ZERO
      END DO
C     -----------------------------------------
C     DDSDDE CALCULATION
      BB= - (TWO*(GG**2)*LAM)/QQ
      CC=(TWO*SIX*(GG**2)*LAM)/QQ
      X= - (((THREE*GCON*QCON)/QQ)-((SIX*(GG**2)*LAM*QCON)/(QQ**2)))
      Y=((VOL*GCON*QCON)/QQ) + ((TWO*(GG**2)*LAM*QCON*VOL)/(QQ**2))
C     -----------------------------------------
C     PHI CALCULATION
      JJ=ZERO
      DO I=1,6
            JJ=JJ + DSIG(I)
      END DO
      W=THREE*JJ
      PHI=SQRT(W) - ZZZ - H*SP
      QW=(LAMBDA+TWO*MU)+(TWO*BB)+CC
      TTTT=QW+((X*STT(1))+Y)*DS(1)
      DO I=1,6
            DO J=1,6
                  DDSDDE(I,J)=ZERO
            END DO
      END DO
C      IF(PHI.LE.ZERO) THEN
C            DO I=1,6
C                  STRESS(I)=DSIG(I)+SIGV(I)
C            END DO
C            STATEV(1)=SP
C      ELSE
C            DDSDDE(1,1)=QW+((X*STT(1))+Y)*DS(1)
C            DDSDDE(2,2)=QW+((X*STT(2))+Y)*DS(2)
C            DDSDDE(3,3)=QW+((X*STT(3))+Y)*DS(3)
C            DO I=1,3
C                  DDSDDE(I,I)=QW+((X*STT(I))+Y)*DS(I)
C            END DO
C            DO I=4,6
C                  DDSDDE(I,I)=MU+BB+X*STT(I)*DS(I)
C            END DO
C            DDSDDE(1,2)=LAMBDA+CC+((X*STT(1))+Y)*DS(2)
C            DDSDDE(1,3)=LAMBDA+CC+((X*STT(1))+Y)*DS(3)
C            DDSDDE(2,3)=LAMBDA+CC+((X*STT(2))+Y)*DS(3)
C            DO I=1,3
C                  DO J=4,6
C                        DDSDDE(I,J)=((X*STT(I))+Y)*DS(J)
C                  END DO
C            END DO
C            DDSDDE(4,5)=X*STT(4)*DS(5)
C            DDSDDE(4,6)=X*STT(4)*DS(6)
C            DDSDDE(5,6)=X*STT(5)*DS(6)
C            DO I=1,6
C                  DO J=1,6
C                        DDSDDE(J,I)=DDSDDE(I,J)
C                  END DO
C            END DO
C            SP=SP+LAM
C            STATEV(1)=SP
            DO I=1,6
                DO J=1,6
                      STRESS(I)=STRESS(I)+DDSDDE(I,J)*DS(J)
                END DO
            END DO
C      END IF
      RETURN
      END
