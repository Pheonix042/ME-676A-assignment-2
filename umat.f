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
      E=1000.0D0
      MU=0.33D0
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
      VOL=STT(1)+STT(2)+STT(3)
      DO I=1,3
            SIGV(I)=KK*VOL/THREE
      END DO
      DO I=4,6
            SIGV(I)=ZERO
      END DO
C     -----------------------------------------
C     DDSDDE CALCULATION
      BB= - (TWO*(GG**2)*LAM)/QQ
      CC=(12*(GG**2)*LAM)/QQ
      X= - ((THREE*GCON*QCON)/QQ)-((SIX*(GG**2)*LAM*QCON)/(QQ**2))
      Y=((VOL*GCON*QCON)/QQ) + ((TWO*(GG**2)*LAM*QCON*VOL)/(QQ**2))
C     -----------------------------------------
C     PHI CALCULATION
      JJ=ZERO
      DO I=1,6
            JJ=JJ + DSIG(I)
      END DO
      W=THREE*JJ
      PHI=SQRT(W) - ONE - H*SP
      IF(PHI.LE.ZERO) THEN
            DO I=1,6
                  STRESS(I)=DSIG(I)+SIGV(I)
            END DO
            STATEV(1)=SP
      ELSE
            DO I=1,3
                  DDSDDE(I,I)=(LAMBDA+TWO*MU)+(TWO*BB)+CC+((X*STT(I))+Y)*DS(I)
            END DO
            DO I=4,6
                  DDSDDE(I,I)=MU+BB+X*STT(I)*DS(I)
            END DO
            DDSDDE(1,2)=LAMBDA+CC+((X*STT(1))+Y)*DS(2)
            DDSDDE(1,3)=LAMBDA+CC+((X*STT(1))+Y)*DS(3)
            DDSDDE(2,3)=LAMBDA+CC+((X*STT(2))+Y)*DS(3)
            DO I=1,3
                  DO J=4,6
                        DDSDDE(I,J)=((X*STT(I))+Y)*DS(J)
                  END DO
            END DO
            DDSDDE(4,5)=X*STT(4)*DS(5)
            DDSDDE(4,6)=X*STT(4)*DS(6)
            DDSDDE(5,6)=X*STT(5)*DS(6)
            DO I=1,6
                  DO J=1,6
                        DDSDDE(J,I)=DDSDDE(I,J)
                  END DO
            END DO
            SP=SP+LAM
            STATEV(1)=SP
            DO I=1,6
                  DO J=1,6
                        STRESS(I)=STRESS(I)+DDSDDE(I,J)*DS(J)
                  END DO
            END DO
      END IF
      RETURN
      END
