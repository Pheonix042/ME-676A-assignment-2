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
C-----------CODE FROM HERE-----------
C-------We need to solve for delta Lambda-----------
C-----First we need to get trail Strain for n+1 step------------
C------Define Some Variables---------------
C
      DIMENSION STNMT(3,3),DSTNMT(3,3),STNTR(3,3), STNVTR(3,3),
     1 STNDTR(3,3), STRVTR(3,3), STRDTR(3,3),
     2 STRDTV(6),STRTRI(6), SNDTRI(6)
C
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     1 SIX=6.0D0)
C
C----------WRITE GIVEN DATA------
C
      E=1000.0D0
      NU=0.3D0
      SIG0=1.0D0
      H=10.0D0
C
C-----GET G AND K-----------
C
      G=E/(TWO*(ONE+NU))
      K=E/(THREE*(1-TWO*NU))
      LAMB=E*NU/((ONE+NU)*(ONE-TWO*NU))
C
C---------GET ACCUMULATED STRAIN---------
C
      EPP=STATEV(1)
C----------GET STRAIN AND DSTRAIN MATRIX FROM THE DATA-------------
C
      DO I=1,3
            STNMT(I,I)=STRAN(I)
            DSTNMT(I,I)=DSTRAN(I)
      END DO
      STNMT(1,2)=STRAN(4)
      STNMT(1,3)=STRAN(5)
      STNMT(2,3)=STRAN(6)
      DSTNMT(1,2)=DSTRAN(4)
      DSTNMT(1,3)=DSTRAN(5)
      DSTNMT(2,3)=DSTRAN(6)
      DO I=1,3
            DO J=1,I-1
                  STNMT(I,J)=STNMT(J,I)
                  DSTNMT(I,J)=DSTNMT(J,I)
            END DO
      END DO
C
C--------CALCULATE TRIAL STRAIN------
C
      DO I=1,3
            DO J=1,3
                  STNTR(I,J)=STNMT(I,J)+DSTNMT(I,J)
            END DO
      END DO
C
C----------CALCULATE VOLUMETRIC AND DEVIATORIC PART OF TRIAL STRAIN------
C
      STNPSS=(STNTR(1,1)+STNTR(2,2)+STNTR(3,3))/THREE
      DO I=1,3
            DO J=1,3
                  IF(I.EQ.J) THEN
                        STNVTR(I,J)=STNPSS
                  ELSE
                        STNVTR(I,J)=ZERO
                  END IF
                  STNDTR(I,J)=STNTR(I,J)-STNVTR(I,J)
            END DO
      END DO
C
C-------------GET TRIAL STRESSES-------
C
      STRVTR=K*STNVTR
      STRDTR=TWO*G*STNDTR
C
C-----NOW WE CAN FIND Q TRIAL AND ||SIGMA'||^2 AS NEEDED LATER------
C
      MSSQ=ZERO
      DO I=1,3
            DO J=1,3
                  MSSQ=MSSQ+STRDTR(I,J)*STRDTR(I,J)
            END DO
      END DO
      QTR=(THREE*MSSQ)/TWO
      QTR=SQRT(QTR)
C
C--------CHECK PHI--------
C
      PHI=QTR-SIG0-H*EPP
      IF(PHI.LE.ZERO) THEN
      STATEV(1)=EPP
      DO I=1,3
            STRESS(I)=STRVTR(I,I)+STRDTR(I,I)
      END DO
      STRESS(4)=STRVTR(1,2)+STRDTR(1,2)
      STRESS(5)=STRVTR(1,3)+STRDTR(1,3)
      STRESS(6)=STRVTR(2,3)+STRDTR(2,3)
C
C----------OR-----------
C
      ELSE
C
C------NOW WE CAN MOVE ON TO DEL LAMBDA----
C
      DLAM=(QTR-SIG0-H*EPP)/(THREE*G+H)
C
C--------UPDATE EPP AND STATEV--------
      EPP=EPP+DLAM
      STATEV(1)=EPP
C
C------WRITE DEV STRESS IN VOIGT NOTA.STRDTV--
C
      DO I=1,3
            STRDTV(I)=STRDTR(I,I)
      END DO
            STRDTV(4)=STRDTR(1,2)
            STRDTV(5)=STRDTR(1,3)
            STRDTV(6)=STRDTR(2,3)
C
C---DEFINE TERMS AS NEEDED IN DDSDDE---
C
      Q=(SIX*DLAM*G*G)/QTR
      R=((SIX*G*G)/MSSQ)*((ONE/(3*G+H))-(DLAM/QTR))
C
C---------DDSDDE---------
C
      DO I=1,3
            DDSDDE(I,I)=LAMB+TWO*G-(TWO*Q/THREE)-R*STRDTV(I)**2
            J=I+3
            DDSDDE(J,J)=G-(Q/TWO)-R*STRDTV(J)**2
      END DO
      DDSDDE(1,2)=LAMB+(Q/THREE)-R*STRDTV(1)*STRDTV(2)
      DDSDDE(1,3)=LAMB+(Q/THREE)-R*STRDTV(1)*STRDTV(3)
      DDSDDE(1,4)=-R*STRDTV(1)*STRDTV(4)
      DDSDDE(1,5)=-R*STRDTV(1)*STRDTV(5)
      DDSDDE(1,6)=-R*STRDTV(1)*STRDTV(6)
      DDSDDE(2,3)=LAMB+(Q/THREE)-R*STRDTV(2)*STRDTV(3)
      DDSDDE(2,4)=-R*STRDTV(2)*STRDTV(4)
      DDSDDE(2,5)=-R*STRDTV(2)*STRDTV(5)
      DDSDDE(2,6)=-R*STRDTV(2)*STRDTV(6)
      DDSDDE(3,4)=-R*STRDTV(3)*STRDTV(4)
      DDSDDE(3,5)=-R*STRDTV(3)*STRDTV(5)
      DDSDDE(4,5)=-R*STRDTV(4)*STRDTV(5)
      DDSDDE(4,6)=-R*STRDTV(4)*STRDTV(6)
      DDSDDE(5,6)=-R*STRDTV(5)*STRDTV(6)
C
C----------FINALLY STRESS UPDATE-----
C
C-----------LET'S WRITE TRIAL STRESS AND DEVI STRAIN TRIAL IN VOIGT---
C
      DO I=1,3
            STRTRI(I)=STRVTR(I,I)+STRDTR(I,I)
            SNDTRI(I)=STNDTR(I,I)
      END DO
      STRTRI(4)=STRVTR(1,2)+STRDTR(1,2)
      STRTRI(5)=STRVTR(1,3)+STRDTR(1,3)
      STRTRI(6)=STRVTR(2,3)+STRDTR(2,3)
      SNDTRI(4)=STNDTR(1,2)
      SNDTRI(5)=STNDTR(1,3)
      SNDTRI(6)=STNDTR(2,3)
C
C---------WRITING STRESS UPDATE EQUATION----
C
      DO I=1,6
            STRESS(I)=STRTRI(I)-((SIX*DLAM*G**2)/QTR)*SNDTRI(I)
      END DO
      END IF
C--------END----
      RETURN
      END
