module fingreen3D_module

  use fingrnextsubs_module

  implicit none

contains

! ==================================================================================
!
!     Purpose: This program computes the three-dimensional free-surface
!             Green function for finite depth and the derivatives of it
!
!             Code Original Author: Yingyi Liu       created on  2013.09.07
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!         Remark:    The Green function can be decomposed into following form,
!                    where the rankine part can be evaluated separately
!                    based on the method of Hess and Smith (1964) or Newman (1986)
!                    in constant panel method due to its high singularity
!
!                             G= 1/r+1/r'+Gw,
!                    where
!                            r:  dsqrt((x-x')^2+(y-y')^2+(z-z')^2);
!                            r': dsqrt((x-x')^2+(y-y')^2+(z+z')^2).
!  Parameters:
!       Input:      RR  --- horizontal distance between the field and the source point
!                   ZF  --- z,  vertical coordinate of the field point
!                   ZP  --- z', vertical coordinate of the source point
!                    V  --- wave number in deep water,i.e., w^2/g
!                  WVN  --- the first elememnt is the positive root of the dispersion
!                        equation V = k*tanh(k*h), the rest elements are the real
!                        roots of the equation um*tanh(um*h)= -V
!                        (w:circular frequency,g:gravity acceleration)
!                   NK  --- number of elements in the array WVN
!                   H   --- water depth (h>0)
!
!       Output:    GRN(1)  --- value of Green's function
!                  GRN(2)  --- derivative of Green's function with respect to R
!                  GRN(3)  --- derivative of Green's function with respect to z
!
!  Contributors list:
!         Yingyi Liu
!         Matthieu Ancellin
!         to be continued...
!
! ==================================================================================
!!-------------------------------------------------------------------------------!!
!                      Level-1 (Top-level) driver subroutine                      !
!!-------------------------------------------------------------------------------!!

      SUBROUTINE fingreen3d_routine(RR,ZF,ZP,V,WVN,NK,H,GRN)

        INTEGER, INTENT(IN) :: NK
        REAL*8,INTENT(IN) :: RR,ZF,ZP,V,WVN(1:NK),H
        COMPLEX*16,INTENT(OUT) :: GRN(3)

        REAL*8 R,R1,PI4,G,RHP,WK
        COMPLEX*16 FG(3)

        DATA  G,PI4 /9.807D0,12.56637061435917D0/

! Initialize parameters based on the inputs

        R =DSQRT(RR**2+(ZP-ZF)**2)
        R1 =DSQRT(RR**2+(ZP+ZF)**2)

        WK=WVN(1)

! Choose an appropriate method to calculate the Green function
! based on the value of RHP (R/h)
!
        RHP=DABS(RR/H)
         IF (RHP.LT.0.0005D0) THEN
          CALL LINTON(RR,ZF,ZP,V,WK,WVN,NK,H,FG)
        ELSEIF (RHP.GE.0.0005D0.AND.RHP.LT.0.05D0) THEN
          CALL PIDCOCK(RR,ZF,ZP,V,WK,WVN,NK,H,FG)
         ELSEIF (RHP.GE.0.05D0.AND.RHP.LT.0.5D0) THEN
          CALL EIGENE(RR,ZF,ZP,WK,WVN,NK,H,FG)
         ELSE
          CALL EIGEN(RR,ZF,ZP,WK,WVN,NK,H,FG)
        ENDIF

! Post-process the Green function's value
! if R=0, the value of its derivative with respect to R should be zero
!
        GRN(:)=DCMPLX(0.D0,0.D0)

        GRN(1)=FG(1)
        IF(RR.NE.0.0D0) THEN
        GRN(2)=FG(2)
        ELSE
        GRN(2)=DCMPLX(0.D0,0.D0)
        ENDIF
        GRN(3)=FG(3)


        RETURN
        END SUBROUTINE fingreen3d_routine

! ==================================================================================
!   Purpose: This subroutine computes roots of the water-wave dispersion equation
!             in finite water depth, by using a higher-order iterative method
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Author:
!
!    Yingyi Liu on Mar.23, 2017
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!    J.N. Newman
!    Numerical solutions of the water-wave dispersion relation
!    Applied Ocean Research 12 (1990) 14-18
!
!  Parameters:
!      Input:   NRT --- Integer, the number of roots required
!                W   --- Real, wave angular frequency
!                H   --- Real, water depth (h>0)
!      Output:  WVN --- Real, an array storing roots of the dispersion equation
! ==================================================================================

    SUBROUTINE DISPERSION(WVN,NRT,W,H)
!
!   Evaluation of the roots of the following equations
!   by higher-order iterative method
!   first root stored in WVN is from Eq. (i)
!   the rest roots are from Eq. (ii)
!   i) w*w/g = k tanh ( kh )
!   ii) -w*w/g = Um tan ( Umh )
!

      INTEGER,INTENT(IN):: NRT
      REAL*8,INTENT(IN):: W,H
      REAL*8,INTENT(OUT):: WVN(1:NRT)
      INTEGER I, M
      REAL*8 T,X,U,Y,DNM,G,PI
      REAL*8 FUN,DFUN,D2FUN,TRIAL,EXX

      DATA G,PI/9.807d0,3.141592653589793d0/

!------------------------------------------------------------------
! I. calculation of wave number (root of Eq. (i))
!------------------------------------------------------------------
!
!   initialize iteration by an accurate Chebyshev approximation
!   if y=x, use the approximation directly insteady of iteration
!   to avoid the singularity in the denomenator of the transcendental
!   function; otherwise, do the iterative procedure.
!
      X=W*W*H/G
      IF (X.GT.0.D0.AND.X.LE.2.D0) THEN
       Y=DSQRT(X)*(0.9994D0+0.1701D0*X+0.0305*X*X)
      ELSE
       T=X*DEXP(-2.D0*X)
       Y=X+2.D0*T-6.D0*T*T
      ENDIF

      IF (DABS(Y-X).LT.1.E-10) THEN
       WVN(1)=X/H
      ELSE
       M=0
       EXX=1.D0
       DO WHILE (EXX.GT.1.0D-10)
        TRIAL=Y
        DNM=TRIAL*TRIAL-X*X
        FUN=DLOG((TRIAL+X)/(TRIAL-X))/2.D0-TRIAL
        DFUN=-X/DNM-1.D0
        D2FUN=2.D0*X*TRIAL/(DNM*DNM)
        Y=TRIAL-FUN/DFUN*(1.D0+(FUN/DFUN)*(D2FUN/DFUN)/2.D0)
        EXX=DABS(Y-TRIAL)
        M=M+1
       ENDDO
       WVN(1)=Y/H
      ENDIF

!------------------------------------------------------------------
! II. calcultion of roots of Eq. (ii), which characterizes
!     the evanescene modes in eigenfunction
!------------------------------------------------------------------
!
!   initialize iteration by a suitable starting approximation
!
      U=3.D0*X/(7.D0+3.D0*X)
      T=0.0159D0+0.1032D0*U+4.3152D0*U*U-2.8768D0*U*U*U
!
!   perform iterative procedure to find exact solution of Um (m=1,..NRT-1)
!   of the transcendental equation Eq. (ii)
!
      DO I=2,NRT
       M=0
       EXX=1.D0
       DO WHILE (EXX.GT.1.0D-10)
        TRIAL=T
        Y=(I-1)*PI-TRIAL
        DNM=Y*Y+X*X
        FUN=ATAN2(X,Y)-TRIAL
        DFUN=X/DNM-1.D0
        D2FUN=2.D0*X*TRIAL/(DNM*DNM)
        T=TRIAL-FUN/DFUN*(1.D0+(FUN/DFUN)*(D2FUN/DFUN)/2.D0)
        EXX=DABS(T-TRIAL)
        M=M+1
       ENDDO
       Y=(I-1)*PI-T
       WVN(I)=Y/H
       T=T-PI*X/(X*X+PI*I*(PI*(I-1)-T))
      ENDDO

      END SUBROUTINE DISPERSION

!!-------------------------------------------------------------------------------!!
!                Belows are Level-2 (intermediate-level) subroutines              !
!!-------------------------------------------------------------------------------!!

! ==================================================================================
!       Purpose: This program computes Green function and its derivatives
!                using series expansion, without epsilon algorithm, in the region
!                R/H.GT.0.5. The terms needed is 2~10 to 10^-6 accuracy,
!                with a maximum number of term 10.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                V --- wave number in deep water , i.e., V= w^2/g
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!                WVN --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                NK --- number of elements in the array WVN
!                H --- water depth (h>0)
!       Output:  GRN --- Green function value and its derivatives
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!.
! ==================================================================================

        SUBROUTINE EIGEN(R,ZF,ZP,WK,WVN,NK,H,GRN)

         INTEGER I,NK,NT
         REAL*8 J0,J1,Y0,Y1,K0,K1
         REAL*8 PI,KM,PM,RM(3)
         REAL*8 SR,SR1,DSRR,DSRZ,DSR1R,DSR1Z
         REAL*8,INTENT(IN)::R,ZF,ZP,H,WK,WVN(1:NK)

         COMPLEX*16,INTENT(OUT):: GRN(3)
         COMPLEX*16 CI,P0,DP0,NM0

         DATA  PI,CI/3.141592653589793d0,(0.D0,1.0D0)/

! initialize some important parameters
!
         NT=10

         NM0=H/2.D0*(1.D0+DSINH(2.D0*WK*H)/(2.D0*WK*H))
         P0=PI/NM0*DCOSH(WK*(ZF+H))*DCOSH(WK*(ZP+H))
         DP0=PI/NM0*DSINH(WK*(ZF+H))*DCOSH(WK*(ZP+H))

! calculate the imaginary part of the Green function
!
         CALL JY01B(WK*R,J0,J1,Y0,Y1)

         GRN(1)=P0*(CI*J0-Y0)
         GRN(2)=-WK*P0*(CI*J1-Y1)
         GRN(3)=WK*DP0*(CI*J0-Y0)

! calculate the real part of the Green function
! the trunction terms number depends on the value of RHP
!
         DO I=2,NT

         KM=WVN(I)

         PM = 4.D0/(1.D0+DSIN(2.D0*KM*H)/(2.D0*KM*H))/H

         CALL IK01A(KM*R,K0,K1)

         RM(1) = PM*K0*DCOS(KM*(ZF+H))*DCOS(KM*(ZP+H))
         RM(2) = -PM*KM*K1*DCOS(KM*(ZF+H))*DCOS(KM*(ZP+H))
         RM(3) = -PM*KM*K0*DSIN(KM*(ZF+H))*DCOS(KM*(ZP+H))

         GRN(1)= GRN(1)+RM(1)
         GRN(2)= GRN(2)+RM(2)
         GRN(3)= GRN(3)+RM(3)

         IF (DABS(RM(1)).LT.1.E-6.AND.DABS(RM(2)).LT.1.E-6.AND.DABS(RM(3)).LT.1.E-6) THEN
          EXIT
         ENDIF

         ENDDO

! exclude the singular terms for the purpose of analytical integration over panels
!
         SR=DSQRT(R**2+(ZF-ZP)**2)
         SR1=DSQRT(R**2+(ZF+ZP)**2)

         DSRR=R/SR
         DSRZ=(ZF-ZP)/SR
         DSR1R=R/SR1
         DSR1Z=(ZF+ZP)/SR1

         GRN(1)=GRN(1)-(1.D0/SR+1.D0/SR1)
         GRN(2)=GRN(2)+(DSRR/SR**2+DSR1R/SR1**2)
         GRN(3)=GRN(3)+(DSRZ/SR**2+DSR1Z/SR1**2)

        RETURN
        END SUBROUTINE EIGEN

! ==================================================================================
!       Purpose: This program computes Green function and its derivatives
!                using series expansion, with epsilon algorithm, in the region
!                R/H.GT.0.01 and R/H.LT.0.5. The terms needed is 10~80 to
!                10^-5 accuracy, with a maximum number of term 80.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- Horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                V --- wave number in deep water , i.e., V= w^2/g
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!                WVN --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                NK --- number of elements in the array WVN
!                H --- water depth (h>0)
!       Output:  GRN --- Green function value and its derivatives
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

        SUBROUTINE EIGENE(R,ZF,ZP,WK,WVN,NK,H,GRN)

         INTEGER I,NK,NT,I_PADE,K_PADE
         REAL*8 J0,J1,Y0,Y1,K0,K1
         REAL*8 PI,KM,PM,RM(3),RL(3),ERR,RHP
         REAL*8 SR,SR1,DSRR,DSRZ,DSR1R,DSR1Z
         REAL*8 G(0:NK-1),GR(0:NK-1),GZ(0:NK-1)
         REAL*8,INTENT(IN)::R,ZF,ZP,H,WK,WVN(1:NK)

         COMPLEX*16,INTENT(OUT):: GRN(3)
         COMPLEX*16 CI,P0,DP0,NM0

         DATA  PI,CI/3.141592653589793d0,(0.D0,1.0D0)/

! initialize some important parameters
!
         RHP=DABS(R/H)
         NT=INT(-88.89*RHP+54.45)

         NM0=H/2.D0*(1.D0+DSINH(2.D0*WK*H)/(2.D0*WK*H))
         P0=PI/NM0*DCOSH(WK*(ZF+H))*DCOSH(WK*(ZP+H))
         DP0=PI/NM0*DSINH(WK*(ZF+H))*DCOSH(WK*(ZP+H))

! calculate the imaginary part of the Green function
!
         CALL JY01B(WK*R,J0,J1,Y0,Y1)

         GRN(1)=P0*(CI*J0-Y0)
         GRN(2)=-WK*P0*(CI*J1-Y1)
         GRN(3)=WK*DP0*(CI*J0-Y0)

! calculate the real part of the Green function
! the trunction terms number depends on the value of RHP
!
         G(0)= 0.D0
         GR(0)= 0.D0
         GZ(0)= 0.D0

         DO I=2,NT

         KM=WVN(I)

         PM = 4.D0/(1.D0+DSIN(2.D0*KM*H)/(2.D0*KM*H))/H

         CALL IK01A(KM*R,K0,K1)

         RM(1) = PM*K0*DCOS(KM*(ZF+H))*DCOS(KM*(ZP+H))
         RM(2) = -PM*KM*K1*DCOS(KM*(ZF+H))*DCOS(KM*(ZP+H))
         RM(3) = -PM*KM*K0*DSIN(KM*(ZF+H))*DCOS(KM*(ZP+H))

         G(I-1)= G(I-2)+RM(1)
         GR(I-1)= GR(I-2)+RM(2)
         GZ(I-1)= GZ(I-2)+RM(3)

         ENDDO

! call the epsilon algorithm to approximate the limit of the infinite series
!
         CALL Limes (NT-1,G(0:NT-1),RL(1),I_PADE,K_PADE,ERR)
         CALL Limes (NT-1,GR(0:NT-1),RL(2),I_PADE,K_PADE,ERR)
         CALL Limes (NT-1,GZ(0:NT-1),RL(3),I_PADE,K_PADE,ERR)

         GRN(1)= GRN(1)+RL(1)
         GRN(2)= GRN(2)+RL(2)
         GRN(3)= GRN(3)+RL(3)

! exclude the singular terms for the purpose of analytical integration over panels
!
         SR=DSQRT(R**2+(ZF-ZP)**2)
         SR1=DSQRT(R**2+(ZF+ZP)**2)

         DSRR=R/SR
         DSRZ=(ZF-ZP)/SR
         DSR1R=R/SR1
         DSR1Z=(ZF+ZP)/SR1

         GRN(1)=GRN(1)-(1.D0/SR+1.D0/SR1)
         GRN(2)=GRN(2)+(DSRR/SR**2+DSR1R/SR1**2)
         GRN(3)=GRN(3)+(DSRZ/SR**2+DSR1Z/SR1**2)

        RETURN
        END SUBROUTINE EIGENE


! ==================================================================================
!       Purpose: This program computes Green function and its derivatives
!                using Pidcock's method, with epsilon algorithm, in the region
!                R/H.GT.0.0005 and R/H.LT.0.01. The terms needed is 50~100 to
!                10^-4 accuracy, with a maximum number of term 100.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                V --- wave number in deep water , i.e., V= w^2/g
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!                WVN --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                NK --- number of elements in the array WVN
!                H --- water depth (h>0)
!       Output:  GRN --- Green function value and its derivatives
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

        SUBROUTINE PIDCOCK(R,ZF,ZP,V,WK,WVN,NK,H,GRN)

         INTEGER NK,NT,I,I_PADE,K_PADE

         REAL*8 J0,J1,Y0,Y1,K0,K0T,K1,K1T
         REAL*8 NM0,P0,DP0,X,PI,GA
         REAL*8 RM(3),SM(3),RL(3),ERR,RHP
         REAL*8 UM,UMT,MP,AC,AM,AMT,DAM,DAMT
         REAL*8 SR,SR1,DSRR,DSRZ,DSR1R,DSR1Z
         REAL*8 G(0:NK-1),GR(0:NK-1),GZ(0:NK-1)
         REAL*8,INTENT(IN)::R,ZF,ZP,H,V,WK,WVN(1:NK)
         !REAL*8,EXTERNAL:: IMGS,DGSR,DGSZ

         COMPLEX*16,INTENT(OUT):: GRN(3)
         COMPLEX*16 CI

         DATA PI,GA /3.141592653589793D0,0.5772156649015328D0/
         DATA CI/(0.0D0, 1.0D0)/

! initialize some important parameters
!
         NM0=H/2.D0*(1.D0+DSINH(2.D0*WK*H)/(2.D0*WK*H))
         P0=PI/NM0*DCOSH(WK*(ZF+H))*DCOSH(WK*(ZP+H))
         DP0=WK*PI/NM0*DSINH(WK*(ZF+H))*DCOSH(WK*(ZP+H))

! calculate the imaginary part of the Green function
!
         CALL JY01B(WK*R,J0,J1,Y0,Y1)

         GRN(1)=P0*(CI*J0-Y0)+2.D0/H*(GA+DLOG(R/4.D0/H))
         GRN(2)=-WK*P0*(CI*J1-Y1)+2.D0/(R*H)
         GRN(3)=DP0*(CI*J0-Y0)

! calculate the real part, for the infinite series of Eq.(13)
! the series applies an acceleration method by subtracting an asymptotic series
!
         RHP=DABS(R/H)
         NT=INT(-1010.10*RHP+100.50)

         G(0)= 0.D0
         GR(0)= 0.D0
         GZ(0)= 0.D0

         DO I=2,NT

          MP=(I-1)*PI
          UMT=MP/H
          UM=WVN(I)
          X=UM*R

          CALL IK01A(X,K0,K1)
          CALL IK01A(UMT*R,K0T,K1T)

          AC=(UM**2+V**2)/((UM**2+V**2)*H-V)
          AM=AC*DCOS(UM*(ZF+H))*DCOS(UM*(ZP+H))
          AMT=DCOS(UMT*(ZF+H))*DCOS(UMT*(ZP+H))/H
          DAM=-UM*AC*DSIN(UM*(ZF+H))*DCOS(UM*(ZP+H))
          DAMT=-UMT*DSIN(UMT*(ZF+H))*DCOS(UMT*(ZP+H))/H

          RM(1)= AM*K0-AMT*K0T
          RM(2)= -UM*AM*K1+UMT*AMT*K1T
          RM(3)= DAM*K0-DAMT*K0T

          G(I-1)= G(I-2)+RM(1)
          GR(I-1)= GR(I-2)+RM(2)
          GZ(I-1)= GZ(I-2)+RM(3)

         ENDDO

! call the epsilon algorithm to approximate the limit of the infinite series
!
         CALL Limes (NT-1,G(0:NT-1),RL(1),I_PADE,K_PADE,ERR)
         CALL Limes (NT-1,GR(0:NT-1),RL(2),I_PADE,K_PADE,ERR)
         CALL Limes (NT-1,GZ(0:NT-1),RL(3),I_PADE,K_PADE,ERR)

         SM(1)= 4.D0*RL(1)
         SM(2)= 4.D0*RL(2)
         SM(3)= 4.D0*RL(3)

! calculate the sum of Rankine terms by Chebyshev approximations
!
         GRN(1)=GRN(1)+SM(1)+IMGS(R,ZF,ZP,H)
         GRN(2)=GRN(2)+SM(2)+DGSR(R,ZF,ZP,H)
         GRN(3)=GRN(3)+SM(3)+DGSZ(R,ZF,ZP,H)

! exclude the singular terms for the purpose of analytical integration over panels
!
         SR=DSQRT(R**2+(ZF-ZP)**2)
         SR1=DSQRT(R**2+(ZF+ZP)**2)

         DSRR=R/SR
         DSRZ=(ZF-ZP)/SR
         DSR1R=R/SR1
         DSR1Z=(ZF+ZP)/SR1

         GRN(1)=GRN(1)-(1.D0/SR+1.D0/SR1)
         GRN(2)=GRN(2)+(DSRR/SR**2+DSR1R/SR1**2)
         GRN(3)=GRN(3)+(DSRZ/SR**2+DSR1Z/SR1**2)

        RETURN
        END SUBROUTINE PIDCOCK


! ==================================================================================
!       Purpose: This program computes Green function and its derivatives
!                using Linton's method, with epsilon algorithm, in the region
!                R/H.LT.0.0005. The terms needed is about 3~10 to
!                10^-6 accuracy, with a maximum number of term 10.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                V --- wave number in deep water , i.e., V= w^2/g
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!                WVNO --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                NK --- number of elements in the array WVNO
!                H --- water depth (h>0)
!       Output:  GRN --- Green function value and its derivatives
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

        SUBROUTINE LINTON(R,ZF,ZP,V,WK,WVNO,NK,H,GRN)

         INTEGER NK,NT,I,FLAG
         REAL*8 PA,A,PI,EPS,SQPI
         REAL*8 SR,SR1,SR2,DSRR,DSRZ,DSR1R,DSR1Z,DSR2R,DSR2Z
         REAL*8 RO(4),DROR(4),DROZ(4),ERRSUM(3),BJ0,BJ1
         REAL*8 COF(NK),DCOF(NK)
         REAL*8,INTENT(IN)::R,ZF,ZP,H,V,WK,WVNO(NK)
         COMPLEX*16,INTENT(OUT):: GRN(3)
         COMPLEX*16 KM,RM(3),CI
         !REAL*8,EXTERNAL:: ERFCC,DERFCC,AQUADF,G3,G1,G2,FAC,F2
         !COMPLEX*16,EXTERNAL::NM

        DATA PI,EPS/3.141592653589793D0,1.E-6/
        DATA CI/(0.0D0, 1.0D0)/

        SQPI=DSQRT(PI)

!-----------------------------------------------------------------------------
!  I. Optimal selection of the parameter A
!-----------------------------------------------------------------------------

        FLAG=2

        IF (FLAG.EQ.1) THEN
          A=0.25D0
        ELSEIF(FLAG.EQ.2) THEN
          IF (V.LT.0.02D0) THEN
            A=3324.2D0*V**3-258.89D0*V**2+8.3958D0*V+1.3208D0
          ELSEIF (V.LT.0.1D0) THEN
            A=-30.303D0*V**3-2.8138D0*V**2+2.7740D0*V+0.05052D0
          ELSEIF (V.LT.1.0D0) THEN
            A=1.5699D0*V**3-2.5075D0*V**2+1.8316D0*V+0.1077D0
          ELSE
            A=-0.000069641D0*V**3-0.0017082D0*V**2+0.45566D0*V+0.80551D0
          ENDIF
           A=DABS(A)/(WK*H)
        ENDIF
!
!       PRINT*,A
!
!-----------------------------------------------------------------------------
!  II. Initialisation of various parameters
!-----------------------------------------------------------------------------

         PA=A*H
         SR=DSQRT(R**2+(ZF-ZP)**2)
         SR1=DSQRT(R**2+(2.D0*H+ZF+ZP)**2)
         SR2=DSQRT(R**2+(ZF+ZP)**2)

         CALL RCHI(R,ZF,ZP,H,RO,DROR,DROZ)

         DSRR=R/SR
         DSRZ=(ZF-ZP)/SR
         DSR1R=R/SR1
         DSR1Z=(2.D0*H+ZF+ZP)/SR1
         DSR2R=R/SR2
         DSR2Z=(ZF+ZP)/SR2
!
!-----------------------------------------------------------------------------
!  III. Calculating the part of evanescent modes
!-----------------------------------------------------------------------------

         GRN(:) = (0.0D0, 0.0D0)
         NT=20

         CALL COEF(R,PA,A,NT,WVNO,WK,COF)
         CALL DCOEF(R,PA,A,NT,WVNO,WK,DCOF)

         DO I=1,NT

         IF (I.EQ.1) THEN
          KM=-CI*WK
         ELSE
          KM=WVNO(I)
         ENDIF

         RM(1)= -COF(I)/NM(H,KM)*CDCOS(KM*(ZF+H))*CDCOS(KM*(ZP+H))
         RM(2)= -DCOF(I)/NM(H,KM)*CDCOS(KM*(ZF+H))*CDCOS(KM*(ZP+H))
         RM(3)= COF(I)*KM/NM(H,KM)*CDSIN(KM*(ZF+H))*CDCOS(KM*(ZP+H))

         GRN(1)= GRN(1)+RM(1)
         GRN(2)= GRN(2)+RM(2)
         GRN(3)= GRN(3)+RM(3)


         IF (I.GT.2.AND.CDABS(RM(1)).LT.EPS.AND.CDABS(RM(2)).LT.EPS.AND.CDABS(RM(3)).LT.EPS) THEN
         EXIT
         ENDIF

         ENDDO
!
!-----------------------------------------------------------------------------
!  IV. Calculating the part of Error functions
!-----------------------------------------------------------------------------

         ERRSUM(:)=0.D0

         DO I=1,4

         ERRSUM(1)=ERRSUM(1)+ERFCC(RO(I)/PA)/(4.D0*PI*RO(I))
         ERRSUM(2)=ERRSUM(2)+(DERFCC(RO(I)/PA)/PA*DROR(I)*RO(I)-DROR(I)*ERFCC(RO(I)/PA))/(4.D0*PI*RO(I)**2)
         ERRSUM(3)=ERRSUM(3)+(DERFCC(RO(I)/PA)/PA*DROZ(I)*RO(I)-DROZ(I)*ERFCC(RO(I)/PA))/(4.D0*PI*RO(I)**2)

         ENDDO


       IF (DABS(SR).GT.1.E-6) THEN

         ERRSUM(1)=ERRSUM(1)+(ERFCC(SR/PA)-1.D0)/(4.D0*PI*SR)+ERFCC(SR1/PA)/(4.D0*PI*SR1)-1.D0/(4.D0*PI*SR2)
         ERRSUM(2)=ERRSUM(2)-(ERFCC(SR/PA)-1.D0)/(4.D0*PI*SR**2)*DSRR+DERFCC(SR/PA)/(4.D0*PI*SR*PA)*DSRR-ERFCC(SR1/PA)    &
                   /(4.D0*PI*SR1**2)*DSR1R+DERFCC(SR1/PA)/(4.D0*PI*SR1*PA)*DSR1R+DSR2R/(4.D0*PI*SR2**2)
         ERRSUM(3)=ERRSUM(3)-(ERFCC(SR/PA)-1.D0)/(4.D0*PI*SR**2)*DSRZ+DERFCC(SR/PA)/(4.D0*PI*SR*PA)*DSRZ-ERFCC(SR1/PA)    &
                   /(4.D0*PI*SR1**2)*DSR1Z+DERFCC(SR1/PA)/(4.D0*PI*SR1*PA)*DSR1Z+DSR2Z/(4.D0*PI*SR2**2)

       ELSE

         ERRSUM(1)=ERRSUM(1)+ERFCC(SR1/PA)/(4.D0*PI*SR1)-1.D0/(4.D0*PI*SR2)
         ERRSUM(2)=ERRSUM(2)-ERFCC(SR1/PA)/(4.D0*PI*SR1**2)*DSR1R+DERFCC(SR1/PA)/(4.D0*PI*SR1*PA)*DSR1R+DSR2R/(4.D0*PI*SR2**2)
         ERRSUM(3)=ERRSUM(3)-ERFCC(SR1/PA)/(4.D0*PI*SR1**2)*DSR1Z+DERFCC(SR1/PA)/(4.D0*PI*SR1*PA)*DSR1Z+DSR2Z/(4.D0*PI*SR2**2)

       ENDIF
!
!-----------------------------------------------------------------------------
!  V. Calculating the part of integration
!-----------------------------------------------------------------------------

         GRN(1)= GRN(1)-ERRSUM(1)-V/2.D0/PI*AQUADF(V,R,ZF,ZP,H,0.D0,PA/2.D0,EPS,1)
         GRN(2)= GRN(2)-ERRSUM(2)+V/2.D0/PI*AQUADF(V,R,ZF,ZP,H,0.D0,PA/2.D0,EPS,2)
         GRN(3)= GRN(3)-ERRSUM(3)-V/2.D0/PI*AQUADF(V,R,ZF,ZP,H,0.D0,PA/2.D0,EPS,3)


         CALL JY01BJ(WK*R,BJ0,BJ1)


         KM=-CI*WK
         GRN(1)= GRN(1)-CI/4.D0/NM(H,KM)*BJ0*DCOSH(WK*(ZF+H))*DCOSH(WK*(ZP+H))
         GRN(2)= GRN(2)+CI*WK/4.D0/NM(H,KM)*BJ1*DCOSH(WK*(ZF+H))*DCOSH(WK*(ZP+H))
         GRN(3)= GRN(3)-CI*WK/4.D0/NM(H,KM)*BJ0*DSINH(WK*(ZF+H))*DCOSH(WK*(ZP+H))


         GRN(1)= -4.D0*PI*GRN(1)
         GRN(2)= -4.D0*PI*GRN(2)
         GRN(3)= -4.D0*PI*GRN(3)


        RETURN
        END SUBROUTINE LINTON

!!-------------------------------------------------------------------------------!!
!             Belows are Level-3 (low level) self-contained subroutines           !
!!-------------------------------------------------------------------------------!!

! ==================================================================================
!       Purpose: This program computes Rankine image source
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                H --- water depth (h>0)
!       Output:  IMGS --- value of the Rankine image source in Pidcock's method
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION IMGS(R,ZF,ZP,H)

       INTEGER L,M,N
       REAL*8 R,ZF,ZP,H,AMN(0:4,0:4)
       REAL*8 T1,T2,T3,SUM1,SUM2
       !REAL*8,EXTERNAL::RS

       DATA AMN /0.D0,-0.02525711D0,0.00086546D0,-0.00004063D0,        &
            0.00000193D0,0.05051418D0,-0.00692380D0,0.00073292D0,     &
            -0.00006636D0,0.00000398D0,0.00230838D0,-0.00097875D0,     &
            0.00020597D0,-0.00003333D0,0.00000524D0,0.00012934D0,      &
            -0.00010879D0,0.00003965D0,-0.00000891D0,0.0D0,          &
            0.00000913D0,-0.00001270D0,0.00000466D0,0.0D0,0.D0/

       T1=R/H
       T2=(ZF-ZP)/H
       T3=(ZF+ZP+2.D0*H)/H

       SUM1=0.D0
       DO L=-1,1
        SUM1=SUM1+1.D0/RS(T1,T2+2.D0*L)+1.D0/RS(T1,T3+2.D0*L)
       ENDDO

       SUM2=0.D0
       DO M=0,4
        DO N=0,4
        SUM2=SUM2+AMN(M,N)*T1**(2.D0*M)*(T2**(2.D0*N)+T3**(2.D0*N))
       ENDDO
       ENDDO

       IMGS=(SUM1-2.D0+SUM2)/H

       RETURN
       END FUNCTION IMGS


! ==================================================================================
!       Purpose: This program computes the derivative of Rankine image source
!                with respect to R
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                H --- water depth (h>0)
!       Output:  DGSR --- derivative of the Rankine image source with respect to R
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION DGSR(R,ZF,ZP,H)

       INTEGER L,M,N
       REAL*8 R,ZF,ZP,H,AMN(0:4,0:4)
       REAL*8 T1,T2,T3,SUM1,SUM2
       !REAL*8,EXTERNAL::RS

       DATA AMN /0.D0,-0.02525711D0,0.00086546D0,-0.00004063D0,        &
            0.00000193D0,0.05051418D0,-0.00692380D0,0.00073292D0,     &
            -0.00006636D0,0.00000398D0,0.00230838D0,-0.00097875D0,     &
            0.00020597D0,-0.00003333D0,0.00000524D0,0.00012934D0,      &
            -0.00010879D0,0.00003965D0,-0.00000891D0,0.0D0,          &
            0.00000913D0,-0.00001270D0,0.00000466D0,0.0D0,0.D0/

       T1=R/H
       T2=(ZF-ZP)/H
       T3=(ZF+ZP+2.D0*H)/H

       SUM1=0.D0
       DO L=-1,1
        SUM1=SUM1+T1/(RS(T1,T2+2.D0*L))**3+T1/(RS(T1,T3+2.D0*L))**3
       ENDDO

       SUM2=0.D0
       DO M=1,4
        DO N=0,4
        SUM2=SUM2+2.D0*M*AMN(M,N)*T1**(2.D0*M-1.D0)*(T2**(2.D0*N)+T3**(2.D0*N))
       ENDDO
       ENDDO

       DGSR=(-SUM1+SUM2)/H**2

       RETURN
       END FUNCTION DGSR

! ==================================================================================
!       Purpose: This program computes the derivative of Rankine image source
!                with respect to z
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                H --- water depth (h>0)
!       Output:  DGSZ --- derivative of the Rankine image source with respect to z
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION DGSZ(R,ZF,ZP,H)

       INTEGER L,M,N
       REAL*8 R,ZF,ZP,H,AMN(0:4,0:4)
       REAL*8 T1,T2,T3,SUM1,SUM2
       !REAL*8,EXTERNAL::RS

       DATA AMN /0.D0,-0.02525711D0,0.00086546D0,-0.00004063D0,        &
            0.00000193D0,0.05051418D0,-0.00692380D0,0.00073292D0,     &
            -0.00006636D0,0.00000398D0,0.00230838D0,-0.00097875D0,     &
            0.00020597D0,-0.00003333D0,0.00000524D0,0.00012934D0,      &
            -0.00010879D0,0.00003965D0,-0.00000891D0,0.0D0,          &
            0.00000913D0,-0.00001270D0,0.00000466D0,0.0D0,0.D0/

       T1=R/H
       T2=(ZF-ZP)/H
       T3=(ZF+ZP+2.D0*H)/H

       SUM1=0.D0
       DO L=-1,1
        SUM1=SUM1+(T2+2.D0*L)/(RS(T1,T2+2.D0*L))**3+(T3+2.D0*L)/(RS(T1,T3+2.D0*L))**3
       ENDDO

       SUM2=0.D0
       DO M=0,4
        DO N=1,4
        SUM2=SUM2+2.D0*N*AMN(M,N)*T1**(2.D0*M)*(T2**(2.D0*N-1.D0)+T3**(2.D0*N-1.D0))
       ENDDO
       ENDDO

       DGSZ=(-SUM1+SUM2)/H**2

       RETURN
       END FUNCTION DGSZ

!
!    ======================================================
!          Function Dsqrt(R^2+Z^2)
!    ======================================================

       REAL*8 FUNCTION RS(R,Z)

       REAL*8 R,Z

       RS=DSQRT(R**2+Z**2)


       RETURN
       END FUNCTION RS


! ==================================================================================
!       Purpose: This program computes the expansion coefficients in
!                Linton's method
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                PA --- a key parameter determining convergence, i.e., PA= a*h
!                A ---  a key parameter related with convergence
!                M --- number of elements in the array WVNO
!                WVNO --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!       Output:  COF --- expansion coefficients in Linton's expression
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       SUBROUTINE COEF(R,PA,A,M,WVNO,WK,COF)

        INTEGER I
        INTEGER,INTENT(IN)::M
        REAL*8,INTENT(IN)::R,PA,A,WVNO(M),WK
        REAL*8,INTENT(OUT)::COF(M)
        !REAL*8,EXTERNAL:: SA0,AQUAD,SAM,F1
        REAL*8 BY0,DNT,EPS,STEP,P1,E1,EI
        REAL*8 PI

        DATA  PI/3.141592653589793D0/

        EPS=1.E-6
        STEP=100.D0
        DNT=PA**2/4.D0

! calculate the first expansion coefficient A0
!
       IF (DABS(R).LT.1.E-6) THEN

         CALL EIX(WK**2*DNT,EI)
         COF(1)=-EI/4.D0/PI

       ELSEIF (R/A.LE.0.5D0) THEN

         COF(1)=SA0(R,PA,WK)

       ELSE

         CALL JY01BY0(WK*R,BY0)
         COF(1)=-BY0/4.D0-AQUAD(R,WK,1,0.0D0,DNT,EPS,1)

       ENDIF

! calculate the rest expansion coefficients Ai (i=2,m)
! the truncation terms number depends on whether the series
! converged to a specified accuracy or not

       DO I=2,M

        IF (DABS(R).LT.1.E-6) THEN

         CALL E1XA(WVNO(I)**2*DNT,E1)
         COF(I)=E1/4.D0/PI

        ELSEIF (R/A.LE.1.D0) THEN

         COF(I)=SAM(R,PA,WVNO(I))

        ELSE

         P1=DNT+STEP
         DO WHILE (DABS(F1(R,WVNO(I),I,P1)).GT.EPS)
         P1=P1+STEP
         ENDDO

         COF(I)=AQUAD(R,WVNO(I),I,DNT,P1,EPS,1)

        ENDIF

       ENDDO


       RETURN
       END SUBROUTINE COEF

! ==================================================================================
!       Purpose: This program computes derivative of the expansion coefficients
!                with respect to R in Linton's method
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                PA --- a key parameter determining convergence, i.e., PA= a*h
!                A ---  a key parameter related with convergence
!                M --- number of elements in the array WVNO
!                WVNO --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!       Output:  DCOEF --- derivative of expansion coefficients with respect to R
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       SUBROUTINE DCOEF(R,PA,A,M,WVNO,WK,DCOF)

        INTEGER I
        INTEGER,INTENT(IN)::M
        REAL*8,INTENT(IN)::R,PA,A,WVNO(M),WK
        REAL*8,INTENT(OUT)::DCOF(M)
        !REAL*8,EXTERNAL:: DSA0,DSAM,AQUAD,F2
        REAL*8 BY1,DNT,EPS,STEP,P1

        EPS=1.E-6
        STEP=100.D0
        DNT=PA**2/4.D0

! calculate the derivative of the first expansion coefficient A0
!
        IF (DABS(R).LT.1.E-6) THEN

          DCOF(1)=0.D0

        ELSEIF (R/A.LE.0.5D0) THEN

          DCOF(1)=DSA0(R,PA,WK)

        ELSE

          CALL JY01BY1(WK*R,BY1)
          DCOF(1)=WK*BY1/4.D0-AQUAD(R,WK,1,0.0D0,DNT,EPS,2)

        ENDIF

! calculate the derivative of the rest expansion coefficients Ai (i=2,m)
! the truncation terms number depends on whether the series
! converged to a specified accuracy or not

       DO I=2,M

        IF (DABS(R).LT.1.E-6) THEN

         DCOF(I)=0.D0

        ELSEIF (R/A.LE.1.D0) THEN

         DCOF(I)=DSAM(R,PA,WVNO(I))

        ELSE

         P1=DNT+STEP
         DO WHILE (DABS(F2(R,WVNO(I),I,P1)).GT.EPS)
         P1=P1+STEP
         ENDDO

         DCOF(I)=AQUAD(R,WVNO(I),I,DNT,P1,EPS,2)

        ENDIF

       ENDDO

       RETURN
       END SUBROUTINE DCOEF

! ==================================================================================
!       Purpose: This program computes the series representations of expansion
!                coefficients (m=0)
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 17, 2018
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                PA --- a key parameter determining convergence, i.e., PA= a*h
!                A ---  a key parameter related with convergence
!                M --- number of elements in the array WVNO
!                WVNO --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!       Output:  SA0 --- value of first expansion coefficient by series representation (m=0)
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================


       REAL*8 FUNCTION SA0(R,PA,WK)

       REAL*8,INTENT(IN)::R,PA,WK
       INTEGER N,P
       REAL*8 GA,PI,EPS,PSI,RM,EN,SGN,FACN
       !REAL*8,EXTERNAL::SYMB,FAC

       DATA GA,PI,EPS/0.5772156649015328D0 ,3.141592653589793D0,1.E-6/

       N=1
       RM=100.D0
       SA0=-GA-2.D0*DLOG(WK*PA/2.D0)
       SGN=-1.D0
       FACN=1.D0

! Loop: the summation stops at the condition when RM meets the tolerance
!
       DO WHILE (DABS(RM).GT.EPS)

        CALL ENXA(N+1,(R/PA)**2,EN)

        PSI=0.D0
        DO P=1,N
         PSI=PSI+1.D0/P
        ENDDO

! calculate the value of each term in the infinite series expansion
!
        RM=SGN*(R/PA)**(2.D0*N)/FACN/N-(WK*PA/2.D0)**(2.D0*N)/FACN*EN-2.D0*SGN     &
            /FACN**2*(WK*R/2.D0)**(2.D0*N)*(DLOG(WK*R/2.D0)+GA-PSI)

        !RM=-(WK*PA/2.D0)**(2.D0*N)/FACN*EN

        SA0=SA0+RM

        N=N+1
        SGN=-SGN   ! Sign function depends on the power n
        FACN=FACN*N   ! factorial of the numbers from 1 to n

       ENDDO

        SA0=SA0/4.D0/PI

       RETURN
       END FUNCTION SA0

! ==================================================================================
!       Purpose: This program computes the series representations of expansion
!                coefficients with respect to R in Linton's method (m=0)
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 17, 2018
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                PA --- a key parameter determining convergence, i.e., PA= a*h
!                A ---  a key parameter related with convergence
!                M --- number of elements in the array WVNO
!                WVNO --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!       Output:  DSA0 --- derivative value of first expansion coefficient by series
!                        representation (m=0)
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION DSA0(R,PA,WK)

       REAL*8,INTENT(IN)::R,PA,WK
       INTEGER N,P
       REAL*8 GA,PI,EPS,PSI,RM,EN,SGN,FACN
       !REAL*8,EXTERNAL::SYMB,FAC

       DATA GA,PI,EPS/0.5772156649015328D0 ,3.141592653589793D0,1.E-6/

       N=1
       RM=100.D0
       DSA0=0.D0
       SGN=-1.D0
       FACN=1.D0

! Loop: the summation stops at the condition when RM meets the tolerance
!
       DO WHILE (DABS(RM).GT.EPS)

        CALL ENXA(N,(R/PA)**2,EN)

        PSI=0.D0
        DO P=1,N
         PSI=PSI+1.D0/P
        ENDDO

! calculate the value of each term in the infinite series expansion
!
        RM=2.D0*SGN*(R/PA)**(2.D0*N-1.D0)/FACN/PA+(WK*PA/2.D0)**(2.D0*N)/FACN*EN*2.D0*R/PA**2-2.D0*SGN*N*WK     &
           /FACN**2*(WK*R/2.D0)**(2.D0*N-1.D0)*(DLOG(WK*R/2.D0)+GA-PSI)-2.D0*SGN/FACN**2*(WK*R/2.D0)**(2.D0*N)/R

        DSA0=DSA0+RM

        N=N+1
        SGN=-SGN   ! Sign function depends on the power n
        FACN=FACN*N   ! factorial of the numbers from 1 to n

       ENDDO

        DSA0=DSA0/4.D0/PI

       RETURN
       END FUNCTION DSA0

! ==================================================================================
!       Purpose: This program computes the series representations of expansion
!                coefficients (m>=1,m=1,2,3,4...)
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 17, 2018
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                PA --- a key parameter determining convergence, i.e., PA= a*h
!                A ---  a key parameter related with convergence
!                M --- number of elements in the array WVNO
!                WVNO --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!       Output:  SAM --- value of first expansion coefficient by series representation
!                       (m>=1,m=1,2,3,4...)
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION SAM(R,PA,WVN)

       REAL*8,INTENT(IN)::R,PA,WVN
       INTEGER N
       REAL*8 GA,PI,EPS,RM,EN,SGN,FACN
       !REAL*8,EXTERNAL::SYMB,FAC

       DATA GA,PI,EPS/0.5772156649015328D0 ,3.141592653589793D0,1.E-6/

       N=0
       RM=100.D0
       SAM=0.D0
       SGN=1.D0
       FACN=1.D0

! Loop: the summation stops at the condition when RM meets the tolerance
!
       DO WHILE (DABS(RM).GT.EPS)

        CALL ENXA(N+1,(WVN*PA/2.D0)**2,EN)

        RM=SGN*(R/PA)**(2.D0*N)/FACN*EN

        SAM=SAM+RM

        N=N+1
        SGN=-SGN   ! Sign function depends on the power n
        FACN=FACN*N   ! factorial of the numbers from 1 to n

       ENDDO

        SAM=SAM/4.D0/PI

       RETURN
       END FUNCTION SAM

! ==================================================================================
!       Purpose: This program computes the series representations of expansion
!                coefficients with respect to R in Linton's method (m>=1,m=1,2,3,4...)
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 17, 2018
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                PA --- a key parameter determining convergence, i.e., PA= a*h
!                A ---  a key parameter related with convergence
!                M --- number of elements in the array WVNO
!                WVNO --- The first elememnt is WK, the rest elements are
!                        the real roots of the equation um*tanh(um*h)= -V
!                WK --- positive root of the dispersion equation k*tanh(k*h)= V
!       Output:  DSAM --- derivative value of first expansion coefficient by series
!                        representation (m>=1,m=1,2,3,4...)
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION DSAM(R,PA,WVN)

       REAL*8,INTENT(IN)::R,PA,WVN
       INTEGER N
       REAL*8 GA,PI,EPS,RM,EN,SGN,FACN
       !REAL*8,EXTERNAL::SYMB,FAC

       DATA GA,PI,EPS/0.5772156649015328D0 ,3.141592653589793D0,1.E-6/

       N=1
       RM=100.D0
       DSAM=0.D0
       SGN=-1.D0
       FACN=1.D0

! Loop: the summation stops at the condition when RM meets the tolerance
!
       DO WHILE (DABS(RM).GT.EPS)

        CALL ENXA(N+1,(WVN*PA/2.D0)**2,EN)

        RM=2.D0*N*SGN*(R/PA)**(2.D0*N-1.D0)/FACN/PA*EN

        DSAM=DSAM+RM

        N=N+1
        SGN=-SGN   ! Sign function depends on the power n
        FACN=FACN*N   ! factorial of the numbers from 1 to n

       ENDDO

       DSAM=DSAM/4.D0/PI

       RETURN
       END FUNCTION DSAM

! ==================================================================================
!       Purpose: This program computes the expansion coefficients and its
!                derivative in Linton's method by integral form.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                WVN --- The input wave number
!                NO --- the term number to be calculated in Linton's expansion
!                T ---  integration variable
!                FUNTAG --- 1: calculate function value; 2: calculate derivative value
!       Output:  FUN --- value of expansion coefficient or its derivative by integral
!                        representation
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

        REAL*8 FUNCTION FUN(R,WVN,NO,T,FUNTAG)

        INTEGER,INTENT(IN)::NO,FUNTAG
        REAL*8,INTENT(IN)::R,WVN,T
        !REAL*8,EXTERNAL::F1,F2

        SELECT CASE(FUNTAG)


        CASE(1)
          FUN=F1(R,WVN,NO,T)

        CASE(2)
          FUN=F2(R,WVN,NO,T)


        END SELECT

        RETURN
        END FUNCTION FUN


! ==================================================================================
!       Purpose: This program computes the expansion coefficients and its
!                derivative in Linton's method by integral form.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                WVN --- The input wave number
!                NO --- the term number to be calculated in Linton's expansion
!                T ---  integration variable
!                FUNTAG --- 1: calculate function value; 2: calculate derivative value
!       Output:  FUN --- value of expansion coefficient or its derivative by integral
!                        representation
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

        REAL*8 FUNCTION FUNF(V,R,ZF,ZP,H,T,FUNTAG)

        INTEGER,INTENT(IN)::FUNTAG
        REAL*8,INTENT(IN)::V,R,ZF,ZP,H,T
        !REAL*8,EXTERNAL::G1,G2,G3

        SELECT CASE(FUNTAG)


!        INTEGRAND OF L1
        CASE(1)
          FUNF=G1(V,R,ZF,ZP,H,T)

!        INTEGRAND OF L1R
        CASE(2)
          FUNF=G2(V,R,ZF,ZP,H,T)

!        INTEGRAND OF L1Z
        CASE(3)
          FUNF=G3(V,R,ZF,ZP,H,T)


        END SELECT

        RETURN
        END FUNCTION FUNF

! ==================================================================================
!       Purpose: This program computes the value of expansion coefficients
!                in Linton's method by integral form.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                WVN --- The input wave number
!                NO --- the term number to be calculated in Linton's expansion
!                T ---  integration variable
!       Output:  F1 --- value of expansion coefficient by integral representation
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION F1(R,WVN,NO,T)

       INTEGER NO
       REAL*8,INTENT(IN)::R,WVN,T
       REAL*8 PI

       DATA  PI/3.141592653589793D0/

       IF (DABS(T).LT.1.D-8) THEN

        F1=0.D0

       ELSE

        IF (NO.EQ.1) THEN

        F1=DEXP(-R**2/4.D0/T+WVN**2*T)/4.D0/PI/T

        ELSE

        F1=DEXP(-R**2/4.D0/T-WVN**2*T)/4.D0/PI/T

        ENDIF

       ENDIF


       RETURN
       END FUNCTION F1

! ==================================================================================
!       Purpose: This program computes derivative of the expansion coefficients
!                in Linton's method by integral form.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   R --- horizontal distance of the source point and the field point
!                WVN --- The input wave number
!                NO --- the term number to be calculated in Linton's expansion
!                T ---  integration variable
!       Output:  F2 --- derivative value of the expansion coefficient by integral
!                        representation
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION F2(R,WVN,NO,T)


       INTEGER NO
       REAL*8,INTENT(IN)::R,WVN,T
       REAL*8 PI

       DATA  PI/3.141592653589793D0/

       IF (DABS(T).LT.1.D-8) THEN

        F2=0.D0

       ELSE

        IF (NO.EQ.1) THEN

         F2=-R*DEXP(-R**2/4.D0/T+WVN**2*T)/8.D0/PI/T**2

        ELSE

         F2=-R*DEXP(-R**2/4.D0/T-WVN**2*T)/8.D0/PI/T**2

        ENDIF

       ENDIF

       RETURN
       END FUNCTION F2

! ==================================================================================
!       Purpose: This program computes the value of the survival integral
!                in L1 function in Linton's method.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   V --- wave number in deep water , i.e., V= w^2/g
!                R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                H --- water depth (h>0)
!                T ---  integration variable
!       Output:  G1 --- value of the survival integral in L1 in Linton's
!                       representation
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION G1(V,R,ZF,ZP,H,T)

       REAL*8,INTENT(IN)::V,R,ZF,ZP,H,T
       INTEGER I
       REAL*8 SUM,VAL(4)
       !REAL*8,EXTERNAL::ERFCC

       IF (DABS(T).LT.1.D-8) THEN

        G1=0.D0

       ELSE

        CALL CHI(ZF,ZP,H,VAL)

        SUM=0.D0
        DO I=1,4
         SUM=SUM+DEXP(-V*VAL(I))*ERFCC(VAL(I)/2.D0/T-V*T)
        ENDDO

        G1=DEXP(V**2*T**2-R**2/4.D0/T**2)*SUM/T

       ENDIF

        RETURN
        END FUNCTION G1

! ==================================================================================
!       Purpose: This program computes derivative of the survival integral
!                in L1 function with respect to R in Linton's method.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   V --- wave number in deep water , i.e., V= w^2/g
!                R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                H --- water depth (h>0)
!                T ---  integration variable
!       Output:  G2 --- derivative value of the survival integral in L1 with respect
!                       to R in Linton's representation
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION G2(V,R,ZF,ZP,H,T)


       REAL*8,INTENT(IN)::V,R,ZF,ZP,H,T
       INTEGER I
       REAL*8 SUM,VAL(4)
       !REAL*8,EXTERNAL::ERFCC

       IF (DABS(T).LT.1.D-8) THEN

        G2=0.D0

       ELSE

        CALL CHI(ZF,ZP,H,VAL)

        SUM=0.D0
        DO I=1,4
         SUM=SUM+DEXP(-V*VAL(I))*ERFCC(VAL(I)/2.D0/T-V*T)
        ENDDO

        G2=R/2.D0/T**3*DEXP(V**2*T**2-R**2/4.D0/T**2)*SUM

       ENDIF

       RETURN
       END FUNCTION G2

! ==================================================================================
!       Purpose: This program computes derivative of the survival integral
!                in L1 function with respect to z in Linton's method.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on Sept 7, 2013
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!  Parameters:
!       Input:   V --- wave number in deep water , i.e., V= w^2/g
!                R --- horizontal distance of the source point and the field point
!                ZF,ZP --- Z coordinates of the source point and the field point
!                H --- water depth (h>0)
!                T ---  integration variable
!       Output:  G3 --- derivative value of the survival integral in L1 with respect
!                       to z in Linton's representation
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       REAL*8 FUNCTION G3(V,R,ZF,ZP,H,T)

       REAL*8,INTENT(IN)::V,R,ZF,ZP,H,T
       INTEGER I
       REAL*8 SUM,VAL(4),DVAL(4)
       !REAL*8,EXTERNAL::ERFCC,DERFCC

       IF (DABS(T).LT.1.D-8) THEN

        G3=0.D0

       ELSE

        CALL CHI(ZF,ZP,H,VAL)
        CALL DCHI(DVAL)

        SUM=0.D0
        DO I=1,4
         SUM=SUM-V*DVAL(I)*DEXP(-V*VAL(I))*ERFCC(VAL(I)/2.D0/T-V*T)    &
             +DEXP(-V*VAL(I))*DERFCC(VAL(I)/2.D0/T-V*T)*DVAL(I)/2.D0/T
        ENDDO

        G3=DEXP(V**2*T**2-R**2/4.D0/T**2)*SUM/T

       ENDIF


       RETURN
       END FUNCTION G3
!
!    ===========================================================
!      Functions 1/dsqrt(R^2+chi(1,i)^2) and their derivatives
!         Code Original Author: Yingyi Liu   2013.09.07
!    ===========================================================

       SUBROUTINE RCHI(R,ZF,ZP,H,RO,DROR,DROZ)

       INTEGER I
       REAL*8 R,ZF,ZP,H
       REAL*8 RO(4),DROR(4),DROZ(4),VAL(4)


       CALL CHI(ZF,ZP,H,VAL)

       DO I=1,4

       RO(I)=DSQRT(R**2+VAL(I)**2)
       DROR(I)=R/DSQRT(R**2+VAL(I)**2)
       DROZ(I)=VAL(I)*(-1.D0)**I/DSQRT(R**2+VAL(I)**2)

       ENDDO

       RETURN
       END SUBROUTINE RCHI
!
!    ===========================================================
!                        Function Chi(1,i)
!         Code Original Author: Yingyi Liu   2013.09.07
!    ===========================================================

       SUBROUTINE CHI(ZF,ZP,H,VAL)

       REAL*8 ZF,ZP,H,VAL(4)

       VAL(1)=-ZP-ZF
       VAL(2)=2.D0*H-ZP+ZF
       VAL(3)=2.D0*H+ZP-ZF
       VAL(4)=4.D0*H+ZP+ZF

       RETURN
       END SUBROUTINE CHI
!
!    ===========================================================
!       The derivative of function Chi(1,i) with respect to z
!         Code Original Author: Yingyi Liu   2013.09.07
!    ===========================================================

       SUBROUTINE DCHI(DVAL)

       REAL*8 DVAL(4)

       DVAL(1)=-1.D0
       DVAL(2)=1.D0
       DVAL(3)=-1.D0
       DVAL(4)=1.D0

       RETURN
       END SUBROUTINE DCHI

!
!    ===========================================================
!                         Function Nm
!         Code Original Author: Yingyi Liu   2013.09.07
!    ===========================================================

       FUNCTION NM(H,WVN)

       REAL*8 H
       COMPLEX*16 NM,WVN

       NM=H/2.D0*(1.D0+CDSIN(2.D0*WVN*H)/(2.D0*WVN*H))

       RETURN
       END FUNCTION NM
!
!    ===========================================================
!          Derivative of the complementary error function
!         Code Original Author: Yingyi Liu   2013.09.07
!    ===========================================================

       FUNCTION DERFCC(X)

       REAL*8 DERFCC,X
       REAL*8 PI

       DATA  PI/3.141592653589793D0/

       DERFCC=-2.D0/DSQRT(PI)*DEXP(-X**2.D0)

       RETURN
       END FUNCTION DERFCC

!
!    ===========================================================
!             The complementary error function erfc(x)
!         Code Original Author: Yingyi Liu   2013.09.07
!    ===========================================================

       FUNCTION ERFCC(X)
       REAL*8 ERFCC,X,ERR

       CALL ERROR(X,ERR)
       ERFCC=1.D0-ERR

       RETURN
       END FUNCTION ERFCC



! ==================================================================================
!       Purpose: This program computes the 1D integral value of an external
!                integrand function in a finite interval (A,B) using standard
!                Kronrod extension of Gauss-Lengendre rule,  here using the
!                commonly used standard pair (G7,K15).
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on June 2, 2012
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!    Yingyi Liu, Ying Gou, Bin Teng, Shigeo Yoshida
!    An Extremely Efficient Boundary Element Method for Wave Interaction with Long
!    Cylindrical Structures Based on Free-Surface Green's Functions.
!    Computation, 4(2016), 36
!
!  Parameters:
!       Input:   V,R,ZF,ZP,H --- variables of the integrated function, should be
!                        changed by the user when he/she uses this program separately
!                        from this library for other purpose.
!                A,B --- the two end-points of the interval (A,B).
!                FUNTAG --- a flag for indicating which function components to be integrated.
!  Ext. Routine: FUN --- the external subroutine of the function integrand to be integrated,
!                        this should be specified by the user when he/she uses this program
!                        separately from this library for other purpose.
!       Output:  GK_INT --- the integral value calculated by Gauss-Kronrod rule,
!                           which is the final output and more accurate.
!                ERR --- the error estimate of this integration within interval (A,B).
!
!  Passed Variables:
!                XGQ(7),WGQ(7) --- the abscissas and weights for the Gauss rule.
!                XKQ(15),WKQ(15) --- the abscissas and weights for the Gauss-Kronrod rule.
!                G_INT --- the integral value calculated by Gauss rule.
!                STF(7) ---  the stored integrand function value calculated by Gauss rule,
!                            which can be reused by the Gauss-Kronrod rule.
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       SUBROUTINE GK_INTEGF(V,R,ZF,ZP,H,A,B,GK_INT,ERR,FUNTAG)

       INTEGER,INTENT(IN)::FUNTAG
       REAL*8,INTENT(IN)::V,R,ZF,ZP,H,A,B
       REAL*8,INTENT(OUT)::GK_INT,ERR
       !REAL*8,EXTERNAL::FUNF
       INTEGER I
       REAL*8 XGQ(7),WGQ(7),XKQ(15),WKQ(15),STF(7),G_INT


! ** matrix XGQ store the Gauss-Legendre sampling points(7)
!
        DATA XGQ/-0.949107912342759D0,-0.741531185599394D0,     &
                 -0.405845151377397D0, 0.000000000000000D0,     &
                  0.405845151377397D0, 0.741531185599394D0,     &
                  0.949107912342759D0/
!
! ** matrix WGQ store the Gauss-Legendre weighting factors
!
        DATA WGQ/ 0.129484966168870D0, 0.279705391489277D0,     &
                  0.381830050505119D0, 0.417959183673469D0,     &
                  0.381830050505119D0, 0.279705391489277D0,     &
                  0.129484966168870D0/
!
! ** matrix XKQ store the Kronrod sampling points(15)
!
        DATA XKQ/-0.991455371120813D0,-0.949107912342759D0,     &
                 -0.864864423359769D0,-0.741531185599394D0,     &
                 -0.586087235467691D0,-0.405845151377397D0,     &
                 -0.207784955007898D0, 0.000000000000000D0,     &
                  0.207784955007898D0, 0.405845151377397D0,     &
                  0.586087235467691D0, 0.741531185599394D0,     &
                  0.864864423359769D0, 0.949107912342759D0,     &
                  0.991455371120813D0/
!
! ** matrix WKQ store the weighting factors for the Kronrod
!    sampling points
!
        DATA WKQ/ 0.022935322010529D0, 0.063092092629979D0,     &
                  0.104790010322250D0, 0.140653259715525D0,     &
                  0.169004726639267D0, 0.190350578064785D0,     &
                  0.204432940075298D0, 0.209482141084728D0,     &
                  0.204432940075298D0, 0.190350578064785D0,     &
                  0.169004726639267D0, 0.140653259715525D0,     &
                  0.104790010322250D0, 0.063092092629979D0,     &
                  0.022935322010529D0/

        G_INT=0.0D0
        DO I=1,7
         STF(I)=FUNF(V,R,ZF,ZP,H,1/2.0D0*(B+A+(B-A)*XGQ(I)),FUNTAG)
         G_INT=G_INT+WGQ(I)*(B-A)/2.0D0*STF(I)
        END DO

        GK_INT=0.0D0
        DO I=1,15
         IF (MOD(I,2)==0) THEN
           GK_INT=GK_INT+WKQ(I)*(B-A)/2.0D0*STF(I/2)
         ELSE
         GK_INT=GK_INT+WKQ(I)*(B-A)/2.0D0*FUNF(V,R,ZF,ZP,H,1/2.0D0*(B+A+(B-A)*XKQ(I)),FUNTAG)
         ENDIF
       END DO
          ERR=(200.0D0*DABS(GK_INT-G_INT))**(1.5D0)

       RETURN
      END SUBROUTINE GK_INTEGF


! ==================================================================================
!       Purpose: This program using an adaptive quadrature method to compute the
!                integral value of an external integrand function in a finite
!                interval (A,B) by automatically dividing the interval into halves
!                and continuously calling the Gauss-Kronrod subroutine in order
!                to finally meet the requested accuracy Eps.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on June 2, 2012
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!    Yingyi Liu, Ying Gou, Bin Teng, Shigeo Yoshida
!    An Extremely Efficient Boundary Element Method for Wave Interaction with Long
!    Cylindrical Structures Based on Free-Surface Green's Functions.
!    Computation, 4(2016), 36
!
!  Parameters:
!       Input:   V,R,ZF,ZP,H --- variables of the integrated function, should be
!                        changed by the user when he/she uses this program separately
!                        from this library for other purpose.
!                A,B --- the two end-points of the interval (A,B).
!                EPS --- the requested tolerance for the integration.
!                FUNTAG --- a flag for indicating which function components to be integrated.
!  Ext. Routine: GK_INTEG --- the external subroutine of Gauss-Kronrod rule, this should be
!                        modified a little by the user when he/she uses this program
!                        separately from this library for other purpose.
!       Output:  ANS --- the result of integration.
!
!  Passed Variables:
!                GK_INT --- the integral value calculated by subroutine GK_INTEG.
!                ERR ---  the error estimate calculated by subroutine GK_INTEG.
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       RECURSIVE REAL*8 FUNCTION AQUADF(V,R,ZF,ZP,H,A,B,EPS,FUNTAG) RESULT(ANS)

       REAL*8,INTENT(IN)::V,R,ZF,ZP,H,A,B,EPS
       INTEGER,INTENT(IN)::FUNTAG
       REAL*8 GK_INT,ERR

       ANS=0.0D0
       CALL GK_INTEGF(V,R,ZF,ZP,H,A,B,GK_INT,ERR,FUNTAG)

       IF (ERR>EPS.AND.DABS(A-B)>EPS) THEN
         ANS=ANS+AQUADF(V,R,ZF,ZP,H,A,(A+B)/2.0D0,EPS,FUNTAG)+AQUADF(V,R,ZF,ZP,H,(A+B)/2.0D0,B,EPS,FUNTAG)
       ELSE
          ANS=ANS+GK_INT
       ENDIF

       RETURN
       END FUNCTION AQUADF



! ==================================================================================
!       Purpose: This program computes the 1D integral value of an external
!                integrand function in a finite interval (A,B) using standard
!                Kronrod extension of Gauss-Lengendre rule,  here using the
!                commonly used standard pair (G7,K15).
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on June 2, 2012
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!    Yingyi Liu, Ying Gou, Bin Teng, Shigeo Yoshida
!    An Extremely Efficient Boundary Element Method for Wave Interaction with Long
!    Cylindrical Structures Based on Free-Surface Green's Functions.
!    Computation, 4(2016), 36
!
!  Parameters:
!       Input:   R,WVN,NO --- variables of the integrated function, should be
!                        changed by the user when he/she uses this program separately
!                        from this library for other purpose.
!                A,B --- the two end-points of the interval (A,B).
!                FUNTAG --- a flag for indicating which function components to be integrated.
!  Ext. Routine: FUN --- the external subroutine of the function integrand to be integrated,
!                        this should be specified by the user when he/she uses this program
!                        separately from this library for other purpose.
!       Output:  GK_INT --- the integral value calculated by Gauss-Kronrod rule,
!                           which is the final output and more accurate.
!                ERR --- the error estimate of this integration within interval (A,B).
!
!  Passed Variables:
!                XGQ(7),WGQ(7) --- the abscissas and weights for the Gauss rule.
!                XKQ(15),WKQ(15) --- the abscissas and weights for the Gauss-Kronrod rule.
!                G_INT --- the integral value calculated by Gauss rule.
!                STF(7) ---  the stored integrand function value calculated by Gauss rule,
!                            which can be reused by the Gauss-Kronrod rule.
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       SUBROUTINE GK_INTEG(R,WVN,NO,A,B,GK_INT,ERR,FUNTAG)

        INTEGER,INTENT(IN)::NO,FUNTAG
        REAL*8,INTENT(IN)::R,WVN,A,B
        REAL*8,INTENT(OUT)::GK_INT,ERR
        !REAL*8,EXTERNAL::FUN
        INTEGER I
        REAL*8 XGQ(7),WGQ(7),XKQ(15),WKQ(15),STF(7),G_INT


! ** matrix XGQ store the Gauss-Legendre sampling points(7)
!
        DATA XGQ/-0.949107912342759D0,-0.741531185599394D0,     &
                 -0.405845151377397D0, 0.000000000000000D0,     &
                  0.405845151377397D0, 0.741531185599394D0,     &
                  0.949107912342759D0/
!
! ** matrix WGQ store the Gauss-Legendre weighting factors
!
        DATA WGQ/ 0.129484966168870D0, 0.279705391489277D0,     &
                  0.381830050505119D0, 0.417959183673469D0,     &
                  0.381830050505119D0, 0.279705391489277D0,     &
                  0.129484966168870D0/
!
! ** matrix XKQ store the Kronrod sampling points(15)
!
        DATA XKQ/-0.991455371120813D0,-0.949107912342759D0,     &
                 -0.864864423359769D0,-0.741531185599394D0,     &
                 -0.586087235467691D0,-0.405845151377397D0,     &
                 -0.207784955007898D0, 0.000000000000000D0,     &
                  0.207784955007898D0, 0.405845151377397D0,     &
                  0.586087235467691D0, 0.741531185599394D0,     &
                  0.864864423359769D0, 0.949107912342759D0,     &
                  0.991455371120813D0/
!
! ** matrix WKQ store the weighting factors for the Kronrod
!    sampling points
!
        DATA WKQ/ 0.022935322010529D0, 0.063092092629979D0,     &
                  0.104790010322250D0, 0.140653259715525D0,     &
                  0.169004726639267D0, 0.190350578064785D0,     &
                  0.204432940075298D0, 0.209482141084728D0,     &
                  0.204432940075298D0, 0.190350578064785D0,     &
                  0.169004726639267D0, 0.140653259715525D0,     &
                  0.104790010322250D0, 0.063092092629979D0,     &
                  0.022935322010529D0/

        G_INT=0.0D0
        DO I=1,7
         STF(I)=FUN(R,WVN,NO,1/2.0D0*(B+A+(B-A)*XGQ(I)),FUNTAG)
         G_INT=G_INT+WGQ(I)*(B-A)/2.0D0*STF(I)
        END DO

        GK_INT=0.0D0
        DO I=1,15
         IF (MOD(I,2)==0) THEN
           GK_INT=GK_INT+WKQ(I)*(B-A)/2.0D0*STF(I/2)
         ELSE
         GK_INT=GK_INT+WKQ(I)*(B-A)/2.0D0*FUN(R,WVN,NO,1/2.0D0*(B+A+(B-A)*XKQ(I)),FUNTAG)
         ENDIF
       END DO

       ERR=(200.0D0*DABS(GK_INT-G_INT))**(1.5D0)

       RETURN
       END SUBROUTINE GK_INTEG


! ==================================================================================
!       Purpose: This program using an adaptive quadrature method to compute the
!                integral value of an external integrand function in a finite
!                interval (A,B) by automatically dividing the interval into halves
!                and continuously calling the Gauss-Kronrod subroutine in order
!                to finally meet the requested accuracy Eps.
!
!  License:
!
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Code Original Author:
!
!    Yingyi Liu on June 2, 2012
!
!  Modified on:
!
!    July 21, 2016
!
!  Reference:
!
!    Yingyi Liu et al. A reliable open-source package for performance evaluation of
!    floating renewable energy systems in coastal and offshore regions. Energy Conversion
!    and management, 2018
!
!    Yingyi Liu, Ying Gou, Bin Teng, Shigeo Yoshida
!    An Extremely Efficient Boundary Element Method for Wave Interaction with Long
!    Cylindrical Structures Based on Free-Surface Green's Functions.
!    Computation, 4(2016), 36
!
!  Parameters:
!       Input:   R,WVN,NO --- variables of the integrated function, should be
!                        changed by the user when he/she uses this program separately
!                        from this library for other purpose.
!                A,B --- the two end-points of the interval (A,B).
!                EPS --- the requested tolerance for the integration.
!                FUNTAG --- a flag for indicating which function components to be integrated.
!  Ext. Routine: GK_INTEG --- the external subroutine of Gauss-Kronrod rule, this should be
!                        modified a little by the user when he/she uses this program
!                        separately from this library for other purpose.
!       Output:  ANS --- the result of integration.
!
!  Passed Variables:
!                GK_INT --- the integral value calculated by subroutine GK_INTEG.
!                ERR ---  the error estimate calculated by subroutine GK_INTEG.
!
!  Contributors list:
!         Yingyi Liu
!         to be continued...
!
! ==================================================================================

       RECURSIVE REAL*8 FUNCTION AQUAD(R,WVN,NO,A,B,EPS,FUNTAG) RESULT(ANS)

        REAL*8,INTENT(IN)::R,WVN,A,B,EPS
        INTEGER,INTENT(IN)::NO,FUNTAG
        REAL*8 GK_INT,ERR

       ANS=0.0D0
       CALL GK_INTEG(R,WVN,NO,A,B,GK_INT,ERR,FUNTAG)

       IF (ERR>EPS.AND.DABS(A-B)>EPS) THEN
         ANS=ANS+AQUAD(R,WVN,NO,A,(A+B)/2.0D0,EPS,FUNTAG)+AQUAD(R,WVN,NO,(A+B)/2.0D0,B,EPS,FUNTAG)
       ELSE
          ANS=ANS+GK_INT
       ENDIF

       RETURN
       END FUNCTION AQUAD

  end module fingreen3D_module
