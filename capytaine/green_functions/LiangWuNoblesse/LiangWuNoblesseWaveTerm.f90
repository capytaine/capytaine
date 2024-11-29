! 
!==================================================================================
! 
!     Purpose: This program evaluates the three-dimensional free-surface
!              Green function and derivatives in deep water
! 
!              Code Original Author: Hui Liang       created on  2016.12.26    
!              Superficially adapted for Capytaine by Matthieu Ancellin, 2022--2024 
!
!  License:
!
!
!    This routine is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option) 
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3  
!    or later), along with this routine. If not, see <http://www.gnu.org/licenses/>.
!
!  Modified on:
! 
!    January 02, 2018
! 
!  Reference:
! 
!
!  [1]  H. Wu, C. Zhang, Y. Zhu, W. Li, D. Wan, F. Noblesse, 
!    A global approximation to the Green function for 
!    diffraction radiation of water waves, 
!    Eur. J. Mech. B Fluids 65 (2017) 54-64.
!
!  [2]  H. Liang, H. Wu, F. Noblesse,
!    Validation of a global approximation for 
!    wave diffraction-radiation in deep water,
!    Appl. Ocean Res. 74 (2018) 80-86.
! 
!   Remarks
! 
!   The local-flow component is approximated by mean of the global 
!   approximations [1]. The computations reported in [2] provides
!   strong evidence that the global approximations are sufficiently 
!   accurate to compute linear and second-order wave loads in practice.
!  
!
!   It should be noted that the Rankine source term -2/d appeared in L_z 
!   given by (8a) in Ref. [2] is not evaluated here (there is a typo in (8a)).  
!   These Rankine source components are required to evaluate in another routine  
!   due to strong singular behavior.
!
!   For any questions, please contact: lianghuistar@gmail.com
!                   
!   We define the flow-field point (x,y,z) and source point (xi,eta,zeta).               
!   Please note that all variables are non-dimensionalized with respect to 
!   the wavenumber k0
!                   
!   The Green function is defiend as                
!                   
!   G = -1/r-1/d+GF,
!
!   where   
!                                          
!   r = sqrt((x-xi)^2+(y-eta)^2+(z-zeta)^2);                  
!   d = sqrt((x-xi)^2+(y-eta)^2+(z+zeta)^2).
!
!   Parameters:
!
!       Input:      hh  --- sqrt((x-xi)^2 + (y-eta)^2)
!                   vv  --- z+zeta
! 
!       Output:     GF  --- free surface itself GF
!                   GFh  --- derivative of GF with respect to hh
!
! ==================================================================================

module LiangWuNoblesseWaveTerm

implicit none

real(8), parameter    ::  pi   =  4.0d0*datan(1.0d0)
real(8), parameter    ::  gama =  0.5772156649d0
complex(8), parameter ::  Im   =  dcmplx(0.0d0,1.0d0)

contains

!===============================================================
subroutine HavelockGF(hh,vv,GF,GFh)
  implicit none
! --- Variables -------------------------------------------
  real(8),intent(in)    ::  hh,vv
  complex(8),intent(out)  ::  GF, GFh

! --- Local variables -------------------------------------
  real(8)                ::  dd, alpha, beta, rho

  dd = sqrt(hh*hh+vv*vv)
  alpha = -vv/dd
  beta = hh/dd
  rho = dd/(1.0d0+dd)
 
  GF = GF_Func_L0(hh, vv, dd, alpha, beta, rho) + GF_Func_W(hh, vv)
  GFh = GF_Func_Ls(hh, vv, dd, alpha, beta, rho) + GF_Func_Wh(hh, vv)

end subroutine HavelockGF

!=============================================================

function GF_Func_L0(hh,vv, dd, alpha, beta, rho)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  hh,vv, dd, alpha, beta, rho

! --- Local variables -------------------------------------
  real(8)          ::  PP,Lp
  real(8)          ::  GF_Func_L0

  PP  =  dlog(0.5d0*(dd-vv))+gama-2.0d0*dd*dd
  PP  =  dexp(vv)*PP
  PP  =  PP+dd*dd-vv
    
  Lp  =  GF_Func_Lp(hh,vv, dd, alpha, beta, rho)
    
  GF_Func_L0  =  2.0d0*PP/(1.0d0+dd**3)+2.0d0*Lp

  return
end function GF_Func_L0

!=============================================================
function GF_Func_Lp(hh,vv, dd, alpha, beta, rho)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  hh,vv, dd, alpha, beta, rho

! --- Local variables -------------------------------------
  real(8)          ::  A,B,C,D,RR
  real(8)          ::  GF_Func_Lp

  A  =  GF_FuncA(rho)
  B  =  GF_FuncB(rho)
  C  =  GF_FuncC(rho)
  D  =  GF_FuncD(rho)

  RR  =  (1.0d0-beta)*A
  RR  =  RR-beta*B
  RR  =  RR-alpha*C/(1.0d0+6.0d0*alpha*rho*(1.0d0-rho))
  RR  =  RR+beta*(1.0d0-beta)*D

  GF_Func_Lp  =  rho*(1.0d0-rho)**3*RR

  return
end function GF_Func_Lp

!=============================================================
function GF_Func_W(hh,vv)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  hh,vv

! --- Local variables -------------------------------------
  real(8)          ::  H0,J0
  complex(8)        ::  GF_Func_W

  H0  =  StruveH0(hh)
  J0  =  BesselJ0(hh)

  GF_Func_W  =  2.0d0*pi*(H0-Im*J0)*dexp(vv)

  return
end function GF_Func_W

!=============================================================
function GF_FuncA(tt)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  tt

! --- Local variables -------------------------------------
  real(8)          ::  A0,A1,A2,A3,A4
  real(8)          ::  A5,A6,A7,A8,A9
  real(8)          ::  GF_FuncA

  data  A0,A1,A2,A3,A4,A5,A6,A7,A8,A9      &
    /  +1.21d0,  -13.328d0,  +215.896d0,    &
      -1763.96d0,  +8418.94d0,  -24314.21d0,  &
      +42002.57d0,-41592.9d0,  21859.0d0,    &
      -4838.6d0  /

  GF_FuncA  =  (A8+A9*tt)*tt
  GF_FuncA  =  (A7+GF_FuncA)*tt
  GF_FuncA  =  (A6+GF_FuncA)*tt
  GF_FuncA  =  (A5+GF_FuncA)*tt
  GF_FuncA  =  (A4+GF_FuncA)*tt
  GF_FuncA  =  (A3+GF_FuncA)*tt
  GF_FuncA  =  (A2+GF_FuncA)*tt
  GF_FuncA  =  (A1+GF_FuncA)*tt
  GF_FuncA  =  A0+GF_FuncA

  return
end function GF_FuncA

!=============================================================
function GF_FuncB(tt)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  tt

! --- Local variables -------------------------------------
  real(8)          ::  B0,B1,B2,B3,B4
  real(8)          ::  B5,B6,B7,B8,B9
  real(8)          ::  GF_FuncB

  data  B0,B1,B2,B3,B4,B5,B6,B7,B8,B9      &
    /  +0.938d0,  +5.737d0,  -67.92d0,    &
      +796.534d0,  -4780.77d0,  +17137.74d0,  &
      -36618.81d0,+44894.06d0,-29030.24d0,  &
      +7671.22d0  /

  GF_FuncB  =  (B8+B9*tt)*tt
  GF_FuncB  =  (B7+GF_FuncB)*tt
  GF_FuncB  =  (B6+GF_FuncB)*tt
  GF_FuncB  =  (B5+GF_FuncB)*tt
  GF_FuncB  =  (B4+GF_FuncB)*tt
  GF_FuncB  =  (B3+GF_FuncB)*tt
  GF_FuncB  =  (B2+GF_FuncB)*tt
  GF_FuncB  =  (B1+GF_FuncB)*tt
  GF_FuncB  =  B0+GF_FuncB

  return
end function GF_FuncB

!=============================================================
function GF_FuncC(tt)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  tt

! --- Local variables -------------------------------------
  real(8)          ::  C0,C1,C2,C3,C4
  real(8)          ::  C5,C6,C7
  real(8)          ::  GF_FuncC

  data  C0,C1,C2,C3,C4,C5,C6,C7          &
    /  +1.268d0,  -9.747d0,  +209.653d0,    &
      -1397.89d0,  +5155.67d0,  -9844.35d0,    &
      +9136.4d0,  -3272.62d0    /

  GF_FuncC  =  (C6+C7*tt)*tt
  GF_FuncC  =  (C5+GF_FuncC)*tt
  GF_FuncC  =  (C4+GF_FuncC)*tt
  GF_FuncC  =  (C3+GF_FuncC)*tt
  GF_FuncC  =  (C2+GF_FuncC)*tt
  GF_FuncC  =  (C1+GF_FuncC)*tt
  GF_FuncC  =  C0+GF_FuncC

  return
end function GF_FuncC

!=============================================================
function GF_FuncD(tt)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  tt

! --- Local variables -------------------------------------
  real(8)          ::  D0,D1,D2,D3,D4
  real(8)          ::  D5,D6,D7,D8,D9
  real(8)          ::  GF_FuncD
      
  data  D0,D1,D2,D3,D4,D5,D6,D7,D8,D9      &
    /  +0.632d0,  -40.97d0,  +667.16d0,    &
      -6072.07d0,  +31127.39d0,-96293.05d0,  &
      +181856.75d0,      -205690.43d0,  &
      +128170.2d0,-33744.6d0  /
    
  GF_FuncD  =  (D8+D9*tt)*tt
  GF_FuncD  =  (D7+GF_FuncD)*tt
  GF_FuncD  =  (D6+GF_FuncD)*tt
  GF_FuncD  =  (D5+GF_FuncD)*tt
  GF_FuncD  =  (D4+GF_FuncD)*tt
  GF_FuncD  =  (D3+GF_FuncD)*tt
  GF_FuncD  =  (D2+GF_FuncD)*tt
  GF_FuncD  =  (D1+GF_FuncD)*tt
  GF_FuncD  =  D0+GF_FuncD

  return
end function GF_FuncD
    
!=============================================================
function GF_Func_Ls(hh,vv, dd, alpha, beta, rho)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  hh,vv, dd, alpha, beta, rho

! --- Local variables -------------------------------------
  real(8)          ::  PS,QS,Lsp
  real(8)          ::  GF_Func_Ls

  PS  =  (beta+hh)/(dd-vv)
  PS  =  PS-2.0d0*beta+2.0d0*dexp(vv)*dd-hh
 
  QS  =  dexp(-dd)*(1.0d0-beta)
  QS  =  QS*(1.0d0+dd/(1.0d0+dd**3))
 
  Lsp  =  GF_Func_Lsp(hh,vv, dd, alpha, beta, rho)
 
  GF_Func_Ls  =  2.0d0*PS/(1.0d0+dd**3)-4.0d0*QS+2.0d0*Lsp

  return
end function GF_Func_Ls

!=============================================================
function GF_Func_Lsp(hh,vv, dd, alpha, beta, rho)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  hh,vv, dd, alpha, beta, rho

! --- Local variables -------------------------------------
  real(8)          ::  A,B,C,RR
  real(8)          ::  GF_Func_Lsp

  A  =  GF_dFuncA(rho)
  B  =  GF_dFuncB(rho)
  C  =  GF_dFuncC(rho)

  RR  =  beta*A
  RR  =  RR-(1.0d0-alpha)*B
  RR  =  RR+beta*(1.0d0-beta)*rho*(1.0d0-2.0d0*rho)*C

  GF_Func_Lsp  =  rho*(1.0d0-rho)**3*RR

  return
end function GF_Func_Lsp

!=============================================================
function GF_Func_Wh(hh,vv)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  hh,vv

! --- Local variables -------------------------------------
  real(8)          ::  H1,J1
  complex(8)        ::  GF_Func_Wh

  H1  =  StruveH1(hh)
  J1  =  BesselJ1(hh)
    
  GF_Func_Wh  =  2.0d0*pi*(2.0d0/pi-H1+Im*J1)*dexp(vv)

  return
end function GF_Func_Wh

!=============================================================
function GF_dFuncA(tt)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  tt

! --- Local variables -------------------------------------
  real(8)          ::  A0,A1,A2,A3,A4
  real(8)          ::  A5,A6,A7,A8,A9
  real(8)          ::  GF_dFuncA
      
  data  A0,A1,A2,A3,A4,A5,A6,A7,A8,A9      &
    /  +2.948d0,  -24.53d0,  +249.69d0,    &
      -754.85d0,  -1187.71d0,  +16370.75d0,  &
      -48811.41d0,+68220.87d0,-46688.0d0,    &
      +12622.25d0  /

  GF_dFuncA  =  (A8+A9*tt)*tt
  GF_dFuncA  =  (A7+GF_dFuncA)*tt
  GF_dFuncA  =  (A6+GF_dFuncA)*tt
  GF_dFuncA  =  (A5+GF_dFuncA)*tt
  GF_dFuncA  =  (A4+GF_dFuncA)*tt
  GF_dFuncA  =  (A3+GF_dFuncA)*tt
  GF_dFuncA  =  (A2+GF_dFuncA)*tt
  GF_dFuncA  =  (A1+GF_dFuncA)*tt
  GF_dFuncA  =  A0+GF_dFuncA

  return
end function GF_dFuncA

!=============================================================
function GF_dFuncB(tt)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  tt

! --- Local variables -------------------------------------
  real(8)          ::  B0,B1,B2,B3,B4
  real(8)          ::  B5,B6,B7,B8,B9
  real(8)          ::  GF_dFuncB
      
  data  B0,B1,B2,B3,B4,B5,B6,B7,B8,B9      &
    /  +1.11d0,  +2.894d0,  -76.765d0,    &
      +1565.35d0,  -11336.19d0,+44270.15d0,  &
      -97014.11d0,+118879.26d0,-76209.82d0,  &
      +19923.28d0  /
    
  GF_dFuncB  =  (B8+B9*tt)*tt
  GF_dFuncB  =  (B7+GF_dFuncB)*tt
  GF_dFuncB  =  (B6+GF_dFuncB)*tt
  GF_dFuncB  =  (B5+GF_dFuncB)*tt
  GF_dFuncB  =  (B4+GF_dFuncB)*tt
  GF_dFuncB  =  (B3+GF_dFuncB)*tt
  GF_dFuncB  =  (B2+GF_dFuncB)*tt
  GF_dFuncB  =  (B1+GF_dFuncB)*tt
  GF_dFuncB  =  B0+GF_dFuncB

  return
end function GF_dFuncB

!=============================================================
function GF_dFuncC(tt)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  tt

! --- Local variables -------------------------------------
  real(8)          ::  C0,C1,C2,C3,C4,C5
  real(8)          ::  GF_dFuncC
      
  data  C0,C1,C2,C3,C4,C5            &
    /  +14.19d0,  -148.24d0,  +847.8d0,    &
      -2318.58d0,  +3168.35d0,  -1590.27d0    /
    
  GF_dFuncC  =  (C4+C5*tt)*tt
  GF_dFuncC  =  (C3+GF_dFuncC)*tt
  GF_dFuncC  =  (C2+GF_dFuncC)*tt
  GF_dFuncC  =  (C1+GF_dFuncC)*tt
  GF_dFuncC  =  C0+GF_dFuncC

  return
end function GF_dFuncC

!**************************************
function BesselJ0(xx)
!
  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  xx

! --- Local variables -------------------------------------
  real(8)          ::  yy,y2,f0,theta0
  real(8)          ::  BesselJ0        
  real(8)          ::  P0,P1,P2,P3,P4,P5,P6
  real(8)          ::  R0,R1,R2,R3,R4,R5
  real(8)          ::  S1,S2,S3,S4,S5
    
  data  P0,P1,P2,P3,P4,P5,P6    &
    /  +0.999999999d0,  -2.249999879d0,  +1.265623060d0,  &
      -0.316394552d0,  +0.044460948d0,  -0.003954479d0,  &
      +0.000212950d0  /
    
  data  R0,R1,R2,R3,R4,R5      &
    /  +0.79788454d0,  -0.00553897d0,  +0.00099336d0,  &
      -0.00044346d0,  +0.00020445d0,  -0.00004959d0  /
    
  data  S1,S2,S3,S4,S5        &
    /  -0.04166592d0,  +0.00239399d0,  -0.00073984d0,  &
      +0.00031099d0,  -0.00007605d0  /


  if(xx <= 3.0d0) then
    yy  =  (xx/3.0d0)**2
    BesselJ0  =   P0+(P1+(P2+(P3+(P4+(P5+P6*yy)*yy)*yy)*yy)*yy)*yy
  else
    yy  =  3.0d0/xx
    y2  =  yy**2
        
    f0      =  R0+(R1+(R2+(R3+(R4+R5*y2)*y2)*y2)*y2)*y2

    theta0    =   xx  -  0.25d0*pi    &
            +(S1+(S2+(S3+(S4+S5*y2)*y2)*y2)*y2)*yy

    BesselJ0  =  f0*dcos(theta0)/dsqrt(xx)
  end if

  return
end function BesselJ0

!===============================================================
function BesselJ1(xx)
!
  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  xx

! --- Local variables -------------------------------------
  real(8)          ::  yy,y2,f1,theta1
  real(8)          ::  BesselJ1
    
  real(8)          ::  P0,P1,P2,P3,P4,P5,P6
  real(8)          ::  R0,R1,R2,R3,R4,R5
  real(8)          ::  S1,S2,S3,S4,S5
    
  data  P0,P1,P2,P3,P4,P5,P6    &
    /  +0.500000000d0,  -0.562499992d0,  +0.210937377d0,  &
      -0.039550040d0,  +0.004447331d0,  -0.000330547d0,  &
      +0.000015525d0  /
    
  data  R0,R1,R2,R3,R4,R5      &
    /  +0.79788459d0,  +0.01662008d0,  -0.00187002d0,  &
      +0.00068519d0,  -0.00029440d0,  +0.00006952d0  /
    
  data  S1,S2,S3,S4,S5        &
    /  +0.12499895d0,  -0.00605240d0,  +0.00135825d0,  &
      -0.00049616d0,  +0.00011531d0  /

  if(xx <= 3.0d0) then
    yy  =  xx/3.0d0
    y2  =  yy*yy
    BesselJ1  =  P0+(P1+(P2+(P3+(P4+(P5+P6*y2)*y2)*y2)*y2)*y2)*y2

    BesselJ1  =  BesselJ1*xx
  else
    yy  =  3.0d0/xx
    y2  =  yy*yy
    f1      =  R0+(R1+(R2+(R3+(R4+R5*y2)*y2)*y2)*y2)*y2

    theta1    =   xx  -  0.75d0*pi    &
            +(S1+(S2+(S3+(S4+S5*y2)*y2)*y2)*y2)*yy

    BesselJ1  =  f1*dcos(theta1)/dsqrt(xx)
  end if

  return
end function BesselJ1

!===============================================================
function BesselY0(xx)
!
  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  xx

! --- Local variables -------------------------------------
  real(8)          ::  f0,theta0
  real(8)          ::  BesselY0        

  if(xx <= 3.0d0) then
    BesselY0  =  (2.0d0/pi)*dlog(xx/2.0d0)*BesselJ0(xx)  &  
            +0.367466907d0              &
            +0.605593797d0*(xx/3.0d0)**2      &                            
            -0.743505078d0*(xx/3.0d0)**4      &                            
            +0.253005481d0*(xx/3.0d0)**6      &                            
            -0.042619616d0*(xx/3.0d0)**8      &                            
            +0.004285691d0*(xx/3.0d0)**10      &                            
            -0.000250716d0*(xx/3.0d0)**12                            
  else
    f0      =   0.79788454d0          &
            -0.00553897d0*(3.0d0/xx)**2    &                            
            +0.00099336d0*(3.0d0/xx)**4    &                            
            -0.00044346d0*(3.0d0/xx)**6    &                            
            +0.00020445d0*(3.0d0/xx)**8    &                            
            -0.00004959d0*(3.0d0/xx)**10                            

    theta0    =   xx    -  0.25d0*pi      &
            -0.04166592d0*(3.0d0/xx)**1    &                            
            +0.00239399d0*(3.0d0/xx)**3    &                            
            -0.00073984d0*(3.0d0/xx)**5    &                            
            +0.00031099d0*(3.0d0/xx)**7    &                            
            -0.00007605d0*(3.0d0/xx)**9                            

    BesselY0  =  f0*dsin(theta0)/dsqrt(xx)                                    

  end if

  return
end function BesselY0

!===============================================================
  function BesselY1(xx)
!
  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  xx

! --- Local variables -------------------------------------
  real(8)          ::  f1,theta1
  real(8)          ::  BesselY1        

  if(xx <= 3.0d0) then
    BesselY1  =  (2.0d0/pi)*(dlog(xx/2.0d0)*BesselJ1(xx)-1.0d0/xx)  &  
            +0.07373571d0*(xx/3.0d0)**1        &                            
            +0.72276433d0*(xx/3.0d0)**3        &                            
            -0.43885620d0*(xx/3.0d0)**5        &                            
            +0.10418264d0*(xx/3.0d0)**7        &                            
            -0.01340825d0*(xx/3.0d0)**9        &                            
            +0.00094249d0*(xx/3.0d0)**11                            
  else
    f1      =   0.79788459d0            &
            +0.01662008d0*(3.0d0/xx)**2    &                            
            -0.00187002d0*(3.0d0/xx)**4    &                            
            +0.00068519d0*(3.0d0/xx)**6    &                            
            -0.00029440d0*(3.0d0/xx)**8    &                            
            +0.00006952d0*(3.0d0/xx)**10                            

    theta1    =   xx    -  3.0d0*pi/4.0d0    &
            +0.12499895d0*(3.0d0/xx)**1    &                            
            -0.00605240d0*(3.0d0/xx)**3    &                            
            +0.00135825d0*(3.0d0/xx)**5    &                            
            -0.00049616d0*(3.0d0/xx)**7    &                            
            +0.00011531d0*(3.0d0/xx)**9                            

    BesselY1  =  f1*dsin(theta1)/dsqrt(xx)                                    

  end if

  return
end function BesselY1


!===============================================================
function StruveH0(xx)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  xx

! --- Local variables -------------------------------------

  real(8)          ::  P0,P1,P2
  real(8)          ::  P3,P4,P5
  real(8)          ::  a0,a1,a2,a3
  real(8)          ::  b1,b2,b3
  real(8)          ::  c1,c2
  real(8)          ::  yy,StruveH0

  if(xx <= 3.0d0) then
        
    yy  =  (xx/3.0d0)**2
        
    P0  =  +1.909859164d0
    P1  =  -1.909855001d0
    P2  =  +0.687514637d0
    P3  =  -0.126164557d0
     P4  =  +0.013828813d0
    P5  =  -0.000876918d0
        
    StruveH0  =  P0+(P1+(P2+(P3+(P4+P5*yy)*yy)*yy)*yy)*yy
        
    StruveH0  =  StruveH0*(xx/3.0d0)
        
    else
        
    yy  =  (3.0d0/xx)**2

    a0  =  0.99999906d0
    a1  =  4.77228920d0
    a2  =  3.85542044d0
    a3  =  0.32303607d0

    b1  =  4.88331068d0
    b2  =  4.28957333d0
    b3  =  0.52120508d0

    c1  =  2.0d0*(a0  +  (a1+(a2+a3*yy)*yy)*yy) 
    c2  =  pi*xx*(1.0d0+  (b1+(b2+b3*yy)*yy)*yy) 
                                
    StruveH0  =  c1/c2+  BesselY0(xx)                                
                                        
  end if

  return
end function StruveH0

!===============================================================
function StruveH1(xx)

  implicit none

! --- Variables -------------------------------------------
  real(8),intent(in)    ::  xx

! --- Local variables -------------------------------------

  real(8)          ::  P1,P2,P3
  real(8)          ::  P4,P5,P6
  real(8)          ::  a0,a1,a2,a3
  real(8)          ::  b1,b2,b3
  real(8)          ::  c1,c2,yy
  real(8)          ::  StruveH1        

  if(xx <= 3.0d0) then
        
    yy  =  (xx/3.0d0)**2
        
    P1  =  +1.909859286d0
    P2  =  -1.145914713d0
    P3  =  +0.294656958d0
    P4  =  -0.042070508d0
     P5  =  +0.003785727d0
    P6  =  -0.000207183d0
      
    StruveH1  =  (P1+(P2+(P3+(P4+(P5+P6*yy)*yy)*yy)*yy)*yy)*yy

    else
        
    yy  =  (3.0d0/xx)**2

    a0  =  1.00000004d0
    a1  =  3.92205313d0
    a2  =  2.64893033d0
    a3  =  0.27450895d0

    b1  =  3.81095112d0
    b2  =  2.26216956d0
    b3  =  0.10885141d0

    c1  =  2.0d0*(a0  +  (a1+(a2+a3*yy)*yy)*yy) 
    c2  =  pi*(1.0d0  +  (b1+(b2+b3*yy)*yy)*yy) 

    StruveH1  =  c1/c2  +  BesselY1(xx)

  end if

  return
end function StruveH1

end module LiangWuNoblesseWaveTerm
