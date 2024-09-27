!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   auto 12 dof aircraft model in radians
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
SUBROUTINE FUNC(NDIM, U, ICP, PAR, IJAC, F, DFDU, DFDP)
    ! USE ThrustModel
    USE SimulationParams
    ! USE CalcCoeffExpRadModule
    ! USE CalcCoeffConstRadModule
    USE CalcCoeffComb
    IMPLICIT NONE

    
    INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

    DOUBLE PRECISION :: V, alpha, beta, p, q, r, phi, theta
    DOUBLE PRECISION, DIMENSION(3) :: CFor, CMom
    DOUBLE PRECISION :: del_e, del_a, del_r, eta
    DOUBLE PRECISION :: CD, CY, CLIFT, Cl, Cm, Cn

    ! Intermediate calculations
    DOUBLE PRECISION :: a11, a12, a13, a21, a22, a23, a24
    DOUBLE PRECISION :: a31, a32, a33, a34, a41, a42, a43, a44
    DOUBLE PRECISION :: a51, a52, a53, a61, a62, a63, a64

    CALL InitializeParameters()
    ! Tm = Thrust_Model(max_thrust)
    ! print *, 'Tm:',Tm
    ! print *, 'rho', rho
    CALL Calc_CoeffComb_exp(U, PAR, CFor, CMom)
    ! CALL Calc_Coeff_Const(U, PAR, CFor, CMom)
      
    ! Force coefficients
    CD = CFor(1)
    CY = CFor(2)
    CLIFT = CFor(3)

    ! Moment coefficients
    Cl = CMom(1)
    Cm = CMom(2)
    Cn = CMom(3)


    ! States
    V = U(1)
    alpha = U(2)
    beta = U(3)
    p = U(4)
    q = U(5)
    r = U(6)
    phi = U(7)
    theta = U(8)

    ! Control Parameters
    del_e = PAR(1)
    del_a = PAR(2)
    del_r = PAR(3)
    eta = PAR(4)

    ! Force equations in body frame
    a11 = (1.0/m) * eta * Tm * COS(alpha) * COS(beta)
    a12 = -(1.0/m) * 0.5 * rho * (v**2) * S * CD
    a13 = -g * (SIN(theta) * COS(alpha) * COS(beta) - COS(theta) * SIN(phi) * SIN(beta) - SIN(alpha) * & 
      & COS(beta) * COS(phi) * COS(theta))
    F(1)  = a11 + a12 + a13

    a21 = q - (p * COS(alpha) + r * SIN(alpha)) * TAN(beta)
    a22 = -(1.0/m) * (1.0/v) * (1.0/COS(beta)) * eta * Tm * SIN(alpha)
    a23 = -(1.0/m) * (1.0/COS(beta)) * 0.5 * rho * v * S * CLIFT
    a24 = (1.0/v) * (1.0/COS(beta)) * g * (SIN(alpha) * SIN(theta) + COS(alpha) * COS(phi) * COS(theta))
    F(2)  = a21 + a22 + a23 + a24

    a31 = p * sin(alpha) - r * cos(alpha)
    a32 = -(1 / m) * (1 / V) * eta * Tm * cos(alpha) * sin(beta)
    a33 = -(1 / m) * 0.5 * rho * V * S * CY
    a34 = (1 / V) * g * (cos(alpha) * sin(beta) * sin(theta) + cos(beta) * sin(phi) * cos(theta) - sin(alpha) * sin(beta) * &
    & cos(phi) * cos(theta))

    F(3) = a31 + a32 + a33 + a34

    ! Moments in body axis
    a41 = ((Iyy - Ixx - Izz) * Ixz / (Ixz**2 - Ixx * Izz)) * p * q
    a42 = ((Ixz**2 + Izz**2 - Izz * Iyy) / (Ixz**2 - Ixx * Izz)) * q * r
    a43 = (Izz / (Ixx * Izz - Ixz**2)) * 0.5 * rho * V**2 * S * b * Cl
    a44 = (Ixz / (Ixx * Izz - Ixz**2)) * 0.5 * rho * V**2 * S * b * Cn

    F(4) = a41 + a42 + a43 + a44

    a51 = ((Izz - Ixx) / Iyy) * r * p
    a52 = -(Ixz / Iyy) * (p**2 - r**2)
    a53 = (1 / Iyy) * 0.5 * rho * V**2 * S * cbar * Cm

    F(5) = a51 + a52 + a53

    a61 = ((Ixx * Iyy - Ixx**2 - Ixz**2) / (Ixz**2 - Ixx * Izz)) * p * q
    a62 = (((Ixx - Iyy + Izz) * Ixz) / (Ixz**2 - Ixx * Izz)) * q * r
    a63 = (Ixz / (Ixx * Izz - Ixz**2)) * 0.5 * rho * V**2 * S * b * Cl
    a64 = (Ixx / (Ixx * Izz - Ixz**2)) * 0.5 * rho * V**2 * S * b * Cn

    F(6) = a61 + a62 + a63 + a64

    ! Kinematic Equations
    F(7) = p + (q * sin(phi) + r * cos(phi)) * (sin(theta) / cos(theta))
    F(8) = q * cos(phi) - r * sin(phi)
    !F(9) = (q * sin(phi) + r * cos(phi)) / cos_theta
    ! PRINT *, 'U', U
    ! PRINT *, 'PAR:', del_e,del_a,del_r, eta

END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----

IMPLICIT NONE
INTEGER :: NDIM
DOUBLE PRECISION :: T
DOUBLE PRECISION :: U(NDIM), PAR(*)


  ! Initial conditions for states (U) and control parameters (PAR)
  U(1) = 250.0     ! V: Velocity (m/s)
  U(2) = 5.35E-002  ! alpha: Angle of attack (rad)
  U(3) = 0.0       ! beta: Sideslip angle (rad)
  U(4) = 0.0       ! p: Roll rate (rad/s)
  U(5) = 0.0       ! q: Pitch rate (rad/s)
  U(6) = 0.0       ! r: Yaw rate (rad/s)
  U(7) = 0.0       ! phi: Roll angle (rad)
  U(8) = 5.35E-002  ! theta: Pitch angle (rad)

  ! Control parameters (for example: elevator, aileron, rudder deflections)
  PAR(1) = 0.090342   ! del_e: Elevator deflection (rad)
  PAR(2) = 0.0    ! del_a: Aileron deflection (rad)
  PAR(3) = 0.0    ! del_r: Rudder deflection (rad)
  PAR(4) = 0.2488    ! eta: Throttle setting (0 to 1)
  ! PRINT *, 'U:', U
  ! PRINT *, 'PAR:', PAR(1)

END SUBROUTINE STPNT


SUBROUTINE BCND 
END SUBROUTINE BCND

SUBROUTINE ICND 
END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT

SUBROUTINE PVLS
END SUBROUTINE PVLS
