MODULE EOMWindRefModule
  USE ThrustModel
  USE SimulationParams
  USE CalcCoeffExpModule
  IMPLICIT NONE
CONTAINS
  SUBROUTINE eom_WindRef(x, U, CFor, CMom, F)
    IMPLICIT NONE
    REAL, DIMENSION(8) :: x, F
    REAL, DIMENSION(4) :: U
    REAL, DIMENSION(3) :: CFor, CMom
    REAL :: del_e, del_a, del_r, eta
    REAL :: CD, CY, CLIFT, Cl, Cm, Cn, Tm

    ! Intermediate calculations
    REAL :: a11, a12, a13, a21, a22, a23, a24
    REAL :: a31, a32, a33, a34, a41, a42, a43, a44
    REAL :: a51, a52, a53, a61, a62, a63, a64
    
    CALL InitializeParameters()
    ! CALL Calc_Coeff_exp(x, U, CFor, CMom)

    Tm = Thrust_Model(max_thrust)

    ! Control inputs
    del_e = U(1)
    del_a = U(2)
    del_r = U(3)
    eta   = U(4)

    ! Force coefficients
    CD = CFor(1)
    CY = CFor(2)
    CLIFT = CFor(3)

    ! Moment coefficients
    Cl = CMom(1)
    Cm = CMom(2)
    Cn = CMom(3)


    ! Print *,'rho:', rho
    ! Print *,'Tm:', Tm
    !  print *, 'cbar:',cbar

    ! Calculate thrust

    ! Force equations in body frame
    a11 = (1.0/m) * eta * Tm * COS(x(2)) * COS(x(3))
    a12 = -(1.0/m) * 0.5 * rho * (x(1)**2) * S * CD
    a13 = -g * (SIN(x(8)) * COS(x(2)) * COS(x(3)) - COS(x(8)) * SIN(x(7)) * SIN(x(3)) - SIN(x(2)) * & 
      & COS(x(3)) * COS(x(7)) * COS(x(8)))
    F(1)  = a11 + a12 + a13

    a21 = x(5) - (x(4) * COS(x(2)) + x(6) * SIN(x(2))) * TAN(x(3))
    a22 = -(1.0/m) * (1.0/x(1)) * (1.0/COS(x(3))) * eta * Tm * SIN(x(2))
    a23 = -(1.0/m) * (1.0/COS(x(3))) * 0.5 * rho * x(1) * S * CLIFT
    a24 = (1.0/x(1)) * (1.0/COS(x(3))) * g * (SIN(x(2)) * SIN(x(8)) + COS(x(2)) * COS(x(7)) * COS(x(8)))
    F(2)  = a21 + a22 + a23 + a24

    a31 = x(4) * sin(x(2)) - x(6) * cos(x(2))
    a32 = -(1 / m) * (1 / X(1)) * eta * Tm * cos(x(2)) * sin(x(3))
    a33 = -(1 / m) * 0.5 * rho * X(1) * S * CY
    a34 = (1 / X(1)) * g * (cos(x(2)) * sin(x(3)) * sin(x(8)) + cos(x(3)) * sin(x(7)) * cos(x(8)) - sin(x(2)) * sin(x(3)) * &
    & cos(x(7)) * cos(x(8)))

    F(3) = a31 + a32 + a33 + a34

    ! Moments in body axis
    a41 = ((Iyy - Ixx - Izz) * Ixz / (Ixz**2 - Ixx * Izz)) * x(4) * x(5)
    a42 = ((Ixz**2 + Izz**2 - Izz * Iyy) / (Ixz**2 - Ixx * Izz)) * x(5) * x(6)
    a43 = (Izz / (Ixx * Izz - Ixz**2)) * 0.5 * rho * X(1)**2 * S * b * Cl
    a44 = (Ixz / (Ixx * Izz - Ixz**2)) * 0.5 * rho * X(1)**2 * S * b * Cn

    F(4) = a41 + a42 + a43 + a44

    a51 = ((Izz - Ixx) / Iyy) * x(6) * x(4)
    a52 = -(Ixz / Iyy) * (x(4)**2 - x(6)**2)
    a53 = (1 / Iyy) * 0.5 * rho * X(1)**2 * S * cbar * Cm

    F(5) = a51 + a52 + a53

    a61 = ((Ixx * Iyy - Ixx**2 - Ixz**2) / (Ixz**2 - Ixx * Izz)) * x(4) * x(5)
    a62 = (((Ixx - Iyy + Izz) * Ixz) / (Ixz**2 - Ixx * Izz)) * x(5) * x(6)
    a63 = (Ixz / (Ixx * Izz - Ixz**2)) * 0.5 * rho * X(1)**2 * S * b * Cl
    a64 = (Ixx / (Ixx * Izz - Ixz**2)) * 0.5 * rho * X(1)**2 * S * b * Cn

    F(6) = a61 + a62 + a63 + a64

    ! Kinematic Equations
    F(7) = x(4) + (x(5) * sin(x(7)) + x(6) * cos(x(7))) * (sin(x(8)) / cos(x(8)))
    F(8) = x(5) * cos(x(7)) - x(6) * sin(x(7))
    !F(9) = (x(5) * sin(x(7)) + x(6) * cos(x(7))) / cos_theta

    ! (More calculations based on the rest of the MATLAB code)
  END SUBROUTINE eom_WindRef
END MODULE EOMWindRefModule
