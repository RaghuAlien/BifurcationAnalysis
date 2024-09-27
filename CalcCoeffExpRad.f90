MODULE CalcCoeffExpRadModule
  USE DataReaderModule    ! Reads aerodynamic coefficients from a file
  USE SimulationParams

  IMPLICIT NONE
  DOUBLE PRECISION :: CL0_a, CLalpha, CLift_a_de, CLalpha_dot, CLift_a_q, CD0_a, CD_a_q, CD_a_de, CDalpha
  DOUBLE PRECISION :: CY_a_b, CY_a_dr, CY_a_da, CY_a_r, CY_a_p, CY_a_de, Cl_a_b, Cl_a_dr, Cl_a_da, Cl_a_r, Cl_a_p
  DOUBLE PRECISION :: Cl_ac_de, Cm0_a, Cmalpha, Cm_a_de, Cmalpha_dot, Cm_a_q, Cn_a_b, Cn_a_dr, Cn_a_da, Cn_a_r, Cn_a_p, Cn_a_de
  
  
  
CONTAINS
  SUBROUTINE Calc_Coeff_exp(U, PAR, CFor, CMom)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(8) :: U
    DOUBLE PRECISION, DIMENSION(4) :: PAR
    DOUBLE PRECISION, DIMENSION(3) :: CFor, CMom
    DOUBLE PRECISION :: alpha_dot, del_e, del_a, del_r, eta
    DOUBLE PRECISION :: CD_AC, CLiFT_ac, CY_AC, Cl_ac, Cm_ac, Cn_ac
    INTEGER :: index
    INTEGER, PARAMETER :: NUM_POINTS = 181500  ! 3001 - 1 (header)
    DOUBLE PRECISION :: CD0_a, CD_a_q, CD_a_de
    DOUBLE PRECISION :: CL0_a, CLift_a_q, CLift_a_de
    DOUBLE PRECISION :: CY_a_b, CY_a_p, CY_a_r, CY_a_de, CY_a_da, CY_a_dr
    DOUBLE PRECISION :: Cl_a_b, Cl_a_p, Cl_a_r, Cl_ac_de, Cl_a_da, Cl_a_dr
    DOUBLE PRECISION :: Cn_a_b, Cn_a_p, Cn_a_r, Cn_a_de, Cn_a_da, Cn_a_dr
    DOUBLE PRECISION :: Cm0_a, Cm_a_q, Cm_a_de
    
    DOUBLE PRECISION :: alpha_data(NUM_POINTS)
    DOUBLE PRECISION :: clift(NUM_POINTS), clift_q(NUM_POINTS), cd0(NUM_POINTS), cd_q(NUM_POINTS), cd_de(NUM_POINTS)
    DOUBLE PRECISION :: cy_b(NUM_POINTS), cy_p(NUM_POINTS), cy_r(NUM_POINTS), clift_de(NUM_POINTS), cy_de(NUM_POINTS)
    DOUBLE PRECISION :: cr_b(NUM_POINTS), cr_p(NUM_POINTS), cr_r(NUM_POINTS), cr_de(NUM_POINTS), cr_da(NUM_POINTS)
    DOUBLE PRECISION :: cn_b(NUM_POINTS), cn_p(NUM_POINTS), cn_r(NUM_POINTS), cn_de(NUM_POINTS), cn_da(NUM_POINTS)
    DOUBLE PRECISION :: cm0(NUM_POINTS), cm_q(NUM_POINTS), cm_de(NUM_POINTS), cy_da(NUM_POINTS), cy_dr(NUM_POINTS),& 
                                cr_dr(NUM_POINTS), cn_dr(NUM_POINTS)


    call InitializeParameters()
    ! Extract control inputs
    del_e = PAR(1)
    del_a = PAR(2)
    del_r = PAR(3)
    eta   = PAR(4)

    ! Read data from the file once
    CALL read_data('aero_derivatives_data.txt', alpha_data, clift, clift_q, cd0, cd_q, cd_de, &
                       cy_b, cr_b, cr_p, cr_r, cm0, cm_q, cn_b, cn_p, cn_r, &
                       cy_p, cy_r, clift_de, cy_da, cy_dr, cr_da, cn_da, &
                       cr_dr, cn_dr, cm_de, cy_de, cr_de, cn_de)

    ! Find the index closest to the current alpha value
    ! 
    index = 1
    DO WHILE (index <= num_points)
        IF (alpha_data(index) >= U(2)) EXIT
        index = index + 1
    END DO

    ! Ensure index is within bounds
    IF (index > num_points) THEN
        index = num_points
    END IF
    print *, 'Index: ', index , 'alpha: ', U(2)
    ! Assign aerodynamic coefficients
    CD0_a = cd0(index)
    CD_a_q = cd_q(index)
    CD_a_de = cd_de(index)
    CL0_a = clift(index)
    CLift_a_q = clift_q(index)
    CLift_a_de = clift_de(index)
    CY_a_b = cy_b(index)
    CY_a_p = cy_p(index)
    CY_a_r = cy_r(index)
    CY_a_de = cy_de(index)
    CY_a_da = cy_da(index)
    CY_a_dr = cy_dr(index)
    Cl_a_b = cr_b(index)
    Cl_a_p = cr_p(index)
    Cl_a_r = cr_r(index)
    Cl_ac_de = cr_de(index)
    Cl_a_da = cr_da(index)
    Cl_a_dr = cr_dr(index)
    Cn_a_b = cn_b(index)
    Cn_a_p = cn_p(index)
    Cn_a_r = cn_r(index)
    Cn_a_de = cn_de(index)
    Cn_a_da = cn_da(index)
    Cn_a_dr = cn_dr(index)
    Cm0_a = cm0(index)
    Cm_a_q = cm_q(index)
    Cm_a_de = cm_de(index)
    
    CLalpha = 0.0
    CLalpha_dot = 0.0
    CDalpha = 0.0
    Cmalpha = 0.0
    Cmalpha_dot = 0.0
    alpha_dot = 0.0

    ! Calculate aerodynamic coefficients
    CD_AC = CD0_a + CDalpha * U(2) + CD_a_q * (cbar / 2) * U(5) / U(1) + CD_a_de * del_e
    CLiFT_ac = CL0_a + CLalpha * U(2) + CLalpha_dot * (cbar / (2 * U(1))) * alpha_dot + CLift_a_q * & 
              & (cbar / 2) * U(5) / U(1) + CLift_a_de * del_e
    CY_AC = CY_a_b * U(3) + CY_a_p * (b / 2) * U(4) / U(1) + CY_a_r * (b / 2) * U(6) / U(1) + CY_a_de & 
            & * del_e + CY_a_da * del_a + CY_a_dr * del_r

    Cl_ac = Cl_a_b * U(3) + Cl_a_p * (b / 2) * U(4) / U(1) + Cl_a_r * (b / 2) * U(6) / U(1) + Cl_ac_de & 
            & * del_e + Cl_a_da * del_a + Cl_a_dr * del_r
    Cm_ac = Cm0_a + Cmalpha * U(2) + Cmalpha_dot * (cbar / (2 * U(1))) * alpha_dot + Cm_a_q * & 
            & (cbar / 2) * U(5) / U(1) + Cm_a_de * del_e
    Cn_ac = Cn_a_b * U(3) + Cn_a_p * (b / 2) * U(4) / U(1) + Cn_a_r * (b / 2) * U(6) / U(1) + & 
            & Cn_a_de * del_e + Cn_a_da * del_a + Cn_a_dr * del_r

    ! Output arrays
    CFor(1) = CD_AC
    CFor(2) = CY_AC
    CFor(3) = CLiFT_ac


    CMom(1) = Cl_ac
    CMom(2) = Cm_ac
    CMom(3) = Cn_ac

  ! PRINT *, "Calculated Aerodynamic Forces (CFor):"
  ! PRINT *, "CD_AC  = ", CFor(1)
  ! PRINT *, "CY_AC  = ", CFor(2)
  ! PRINT *, "CLiFT_ac = ", CFor(3)

  ! PRINT *, "Calculated Aerodynamic Moments (CMom):"
  ! PRINT *, "Cl_ac  = ", CMom(1)
  ! PRINT *, "Cm_ac  = ", CMom(2)
  ! PRINT *, "Cn_ac  = ", CMom(3)

  END SUBROUTINE Calc_Coeff_exp
END MODULE CalcCoeffExpRadModule
