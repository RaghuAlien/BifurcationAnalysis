MODULE CalcCoeffConstRadModule
  USE DataReaderModule    ! Reads aerodynamic coefficients from a file
  USE SimulationParams

  IMPLICIT NONE
  DOUBLE PRECISION :: CL0_a, CLalpha, CLift_a_de, CLalpha_dot, CLift_a_q, CD0_a, CD_a_q, CD_a_de, CDalpha
  DOUBLE PRECISION :: CY_a_b, CY_a_dr, CY_a_da, CY_a_r, CY_a_p, CY_a_de, Cl_a_b, Cl_a_dr, Cl_a_da, Cl_a_r, Cl_a_p
  DOUBLE PRECISION :: Cl_ac_de, Cm0_a, Cmalpha, Cm_a_de, Cmalpha_dot, Cm_a_q, Cn_a_b, Cn_a_dr, Cn_a_da, Cn_a_r, Cn_a_p, Cn_a_de
  
  
  
CONTAINS
  SUBROUTINE Calc_Coeff_Const(U, PAR, CFor, CMom)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(8) :: U
    DOUBLE PRECISION, DIMENSION(4) :: PAR
    DOUBLE PRECISION, DIMENSION(3) :: CFor, CMom
    DOUBLE PRECISION :: alpha_dot, del_e, del_a, del_r, eta
    DOUBLE PRECISION :: CD_AC, CLiFT_ac, CY_AC, Cl_ac, Cm_ac, Cn_ac
    INTEGER :: index
    INTEGER, PARAMETER :: NUM_POINTS = 3000  ! 3001 - 1 (header)
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
    ! Assign aerodynamic coefficients 200mps

    CD0_a =   0.0287
	CD_a_q =  -0.0008
	CD_a_de =  -0.0021
	CL0_a =   0.1805
	CLift_a_q =   0.0761
	CLift_a_de =   0.0136
	CY_a_b =  -0.0184
	CY_a_p =  -0.0002
	CY_a_r =   0.0034
	CY_a_de =   0.0000
	CY_a_da =  -0.0005
	CY_a_dr =   0.0037
	Cl_a_b =  -0.0013
	Cl_a_p =  -0.0071
	Cl_a_r =   0.0016
	Cl_ac_de =   0.0000
	Cl_a_da =   0.0012
	Cl_a_dr =   0.0002
	Cn_a_b =   0.0016
	Cn_a_p =  -0.0012
	Cn_a_r =  -0.0031
	Cn_a_de =  -0.0000
	Cn_a_da =  -0.0000
	Cn_a_dr =  -0.0012
	Cm0_a =   0.0041
	Cm_a_q =  -0.0852
	Cm_a_de =  -0.0150

    
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

  END SUBROUTINE Calc_Coeff_Const
END MODULE CalcCoeffConstRadModule
