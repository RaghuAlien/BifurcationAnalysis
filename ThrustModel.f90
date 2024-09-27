MODULE ThrustModel
  USE AtmosphereModule
  IMPLICIT NONE
  DOUBLE PRECISION :: rho
CONTAINS
  FUNCTION Thrust_Model(max_thrust) RESULT(Tm)
    IMPLICIT NONE
    DOUBLE PRECISION :: Tm, max_thrust
    
    ! Call the standard atmosphere function to get rho
    call AtmModel(rho)
    ! Calculate the thrust
    Tm = ((rho / rho_sl)**0.7) * max_thrust * (1.0 - EXP((alt - 17000.0) / 2000.0))
    
  END FUNCTION Thrust_Model
END MODULE ThrustModel
