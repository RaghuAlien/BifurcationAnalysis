module AtmosphereModule
  use SimulationParams   ! Import the SimulationParams module
  implicit none

  ! Parameters for the atmosphere model
  DOUBLE PRECISION, parameter :: gas_const = 287.0     ! Specific gas constant for air (J/(kg·K))
  DOUBLE PRECISION, parameter :: T0 = 288.15   ! Sea-level standard temperature (K)
  DOUBLE PRECISION, parameter :: P0 = 101325.0 ! Sea-level standard pressure (Pa)
  DOUBLE PRECISION, parameter :: L = 0.0065    ! Temperature lapse rate (K/m)
  DOUBLE PRECISION, parameter :: T_strat = 216.65 ! Temperature in the stratosphere (K)
  DOUBLE PRECISION, parameter :: P_strat = 22632.0 ! Pressure at 11 km in the stratosphere (Pa)

contains

  subroutine AtmModel(rho)
    implicit none
    DOUBLE PRECISION, intent(out) :: rho               ! Air density (kg/m^3)
    DOUBLE PRECISION :: T_atm, p_atm
    DOUBLE PRECISION :: altitude

    ! Directly access the altitude from SimulationParams
    altitude = alt  ! Access the public variable directly

    if (altitude <= 11000.0) then
      T_atm = T0 - L * altitude
      p_atm = P0 * (T_atm / T0)**(g / (gas_const * L))
    else
      T_atm = T_strat
      p_atm = P_strat * exp(-g * (altitude - 11000.0) / (gas_const * T_atm))
    end if

    rho = p_atm / (gas_const * T_atm)

  end subroutine AtmModel

end module AtmosphereModule
