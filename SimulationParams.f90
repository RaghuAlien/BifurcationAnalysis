module SimulationParams
  implicit none

  ! Physical constants
  DOUBLE PRECISION :: m, g, Ixx, Iyy, Izz, Ixz, rho_sl, S, b, cbar, max_thrust, alt, Tm, rho
  public :: m, g, Ixx, Iyy, Izz, Ixz, rho_sl, S, b, cbar, max_thrust, alt, Tm, rho

  ! Initial conditions and conversion factors
  DOUBLE PRECISION, parameter :: d2r = 3.141592653589793 / 180.0 ! degrees to radians
  DOUBLE PRECISION, parameter :: r2d = 180.0 / 3.141592653589793 ! radians to degrees

contains

  subroutine InitializeParameters()
    ! Initialize the geometrical and physical parameters
    m = 16374.685
    g = 9.81         ! m/s^2
    Ixx = 30889.2149
    Iyy = 239654.7110
    Izz = 259898.2683
    Ixz = -3124.2985
    rho_sl = 1.225   ! kg/m^3, density at sea level
    S = 37.1612
    b = 11.3995
    cbar = 3.5113
    max_thrust = 142393.1779
    alt = 8000.0     ! Altitude in meters
    Tm =   77810.706875556760     
    rho = 0.52497499767575950

  end subroutine InitializeParameters

end module SimulationParams
