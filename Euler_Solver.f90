PROGRAM Euler_Solver
  IMPLICIT NONE

  ! Declare variables
  INTEGER :: i, NDIM
  DOUBLE PRECISION :: t, t_end, dt
  DOUBLE PRECISION, ALLOCATABLE :: U(:), U_next(:), F(:)
  DOUBLE PRECISION :: ICP(4), PAR(39)
  integer :: unit_num
  character(len=20) :: file_name
  ! Constants
  NDIM = 8 ! Number of state variables

  ! Allocate arrays
  ALLOCATE(U(NDIM), U_next(NDIM), F(NDIM))

  ! Initial conditions for states (U) and control parameters (PAR)
  U(1) = 300.00629290388014     ! V: Velocity (m/s)
  U(2) = 3.9005133574702887E-002  ! alpha: Angle of attack (rad)
  U(3) = 0.0       ! beta: Sideslip angle (rad)
  U(4) = 0.0       ! p: Roll rate (rad/s)
  U(5) = 0.0       ! q: Pitch rate (rad/s)
  U(6) = 0.0       ! r: Yaw rate (rad/s)
  U(7) = 0.0       ! phi: Roll angle (rad)
  U(8) = 3.9005133574702887E-002  ! theta: Pitch angle (rad)

  ! Control parameters (for example: elevator, aileron, rudder deflections)
  PAR(1) = 0.27731001377105713   ! del_e: Elevator deflection (rad)
  PAR(2) = 0.0    ! del_a: Aileron deflection (rad)
  PAR(3) = 0.0    ! del_r: Rudder deflection (rad)
  PAR(4) = 0.31626999378204346    ! eta: Throttle setting (0 to 1)
  
  ! Time settings
  t = 0.0         ! Initial time
  t_end = 50.0    ! End time
  dt = 0.1       ! Time step

! Open file for writing
  file_name = 'state_output.txt'
  open(unit=unit_num, file=file_name, status='replace', action='write')

  ! Main integration loop
  DO WHILE (t < t_end)
    ! Compute F at the current state
    CALL FUNC(NDIM, U, ICP, PAR, 0, F, F, F)

    ! Update state variables using Euler method
    U_next = U + dt * F
    U = U_next

    ! Update time
    t = t + dt

    ! Write the results to the file
    ! Write the time and state to the file
    WRITE(unit_num, '(F10.2, 8F10.4)') t, U
    ! Print progress to command window
    IF (MOD(INT(t / dt), 100) == 0) THEN
      PRINT *, ' Progress: ', t / t_end * 100, '%'
    END IF

  END DO

  ! Close the output file
  CLOSE(unit_num)

  ! Deallocate arrays
  DEALLOCATE(U, U_next, F)

END PROGRAM Euler_Solver
