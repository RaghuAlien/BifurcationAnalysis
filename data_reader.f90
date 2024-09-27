MODULE DataReaderModule
  IMPLICIT NONE
  INTEGER, PARAMETER :: NUM_VARS = 27
  INTEGER, PARAMETER :: NUM_POINTS = 181500  ! 3001 - 1 (header)

  CONTAINS

  SUBROUTINE read_data(filename, alpha, clift, clift_q, cd0, cd_q, cd_de, &
                       cy_b, cr_b, cr_p, cr_r, cm0, cm_q, cn_b, cn_p, cn_r, &
                       cy_p, cy_r, clift_de, cy_da, cy_dr, cr_da, cn_da, &
                       cr_dr, cn_dr, cm_de, cy_de, cr_de, cn_de)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(OUT) :: alpha(num_points), clift(num_points), clift_q(num_points), cd0(num_points), &
              cd_q(num_points), cd_de(num_points), &
             cy_b(num_points), cr_b(num_points), cr_p(num_points), cr_r(num_points), cm0(num_points), cm_q(num_points), &
             cn_b(num_points), cn_p(num_points), cn_r(num_points), cy_p(num_points), cy_r(num_points), clift_de(num_points), &
             cy_da(num_points), cy_dr(num_points), cr_da(num_points), cn_da(num_points), cr_dr(num_points), cn_dr(num_points), &
             cm_de(num_points), cy_de(num_points), cr_de(num_points), cn_de(num_points)
    INTEGER :: i, ios
    CHARACTER(LEN=500) :: line

    ! Open the data file
    OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
      PRINT *, 'Error opening file!'
      STOP
    END IF

    ! Read and discard the header line
    READ(10, '(A)') line

    ! Read the data
    DO i = 1, NUM_POINTS
      READ(10, *) alpha(i), clift(i), clift_q(i), cd0(i), cd_q(i), cd_de(i), &
                    cy_b(i), cr_b(i), cr_p(i), cr_r(i), cm0(i), cm_q(i), cn_b(i), &
                    cn_p(i), cn_r(i), cy_p(i), cy_r(i), clift_de(i), cy_da(i), &
                    cy_dr(i), cr_da(i), cn_da(i), cr_dr(i), cn_dr(i), cm_de(i), &
                    cy_de(i), cr_de(i), cn_de(i)
    END DO

    CLOSE(10)
  END SUBROUTINE read_data

END MODULE DataReaderModule
