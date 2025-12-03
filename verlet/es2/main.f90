PROGRAM main
  USE kinds, ONLY: wp => dp                                           
  USE force_module                                                        
  IMPLICIT NONE
  
  INTEGER :: nk, n, i, k, step
  REAL(KIND=wp) :: tau, sigma, epsilon

  REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: pos, vel, force, force_new
  REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: mass
  
  REAL(KIND=wp) :: inp_m, inp_x, inp_y, inp_z, inp_vx, inp_vy, inp_vz

  ! Read Input
  OPEN(UNIT=10, FILE='input.txt', STATUS='old')
  READ(10, *) nk, tau
  READ(10, *) sigma, epsilon
  READ(10, *) n

  ALLOCATE(pos(n, 3), vel(n, 3), force(n, 3), force_new(n, 3))
  ALLOCATE (mass(n))
  
  DO i = 1, n
    READ(10, *) inp_m, inp_x, inp_y, inp_z, inp_vx, inp_vy, inp_vz
    mass(i)   = inp_m
    pos(i, 1) = inp_x
    pos(i, 2) = inp_y
    pos(i, 3) = inp_z
    vel(i, 1) = inp_vx
    vel(i, 2) = inp_vy
    vel(i, 3) = inp_vz
  END DO

  CLOSE(10)
  
  PRINT *, "------------------------------------------------"
  PRINT *, "System initialized with ", n, " particles."
  PRINT *, "Steps:", nk, " dt:", tau
  PRINT *, "------------------------------------------------"

  ! Calculate the forces at t = 0
  CALL force_calculation(n, sigma, epsilon, pos, force) 

  ! Verlet Algorithm
  OPEN(UNIT=20, FILE='trajectory.xyz', STATUS='replace')
  
  ! Step 1
  DO step = 1, nk

    DO i = 1, n
      DO k = 1, 3
        pos(i, k) = pos(i, k) + (vel(i, k) * tau) + (force(i, k) / (2.0_wp * mass(i))) * (tau**2)
      END DO
    END DO
    
    ! Step 2
    CALL force_calculation(n, sigma, epsilon, pos, force)
    
    ! Step 3
    DO i = 1, n
      DO k = 1, 3
        vel(i, k) = vel(i, k) + (tau / (2.0_wp * mass(i))) * (force(i, k) + force_new(i, k))
      END DO
    END DO

    ! Step 4
    force = force_new

    ! Output Results (e.g., every 100 steps)
    IF (mod(step, 100) == 0) THEN

      ! Print Distance between atom 1 and 2 to check oscillation
      ! WRITE(*, '(A, I5, A, F15.8)') "Step: ", step, " | Distance 1-2: ", &
      !      sqrt((pos(1,1)-pos(2,1))**2 + (pos(1,2)-pos(2,2))**2 + (pos(1,3)-pos(2,3))**2)
             
      ! VMD XYZ format
      WRITE(20, *) n
      WRITE(20, '(A, I8)') "Step: ", step
      DO i = 1, n
        WRITE(20, '(A, 3F15.8)') "Ne ", pos(i, 1), pos(i, 2), pos(i, 3)
      END DO 
    END IF

  END DO

  CLOSE(20)
  
  PRINT *, "------------------------------------------------"
  PRINT *, "Simulation completed."
  PRINT *, "Final Z coordinate of Atom 1: ", pos(1, 3)
  PRINT *, "------------------------------------------------"
  
  DEALLOCATE(pos, vel, force, force_new, mass)

END PROGRAM main
