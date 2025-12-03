MODULE force_module                                                       
  USE kinds, ONLY: wp => dp                                           
  IMPLICIT NONE

  PUBLIC :: force_calculation                                                                                                                                           
  CONTAINS                                                           
                                                                      
    SUBROUTINE force_calculation(n, sigma, epsilon, pos, forces, V)
      INTEGER, INTENT(IN) :: n
      REAL (KIND=wp), INTENT(IN):: sigma, epsilon
      REAL (KIND=wp), DIMENSION(n,3), INTENT(IN) :: pos
      REAL (KIND=wp), DIMENSION (n,3), INTENT(OUT) :: forces
      REAL (KIND=wp), INTENT(OUT) :: V

      INTEGER :: a, b, k
      REAL (KIND=wp) :: r2, r, r_inv, sigma_r
      REAL (KIND=wp) :: V_prime, f_component
      REAL (KIND=wp) :: term6, term12
      REAL (KIND=wp) :: diff(3) ! To store x_ab, y_ab, z_ab

      ! Initialize outputs
      forces = 0.0_wp
      V = 0.0_wp

      ! Loop over pairs of atoms
      DO a = 1, n - 1
        DO b = a + 1, n
          
          ! Calculate distance components (x_ab, y_ab, z_ab)
          r2 = 0.0_wp 
          DO k = 1, 3
            diff(k) = pos(a, k) - pos(b, k)
            r2 = r2 + diff(k)**2
          END DO

          r = sqrt(r2)
          r_inv = 1.0_wp / r

          ! Calculate V'_LJ
          sigma_r = sigma / r
          
          term6 = sigma_r**6
          term12 = term6**2
          
          V = 4 * epsilon * (term12 - term6)

          V_prime = 4 * epsilon * (- 12.0_wp * (sigma_r**12) * r_inv + 6.0_wp * (sigma_r**6) * r_inv)

          ! Calculate force components
          DO k = 1, 3
            f_component = - (diff(k) * r_inv) * V_prime

            forces(a, k) = forces(a, k) + f_component

            forces(b, k) = forces(b, k) - f_component ! Newton's 3rd Law 
          END DO

        END DO
      END DO

    END SUBROUTINE force_calculation     
                                                                      
END MODULE force_module
