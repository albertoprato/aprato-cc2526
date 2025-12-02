PROGRAM verlet
  IMPLICIT NONE
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(p=13, r=300)
  INTEGER :: n, i, j
  REAL(KIND=wp) :: t, tau, m, fx, fy, fz, K, U, E
  REAL(KIND=wp) :: K_init, U_init, E_init, K_final, U_final, E_final, dE
  REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: x, vx, y, vy, z, vz
  
  t = 120.0_wp
  m = 1.0_wp
  fx = 0.0_wp
  fy = 0.1_wp
  fz = 0.0_wp

  ! Tabella Output
  PRINT '(A)', REPEAT("=", 95)
  PRINT '(A)', "   Verlet Algorithm - Classical Trajectory with Constant Force"
  PRINT '(A)', REPEAT("=", 95)
  PRINT '(A6, A10, A15, A15, A15, A20)', &
        "Steps", "Tau", "K (final)", "U (final)", "y(final)", "dE (E_fin-E_init)"
  PRINT '(A)', REPEAT("-", 95)
  
  ! Loop sui timestep
  DO j = 0, 20, 1
    tau = 1.0_wp + REAL(j, KIND = wp)
    n = NINT(t / tau)
    
    ALLOCATE(x(n+1), y(n+1), z(n+1), vx(n+1), vy(n+1), vz(n+1))
    
    ! Condizioni iniziali
    x(1) = 0.0_wp
    y(1) = 0.0_wp
    z(1) = 0.0_wp
    vx(1) = 0.0_wp
    vy(1) = 0.0_wp
    vz(1) = 0.0_wp
    
    ! Energia iniziale (t=0)
    K_init = kinetic_energy(m, vx(1), vy(1), vz(1))
    U_init = 0.0_wp
    E_init = K_init + U_init
    
    U = U_init

    ! Algoritmo Verlet
    DO i = 1, n
      x(i+1) = x(i) + tau * vx(i) + tau**2 * fx / (2.0_wp * m)
      vx(i+1) = vx(i) + tau / (2.0_wp * m) * (fx + fx)
      
      y(i+1) = y(i) + tau * vy(i) + tau**2 * fy / (2.0_wp * m)
      vy(i+1) = vy(i) + tau / (2.0_wp * m) * (fy + fy)
      
      z(i+1) = z(i) + tau * vz(i) + tau**2 * fz / (2.0_wp * m)
      vz(i+1) = vz(i) + tau / (2.0_wp * m) * (fz + fz)

      U = U + potential_energy(fx, fy, fz, x(i), y(i), z(i), x(i+1), y(i+1), z(i+1))
    END DO

    ! Energie finali (t=T)
    K_final = kinetic_energy(m, vx(n+1), vy(n+1), vz(n+1))
    U_final = U
    E_final = K_final + U_final

    ! Differenza di energia totale
    dE = E_final - E_init
    
    ! Output formattato
    PRINT '(I6, 1X, F9.2, 1X, ES15.6, 1X, ES15.6, 1X, ES15.6, 1X, ES20.6)', &
           n, tau, K_final, U_final, y(n+1), dE

    DEALLOCATE(x, y, z, vx, vy, vz)

  END DO
  
  PRINT '(A)', REPEAT("=", 95)
  PRINT '(A)', ""


CONTAINS

   FUNCTION kinetic_energy(mass, vx_k, vy_k, vz_k) RESULT(T_k)
       IMPLICIT NONE
       INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
       REAL(KIND=wp), INTENT(IN) :: mass, vx_k, vy_k, vz_k
       REAL(KIND=wp) :: T_k
    
       T_k = 0.5_wp * mass * (vx_k**2 + vy_k**2 + vz_k**2)
   END FUNCTION kinetic_energy

   FUNCTION potential_energy(fx, fy, fz, x_k, y_k, z_k, x_kp1, y_kp1, z_kp1) RESULT(dU)
       IMPLICIT NONE
       INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
       REAL(KIND=wp), INTENT(IN) :: fx, fy, fz
       REAL(KIND=wp), INTENT(IN) :: x_k, y_k, z_k, x_kp1, y_kp1, z_kp1
       REAL(KIND=wp) :: dU
    
       dU = fx * (x_kp1 - x_k) + fy * (y_kp1 - y_k) + fz * (z_kp1 - z_k)
   END FUNCTION potential_energy

END PROGRAM verlet
