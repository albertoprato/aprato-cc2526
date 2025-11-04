PROGRAM verlet
  IMPLICIT NONE
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(p=13, r=300)
  INTEGER :: n, i, j
  REAL(KIND=wp) :: t, tau, m, fx, fy, fz, K, U, E
  REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: x, vx, y, vy, z, vz
  
  t = 120.0_wp
  m = 1.0_wp
  fx = 0.0_wp
  fy = 0.1_wp
  fz = 0.0_wp
  
  ! Tabella per Output
  PRINT '(A)', REPEAT("=", 90)
  PRINT '(A)', "   Simulazione Verlet - Particella sotto forza costante"
  PRINT '(A)', REPEAT("=", 90)
  PRINT '(A6, A10, A15, A15, A15, A15)', &
        "Steps", "Tau", "K (kin)", "U (pot)", "E (tot)", "y(final)"
  PRINT '(A)', REPEAT("-", 90)
  
  ! Loop sui timestep
  DO j = 0, 30, 1
    tau = 0.1_wp + REAL(j, KIND=wp) / 2
    n = NINT(t / tau)
    
    ! Allocazione dinamica
    IF (ALLOCATED(x)) DEALLOCATE(x, y, z, vx, vy, vz)
    ALLOCATE(x(n+1), y(n+1), z(n+1), vx(n+1), vy(n+1), vz(n+1))
    
    ! Condizioni iniziali
    x(1) = 0.0_wp
    y(1) = 0.0_wp
    z(1) = 0.0_wp
    vx(1) = 0.0_wp
    vy(1) = 0.0_wp
    vz(1) = 0.0_wp
    
    ! Algoritmo Verlet
    DO i = 1, n
      x(i+1) = x(i) + tau * vx(i) + tau**2 * fx / (2.0_wp * m)
      vx(i+1) = vx(i) + tau / (2.0_wp * m) * (fx + fx)
      
      y(i+1) = y(i) + tau * vy(i) + tau**2 * fy / (2.0_wp * m)
      vy(i+1) = vy(i) + tau / (2.0_wp * m) * (fy + fy)
      
      z(i+1) = z(i) + tau * vz(i) + tau**2 * fz / (2.0_wp * m)
      vz(i+1) = vz(i) + tau / (2.0_wp * m) * (fz + fz)
    END DO
    
    ! Calcolo energie
    K = (m * vy(n+1)**2) / 2.0_wp
    U = -fy * y(n+1)
    E = K + U
    
    ! Output formattato
    PRINT '(I6, F10.2, 4ES15.6)', n, tau, K, U, E, y(n+1)
  END DO
  
  PRINT '(A)', REPEAT("=", 90)
  
  ! Deallocazione
  DEALLOCATE(x, y, z, vx, vy, vz)
  
END PROGRAM verlet
