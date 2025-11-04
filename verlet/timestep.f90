PROGRAM verlet
  IMPLICIT NONE
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  INTEGER :: n, i, j
  REAL (KIND=wp) :: t, tau, m, fx, fy, fz, K, U, E
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: x, vx, y, vy, z, vz
  ALLOCATE ( x(n) , y(n), z(n) , vx(n) , vy(n) , vz(n) )
  t = 120_wp
  m = 1.0_wp
  fx = 0.0_wp
  fy = 0.1_wp
  fz = 0.0_wp
  vx(1) = 0.0_wp
  vy(1) = 0.0_wp
  vz(1) = 0.0_wp
  x(1) = 0.0_wp
  y(1) = 0.0_wp
  z(1) = 0.0_wp

  DO j = 0, 70, 2
    tau = 0.1_wp + j
    n = t / tau 
    DO i = 1, n
      x(i+1) = x(i)+ tau * vx(i) + tau**2 * fx / 2 * m
      vx(i+1) = vx(i)+ tau / (2 * m) * (fx+fx)    
      y(i+1) = y(i) + tau * vy(i) + tau**2 * fy / 2 * m
      vy(i+1) = vy(i) + tau / (2 * m) * (fy+fy)
      z(i+1) = z(i) + tau * vz(i) + tau**2 * fz / 2 * m
      vz(i+1) = vz(i) + tau / (2 * m) * (fz + fz)
    ENDDO
    K = (m * vy(n)**2) / 2
    U = -0.1_wp * y(n)
    E = K + U
    PRINT *, n, K, U, E, y(n)
  ENDDO
END PROGRAM verlet
