PROGRAM cda
  
  USE kinds, ONLY: wp => dp
  USE cubes

  IMPLICIT NONE
 
  TYPE (cube) :: cubeA, cubeB, cubeAB, cubeREF, DELTAcube
  REAL (KIND=wp), DIMENSION (:), ALLOCATABLE ::cube_integral, cubeCDZ
  REAL (KIND=wp) :: z_axis
  INTEGER :: i
  
  CHARACTER (LEN=*), PARAMETER :: file_a  = '../test/CuCO+/a.cube'
  CHARACTER (LEN=*), PARAMETER :: file_b  = '../test/CuCO+/b.cube'
  CHARACTER (LEN=*), PARAMETER :: file_ab = '../test/CuCO+/ab.cube'

  CALL cube_get(cubeA, file_a)
  CALL cube_get(cubeB, file_b)
  CALL cube_get(cubeAB, file_ab)

  cubeREF = cubeA + cubeB
  
  DELTAcube = cubeAB - cubeREF
  
  cube_integral = cube_int(DELTAcube)

  cubeCDZ = cube_cdz(DELTAcube)

  OPEN(UNIT=12, FILE="output.txt", STATUS='replace', ACTION="WRITE")
  
  DO i=1, DELTAcube%nz
    
    z_axis = DELTAcube%zmin + (i-1) * DELTAcube%dz

    WRITE(UNIT=12, FMT='(3(ES14.6, 2X))') z_axis, cube_integral(i), cubeCDZ(i)

  ENDDO  

  CLOSE(12)

END PROGRAM cda
