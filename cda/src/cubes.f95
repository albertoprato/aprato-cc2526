MODULE cubes
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE


  TYPE, PUBLIC :: cube
    ! PRIVATE
    CHARACTER (LEN=72) :: str1
    CHARACTER (LEN=72) :: str2
    REAL (KIND=wp) :: xmin, ymin, zmin, dx, dy, dz
    INTEGER :: nx, ny, nz, natom
    INTEGER, DIMENSION(:), POINTER :: zahl
    REAL (KIND=wp), DIMENSION(:), POINTER :: chrg, x, y, z
    REAL (KIND=wp), DIMENSION(:), POINTER :: array
  END TYPE cube

  PUBLIC :: cube_get, &
            cube_add, &
            cube_sub, &
            cube_int, &
            cube_cdz, &
            cube_del


  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cube_add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cube_sub
  END INTERFACE



  CONTAINS



  SUBROUTINE cube_get (mycube, infile)
    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE (cube), INTENT(OUT) :: mycube
     
    INTEGER :: i
    
    REAL (KIND=wp) :: a !dummy variable to store zeros in the dx, dy, dz matrix
    

    OPEN(UNIT = 10, FILE=infile, STATUS='old')

    READ(10, *) mycube%str1
    READ(10, *) mycube%str2
    READ(10, *) mycube%natom, mycube%xmin, mycube%ymin, mycube%zmin
    READ(10, *) mycube%nx, mycube%dx, a, a
    READ(10, *) mycube%ny, a, mycube%dy, a
    READ(10, *) mycube%nz, a, a, mycube%dz
    
    ALLOCATE (mycube%zahl(mycube%natom), mycube%chrg(mycube%natom), mycube%x(mycube%natom), &
               mycube%y(mycube%natom), mycube%z(mycube%natom))    
    
    DO i=1, mycube%natom
    
      READ(10,*) mycube%zahl(i), mycube%chrg(i), mycube%x(i), mycube%y(i), mycube%z(i)
    
    ENDDO
    
    ALLOCATE ( mycube%array(mycube%nx * mycube%ny * mycube%nz) )

    READ(10, *) mycube%array

    CLOSE(10)
 
  END SUBROUTINE cube_get



!==========================================================================================================



  FUNCTION cube_add (mycube1, mycube2)
    TYPE(cube) :: cube_add
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    
    INTEGER :: i
    
    ! Assegno le caratteristiche irrilevanti del cubo 1
 
    cube_add%str1 = mycube1%str1
    cube_add%str2 = mycube1%str2

    cube_add%xmin = mycube1%xmin
    cube_add%ymin = mycube1%ymin
    cube_add%zmin = mycube1%zmin
    
    cube_add%nx = mycube1%nx
    cube_add%ny = mycube1%ny
    cube_add%nz = mycube1%nz
    
    cube_add%dx = mycube1%dx
    cube_add%dy = mycube1%dy
    cube_add%dz = mycube1%dz
    
    cube_add%natom = mycube1%natom + mycube2%natom

    ALLOCATE ( cube_add%zahl(cube_add%natom), cube_add%chrg(cube_add%natom), &
               cube_add%x(cube_add%natom), cube_add%y(cube_add%natom), cube_add%z(cube_add%natom) )
    
    DO i=1, mycube1%natom
      cube_add%zahl(i) = mycube1%zahl(i)
      cube_add%chrg(i) = mycube1%chrg(i)
      cube_add%x(i) = mycube1%x(i)
      cube_add%y(i) = mycube1%y(i)
      cube_add%z(i) = mycube1%z(i)
    ENDDO

    DO i=1, mycube2%natom
      cube_add%zahl(mycube1%natom + i) = mycube2%zahl(i)
      cube_add%chrg(mycube1%natom + i) = mycube2%chrg(i)
      cube_add%x(mycube1%natom + i) = mycube2%x(i)
      cube_add%y(mycube1%natom + i) = mycube2%y(i)
      cube_add%z(mycube1%natom + i) = mycube2%z(i)
    ENDDO

    ALLOCATE ( cube_add%array(cube_add%nx * cube_add%ny * cube_add%nz) )

    cube_add%array = mycube1%array + mycube2%array 

  END FUNCTION cube_add



!==========================================================================================================



  FUNCTION cube_sub (mycube1, mycube2)
    TYPE(cube) :: cube_sub
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
        
    cube_sub%str1 = mycube1%str1
    cube_sub%str2 = mycube1%str2

    cube_sub%xmin = mycube1%xmin
    cube_sub%ymin = mycube1%ymin
    cube_sub%zmin = mycube1%zmin
    
    cube_sub%nx = mycube1%nx
    cube_sub%ny = mycube1%ny
    cube_sub%nz = mycube1%nz
    
    cube_sub%dx = mycube1%dx
    cube_sub%dy = mycube1%dy
    cube_sub%dz = mycube1%dz
    
    cube_sub%natom = mycube1%natom 

    ALLOCATE ( cube_sub%zahl(cube_sub%natom), cube_sub%chrg(cube_sub%natom), &
               cube_sub%x(cube_sub%natom), cube_sub%y(cube_sub%natom), cube_sub%z(cube_sub%natom) )
  
      cube_sub%zahl = mycube1%zahl
      cube_sub%chrg = mycube1%chrg
      cube_sub%x = mycube1%x
      cube_sub%y = mycube1%y
      cube_sub%z = mycube1%z

    ALLOCATE ( cube_sub%array(cube_sub%nx * cube_sub%ny * cube_sub%nz) )

    cube_sub%array = mycube1%array - mycube2%array 
  
  END FUNCTION cube_sub



!==========================================================================================================


  
  SUBROUTINE  cube_unroll (mycube, array3d)

    TYPE (cube), INTENT(IN) :: mycube
    REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: array3d
    INTEGER :: i, i_x, i_y, i_z

    ALLOCATE(array3d(mycube%nx, mycube%ny, mycube%nz))
    
    i = 1

    DO i_x = 1, mycube%nx
        DO i_y = 1, mycube%ny
            DO i_z = 1, mycube%nz

                array3d(i_x, i_y, i_z) = mycube%array(i)

                i = i + 1

            ENDDO
        ENDDO
    ENDDO 

  END SUBROUTINE cube_unroll



!==========================================================================================================
  


  FUNCTION cube_int (mycube) ! Integral in xy space
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: cube_int ! Dp(z)
    REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: array3d
    TYPE (cube), INTENT(IN) :: mycube
    INTEGER :: i_x, i_y, i_z
    
    ALLOCATE(cube_int(mycube%nz))

    CALL cube_unroll(mycube, array3d)
        
    DO i_z = 1, mycube%nz

      cube_int(i_z) = 0.0_wp

        DO i_x = 1, mycube%nx
          DO i_y = 1, mycube%ny

            cube_int(i_z) = cube_int(i_z) + (array3d(i_x, i_y, i_z) * mycube%dx * mycube%dy)

          ENDDO
        ENDDO
    ENDDO
  END FUNCTION cube_int

  

!==========================================================================================================


  FUNCTION cube_cdz (mycube)
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: cube_integral, cube_cdz
    TYPE (cube), INTENT(IN) :: mycube

    INTEGER :: i, j

    cube_integral = cube_int(mycube)

    ALLOCATE(cube_cdz(mycube%nz))

    DO i = 1, mycube%nz
        DO j = 1, i
                cube_cdz(i) = cube_cdz(i) + (cube_integral(j) * mycube%dz) 
        ENDDO
    ENDDO
 END FUNCTION cube_cdz
  


!==========================================================================================================



  SUBROUTINE cube_del (mycube)
    TYPE (cube), INTENT(IN) :: mycube
    ! ...
  END SUBROUTINE cube_del

END MODULE cubes
