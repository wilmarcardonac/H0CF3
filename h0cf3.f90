Program h0cf3 

!####################
! LOAD NEEDED MODULES 
!####################

  use healpix_types
  use healpix_modules
  use udgrade_nr, only: udgrade_ring, udgrade_nest
  use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring, remove_dipole, ang2pix_ring
  use fitstools, only: getsize_fits, input_map, output_map, fits2cl
  use head_fits
!  use fgsl
  use fiducial
  use arrays
  use functions 

!#################################
! DECLARE VARIABLES AND PARAMETERS
!#################################

    Implicit none

    Integer*4 :: m,n,i                      ! INTERGER FOR SHORT LOOPS 

    Real(kind=DP),dimension(0:DEGREE_REMOVE_DIPOLE*DEGREE_REMOVE_DIPOLE-1) :: multipoles ! SAVES MONOPOLE AND DIPOLE OF CMB MAP 
    Real(kind=DP),dimension(1:2) :: zbounds ! BOUNDS TO COMPUTE DIPOLE AND MONOPOLE

    Logical :: exist

    Character(len=4) :: x

    Character(len=80),dimension(1:60) :: header

!######################################
! INITIALIZATION OF VARIABLES AND FILES
!######################################

    open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)

    zbounds(:) = 0.d0 

    npixC = nside2npix(nsmax)

    call read_data_CF3(PATH_TO_COSMICFLOWS_DATA)

    write(UNIT_EXE_FILE,*) 'THIS PROJECT AIMS AT STUDYING THE DISCREPANCY BETWEEN LOCAL AND GLOBAL MEASUREMENTS '
    write(UNIT_EXE_FILE,*) 'OF THE HUBBLE EXPANSION RATE'

    write(UNIT_EXE_FILE,*) ' '

    write(UNIT_EXE_FILE,*) 'DATA FROM THE COSMICFLOWS-3 COMPILATION SUCCESFULLY READ '

    write(UNIT_EXE_FILE,*) ' '

    write(UNIT_EXE_FILE,*) 'N_side = ', nsmax, ' IN THE CURRENT ANALYSIS; THIS CORRESPONDS TO A NUMBER '
    write(UNIT_EXE_FILE,*) 'OF PIXELS IN THE NUMBER COUNTS MAPS EQUAL TO ', npixC

!################
! ANALYSIS STARTS 
!################

    call write_python_script_angular_distribution_galaxies()

    call write_mpi_file('angular_distribution_galaxies_CF3')

    call system('cd scripts; sbatch angular_distribution_galaxies_CF3.mpi')
!    call ang2pix_ring(1,0.d0,0.d0,i)


    close(UNIT_EXE_FILE)

End Program h0cf3




