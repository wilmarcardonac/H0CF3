Program h0cf3 

!####################
! LOAD NEEDED MODULES 
!####################

  use healpix_types
  use healpix_modules
  use udgrade_nr, only: udgrade_ring, udgrade_nest
  use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring, remove_dipole, ang2pix_ring, vec2ang
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

    Integer*8 :: i

    Real*8 :: testdipamp,test1,test2

!######################################
! INITIALIZATION OF VARIABLES AND FILES
!######################################

    i = -1

    open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)

    npixC = nside2npix(nsmax)

    call read_data_CF3(PATH_TO_COSMICFLOWS_DATA)

    write(UNIT_EXE_FILE,*) 'THIS PROJECT AIMS AT STUDYING THE DISCREPANCY BETWEEN LOCAL AND GLOBAL MEASUREMENTS '
    write(UNIT_EXE_FILE,*) 'OF THE HUBBLE EXPANSION RATE'

    write(UNIT_EXE_FILE,*) ' '

    write(UNIT_EXE_FILE,*) 'DATA FROM THE COSMICFLOWS-3 COMPILATION SUCCESFULLY READ '

    write(UNIT_EXE_FILE,*) ' '

    write(UNIT_EXE_FILE,*) 'N_side = ', nsmax, ' IN THE CURRENT ANALYSIS; THIS CORRESPONDS TO A NUMBER '
    write(UNIT_EXE_FILE,*) 'OF PIXELS IN THE NUMBER COUNTS MAPS EQUAL TO ', npixC

    write(UNIT_EXE_FILE,*) ' '

!################
! ANALYSIS STARTS 
!################

    If (do_galaxy_distribution_plots) then

       write(UNIT_EXE_FILE,*) 'FIGURES ANGULAR GALAXY DISTRIBUTION BEING CREATED'

       call write_python_script_angular_distribution_galaxies()

       call write_mpi_file('angular_distribution_galaxies_CF3')

       call system('cd scripts; sbatch angular_distribution_galaxies_CF3.mpi')

       write(UNIT_EXE_FILE,*) ' '

    Else

       write(UNIT_EXE_FILE,*) 'NOT DOING ANGULAR GALAXY DISTRIBUTION FIGURES'

       write(UNIT_EXE_FILE,*) ' '

    End If

    stop

    call compute_number_counts_map(1.d-2,5.d-2,.false.,i,testdipamp,test1,test2)

    print *, testdipamp, test1, test2

    close(UNIT_EXE_FILE)

End Program h0cf3




