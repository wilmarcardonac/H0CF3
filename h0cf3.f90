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
  use omp_lib

!#################################
! DECLARE VARIABLES AND PARAMETERS
!#################################

    Implicit none

    Integer*8 :: i
    Integer*4 :: index_options
    
    Real*8 :: testdipamp,test1,test2,wtime

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

       Do index_options=1,3

          call write_python_script_angular_distribution_galaxies(index_options)

       End Do

       call write_mpi_file('ang_dist_gal_whole_CF3_dataset',1)

       call write_mpi_file('ang_dist_gal_redshift_bins_CF3_dataset',2)

       call write_mpi_file('ang_dist_gal_redshift_bins_cumulative_CF3_dataset',3)

       write(UNIT_EXE_FILE,*) ' '

       call system('cd scripts; sbatch ang_dist_gal_whole_CF3_dataset.mpi')

       write(UNIT_EXE_FILE,*) ' '

       call system('cd scripts; sbatch ang_dist_gal_redshift_bins_CF3_dataset.mpi')

       write(UNIT_EXE_FILE,*) ' '

       call system('cd scripts; sbatch ang_dist_gal_redshift_bins_cumulative_CF3_dataset.mpi')

       write(UNIT_EXE_FILE,*) ' '

    Else

       write(UNIT_EXE_FILE,*) 'NOT DOING ANGULAR GALAXY DISTRIBUTION FIGURES'

       write(UNIT_EXE_FILE,*) ' '

    End If

    write(UNIT_EXE_FILE,*) 'COMPUTING NUMBER COUNTS MAPS FOR DIFFERENT REDSHIFT BINS'

    write(UNIT_EXE_FILE,*) ' '

    write(UNIT_EXE_FILE,*) 'AND DOING CORRESPONDING JACKKNIFE ANALYSIS'

    write(UNIT_EXE_FILE,*) ' '

    wtime = omp_get_wtime() ! SETTING STARTING TIME

    !$omp Parallel Do
    Do index_options=1,number_redshift_bins

       call compute_number_counts_map(redshift_min,redshift_min+index_options*redshift_step,.false.,&
            i,testdipamp,test1,test2)

       call jackknife_analysis(redshift_min,redshift_min+index_options*redshift_step)

    End Do
    !$omp End Parallel Do 

    write(UNIT_EXE_FILE,*) 'THE ANALYSIS WAS PERFORMED IN ',(omp_get_wtime()-wtime)/3.6d3,' HOURS'

    close(UNIT_EXE_FILE)

End Program h0cf3




