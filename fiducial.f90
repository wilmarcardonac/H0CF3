Module fiducial

  use healpix_types

  Implicit none

  save 

  !#################################
  ! PARAMETERS OF NUMBER COUNTS MAPS
  !#################################
  
  Real*8,parameter    :: miss = -1.63750d-13
  Integer(kind=I4B),parameter :: nsmax = 1  ! Nside FOR NUMBER COUNTS MAP
  Integer*4,parameter :: nlmax = 2*nsmax    ! HIGHEST MULTIPOLE
  Integer*4,parameter :: UNIT_EXE_FILE = 90           ! UNIT NUMBER FOR EXECUTION INFORMATION FILE
  Integer*4,parameter :: UNIT_ANAFAST_PAR_FILE = 92   ! UNIT NUMBER FOR FILE
  Integer(kind=I4B), parameter :: RING_ORDERING = 1 
  Integer(kind=I4B), parameter :: DEGREE_REMOVE_DIPOLE = 2
  Real*8, parameter :: redshift_min = 1.d-2 ! SET LOWEST REDSHIFT FOR ANALYSIS 
  Real*8, parameter :: redshift_max = 5.d-2 ! SET HIGHEST REDSHIFT FOR ANALYSIS 
  Real*8, parameter :: redshift_step = 5.d-3 ! SET STEP IN REDSHIFT FOR ANALYSIS  
  Character(len=*),parameter :: ORDERING_NC_MAPS = 'RING'!'NESTED'    ! ORDERING NUMBER COUNTS MAPS
  Character(len=*),parameter :: SYS_COORD = 'G' ! MAP COORDINATE SYSTEM 



  !################################
  ! VARIABLES IN NUMBER COUNTS MAPS 
  !################################

  Integer(kind=I8B) :: npixC                ! NUMBER OF PIXELS IN NUMBER COUNTS MAPS

  !################################
  ! SPECIFICATIONS FOR THE ANALYSIS
  !################################

  Logical,parameter   :: do_galaxy_distribution_plots = .false. ! DO PLOTS OF ANGULAR DISTRIBUTION OF GALAXIES IF SET IT TRUE
  Logical,parameter   :: do_jackknife_analysis = .false.        ! DO JACKKNIFE ANALYSIS IF SET IT TRUE

  Character(len=*),parameter :: beam_file = " '' "    ! PATH TO BEAM FILE
  Character(len=*),parameter :: almsfile = " '' "     ! PATH ALMS FILE
  Character(len=*),parameter :: plmfile = " '' "      ! PATH TO PLM FILE
  Character(len=*),parameter :: outfile_alms = " '' "    ! PATH TO ALMS OUTPUT FILE 
!  Character(len=*),parameter :: PATH_TO_ANAFAST_PARAMETER_FILE = './anafast_parameter_files/g_cmb_' ! PATH TO SYNFAST PARAMETERS FILES
!  Character(len=*),parameter :: PATH_TO_CMB_MAPS = './cmb_maps/map_g_' ! PATH TO SIMULATED CMB MAPS
  Character(len=*),parameter :: PATH_TO_HEALPIX_DATA = '/home/wcardona/software/Healpix/Healpix_3.31/data' ! PATH TO HEALPIX DATA

  Character(len=*),parameter :: EXECUTION_INFORMATION = './output/execution_information.txt' ! PATH TO EXECUTION INFORMATION FILE
  character(len=*),parameter :: fmt = '(I4.4)'  ! FORMAT NUMBERS 1 - 1000
!  Character(len=*),parameter :: PATH_TO_CMB_MASK = './data/COM_CMB_IQU-common-field-MaskInt_2048_R2.01.fits' ! PLANCK MASK TO BE USED (NESTED ORDERING)
!  Character(len=*),parameter :: PATH_TO_PLANCK_CMB_MAP = './data/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits' ! PLANCK CMB MAP TO BE USED (NESTED ORDERING)
!  Character(len=*),parameter :: PATH_TO_PLANCK_CMB_MAP = './data/COM_CompMap_CMB-smica_2048_R1.20.fits' ! PLANCK CMB MAP TO BE USED (NESTED ORDERING). 2013 release including inpainted SMICA
!  Character(len=*),parameter :: PATH_TO_VSK_MASK = './vsk_maps/vsk_mask.fits' ! VSK MASK TO BE USED (RING ORDERING)
!  Character(len=*),parameter :: PATH_TO_VSK_SPECTRA = './vsk_angular_power_spectrum/' 
!  Character(len=*),parameter :: PATH_TO_CMB_FREQUENCY_MAPS = './cmb_maps/frequency-maps/' ! PATH TO CMB FREQUENCY MAPS (DATA AND FFP8.1 SIMULATIONS)

  !################
  ! MCMC PARAMETERS
  !################


End Module fiducial
