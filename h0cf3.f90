Program h0cf3 

!####################
! LOAD NEEDED MODULES 
!####################

    use fiducial
    use arrays
    use functions 

!#################################
! DECLARE VARIABLES AND PARAMETERS
!#################################

    Implicit none

!!$    Integer*4 :: m,n,i                                      ! INTERGER FOR SHORT LOOPS 
!!$    Integer*4 :: seed1,seed2                                      ! SEEDS FOR RANDOM NUMBER GENERATOR 
!!$    Integer*4 :: number_accepted_points,number_rejected_points    ! MCMC PARAMETERS
!!$    Integer*4 :: weight                                           ! IT COUNTS THE NUMBER OF STEPS TAKEN BEFORE MOVING TO A NEW POINT IN MCMC 
!!$
!!$    Real*4 :: average_acceptance_probability           ! SAVES ACCEPTANCE PROBABILITY 
!!$    Real*4 :: genunf                            ! RANDOM UNIFOR DEVIATES 
!!$
!!$    Real*8 :: random_uniform                           ! SAVES RANDOM UNIFORM DEVIATE BETWEEN 0 AND 1 
!!$    Real*8 :: old_loglikelihood,current_loglikelihood  ! TEMPORALY SAVES LIKELIHOOD VALUES 
!!$!    Real*8 :: chi2SNIabestfit                          ! SAVES CHI2 OF SNIA AT THE BEST FIT
!!$    Real*4,dimension(number_of_parameters*(number_of_parameters+3)/2 + 1) :: parm     ! ARRAY NEEDED BY RANDON NUMBER GENERATOR 
!!$    Real*4,dimension(number_of_parameters) :: work,x_old,x_new                        ! ARRAYS NEEDED BY RANDOM NUMBER GENERATOR 
!!$    Real*8,dimension(number_of_parameters) :: bestfit,means                           ! SAVES BESTFIT AND MEANS OF PARAMETERS 
!!$    Real*8,dimension(number_of_hosts_galaxies) :: sigmainthost
!!$
!!$    Logical :: not_good_aap,non_plausible_parameters   ! IT CONTROLS PLAUSIBLE VALUES OF COSMOLOGICAL PARAMETERS
!!$    Logical,dimension(number_of_parameters) :: plausibility  
!!$
!!$    Character(len=10) :: string ! STORES STRINGS FOR INTEGERS
!!$    Character(len=12),dimension(number_hyperparameters) :: alpha_string
!!$    Character(len=30),dimension(number_of_parameters) :: paramnames,latexname
!!$    Character(len=5) :: galaxy

!##########################################################
! ASSIGNMENTS AND INITIALIZATION OF RANDOM NUMBER GENERATOR
!##########################################################

!!$    galaxy = host(8)
!!$
!!$    sigmainthost(:) = 0.d0
!!$
!!$    weight = 1
!!$
!!$    number_rejected_points = 0
!!$
!!$    number_accepted_points = 0
!!$
!!$    call initialize()                               ! INITIALIZE RANDOM NUMBER GENERATOR
!!$
!!$    call phrtsd(phrase,seed1,seed2)                 ! GENERATE SEEDS FOR RANDOM NUMBERS FROM PHRASE
!!$
!!$    call set_initial_seed(seed1,seed2)              ! SET INITIAL SEEDS FOR RANDOM NUMBER GENERATOR
!!$
!!$    ! ALLOCATING MEMORY FOR POINTS IN PARAMETER SPACE AND ACCEPTANCE PROBABILITY
!!$    allocate (old_point(1:number_of_parameters),current_point(1:number_of_parameters),&
!!$    acceptance_probability(number_iterations),Covgauss(number_of_parameters,number_of_parameters),&
!!$    Covguess(number_of_parameters,number_of_parameters),stat = status1)
!!$
!!$    call set_covariance_matrix()
!!$
!!$    call read_data()

!##################################
! MARKOV CHAIN MONTE CARLO ANALYSIS
!##################################

 

End Program h0cf3




