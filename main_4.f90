program main_progression
  use progression_model, only: dp
  use Operaciones, only: generate_pr_otgr_data, &
    generate_full_limits, generate_zero_limits, generate_uniform_curve, &
    generate_intermediate_limits
  implicit none

  !=====================================================================
  ! INPUT PARAMETERS
  ! Modify only this block to define a simulation.
  !=====================================================================

  ! Nominal number of formative stages.
  integer, parameter :: d_input = 4

  ! Operation to perform:
  !   1 = random sampling of PR and OTGR(d+1)
  !   2 = analytical limits for full persistence, lambda = 1 (falsa)
  !   3 = analytical limits for zero persistence, lambda = 0
  !   4 = homogeneous p and homogeneous lambda curve
  !   5 = numerical limits for a fixed intermediate persistence profile
  integer, parameter :: operation_case = 5

  ! Number of random configurations of the stage pass probabilities.
  integer, parameter :: number_of_samples = 10000

  ! Persistence model:
  !   1 = homogeneous persistence
  !   2 = heterogeneous persistence specified stage by stage
  !   3 = persistence increasing linearly from lambda_1 to lambda_2
  integer, parameter :: persistence_case = 1

  ! Input used only when persistence_case = 1.
  real(kind=dp), parameter :: lambda_homogeneous = 0.10_dp

  ! Input used only when persistence_case = 2.
  ! The number of elements must equal d_input.
  real(kind=dp), parameter :: lambda_heterogeneous(d_input) = &
    (/ 0.70_dp, 0.80_dp, 0.90_dp , 0.95_dp  /)!, 1.00_dp, 1.00_dp, 1.00_dp  /)!

  ! Inputs used only when persistence_case = 3.
  ! lambda_1 applies to the first stage and lambda_2 to the last stage.
  real(kind=dp), parameter :: lambda_1 = 0.00_dp
  real(kind=dp), parameter :: lambda_2 = 1.00_dp

  ! Name of the random-sampling output file.
  character(len=*), parameter :: output_file = 'PR_OTGR_data.dat'

  ! Number of points used to draw the analytical limit curves.
  integer, parameter :: number_of_curve_points = 1001

  ! File containing the analytical full-persistence limit curves.
  character(len=*), parameter :: limits_output_file = &
    'full_persistence_limits.dat'

  ! File containing the analytical zero-persistence limit curves.
  character(len=*), parameter :: zero_limits_output_file = &
    'zero_persistence_limits.dat'

  ! File containing the homogeneous-p and homogeneous-lambda curve.
  character(len=*), parameter :: uniform_output_file = &
    'uniform_case.dat'

  ! Numerical optimization grid and Differential Evolution parameters.
  integer, parameter :: optimization_pr_points = 201
  integer, parameter :: de_population_size = 200
  integer, parameter :: de_max_generations = 1000
  integer, parameter :: de_restarts = 5
  real(kind=dp), parameter :: de_mutation_factor = 0.80_dp
  real(kind=dp), parameter :: de_crossover_rate = 0.90_dp
  real(kind=dp), parameter :: optimization_pr_tolerance = 1.0e-10_dp
  real(kind=dp), parameter :: probability_epsilon = 1.0e-10_dp
  real(kind=dp), parameter :: local_search_tolerance = 1.0e-8_dp

  ! Files containing the numerical limits and their optimal configurations.
  character(len=*), parameter :: intermediate_limits_file = &
    'intermediate_persistence_limits.dat'
  character(len=*), parameter :: intermediate_optima_file = &
    'intermediate_persistence_optima.dat'

  !=====================================================================
  ! END OF INPUT PARAMETERS
  !=====================================================================

  select case (operation_case)

  case (1)

    call generate_pr_otgr_data(d_input, number_of_samples, &
      persistence_case, lambda_homogeneous, lambda_heterogeneous, &
      lambda_1, lambda_2, output_file)

  case (2)

    call generate_full_limits(d_input, &
      number_of_curve_points, limits_output_file)

  case (3)

    call generate_zero_limits(d_input, &
      number_of_curve_points, zero_limits_output_file)

  case (4)

    call generate_uniform_curve(d_input, number_of_curve_points, &
      lambda_homogeneous, uniform_output_file)

  case (5)

    call generate_intermediate_limits(d_input,optimization_pr_points, &
      persistence_case,lambda_homogeneous,lambda_heterogeneous, &
      lambda_1,lambda_2,de_population_size,de_max_generations, &
      de_restarts,de_mutation_factor,de_crossover_rate, &
      optimization_pr_tolerance,probability_epsilon, &
      local_search_tolerance,intermediate_limits_file, &
      intermediate_optima_file)

  case default

    write(*,*) 'Invalid operation_case. Use 1, 2, 3, 4, or 5.'

  end select

end program main_progression
