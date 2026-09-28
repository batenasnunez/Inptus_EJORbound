module Operaciones
  use progression_model
  use Optimizacion, only: optimize_pr_bound
  use Funciones, only: full_persistence_lower, full_persistence_upper, &
    full_persistence_first_threshold, full_persistence_critical_threshold, &
    zero_persistence_lower, zero_persistence_upper, &
    zero_persistence_threshold, uniform_otgr_dplus1
  implicit none
  private

  public :: generate_pr_otgr_data
  public :: generate_full_limits
  public :: generate_zero_limits
  public :: generate_uniform_curve
  public :: generate_intermediate_limits

contains

  subroutine generate_intermediate_limits(number_of_stages,number_of_points, &
      persistence_case,lambda_homogeneous_input, &
      lambda_heterogeneous_input,lambda_first,lambda_last, &
      population_size,maximum_generations,number_of_restarts, &
      mutation_factor,crossover_rate,pr_tolerance,probability_epsilon, &
      local_tolerance,limits_file,optima_file)

    integer, intent(in) :: number_of_stages,number_of_points
    integer, intent(in) :: persistence_case
    real(kind=dp), intent(in) :: lambda_homogeneous_input
    real(kind=dp), intent(in) :: lambda_heterogeneous_input(:)
    real(kind=dp), intent(in) :: lambda_first,lambda_last
    integer, intent(in) :: population_size,maximum_generations
    integer, intent(in) :: number_of_restarts
    real(kind=dp), intent(in) :: mutation_factor,crossover_rate
    real(kind=dp), intent(in) :: pr_tolerance,probability_epsilon
    real(kind=dp), intent(in) :: local_tolerance
    character(len=*), intent(in) :: limits_file,optima_file

    integer, parameter :: limits_unit = 24
    integer, parameter :: optima_unit = 25
    integer :: i,k,limits_status,optima_status
    real(kind=dp) :: target_pr,gmin,gmax,pr_min,pr_max
    real(kind=dp), allocatable :: pmin(:),pmax(:),warm_min(:),warm_max(:)
    logical :: min_success,max_success,have_warm_min,have_warm_max

    if (number_of_stages < 1 .or. number_of_points < 2) then
      write(*,*) 'Invalid stage count or PR-grid size.'
      return
    end if

    if (population_size < 4 .or. maximum_generations < 1 .or. &
        number_of_restarts < 1) then
      write(*,*) 'Invalid Differential Evolution parameters.'
      return
    end if

    if (mutation_factor <= 0.0_dp .or. crossover_rate < 0.0_dp .or. &
        crossover_rate > 1.0_dp) then
      write(*,*) 'Invalid mutation factor or crossover rate.'
      return
    end if

    call initialize_model(number_of_stages)
    if (.not. define_persistence(persistence_case, &
        lambda_homogeneous_input,lambda_heterogeneous_input, &
        lambda_first,lambda_last)) then
      call finalize_model()
      return
    end if

    call initialize_random_generator()
    allocate(pmin(d),pmax(d),warm_min(d),warm_max(d))
    warm_min = 0.5_dp
    warm_max = 0.5_dp
    have_warm_min = .false.
    have_warm_max = .false.

    open(unit=limits_unit,file=limits_file,status='unknown', &
      action='write',iostat=limits_status)
    if (limits_status /= 0) then
      write(*,*) 'Error opening intermediate limits file: ',limits_file
      deallocate(pmin,pmax,warm_min,warm_max)
      call finalize_model()
      return
    end if

    open(unit=optima_unit,file=optima_file,status='unknown', &
      action='write',iostat=optima_status)
    if (optima_status /= 0) then
      write(*,*) 'Error opening intermediate optima file: ',optima_file
      close(limits_unit)
      deallocate(pmin,pmax,warm_min,warm_max)
      call finalize_model()
      return
    end if

    write(limits_unit,'(A)') 'PR GMIN GMAX'
    write(optima_unit,'(A)',advance='no') &
      'PR BOUND OTGR PR_CHECK RESIDUAL'
    do k = 1,d
      call write_probability_header(optima_unit,k)
    end do
    write(optima_unit,*)

    do i = 0,number_of_points-1
      target_pr = real(i,kind=dp)/real(number_of_points-1,kind=dp)

      call optimize_pr_bound(target_pr,.true.,population_size, &
        maximum_generations,number_of_restarts,mutation_factor, &
        crossover_rate,pr_tolerance,probability_epsilon,local_tolerance, &
        have_warm_min,warm_min,gmin,pmin,pr_min,min_success)

      call optimize_pr_bound(target_pr,.false.,population_size, &
        maximum_generations,number_of_restarts,mutation_factor, &
        crossover_rate,pr_tolerance,probability_epsilon,local_tolerance, &
        have_warm_max,warm_max,gmax,pmax,pr_max,max_success)

      if (min_success) then
        warm_min = pmin
        have_warm_min = .true.
      else
        gmin = -1.0_dp
        pr_min = -1.0_dp
        pmin = -1.0_dp
      end if

      if (max_success) then
        warm_max = pmax
        have_warm_max = .true.
      else
        gmax = -1.0_dp
        pr_max = -1.0_dp
        pmax = -1.0_dp
      end if

      write(limits_unit,'(ES24.16,1X,ES24.16,1X,ES24.16)') &
        target_pr,gmin,gmax
      call write_optimum_row(optima_unit,target_pr,-1,gmin,pr_min,pmin)
      call write_optimum_row(optima_unit,target_pr,1,gmax,pr_max,pmax)

      if (mod(i,max(1,(number_of_points-1)/10)) == 0) then
        write(*,*) 'Intermediate-frontier progress: ',i+1,' / ',number_of_points
      end if
    end do

    close(limits_unit)
    close(optima_unit)
    deallocate(pmin,pmax,warm_min,warm_max)
    call finalize_model()

    write(*,*) 'Intermediate-persistence optimization completed.'
    write(*,*) 'Limits file: ',limits_file
    write(*,*) 'Optima file: ',optima_file
  end subroutine generate_intermediate_limits


  subroutine write_probability_header(unit_number,index_value)
    integer, intent(in) :: unit_number,index_value
    character(len=16) :: index_text

    write(index_text,'(I10)') index_value
    write(unit_number,'(A)',advance='no') &
      ' P'//trim(adjustl(index_text))
  end subroutine write_probability_header


  subroutine write_optimum_row(unit_number,target_pr,bound_type,otgr_value, &
      calculated_pr,probabilities)

    integer, intent(in) :: unit_number,bound_type
    real(kind=dp), intent(in) :: target_pr,otgr_value,calculated_pr
    real(kind=dp), intent(in) :: probabilities(:)

    integer :: k
    real(kind=dp) :: residual

    if (calculated_pr >= 0.0_dp) then
      residual = abs(calculated_pr-target_pr)
    else
      residual = -1.0_dp
    end if

    write(unit_number,'(ES24.16,1X,I3,1X,ES24.16,1X,ES24.16,1X,ES24.16)', &
      advance='no') target_pr,bound_type,otgr_value,calculated_pr,residual
    do k = 1,size(probabilities)
      write(unit_number,'(1X,ES24.16)',advance='no') probabilities(k)
    end do
    write(unit_number,*)
  end subroutine write_optimum_row

  subroutine generate_uniform_curve(number_of_stages, number_of_points, &
      lambda_value, output_file)

    integer, intent(in) :: number_of_stages
    integer, intent(in) :: number_of_points
    real(kind=dp), intent(in) :: lambda_value
    character(len=*), intent(in) :: output_file

    integer, parameter :: output_unit = 23
    integer :: i, io_status
    real(kind=dp) :: p_value, pr_value, otgr_value

    if (number_of_stages < 1) then
      write(*,*) 'The number of stages must be positive.'
      return
    end if

    if (number_of_points < 2) then
      write(*,*) 'At least two points are required for the uniform curve.'
      return
    end if

    if (lambda_value < 0.0_dp .or. lambda_value > 1.0_dp) then
      write(*,*) 'Uniform persistence must lie in the interval [0,1].'
      return
    end if

    open(unit=output_unit, file=output_file, status='unknown', &
         action='write', iostat=io_status)

    if (io_status /= 0) then
      write(*,*) 'Error opening uniform-case output file: ', output_file
      return
    end if

    write(output_unit,'(A)') 'PR OTGR_DPLUS1'

    do i = 0, number_of_points-1
      p_value = real(i,kind=dp) / real(number_of_points-1,kind=dp)
      pr_value = p_value
      otgr_value = uniform_otgr_dplus1(p_value, lambda_value, &
        number_of_stages)
      write(output_unit,'(ES24.16,1X,ES24.16)') pr_value, otgr_value
    end do

    close(output_unit)

    write(*,*) 'Uniform-case curve completed.'
    write(*,*) 'Uniform persistence: ', lambda_value
    write(*,*) 'Output file: ', output_file
  end subroutine generate_uniform_curve

  subroutine generate_zero_limits(number_of_stages, &
      number_of_points, output_file)

    integer, intent(in) :: number_of_stages
    integer, intent(in) :: number_of_points
    character(len=*), intent(in) :: output_file

    integer, parameter :: output_unit = 22
    integer :: i, io_status
    real(kind=dp) :: a, gmin, gmax, threshold

    if (number_of_stages < 1) then
      write(*,*) 'The number of stages must be positive.'
      return
    end if

    if (number_of_points < 2) then
      write(*,*) 'At least two points are required for the limits.'
      return
    end if

    threshold = zero_persistence_threshold(number_of_stages)

    open(unit=output_unit, file=output_file, status='unknown', &
         action='write', iostat=io_status)

    if (io_status /= 0) then
      write(*,*) 'Error opening zero-persistence limits file: ', output_file
      return
    end if

    write(output_unit,'(A)') 'A GMIN GMAX'

    do i = 0, number_of_points-1
      a = real(i,kind=dp) / real(number_of_points-1,kind=dp)
      gmin = zero_persistence_lower(a, number_of_stages)
      gmax = zero_persistence_upper(a, number_of_stages)
      write(output_unit,'(ES24.16,1X,ES24.16,1X,ES24.16)') a, gmin, gmax
    end do

    close(output_unit)

    write(*,*) 'Zero-persistence limits completed.'
    write(*,*) 'Lower-frontier threshold: ', threshold
    write(*,*) 'Limits file: ', output_file
  end subroutine generate_zero_limits

  subroutine generate_full_limits(number_of_stages, &
      number_of_points, output_file)

    integer, intent(in) :: number_of_stages
    integer, intent(in) :: number_of_points
    character(len=*), intent(in) :: output_file

    integer, parameter :: output_unit = 21
    integer :: i, io_status
    real(kind=dp) :: a, gmin, gmax, a1, ac

    if (number_of_stages < 1) then
      write(*,*) 'The number of stages must be positive.'
      return
    end if

    if (number_of_points < 2) then
      write(*,*) 'At least two points are required for the limits.'
      return
    end if

    a1 = full_persistence_first_threshold(number_of_stages)
    ac = full_persistence_critical_threshold(number_of_stages)

    if (a1 < 0.0_dp .or. ac < 0.0_dp) then
      write(*,*) 'The full-persistence thresholds could not be calculated.'
      return
    end if

    open(unit=output_unit, file=output_file, status='unknown', &
         action='write', iostat=io_status)

    if (io_status /= 0) then
      write(*,*) 'Error opening limits file: ', output_file
      return
    end if

    write(output_unit,'(A)') 'A GMIN GMAX'

    do i = 0, number_of_points-1
      a = real(i,kind=dp) / real(number_of_points-1,kind=dp)
      gmin = full_persistence_lower(a, number_of_stages)
      gmax = full_persistence_upper(a, number_of_stages)
      write(output_unit,'(ES24.16,1X,ES24.16,1X,ES24.16)') a, gmin, gmax
    end do

    close(output_unit)

    write(*,*) 'Full-persistence limits completed.'
    write(*,*) 'First threshold a1: ', a1
    write(*,*) 'Critical threshold ac: ', ac
    write(*,*) 'Limits file: ', output_file
  end subroutine generate_full_limits

  subroutine generate_pr_otgr_data(number_of_stages, number_of_samples, &
      persistence_case, lambda_homogeneous_input, &
      lambda_heterogeneous_input, lambda_first, lambda_last, &
      output_file)

    integer, intent(in) :: number_of_stages
    integer, intent(in) :: number_of_samples
    integer, intent(in) :: persistence_case
    real(kind=dp), intent(in) :: lambda_homogeneous_input
    real(kind=dp), intent(in) :: lambda_heterogeneous_input(:)
    real(kind=dp), intent(in) :: lambda_first, lambda_last
    character(len=*), intent(in) :: output_file

    integer, parameter :: output_unit = 20
    integer :: i, io_status
    real(kind=dp) :: pr_value, otgr_dplus1_value
    logical :: updated

    if (number_of_stages < 1) then
      write(*,*) 'The number of stages must be positive.'
      return
    end if

    if (number_of_samples < 1) then
      write(*,*) 'The number of samples must be positive.'
      return
    end if

    call initialize_model(number_of_stages)

    if (.not. define_persistence(persistence_case, &
        lambda_homogeneous_input, lambda_heterogeneous_input, &
        lambda_first, lambda_last)) then
      call finalize_model()
      return
    end if

    call initialize_random_generator()

    open(unit=output_unit, file=output_file, status='unknown', &
         action='write', iostat=io_status)

    if (io_status /= 0) then
      write(*,*) 'Error opening output file: ', output_file
      call finalize_model()
      return
    end if

    ! Plain column names facilitate direct import into OriginPro.
    write(output_unit,'(A)') 'PR OTGR_DPLUS1'

    do i = 1, number_of_samples
      call random_number(p)

      updated = update_model()
      if (.not. updated) then
        write(*,*) 'Invalid model parameters in sample ', i
        close(output_unit)
        call finalize_model()
        return
      end if

      pr_value = performance_rate()
      otgr_dplus1_value = graduation_rate_dplus1()

      write(output_unit,'(ES24.16,1X,ES24.16)') &
        pr_value, otgr_dplus1_value
    end do

    close(output_unit)
    call finalize_model()

    write(*,*) 'Simulation completed.'
    write(*,*) 'Number of samples: ', number_of_samples
    write(*,*) 'Output file: ', output_file
  end subroutine generate_pr_otgr_data


  logical function define_persistence(persistence_case, &
      lambda_homogeneous_input, lambda_heterogeneous_input, &
      lambda_first, lambda_last)

    integer, intent(in) :: persistence_case
    real(kind=dp), intent(in) :: lambda_homogeneous_input
    real(kind=dp), intent(in) :: lambda_heterogeneous_input(:)
    real(kind=dp), intent(in) :: lambda_first, lambda_last

    integer :: k
    real(kind=dp) :: fraction

    define_persistence = .false.

    select case (persistence_case)

    case (1)
      call set_homogeneous_persistence(lambda_homogeneous_input)

    case (2)
      if (size(lambda_heterogeneous_input) /= d) then
        write(*,*) 'The heterogeneous lambda array must have d elements.'
        return
      end if
      lambda = lambda_heterogeneous_input

    case (3)
      if (d == 1) then
        lambda(1) = lambda_first
      else
        do k = 1, d
          fraction = real(k-1, kind=dp) / real(d-1, kind=dp)
          lambda(k) = lambda_first + &
            fraction * (lambda_last - lambda_first)
        end do
      end if

    case default
      write(*,*) 'Invalid persistence_case. Use 1, 2, or 3.'
      return

    end select

    if (any(lambda < 0.0_dp) .or. any(lambda > 1.0_dp)) then
      write(*,*) 'Persistence values must lie in the interval [0,1].'
      return
    end if

    define_persistence = .true.
  end function define_persistence


  subroutine initialize_random_generator()
    integer, parameter :: seed_modulus = 1000000007
    integer :: seed_size, j, clock_count
    integer :: date_values(8)
    integer, allocatable :: seed(:)

    call system_clock(count=clock_count)
    call date_and_time(values=date_values)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))

    do j = 1, seed_size
      seed(j) = modulo(modulo(clock_count, seed_modulus) + 104729*j + &
        1009*date_values(8) + 101*date_values(7) + &
        17*date_values(6) + date_values(5), seed_modulus)
      if (seed(j) == 0) seed(j) = j
    end do

    call random_seed(put=seed)
    deallocate(seed)
  end subroutine initialize_random_generator

end module Operaciones
