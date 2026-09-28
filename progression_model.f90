module progression_model
  implicit none

  ! Double precision requested for the numerical implementation.
  integer, parameter, public :: dp = 8
  real(kind=dp), parameter :: invalid_value = -1.0_dp
  real(kind=dp), parameter :: tolerance = 100.0_dp * epsilon(1.0_dp)

  ! Global model variables. They are allocated by initialize_model.
  integer, public :: d = 0
  real(kind=dp), allocatable, public :: p(:)
  real(kind=dp), allocatable, public :: lambda(:)

  ! Internally calculated stage quantities.
  real(kind=dp), allocatable, private :: q(:)
  real(kind=dp), allocatable, private :: stage_pass(:)
  real(kind=dp), allocatable, private :: survival(:)
  real(kind=dp), allocatable, private :: stage_dropout(:)
  real(kind=dp), allocatable, private :: stage_duration(:)

  logical, private :: model_initialized = .false.
  logical, private :: parameters_valid = .false.

  public :: initialize_model
  public :: finalize_model
  public :: update_model
  public :: validate_parameters
  public :: set_homogeneous_persistence
  public :: P_k
  public :: Pi_k
  public :: Q_k
  public :: stage_time
  public :: performance_rate
  public :: graduation_rate_d
  public :: graduation_rate_dplus1
  public :: total_graduation_probability
  public :: dropout_probability
  public :: normalized_graduation_time
  public :: mean_dropout_time

contains

  subroutine initialize_model(number_of_stages)
    integer, intent(in) :: number_of_stages

    call finalize_model()

    if (number_of_stages < 1) then
      d = 0
      return
    end if

    d = number_of_stages
    allocate(p(1:d), lambda(1:d))
    allocate(q(1:d), stage_pass(1:d), stage_dropout(1:d))
    allocate(stage_duration(1:d), survival(0:d))

    p = 0.0_dp
    lambda = 0.0_dp
    q = 0.0_dp
    stage_pass = 0.0_dp
    stage_dropout = 0.0_dp
    stage_duration = 0.0_dp
    survival = 0.0_dp
    survival(0) = 1.0_dp

    model_initialized = .true.
    parameters_valid = .false.
  end subroutine initialize_model


  subroutine finalize_model()
    if (allocated(p)) deallocate(p)
    if (allocated(lambda)) deallocate(lambda)
    if (allocated(q)) deallocate(q)
    if (allocated(stage_pass)) deallocate(stage_pass)
    if (allocated(stage_dropout)) deallocate(stage_dropout)
    if (allocated(stage_duration)) deallocate(stage_duration)
    if (allocated(survival)) deallocate(survival)

    d = 0
    model_initialized = .false.
    parameters_valid = .false.
  end subroutine finalize_model


  subroutine set_homogeneous_persistence(lambda_value)
    real(kind=dp), intent(in) :: lambda_value

    if (.not. model_initialized) return
    lambda = lambda_value
    parameters_valid = .false.
  end subroutine set_homogeneous_persistence


  logical function validate_parameters()
    if (.not. model_initialized) then
      validate_parameters = .false.
      return
    end if

    validate_parameters = all(p >= 0.0_dp .and. p <= 1.0_dp) .and. &
                          all(lambda >= 0.0_dp .and. lambda <= 1.0_dp)
  end function validate_parameters


  logical function update_model()
    integer :: k
    real(kind=dp) :: denominator

    parameters_valid = validate_parameters()
    update_model = parameters_valid
    if (.not. parameters_valid) return

    q = 1.0_dp - p
    survival(0) = 1.0_dp

    do k = 1, d
      denominator = 1.0_dp - lambda(k) * q(k)

      if (denominator > tolerance) then
        stage_pass(k) = p(k) / denominator
        stage_dropout(k) = q(k) * (1.0_dp - lambda(k)) / denominator
        stage_duration(k) = 1.0_dp / denominator
      else
        ! Singular trap: p_k=0 and lambda_k=1. The student neither
        ! advances nor withdraws, and the residence time is infinite.
        stage_pass(k) = 0.0_dp
        stage_dropout(k) = 0.0_dp
        stage_duration(k) = huge(1.0_dp)
      end if

      survival(k) = survival(k-1) * stage_pass(k)
    end do
  end function update_model


  real(kind=dp) function P_k(k)
    integer, intent(in) :: k
    logical :: updated

    P_k = invalid_value
    if (k < 1 .or. k > d) return
    updated = update_model()
    if (.not. updated) return
    P_k = stage_pass(k)
  end function P_k


  real(kind=dp) function Pi_k(k)
    integer, intent(in) :: k
    logical :: updated

    Pi_k = invalid_value
    if (k < 0 .or. k > d) return
    updated = update_model()
    if (.not. updated) return
    Pi_k = survival(k)
  end function Pi_k


  real(kind=dp) function Q_k(k)
    integer, intent(in) :: k
    logical :: updated

    Q_k = invalid_value
    if (k < 1 .or. k > d) return
    updated = update_model()
    if (.not. updated) return
    Q_k = stage_dropout(k)
  end function Q_k


  real(kind=dp) function stage_time(k)
    integer, intent(in) :: k
    logical :: updated

    stage_time = invalid_value
    if (k < 1 .or. k > d) return
    updated = update_model()
    if (.not. updated) return
    stage_time = stage_duration(k)
  end function stage_time


  real(kind=dp) function performance_rate()
    integer :: k
    real(kind=dp) :: numerator, denominator, enrolled_stage
    logical :: updated

    performance_rate = invalid_value
    updated = update_model()
    if (.not. updated) return

    numerator = sum(survival(1:d))
    denominator = 0.0_dp

    do k = 1, d
      if (stage_duration(k) >= huge(1.0_dp)) then
        if (survival(k-1) > tolerance) then
          performance_rate = 0.0_dp
          return
        end if
      else
        ! Pi_k/p_k = Pi_(k-1)/(1-lambda_k*q_k). This form remains
        ! regular when p_k=0, provided the stage is not a singular trap.
        enrolled_stage = survival(k-1) * stage_duration(k)
        denominator = denominator + enrolled_stage
      end if
    end do

    if (denominator > tolerance) performance_rate = numerator / denominator
  end function performance_rate


  real(kind=dp) function graduation_rate_d()
    logical :: updated

    graduation_rate_d = invalid_value
    updated = update_model()
    if (.not. updated) return
    graduation_rate_d = product(p)
  end function graduation_rate_d


  real(kind=dp) function graduation_rate_dplus1()
    logical :: updated

    graduation_rate_dplus1 = invalid_value
    updated = update_model()
    if (.not. updated) return
    graduation_rate_dplus1 = product(p) * &
      (1.0_dp + sum(lambda * q))
  end function graduation_rate_dplus1


  real(kind=dp) function total_graduation_probability()
    logical :: updated

    total_graduation_probability = invalid_value
    updated = update_model()
    if (.not. updated) return
    total_graduation_probability = survival(d)
  end function total_graduation_probability


  real(kind=dp) function dropout_probability()
    integer :: k
    logical :: updated

    dropout_probability = invalid_value
    updated = update_model()
    if (.not. updated) return

    dropout_probability = 0.0_dp
    do k = 1, d
      dropout_probability = dropout_probability + &
        survival(k-1) * stage_dropout(k)
    end do
  end function dropout_probability


  real(kind=dp) function normalized_graduation_time()
    integer :: k
    real(kind=dp) :: total_time
    logical :: updated

    normalized_graduation_time = invalid_value
    updated = update_model()
    if (.not. updated) return

    total_time = 0.0_dp
    do k = 1, d
      if (stage_duration(k) >= huge(1.0_dp)) then
        normalized_graduation_time = 0.0_dp
        return
      end if
      total_time = total_time + stage_duration(k)
    end do

    if (total_time > tolerance) then
      normalized_graduation_time = real(d, kind=dp) / total_time
    end if
  end function normalized_graduation_time


  real(kind=dp) function mean_dropout_time()
    integer :: k
    real(kind=dp) :: dropout_total, cumulative_time, weighted_time
    real(kind=dp) :: dropout_at_stage
    logical :: updated

    mean_dropout_time = invalid_value
    updated = update_model()
    if (.not. updated) return

    dropout_total = 0.0_dp
    cumulative_time = 0.0_dp
    weighted_time = 0.0_dp

    do k = 1, d
      if (stage_duration(k) >= huge(1.0_dp)) then
        cumulative_time = huge(1.0_dp)
      else if (cumulative_time < huge(1.0_dp)) then
        cumulative_time = cumulative_time + stage_duration(k)
      end if

      dropout_at_stage = survival(k-1) * stage_dropout(k)
      dropout_total = dropout_total + dropout_at_stage

      ! Avoid the indeterminate numerical product 0*HUGE at an
      ! unreachable stage located after a singular trap.
      if (dropout_at_stage > tolerance) then
        weighted_time = weighted_time + dropout_at_stage * cumulative_time
      end if
    end do

    ! TTA is conditional on dropout and is undefined when dropout is zero.
    if (dropout_total > tolerance) then
      mean_dropout_time = weighted_time / dropout_total
    end if
  end function mean_dropout_time

end module progression_model
