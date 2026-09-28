module Funciones
  use progression_model, only: dp
  implicit none
  private

  real(kind=dp), parameter :: invalid_value = -1.0_dp
  real(kind=dp), parameter :: root_tolerance = 1.0e-12_dp

  integer, save :: cached_d = -1
  real(kind=dp), save :: cached_ac = invalid_value

  public :: full_persistence_lower
  public :: full_persistence_upper
  public :: full_persistence_corner
  public :: full_persistence_interior
  public :: full_persistence_homogeneous
  public :: full_persistence_first_threshold
  public :: full_persistence_critical_threshold
  public :: zero_persistence_lower
  public :: zero_persistence_upper
  public :: zero_persistence_threshold
  public :: uniform_otgr_dplus1

contains

  real(kind=dp) function uniform_otgr_dplus1(p_value, lambda_value, stages)
    real(kind=dp), intent(in) :: p_value, lambda_value
    integer, intent(in) :: stages

    uniform_otgr_dplus1 = invalid_value

    if (stages < 1) return
    if (p_value < 0.0_dp .or. p_value > 1.0_dp) return
    if (lambda_value < 0.0_dp .or. lambda_value > 1.0_dp) return

    uniform_otgr_dplus1 = p_value**stages * &
      (1.0_dp + real(stages,kind=dp) * &
      lambda_value * (1.0_dp-p_value))
  end function uniform_otgr_dplus1

  real(kind=dp) function full_persistence_homogeneous(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    full_persistence_homogeneous = invalid_value
    if (.not. valid_arguments(a, stages)) return

    full_persistence_homogeneous = a**stages * &
      (1.0_dp + real(stages,kind=dp) * (1.0_dp-a))
  end function full_persistence_homogeneous


  real(kind=dp) function full_persistence_lower(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    full_persistence_lower = full_persistence_homogeneous(a, stages)
  end function full_persistence_lower


  real(kind=dp) function full_persistence_corner(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    real(kind=dp) :: r, denominator

    full_persistence_corner = invalid_value
    if (.not. valid_arguments(a, stages)) return

    if (stages == 1) then
      full_persistence_corner = a * (2.0_dp-a)
      return
    end if

    denominator = real(stages,kind=dp) - &
      real(stages-1,kind=dp) * a
    r = a / denominator
    full_persistence_corner = r * (2.0_dp-r)
  end function full_persistence_corner


  real(kind=dp) function full_persistence_interior(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    real(kind=dp) :: delta_squared, delta, u, v
    real(kind=dp) :: dreal

    full_persistence_interior = invalid_value
    if (.not. valid_arguments(a, stages)) return
    if (stages < 2) return

    dreal = real(stages,kind=dp)
    delta_squared = (dreal+1.0_dp)**2 * (dreal-2.0_dp)**2 * a**2 &
      - 2.0_dp*dreal*(dreal+1.0_dp) * &
        (dreal**2-dreal+2.0_dp) * a &
      + dreal**2 * (dreal+1.0_dp)**2

    if (delta_squared < -root_tolerance) return
    delta = sqrt(max(0.0_dp,delta_squared))

    v = ((dreal+1.0_dp) * (a*dreal-2.0_dp*a+dreal) + delta) / &
      (2.0_dp*dreal**2)
    u = (1.0_dp+dreal-dreal*v) / 2.0_dp

    if (u < 0.0_dp .or. u > 1.0_dp .or. &
        v < 0.0_dp .or. v > 1.0_dp) return

    full_persistence_interior = u * v**(stages-1) * &
      (1.0_dp+dreal-u-real(stages-1,kind=dp)*v)
  end function full_persistence_interior


  real(kind=dp) function full_persistence_first_threshold(stages)
    integer, intent(in) :: stages

    full_persistence_first_threshold = invalid_value
    if (stages < 1) return

    if (stages == 1) then
      full_persistence_first_threshold = 1.0_dp
    else
      full_persistence_first_threshold = real(stages,kind=dp) / &
        real(stages+1,kind=dp)
    end if
  end function full_persistence_first_threshold


  real(kind=dp) function full_persistence_critical_threshold(stages)
    integer, intent(in) :: stages

    integer, parameter :: scan_points = 10000
    integer, parameter :: maximum_iterations = 200
    integer :: i, iteration
    real(kind=dp) :: a1, left, right, middle
    real(kind=dp) :: f_left, f_right, f_middle

    full_persistence_critical_threshold = invalid_value
    if (stages < 1) return

    if (stages == 1) then
      full_persistence_critical_threshold = 1.0_dp
      return
    end if

    if (stages == cached_d) then
      full_persistence_critical_threshold = cached_ac
      return
    end if

    a1 = full_persistence_first_threshold(stages)
    left = a1
    f_left = crossing_function(left, stages)
    right = left
    f_right = f_left

    do i = 1, scan_points
      right = a1 + (1.0_dp-a1) * &
        real(i,kind=dp) / real(scan_points,kind=dp)
      f_right = crossing_function(right, stages)

      if (f_left > 0.0_dp .and. f_right <= 0.0_dp) exit

      left = right
      f_left = f_right
    end do

    if (i > scan_points) return

    do iteration = 1, maximum_iterations
      middle = 0.5_dp * (left+right)
      f_middle = crossing_function(middle, stages)

      if (abs(f_middle) <= root_tolerance) then
        left = middle
        right = middle
        exit
      end if

      if (abs(right-left) <= root_tolerance) exit

      if (f_middle > 0.0_dp) then
        left = middle
      else
        right = middle
      end if
    end do

    cached_d = stages
    cached_ac = 0.5_dp * (left+right)
    full_persistence_critical_threshold = cached_ac
  end function full_persistence_critical_threshold


  real(kind=dp) function full_persistence_upper(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    real(kind=dp) :: a1, ac

    full_persistence_upper = invalid_value
    if (.not. valid_arguments(a, stages)) return

    if (stages == 1) then
      full_persistence_upper = a * (2.0_dp-a)
      return
    end if

    a1 = full_persistence_first_threshold(stages)
    ac = full_persistence_critical_threshold(stages)
    if (ac < 0.0_dp) return

    if (a <= a1) then
      full_persistence_upper = full_persistence_corner(a, stages)
    else if (a < ac) then
      full_persistence_upper = full_persistence_interior(a, stages)
    else
      full_persistence_upper = full_persistence_homogeneous(a, stages)
    end if
  end function full_persistence_upper


  real(kind=dp) function zero_persistence_threshold(stages)
    integer, intent(in) :: stages

    zero_persistence_threshold = invalid_value
    if (stages < 1) return

    zero_persistence_threshold = real(stages-1,kind=dp) / &
      real(stages,kind=dp)
  end function zero_persistence_threshold


  real(kind=dp) function zero_persistence_lower(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    real(kind=dp) :: threshold

    zero_persistence_lower = invalid_value
    if (.not. valid_arguments(a, stages)) return

    threshold = zero_persistence_threshold(stages)

    if (a <= threshold) then
      zero_persistence_lower = 0.0_dp
    else
      zero_persistence_lower = real(stages,kind=dp)*a - &
        real(stages-1,kind=dp)
    end if
  end function zero_persistence_lower


  real(kind=dp) function zero_persistence_upper(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    real(kind=dp) :: denominator

    zero_persistence_upper = invalid_value
    if (.not. valid_arguments(a, stages)) return

    denominator = real(stages,kind=dp) - &
      real(stages-1,kind=dp)*a
    zero_persistence_upper = a / denominator
  end function zero_persistence_upper


  logical function valid_arguments(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    valid_arguments = stages >= 1 .and. a >= 0.0_dp .and. a <= 1.0_dp
  end function valid_arguments


  real(kind=dp) function crossing_function(a, stages)
    real(kind=dp), intent(in) :: a
    integer, intent(in) :: stages

    crossing_function = full_persistence_interior(a, stages) - &
      full_persistence_homogeneous(a, stages)
  end function crossing_function

end module Funciones
