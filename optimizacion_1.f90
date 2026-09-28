module Optimizacion
  use progression_model
  implicit none
  private

  public :: optimize_pr_bound

contains

  subroutine optimize_pr_bound(target_pr, minimize_bound, population_size, &
      maximum_generations, number_of_restarts, mutation_factor, &
      crossover_rate, pr_tolerance, probability_epsilon, local_tolerance, &
      use_warm_start, warm_start, best_otgr, best_probabilities, &
      pr_check, success)

    real(kind=dp), intent(in) :: target_pr
    logical, intent(in) :: minimize_bound
    integer, intent(in) :: population_size, maximum_generations
    integer, intent(in) :: number_of_restarts
    real(kind=dp), intent(in) :: mutation_factor, crossover_rate
    real(kind=dp), intent(in) :: pr_tolerance, probability_epsilon
    real(kind=dp), intent(in) :: local_tolerance
    logical, intent(in) :: use_warm_start
    real(kind=dp), intent(in) :: warm_start(:)
    real(kind=dp), intent(out) :: best_otgr
    real(kind=dp), intent(out) :: best_probabilities(:)
    real(kind=dp), intent(out) :: pr_check
    logical, intent(out) :: success

    integer :: nvar, restart, generation, member, component
    integer :: r1, r2, r3, forced_component
    real(kind=dp), allocatable :: population(:,:), fitness(:)
    real(kind=dp), allocatable :: trial(:), mutant(:), full_candidate(:)
    real(kind=dp), allocatable :: restart_best(:)
    logical, allocatable :: feasible(:)
    real(kind=dp) :: random_value, trial_fitness, restart_fitness
    real(kind=dp) :: candidate_pr
    logical :: trial_feasible, restart_success

    success = .false.
    best_otgr = -1.0_dp
    pr_check = -1.0_dp
    best_probabilities = -1.0_dp

    if (d < 1) return
    if (size(best_probabilities) /= d) return
    if (size(warm_start) /= d) return
    if (target_pr < 0.0_dp .or. target_pr > 1.0_dp) return

    if (target_pr <= pr_tolerance) then
      best_probabilities = 0.0_dp
      p = best_probabilities
      best_otgr = 0.0_dp
      pr_check = 0.0_dp
      success = .true.
      return
    end if

    if (1.0_dp-target_pr <= pr_tolerance) then
      best_probabilities = 1.0_dp
      p = best_probabilities
      best_otgr = 1.0_dp
      pr_check = 1.0_dp
      success = .true.
      return
    end if

    if (d == 1) then
      best_probabilities(1) = target_pr
      p = best_probabilities
      best_otgr = graduation_rate_dplus1()
      pr_check = performance_rate()
      success = abs(pr_check-target_pr) <= pr_tolerance
      return
    end if

    if (population_size < 4 .or. maximum_generations < 1 .or. &
        number_of_restarts < 1) return

    nvar = d-1
    allocate(population(nvar,population_size))
    allocate(fitness(population_size), feasible(population_size))
    allocate(trial(nvar), mutant(nvar), full_candidate(d))
    allocate(restart_best(d))

    do restart = 1, number_of_restarts
      call initialize_population(population, probability_epsilon)

      ! The homogeneous configuration p_k=target_pr is feasible for any
      ! fixed persistence profile and guarantees a feasible initial member.
      population(:,1) = max(probability_epsilon,min(1.0_dp,target_pr))

      if (restart == 1 .and. use_warm_start .and. population_size >= 2) then
        population(:,2) = max(probability_epsilon, &
          min(1.0_dp,warm_start(1:nvar)))
      end if

      do member = 1, population_size
        call evaluate_reduced_candidate(population(:,member), target_pr, &
          pr_tolerance, probability_epsilon, fitness(member), &
          candidate_pr, full_candidate, feasible(member))
      end do

      do generation = 1, maximum_generations
        do member = 1, population_size
          call choose_distinct_indices(population_size, member, r1, r2, r3)

          mutant = population(:,r1) + mutation_factor * &
            (population(:,r2)-population(:,r3))
          mutant = max(probability_epsilon,min(1.0_dp,mutant))

          call random_number(random_value)
          forced_component = 1 + int(random_value*real(nvar,kind=dp))
          if (forced_component > nvar) forced_component = nvar

          trial = population(:,member)
          do component = 1, nvar
            call random_number(random_value)
            if (random_value <= crossover_rate .or. &
                component == forced_component) then
              trial(component) = mutant(component)
            end if
          end do

          call evaluate_reduced_candidate(trial, target_pr, pr_tolerance, &
            probability_epsilon, trial_fitness, candidate_pr, &
            full_candidate, trial_feasible)

          if (trial_feasible) then
            if (.not. feasible(member) .or. &
                is_better(trial_fitness,fitness(member),minimize_bound)) then
              population(:,member) = trial
              fitness(member) = trial_fitness
              feasible(member) = .true.
            end if
          end if
        end do
      end do

      call best_population_member(population,fitness,feasible,target_pr, &
        pr_tolerance,probability_epsilon,minimize_bound,restart_fitness, &
        restart_best,restart_success)

      if (restart_success) then
        if (.not. success .or. &
            is_better(restart_fitness,best_otgr,minimize_bound)) then
          best_otgr = restart_fitness
          best_probabilities = restart_best
          success = .true.
        end if
      end if
    end do

    if (success) then
      call coordinate_refinement(target_pr,minimize_bound,pr_tolerance, &
        probability_epsilon,local_tolerance,best_otgr,best_probabilities)
      p = best_probabilities
      pr_check = performance_rate()
      success = abs(pr_check-target_pr) <= pr_tolerance
    end if

    deallocate(population,fitness,feasible,trial,mutant)
    deallocate(full_candidate,restart_best)
  end subroutine optimize_pr_bound


  subroutine initialize_population(population, probability_epsilon)
    real(kind=dp), intent(out) :: population(:,:)
    real(kind=dp), intent(in) :: probability_epsilon

    call random_number(population)
    population = probability_epsilon + &
      (1.0_dp-probability_epsilon)*population
  end subroutine initialize_population


  subroutine evaluate_reduced_candidate(reduced, target_pr, pr_tolerance, &
      probability_epsilon, objective, calculated_pr, full_candidate, &
      feasible)

    real(kind=dp), intent(in) :: reduced(:), target_pr
    real(kind=dp), intent(in) :: pr_tolerance, probability_epsilon
    real(kind=dp), intent(out) :: objective, calculated_pr
    real(kind=dp), intent(out) :: full_candidate(:)
    logical, intent(out) :: feasible

    real(kind=dp) :: last_probability
    logical :: reconstructed

    feasible = .false.
    objective = -1.0_dp
    calculated_pr = -1.0_dp
    full_candidate = -1.0_dp

    call reconstruct_last_probability(reduced,target_pr, &
      probability_epsilon,last_probability,reconstructed)
    if (.not. reconstructed) return

    full_candidate(1:d-1) = reduced
    full_candidate(d) = last_probability
    p = full_candidate

    calculated_pr = performance_rate()
    if (calculated_pr < 0.0_dp) return
    if (abs(calculated_pr-target_pr) > pr_tolerance) return

    objective = graduation_rate_dplus1()
    if (objective < 0.0_dp) return
    feasible = .true.
  end subroutine evaluate_reduced_candidate


  subroutine reconstruct_last_probability(reduced,target_pr, &
      probability_epsilon,last_probability,success)

    real(kind=dp), intent(in) :: reduced(:), target_pr
    real(kind=dp), intent(in) :: probability_epsilon
    real(kind=dp), intent(out) :: last_probability
    logical, intent(out) :: success

    integer :: k
    real(kind=dp) :: sum_numerator, sum_denominator, survival_value
    real(kind=dp) :: denominator, c_value, coefficient, constant_term
    real(kind=dp) :: scale

    success = .false.
    last_probability = -1.0_dp
    sum_numerator = 0.0_dp
    sum_denominator = 0.0_dp
    survival_value = 1.0_dp

    do k = 1, d-1
      if (reduced(k) < probability_epsilon .or. &
          reduced(k) > 1.0_dp) return

      denominator = 1.0_dp-lambda(k)*(1.0_dp-reduced(k))
      if (denominator <= 0.0_dp) return

      survival_value = survival_value*reduced(k)/denominator
      sum_numerator = sum_numerator+survival_value
      sum_denominator = sum_denominator+survival_value/reduced(k)
    end do

    c_value = 1.0_dp-lambda(d)
    coefficient = target_pr*sum_denominator*lambda(d) - &
      sum_numerator*lambda(d)-survival_value
    constant_term = target_pr*sum_denominator*c_value + &
      target_pr*survival_value-sum_numerator*c_value

    scale = max(1.0_dp,abs(coefficient),abs(constant_term))

    if (abs(coefficient) > 100.0_dp*epsilon(1.0_dp)*scale) then
      last_probability = -constant_term/coefficient
    else
      call bisect_last_probability(sum_numerator,sum_denominator, &
        survival_value,target_pr,probability_epsilon,last_probability,success)
      return
    end if

    if (last_probability < -probability_epsilon .or. &
        last_probability > 1.0_dp+probability_epsilon) return

    last_probability = max(0.0_dp,min(1.0_dp,last_probability))
    success = .true.
  end subroutine reconstruct_last_probability


  subroutine bisect_last_probability(sum_numerator,sum_denominator, &
      survival_value,target_pr,probability_epsilon,last_probability,success)

    real(kind=dp), intent(in) :: sum_numerator, sum_denominator
    real(kind=dp), intent(in) :: survival_value, target_pr
    real(kind=dp), intent(in) :: probability_epsilon
    real(kind=dp), intent(out) :: last_probability
    logical, intent(out) :: success

    integer, parameter :: maximum_iterations = 200
    integer :: iteration
    real(kind=dp) :: left, right, middle, f_left, f_right, f_middle

    success = .false.
    last_probability = -1.0_dp
    left = probability_epsilon
    right = 1.0_dp
    f_left = reduced_pr(left,sum_numerator,sum_denominator,survival_value) - &
      target_pr
    f_right = reduced_pr(right,sum_numerator,sum_denominator,survival_value) - &
      target_pr

    if (f_left*f_right > 0.0_dp) return

    do iteration = 1, maximum_iterations
      middle = 0.5_dp*(left+right)
      f_middle = reduced_pr(middle,sum_numerator,sum_denominator, &
        survival_value)-target_pr

      if (abs(f_middle) <= 1.0e-13_dp .or. &
          abs(right-left) <= 1.0e-13_dp) exit

      if (f_left*f_middle <= 0.0_dp) then
        right = middle
        f_right = f_middle
      else
        left = middle
        f_left = f_middle
      end if
    end do

    last_probability = 0.5_dp*(left+right)
    success = .true.
  end subroutine bisect_last_probability


  real(kind=dp) function reduced_pr(last_probability,sum_numerator, &
      sum_denominator,survival_value)

    real(kind=dp), intent(in) :: last_probability
    real(kind=dp), intent(in) :: sum_numerator, sum_denominator
    real(kind=dp), intent(in) :: survival_value

    real(kind=dp) :: denominator

    denominator = 1.0_dp-lambda(d)*(1.0_dp-last_probability)
    if (denominator <= 0.0_dp) then
      reduced_pr = -1.0_dp
      return
    end if

    reduced_pr = (sum_numerator + &
      survival_value*last_probability/denominator) / &
      (sum_denominator+survival_value/denominator)
  end function reduced_pr


  subroutine choose_distinct_indices(population_size,excluded,r1,r2,r3)
    integer, intent(in) :: population_size, excluded
    integer, intent(out) :: r1, r2, r3

    r1 = random_index(population_size)
    do while (r1 == excluded)
      r1 = random_index(population_size)
    end do

    r2 = random_index(population_size)
    do while (r2 == excluded .or. r2 == r1)
      r2 = random_index(population_size)
    end do

    r3 = random_index(population_size)
    do while (r3 == excluded .or. r3 == r1 .or. r3 == r2)
      r3 = random_index(population_size)
    end do
  end subroutine choose_distinct_indices


  integer function random_index(upper_bound)
    integer, intent(in) :: upper_bound
    real(kind=dp) :: random_value

    call random_number(random_value)
    random_index = 1+int(random_value*real(upper_bound,kind=dp))
    if (random_index > upper_bound) random_index = upper_bound
  end function random_index


  logical function is_better(candidate,current,minimize_bound)
    real(kind=dp), intent(in) :: candidate,current
    logical, intent(in) :: minimize_bound

    if (minimize_bound) then
      is_better = candidate < current
    else
      is_better = candidate > current
    end if
  end function is_better


  subroutine best_population_member(population,fitness,feasible,target_pr, &
      pr_tolerance,probability_epsilon,minimize_bound,best_fitness, &
      best_candidate,success)

    real(kind=dp), intent(in) :: population(:,:),fitness(:)
    logical, intent(in) :: feasible(:),minimize_bound
    real(kind=dp), intent(in) :: target_pr,pr_tolerance,probability_epsilon
    real(kind=dp), intent(out) :: best_fitness,best_candidate(:)
    logical, intent(out) :: success

    integer :: member
    real(kind=dp) :: objective,calculated_pr
    real(kind=dp), allocatable :: candidate(:)
    logical :: candidate_feasible

    success = .false.
    best_fitness = -1.0_dp
    best_candidate = -1.0_dp
    allocate(candidate(d))

    do member = 1, size(fitness)
      if (.not. feasible(member)) cycle
      call evaluate_reduced_candidate(population(:,member),target_pr, &
        pr_tolerance,probability_epsilon,objective,calculated_pr, &
        candidate,candidate_feasible)
      if (.not. candidate_feasible) cycle

      if (.not. success .or. &
          is_better(objective,best_fitness,minimize_bound)) then
        best_fitness = objective
        best_candidate = candidate
        success = .true.
      end if
    end do

    deallocate(candidate)
  end subroutine best_population_member


  subroutine coordinate_refinement(target_pr,minimize_bound,pr_tolerance, &
      probability_epsilon,local_tolerance,best_fitness,best_candidate)

    real(kind=dp), intent(in) :: target_pr,pr_tolerance
    real(kind=dp), intent(in) :: probability_epsilon,local_tolerance
    logical, intent(in) :: minimize_bound
    real(kind=dp), intent(inout) :: best_fitness,best_candidate(:)

    integer, parameter :: maximum_sweeps = 500
    integer :: sweep,component,direction
    real(kind=dp) :: step,objective,calculated_pr
    real(kind=dp), allocatable :: reduced(:),trial(:),full_candidate(:)
    logical :: feasible,improved

    if (d <= 1) return
    allocate(reduced(d-1),trial(d-1),full_candidate(d))
    reduced = best_candidate(1:d-1)
    step = 0.05_dp

    do sweep = 1, maximum_sweeps
      improved = .false.

      do component = 1, d-1
        do direction = -1, 1, 2
          trial = reduced
          trial(component) = trial(component) + &
            real(direction,kind=dp)*step
          trial(component) = max(probability_epsilon, &
            min(1.0_dp,trial(component)))

          call evaluate_reduced_candidate(trial,target_pr,pr_tolerance, &
            probability_epsilon,objective,calculated_pr,full_candidate, &
            feasible)

          if (feasible .and. &
              is_better(objective,best_fitness,minimize_bound)) then
            reduced = trial
            best_candidate = full_candidate
            best_fitness = objective
            improved = .true.
          end if
        end do
      end do

      if (.not. improved) step = 0.5_dp*step
      if (step <= local_tolerance) exit
    end do

    deallocate(reduced,trial,full_candidate)
  end subroutine coordinate_refinement

end module Optimizacion
