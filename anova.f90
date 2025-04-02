module anova_mod
  use iso_fortran_env
  use kind_mod, only: dp
  implicit none
  private
  public :: aonew
  private :: log_gamma, betacf
  real(real64), parameter :: nan =  transfer(-2251799813685248_int64, 1._real64)
contains

  !--------------------------------------------------------------------
  ! function log_gamma: computes the natural logarithm of the gamma
  ! function using a Lanczos approximation.
  !--------------------------------------------------------------------
  elemental real(kind=dp) function log_gamma(xx)
    real(kind=dp), intent(in) :: xx
    real(kind=dp) :: x, tmp, ser
    integer :: j
    real(kind=dp), parameter :: cof(6) = (/ &
         76.18009172947146d0, -86.50532032941677d0, 24.01409824083091d0, &
         -1.231739572450155d0, 0.1208650973866179d0, -0.5395239384953d0 /)
    x = xx - 1.0d0
    tmp = x + 5.5d0
    tmp = tmp - (x + 0.5d0) * dlog(tmp)
    ser = 1.000000000190015d0
    do j = 1, 6
       x = x + 1.0d0
       ser = ser + cof(j) / x
    end do
    log_gamma = -tmp + dlog(2.5066282746310005d0 * ser)
  end function log_gamma

  !--------------------------------------------------------------------
  ! function betacf: evaluates the continued fraction for the incomplete
  ! beta function by modified Lentz's method.
  !--------------------------------------------------------------------
  elemental real(kind=dp) function betacf(a, b, x)
    implicit none
    real(kind=dp), intent(in) :: a, b, x
    real(kind=dp) :: qab, qap, qam, c, d, h, aa, del
    integer :: m_int, itmax
    real(kind=dp), parameter :: eps = 3.0d-7
    itmax = 100

    qab = a + b
    qap = a + 1.0d0
    qam = a - 1.0d0
    c = 1.0d0
    d = 1.0d0 - qab * x / qap
    if (d == 0.0d0) d = 1.0d-30
    d = 1.0d0 / d
    h = d
    do m_int = 1, itmax
       aa = m_int * (b - m_int) * x / ((qam + 2.0d0 * m_int) * (a + 2.0d0 * m_int))
       d = 1.0d0 + aa * d
       if (d == 0.0d0) d = 1.0d-30
       c = 1.0d0 + aa / c
       if (c == 0.0d0) c = 1.0d-30
       d = 1.0d0 / d
       h = h * d * c
       aa = -(a + m_int) * (qab + m_int) * x / ((a + 2.0d0 * m_int) * (qap + 2.0d0 * m_int))
       d = 1.0d0 + aa * d
       if (d == 0.0d0) d = 1.0d-30
       c = 1.0d0 + aa / c
       if (c == 0.0d0) c = 1.0d-30
       d = 1.0d0 / d
       del = d * c
       h = h * del
       if (dabs(del - 1.0d0) < eps) exit
    end do
    betacf = h
  end function betacf

  !--------------------------------------------------------------------
  ! function dincbet: computes the regularized incomplete beta function.
  !--------------------------------------------------------------------
  real(kind=dp) function dincbet(a, b, x)
    implicit none
    real(kind=dp), intent(in) :: a, b, x
    real(kind=dp) :: bt
    if (x == 0.0d0 .or. x == 1.0d0) then
       bt = 0.0d0
    else
       bt = dexp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * dlog(x) + b * dlog(1.0d0 - x))
    end if
    if (x < (a + 1.0d0) / (a + b + 2.0d0)) then
       dincbet = bt * betacf(a, b, x) / a
    else
       dincbet = 1.0d0 - bt * betacf(b, a, 1.0d0 - x) / b
    end if
  end function dincbet

  !--------------------------------------------------------------------
  ! subroutine aonew: performs one-way analysis of variance.
  ! follows interface of IMSL aonew subroutine
  !
  ! Required arguments:
  !   ni   - integer vector (length ngroup) with number of responses per group.
  !   y    - real(kind=dp) vector containing responses (grouped contiguously).
  !   aov  - real(kind=dp) vector (length 15) that will hold the analysis
  !          of variance statistics.
  !
  ! Optional arguments:
  !   ngroup - number of groups (default: size(ni)).
  !   iprint - printing option (0: none; 1: print aov; 2: print stat; 3: print both)
  !   stat   - real(kind=dp) matrix (ngroup by 4) for group statistics.
  !   ldstat - leading dimension of stat (default: size(stat,1)).
  !   nmiss  - integer output; number of missing values (elements of y that are nan).
  !
  ! The aov vector elements are:
  !   1: degrees of freedom among groups
  !   2: degrees of freedom within groups
  !   3: total (corrected) degrees of freedom
  !   4: sum of squares among groups
  !   5: sum of squares within groups
  !   6: total (corrected) sum of squares
  !   7: among-groups mean square
  !   8: within-groups mean square
  !   9: F-statistic
  !  10: p-value
  !  11: R2 (in percent)
  !  12: Adjusted R2 (in percent)
  !  13: Estimated standard deviation of the error within groups
  !  14: Overall mean of y
  !  15: Coefficient of variation (in percent)
  !--------------------------------------------------------------------
  subroutine aonew(ni, y, aov, ngroup, iprint, stat, ldstat, nmiss)
    integer, intent(in) :: ni(:)
    real(kind=dp), intent(in) :: y(:)
    real(kind=dp), intent(out) :: aov(15)
    integer, intent(in), optional :: ngroup
    integer, intent(in), optional :: iprint
    real(kind=dp), intent(out), optional :: stat(:,:)
    integer, intent(in), optional :: ldstat
    integer, intent(out), optional :: nmiss

    integer :: ngrp, prnt, ldst
    integer :: i, j, count, total, nm
    integer, allocatable :: group_n(:)
    integer, allocatable :: group_start(:)
    real(kind=dp), allocatable :: group_mean(:), group_sd(:)
    real(kind=dp) :: sum, ssum, mean, overall_mean, overall_sum
    real(kind=dp) :: ss_among, ss_within, ss_total
    real(kind=dp) :: ms_among, ms_within, fstat, pval, r2, adj_r2, s_error, cv

    nm = 0
    total = 0
    overall_sum = 0.0d0
    ss_within = 0.0d0

    ! determine number of groups (default: size(ni))
    if (present(ngroup)) then
       ngrp = ngroup
    else
       ngrp = size(ni)
    end if

    if (present(iprint)) then
       prnt = iprint
    else
       prnt = 0
    end if

    if (present(stat)) then
       if (present(ldstat)) then
          ldst = ldstat
       else
          ldst = size(stat, 1)
       end if
    end if

    allocate(group_n(ngrp))
    allocate(group_mean(ngrp))
    allocate(group_sd(ngrp))
    allocate(group_start(ngrp))

    ! compute starting indices for each group in y
    group_start(1) = 1
    do i = 2, ngrp
       group_start(i) = group_start(i - 1) + ni(i - 1)
    end do

    ! compute group means and count nonmissing values
    do i = 1, ngrp
       count = 0
       sum = 0.0d0
       do j = group_start(i), group_start(i) + ni(i) - 1
          ! test for nan: a value is nan if it does not equal itself
          if (y(j) == y(j)) then
             count = count + 1
             sum = sum + y(j)
          else
             nm = nm + 1
          end if
       end do
       if (count > 0) then
          mean = sum / dble(count)
       else
          mean = 0.0d0
       end if
       group_n(i) = count
       group_mean(i) = mean
       overall_sum = overall_sum + sum
       total = total + count
    end do

    if (total > 0) then
       overall_mean = overall_sum / dble(total)
    else
       overall_mean = 0.0d0
    end if

    ! compute within-group sum of squares and group standard deviations
    do i = 1, ngrp
       ssum = 0.0d0
       do j = group_start(i), group_start(i) + ni(i) - 1
          if (y(j) == y(j)) then
             ssum = ssum + (y(j) - group_mean(i))**2
          end if
       end do
       if (group_n(i) > 1) then
          group_sd(i) = dsqrt(ssum / (dble(group_n(i)) - 1.0d0))
       else
          group_sd(i) = nan  ! set to nan when undefined
       end if
       ss_within = ss_within + ssum
    end do

    ! compute among-group sum of squares
    ss_among = 0.0d0
    do i = 1, ngrp
       ss_among = ss_among + dble(group_n(i)) * (group_mean(i) - overall_mean)**2
    end do

    ss_total = ss_among + ss_within

    ! degrees of freedom
    aov(1) = dble(ngrp - 1)         ! among groups
    aov(2) = dble(total - ngrp)      ! within groups
    aov(3) = dble(total - 1)         ! total (corrected)

    ! sum of squares
    aov(4) = ss_among
    aov(5) = ss_within
    aov(6) = ss_total

    if (aov(1) > 0.0d0) then
       ms_among = ss_among / aov(1)
    else
       ms_among = nan
    end if

    if (aov(2) > 0.0d0) then
       ms_within = ss_within / aov(2)
    else
       ms_within = nan
    end if

    aov(7) = ms_among
    aov(8) = ms_within

    if (ms_within > 0.0d0) then
       fstat = ms_among / ms_within
    else
       fstat = nan
    end if
    aov(9) = fstat

    ! compute p-value using the incomplete beta function:
    ! p = dincbet(d2/2, d1/2, d2/(d2+d1*f))
    if ((ms_within > 0.0d0) .and. (aov(1) > 0.0d0) .and. (aov(2) > 0.0d0)) then
       pval = dincbet(aov(2) / 2.0d0, aov(1) / 2.0d0, aov(2) / (aov(2) + aov(1) * fstat))
    else
       pval = nan
    end if
    aov(10) = pval

    if (ss_total > 0.0d0) then
       r2 = (ss_among / ss_total) * 100.0d0
    else
       r2 = nan
    end if
    aov(11) = r2

    if (total - ngrp > 0) then
       adj_r2 = (1.0d0 - (dble(total - 1) / dble(total - ngrp)) * (1.0d0 - ss_among / ss_total)) * 100.0d0
    else
       adj_r2 = nan
    end if
    aov(12) = adj_r2

    if (ms_within >= 0.0d0) then
       s_error = dsqrt(ms_within)
    else
       s_error = nan
    end if
    aov(13) = s_error

    aov(14) = overall_mean

    if (overall_mean /= 0.0d0) then
       cv = (s_error / overall_mean) * 100.0d0
    else
       cv = nan
    end if
    aov(15) = cv

    ! fill the stat matrix with group statistics if provided
    if (present(stat)) then
       do i = 1, ngrp
          if (i <= ldst) then
             stat(i, 1) = dble(i)           ! group number
             stat(i, 2) = dble(group_n(i))    ! number of nonmissing observations
             stat(i, 3) = group_mean(i)       ! group mean
             stat(i, 4) = group_sd(i)         ! group standard deviation
          end if
       end do
    end if

    ! printing based on iprint option:
    !   1: print aov only, 2: print stat only, 3: print both
    if (prnt == 1 .or. prnt == 3) then
       write(*,*) '*** analysis of variance ***'
       write(*,*) 'df among groups = ', aov(1)
       write(*,*) 'df within groups = ', aov(2)
       write(*,*) 'total df = ', aov(3)
       write(*,*) 'ss among groups = ', aov(4)
       write(*,*) 'ss within groups = ', aov(5)
       write(*,*) 'total ss = ', aov(6)
       write(*,*) 'ms among groups = ', aov(7)
       write(*,*) 'ms within groups = ', aov(8)
       write(*,*) 'f statistic = ', aov(9)
       write(*,*) 'p value = ', aov(10)
       write(*,*) 'r2 (percent) = ', aov(11)
       write(*,*) 'adjusted r2 (percent) = ', aov(12)
       write(*,*) 'estimated std dev error = ', aov(13)
       write(*,*) 'overall mean = ', aov(14)
       write(*,*) 'coefficient of variation (percent) = ', aov(15)
    end if

    if (prnt == 2 .or. prnt == 3) then
       write(*,*) '*** group statistics ***'
       write(*,*) 'group   n   mean   std deviation'
       do i = 1, ngrp
          write(*,*) i, group_n(i), group_mean(i), group_sd(i)
       end do
    end if

    if (present(nmiss)) then
       nmiss = nm
    end if

    deallocate(group_n)
    deallocate(group_mean)
    deallocate(group_sd)
    deallocate(group_start)
  end subroutine aonew

end module anova_mod