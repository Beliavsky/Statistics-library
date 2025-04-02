module time_series_mod
  use kind_mod, only: dp
  implicit none
contains

  subroutine estimate_missing(itime_points, w, z, nobsw, imeth, maxlag, &
    maxbc, tolbc, tolss, relerr, maxit, wmean, nobsz, itime_points_full, miss_index)
    ! implementation of missing value estimation with AR(1) method implemented
    integer, intent(in) :: itime_points(:)
    real(kind=dp), intent(in) :: w(:)
    real(kind=dp), allocatable, intent(out) :: z(:)
    integer, intent(in), optional :: nobsw, imeth, maxlag, maxbc, maxit
    real(kind=dp), intent(in), optional :: tolbc, tolss, relerr, wmean
    integer, intent(out), optional :: nobsz
    integer, allocatable, intent(out), optional :: itime_points_full(:)
    integer, allocatable, intent(out), optional :: miss_index(:)

    integer :: nobs, i, full_length, ip, j, gap_start, gap_end
    integer :: method
    integer :: start_time
    real(kind=dp) :: median_val, dt, slope
    real(kind=dp), allocatable :: median_array(:)
    integer :: n_before, n_after, k, idx_prev, idx_next
    integer :: count_miss
    real(kind=dp) :: mu, phi, forecast, prev_obs, last_obs
    real(kind=dp), parameter :: eps = 1e-8_dp

    ! set defaults
    if (present(nobsw)) then
       nobs = nobsw
    else
       nobs = size(w)
    end if
    if (present(imeth)) then
       method = imeth
    else
       method = 3
    end if

    if (present(wmean)) then
       mu = wmean
    else
       mu = 0.0_dp
    end if

    start_time = itime_points(1)
    full_length = itime_points(nobs) - start_time + 1

    allocate(z(full_length))

    ! fill full time axis if requested
    if (present(itime_points_full)) then
       allocate(itime_points_full(full_length))
       do i = 1, full_length
          itime_points_full(i) = start_time + i - 1
       end do
    end if

    ! initialize z with a flag value
    z = -9999.0_dp

    ! place observed values into z at correct positions
    do i = 1, nobs
       ip = itime_points(i) - start_time + 1
       z(ip) = w(i)
    end do

    ! record missing indices if requested
    count_miss = 0
    do i = 1, full_length
       if (z(i) == -9999.0_dp) then
          count_miss = count_miss + 1
       end if
    end do
    if (present(miss_index)) then
       allocate(miss_index(count_miss))
       count_miss = 0
       do i = 1, full_length
          if (z(i) == -9999.0_dp) then
             count_miss = count_miss + 1
             miss_index(count_miss) = i
          end if
       end do
    end if

    ! process gaps of missing values in increasing time order
    i = 1
    do while (i <= full_length)
       if (z(i) == -9999.0_dp) then
          gap_start = i
          j = i
          do while (j <= full_length .and. z(j) == -9999.0_dp)
             j = j + 1
          end do
          gap_end = j - 1
          ! if gap is at beginning or end, assign a default value (here zero)
          if (gap_start == 1 .or. j > full_length) then
             do ip = gap_start, gap_end
                z(ip) = 0.0_dp
             end do
          else
             idx_prev = gap_start - 1
             idx_next = j
             select case (method)
             case (0)  ! median method
                n_before = min(4, idx_prev)
                n_after = min(4, full_length - idx_next + 1)
                allocate(median_array(n_before+n_after))
                k = 0
                do ip = idx_prev - n_before + 1, idx_prev
                   k = k + 1
                   median_array(k) = z(ip)
                end do
                do ip = idx_next, idx_next + n_after - 1
                   k = k + 1
                   median_array(k) = z(ip)
                end do
                call sort_array(median_array)
                median_val = compute_median(median_array)
                deallocate(median_array)
                do ip = gap_start, gap_end
                   z(ip) = median_val
                end do
             case (1)  ! cubic spline interpolation (dummy: linear interpolation)
                do ip = gap_start, gap_end
                   dt = real(ip - idx_prev, kind=dp)
                   slope = (z(idx_next) - z(idx_prev)) / real(idx_next - idx_prev, kind=dp)
                   z(ip) = z(idx_prev) + slope * dt
                end do
             case (2)  ! ar(1) forecast
                ! use the last two values before the gap; if available, use indices gap_start-2 and gap_start-1
                if (gap_start >= 3) then
                   prev_obs = z(gap_start - 2)
                   last_obs = z(gap_start - 1)
                else
                   prev_obs = mu
                   last_obs = z(gap_start - 1)
                end if
                do ip = gap_start, gap_end
                   if (abs(prev_obs - mu) < eps) then
                      phi = 0.0_dp
                   else
                      phi = (last_obs - mu) / (prev_obs - mu)
                   end if
                   forecast = mu + phi * (last_obs - mu)
                   z(ip) = forecast
                   ! update previous observations for next forecast in the gap
                   prev_obs = last_obs
                   last_obs = forecast
                end do
             case (3)  ! ar(p) forecast (dummy: use last observed value)
                do ip = gap_start, gap_end
                   z(ip) = z(idx_prev)
                end do
             end select
          end if
          i = gap_end + 1
       else
          i = i + 1
       end if
    end do

    if (present(nobsz)) then
       nobsz = full_length
    end if

  end subroutine estimate_missing

  subroutine sort_array(a)
    real(kind=dp), intent(inout) :: a(:)
    integer :: i, j, n
    real(kind=dp) :: temp
    n = size(a)
    do i = 1, n - 1
       do j = i + 1, n
          if (a(i) > a(j)) then
             temp = a(i)
             a(i) = a(j)
             a(j) = temp
          end if
       end do
    end do
  end subroutine sort_array

  function compute_median(a) result(med)
    real(kind=dp), intent(in) :: a(:)
    real(kind=dp) :: med
    integer :: n
    n = size(a)
    if (mod(n,2) == 1) then
       med = a((n+1)/2)
    else
       med = 0.5_dp * (a(n/2) + a(n/2+1))
    end if
  end function compute_median

end module time_series_mod
