program test_estimate_missing
  use kind_mod
  use time_series_mod
  implicit none
  integer, parameter :: n = 200
  integer :: i, j, k, n_obs, n_miss, ntimes, miss_ind
  integer, dimension(n) :: tpoints
  integer, allocatable :: times(:), missing_index(:), times2(:)
  real(kind=dp), dimension(n) :: y, x1, x2
  real(kind=dp), allocatable :: result(:)

  ! initialize tpoints with values 1,2,...,200
  tpoints = (/ (i, i = 1, n) /)

  ! initialize y with the given data
  y = (/ &
    1.30540_dp, -1.37166_dp, 1.47905_dp, -0.91059_dp, 1.36191_dp, -2.16966_dp, 3.11254_dp, -1.99536_dp, 2.29740_dp, -1.82474_dp, &
   -0.25445_dp,  0.33519_dp, -0.25480_dp, -0.50574_dp, -0.21429_dp, -0.45932_dp, -0.63813_dp,  0.25646_dp, -0.46243_dp, -0.44104_dp, &
    0.42733_dp,  0.61102_dp, -0.82417_dp,  1.48537_dp, -1.57733_dp, -0.09846_dp,  0.46311_dp,  0.49156_dp, -1.66090_dp,  2.02808_dp, &
   -1.45768_dp,  1.36115_dp, -0.65973_dp,  1.13332_dp, -0.86285_dp,  1.23848_dp, -0.57301_dp, -0.28210_dp,  0.20195_dp,  0.06981_dp, &
    0.28454_dp,  0.19745_dp, -0.16490_dp, -1.05019_dp,  0.78652_dp, -0.40447_dp,  0.71514_dp, -0.90003_dp,  1.83604_dp, -2.51205_dp, &
    1.00526_dp, -1.01683_dp,  1.70691_dp, -1.86564_dp,  1.84912_dp, -1.33120_dp,  2.35105_dp, -0.45579_dp, -0.57773_dp, -0.55226_dp, &
    0.88371_dp,  0.23138_dp,  0.59984_dp,  0.31971_dp,  0.59849_dp,  0.41873_dp, -0.46955_dp,  0.53003_dp, -1.17203_dp,  1.52937_dp, &
   -0.48017_dp, -0.93830_dp,  1.00651_dp, -1.41493_dp, -0.42188_dp, -0.67010_dp,  0.58079_dp, -0.96193_dp,  0.22763_dp, -0.92214_dp, &
    1.35697_dp, -1.47008_dp,  2.47841_dp, -1.50522_dp,  0.41650_dp, -0.21669_dp, -0.90297_dp,  0.00274_dp, -1.04863_dp,  0.66192_dp, &
   -0.39143_dp,  0.40779_dp, -0.68174_dp, -0.04700_dp, -0.84469_dp,  0.30735_dp, -0.68412_dp,  0.25888_dp, -1.08642_dp,  0.52928_dp, &
    0.72168_dp, -0.18199_dp, -0.09499_dp,  0.67610_dp,  0.14636_dp,  0.46846_dp, -0.13989_dp,  0.50856_dp, -0.22268_dp,  0.92756_dp, &
    0.73069_dp,  0.78998_dp, -1.01650_dp,  1.25637_dp, -2.36179_dp,  1.99616_dp, -1.54326_dp,  1.38220_dp,  0.19674_dp, -0.85241_dp, &
    0.40463_dp,  0.39523_dp, -0.60721_dp,  0.25041_dp, -1.24967_dp,  0.26727_dp,  1.40042_dp, -0.66963_dp,  1.26049_dp, -0.92074_dp, &
    0.05909_dp, -0.61926_dp,  1.41550_dp,  0.25537_dp, -0.13240_dp, -0.07543_dp,  0.10413_dp,  1.42445_dp, -1.37379_dp,  0.44382_dp, &
   -1.57210_dp,  2.04702_dp, -2.22450_dp,  1.27698_dp,  0.01073_dp, -0.88459_dp,  0.88194_dp, -0.25019_dp,  0.70224_dp, -0.41855_dp, &
    0.93850_dp,  0.36007_dp, -0.46043_dp,  0.18645_dp,  0.06337_dp,  0.29414_dp, -0.20054_dp,  0.83078_dp, -1.62530_dp,  2.64925_dp, &
   -1.25355_dp,  1.59094_dp, -1.00684_dp,  1.03196_dp, -1.58045_dp,  2.04295_dp, -2.38264_dp,  1.65095_dp, -0.33273_dp, -1.29092_dp, &
    0.14020_dp, -0.11434_dp,  0.04392_dp,  0.05293_dp, -0.42277_dp,  0.59143_dp, -0.03347_dp, -0.58457_dp,  0.87030_dp,  0.19985_dp, &
   -0.73500_dp,  0.73640_dp,  0.29531_dp,  0.22325_dp, -0.60035_dp,  1.42253_dp, -1.11278_dp,  1.30468_dp, -0.41923_dp, -0.38019_dp, &
    0.50937_dp,  0.23051_dp,  0.46496_dp,  0.02459_dp, -0.68478_dp,  0.25821_dp,  1.17655_dp, -2.26629_dp,  1.41173_dp, -0.68331_dp /)

  ! full series x1 is just y
  x1 = y

  ! allocate times2 to build the series with missing observations
  allocate(times2(n))
  k = 0
  ! include first observation always
  times2(1) = tpoints(1)
  x2(1) = y(1)
  do i = 1, n - 1
     ! remove observations at indices 130, 140, 141, 160, 175, and 176
     if (i /= 129 .and. i /= 139 .and. i /= 140 .and. i /= 159 .and. i /= 174 .and. i /= 175) then
        k = k + 1
        times2(k+1) = tpoints(i+1)
        x2(k+1) = y(i+1)
     end if
  end do
  n_obs = k + 1
  ntimes = tpoints(n) - tpoints(1) + 1
  n_miss = ntimes - n_obs

  ! allocate full time axis and missing index arrays
  allocate(times(ntimes))
  allocate(missing_index(n_miss))

  do j = 0, 3
     if (j <= 2) then
        call estimate_missing(times2(1:n_obs), x2(1:n_obs), result, nobsw=n_obs, &
             imeth=j, itime_points_full=times, miss_index=missing_index)
     else
        call estimate_missing(times2(1:n_obs), x2(1:n_obs), result, nobsw=n_obs, &
             imeth=j, maxlag=20, itime_points_full=times, miss_index=missing_index)
     end if
     write(*,*) "method: ", j
     write(*,*) "time      actual      predicted      difference"
     do i = 1, size(missing_index)
        miss_ind = missing_index(i)
        write(*, '(i4,6x,f6.3,8x,f6.3,7x,f6.3)') times(miss_ind), x1(miss_ind), result(miss_ind), &
             abs(x1(miss_ind)-result(miss_ind))
     end do
     write(*,*)
     if (allocated(result)) then
        deallocate(result)
     end if
  end do

  deallocate(times)
  deallocate(missing_index)
  deallocate(times2)

end program test_estimate_missing
