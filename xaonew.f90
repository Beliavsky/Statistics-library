! https://help.imsl.com/fortran/7.0/stat/stat.htm?turl=docs%2Faonew.htm
use kind_mod, only: dp
use anova_mod, only: aonew
implicit none
integer, parameter :: ngroup = 3, nobs = 6, ni(ngroup) = [3, 2, 1], iprint = 3
real(kind=dp) :: aov(15), y(nobs)
data y/101.0, 105.0, 94.0, 84.0, 88.0, 32.0/
call aonew (ni, y, aov, iprint=iprint)
end