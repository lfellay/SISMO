# New zone table on the r21 grid
# class convention: 1=radiative, 2=convective, 3=atmosphere
# density factors radiative convective atmosphere:  1.000000E+00  1.250000E+00  1.500000E+00
# MULTI_P, Ntarget, N5, nz_new
          16       25000       25005           4
# k  i_start  i_end  class
         1         1      5376         2
         2      5376     15839         1
         3     15839     22062         2
         4     22062     25005         3
