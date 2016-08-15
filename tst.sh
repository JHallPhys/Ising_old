#!/bin/bash

gfortran con.f90 rng.f90 flippy.f90 spin_en.f90 neighbour.f90 probability.f90  test.f90 -o test.exe

./test.exe


#xmgrace mt.dat &
#xmgrace mtsqr.dat &
#xmgrace en.dat &
#xmgrace ensqr.dat &
#xmgrace cv.dat &
#xmgrace msus.dat & 
#xmgrace cvp.dat & 
