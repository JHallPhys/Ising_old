module constants

implicit none 

integer,parameter :: dp = selected_real_kind(15,300)
real(kind=dp),parameter :: k_boltz = 1.3806503E-23
real(kind=dp),parameter :: t_curie = 2.269185
real(kind=dp),parameter :: rt_2 = 1.41421356237

end module 
