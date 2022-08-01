! "Nonadiabatic transition paths from quantum jump trajectories"
! Michelle C Anderson, Addison J. Schile, and David T. Limmer, University of California Berkeley
! Copyright 2022

! This file is part of TPT_CI 
! TPT_CI is free software: you can redistribute it and/or modify it under the terms of the GNU General 
! Public License as published by the Free Software Foundation, either version 3 of the License, 
! or (at your option) any later version.
! TPT_CI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with TPT_CI. 
! If not, see <https://www.gnu.org/licenses/>.


! All parameters should be speicifed in atomic units to avoid any
! complications.
! Parameters to specify how the run proceeds
MODULE parameters
IMPLICIT NONE

! In all cases save where otherwise noted, 1 is a 'true' or 'yes' on a binary
! switch and '0' is 'no' or 'false' but usually any number except 1 is ignored.

! Option to print extra information to deal with problems
LOGICAL, PARAMETER :: DEBUG = .FALSE.

! Basic system setup parameters
INTEGER, PARAMETER :: CI = 0  ! Conical intersection

! What kind of run to perform (CI is the only setup included in this code, so leave unchanged)
INTEGER, PARAMETER :: run_type = CI  
! Should the CI off-diagonal coupling be \lambda*|Q_c| rather than \lambda*Q_c 
! Not covered in the paper but an intersting avenue of future exploration
INTEGER, PARAMETER :: abs_coupling = 1

! When creating a markov state model, how many timesteps to run (tau = dt*msm_steps)
INTEGER, PARAMETER :: msm_steps = 50

! Only set one of these options to 1 at any time
! Run the secular Redfield approximation (1)
INTEGER, PARAMETER :: secular = 0
! Run lindblad jump dynamics  (1)
INTEGER, PARAMETER :: lb_switch = 1

! Which eigenstates LB propogation should end with?
! Propagation will continue until reaching one of these states.
INTEGER, PARAMETER :: flag_estate_1 = 1  ! First eigenstate to stop propogation (should be 1)
INTEGER, PARAMETER :: flag_estate_2 = 3  ! Second eigenstate to stop propogation (should not be 1)

! Number of trajectories for a lindblad run to carry out
INTEGER, PARAMETER :: ntraj = 10000
! Duration for the Redfield propagation (au)
REAL*8, PARAMETER :: duration = 10000.0d0
! Time step for rk4 routine
REAL*8, PARAMETER :: dt = 5.000d0 !duration / (ntimes - 1)
! How many iterations of redfield to do between printings of position/diabatic state
INTEGER, PARAMETER :: print_num = 10

! Constants and conversions
REAL*8, PARAMETER :: pi = 3.1415926d0
REAL*8, PARAMETER :: hbar = 1.0d0 

! Bath parameters
! Number of baths which will be coupled to some coordinate; should be 2
INTEGER, PARAMETER :: nbath = 2
! Cutoff parameter in spectral density for each bath
REAL*8, PARAMETER, DIMENSION(nbath) :: omegac = (/0.00048378909d0, 0.00048378909d0/)
! Strength parameter in spectral density
REAL*8, PARAMETER, DIMENSION(nbath) :: eta = (/9.3453302887d-9, 9.3453302887d-9/) 
! Debye bath is 1, Ohmic is 0
! ohmic is the default if no valid type is given
INTEGER, PARAMETER, DIMENSION(nbath) :: bath_type = (/1,1/)
! Temperature 1/(kb*T). To indicate zero temperature set beta = -1.0d0;
! any negative value will do but please use -1.0d0
REAL*8, PARAMETER :: beta = 1052.584412992859d0

! IMPORTANT!!! Note that internally and in output files the HIGHER ENERGY WELL IS WELL 2, 
! NOT WELL 1 AS IN THE PAPER.
! The order of well 1 and well 2 was changed in the writeup but the original order has not
! been altered in the setup code to avoid potentially introducing errors in the code or plots.

! Conical Intersection Parameters
INTEGER, PARAMETER :: dim_nc = 40   ! Basis set size, coupling coodinate
INTEGER, PARAMETER :: dim_nt = 110   ! Basis set size, tuning coordiante
REAL*8, PARAMETER :: omegas_c = 1.00*0.0043364174d0 ! Coupling coordinate frequency
REAL*8, PARAMETER :: omegas_t = 1.00*0.0027194482d0 ! Tuning coordinate frequency
REAL*8, PARAMETER, DIMENSION(2) :: kappa_t = (/-0.0038586765d0*3.0d0, 0.0054756457d0*2.4d0/) ! Oscillator displacements
REAL*8, PARAMETER, DIMENSION(2) :: exe = (/0.144792242d0, 0.177866612d0*0.87d0/) ! Base energies of both wells
REAL*8, PARAMETER :: lambda_s = 0.0096283166d0*1.00d0 ! Diabatic coupling strength

! Size parameters
! Dimension of density matrix; must be an even number
INTEGER, PARAMETER :: n = 700   ! Truncated basis size, must be less than or equal to m
INTEGER, PARAMETER :: m = dim_nc*dim_nt*2 ! Basis size for initial diagonalization (dim_nc*dim_nt*2 for a 2-surface CI)

! The default initialization is vertical excitation
! Set init_stat_therm to 1 to initialize in an eigenstate
! Populate a single eigenstate
INTEGER, PARAMETER :: init_stat_therm = 0
! Eigenstate to initialize all density in if init_stat_therm = 1
INTEGER, PARAMETER :: init_es = 4
! Initialize in either well 1 of the CI or well 2
! Well 2 is the one used for vertical excitation in the paper. Again
! NOTE THAT WELL 1 AND WELL 2 SWITCHED NAMES in the paper writeup
INTEGER, PARAMETER :: init_well = 2

! Lindblad jumps tolerance for wfn norm to match a random number
! during a binary search
REAL*8, PARAMETER :: tol = 1.0d-5

! Plotting bounds
INTEGER, PARAMETER :: bound1 = 8
INTEGER, PARAMETER :: bound2 = 14
INTEGER, PARAMETER :: bound3 = 16
INTEGER, PARAMETER :: bound4 = 21

END MODULE parameters
