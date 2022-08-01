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

! Sorts an ensemble of trajectories into two 
! files based on what the last eigenstate in 
! the trajectory is.

! Will give average and block averging standard deviation information
! after the run.


PROGRAM sort_trajectories
USE parameters
USE prequel
USE lindblad_coherent
USE lindblad_jumps
USE rk_ns_utils
USE set_up_H

IMPLICIT NONE

! The flag eigenstates which indicate the two different trajectory fates to distinguish
INTEGER :: end_state_1 = flag_estate_1
INTEGER :: end_state_2 = flag_estate_2
! How large to make "ending"
INTEGER :: capacity = ntraj
! Which trajectory has which ending; the read through twice strategy is employed
INTEGER, DIMENSION(:), ALLOCATABLE :: ending
! Index, not used, eigenstate, previous eigenstate seen, indicator of read success
INTEGER :: i, id, estate, last_estate, success
REAL*8 :: time
TYPE(wfunc) :: wfn
TYPE(jump_data) :: jump_info
INTEGER :: total_ending_1, total_ending_2  ! Total number of each ending

! I would like an estimate of the standard deviation of the yield
! so I will procure some subsamples
INTEGER, PARAMETER :: subs = 5
! Name is due to my original interest in having 5 subsamples
! real*8 because the sdev functions take reals so why not?
REAL*8, DIMENSION(subs) :: total_ending_1_5, total_ending_2_5

OPEN(15,FILE="jump_data.txt")
CALL init_wfn(wfn)
last_estate = 0
total_ending_1_5 = 0.0d0
total_ending_2_5 = 0.0d0

ALLOCATE(ending(capacity))
ending = -1

total_ending_1 = 0
total_ending_2 = 0

! First, read through the file once to determine where everything goes
DO i = 1, capacity
  DO 
    success = read_dj_line(15,id,estate,time,wfn,jump_info)
    IF (success .NE. 0) THEN
       WRITE(*,*) "Read error. Exiting on ", i
       EXIT
    END IF
    IF (estate .EQ. -1) THEN
      ending(i) = last_estate
      WRITE(*,*) "Ending ", i, "is", ending(i)
      EXIT 
    END IF
    last_estate = estate
  END DO
  IF (success .NE. 0) THEN  ! Get out of the outer loop, too
    EXIT
  END IF
END DO

CLOSE(15)


OPEN(15,FILE="jump_data.txt")
OPEN(16,FILE="ending1.txt")
OPEN(17,FILE="ending2.txt")

! Read through again and write the trajectories to their
! proper files
DO i = 1, capacity
  IF (ending(i) .EQ. -1) THEN
    EXIT  ! Flag value; we're finished
  END IF

  IF (ending(i) .EQ. end_state_1) THEN
    WRITE(*,*) "ending 1 found", ending(i)
    total_ending_1 = total_ending_1 + 1

    total_ending_1_5(MOD(i,subs) + 1) = total_ending_1_5(MOD(i,subs) + 1) + 1

    DO 
      success = read_dj_line(15,id,estate,time,wfn,jump_info)
      IF (success .NE. 0) THEN
        WRITE(*,*) "Read error ", success
        STOP
      END IF
      WRITE(16,*) id, time, estate
      IF (estate .EQ. 0) THEN
        WRITE(16,*) org_out(wfn%fe) 
      ELSE IF (estate .EQ. -1) THEN
        EXIT
      ELSE
        WRITE(16,*) REAL(REAL(wfn%fe(estate,1))), AIMAG(wfn%fe(estate,1))
      END IF
    END DO
  ELSE IF (ending(i) .EQ. end_state_2) THEN
    WRITE(*,*) "ending 2 found", ending(i)
    total_ending_2 = total_ending_2 + 1

    total_ending_2_5(MOD(i,subs) + 1) = total_ending_2_5(MOD(i,subs) + 1) + 1

    DO 
      success = read_dj_line(15,id,estate,time,wfn,jump_info)
      IF (success .NE. 0) THEN
        WRITE(*,*) "Read error ", success
        STOP
      END IF
      WRITE(17,*) id, time, estate
      IF (estate .EQ. 0) THEN
        WRITE(17,*) org_out(wfn%fe) 
      ELSE IF (estate .EQ. -1) THEN
        EXIT
      ELSE
        WRITE(17,*) REAL(REAL(wfn%fe(estate,1))), AIMAG(wfn%fe(estate,1))
      END IF
    END DO
  ELSE
    WRITE(*,*) "neither ending 1 nor 2 found", ending(i)
  END IF

END DO

WRITE(16,*) "TOTAL", total_ending_1, "STANDARD DEVIATION FROM ", subs, "SAMPLES", &
            sdev(subs,total_ending_1_5)
WRITE(17,*) "TOTAL", total_ending_2, "STANDARD DEVIATION FROM ", subs, "SAMPLES", &
            sdev(subs,total_ending_2_5)

CLOSE(15)
CLOSE(16)
CLOSE(17)

END PROGRAM
