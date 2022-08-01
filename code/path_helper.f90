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

! Reads through an ensemble of jump trajectories and counts
! how many jumps occur between all distinct eigenstates then
! prepares a jump graph for analysis in the same manner
! as a Markov state model. Block averaging error estimates
! are accomodated by splitting the whole ensemble into five
! groups which can then be analyzed separately.

PROGRAM path_helper
USE parameters
USE prequel
USE lindblad_coherent
USE lindblad_jumps
USE rk_ns_utils
USE set_up_H

IMPLICIT NONE

! When do things collapse to eigenstates that are not zero? When a trajectory collapses
! within the first 'resolution' time, increment collapse_time_array(1)
INTEGER, DIMENSION(100000) :: collapse_time_array
REAL*8 :: resolution = 100.0d0  ! au time between each entry in "collapse_array"
INTEGER :: max_index  ! Keeping track of the maximum occupied entry of collapse_time_array
! When do trajectories end once and for all?
INTEGER, DIMENSION(100000) :: end_time_array
INTEGER :: max_end_index  ! end of that array
! Binning for how many trajectories have a certain number of jumps; if 1500 jumps ! are seen in a trajectory, increment jump_num_array(1500)
INTEGER, DIMENSION(10000) :: jump_num_array
! Resolution 10 on the number of jumps
INTEGER, DIMENSION(1000) :: jump_num_array_10
INTEGER :: running_jump_total  ! temporary to count up for any given trajectory
INTEGER :: max_jump_num_index  ! maximum number of jumps BETWEEN EIGENSTATES seen 

! index, indicator, previous eigenstate, id (unused), current eigenstate, total number
! of jumps, number of jumps between different states, number of trajectories read through
INTEGER :: i, success, prev_es, id, estate, num_jumps_seen, changes_seen, read_over
REAL*8 :: time
TYPE(wfunc) :: wfn
TYPE(jump_data) :: jump_info
! If one sees jump i-->j, increment entry (i,j), the second is incremented only if the eigenstate changes
INTEGER, DIMENSION(0:n,0:n) :: jumps_seen, jumps_seen_changes  
! The same idea as jumps_seen_changes but for block averaging standard deviation estimates (1/5 data in each)
INTEGER, DIMENSION(0:n,0:n) :: jsc_1, jsc_2, jsc_3, jsc_4, jsc_5  

INTEGER :: flag
INTEGER :: temp_int


CALL init_wfn(wfn)
OPEN(21,FILE="ph_jumps.txt")
num_jumps_seen = 0
jumps_seen_changes = 0
jsc_1 = 0
jsc_2 = 0
jsc_3 = 0
jsc_4 = 0
jsc_5 = 0
changes_seen = 0
jumps_seen = 0

end_time_array = 0
jump_num_array = 0
jump_num_array_10 = 0
running_jump_total = 0
max_jump_num_index = 0

collapse_time_array = 0
max_index = 0

DO i = 1, ntraj
  flag = 0
  success = read_dj_line(21,id,estate,time,wfn,jump_info)
  estate = identify_estate(wfn)  ! Might not be recorded properly on first entry of a traj file 
  prev_es = estate
  IF (success .NE. 0) THEN
    WRITE(*,*) "Read failed. Out of trajectories on error ", success, "trajectories read ", i - 1
    EXIT
  END IF
  DO
    success = read_dj_line(21,id,estate,time,wfn,jump_info)
    IF (estate .EQ. -1 .OR. success .NE. 0) THEN  ! Finished with that trajectory or error; start outer loop over
      
      temp_int = NINT(time/resolution) + 1
      IF (temp_int .GT. 100000) THEN
        temp_int = 100000
      END IF
      end_time_array(temp_int) = end_time_array(temp_int) + 1
      IF (temp_int .GT. max_end_index) THEN
        max_end_index = temp_int
      END IF
   
      IF (running_jump_total .GT. 10000) THEN
        running_jump_total = 10000
      END IF
      IF (running_jump_total .GT. max_jump_num_index) THEN
        max_jump_num_index = running_jump_total
      END IF
      jump_num_array(running_jump_total) = jump_num_array(running_jump_total) + 1
      temp_int = NINT(running_jump_total/10.0d0) + 1
      jump_num_array_10(temp_int) = jump_num_array_10(temp_int) + 1
      running_jump_total = 0
 
      EXIT
    END IF

    ! Did we just collapse to an eigenstate from an uncollapsed trajectory? If so, record some information
    ! about where it happens
    IF (estate .NE. 0 .AND. prev_es .EQ. 0) THEN
      temp_int = NINT(time/resolution) + 1
      IF (temp_int .GT. 100000) THEN
        temp_int = 100000
      END IF
      IF (temp_int .GT. max_index) THEN
        max_index = temp_int  ! Keep track of how much of this array is occupied
      END IF
      collapse_time_array(temp_int) = collapse_time_array(temp_int) + 1
    END IF


    IF (estate .EQ. prev_es) THEN ! Jumped into same state as started in
      jumps_seen(prev_es,estate) = jumps_seen(prev_es,estate) + 1
      num_jumps_seen = num_jumps_seen + 1
    ELSE  ! Changed identity by jumping
      running_jump_total = running_jump_total + 1

      ! Update sub jump matrices for an STDEV estimate later
      IF (MOD(i,5) .EQ. 0) THEN
        jsc_1(prev_es,estate) = jsc_1(prev_es,estate) + 1
      ELSE IF (MOD(i,5) .EQ. 1) THEN
        jsc_2(prev_es,estate) = jsc_2(prev_es,estate) + 1
      ELSE IF (MOD(i,5) .EQ. 2) THEN
        jsc_3(prev_es,estate) = jsc_3(prev_es,estate) + 1
      ELSE IF (MOD(i,5) .EQ. 3) THEN
        jsc_4(prev_es,estate) = jsc_4(prev_es,estate) + 1
      ELSE IF (MOD(i,5) .EQ. 4) THEN
        jsc_5(prev_es,estate) = jsc_5(prev_es,estate) + 1
      END IF

      jumps_seen_changes(prev_es,estate) = jumps_seen_changes(prev_es,estate) + 1
      jumps_seen(prev_es,estate) = jumps_seen(prev_es,estate) + 1
      num_jumps_seen = num_jumps_seen + 1
      changes_seen = changes_seen + 1
    END IF
    prev_es = estate
  END DO

END DO
CLOSE(21)

read_over = i - 1  ! How many trajectories were actually read in

WRITE(*,*) "Read over", i-1, "trajectories. Average jumps per trajectory", (1.0d0*num_jumps_seen)/(i-1), &
           "average number of state change jumps per trajectory", (1.0d0*changes_seen)/(i-1)


OPEN(25,FILE="jumps_markov_mat.txt")
WRITE(25,*) n+1
DO i = 1, n
  WRITE(25,*) jumps_seen_changes(i,1:n), jumps_seen_changes(i,0)
END DO
WRITE(25,*) jumps_seen_changes(0,1:n), jumps_seen_changes(0,0)
CLOSE(25)

! Sub matrices for analysis/looking at variation within the sample
! block averaging for standard deviation checks
OPEN(26,FILE="jumps_markov_mat_s1.txt")
WRITE(26,*) n+1
DO i = 1, n
  WRITE(26,*) jsc_1(i,1:n), jsc_1(i,0)
END DO
WRITE(26,*) jsc_1(0,1:n), jsc_1(0,0)
CLOSE(26)

OPEN(26,FILE="jumps_markov_mat_s2.txt")
WRITE(26,*) n+1
DO i = 1, n
  WRITE(26,*) jsc_2(i,1:n), jsc_2(i,0)
END DO
WRITE(26,*) jsc_2(0,1:n), jsc_2(0,0)
CLOSE(26)

OPEN(26,FILE="jumps_markov_mat_s3.txt")
WRITE(26,*) n+1
DO i = 1, n
  WRITE(26,*) jsc_3(i,1:n), jsc_3(i,0)
END DO
WRITE(26,*) jsc_3(0,1:n), jsc_3(0,0)
CLOSE(26)

OPEN(26,FILE="jumps_markov_mat_s4.txt")
WRITE(26,*) n+1
DO i = 1, n
  WRITE(26,*) jsc_4(i,1:n), jsc_4(i,0)
END DO
WRITE(26,*) jsc_4(0,1:n), jsc_4(0,0)
CLOSE(26)

OPEN(26,FILE="jumps_markov_mat_s5.txt")
WRITE(26,*) n+1
DO i = 1, n
  WRITE(26,*) jsc_5(i,1:n), jsc_5(i,0)
END DO
WRITE(26,*) jsc_5(0,1:n), jsc_5(0,0)
CLOSE(26)

! Write out the collapse time array for histogramming
OPEN(29,FILE="collapse_time.txt")
DO i = 1, max_index
  WRITE(29,*) resolution*(i-1), collapse_time_array(i), SUM(collapse_time_array(1:i))
END DO
CLOSE(29)

! Write out the ending time array
OPEN(29,FILE="end_time.txt")
DO i = 1, max_end_index
  WRITE(29,*) resolution*(i-1), end_time_array(i), SUM(end_time_array(1:i))
END DO
CLOSE(29)

! Write out the number of jumps array
WRITE(*,*) "Maximum number of inter-eigenstate jumps observed", max_jump_num_index
OPEN(29,FILE="num_jumps_histogram.txt")
temp_int = NINT(max_jump_num_index/10.0) + 1
DO i = 1, temp_int
  WRITE(29,*) (i-1)*10, i*10, jump_num_array_10(i), SUM(jump_num_array_10(1:i))
END DO
CLOSE(29)


END PROGRAM path_helper


