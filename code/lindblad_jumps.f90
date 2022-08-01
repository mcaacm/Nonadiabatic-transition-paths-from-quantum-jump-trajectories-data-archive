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

! Performs Lindblad jump dynamics. See lindblad_coherent.f90
! for the machinery it depends on.

! The implementation of Lindblad jumps here is largely based on 
! "Vogt, N., Jeske, J. and Cole, J. H. "Stochastic Bloch-Redfield theory:
! Quantum jumps in a solid-state environment" Phys. Rev. B 88, 2013.

MODULE lindblad_jumps
USE parameters
USE prequel
USE rk_ns_utils
USE lindblad_coherent
USE OMP_LIB
USE set_up_H
IMPLICIT NONE

CONTAINS

! Determine which quantum jump has occurred and 
! perform that jump. Requires the jump operators in 
! spaces, the wavefunction about to jump (wfn), a place
! to return the output wavefunction (wfno) and can take 
! the random number associated with the jump for help
! debugging should an error occur
SUBROUTINE perform_jump(d,wfn,wfno,spaces,time,prev_rn)
  INTEGER, INTENT (IN) :: d  ! Number of Lindblad subspaces
  TYPE(wfunc), INTENT(IN) :: wfn ! Input wavefuntion
  TYPE(wfunc), INTENT(INOUT) :: wfno  ! Output wavefunction (already allocated)
  TYPE(sspace), DIMENSION(d), INTENT(IN) :: spaces ! Lindblad subspaces
  REAL*8, INTENT(IN) :: prev_rn ! Random number associated with when to make this jump (for debugging)
  REAL*8, INTENT(IN) :: time  ! When the jump occurred
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: temp
  REAL*8, DIMENSION(d,nbath) :: probs ! Probabilities to jump to each other subspace
  REAL*8 :: pt, rn
  INTEGER :: i, j
  ALLOCATE(temp(n,n))
  ! Get the probabilities of all jumps which may occur
  CALL get_p_jumps(wfn,spaces,dt,d,probs)
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Jump probabilities (normalized) ", probs, " dimension ", d
  END IF
  ! Get a random number
  rn = get_rand_num()
  pt = 0.0d0
  ! Choose which jump to take
  DO i=1,d
    DO j = 1, nbath
      pt =  pt + probs(i,j)
      IF (DEBUG .EQV. .TRUE.) THEN
        WRITE(*,*) "Should I jump on i,j", i, j, "? Prob ", probs(i,j), "total now", probs(i,j)
      END IF
      IF (pt .GE. rn) THEN
        EXIT
      END IF
    END DO
    IF (pt .GE. rn) THEN
      EXIT
    END IF
  END DO 

  ! Error has occurred if this goes through; this can happen sometimes if the selected
  ! random number is zero to propogation precision. 
  IF (i .GT. d) THEN
    WRITE(*,*) "No jump selected! for wfn ", wfn%fe, "probs", probs
    WRITE(*,*) "Rn ", rn, "collapse time rn", prev_rn, " pt ", pt, " last operator ", i !, " is "
    WRITE(*,*) "The most likely problem is a bad system setup or unusually small random number seletion"
    WRITE(*,*) "causing the wavefunction norm to round to zero. Try reruning with a different three"
    WRITE(*,*) "digit random seed."
    STOP
  END IF
  CALL spam_times_vec(spaces(i)%Aomega,wfn%fe,wfno%fe,j)
  CALL normalize(wfno%fe)
  wfno%sn = identify_estate(wfno) 
  ! It would be better to pass in the file to output to, but just knowing that 42 is always
  ! the output file for the jump information works, too. Output the whole wavefunction only
  ! if the wavefunction is uncollapsed, else record the eigenstate and one relevant entry
  IF (wfno%sn .EQ. 0) THEN
    WRITE(42,*) 0, time, wfno%sn 
    WRITE(42,*) org_out(wfno%fe) 
  ELSE
    WRITE(42,*) 0, time, wfno%sn 
    WRITE(42,*) REAL(REAL(wfno%fe(wfno%sn,1))), AIMAG(wfno%fe(wfno%sn,1))
  END IF
  DEALLOCATE(temp)
END SUBROUTINE


! Get the probabilities of all the jumps from operators in spaces,
! assumes the wavefunction is normalized.
! Returns a normalized vector of probabilities to make the jump
! associated with each possible jump operator
SUBROUTINE get_p_jumps(wfn,spaces,dt,d,probs)
  TYPE(wfunc), INTENT(IN) :: wfn ! Input wavefunction
  REAL*8, DIMENSION(d,nbath), INTENT(OUT) :: probs  ! Jump probabilities to output
  INTEGER, INTENT(IN) :: d  ! Number of Lindblad spaces
  REAL*8, INTENT(IN) :: dt  ! Time step
  REAL*8 :: total ! For normalization
  TYPE(sspace), DIMENSION(d), INTENT(IN) :: spaces ! Lindblad subspaces
  TYPE(wfunc) :: nwfn
  INTEGER i, j, jpn
  total = 0.0d0
  CALL init_wfn(nwfn)  ! That wfn is read only so have a copy with its eigenstate attached
  nwfn%fe = wfn%fe
  nwfn%sn = identify_estate(wfn)
  probs = 0.0d0
  
  IF (nwfn%sn .EQ. 0) THEN  ! Uncollapsed treatment; all operators are relevant
    DO i = 1, d
      DO j = 1, nbath
        probs(i,j) = get_pa(nwfn,spaces(i)%thpl(j),spaces(i)%Aomega,dt,j)
        total = probs(i,j) + total
      END DO
    END DO 
  ELSE  ! Collapsed treatment; look at the organize_operators indexes
    i = 1
    jpn = jump_pointers(nwfn%sn,i)  ! Index of next jump operator that may work
    DO WHILE (jpn .NE. -1)
      IF (DEBUG .EQV. .TRUE.) THEN
        WRITE(*,*) "Selecting jump pointer ", jpn
      END IF
      DO j = 1, nbath
        probs(jpn,j) = get_pa(nwfn,spaces(jpn)%thpl(j),spaces(jpn)%Aomega,dt,j)
        IF (DEBUG .EQV. .TRUE.) THEN
          WRITE(*,*) "Adding probability from operator ", jpn, "of", probs(jpn,j)
        END IF
        total = probs(jpn,j) + total
      END DO
      i = i + 1
      jpn = jump_pointers(nwfn%sn,i)  ! Update jump pointer index
    END DO 
  END IF

  probs = probs/total   
END SUBROUTINE

! Given a wavefunction and a jump operator and
! associated jump rate, this is the 
! probability of jump operator lbo occuring
REAL*8 FUNCTION get_pa(wfn,thpl,lbo,dt,j)
  INTEGER, INTENT(IN) :: j
  TYPE(wfunc), INTENT(IN) :: wfn  ! Input wavefunction
  COMPLEX*16, INTENT(IN) :: thpl  ! Lindblad rate
  REAL*8, INTENT(IN) :: dt ! Time step
  TYPE(SPAM), INTENT(IN) :: lbo  ! Proposed jump operator
  COMPLEX*16, DIMENSION(n,1) :: rvec  ! A right vector
  COMPLEX*16 :: res  ! Temporary result
  REAL*8 :: norm  ! norm of sparce multiplication result

  IF (lbo%nentries .EQ. 1 .AND. wfn%sn .NE. 0) THEN
    ! Exit early since this obviously will not have any density
    IF (lbo%cols(1) .NE. wfn%sn) THEN
      get_pa = 0.0d0
      RETURN
    END IF
  END IF
  ! Perform the sparce matrix/vector multiplication 
  CALL spam_times_vec(lbo,wfn%fe,rvec,j,norm)
  res = norm*thpl*dt/(hbar**2.0d0)
  get_pa = res
END FUNCTION

! Driver function for perparing the lindblad operators and then
! running propagation from the initial conditions
! Translated_evecs are the wavefunctions to start from with
! associated start probability estate_probs. This is no longer
! relevant but back when true thermal initialization was an 
! option multiple start vectors were used.
SUBROUTINE drive_lindblad(Hs,evals,evecs_inv,translated_evecs,estate_probs)
  ! Eigenvectors and inverses, diagonalized Hamiltonian, position operator, 
  ! the vectors used to initialize the jumps
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: Hs, translated_evecs
  COMPLEX*16, DIMENSION(n,m), INTENT(IN) :: evecs_inv
  ! Probability to initialize from each jump starting vector
  REAL*8, INTENT(IN), DIMENSION(n) :: estate_probs
  ! Eigenvalues of the hamiltonian
  COMPLEX*16, DIMENSION(n), INTENT(IN) :: evals
  ! Lindblad jump operator subspaces
  TYPE(sspace), DIMENSION(:), ALLOCATABLE :: spaces 
  ! Number of Lindblad jump operators
  INTEGER :: d
  ! Lambshift Hamiltonian
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: Hls
  ! Real type for eigenvalues
  REAL*8, DIMENSION(n) :: revals

  ALLOCATE(Hls(n,n))

  OPEN(21,file="rkevolution.txt")
  OPEN(25,file="p_in_b_data.txt")
  revals = REAL(REAL(evals))
  WRITE(*,*) "Sorting Lindblad operators"
  CALL new_sort(revals,spaces,d) 
  WRITE(*,*) "Sorted, dimension of Lindblad operator space is : ", d
  CALL make_sparce_operator(d,spaces,G_global)
  WRITE(*,*) "Made sparce matrices for jump calculations."
  CALL get_Hls(d,spaces,Hls)
  WRITE(*,*) "Got Lamb Shift H. Organizing operators"
  CALL organize_operators(d, spaces)
  ! Call the jump driver routine.  
  CALL drive_lb_jumps(d,Hs,evecs_inv,spaces,translated_evecs,estate_probs,Hls)

  CLOSE(21)
  CLOSE(25)
  DEALLOCATE(Hls)
END SUBROUTINE


! Propogate a wavefunction forward in time from the
! coherent Hamiltonian for the dynamics.
! hco MUST be diagonal in the basis passed in. 
! This should always be the case under the strict
! secular approximation assumptions made in which
! all subspaces are presumed to have only one entry.
SUBROUTINE rk4_coherent(hco,wfn,wfno,dt)
  ! Coherent Hamiltonian
  COMPLEX*16, DIMENSION(n), INTENT(IN) :: hco 
  ! Wavefunction in
  TYPE(wfunc), INTENT(IN) :: wfn
  ! Wavefunction out
  TYPE(wfunc), INTENT(OUT) :: wfno
  ! Time step
  REAL*8, INTENT(IN) :: dt
  INTEGER :: i
  DO i = 1, n
    wfno%fe(i,1) = EXP(CMPLX(0.0d0,-dt)*hco(i))*wfn%fe(i,1)
  END DO
END SUBROUTINE

! Use a bisection search to find an exact time of jump
! return the time when the jump occurs and the value of the
! wavefuction (wfi) where the jump occurs
! dt is the width of the interval at start.
! wfl is the start of the interval, wfi is the end
! they are separated by dt in time and t is the time
! at which wfi exists (end time of the interval)
! Note that changing "dt" can change the exact time at which
! the jump occurs, biasing high or low, but usually not by enough to affect
! results significantly.
SUBROUTINE find_jump_time(hco,hco_ev,hco_evin,p,wfl,wfi,dt,t)
  COMPLEX*16, DIMENSION(n), INTENT(IN) :: hco  ! Coherent hamiltonian
  ! Wavefunction at middle (wfi) and start (wfl) of search interval
  TYPE(wfunc), INTENT(INOUT) :: wfi, wfl
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: hco_ev, hco_evin  ! Hco eigenvectors
  ! Timestep, search tolerance
  REAL*8, INTENT(IN) :: dt, p
  REAL*8, INTENT(INOUT) :: t ! current time
  REAL*8 :: cdt  ! Time at which a jump occurred
  ! Ellapsed iterations 
  INTEGER :: iter
  cdt = dt
  iter = 0
  IF (non_diag_hco .EQ. 1) THEN
    wfi%fe = MATMUL(hco_ev,wfi%fe)
  END IF
  DO WHILE (ABS((normsq(wfi%fe) - p)) .GE. tol .AND. iter .LT. 50) 
    ! If it is not, prepare to take a 1/2 dt step 
    cdt = cdt/2.0d0
    ! See whether the norm is too high or too low
    ! if the norm is too low, use the left point and take a 1/2 dt step forward
    IF ((normsq(wfi%fe) - p) .LT. 0) THEN ! The norm is too low, back up, half step
      t = t - cdt
      CALL rk4_coherent(hco,wfl,wfi,cdt) 
      IF (non_diag_hco .EQ. 1) THEN
        wfi%fe = MATMUL(hco_ev,wfi%fe)
      END IF
    ! otherwise, use the middle point and take a half step forward
    ELSE ! Norm is too high; take half step forward
      t = t + cdt
      IF (non_diag_hco .EQ. 1) THEN
        wfi%fe = MATMUL(hco_evin,wfi%fe)
      END IF
      wfl = wfi
      CALL rk4_coherent(hco,wfl,wfi,cdt) 
      IF (non_diag_hco .EQ. 1) THEN
        wfi%fe = MATMUL(hco_ev,wfi%fe)
      END IF
    END IF
    iter = iter + 1
  END DO
  IF (non_diag_hco .EQ. 1) THEN
    wfi%fe = MATMUL(hco_evin,wfi%fe)  ! Used by previous program, pass back in hco basis
  END IF
  IF (ABS(normsq(wfi%fe)) .LT. 1.0d-15) THEN
    WRITE(*,*) "WARNING: Wavefunction norm following jump time search extremely small: ", &
               ABS(normsq(wfi%fe)), "wfn", wfi%fe
  END IF

  ! Usually results in a crash but not always; unusually small random numbers can
  ! sometimes cause this, or badly behaved dynamics
  IF (iter .GE. 50) THEN
    WRITE(*,*) "Iterations exceeded norms ", p, normsq(MATMUL(hco_ev,wfi%fe))
    WRITE(*,*) "hco ", hco
    WRITE(*,*) "wfi ", wfi%fe
  END IF
END SUBROUTINE

! Drive the jump version of Lindblad dynamics, calling trajectory
! driving functions the desired number of times
SUBROUTINE drive_lb_jumps(d,Hs,evec_inv,spaces,translated_evecs,estate_probs,Hls)
  INTEGER, INTENT(IN) :: d ! Number of Lindblad supscaces
  ! System hamiltonian, starting vectors (used to allow more than possibility),
  ! lamb shift hamiltonian
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: Hs, translated_evecs, Hls
  COMPLEX*16, DIMENSION(n,m), INTENT(IN) :: evec_inv  ! Used to be used for printing; deprecated
  ! Lindblad subspaces
  TYPE(sspace), DIMENSION(d), INTENT(IN) :: spaces 
  ! Probabilities to use each translated_evecs as a starting point
  REAL*8, INTENT(IN), DIMENSION(n) :: estate_probs
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: temp  ! Temporary matrix
  ! Coherent Hamiltonian, its eigenvectors, and eigenvector inverse
  ! which used to be important when allowing less stric secular approximations
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: hco, hco_ev, hco_evin
  ! Coherent Hamitlonian eigenvalues
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: hco_evals
  INTEGER :: i

  ALLOCATE(temp(n,n)) 
  ALLOCATE(hco(n,n))
  ALLOCATE(hco_evin(n,n))
  ALLOCATE(hco_ev(n,n))
  ALLOCATE(hco_evals(n))
  OPEN(42,FILE="jump_data.txt")

  CALL init_rand(get_seed())
  CALL get_hco(d,Hs,spaces,hco,Hls)

  ! Diagonalizing Hco should never be necessary anymore; this is 
  ! a holdover from a more complicated software version.
  ! Diagonalize Hco and record left and right eigenvectors
  ! if necessary; if not avoid potential reshuffling of eigenvector
  ! order and any associated problems by just copying in directly;
  ! make the eigenvectors identity incase they get used by accident
  IF (non_diag_hco .EQ. 1) THEN
    CALL diagonalize_g(hco,hco_evals,hco_ev)
    hco_evin = hco_ev
    CALL invert(hco_evin)
  ELSE
    hco_ev = 0.0d0
    hco_evin = 0.0d0
    DO i = 1, n
      hco_ev(i,i) = 1.0d0
      hco_evin(i,i) = 1.0d0  
      hco_evals(i) = hco(i,i)
    END DO
  END IF
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Diagonalizing Hco; evals are ", hco_evals
    WRITE(*,*) "hco_ev (Hco eigenvectors)"
    CALL print_mat(hco_ev)
    WRITE(*,*) "hco_evin (inverse of Hco eigenvectors)"
    CALL print_mat(hco_evin)
  END IF
  

  DO i = 1, ntraj
     WRITE(*,*) "Starting trajectory ", i
     CALL drive_trajectory_early_end(d,spaces,hco_evals,translated_evecs,estate_probs, &
       hco_ev,hco_evin)
     WRITE(*,*) "Ending trajectory ", i
  END DO

  CLOSE(42)
END SUBROUTINE


! Ends only when a collapse to one of two given eigenstates is recorded
! Drives single lindblad jumps trajectory
! Needs access to the translated eigenvectors from the first set up
! (translated_evecs) which are already in the eigenbasis and specify
! states from which to begin wavefunction propagation
! Conversions to and from the Hco eigenbasis remain despite the fact that
! all the cases in which Hco might not be diagonal have been removed. They
! might be reimplemented later.
SUBROUTINE drive_trajectory_early_end(d,spaces,hco,translated_evecs,estate_probs,&
                            hco_ev,hco_evin)
  TYPE(sspace), DIMENSION(d), INTENT(IN) :: spaces ! Lindblad subspaces
  ! inverse of eigenvectors, Q in proper basis, starting eigenvectors
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) ::  translated_evecs
  ! Coherent hamiltonian and its eigenvectors
  COMPLEX*16, DIMENSION(n), INTENT(IN) :: hco
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: hco_ev, hco_evin  ! Hco eigenvectors (Hco could have be non-diagonal in early cases)
  ! Probability to start in each of the translated_evecs respectively
  REAL*8, INTENT(IN), DIMENSION(n) :: estate_probs
  INTEGER, INTENT(IN) :: d
  ! Wavefunctions for various purposes before/after evolution
  TYPE(wfunc) :: wfn, res, nwfn
  REAL*8 :: ct, jt, rn, acc, gap ! Current time, jump time, random number, prob. accumulation, temp for multiple jumps/timestep
  INTEGER i, temp_count

  rn = get_rand_num()  ! Historical reasons; there used to be a need for this here and
  ! if it's removed it will mess up the ensembles
  ct = 0.0d0
  acc = 0.0d0
  rn = get_rand_num()

  ! Iterate through the start options and, choose the state
  ! to start this trajectory in
  DO i = 1, n/2
    acc = acc + estate_probs(i) 
    IF (acc .GE. rn) THEN
      EXIT
    END IF
  END DO
  CALL init_wfn(res)
  CALL init_wfn(nwfn)

  ! Initialize wavefunction with proper entries filled
  CALL init_wfn(wfn)
  wfn%fe = 0.0d0
  wfn%fe(1:n,1) = translated_evecs(1:n,i)
  CALL normalize(wfn%fe)  ! For safety

  ! Write init data to jump file
  WRITE(42,*) 0, 0.0d0, 0
  WRITE(42,*) org_out(wfn%fe)

  ! Hco has been diagonalized in order to propogate forward
  ! in time; transfer the wavefunction to the diagonal HCO
  ! basis for propogation
  IF (non_diag_hco .EQ. 1) THEN
    wfn%fe = MATMUL(hco_evin,wfn%fe)
  END IF
  res = wfn
  rn = get_rand_num()

  ct = 0.0d0
  DO 
    CALL rk4_coherent(hco,wfn,res,dt)
    IF (non_diag_hco .EQ. 1) THEN
      nwfn%fe = MATMUL(hco_ev,res%fe)
    ELSE
      nwfn%fe = res%fe
    END IF
    ! Time to perform a quantum jump
    gap = dt
    temp_count = 0
    
    ! Jump has occurred
    DO WHILE (normsq(nwfn%fe) .LE. rn)  ! Has been basis converted iff necessary
      jt = ct + dt
      CALL find_jump_time(hco,hco_ev,hco_evin,rn,wfn,res,gap,jt)
      wfn = res
      ! Move WFN out of the basis in which Hco is diagonal in order
      ! to complete a jump
      IF (non_diag_hco .EQ. 1) THEN
        wfn%fe = MATMUL(hco_ev,wfn%fe)
      END IF
      CALL normalize(wfn%fe)
      CALL perform_jump(d,wfn,res,spaces,jt,rn)
      temp_count = temp_count + 1


      wfn = res
      CALL normalize(wfn%fe)  ! For safety's sake
      wfn%sn = identify_estate(wfn) 
      ! The wavefunction has reached one of the ending states. End propogation
      IF (wfn%sn .EQ. flag_estate_1 .OR. wfn%sn .EQ. flag_estate_2) THEN 
        CALL destroy_wfn(wfn)
        CALL destroy_wfn(res)
        CALL destroy_wfn(nwfn)
        WRITE(42,*) 0, jt, -1
        RETURN
      END IF

      ! Put res back in the basis in which Hco is diagonal to propogate
      IF (non_diag_hco .EQ. 1) THEN
        wfn%fe = MATMUL(hco_evin,wfn%fe)
      END IF
      rn = get_rand_num()

      ! After performing the jump, propagate forward
      gap = ct + gap - jt
      CALL rk4_coherent(hco,wfn,res,gap)

      IF (non_diag_hco .EQ. 1) THEN
        nwfn%fe = MATMUL(hco_ev,res%fe)
      ELSE
        nwfn%fe = res%fe
      END IF

    END DO
    ! Move forward the amount of a full dt 
    ct = ct + dt
    wfn = res
  END DO

END SUBROUTINE

END MODULE lindblad_jumps
