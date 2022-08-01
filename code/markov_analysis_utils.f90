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

! Useful functions for analyzing a Markov model to try to find
! transition paths, committor, etc
MODULE analyze_markov_utils

! Equations are implemented from 
! Noe, F., Schutte, C., Vanden-Eijnden, E. and Weikl, T. R. 
! "Constructing the equilibrium ensemble of folding pathways 
! from short off-equilibrium simulations", PNAS 106, 2009.

IMPLICIT NONE

! Print lots of extra information
LOGICAL, PARAMETER :: DEBUG = .FALSE.
! Reverse committors can be done without detailed balance imposed (.FALSE.)
! or with it (.TRUE.). In detailed balance systems either should be equivalent.
! Testing that results are invariatnt to this is one way I set my mind at ease
LOGICAL, PARAMETER :: detailed_balance = .FALSE.  

! Eigenstate B (the product eigenstate, bottom of the higher energy well usually)
INTEGER, PARAMETER :: B_estate = 3

CONTAINS

! Quick sort an array of dimension d with indexes returned in 'ind'
! which show which entry has ended up where in the list
! low and high are where to sort on the array (usually 1, d)
! Returns the new order in 'ind' and the sorted values in 'arr'
RECURSIVE SUBROUTINE quick_sort(d,arr,ind,low,high)
  INTEGER, INTENT(IN) :: d, low, high
  REAL*8, DIMENSION(d), INTENT(INOUT) :: arr
  INTEGER, DIMENSION(d),  INTENT(INOUT) :: ind
  INTEGER :: pi

  IF (low .LT. high) THEN
    pi = partition(d,arr,ind,low,high)
    CALL quick_sort(d,arr,ind,low,pi-1)
    CALL quick_sort(d,arr,ind,pi+1,high)
  END IF
END SUBROUTINE

! Quick sort helper; arr is values, ind is their associated indices
INTEGER FUNCTION partition(d,arr,ind,low,high)
  INTEGER, INTENT(IN) :: d, high, low ! array dimension, sort bounds
  REAL*8, DIMENSION(d), INTENT(INOUT) :: arr
  INTEGER, DIMENSION(d),  INTENT(INOUT) :: ind
  REAL*8 :: pivot, tr  ! pivot, temp
  INTEGER :: i, j, ti  ! two indices and a temp

  i = low - 1
  pivot = arr(high)
  
  DO j = low, high -1 
    IF (arr(j) .LT. pivot) THEN
      i = i + 1
      tr = arr(i) 
      ti = ind(i)
      arr(i) = arr(j)
      ind(i) = ind(j)
      arr(j) = tr
      ind(j) = ti
    END IF
  END DO

  tr = arr(i+1)
  ti = ind(i+1)
  arr(i+1) = arr(high)
  ind(i+1) = ind(high)
  arr(high) = tr
  ind(high) = ti

  partition = i + 1  ! Return the boundary
END FUNCTION

! Pass in 'vals' as a set of percentages. They will be sorted and then
! passed back in running summation form
SUBROUTINE running_sum(d,vals)
  INTEGER, INTENT(IN) :: d
  REAL*8, DIMENSION(d), INTENT(INOUT) :: vals
  INTEGER, DIMENSION(d) :: ind
  REAL*8, DIMENSION(d) :: temp
  INTEGER :: low, high

  DO high = 1, d
    ind(high) = high
  END DO
  high = d
  low = 1
  CALL quick_sort(d,vals,ind,low,high)
  temp = vals

  ! Running sum conversion
  vals(1) = vals(d)
  DO high = 2, d
    vals(high) = temp(d - high + 1) + vals(high - 1)
  END DO
END SUBROUTINE

! Quick Dijkstra widest path implementation; uses slow O(n) minsearch on
! each implementation rather than a fast priority queue because it probably
! doesn't matter much for these small graphs and the more complicated
! method is easier to mess up.
! Can work on general connection matrices. Pass graph in as V. Widths
! and backtracking information in A_width and A_previous are passed in
! and overwritten with information to follow back along the path
SUBROUTINE dijkstra_wp(d,V,A_width,A_previous,init_node)
  INTEGER, INTENT(IN) :: d, init_node ! Graph dimension, the node where we start
  REAL*8, INTENT(IN), DIMENSION(d,d) :: V  ! Edge connection and weights
  REAL*8, INTENT(OUT), DIMENSION(d) :: A_width  ! Width to you
  INTEGER, INTENT(OUT), DIMENSION(d) :: A_previous  ! Backtrack information
  LOGICAL, DIMENSION(d) :: A_present  ! List of those present in working set
  REAL*8, DIMENSION(d) :: B_width  ! Auxilliary weidth
  REAL*8 :: temp_width
  INTEGER :: i, x, num_A

  num_A = d
  A_width = -1.0d0  ! Flag value; unreachable
  B_width = -1.0d0  ! Auxiliary width info for those removed from queue
  A_present = .TRUE.  ! Yes, all nodes are present in the queue at start
  A_previous = -1
  A_width(init_node) = MAXVAL(V)*2.0d0  ! Thou shalt not be a bottle neck
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "init_node is", init_node
  END IF

  DO WHILE (num_A .NE. 0)
    x = MAXLOC(A_width,1)
    IF (DEBUG .EQV. .TRUE.) THEN 
      WRITE(*,*) "MAXLOC found", x, "width", A_width(x), "from", A_width
    END IF
    A_present(x) = .FALSE.  ! Removed from queue
    B_width(x) = A_width(x)
    num_A = num_A - 1
    
    DO i = 1, d
      ! There is a connection between the nodes in the weight matrix and destination
      ! node is still present in the queue
      IF (V(x,i) .GT. 1.0d-15 .AND. A_present(i) .EQV. .TRUE.) THEN
        temp_width = MIN(A_width(x),V(x,i))  ! Maximum flux allowed from x to i
        IF (DEBUG .EQV. .TRUE.) THEN
          WRITE(*,*) "temp_width ", temp_width, "A_width(x)", A_width(x), "V(x,i)", V(x,i)
        END IF
        IF (temp_width .GT. A_width(i)) THEN
          IF (DEBUG .EQV. .TRUE.) THEN
            WRITE(*,*) "temp_width", temp_width, "A_width", A_width(i), "i =", i, "V", x, i, "=", V(x,i)
          END IF
          A_previous(i) = x
          A_width(i) = temp_width 
        END IF
      END IF
    END DO
    A_width(x) = -2.0d0  ! Will no longer be found by MAXLOC
  END DO
  A_width = B_width
END SUBROUTINE

! Backtrack along a path from "sink" to "source" and subtract that
! flux from V so that it can no longer be used in computing a new
! widest path; also record the path and the flux. Return -1 on failure,
! else the path length
INTEGER FUNCTION back_track_sub(d,V,sink,source,A_width,A_previous,path)
  INTEGER, INTENT(IN) :: d, sink, source
  REAL*8, INTENT(INOUT), DIMENSION(d,d) :: V ! Flux matrix to have the path removed
  INTEGER, INTENT(OUT), DIMENSION(d) :: path  ! traversed path
  REAL*8, INTENT(IN), DIMENSION(d) :: A_width  ! Width to you
  INTEGER, INTENT(IN), DIMENSION(d) :: A_previous  ! Backtrack information
  INTEGER :: i, j  ! Where we are (i) where we are headed back (j)
  REAL*8  :: flux
  back_track_sub = 1

  flux = A_width(sink)
  i = sink
  j = A_previous(sink)
  path = source  ! Convenient for plotting
  path(1) = sink
  ! Backtrack until we reach the source or we can't anymore
  DO WHILE (i .NE. source .AND. j .NE. -1) 
    V(j,i) = V(j,i) - flux
    IF (V(j,i) .LT. -1.0d-15) THEN
      WRITE(*,*) "WARNING: negative flux now present in V during backtrack ", j, i, "value is", V(j,i)
    END IF 
    i = j
    j = A_previous(i)
    back_track_sub = back_track_sub + 1
    path(back_track_sub) = i
  END DO

  IF (i .NE. source) THEN  ! Success
    back_track_sub = -1  ! Failure
  END IF
END FUNCTION

! Print a matrix row by row
SUBROUTINE print_mat_r8(M,fnum)
  REAL*8, DIMENSION(:,:), INTENT(IN) :: M
  INTEGER, OPTIONAL, INTENT(IN) :: fnum
  INTEGER i
  INTEGER, DIMENSION(2) :: info
  info = SHAPE(M)
  IF (PRESENT(fnum)) THEN
    DO i=1,info(1)
      WRITE(fnum,*) M(i,:)
    END DO
  ELSE
    DO i=1,info(1)
      WRITE(*,*) M(i,:)
    END DO
  END IF
END SUBROUTINE

! Solve a double precision A*X=B problem via Lapack
! I think this has more variables than it really needed
! but this isn't Fortran 66 so it's fine.
SUBROUTINE solve_system(d,A,B,X)
  INTEGER, INTENT(IN) :: d  ! Dimension 
  REAL*8, DIMENSION(d,d), INTENT(IN) :: A
  REAL*8, DIMENSION(d), INTENT(IN) :: B
  REAL*8, DIMENSION(d), INTENT(OUT) :: X
  INTEGER, DIMENSION(d) :: ipiv
  REAL*8, DIMENSION(d) :: r, c
  REAL*8, DIMENSION(d,1) :: B_temp, X_temp
  REAL*8, DIMENSION(4*d) :: work
  REAL*8, DIMENSION(1) :: berr, ferr
  REAL*8 :: rcond
  INTEGER, DIMENSION(d) :: iwork
  INTEGER :: info, i
  REAL*8, DIMENSION(d,d) :: A_mod, Af
  CHARACTER(1) :: equed

  Af = A
  A_mod = A
  X = B
  B_temp(1:d,1) = B(1:d)
  X_temp(1:d,1) = B(1:d)
  ipiv = 1
  i = 1
  CALL dgesvx('E','N',d,i,A_mod,d,Af,d,ipiv,equed,r,c,B_temp,d,X_temp,d,rcond, &
               ferr,berr,work,iwork,info)
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Rcond", rcond, "equed", equed
  END IF
  IF (info .NE. 0) THEN
    WRITE(*,*) "System not successfully solved"
  END IF
  X(1:d) = X_temp(1:d,1)
END SUBROUTINE

! Read an integer matrix graph of a given size, i.e. entry
! i,j of the graph says how many jumps were observed from eigenstate
! i to eigenstate j over the course of an ensemble
! The first entry read in is the size of the graph; entry d is 
! assumed to be the source.
SUBROUTINE read_graph(graph,d,fnum)
  INTEGER, INTENT(OUT) :: d  ! graph dimension
  INTEGER, INTENT(IN) :: fnum  ! File to read from
  INTEGER, INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: graph  ! Graph read in
  INTEGER :: reason, i, j

  READ(fnum,*,IOSTAT=reason) d

  ! Well, that's a problem
  IF (reason .NE. 0 .OR. d .LE. 0) THEN
    WRITE(*,*) "Error on MSM specifications, iostat, dimension or timescale", reason, d
    STOP
  END IF

  ALLOCATE(graph(d,d))
  DO i = 1, d
    READ(fnum,*,IOSTAT=reason) (graph(i,j), j=1, d)
    IF (reason .NE. 0) THEN
      WRITE(*,*) "Failed read in graph on ", i, "reason", reason, "graph", graph(i,1:d)
      STOP
    END IF
  END DO 

  WRITE(*,*) "Read graph."
END SUBROUTINE

! Read the eigenvector information from a file 
! Format: pos1, pos2, energy, probability to be in diabatic state 1, probability to be in diabatic state 2
SUBROUTINE read_evi(d,fnum,even,evpc,evpt,ps1,ps2)
  INTEGER, INTENT(IN) :: d, fnum
  REAL*8, DIMENSION(d), INTENT(OUT) :: even, evpc, evpt, ps1, ps2  ! Eigenvector energy, position, elec state character
  REAL*8 :: pos1, pos2, energy, pin1, pin2
  INTEGER :: i, success

  DO i = 1, d
    READ(fnum,*,IOSTAT=success) pos1, pos2, energy, pin1, pin2
    even(i) = energy
    evpc(i) = pos1
    evpt(i) = pos2
    ps1(i) = pin1
    ps2(i) = pin2
    IF (success .NE. 0) THEN  ! The program can continue without this although some numbers will be garbage
      WRITE(*,*) "Failed to read all eigenvector information", success
    END IF
  END DO
END SUBROUTINE

! Read a Markov State Model out of a file in the form
! n, dt 
! equilibrium probabilities (n of them)
! lines of numbers, n per line, each line normalized
! The file specified by fnum is presumed to already be open
SUBROUTINE read_msm(msm,pis,neq_pis,cprobs,evals,d,dt,fnum)
  REAL*8, INTENT(OUT) :: dt  ! tau
  INTEGER, INTENT(OUT) :: d  ! Dimension of model
  INTEGER, INTENT(IN) :: fnum  ! Open file
  REAL*8, INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: msm  ! The Markov state model
  ! neq_pis and cprobs are relics from an idea that did not work
  ! evals are the energies, pis are the equilibrium probabilities
  REAL*8, INTENT(OUT), DIMENSION(:), ALLOCATABLE :: pis, neq_pis, cprobs, evals  
  INTEGER :: reason, i, j
  REAL*8 :: temp

  READ(fnum,*,IOSTAT=reason) d, dt

  ! Well, that's a problem
  IF (reason .NE. 0 .OR. d .LE. 0 .OR. dt .LT. 0) THEN
    WRITE(*,*) "Error on MSM specifications, iostat, dimension or timescale", reason, d, dt
    STOP
  END IF

  ALLOCATE(msm(d,d))
  ALLOCATE(pis(d))
  pis = 0.0d0
  READ(fnum,*,IOSTAT=reason) (pis(j), j=1, d)
  IF (reason .NE. 0) THEN
    WRITE(*,*) "Failed read in MSM", reason, "read", pis
    STOP
  END IF

  DO i = 1, d
    READ(fnum,*,IOSTAT=reason) (msm(i,j), j=1, d)
    IF (reason .NE. 0) THEN
      WRITE(*,*) "Failed read in MSM on ", i, "reason", reason, "msm", msm(i,1:d)
      STOP
    END IF
  END DO 

  ! Read in the eigenvalues
  ALLOCATE(evals(d)) ! The uncollapsed state can stay uncollapsed, so it has n+1 options
  evals = 0.0d0
  READ(fnum,*,IOSTAT=reason) (evals(j), j = 1, d) 
  IF (reason .NE. 0) THEN
    WRITE(*,*) "Failed read in evals", reason, "read", evals
  END IF

  ! Check that the MSM is normalized
  DO i = 1, d
    temp = 0.0d0
    DO j = 1, d
      temp = temp + msm(i,j)
    END DO
    IF (ABS(temp - 1.0d0) .GT. 1.0d-5) THEN
      WRITE(*,*) "WARNING: Normalization doesn't hold on line", i, "value", temp  
    END IF
  END DO
  WRITE(*,*) "Checked MSM. Normalized."

  ! Check that detailed balance is maintained
  WRITE(*,*) "Checking detailed balance"
  DO i = 1, d
    DO j = 1, d
      IF (ABS(msm(i,j)*pis(i)/pis(j) - msm(j,i)) .GT. 1.0d-5) THEN
        WRITE(*,*) "Detailed balance not holding on ", i, j, "val", msm(i,j)*pis(i)/pis(j), "vs", &
                   msm(j,i)
      END IF
    END DO
  END DO
END SUBROUTINE

! Analysis of MSM. Calculate the committors according to the formula presented in
! www.pnas.org/cgi/doi/10.1073/pnas.0905466106
! i_def, a_def, and b_def contain the indexes of entries in I, B, and A
! They should be in increasing order but it will probably work if they aren't.
SUBROUTINE calc_comm(d,msm,comm_f,comm_b,ni,nb,na,i_def,a_def,b_def,pis)
  REAL*8, DIMENSION(d,d), INTENT(IN) :: msm ! Markov state model
  REAL*8, DIMENSION(d), INTENT(OUT) :: comm_f, comm_b  ! Forward committor, backward committor
  REAL*8, DIMENSION(d), INTENT(IN) :: pis  ! Equilibrium probabilities
  INTEGER, DIMENSION(ni), INTENT(IN) :: i_def  ! Indexes not in reactants or products
  INTEGER, DIMENSION(na), INTENT(IN) :: a_def  ! Reactant indexes
  INTEGER, DIMENSION(nb), INTENT(IN) :: b_def  ! Product indexes
  INTEGER, INTENT(IN) :: d, ni, nb, na  ! Dimension of n, Numbers of eigenstate in I, A, B
  INTEGER :: i, j
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: comm_mat
  REAL*8, DIMENSION(ni) :: B_mat, solu
  REAL*8, DIMENSION(ni,1) :: solu_mat, B_check
  REAL*8 :: check_1, check_2

  ALLOCATE(comm_mat(ni,ni))
  ! \hat{T} = T_i,j, if i in I and j != i
  ! \hat{T} = T_i,j - 1 if i in I j = i
  ! Build the committor matrix to solve the system of equations for the intermediate state
  ! A in A*X = B system to be solved
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "I def ", i_def, "a def", a_def, "b def", b_def
  END IF

  comm_mat = 0.0d0
  comm_f = 0.0d0
  comm_b = 0.0d0
  B_mat = 0.0d0
  solu = 0.0d0
  solu_mat = 0.0d0
  B_check = 0.0d0

  DO i = 1, ni
    DO j = 1, ni
      ! Indexes of the I matrix extracted form the array of its definition
      IF (i .EQ. j) THEN
        comm_mat(i,j) = msm(i_def(i),i_def(j)) - 1.0d0
      ELSE
        comm_mat(i,j) = msm(i_def(i),i_def(j))
      END IF
    END DO
  END DO
  ! Assemble the B of the A*X = B system; B_i = summation over -T_{i,j} where j is in B
  B_mat = 0.0d0
  DO i = 1, ni
    DO j = 1, nb
      B_mat(i) = B_mat(i) - msm(i_def(i),b_def(j)) 
    END DO
  END DO 
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Assembled B_mat", B_mat
    WRITE(*,*) "Comm_mat "
    CALL print_mat_r8(comm_mat,6)
    WRITE(*,*) "Sending to solve"
  END IF

  ! Solve the system to get the forward committors.
  CALL solve_system(ni,comm_mat,B_mat,solu)
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Solved. Solution ", solu
  END IF
  solu_mat(1:ni,1) = solu(1:ni)

  B_check = MATMUL(comm_mat,solu_mat)
  IF (DEBUG .EQV. .TRUE.) THEN
    check_1 = SUM(ABS(B_check(1:ni,1) - B_mat(1:ni)))
    WRITE(*,*) "CHECK: ", B_check(1:ni,1) - B_mat(1:ni), SUM(ABS(B_check(1:ni,1) - B_mat(1:ni)))
  END IF

  ! Set up the commutator Vector putting I then A then B entries in place by indexes given
  DO i = 1, ni
    comm_f(i_def(i)) = solu(i)
  END DO
  DO i = 1, na
    comm_f(A_def(i)) = 0.0d0
  END DO
  DO i = 1, nb
    comm_f(B_def(i)) = 1.0d0
  END DO

  ! Backwards committor is 1 - forwards committor for detailed balance systems; but those aren't all I might care about

  ! backwards ~T_{i,j) = pi(j)/pi(i)*T_{j,i}
  DO i = 1, ni
    DO j = 1, ni
      ! Indexes of the I matrix extracted form the array of its definition
      IF (i .EQ. j) THEN
        comm_mat(i,j) = msm(i_def(j),i_def(i))*pis(i_def(j))/pis(i_def(i)) - 1.0d0
      ELSE
        comm_mat(i,j) = msm(i_def(j),i_def(i))*pis(i_def(j))/pis(i_def(i))
      END IF
    END DO
  END DO
  ! Assemble the B of the A*X = B system; B_i = summation over -T_{i,j} where j is in B
  B_mat = 0.0d0
  DO i = 1, ni
    DO j = 1, na
      B_mat(i) = B_mat(i) - msm(a_def(j),i_def(i))*pis(a_def(j))/pis(i_def(i))
    END DO
  END DO 
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Assembled B_mat", B_mat
    WRITE(*,*) "Comm_mat "
    CALL print_mat_r8(comm_mat,6)
    WRITE(*,*) "Sending to solve"
  END IF

  ! Solve the system to get the backwards committors
  CALL solve_system(ni,comm_mat,B_mat,solu)
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Solved. Solution ", solu
  END IF
  solu_mat(1:ni,1) = solu(1:ni)

  !B_check = MATMUL(comm_mat,solu_mat)
  IF (DEBUG .EQV. .TRUE.) THEN
    check_2 = SUM(ABS(B_check(1:ni,1) - B_mat(1:ni)))
    WRITE(*,*) "CHECK: ", B_check(1:ni,1) - B_mat(1:ni), SUM(ABS(B_check(1:ni,1) - B_mat(1:ni)))
  END IF

  ! Set up the commutator Vector putting I then A then B entries in place by indexes given
  DO i = 1, ni
    comm_b(i_def(i)) = solu(i)
  END DO
  DO i = 1, na
    comm_b(A_def(i)) = 1.0d0
  END DO
  DO i = 1, nb
    comm_b(B_def(i)) = 0.0d0
  END DO

  IF (detailed_balance .EQV. .TRUE.) THEN  ! Overwrite with detailed balance information
    IF (check_2 .LT. check_1) THEN
      comm_f = 1.0d0 - comm_b
    ELSE
      comm_b = 1.0d0 - comm_f  ! Detailed balance solution
    END IF
  END IF


  ! No committor may be larger than 1.0 or smaller than 0.0. Brute force fix for rounding errors.
  DO i = 1, d
    IF (comm_f(i) .GT. 1.0d0) THEN
      comm_f(i) = 1.0d0
    END IF
    IF (comm_f(i) .LT. 0.0d0) THEN
      comm_f(i) = 0.0d0
    END IF
    IF (comm_b(i) .GT. 1.0d0) THEN
      comm_b(i) = 1.0d0
    END IF
    IF (comm_b(i) .LT. 0.0d0) THEN
      comm_b(i) = 0.0d0
    END IF
  END DO
  
END SUBROUTINE

! Go through a MSM and change all fluxes into net A to B fluxes
! MSM(i,j) = MAX[0,MSM(i,j)-MSM(j,i)]
SUBROUTINE msm_to_net_flux(d,msm)
  INTEGER, INTENT(IN) :: d
  REAL*8, DIMENSION(d,d), INTENT(INOUT) :: msm
  INTEGER :: i, j
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "msm to net flux start flux "
    CALL print_mat_r8(msm)
  END IF
  DO i = 1, d
    DO j = i + 1, d
      IF (msm(i,j) .GT. msm(j,i)) THEN
        msm(i,j) = msm(i,j) - msm(j,i)
        msm(j,i) = 0.0d0
      ELSE
        msm(j,i) = msm(j,i) - msm(i,j)
        msm(i,j) = 0.0d0
      END IF
    END DO
  END DO
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "msm to net flux end flux "
    CALL print_mat_r8(msm)
  END IF
END SUBROUTINE

! Given a MSM, calculate for each f_ij entry
! \Pi_i q_i^- Tij q_j^+ (or 0 if i = j) meaning the net
! flux heading to the products along that line during dt
SUBROUTINE calc_fij_AB_flux(d,msm,flux_msm,comm_f,comm_b,pis,factor)
  INTEGER, INTENT(IN) :: d
  REAL*8, DIMENSION(d,d), INTENT(IN) :: msm  ! The Markov state model
  REAL*8, DIMENSION(d,d), INTENT(OUT) :: flux_msm  ! Flux matrix, not net flux
  REAL*8, DIMENSION(d), INTENT(IN) :: comm_f, comm_b, pis  ! Forward committor, backward committor, equilib. prob to be in state
  REAL*8, INTENT(IN) :: factor  ! Constant multiple to apply to all fluxes; use 1 if not needed; deals with machine precision
  INTEGER :: i, j

  flux_msm = 0.0d0
  DO i = 1, d
    DO j = 1, d
      IF (i .NE. j) THEN
        flux_msm(i,j) = pis(i)*factor
        IF (DEBUG .EQV. .TRUE.) THEN
          WRITE(*,*) "pis", pis(i), "factor", factor
        END IF
        flux_msm(i,j) = flux_msm(i,j)*comm_b(i)*msm(i,j)*comm_f(j)  
        IF (DEBUG .EQV. .TRUE.) THEN
          WRITE(*,*) "comm_b", comm_b(i), "comm_f", comm_f(j), "msm(i,j)", msm(i,j), "flux", i, j, "=", flux_msm(i,j)
        END IF
        IF (flux_msm(i,j) .LT. 0) THEN 
          WRITE(*,*) "WARNING: NEGATIVE FLUX i j", i, j, flux_msm(i,j)  ! Something is very wrong
        END IF
      END IF
    END DO
  END DO 
END SUBROUTINE

! Calculate the reactive flux by both \sum a in A, j not in A f(a,j)
! and \sum j not in B, b in B f(j,B)
REAL*8 FUNCTION reactive_flux(d,flux_msm,ni,na,nb,i_def,a_def,b_def)
  INTEGER, INTENT(IN) :: d, ni, nb, na  ! Dimension of n, Numbers of eigenstate in I, A, B
  REAL*8, DIMENSION(d,d), INTENT(IN) :: flux_msm  ! The flux matrix flux_msm(i,j) = f(i-->j)
  INTEGER, DIMENSION(ni), INTENT(IN) :: i_def
  INTEGER, DIMENSION(na), INTENT(IN) :: a_def
  INTEGER, DIMENSION(nb), INTENT(IN) :: b_def
  INTEGER :: i, j
  REAL*8 :: check
  reactive_flux = 0.0d0
  check = 0.0d0
  DO i = 1, ni
    DO j = 1, nb
      reactive_flux = reactive_flux + flux_msm(i_def(i),b_def(j))
      IF (DEBUG .EQV. .TRUE.) THEN
        WRITE(*,*) "Adding reactive flux from ", i_def(i), "to", b_def(j), "value", flux_msm(i_def(i),b_def(j))
      END IF
    END DO
    DO j = 1, na
      check = check + flux_msm(a_def(j),i_def(i))
      IF (DEBUG .EQV. .TRUE.) THEN
        WRITE(*,*) "Adding reactive flux from ", a_def(j), "to", i_def(i), "value", flux_msm(a_def(j),i_def(i))
      END IF
    END DO
  END DO
END FUNCTION

! Continuously extract paths through MSM vertices; pass them back in lengths, paths
! Extracts paths in order of weight (max min flux path/bottleneck selection).
! In the case that there is more than one source or destination node, this gets
! complicated and this routine can't help you there.
SUBROUTINE extract_all_paths(d,d2,flux,sink,source,A_width,A_previous,npaths,paths,lpaths,wpaths,do_scaling)
  INTEGER, INTENT(INOUT), DIMENSION(d) :: A_previous  ! Backtrack information for graph
  REAL*8, INTENT(INOUT), DIMENSION(d,d) :: flux ! Original net flux matrix; will be modified heavily
  INTEGER, INTENT(IN) :: d, d2, source, sink  ! MSM dimension, max allowed number of paths, Origin and destination nodes in graph
  INTEGER, INTENT(OUT) :: npaths  ! Number of paths found
  REAL*8, INTENT(INOUT), DIMENSION(d) :: A_width  ! Widths of original graph, will be modified
  INTEGER, INTENT(OUT), DIMENSION(d2) :: lpaths  ! Length of all the paths in paths
  REAL*8, INTENT(OUT), DIMENSION(d2) :: wpaths  ! The weight of all the paths in paths
  INTEGER, INTENT(OUT), DIMENSION(d,d2) :: paths  ! Each column is a path through the graph
  LOGICAL, INTENT(IN) :: do_scaling  ! Whether to apply the normalizer; don't do this for vertical paths
  INTEGER :: bt_len, i, temp
  REAL*8 :: normalizer

  IF (do_scaling .EQV. .TRUE.) THEN 
    normalizer = MAX(flux(sink,source),flux(source,sink))
    IF (normalizer .LE. 0) THEN  ! Don't let that be a thing; fix it
      normalizer = 1.0d-10
    END IF
    flux = flux/normalizer  ! Do away with errors due to tiny numbers
    IF (DEBUG .EQV. .TRUE.) THEN
      WRITE(*,*) "Normalizing: sink ", sink, "source", source, "d is", d, "d2 is", d2, "normalizer", normalizer
    END IF
  END IF


  paths = sink  ! Convenient choice for plotting 
  lpaths = 0.0d0
  wpaths = 0.0d0
  npaths = 0
  bt_len = 1  ! Arbitrary to get the loop started
  CALL dijkstra_wp(d,flux,A_width,A_previous,source)  ! Rerun Dijkstra to find the next widest path
  DO WHILE (bt_len .GT. 0 .AND. npaths .LT. d2)  ! Until a path subtraction failure or we overflow the path matrix 
    IF (DEBUG .EQV. .TRUE.) THEN
      WRITE(*,*) "Current A_previous", A_previous
    END IF
    npaths = npaths + 1 ! Cool, path has been found
    bt_len = back_track_sub(d,flux,sink,source,A_width,A_previous,paths(1:d,npaths))  ! Follow path backwards and subtract
    IF (DEBUG .EQV. .TRUE.) THEN
      WRITE(*,*) "Back tracked. bt_len ", bt_len
    END IF
    lpaths(npaths) = bt_len
    wpaths(npaths) = A_width(sink)  ! Current width is width to the sink node

    ! Reverse the order of the path
    DO i = 1, FLOOR(bt_len/2.0d0)
      temp = paths(i,npaths)
      paths(i,npaths) = paths(bt_len - i + 1,npaths)
      paths(bt_len - i + 1,npaths) = temp 
    END DO

    IF (DEBUG .EQV. .TRUE.) THEN
      WRITE(*,*) "flux matrix "!, flux
      CALL print_mat_r8(flux)
    END IF

    CALL dijkstra_wp(d,flux,A_width,A_previous,source)  ! Rerun Dijkstra to find the next widest path
    IF (A_width(sink) .LT. 0) THEN  ! Well, that's not going to work
      EXIT
    END IF
  END DO

  IF (do_scaling .EQV. .TRUE.) THEN
    flux = flux*normalizer
    wpaths = wpaths*normalizer  ! Renormalize to the proper width after multiplying for the course of extraction
  END IF
END SUBROUTINE

! Subroutine that, given A, B, I, source and sink will calculate committors, fluxes, paths
! through the graph etc. Requires a lot of information.
! This is the routine for thermal data and an equilibrium MSM. Analyzing a graph representing 
! nonequilibrium dynamics is a different routine
! The 6*d is a totally arbitrarily chosen length allowing for 6*d paths in total to be
! analyzed by the routine (there has to be some upper limit). This has seemed sufficient for
! all work done. If the reactive flux found by the path backtracking doesn't match that found
! by summation of flux leaving the reactants, that 6 may need to be increased to find additional paths
SUBROUTINE perform_analysis(d,ni,na,nb,A_def,B_def,I_def,dt,evals,factor,source,sink,comm_file,&
                            flux_file,p_file1,p_file2,mmat,pis, comm_f_out,comm_b_out, &
                            evec_pc,evec_pt,evec_p1,evec_p2)
  INTEGER, DIMENSION(d), INTENT(IN) :: A_def, B_def, I_def  ! Definitions for wells and intermediate states
  INTEGER, INTENT(IN) :: ni, na, nb  ! Sizes of intermediate, reactant and product definitions
  INTEGER, INTENT(IN) :: d, source, sink  ! Dimension, beginning and end of paths to trace
  REAL*8, INTENT(IN) :: factor, dt  ! Factor by which to multiply reactive fluxes during calculation MSM time
  REAL*8, DIMENSION(d,d), INTENT(IN) :: mmat ! Markov model
  REAL*8, DIMENSION(d), INTENT(IN) :: pis, evals  ! Equilibrium probabilities
  INTEGER, INTENT(IN) :: comm_file, flux_file, p_file1, p_file2  ! Write files
  REAL*8, DIMENSION(d), INTENT(OUT), OPTIONAL :: comm_f_out, comm_b_out ! Optional place to pass back committors to calling program
  REAL*8, DIMENSION(d), INTENT(IN) :: evec_pc, evec_pt, evec_p1, evec_p2  ! Average positions, characters of evecs

  REAL*8, DIMENSION(d) :: comm_f, comm_b  ! Forward and backwards committors
  INTEGER, DIMENSION(d) :: A_previous  ! Backtracking information
  REAL*8, DIMENSION(d) :: A_width  ! Width information 
  REAL*8, DIMENSION(d,d) :: flux, net_flux, net_flux_temp  ! Reactive flux calculations
  INTEGER, DIMENSION(6*d) :: lpath  ! Path backtracking
  REAL*8, DIMENSION(6*d) :: wpath  ! Path width information
  REAL*8, DIMENSION(d,6*d) :: eval_paths  ! Path through energy space
  INTEGER, DIMENSION(d,6*d) :: paths  ! Paths (through indices)
  INTEGER :: i, j, npaths, bt_len
  REAL*8 :: rflux, temp, pi_A

  REAL*8, DIMENSION(d+1,d+1) :: comm_jumps  ! comm_jumps(i,j) contains number of times i-->j was the committor hop
  REAL*8 :: avg_comm_en, avg_comm_pc, avg_comm_pt, avg_comm_p1, avg_comm_p2, avg_comm  ! Average values of the committor states
  REAL*8 :: avg_pcomm_en, avg_pcomm_pc, avg_pcomm_pt, avg_pcomm_p1, avg_pcomm_p2, avg_pcomm ! Average values of the pre-committor states
  REAL*8 :: avg_diff_en, avg_diff_pc, avg_diff_pt, avg_diff_p1, avg_diff_p2, avg_diff  ! Average change at committor jump
  REAL*8 :: total_jumps

  comm_f = 0.0d0
  comm_b = 0.0d0
  CALL calc_comm(d,mmat,comm_f,comm_b,ni,nb,na,I_def,A_def,B_def,pis)
  WRITE(comm_file,*) 0, 0.0, 0.0
  DO i = 1, d
    WRITE(comm_file,*) i, comm_f(i), comm_b(i)
  END DO

  IF (PRESENT(comm_f_out)) THEN
    comm_f_out = comm_f
  END IF
  IF (PRESENT(comm_b_out)) THEN
    comm_b_out = comm_b
  END IF

  ! Total probability to be on a forward move
  pi_A = 0.0d0
  DO i = 1, d
    pi_A = pi_A + pis(i)*comm_b(i)
  END DO

  flux = 0.0d0
  net_flux = 0.0d0
  net_flux_temp = 0.0d0
  CALL calc_fij_AB_flux(d,mmat,flux,comm_f,comm_b,pis,factor)
  rflux = reactive_flux(d,flux,ni,na,nb,I_def,A_def,B_def)
  net_flux = flux
  CALL msm_to_net_flux(d,net_flux) 
  rflux = reactive_flux(d,net_flux,ni,na,nb,I_def,A_def,B_def)
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Calculating rate rflux", rflux, "factor", factor, "pi_A", pi_A, "dt", dt
  END IF
  WRITE(flux_file,*) dt, "NET FLUX ", rflux, "RATE", rflux/(factor*pi_A*dt)
  net_flux_temp = net_flux

  CALL extract_all_paths(d,6*d,net_flux_temp,sink,source,A_width,A_previous,npaths,paths,lpath,wpath,.TRUE.)
  temp = SUM(wpath(1:npaths)) ! In most cases we won't be truncating the number of paths
  WRITE(p_file1,*) "Extracted paths with total flux", temp
  DO i = 1, npaths
    WRITE(p_file1,*) wpath(i), wpath(i)/temp, paths(1:lpath(i),i)
  END DO
  WRITE(p_file1,*) ""
  WRITE(p_file1,*) ""
 
  bt_len = MAXVAL(lpath(1:npaths),1)
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "bt_len is", bt_len
  END IF
  DO i = 1, bt_len
    WRITE(p_file2,*) i, paths(i,1:d), wpath(i), wpath(i)/temp
  END DO  
  WRITE(p_file2,*) ""
  WRITE(p_file2,*) ""

  DO j = 1, d
    ! Translate paths into energy space
    DO i = 1, d
      eval_paths(j,i) = evals(paths(j,i))
    END DO
  END DO

  DO i = 1, npaths
    WRITE(p_file1,*) wpath(i), wpath(i)/temp, eval_paths(1:lpath(i),i)
    WRITE(p_file2,*) i, eval_paths(i,1:d), wpath(i), wpath(i)/temp
  END DO

  WRITE(p_file1,*) ""
  WRITE(p_file1,*) ""
  ! Use committors to find out at which jump the committor passes 0.5
  comm_jumps = 0
  total_jumps = 0
  DO i = 1, npaths
    DO j = 2, lpath(i) ! First entry guaranteed not to have over 50% comm 
      IF (comm_f(paths(j,i)) .GT. 0.5) THEN  ! Committor jump has been made; increment index
        comm_jumps(paths(j-1,i), paths(j,i)) = comm_jumps(paths(j-1,i), paths(j,i)) + wpath(i)
        EXIT  ! Break from the loop and see to the next path
      END IF
    END DO
    ! Increment the total number of jumps
    total_jumps = total_jumps + wpath(i)
  END DO

  avg_pcomm = 0.d0
  avg_pcomm_en = 0.0d0
  avg_pcomm_pc = 0.0d0
  avg_pcomm_pt = 0.0d0
  avg_pcomm_p1 = 0.0d0
  avg_pcomm_p2 = 0.0d0
  avg_diff = 0.0d0
  avg_diff_en = 0.0d0
  avg_diff_pc = 0.0d0
  avg_diff_pt = 0.0d0
  avg_diff_p1 = 0.0d0
  avg_diff_p2 = 0.0d0
  avg_comm = 0.0d0
  avg_comm_en = 0.0d0
  avg_comm_pc = 0.0d0
  avg_comm_pt = 0.0d0
  avg_comm_p1 = 0.0d0
  avg_comm_p2 = 0.0d0
  ! Print out where the committor jump happens
  WRITE(p_file1,*) "Number of paths total", total_jumps
  DO i = 1, d+1
    DO j = 1, d+1
      IF (comm_jumps(i,j) .GT. 0) THEN
        WRITE(p_file1,*) i, "----->", j, "is", comm_jumps(i,j), "percent", 1.0d0*comm_jumps(i,j)/total_jumps, &
                       "from eval", evals(i), "to eval", evals(j), "difference", evals(j) - evals(i)
        ! Running average information about the post committor state
        avg_comm = avg_comm + comm_f(j)*comm_jumps(i,j)/total_jumps
        avg_comm_en = avg_comm_en + evals(j)*comm_jumps(i,j)/total_jumps
        avg_comm_pc = avg_comm_pc + evec_pc(j)*comm_jumps(i,j)/total_jumps
        avg_comm_pt = avg_comm_pt + evec_pt(j)*comm_jumps(i,j)/total_jumps
        avg_comm_p1 = avg_comm_p1 + evec_p1(j)*comm_jumps(i,j)/total_jumps
        avg_comm_p2 = avg_comm_p2 + evec_p2(j)*comm_jumps(i,j)/total_jumps
        ! Running average information about the pre committor state
        avg_pcomm = avg_pcomm + comm_f(i)*comm_jumps(i,j)/total_jumps
        avg_pcomm_en = avg_pcomm_en + evals(i)*comm_jumps(i,j)/total_jumps
        avg_pcomm_pc = avg_pcomm_pc + evec_pc(i)*comm_jumps(i,j)/total_jumps
        avg_pcomm_pt = avg_pcomm_pt + evec_pt(i)*comm_jumps(i,j)/total_jumps
        avg_pcomm_p1 = avg_pcomm_p1 + evec_p1(i)*comm_jumps(i,j)/total_jumps
        avg_pcomm_p2 = avg_pcomm_p2 + evec_p2(i)*comm_jumps(i,j)/total_jumps
        ! Running average information about the difference between pre/post
        avg_diff = avg_diff + (comm_f(j) - comm_f(i))*comm_jumps(i,j)/total_jumps
        avg_diff_en = avg_diff_en + (evals(j) - evals(i))*comm_jumps(i,j)/total_jumps
        avg_diff_pc = avg_diff_pc + (evec_pc(j) - evec_pc(i))*comm_jumps(i,j)/total_jumps
        avg_diff_pt = avg_diff_pt + (evec_pt(j) - evec_pt(i))*comm_jumps(i,j)/total_jumps
        avg_diff_p1 = avg_diff_p1 + (evec_p1(j) - evec_p1(i))*comm_jumps(i,j)/total_jumps
        avg_diff_p2 = avg_diff_p2 + (evec_p2(j) - evec_p2(i))*comm_jumps(i,j)/total_jumps
      END IF
    END DO
  END DO
  WRITE(p_file1,*) ""
  WRITE(p_file1,*) ""
  WRITE(p_file1,*) "AVERAGES FOR THE PRECOMMITTOR JUMP (ENERGY, QC, QT, P1, P2, COMM_F)"
  WRITE(p_file1,*) avg_pcomm_en, avg_pcomm_pc, avg_pcomm_pt, avg_pcomm_p1, avg_pcomm_p2, avg_pcomm
  WRITE(p_file1,*) "AVERAGES FOR THE POSTCOMMITTOR JUMP (ENERGY, QC, QT, P1, P2, COMM_F)"
  WRITE(p_file1,*) avg_comm_en, avg_comm_pc, avg_comm_pt, avg_comm_p1, avg_comm_p2, avg_comm
  WRITE(p_file1,*) "AVERAGE DIFFERENCE THROUGH THE JUMP (ENERGY, QC, QT, P1, P2, COMM_F)"
  WRITE(p_file1,*) avg_diff_en, avg_diff_pc, avg_diff_pt, avg_diff_p1, avg_diff_p2, avg_diff
  WRITE(p_file1,*) ""
  WRITE(p_file1,*) ""

  WRITE(p_file1,*) "RUNNING SUM OF FLUX VS PATHWAYS INCLUDED"
  WRITE(p_file1,*) ""
  WRITE(p_file1,*) ""

  CALL running_sum(npaths,wpath)  ! Making a running summation out of wpath; this is redundant since it's already sorted but whatever
  DO i = 1, npaths
    WRITE(p_file1,*) i, wpath(i), wpath(i)/total_jumps
  END DO
  WRITE(p_file1,*) ""
  WRITE(p_file1,*) ""

END SUBROUTINE


! Subroutine for analyzing a vertical deexcitation graph; goes through most of the same steps
! with most of the same components as the thermal version above.
SUBROUTINE graph_analysis(d,source,sink,file1,file2,graph,evals,comm,evec_pc,evec_pt,evec_p1,evec_p2)
  INTEGER, INTENT(IN) :: d, source, sink  ! Dimension, beginning and end of paths to trace
  INTEGER, DIMENSION(d,d), INTENT(IN) :: graph ! Graph for path tracing
  INTEGER, INTENT(IN) :: file1, file2 ! Write files
  REAL*8, DIMENSION(d), INTENT(IN) :: evals

  REAL*8, DIMENSION(d), INTENT(IN) :: evec_pc, evec_pt, evec_p1, evec_p2  ! Average positions, characters of evecs
  REAL*8, DIMENSION(d), INTENT(IN) :: comm

  INTEGER, DIMENSION(d) :: A_previous  ! Dijkstra working arrays and backtracking information
  REAL*8, DIMENSION(d) :: A_width  ! Dijkstra working array
  REAL*8, DIMENSION(d,d) :: flux  ! Reactive flux calculations
  INTEGER, DIMENSION(6*d) :: lpath  ! Path backtracking length
  REAL*8, DIMENSION(6*d) :: wpath  ! Path weight
  INTEGER, DIMENSION(d,6*d) :: paths  ! Path jump record
  REAL*8, DIMENSION(d,6*d) :: energy_paths  ! Path jump record in terms of energy
  INTEGER :: i, j, npaths, bt_len
  REAL*8 :: temp

  REAL*8, DIMENSION(d+1,d+1) :: comm_jumps  ! comm_jumps(i,j) contains number of times i-->j was the committor hop
  REAL*8 :: avg_comm_en, avg_comm_pc, avg_comm_pt, avg_comm_p1, avg_comm_p2, avg_comm  ! Average values of the committor states
  REAL*8 :: avg_pcomm_en, avg_pcomm_pc, avg_pcomm_pt, avg_pcomm_p1, avg_pcomm_p2, avg_pcomm ! Average values of the pre-committor states
  REAL*8 :: avg_diff_en, avg_diff_pc, avg_diff_pt, avg_diff_p1, avg_diff_p2, avg_diff  ! Average change at committor jump
  INTEGER :: total_jumps
  INTEGER :: cj_from_0  ! Committor jumps which occur straight from zero to x state; no energy knowledge for them

  ! Convert raw graph to net jumps between eigenstates
  flux = 1.0d0*graph
  DO i = 1, d
    DO j = 1, d
      IF (flux(i,j) .GT. flux(j,i)) THEN
        flux(i,j) = flux(i,j) - flux(j,i)
        flux(j,i) = 0.0d0
      ELSE IF (flux(i,j) .LT. flux(j,i)) THEN
        flux(j,i) = flux(j,i) - flux(i,j)
        flux(i,j) = 0.0d0
      ELSE  ! (i,i) entries are 0
        flux(i,j) = 0.0d0
        flux(j,i) = 0.0d0
      END IF
    END DO
  END DO

  OPEN(89,FILE="net_graph.txt")
  WRITE(89,*) "Net jump flux"
  CALL print_mat_r8(flux,89)
  CLOSE(89)

  ! Get all the paths in order of width from Dijkstra's algorithm
  CALL extract_all_paths(d,6*d,flux,sink,source,A_width,A_previous,npaths,paths,lpath,wpath,.FALSE.)

  temp = SUM(wpath(1:npaths))
  DO i = 1, npaths
    WRITE(file1,*) wpath(i), wpath(i)/temp, paths(1:lpath(i),i)
  END DO
  WRITE(file1,*) ""
  WRITE(file1,*) ""
 
  bt_len = MAXVAL(lpath(1:npaths),1)
  DO i = 1, bt_len
    WRITE(file2,*) i, paths(i,1:npaths), wpath(i), wpath(i)/temp
  END DO  
  WRITE(file2,*) ""
  WRITE(file2,*) ""

  ! Translate jumps into energy space
  DO i = 1, d
    DO j = 1, 6*d
      IF (paths(i,j) .GT. d) THEN
        energy_paths(i,j) = evals(d)*2.0d0  ! Arbitrary number; we don't know what energy prior to collapse was
      ELSE
        energy_paths(i,j) = evals(paths(i,j))
      END IF
    END DO
  END DO

  ! Write out the energies of the jump paths rather than eigenstates
  temp = SUM(wpath(1:npaths))
  DO i = 1, npaths
    WRITE(file1,*) wpath(i), wpath(i)/temp, energy_paths(1:lpath(i),i)
  END DO
  WRITE(file1,*) ""
  WRITE(file1,*) ""
 
  bt_len = MAXVAL(lpath(1:npaths),1)
  DO i = 1, bt_len
    WRITE(file2,*) i, energy_paths(i,1:npaths), wpath(i), wpath(i)/temp
  END DO  
  WRITE(file2,*) ""
  WRITE(file2,*) ""

  ! Use committors to find out at which jump the committor passes 0.5
  comm_jumps = 0
  total_jumps = 0
  DO i = 1, npaths
    DO j = 2, lpath(i)  ! Skip over the vertical entry; it does not have a committor
      IF (comm(paths(j,i)) .GT. 0.5) THEN  ! Committor jump has been made; increment index
        comm_jumps(paths(j-1,i), paths(j,i)) = comm_jumps(paths(j-1,i), paths(j,i)) + wpath(i)
        EXIT  ! Break from the loop and see to the next path
      END IF
    END DO
    ! Increment the total number of jumps
    total_jumps = total_jumps + wpath(i)
  END DO

  avg_pcomm = 0.d0
  avg_pcomm_en = 0.0d0
  avg_pcomm_pc = 0.0d0
  avg_pcomm_pt = 0.0d0
  avg_pcomm_p1 = 0.0d0
  avg_pcomm_p2 = 0.0d0
  avg_diff = 0.0d0
  avg_diff_en = 0.0d0
  avg_diff_pc = 0.0d0
  avg_diff_pt = 0.0d0
  avg_diff_p1 = 0.0d0
  avg_diff_p2 = 0.0d0
  avg_comm = 0.0d0
  avg_comm_en = 0.0d0
  avg_comm_pc = 0.0d0
  avg_comm_pt = 0.0d0
  avg_comm_p1 = 0.0d0
  avg_comm_p2 = 0.0d0
  cj_from_0 = 0
  ! Print out where the committor jump happens
  WRITE(file1,*) "Number of paths total", total_jumps
  DO i = 1, d+1
    DO j = 1, d+1
      IF (comm_jumps(i,j) .GT. 0) THEN
        WRITE(file1,*) i, "----->", j, "is", comm_jumps(i,j), "percent", 1.0d0*comm_jumps(i,j)/total_jumps, &
                       "from eval", evals(i), "to eval", evals(j), "difference", evals(j) - evals(i)
        ! Running average information about the post committor state
        avg_comm = avg_comm + comm(j)*comm_jumps(i,j)/total_jumps
        avg_comm_en = avg_comm_en + evals(j)*comm_jumps(i,j)/total_jumps
        avg_comm_pc = avg_comm_pc + evec_pc(j)*comm_jumps(i,j)/total_jumps
        avg_comm_pt = avg_comm_pt + evec_pt(j)*comm_jumps(i,j)/total_jumps
        avg_comm_p1 = avg_comm_p1 + evec_p1(j)*comm_jumps(i,j)/total_jumps
        avg_comm_p2 = avg_comm_p2 + evec_p2(j)*comm_jumps(i,j)/total_jumps
        IF (i .EQ. source) THEN  ! Jump from uncollapsed state about which we know nothing
          cj_from_0 = cj_from_0 + comm_jumps(i,j)
        ELSE
          ! Running average information about the pre committor state
          avg_pcomm = avg_pcomm + comm(i)*comm_jumps(i,j)/total_jumps
          avg_pcomm_en = avg_pcomm_en + evals(i)*comm_jumps(i,j)/total_jumps
          avg_pcomm_pc = avg_pcomm_pc + evec_pc(i)*comm_jumps(i,j)/total_jumps
          avg_pcomm_pt = avg_pcomm_pt + evec_pt(i)*comm_jumps(i,j)/total_jumps
          avg_pcomm_p1 = avg_pcomm_p1 + evec_p1(i)*comm_jumps(i,j)/total_jumps
          avg_pcomm_p2 = avg_pcomm_p2 + evec_p2(i)*comm_jumps(i,j)/total_jumps
          ! Running average information about the difference between pre/post
          avg_diff = avg_diff + (comm(j) - comm(i))*comm_jumps(i,j)/total_jumps
          avg_diff_en = avg_diff_en + (evals(j) - evals(i))*comm_jumps(i,j)/total_jumps
          avg_diff_pc = avg_diff_pc + (evec_pc(j) - evec_pc(i))*comm_jumps(i,j)/total_jumps
          avg_diff_pt = avg_diff_pt + (evec_pt(j) - evec_pt(i))*comm_jumps(i,j)/total_jumps
          avg_diff_p1 = avg_diff_p1 + (evec_p1(j) - evec_p1(i))*comm_jumps(i,j)/total_jumps
          avg_diff_p2 = avg_diff_p2 + (evec_p2(j) - evec_p2(i))*comm_jumps(i,j)/total_jumps
        END IF
      END IF
    END DO
  END DO
  ! Renormalize the averages for the position and energy for the pre and difference jumps due to 
  ! the presence of committor jumps from 0, the uncollapsed state about which nothing is known
  avg_pcomm = avg_pcomm*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_pcomm_en = avg_pcomm_en*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_pcomm_pc = avg_pcomm_pc*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_pcomm_pt = avg_pcomm_pt*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_pcomm_p1 = avg_pcomm_p1*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_pcomm_p2 = avg_pcomm_p2*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_diff = avg_diff*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_diff_en = avg_diff_en*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_diff_pc = avg_diff_pc*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_diff_pt = avg_diff_pt*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_diff_p1 = avg_diff_p1*1.0d0*total_jumps/(total_jumps - cj_from_0)
  avg_diff_p2 = avg_diff_p2*1.0d0*total_jumps/(total_jumps - cj_from_0)

  WRITE(file1,*) ""
  WRITE(file1,*) ""
  WRITE(file1,*) "TOTAL COMMITTOR JUMPS MADE FROM 0 ", cj_from_0!, "RENORM VALUE",&
  !                1.0d0*total_jumps/(total_jumps - cj_from_0)
  ! Once upon a time, I got some information from looking at these commented out numbers, but they
  ! don't really mean anything.
  !WRITE(file1,*) "AVERAGES FOR THE PRECOMMITTOR JUMP (ENERGY, QC, QT, P1, P2, COMM_F)"
  !WRITE(file1,*) avg_pcomm_en, avg_pcomm_pc, avg_pcomm_pt, avg_pcomm_p1, avg_pcomm_p2, avg_pcomm
  WRITE(file1,*) "AVERAGES FOR THE POSTCOMMITTOR JUMP (ENERGY, QC, QT, P1, P2, COMM_F)"
  WRITE(file1,*) avg_comm_en, avg_comm_pc, avg_comm_pt, avg_comm_p1, avg_comm_p2, avg_comm
  !WRITE(file1,*) "AVERAGE DIFFERENCE THROUGH THE JUMP (ENERGY, QC, QT, P1, P2, COMM_F)"
  !WRITE(file1,*) avg_diff_en, avg_diff_pc, avg_diff_pt, avg_diff_p1, avg_diff_p2, avg_diff
  WRITE(file1,*) ""
  WRITE(file1,*) ""



  CALL running_sum(npaths,wpath)  ! Making a running summation out of wpath; this is redundant I think since it's already sorted but whatever
  DO i = 1, npaths
    WRITE(file1,*) i, wpath(i), wpath(i)/total_jumps
  END DO
  WRITE(file1,*) ""
  WRITE(file1,*) ""
END SUBROUTINE



END MODULE



