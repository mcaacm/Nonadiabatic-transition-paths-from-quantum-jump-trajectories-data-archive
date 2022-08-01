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

! Contains functions for the bulk of the machinery of Lindblad
! jump dynamics

MODULE lindblad_coherent
USE parameters
USE prequel
USE rk_ns_utils
IMPLICIT NONE


! The implementation of Lindblad jumps here is largely based on 
! "Vogt, N., Jeske, J. and Cole, J. H. "Stochastic Bloch-Redfield theory:
! Quantum jumps in a solid-state environment" Phys. Rev. B 88, 2013.

CONTAINS

! Returns how many entries in a wavefunction are above
! 1e-15, which is a good indication of nonzero numbers
INTEGER FUNCTION non_zero_num(wfn)
  COMPLEX*16, INTENT(IN), DIMENSION(n,1) :: wfn 
  INTEGER i
  non_zero_num = 0
  DO i = 1, n
    IF (ABS(wfn(i,1)) .GE. 1.0d-15) THEN
      non_zero_num = non_zero_num + 1
    END IF
  END DO
END FUNCTIOn

! Organize operators by index so that jump_pointers (i,j) contains the index of
! the jth operator that will act on eigenstate i. Then, to find only the relevant
! jump operators for eigenstate i, just go through the indexes rather than 
! searching everything. This information is stored globally. It greatly
! cuts down on the amount of time spent checking which operator is relevant.
SUBROUTINE organize_operators(d, spaces)
  INTEGER, INTENT(IN) :: d  ! Number of operators
  TYPE(sspace), DIMENSION(d), INTENT(IN) :: spaces  ! Jump operators
  INTEGER, DIMENSION(n) :: jp_indexes   
  INTEGER :: nent, i, j, ind

  jp_indexes = 1
  ALLOCATE(jump_pointers(n,jump_pointer_dim))
  jump_pointers = -1

  ! For all subspaces
  DO i = 1, d
    nent = spaces(i)%Aomega%nentries
    ! For all the indexes it operators on
    DO j = 1, nent
      ind = spaces(i)%Aomega%cols(j)
      jump_pointers(ind,jp_indexes(ind)) = i   ! Input index of relevant operator
      IF (DEBUG .EQV. .TRUE.) THEN
        WRITE(*,*) "jp_pointers ", ind, "now", jump_pointers(ind,jp_indexes(ind))
      END IF
      jp_indexes(ind) = jp_indexes(ind) + 1
      IF (jp_indexes(ind) .GT. jump_pointer_dim) THEN
        WRITE(*,*) "ERROR: jp_indexes exceeded size of jump_pointer_dim. Raise jump_pointer_dim in prequel.f90"
        STOP
      END IF
    END DO
  END DO
  IF (DEBUG .EQV. .TRUE.) THEN
    WRITE(*,*) "Sorted, here is the result indices", jp_indexes
    WRITE(*,*) jump_pointers
  END IF
END SUBROUTINE


! Prepare jump operators assuming that every delta is different,
! in other words that there are n**2 - n + 1 total jump operators (with
! that 1 being the identity). This assumes that no eigenvalues are
! the same within tolerance and all subspaces have only one entry.
! In other words, the total secular approximation with fully independent
! subspaces.
SUBROUTINE new_sort(evals, spaces, t_len)
  REAL*8, INTENT(IN), DIMENSION(n) :: evals  ! Energy eigenvalues of lambshift H
  INTEGER, INTENT(OUT) :: t_len  ! Number of subspaces
  TYPE(sspace), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: spaces  ! Jump operators for system
  INTEGER :: i, j, k, l
  t_len = n**2 - n + 1

  ! Check for degeneracies 
  DO i = 2, n
    IF (ABS(evals(i-1) - evals(i)) .LT. 1.0d-10) THEN
      ! Problem--eigenvalue degeneracy
      WRITE(*,*) "Eigenvalue degeneracy: ", evals(i-1), evals(i)
      !STOP
    END IF
  END DO

  ALLOCATE(spaces(t_len))

  ! Allocate identity jump operator
  ALLOCATE(spaces(1)%evals(n)) 
  ALLOCATE(spaces(1)%evi1(n))
  ALLOCATE(spaces(1)%evi2(n))
  spaces(1)%occupation = n
  spaces(1)%dims = n
  spaces(1)%deleps = 0.0d0
  DO i = 1, n
    spaces(1)%evi1(i) = i
    spaces(1)%evi2(i) = i
  END DO
  DO l = 1, nbath
    spaces(1)%thpl(l) = get_thpl_element(0.0d0,l,0.0d0,0.0d0)
  END DO

  ! Iterate through every single difference
  k = 2  ! Index into 'spaces'; first space is the identity
  DO i = 1, n
    DO j = 1, n
      ! Not part of identity operator
      IF (i .NE. j) THEN
        ! Allocate space
        ALLOCATE(spaces(k)%evals(1))
        ALLOCATE(spaces(k)%evi1(1))
        ALLOCATE(spaces(k)%evi2(1))
        spaces(k)%dims = 1
        ! Fill information
        spaces(k)%evi1(1) = i 
        spaces(k)%evi2(1) = j 
        spaces(k)%occupation = 1
        spaces(k)%deleps = evals(i) - evals(j)
        DO l = 1, nbath
          spaces(k)%thpl(l) = get_thpl_element(evals(i)-evals(j),l,0.0d0,0.0d0)
        END DO
        k = k + 1

      END IF 
    END DO
  END DO
END SUBROUTINE

! Calculate the Lamb Shift Hamiltonian, sum_omega
! sum_ab S A*(omega) A(omega) as seen in the reference
! paper from the jump operators for the system
SUBROUTINE get_Hls(d,spaces,Hls)
  INTEGER, INTENT(IN) :: d  ! Number of jump operators
  TYPE(sspace), DIMENSION(d), INTENT(IN) :: spaces  ! Jump operators
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: Hls  ! Output matrix
  COMPLEX*16 :: factor
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: temp1, temp2, temp
  INTEGER :: i, j, ind

  ALLOCATE(temp1(n,n))
  ALLOCATE(temp2(n,n))
  ALLOCATE(temp(n,n))

  Hls = 0.0d0
  DO i=1,d
    DO j = 1, nbath
      factor = 1.0d0*AIMAG(spaces(i)%thpl(j))

      ! Well, if there is only one entry in the matrix, as is 
      ! usually the case, then we know what it will look like times
      ! its conjugate, just one diagonal entry at (X_ij* * X_ji) => (i,i)
      IF (spaces(i)%Aomega%nentries .EQ. 1) THEN
        ind = spaces(i)%Aomega%cols(1)
        Hls(ind,ind) = Hls(ind,ind) + factor*(ABS(spaces(i)%Aomega%values(1,j))**2.0d0)
      ELSE  ! Multiple entries require manual multiplication
        CALL sparce_to_full(spaces(i)%Aomega,temp1,j)
        temp2 = CONJG(TRANSPOSE(temp1))
        temp = MATMUL(temp2,temp1)
        ! Multiple baths are all summed up together to create Hls
        Hls = Hls + factor*temp
      END IF  
    END DO
  END DO
  DEALLOCATE(temp1)
  DEALLOCATE(temp2)
  DEALLOCATE(temp)
END SUBROUTINE

! Get Hco (coherent/decay Hamiltonian) for quantum jumps
! When the full secular approximation is applied, Hco
! should be diagonal
SUBROUTINE get_hco(d,Hs,spaces,hco,Hls)
  INTEGER, INTENT(IN) :: d  ! Number of jump operators
  TYPE(sspace), DIMENSION(d), INTENT(IN) :: spaces  ! Jump operators
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: Hs, Hls  ! Original Hamiltonian, Lambshift Hamiltonian
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: temp, temp2 ! Intermediate factors
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: hco  ! Coherent decay hamiltonian
  COMPLEX*16 factor
  INTEGER :: i, j, ind

  ALLOCATE(temp(n,n))
  ALLOCATE(temp2(n,n))
  hco = 0.0d0
  DO i=1,d
    DO j = 1, nbath
      factor = REAL(REAL(spaces(i)%thpl(j)))
      ! Same deal as in Hls calculation; if there is only
      ! one entry, skip MATMUL and put it in the appropriate
      ! diagonal
      IF (spaces(i)%Aomega%nentries .EQ. 1) THEN
        ind = spaces(i)%Aomega%cols(1)
        hco(ind,ind) = hco(ind,ind) - CMPLX(0.0d0,0.5d0)*factor*2.0d0*&
                       (ABS(spaces(i)%Aomega%values(1,j))**2.0d0)
      ELSE
        CALL sparce_to_full(spaces(i)%Aomega,temp,j)
        temp2 = MAtMUL(CONJG(TRANSPOSE(temp)),temp)
        hco = hco - CMPLX(0.0d0,0.5d0)*factor*2.0d0*temp2
      END IF
    END DO
  END DO
  hco = (hco + Hls) /(hbar**2.0d0)
  hco = hco + Hs/hbar

  DO i = 1, n  ! Could be non-diagonal when partial secular approximations were made
    DO j = 1, n  ! Now it is just a problem indicating assumption failure.
      IF ((ABS(hco(i,j)) .GE. 1.0d-10) .AND. (i .NE. j)) THEN
        non_diag_hco = 1  ! Global flag for basis change needed
        WRITE(*,*) "NON-DIAGAONAL HCO"
        OPEN(89,FILE="non_diag_hco.txt")
        WRITE(89,*) hco(i,j), i, j
        STOP
      END IF
    END DO
  END DO

  DEALLOCATE(temp)
  DEALLOCATE(temp2)
END SUBROUTINE

!  Normalize a wavefunction; a wfn of 0 gets 0 back
SUBROUTINE normalize(wfn)
  COMPLEX*16, DIMENSION(n,1), INTENT(INOUT) :: wfn
  IF (norm(wfn) .LT. 1.0d-18) THEN
    wfn = 0.0d0  ! Causes less trouble than a bunch of NaNs
  ELSE
    wfn = wfn/norm(wfn)
  END IF
END SUBROUTINE

! Return <\psi|\psi> of a wavefunction \psi
REAL*8 FUNCTION normsq(wfn)
  COMPLEX*16, DIMENSION(n,1), INTENT(IN) :: wfn
  INTEGER i
  normsq = 0.0d0
  DO i = 1, n
    normsq = wfn(i,1)*CONJG(wfn(i,1)) + normsq
  END DO
END FUNCTION

! Give a wavefunction vector norm
REAL*8 FUNCTION norm(wfn)
  COMPLEX*16, DIMENSION(n,1), INTENT(IN) :: wfn
  INTEGER i
  norm = 0.0d0
  DO i = 1, n
    norm = wfn(i,1)*CONJG(wfn(i,1)) + norm
  END DO
  norm = SQRT(norm)
END FUNCTION

! Sparce matrix multiplication
! SPAM*column vector; select bath
! number via k
SUBROUTINE spam_times_vec(spm,vec,veco,k,norm)
  INTEGER, INTENT(IN) :: k
  TYPE(SPAM), INTENT(IN) :: spm  ! Sparce matrix
  REAL*8, INTENT(OUT), OPTIONAL :: norm  ! Result norm
  COMPLEX*16, DIMENSION(n,1), INTENT(IN) :: vec  ! In vector
  COMPLEX*16, DIMENSION(n,1), INTENT(OUT) :: veco  ! Out vector
  REAL*8 :: total
  INTEGER i
  veco = 0.0d0
  total = 0.0d0
  DO i = 1, spm%nentries
    veco(spm%rows(i),1) = veco(spm%rows(i),1) + vec(spm%cols(i),1)*spm%values(i,k) 
    total = total + (ABS(veco(spm%rows(i),1)))**2.0d0
  END DO
  IF (PRESENT(norm)) THEN
    norm = total
  END IF
END SUBROUTINE

! Once upon a time, full secular approximations weren't
! made and subspaces were of varying size. This machinery
! is still in use.
! Inflate spam storage in each subspace and calculate
! the proper sparse operator to hold within based
! on eigenvectors we have; this could be rewritten
! to happen earlier and thus result in less work...
! but oh well. If I change it now I might break
! something so it's done in this second step.
SUBROUTINE make_sparce_operator(d,spaces,coupling)
  INTEGER, INTENT(IN) :: d  ! Number of spaces
  TYPE(sspace), DIMENSION(d), INTENT(INOUT) :: spaces  ! Spaces
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: coupling  ! Coupling operator
  INTEGER :: i, j, k
  ! Give each subspace the number of row/column/bath entries it needs
  ! ie one for every entry in that energy subspace
  DO i = 1, d
    ALLOCATE(spaces(i)%Aomega%cols(spaces(i)%dims))
    ALLOCATE(spaces(i)%Aomega%rows(spaces(i)%dims))
    ALLOCATE(spaces(i)%Aomega%values(spaces(i)%dims,nbath))
    spaces(i)%Aomega%nentries = spaces(i)%dims
    IF (DEBUG .EQV. .TRUE.) THEN
      WRITE(*,*) "Giving dimension", spaces(i)%dims, "to space", i
    END IF
    DO j = 1, spaces(i)%dims
      spaces(i)%Aomega%cols(j) = spaces(i)%evi2(j)
      spaces(i)%Aomega%rows(j) = spaces(i)%evi1(j)
      DO k = 1, nbath
        spaces(i)%Aomega%values(j,k) = coupling(spaces(i)%Aomega%cols(j),&
          spaces(i)%Aomega%rows(j),k)
      END DO
    END DO
  END DO
END SUBROUTINE

! Take a sparce matix, sm, and return a full matrix, fm
! Integer l selects the matrix values desired, i.e. it
! is the bath index, essentially
SUBROUTINE sparce_to_full(sm,fm, l)
  INTEGER, INTENT(IN) :: l
  TYPE(spam), INTENT(IN) :: sm
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: fm
  INTEGER :: i, j, k
  fm = 0.0d0
  DO i = 1, sm%nentries
    j = sm%rows(i)
    k = sm%cols(i)
    fm(j,k) = sm%values(i,l)
  END DO
END SUBROUTINE

END MODULE lindblad_coherent
