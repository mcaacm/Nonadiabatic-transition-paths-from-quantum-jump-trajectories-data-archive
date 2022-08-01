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

! Runs secular Redfield or Lindblad Jump dynamics as requested in params.f90

PROGRAM testing
USE parameters
USE prequel
USE rk_ns_utils
USE lindblad_coherent
USE lindblad_jumps
USE set_up_H
IMPLICIT NONE

! Density matrix, Hamiltonian which will be diagonalized
! right eigenvectors of H, inverse of evec matrix
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma, H
! Hamiltonian eigenvectors and inverse, matrices for secular Redfield propagation
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: evec, evec_inv, diag_mat, odiag_mat, diag_sig
! Starting vectors for the lindblad jumps routine
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: translated_evecs
! Fourier transform of correlation functions, coupling operator
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:,:) :: theta_plus, theta_minus, couple_op
! Right eigenvalues of Hamiltonian 
COMPLEX*16, DIMENSION(n) :: eval, H_d
! Probability to begin jump run from a given vector when more than one is possible
REAL*8, DIMENSION(n) :: estate_probs
! Indexes, number of steps to run in secular redfield
INTEGER :: i, j, k, stop_point

! Allocate memory
ALLOCATE(sigma(n,n))
ALLOCATE(H(n,n))
ALLOCATE(evec_inv(n,m))
ALLOCATE(evec(m,n))
ALLOCATE(translated_evecs(n,n))
ALLOCATE(theta_plus(n,n,nbath))
ALLOCATE(theta_minus(n,n,nbath))
ALLOCATE(couple_op(n,n,nbath))
! Allocate global storage found in prequel
ALLOCATE(G_min_global(n,n,nbath))
ALLOCATE(G_pls_global(n,n,nbath))
ALLOCATE(G_global(n,n,nbath))
! Allocate global storage in Setup_H
ALLOCATE(Q_x_g(n,n))
ALLOCATE(Q_y_g(n,n))
! These are never deallocated but it doesn't matter;
! They're important until the end


IF (DEBUG .EQV. .TRUE.) THEN
  WRITE(*,*) "Running with following parameters: omega_c",  omegac, "eta", eta, "bath type", bath_type, &
           "dim_nc", dim_nc, "dim_nt", dim_nt, "m", m, "n", n, "omegas_c", omegas_c, "omegas_t", omegas_t, &
           "kappa_t", kappa_t, "exe", exe, "lambda_s", lambda_s
END IF


H = (0.0d0, 0.0d0)
evec = (0.0d0, 0.0d0)
theta_plus = (0.0d0, 0.0d0)
sigma = (0.0d0,0.0d0)
eval = (0.0d0,0.0d0)

CALL setup(H,sigma,evec,eval,translated_evecs,estate_probs,couple_op)
CALL set_G(couple_op)
! Copy and invert the eigenvector matrix
evec_inv = CONJG(TRANSPOSE(evec))

! Number of propagation steps to run
stop_point = IDINT(((duration / dt) + 0.5d0))
DO i = 1, n
  H_d(i) = H(i,i)
END DO

OPEN(21,file="evolution.txt")

IF (lb_switch .EQ. 1) THEN
  WRITE(*,*) "Calling Lindblad driver"
  CALL drive_lindblad(H,eval,evec_inv,translated_evecs,estate_probs)

ELSE  ! Otherwise perform Redfield dynamics

  WRITE(*,*) "Getting theta plus."
  CALL get_Theta_plus(theta_plus, H, 0.0d0, 100.0d0)
  CALL get_Theta_minus(theta_minus, theta_plus)
  CALL set_G_plus(theta_plus)
  CALL set_G_minus(theta_minus)

  ! Perform secular evolution. Non-secular used to be an option, hence this switch remaining
  IF (secular .EQ. 1) THEN 
    OPEN(72,FILE="state_prop.txt")
    IF (n .GT. 16) THEN  ! Write out some eigenbasis populations
      WRITE(72,*) 0, 0, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
                    REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4))), REAL(REAL(sigma(5,5))), &
                    REAL(REAL(sigma(6,6))), REAL(REAL(sigma(7,7))), REAL(REAL(sigma(8,8))), &
                    REAL(REAL(sigma(9,9))), REAL(REAL(sigma(10,10))), REAL(REAL(sigma(11,11))), &
                    REAL(REAL(sigma(12,12))), REAL(REAL(sigma(13,13))), REAL(REAL(sigma(14,14))), &
                    REAL(REAL(sigma(15,15))), REAL(REAL(sigma(16,16))), REAL(REAL(sigma(17,17)))
    ELSE IF (n .GT. 2) THEN
      WRITE(72,*) 0, 0, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), REAL(REAL(sigma(3,3)))
    ELSE
      WRITE(72,*) 0, 0, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2)))
    END IF

    ALLOCATE(diag_mat(n,n))
    ALLOCATE(odiag_mat(n,n))
    ALLOCATE(diag_sig(n,1))

    CALL secular_rf_setup(diag_mat,odiag_mat,theta_plus,theta_minus,H_d)
    ! Diagonalize the diag_mat
    WRITE(*,*) "Setting up secular Redfield: diagonal pieces"
    WRITE(*,*) "Setting up secular Redfield: off diagonal pieces"
    DO k = 1, n
      diag_sig(k,1) = sigma(k,k)
    END DO
  
    DO i = 1, stop_point
      ! Updates the diagonal components
      CALL drho_dt_secular(diag_sig,diag_mat,dt)     
      ! Updates the non-diagonal components
      DO j = 1, n
        DO k = 1, n
          sigma(j,k) = sigma(j,k)*EXP(odiag_mat(j,k)*dt)
        END DO
      END DO
 
      ! Output
      IF (MOD(i,print_num) == 0) THEN
        DO j = 1, n
          sigma(j,j) = diag_sig(j,1)
        END DO
        ! Call output routine for the density matrix
        CALL output_info_CI_dmat(i*dt, 21, sigma, evec, evec_inv)
        IF (n .GT. 16) THEN
          ! Write out the energy eigenbasis populations for a few eigenvalues
          WRITE(72,*) i, i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), &
                      REAL(REAL(sigma(3,3))), REAL(REAL(sigma(4,4))), REAL(REAL(sigma(5,5))), &
                      REAL(REAL(sigma(6,6))), REAL(REAL(sigma(7,7))), REAL(REAL(sigma(8,8))), &
                      REAL(REAL(sigma(9,9))), REAL(REAL(sigma(10,10))), REAL(REAL(sigma(11,11))), &
                      REAL(REAL(sigma(12,12))), REAL(REAL(sigma(13,13))), REAL(REAL(sigma(14,14))), &
                      REAL(REAL(sigma(15,15))), REAL(REAL(sigma(16,16))), REAL(REAL(sigma(17,17)))
        ELSE IF (n .GE. 3) THEN
          WRITE(72,*) i, i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2))), REAL(REAL(sigma(3,3)))
        ELSE
          WRITE(72,*) i, i*dt, REAL(REAL(sigma(1,1))), REAL(REAL(sigma(2,2)))
        END IF
      END IF
    END DO

    CLOSE(72) 
  END IF

END IF  ! End of Redfield section as opposed to Lindblad section

END PROGRAM testing
