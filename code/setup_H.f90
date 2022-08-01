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

! Sets up the Hamiltonian, initial vector to propagate from in
! Lindblad, coupling operators, etc.

! Conical intersection setup follows Chen, Lipeng, Maxim F. Gelin, Vladimir Y. Chernyak, 
! Wolfgang Domcke, and Yang Zhao. "Dissipative dynamics at conical intersections: simulations 
! with the hierarchy equations of motion method." Faraday discussions 194 (2016): 61-80.

MODULE set_up_H
USE parameters
USE prequel
USE rk_ns_utils

IMPLICIT NONE

CONTAINS 

! Output information about a conical intersection
! density matrix; gives the time, the diabatic popluations
! and the coupling/tuning positions
SUBROUTINE output_info_CI_dmat(time,fd,dmat,evec,evec_inv)
  INTEGER, INTENT(IN) :: fd
  COMPLEX*16, INTENT(IN), DIMENSION(n,n) :: dmat  ! Density matrix
  COMPLEX*16, INTENT(IN), DIMENSION(m,n) :: evec  ! Energy eigenvectors
  COMPLEX*16, INTENT(IN), DIMENSION(n,m) :: evec_inv
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: dmat_temp
  REAL*8, INTENT(IN) :: time
  REAL*8 :: Q_x, Q_y, temp, tm1, tm2

  ALLOCATE(dmat_temp(m,m))
  Q_x = trace_mat(n,MATMUL(Q_x_g,dmat),1,n)
  Q_y = trace_mat(n,MATMUL(Q_y_g,dmat),1,n)
  CALL from_basis(dmat,dmat_temp,evec,evec_inv)
  tm1 = trace_mat(m,dmat_temp,1,m/2)
  tm2 = trace_mat(m,dmat_temp,m/2+1,m)
  temp = 1.0d0/(tm1 + tm2)
  WRITE(fd,*) time, trace_mat(m,dmat_temp,1,m/2), trace_mat(m,dmat_temp,m/2+1,m), Q_x, Q_y 
  DEALLOCATE(dmat_temp)
END SUBROUTINE


! Call the setup routines for specific types; when multiple different models were supported,
! this had a much bigger role
SUBROUTINE setup(H, sigma, evec, eval, H_trans, estate_probs,couple_op,Qt_full)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: H, sigma  ! Hamiltonian, starting denisty matrix
  COMPLEX*16, DIMENSION(m,n), INTENT(OUT) :: evec  ! Energy eigenvectors
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: couple_op  ! System Bath coupling operator
  COMPLEX*16, DIMENSION(m,m), INTENT(OUT), OPTIONAL :: Qt_full  ! Full basis position operator
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval  ! Energy eigenvalues
  ! Boltzman likelihood to start in oscillator eigenstate; used by thermal init and DEPRECATED
  REAL*8, INTENT(OUT), DIMENSION(n) :: estate_probs
  ! Whose translation to the unshifted oscillator basis
  ! looks like this; each vector is its own column; used by thermal init and DEPRECATED
  COMPLEX*16, INTENT(OUT), DIMENSION(n,n) :: H_trans

  IF (PRESENT(Qt_full) .AND. run_type .NE. CI) THEN
    WRITE(*,*) "Error: optional argument Qt_full not implemented for this system type (function call to setup)."
    STOP
  END IF 

  IF (run_type .EQ. CI) THEN
    IF (nbath .EQ. 2) THEN
      IF (PRESENT(Qt_full)) THEN
        CALL setup_H_CI_2bath(H,sigma,evec,eval,H_trans,estate_probs,couple_op,Qt_full)
      ELSE
        CALL setup_H_CI_2bath(H,sigma,evec,eval,H_trans,estate_probs,couple_op)
      END IF
    ELSE
      WRITE(*,*) "Error. One bath no longer supported. Use two baths; set second eta to zero."
      STOP
    END IF
  ELSE
    WRITE(*,*) "Unrecognized runtype"
    STOP
  END IF
END SUBROUTINE

! Set up a conical intersection Hamiltonian.
! Used to setup two baths, one on the coupling mode and 
! one on the tuning mode.
SUBROUTINE setup_H_CI_2bath(H, sigma, evec, eval, H_trans, estate_probs,couple_op,Qt_full)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: H, sigma  ! Hamiltonian, initial density matrix
  COMPLEX*16, DIMENSION(m,n), INTENT(OUT) :: evec  ! H eigenvectors
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(OUT) :: couple_op  ! Coupling operator for each bath
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval  ! Hamiltonian eigenvalues
  ! Boltzman likelihood to start in oscillator eigenstate (no longer used)
  REAL*8, INTENT(OUT), DIMENSION(n) :: estate_probs  ! DEPRECATED; first is always 1
  COMPLEX*16, OPTIONAL, INTENT(OUT), DIMENSION(m,m) :: Qt_full  ! An optional full Q operator
  ! Whose translation to the unshifted oscillator basis
  ! looks like this; each vector is its own column
  COMPLEX*16, INTENT(OUT), DIMENSION(n,n) :: H_trans  ! DEPRECATED Only the first vector is used now.
  ! Operators used during setup
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: Q_c, Q_t, H_temp, proj1m, proj2m
  ! For diagonalization of Q_c to get |Q_c| if needed
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: Qc_evec, Qc_sqr, Qc_abs, Q_temp
  COMPLEX*16, DIMENSION(m/2) :: Qc_eval
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: proj1n, proj2n
  ! Bounds and loop indexes to fill the matrices
  INTEGER :: k, k_p, nc, nt, nc_p, nt_p, index_1, index_2, i, bound1, bound2
  REAL*8 :: total
  TYPE(wfunc) :: wfn

  ! For building in the larger Hilbert space before reduction
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: H_m, sigma_m, evec_m
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: evec_inv
  COMPLEX*16, DIMENSION(m) :: evals_m
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: H_trans_m, tvecs
  REAL*8, DIMENSION(:), ALLOCATABLE :: tprobs

  ! Set up a one dimensional position operator 
  ALLOCATE(Q_temp(m,m))
  ALLOCATE(tprobs(m/2))
  ALLOCATE(tvecs(m/2,m/2))
  ALLOCATE(H_m(m,m))
  ALLOCATE(sigma_m(m,m))
  ALLOCATE(evec_m(m,m))
  ALLOCATE(evec_inv(n,m))
  ALLOCATE(H_trans_m(m,m))
  ALLOCATE(proj1n(n,n))
  ALLOCATE(proj2n(n,n))
  ALLOCATE(Q_c(m,m))
  ALLOCATE(Q_t(m,m))
  ALLOCATE(H_temp(m,m))
  ALLOCATE(proj1m(m,m))
  ALLOCATE(proj2m(m,m))


  ! Setup a Q matrix
  DO nc = 1, m 
    DO nc_p = 1, m 
      Q_temp(nc,nc_p) = del(nc,nc_p+1)*SQRT((nc_p)/2.0d0) + &
        del(nc,nc_p-1)*SQRT((nc_p - 1.0)/2.0d0)
     END DO
  END DO

  ! Set up the |Q_c| matrix if absolute coupling was requested
  IF (abs_coupling .EQ. 1) THEN
    ALLOCATE(Qc_evec(m/2,m/2))
    ALLOCATE(Qc_abs(m/2,m/2))
    ALLOCATE(Qc_sqr(m/2,m/2))

    Qc_sqr = Q_temp(1:m/2,1:m/2)


    CALL diagonalize_m(m/2,Qc_sqr,Qc_eval,Qc_evec)
    Qc_sqr = 0.0d0

    DO nc = 1, m/2
      Qc_sqr(nc,nc) = ABS(REAL(REAL(Qc_eval(nc))))
    END DO
    ! Reassemble into the absolute value, named Qc_sqr because it
    ! used to be the square root of the square but this is a neater construction
    Qc_sqr = MATMUL(Qc_evec,MATMUL(Qc_sqr,CONJG(TRANSPOSE(Qc_evec))))
  END IF

  H_m = 0.0d0
  H_trans_m = 0.0d0
  H = 0.0d0
  Q_c = 0.0d0
  Q_t = 0.0d0
  ! Set up full H as in Faraday. Discus. 
  DO k = 1, 2
    DO k_p = 1, 2
      DO nc = 1, dim_nc
        DO nc_p = 1, dim_nc
          DO nt = 1, dim_nt
            DO nt_p = 1, dim_nt
              index_1 = (k-1)*(dim_nc*dim_nt) + (nc-1)*dim_nt + nt
              index_2 = (k_p-1)*(dim_nc*dim_nt) + (nc_p-1)*dim_nt + nt_p
              IF (abs_coupling .EQ. 1) THEN   ! Couple to |Q_c|
                H_m(index_1,index_2) = del(k,k_p)*del(nc,nc_p)*del(nt,nt_p)*  &
                     (exe(k) + hbar*omegas_c*(nc - 0.5) + hbar*omegas_t*(nt - 0.5)) & 
                      + del(k,k_p)*del(nc,nc_p)*kappa_t(k)*(del(nt,nt_p+1)*SQRT((nt_p)/2.0d0) + &
                      del(nt,nt_p-1)*SQRT((nt_p - 1.0)/2.0d0)) + (del(k,k_p-1) + del(k,k_p+1))*&
                      del(nt,nt_p)*lambda_s*Qc_sqr(nc,nc_p)  ! Get the absolute coupling matrix entry
              ELSE  ! Couple to Q_c
                H_m(index_1,index_2) = del(k,k_p)*del(nc,nc_p)*del(nt,nt_p)*  &
                     (exe(k) + hbar*omegas_c*(nc - 0.5) + hbar*omegas_t*(nt - 0.5)) & 
                      + del(k,k_p)*del(nc,nc_p)*kappa_t(k)*(del(nt,nt_p+1)*SQRT((nt_p)/2.0d0) + &
                      del(nt,nt_p-1)*SQRT((nt_p - 1.0)/2.0d0)) + (del(k,k_p-1) + del(k,k_p+1))*&
                      del(nt,nt_p)*lambda_s*Q_temp(nc,nc_p)

              END IF

              ! Coupling mode position
              Q_c(index_1,index_2) = del(k,k_p)*del(nt,nt_p)*Q_temp(nc,nc_p) 
              ! Tuning mode position
              Q_t(index_1,index_2) = del(k,k_p)*del(nc,nc_p)*Q_temp(nt,nt_p)
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO


  H_temp = H_m


  IF (PRESENT(Qt_full)) THEN
    Qt_full = Q_t
  END IF

  CALL diagonalize_m(m,H_m,evals_m,evec_m)

  evec = evec_m(1:m,1:n)
  eval = evals_m(1:n)
  H = 0.0d0
  DO k = 1, n
    H(k,k) = eval(k)
  END DO
  evec_inv = CONJG(TRANSPOSE(evec))
  ! Put global Q matrices in the eigenbasis

  CALL to_basis(Q_c,Q_x_g,evec,evec_inv)
  CALL to_basis(Q_t,Q_y_g,evec,evec_inv)

  ! Setup initial density matrix with either thermal or vertical
  ! state in either well 2 or well 1
  IF (init_well .EQ. 1) THEN
    bound1 = 1
    bound2 = m/2
  ELSE ! Initialize in well 2
    bound1 = m/2 + 1
    bound2 = m
  END IF

  sigma_m = 0.0d0
  total = 0.0d0
  estate_probs = 0.0d0
  H_trans = 0.0d0
  ! Stationary thermal initalization initializes lowest energy eigenstate at the moment
  IF (init_stat_therm .EQ. 1) THEN
    H_trans = 0.0d0
    H_trans(init_es,1) = 1.0d0
    estate_probs = 0.0d0
    estate_probs(1) = 1.0d0
    sigma(init_es,init_es) = 1.0d0
  ELSE
    estate_probs = 0.0d0
    estate_probs(1) = 1.0d0
    H_trans_m(bound1,1) = 1.0d0
    H_trans(1:n,1) = MATMUL(evec_inv,H_trans_m(1:m,1))
    CALL wfn_to_rho_2(H_trans(1:n,1),sigma)
    sigma_m(bound1,bound1) = 1.0d0
  END IF
   

  IF (init_stat_therm .EQ. 0) THEN
    CALL to_basis(sigma_m,sigma,evec,evec_inv)
  END IF

  DO k = 1, n
    IF (REAL(REAL(sigma(k,k))) .LT. -1.0d-8) THEN
      WRITE(*,*) "Error: negative diagonal entry in density matrix."
      STOP
    END IF
  END DO

  CALL from_basis(sigma,sigma_m,evec,evec_inv)
  CALL to_basis(sigma_m,sigma,evec,evec_inv)

  ! Set two separate baths
  couple_op(1:n,1:n,1) = Q_x_g
  couple_op(1:n,1:n,2) = Q_y_g

  proj1m = 0.0d0
  proj2m = 0.0d0
  DO i = 1, m/2
    proj1m(i,i) = 1.0d0 
  END DO
  DO i = m/2 + 1, m
    proj2m(i,i) = 1.0d0
  END DO
  CALL to_basis(proj1m,proj1n,evec,evec_inv)
  CALL to_basis(proj2m,proj2n,evec,evec_inv)

  OPEN(62,FILE="eigenvector_info.txt")
  OPEN(63,FILE="plot_ev.txt")
  OPEN(64,FILE="evic.txt")
  ! Write information about the eigenvectors
  CALL init_wfn(wfn)
  DO i = 1, n
    wfn%fe = 0.0d0
    wfn%fe(i,1) = 1.0d0
    WRITE(64,*) op_avg(wfn,couple_op(1:n,1:n,1)), op_avg(wfn,couple_op(1:n,1:n,2)), op_avg(wfn,H), &
                op_avg(wfn,proj1n),  op_avg(wfn,proj2n)

    WRITE(62,*) "Eigenvector ", i, "pos1", op_avg(wfn,couple_op(1:n,1:n,1)), "pos2", op_avg(wfn,couple_op(1:n,1:n,2)),&
               "energy", op_avg(wfn,H), "p in 1", op_avg(wfn,proj1n), "p in 2", op_avg(wfn,proj2n)
    WRITE(63,*) i, -100.0, op_avg(wfn,H)
    WRITE(63,*) i, 100.0, op_avg(wfn,H)
    WRITE(63,*) ""
    WRITE(63,*) ""

  END DO
  CALL destroy_wfn(wfn)
  CLOSE(62)
  CLOSE(63)
  CLOSE(64)


  DEALLOCATE(H_m)
  DEALLOCATE(sigma_m)
  DEALLOCATE(evec_m)
  DEALLOCATE(evec_inv)
  DEALLOCATE(H_trans_m)
  DEALLOCATE(proj1n)
  DEALLOCATE(proj2n)
  DEALLOCATE(Q_c)
  DEALLOCATE(Q_t)
  DEALLOCATE(H_temp)
  DEALLOCATE(proj1m)
  DEALLOCATE(tprobs)
  DEALLOCATE(tvecs)
  DEALLOCATE(proj2m)

  IF (ALLOCATED(Qc_sqr)) THEN
    DEALLOCATE(Qc_sqr)
  END IF
  IF (ALLOCATED(Qc_evec)) THEN
    DEALLOCATE(Qc_evec)
  END IF
  IF (ALLOCATED(Qc_abs)) THEN
    DEALLOCATE(Qc_abs)
  END IF

END SUBROUTINE



END MODULE
