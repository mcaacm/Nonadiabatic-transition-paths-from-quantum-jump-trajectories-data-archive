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

! Reads specifications in the file in_wfns.txt and eigenstate
! wavefunction specifications in params.f90 and plots the wavefunctions
! in CI_data/wfn_plots/

PROGRAM plot_wf

USE parameters
USE prequel
USE set_up_H
USE rk_ns_utils

IMPLICIT NONE
INTEGER :: i, j, k, reason  ! Loop inexes, read error indicator
CHARACTER(LEN=100) :: fname   ! File name
CHARACTER(LEN=20) :: fstring  ! More file name stuff
REAL*8, DIMENSION(sz,sz) :: hf  ! Hermite polynomial factors
REAL*8, DIMENSION(sz,res) :: grid  ! Hermite evaluation on grid points
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: sigma, H, evec  ! Denisty matrix, Hamiltonian, energy eigenvectors
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:,:) :: couple_op  ! Coupling operator
REAL*8, DIMENSION(n) :: estate_probs  ! Unneeded setup
COMPLEX*16, DIMENSION(n) :: eval  ! Energy eigenvalues
COMPLEX*16, DIMENSION(m,1) :: wfn_m  ! Wavefunctions in full vs reduced basis
COMPLEX*16, DIMENSION(n,1) :: wfn
COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: translated_evecs  ! Unneeded setup
REAL*8, DIMENSION(res,res,6) :: plot, plot9, plot10, plot5, plot14, plot15  ! For wavefunction plotting
REAL*8 :: low, high, integral  ! Plot grid limits and an integration temporary
REAL*8, DIMENSION(res,2) :: out_grid  
REAL*8, DIMENSION(2*n) :: temp

IF (bound1 .LT. 0 .OR. bound2 .LT. 0 .OR. bound2 .LT. bound1) THEN
  WRITE(*,*) "Invalid plot bounds:", bound1, bound2
  STOP
END IF
IF (bound3 .LT. 0 .OR. bound4 .LT. 0 .OR. bound4 .LT. bound3) THEN
  WRITE(*,*) "Invalid plot bounds:", bound3, bound4
  STOP
END IF


ALLOCATE(sigma(n,n))
ALLOCATE(H(n,n))
ALLOCATE(evec(m,n))
ALLOCATE(couple_op(n,n,nbath))
ALLOCATE(translated_evecs(n,n))
! Allocate global storage found in prequel
ALLOCATE(G_min_global(n,n,nbath))
ALLOCATE(G_pls_global(n,n,nbath))
ALLOCATE(G_global(n,n,nbath))
! Allocate global storage in Setup_H
ALLOCATE(Q_x_g(n,n))
ALLOCATE(Q_y_g(n,n))


! Setup system
CALL setup(H,sigma,evec,eval,translated_evecs,estate_probs,couple_op)

! I want to plot the initially excited wavefunction.
OPEN(11,FILE="in_wfns.txt")
WRITE(11,*) org_out(translated_evecs(1:n,1:1))
CLOSE(11)

! Fill a grid of Hermite polynomials for the basis
low = -10.0d0
high = 10.0d0
CALL fill_hermite(sz,hf)
CALL fill_wf_grid(sz,res,hf,grid,low,high)
WRITE(*,*) "Filled hermite grid"
WRITE(*,*) hf

! Plot lower energy eigenstates of relevance
DO i = bound1, bound2
  wfn = 0.0d0
  wfn(i,1) = 1.0d0
  wfn_m = MATMUL(evec,wfn)
  ! Fil and integrate a grid for each eigenstate
  CALL fill_wf_plot(sz,res,grid,low,high,wfn_m,plot)
  integral = integrate_grid(res,plot,low,high)
  WRITE(*,*) "Integrated density eigen fn ", i, "is", integral

  IF (i .EQ. 5) THEN 
    plot5 = plot
  ELSE IF (i .EQ. 9) THEN 
    plot9 = plot
  ELSE IF (i .EQ. 10) THEN
    plot10 = plot
  ELSE IF (i .EQ. 14) THEN
    plot14 = plot
  ELSE IF (i .EQ. 15) THEN
    plot15 = plot
  END IF

  IF (i .LT. 10) THEN
    fstring = "(A24,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A24,I2,A4)"
  ELSE 
    fstring = "(A24,I3,A4)"
  END IF 
   

  WRITE(fname,fstring) "CI_data/wfn_plots/ef_2d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL print_wf_plot(43,plot,low,high,res)
  CLOSE(43)
  WRITE(fname,fstring) "CI_data/wfn_plots/c2_1d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL integrate_out_dimension(2,res,plot,low,high,out_grid)
  DO j = 1, res
    WRITE(43,*) j, low + ((high-low)/res)*(j-1), out_grid(j,1), out_grid(j,2)
  END DO
  CLOSE(43)
  WRITE(fname,fstring) "CI_data/wfn_plots/c1_1d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL integrate_out_dimension(1,res,plot,low,high,out_grid)
  DO j = 1, res
    WRITE(43,*) j, low + ((high-low)/res)*(j-1), out_grid(j,1), out_grid(j,2)
  END DO
  CLOSE(43)

END DO


! Check normalization of individual eigenstates
! Plot higher energy eigenstates of relevance
DO i = bound3, bound4
  wfn = 0.0d0
  wfn(i,1) = 1.0d0
  wfn_m = MATMUL(evec,wfn)
  CALL fill_wf_plot(sz,res,grid,low,high,wfn_m,plot)
  integral = integrate_grid(res,plot,low,high)
  WRITE(*,*) "Integrated density eigen fn ", i, "is", integral

  IF (i .EQ. 5) THEN 
    plot5 = plot
  ELSE IF (i .EQ. 9) THEN 
    plot9 = plot
  ELSE IF (i .EQ. 10) THEN
    plot10 = plot
  ELSE IF (i .EQ. 14) THEN
    plot14 = plot
  ELSE IF (i .EQ. 15) THEN
    plot15 = plot
  END IF

  IF (i .LT. 10) THEN
    fstring = "(A24,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A24,I2,A4)"
  ELSE 
    fstring = "(A24,I3,A4)"
  END IF 
   

  WRITE(fname,fstring) "CI_data/wfn_plots/ef_2d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL print_wf_plot(43,plot,low,high,res)
  CLOSE(43)
  WRITE(fname,fstring) "CI_data/wfn_plots/c2_1d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL integrate_out_dimension(2,res,plot,low,high,out_grid)
  DO j = 1, res
    WRITE(43,*) j, low + ((high-low)/res)*(j-1), out_grid(j,1), out_grid(j,2)
  END DO
  CLOSE(43)
  WRITE(fname,fstring) "CI_data/wfn_plots/c1_1d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  CALL integrate_out_dimension(1,res,plot,low,high,out_grid)
  DO j = 1, res
    WRITE(43,*) j, low + ((high-low)/res)*(j-1), out_grid(j,1), out_grid(j,2)
  END DO
  CLOSE(43)

END DO

! Read specific wavefunction data from "in_wfns.txt" and plot as
! CI_data/wfn_plots/wfn_2d_i.txt
i = 0
OPEN(11,FILE="in_wfns.txt")
DO
  READ(11,*,IOSTAT=reason) (temp(k), k=1, 2*n)

  IF (reason .NE. 0) THEN
    WRITE(*,*) "Read error:", reason, "exiting loop after", i, "reads"
    EXIT
  END IF

  DO j = 1, n
    wfn(j,1) = CMPLX(temp(2*j - 1),temp(2*j))
  END DO
  wfn_m = MATMUL(evec,wfn)
  !wfn_m = 0.0d0
  !wfn_m(m/2 + 1,1) = 1.0d0
  CALL fill_wf_plot(sz,res,grid,low,high,wfn_m,plot)

  IF (i .LT. 10) THEN
    fstring = "(A24,I1,A4)"
  ELSE IF (i .LT. 100) THEN
    fstring = "(A24,I2,A4)"
  ELSE 
    fstring = "(A24,I3,A4)"
  END IF 

  WRITE(fname,fstring) "CI_data/wfn_plots/wf_2d_", i, ".txt"
  OPEN(43,FILE=TRIM(fname))
  integral = integrate_grid(res,plot,low,high)
  CALL print_wf_plot(43,plot,low,high,res)
  WRITE(43,*) "Integrated read wfn to", integral
  CLOSE(43)

  i = i + 1
END DO


CLOSE(11)


END PROGRAM

