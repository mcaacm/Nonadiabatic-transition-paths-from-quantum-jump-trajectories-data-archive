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

! These are functions that all the programs need to call, either
! for integration or linear algebra etc. Global values are also
! defined in this file along with derived types. The
! Markov model analysis doesn't use this file, but everything
! else does.

MODULE prequel
USE parameters
IMPLICIT NONE

! Theta, G and Theta+/Theta- notation in this file follows the derivation in 
! Pollard, W. T. and Frisner, R. A. "Solutions of the Redfield equation for the dissipative quantum dynamics of 
! multilevel systems", J. Chem. Phys. 100, 1994.

! The derivation of the fast, secular Redfield as implemented here is based on
! Balzer, B. and Stock, G. "Modeling of decoherence and dissipation in nonadiabatic
! photoreactions by an effective-scaling nonsecular Redfield algorithm", Chem. Phys.
! 310, 2005.


REAL(KIND=4) :: qp_t_g ! Internal parameter for quadpack correlation integration
INTEGER :: cbi = 1  ! Internal parameter for quadpack integration; bath index
REAL(KIND=4) :: omega_cbi ! Internal parameter for quadpack correlation integration; set before calling integration functions
INTEGER, PARAMETER :: res = 501  ! Resolution of wavefunction plotting grids
INTEGER, PARAMETER :: sz = dim_nt  ! Degree of hermite polynomial to calculate for graphing wavefunctions

! Wavefunction; uncollapsed may have density in any number of eigenstates (sn = 0)
! but once collapsed there is only density in the state given by sn. Only that
! fe entry matters.
TYPE wfunc
  LOGICAL :: collapsed   ! Ended up not really being used in favor of checking sn but remains for historical reasons
  INTEGER :: sn  ! State number 0 for not collapsed, otherwise number of nonzero eigenstate
  COMPLEX*16, DIMENSION(n,1) :: fe  ! Full wavefunction entries
END TYPE

! Standard data for jumping along a trajectory
TYPE jump_data
  TYPE(wfunc) :: wfn  ! Wavefunction after jump
  REAL*8 :: time  ! Time at which jump occurred
  INTEGER :: id  ! Extra information not currently used but used to be important
  INTEGER :: estate  ! Eigenstate jumped to, or 0 if uncollapsed
END TYPE

! Global position operators Q_x_g is Q_c and Q_y_g is Q_t in the current setup
! Passing them around all the time would be
! inconvenient. These must be set up in the energy
! eigenbasis
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: Q_x_g, Q_y_g

! Flags indicating if the Redfield globals have been filled yet
! This is less clumsy than it was; I think it may be
! better than passing them all over the place like
! hot potatoes, even though globals are kind of disreputable
! Used for secular propagation
REAL*8, DIMENSION(nbath) :: g_min_flag = -1  ! g_min_global filled
REAL*8, DIMENSION(nbath) :: g_pls_flag = -1  ! g_pls_global filled
REAL*8, DIMENSION(nbath) :: g_flag = -1  ! g_global filled
! Global storage for G matrices in Redfield calculations
COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: G_min_global, G_pls_global, G_global

! Allocatable storage pointing to which jump operators are relevant
! for the given eigenstate. -1 indicates an unfilled location
! Used in LB jumps
INTEGER, DIMENSION(:,:), ALLOCATABLE :: jump_pointers
INTEGER, PARAMETER :: jump_pointer_dim = 2*n ! How large to make jump_pointers
INTEGER :: non_diag_hco = 0  ! Set to 1 to indicate that Hco is not diagonal
! With the current implementation, it always should be but it remains
! for historical reasons

! Sparce matrix structure for Lindblad operators which
! are either the identity or have a single entry with the
! current setup but that wasn't always the case
TYPE spam
  INTEGER :: nentries  ! Number of entries
  INTEGER, DIMENSION(:), ALLOCATABLE :: cols  ! Column index for each value
  INTEGER, DIMENSION(:), ALLOCATABLE :: rows  ! Row index for each value
  ! The second dimension allows for the G operator from each
  ! bath to have an entry in the matrix in order to save space
  ! because the column/row is the same but the value isn't
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: values
END TYPE

! General dimensional Lindblad subspace corresponding to one jump operator
TYPE sspace 
  INTEGER :: dims  ! Dimension of the subspace (should be 1 in full secular except for identity operators)
  ! Occupation and dims are now redundant but they were important in previous versions
  INTEGER :: occupation  ! How many entries have been added so far (used when adding information)
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: evals  ! Relevant energy eigenvalues involved (not used anymore)
  COMPLEX*16, DIMENSION(nbath) :: thpl  ! FT of correlation function value for that operator and each bath
  INTEGER, DIMENSION(:), ALLOCATABLE :: evi1  ! Eigenvector index 1 for the omega=(E2 - E1) calculation
  INTEGER, DIMENSION(:), ALLOCATABLE :: evi2  ! Eigenvector index 2
  REAL*8 :: deleps  ! Energy difference
  TYPE(spam) :: Aomega  ! The assembled jump operator as a sparce matrix
END TYPE sspace

CONTAINS

! Get the random seed out of the command argument
! for intializing the random number generator
INTEGER FUNCTION get_seed()
  INTEGER :: iarg, i
  CHARACTER(LEN=30) :: arg
 
  i = COMMAND_ARGUMENT_COUNT()
  IF (i .NE. 1) THEN
    get_seed = 1
  ELSE
    CALL GET_COMMAND_ARGUMENT(1, arg) 
    IF (LEN_TRIM(arg) .NE. 3) THEN
      WRITE(*,*) "Length of argument is equal to ", LEN_TRIM(arg)
      get_seed = 1
    ELSE 
      READ(arg,"(I3)") iarg
      get_seed = iarg
    END IF 
  END IF
  WRITE(*,*) "Using random seed base: ", get_seed
END FUNCTION

! Initialize Fortran built in random number generator
! with given seed seed_base_in
SUBROUTINE init_rand(seed_base_in)
  INTEGER, INTENT(IN) :: seed_base_in
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  integer :: n, i
  INTEGER :: seed_base

  CALL RANDOM_SEED(size = n)
  seed_base = seed_base_in
  IF (seed_base .EQ. 0) THEN
    seed_base = 1
  END IF
  ALLOCATE(seed(n))
  DO i = 1, n
    seed(i) = (i**3.0d0)*seed_base + 3*seed_base + 5
  END DO
  CALL RANDOM_SEED(put=seed)
END SUBROUTINE

! Returns a random number
! Call init_rand first
! Not currently using Mersenne Twister; xorshifts seem to
! be almost as good
REAL*8 FUNCTION get_rand_num()
  REAL :: rn 
  CALL RANDOM_NUMBER(rn)
  get_rand_num = rn
END FUNCTION

! Integer delta function delta(i,j)
INTEGER FUNCTION del(i, j)
  INTEGER, INTENT(IN) :: i, j
  IF (i .EQ. j) THEN
    del = 1
  ELSE
    del = 0
  END IF
END FUNCTION

! Print a matrix row by row, real matrix
! of dimension m
SUBROUTINE print_mat_r8(M,fnum)
  REAL*8, DIMENSION(:,:), INTENT(IN) :: M
  INTEGER, OPTIONAL, INTENT(IN) :: fnum  ! Optional file to print to
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

! Print a matrix row by row, complex matrix
! of dimension M
SUBROUTINE print_mat(M,fnum)
  COMPLEX*16, DIMENSION(:,:), INTENT(IN) :: M
  INTEGER, OPTIONAL, INTENT(IN) :: fnum  ! Optional file to print to
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

! Set global matrix required for Redfield tensor assembly
! G+ = G.*Theta+
SUBROUTINE set_G_plus(theta_plus)
  ! Fourier transform of correlation function
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_plus
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: G_pls
  INTEGER :: i
  ALLOCATE(G_pls(n,n))
  DO i = 1, nbath
    IF (G_flag(i) .NE. 1) THEN
      WRITE(*,*) "ERROR: G_flag", i, " not set."
      STOP
    END IF
    g_pls_flag(i) = 1
    G_pls = G_global(1:n,1:n,i)*theta_plus(1:n,1:n,i)
    G_pls_global(1:n,1:n,i) = G_pls
  END DO
  DEALLOCATE(G_pls)
END SUBROUTINE 

! Set global access to the coupling operators 
! all operators passed in (for each bath) must
! be in the eigenbasis already
SUBROUTINE set_G(oper)
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: oper
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: G
  INTEGER :: i
  ALLOCATE(G(n,n))
  DO i = 1, nbath
    G = oper(1:n,1:n,i)
    G_global(1:n,1:n,i) = G
    G_flag(i) = 1
  END DO
  DEALLOCATE(G)
END SUBROUTINE 

! Set global matrix required for Redfield tensor assembly
! G- = G.*theta-
SUBROUTINE set_G_minus(theta_minus)
  COMPLEX*16, DIMENSION(n,n,nbath), INTENT(IN) :: theta_minus
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: G_min
  INTEGER :: i
  ALLOCATE(G_min(n,n))
  DO i = 1, nbath
    IF (G_flag(i) .NE. 1) THEN
      WRITE(*,*) "ERROR: G_flag", i, " not set."
      STOP
    END IF
    G_min = G_global(1:n,1:n,i)*theta_minus(1:n,1:n,i)
    G_min_global(1:n,1:n,i) = G_min
    g_min_flag(i) = 1
  END DO
  DEALLOCATE(G_min)
END SUBROUTINE

! Diagonalize a real matrix, not symmetric, of dimension m
! Eigenvalues could be complex. They are returned in eval
! L is overwritten, no longer original matrix. Eigenvectors in evec
SUBROUTINE diagonalize_real_m(m, L, eval, evec)
  INTEGER, INTENT(IN) :: m
  REAL*8, DIMENSION(m,m), INTENT(INOUT) :: L
  COMPLEX*16, DIMENSION(m), INTENT(OUT) :: eval
  REAL*8, DIMENSION(m,m), INTENT(OUT) :: evec

  REAL*8, DIMENSION(m) :: wr, wi
  REAL*8, ALLOCATABLE, DIMENSION(:) :: work
  INTEGER :: info, i

  ALLOCATE(work(m*10))
  CALL dgeev('V','N',m,L,m,wr,wi,evec,m,evec,m,work,m*10,info)
  WRITE(*,*) "wr is ", wr, "wi is ", wi
  IF (info .NE. 0) THEN
    WRITE(*,*) "Diagonalization error: ", info
  END IF
  DEALLOCATE(work)
  DO i = 1, m
    eval(i) = CMPLX(wr(i),wi(i))
  END DO
END SUBROUTINE

! Diagonalize a matrix of dimension m; return right eval/evecs
! The input matirx must be hermitian. The L matrix
! is overwritten. 
SUBROUTINE diagonalize_m(m, L, eval, evec)
  ! how big is this matrix?
  INTEGER, INTENT(IN) :: m
  ! Matrix to be diagonalized
  COMPLEX*16, DIMENSION(m,m), INTENT(INOUT) :: L
  ! Right eigenvectors
  COMPLEX*16, DIMENSION(m,m), INTENT(OUT) :: evec
  ! Right eigenvalues
  COMPLEX*16, DIMENSION(m), INTENT(OUT) :: eval
  DOUBLE PRECISION, DIMENSION(m) :: eval_r
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: work
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rwork
  INTEGER :: i
  INTEGER :: info = 0
  ALLOCATE(work(m*10))
  ALLOCATE(rwork(3*m -2))
!  INTEGER :: LWORK = -1
  CALL zheev('V','L',m,L,m,eval_r,work,m*10,rwork,info)
  IF (info .NE. 0) THEN
    WRITE(*,*) "Diagonalization error: ", info
  END IF
  DO i = 1, m
    eval(i) = eval_r(i)
  END DO
  evec = L
  DEALLOCATE(work)
  DEALLOCATE(rwork)
END SUBROUTINE diagonalize_m

! Diagonalizes a matrix (hermitian or otherwise)
! of dimension m
! Generalized; works with hermitian or not hermitian
! Diagonalize a matrix; return right eval/evecs
SUBROUTINE diagonalize_m_g(m, L, eval, evec)
  ! how big is this matrix?
  INTEGER, INTENT(IN) :: m
  ! Matrix to be diagonalized
  COMPLEX*16, DIMENSION(m,m), INTENT(INOUT) :: L
  ! Right eigenvectors
  COMPLEX*16, DIMENSION(m,m), INTENT(OUT) :: evec
  ! Right eigenvalues
  COMPLEX*16, DIMENSION(m), INTENT(OUT) :: eval
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: work
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rwork
  COMPLEX*16, DIMENSION(m,1) :: temp_v
  COMPLEX*16 :: temp
  INTEGER :: i, j
  INTEGER :: info = 0
  ALLOCATE(work(m*10))
  ALLOCATE(rwork(m*2))
!  INTEGER :: LWORK = -1
  CALL zgeev('N','V',m,L,m,eval,evec,m,evec,m,work,m*10,rwork,info)
  IF (info .NE. 0) THEN
    WRITE(*,*) "Diagonalization error: ", info
  END IF
  temp = 0.0d0
  temp_v = 0.0d0
  DO j = 1, m
    DO i = 1, m - 1
      IF (ABS(eval(i)) > ABS(eval(i+1))) THEN
        temp_v(1:m,1) = evec(1:m,i)
        evec(1:m,i) = evec(1:m,i+1)
        evec(1:m,i+1) = temp_v(1:m,1)
        temp = eval(i)
        eval(i) = eval(i+1)
        eval(i+1) = temp
      END IF
    END DO
  END DO
  DEALLOCATE(work)
  DEALLOCATE(rwork)
END SUBROUTINE diagonalize_m_g

! Diagonalize a matrix; return right eval/evecs
! only works for hermitian matrices which are
! nxn with n specified in params.f90
SUBROUTINE diagonalize(L, eval, evec)
  ! Matrix to be diagonalized
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: L
  ! Right eigenvectors
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: evec
  ! Right eigenvalues
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval
  DOUBLE PRECISION, DIMENSION(n) :: evals_r
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: work
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rwork
  INTEGER :: i
  INTEGER :: info = 0
!  INTEGER :: LWORK = -1
  ALLOCATE(rwork(3*n - 2))
  ALLOCATE(work(n*10))
  CALL zheev('V','L',n,L,n,evals_r,work,n*10,rwork,info)
  IF (info .NE. 0) THEN
    WRITE(*,*) "Diagonalization error: ", info
  END IF
  DO i = 1, n
    eval(i) = evals_r(i)
  END DO
  evec = L
  DEALLOCATE(rwork)
  DEALLOCATE(work)
END SUBROUTINE diagonalize

! Diagonalize a matrix; return right eval/evecs
! The matrix need not be hermitian but needs to be
! nxn with n given in params.f90
SUBROUTINE diagonalize_g(L, eval, evec)
  ! Matrix to be diagonalized
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: L
  ! Right eigenvectors
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: evec
  ! Right eigenvalues
  COMPLEX*16, DIMENSION(n), INTENT(OUT) :: eval
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: work
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rwork
  COMPLEX*16, DIMENSION(n,1) :: temp_v
  COMPLEX*16 :: temp
  INTEGER :: i, j
  INTEGER :: info = 0
  ALLOCATE(work(n*10))
  ALLOCATE(rwork(n*2))
!  INTEGER :: LWORK = -1
  CALL zgeev('N','V',n,L,n,eval,evec,n,evec,n,work,n*10,rwork,info)
  IF (info .NE. 0) THEN
    WRITE(*,*) "Diagonalization error: ", info
  END IF
  temp = 0.0d0
  temp_v = 0.0d0
  DO j = 1, n
    DO i = 1, n - 1
      IF (ABS(eval(i)) > ABS(eval(i+1))) THEN
        temp_v(1:n,1) = evec(1:n,i)
        evec(1:n,i) = evec(1:n,i+1)
        evec(1:n,i+1) = temp_v(1:n,1)
        temp = eval(i)
        eval(i) = eval(i+1)
        eval(i+1) = temp
      END IF
    END DO
  END DO
  DEALLOCATE(work)
  DEALLOCATE(rwork)
END SUBROUTINE diagonalize_g

! Invert a matrix; return in L 
! Also not really needed because all the vectors I
! care about are orthonormal.
! Makes no assumption about Hermiticity.
! L is overwritten with the inverse matrix
SUBROUTINE invert(L)
  ! Matrix to be inverted
  COMPLEX*16, DIMENSION(n,n), INTENT(INOUT) :: L
  INTEGER, DIMENSION(n) :: pivots
  INTEGER:: info = 0
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: work
  ALLOCATE(work(n*10))
  ! LU decomposition
  CALL zgetrf(n,n,L,n,pivots,info)
  info = 0
  IF (info .NE. 0) THEN
    WRITE(*,*) "Inversion error: ", info
  END IF
  ! Inversion
  CALL zgetri(n,L,n,pivots,work,n*10,info)
  IF (info .NE. 0) THEN
    WRITE(*,*) "Inversion error: ", info
  END IF
  DEALLOCATE(work)
END SUBROUTINE

! Convert a given matrix to basis described by the
! Assumes that the output basis is dimension n and the
! input matix is dimension m
! eigenvector matrix evec and its inverse evec_inv
SUBROUTINE to_basis(matrix_in, matrix_out, evec, evec_inv)
  COMPLEX*16, DIMENSION(m,n), INTENT(IN) :: evec
  COMPLEX*16, DIMENSION(n,m), INTENT(IN) :: evec_inv
  COMPLEX*16, DIMENSION(m,m), INTENT(IN) :: matrix_in
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: matrix_out
  COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: temp
  ALLOCATE(temp(m,n))
  temp = MATMUL(matrix_in,evec)
  matrix_out = MATMUL(evec_inv,temp)
  DEALLOCATE(temp)
END SUBROUTINE

! Reverse the transformation given by to-basis
! Assumes that the input basis is dimension n and the 
! output matirx is dimension m
! Requires right eigenvector basis and their inverse
SUBROUTINE from_basis(matrix_in, matrix_out, evec, evec_inv)
  ! Eigenvectors and their inverse
  COMPLEX*16, DIMENSION(m,n), INTENT(IN) :: evec
  COMPLEX*16, DIMENSION(n,m), INTENT(IN) :: evec_inv
  ! Matrix to be transformed
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: matrix_in
  COMPLEX*16, DIMENSION(m,m), INTENT(OUT) :: matrix_out
  matrix_out = MATMUL(evec,MATMUL(matrix_in,evec_inv))
END SUBROUTINE

! Basis conversion when all pieces are of dimension m
! Requires right eigenvector basis and their inverse
SUBROUTINE to_basis_m(m, matrix, evec, evec_inv)
  INTEGER, INTENT(IN) :: m 
  ! Eigenvectors and their inverse
  COMPLEX*16, DIMENSION(m,m), INTENT(IN) :: evec
  COMPLEX*16, DIMENSION(m,m), INTENT(IN) :: evec_inv
  ! Matrix to be transformed
  COMPLEX*16, DIMENSION(m,m), INTENT(INOUT) :: matrix
  matrix = MATMUL(evec_inv,MATMUL(matrix,evec))
END SUBROUTINE

! Reverse the transformation given by to-basis
! Requires right eigenvector basis and their inverse
! Everything is of dimension m
SUBROUTINE from_basis_m(m, matrix, evec, evec_inv)
  ! Dimension of square matrix
  INTEGER, INTENT(IN) :: m
  ! Eigenvectors and their inverse
  COMPLEX*16, DIMENSION(m,m), INTENT(IN) :: evec, evec_inv
  ! Matrix to be transformed
  COMPLEX*16, DIMENSION(m,m), INTENT(INOUT) :: matrix
  matrix = MATMUL(evec,MATMUL(matrix,evec_inv))
END SUBROUTINE

! Get single element of + rate tensor gamma+(l,j,i,k)
! for Redfield tensor assembly
COMPLEX*16 FUNCTION gamma_plus_element(l,j,i,k,theta,Ga,Gb)
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: theta, Ga, Gb
  INTEGER, INTENT(IN) :: l, j, i, k
  gamma_plus_element = Ga(l,j)*Gb(i,k)*theta(i,k)/(hbar**2.0d0)
END FUNCTION

! Get single element of the - rate tensor gamma-(l,j,i,k)
! for Redfield tensor assembly
COMPLEX*16 FUNCTION gamma_minus_element(l,j,i,k,theta,Ga,Gb)
  COMPLEX*16, DIMENSION(n,n), INTENT(IN) :: theta, Ga, Gb
  INTEGER, INTENT(IN) :: l, j, i, k
  gamma_minus_element = Ga(l,j)*Gb(i,k)*theta(l,j)/(hbar**2.0d0)
END FUNCTION

! Set the index that lets the quadpack functions know
! which bath they are coupling to. Usually baths are
! identical and this is irrelevant.
SUBROUTINE set_cbi(i)
  INTEGER, INTENT(IN) :: i
  cbi = i
END SUBROUTINE

! Imaginary portion of the correlation function for quadpack
! Follows Nimitz quantum bath correlation function definition
! DO NOT CALL for 0 temperature; that is not this function
! and is no longer implemented
REAL(KIND=4) FUNCTION helper_corr_real(omegai)
  REAL(KIND=4), INTENT(IN) :: omegai
  ! omegai is input, but omega is actually used
  REAL(KIND=4) :: omega, t
  ! Helper terms for the calculation
  COMPLEX*16 :: n_of_omega!, e_pow, e_pow_n
  t = qp_t_g  ! get t from the global parameter
  ! We're approximating; get rid of potential singularity issues
  omega = omegai
  IF (omegai < 1.0d-14) THEN
    omega = 1.0d-14
  END IF

  n_of_omega = 1.0/(EXP(beta*omega*hbar) - 1.0)  
  helper_corr_real = spectral_density(1.0d0*omega,cbi)*(2.0*n_of_omega + 1.0)
END FUNCTION

! Real portion of the correlation function for quadpack
! Follows Nimitz quantum bath correlation function definition
! DO NOT CALL for 0 temperature; that is no longer implemented
REAL(KIND=4) FUNCTION helper_corr_imag(omegai)
  REAL(KIND=4), INTENT(IN) :: omegai
  ! omegai is input, but omega is actually used
  REAL(KIND=4) :: omega, t
  ! Helper terms for the calculation
  t = qp_t_g  ! get t from the global parameter
  ! We're approximating; get rid of potential singularity issues
  omega = omegai
  IF (omegai < 1.0d-14) THEN
    omega = 1.0d-14
  END IF

  helper_corr_imag = -spectral_density(1.0d0*omega,cbi)
END FUNCTION


! Spectral density function for bath i and frequency omega
REAL*8 FUNCTION spectral_density(omega, i)
  INTEGER, INTENT(IN) :: i
  REAL*8, INTENT(IN) :: omega
  ! Debye
  IF (bath_type(i) .EQ. 1) THEN
    spectral_density=eta(i)*omega/(omega**2.0d0 + omegac(i)**2.0d0) ! Debye
  ELSE 
  ! Ohmic
    spectral_density = eta(i)*omega*EXP(-ABS(omega)/omegac(i)) ! Ohmic
  END IF
END FUNCTION

! Integer factorial; use reals unless there's a good reason not to
INTEGER*8 FUNCTION factorial(x)
  INTEGER, INTENT(IN) :: x
  INTEGER i
  factorial = 1
  DO i=1,x
    factorial = i*factorial
  END DO
END FUNCTION

! Use recursion relations to fill a table of hermite polynomial coefficients
! for graphing the wavefunctions; stored in hf
SUBROUTINE fill_hermite(sz,hf)
  INTEGER, INTENT(IN) :: sz  ! How many polynomials to prepare
  REAL*8, DIMENSION(sz,sz), INTENT(OUT) :: hf  ! The polynomial coefficients
  INTEGER :: j, k 
  hf = 0
  hf(1,1) = 1.0d0
  hf(2,2) = 2.0d0
  DO j = 3, sz
    hf(j,1) = -hf(j-1,2)
    DO k = 2, j
      hf(j,k) = 2.0d0*hf(j-1,k-1)
      IF (k .LT. sz) THEN
        hf(j,k) = hf(j,k) - k*hf(j-1,k+1) 
      END IF
    END DO
  END DO
END SUBROUTINE

! Factorial in terms of real*8
REAL*16 FUNCTION fact_rn(n)
  INTEGER, INTENT(IN) :: n
  INTEGER :: i
  fact_rn = 1.0d0
  DO i = 1, n
    fact_rn = fact_rn*REAL(REAL(i))
  END DO
END FUNCTION

! Evaluate harmonic oscillator wavefunction density on a grid
! Fills grid of size szxres at points from 'low' to 'high'
! evenly spaced. There are sz functions to be evaluated
SUBROUTINE fill_wf_grid(sz,res,hf,grid,low,high)
  REAL*8, INTENT(IN) :: low, high  ! Start and end
  INTEGER, INTENT(IN) :: sz, res  ! size and resolution
  REAL*8, INTENT(OUT), DIMENSION(sz,res) :: grid  ! Grid of density
  REAL*8, INTENT(IN), DIMENSION(sz,sz) :: hf  ! Hermite factors
  INTEGER :: i, j
  REAL*8 :: gap, x

  gap = (high - low) / (res - 1)
   
  DO j = 1, sz
    DO i = 1, res
      x = (i-1)*gap + low
      grid(j,i) = eval_ho_wf(j,x,sz,hf) 
    END DO
  END DO
END SUBROUTINE


! For help determining the coupling and tuning basis
! Determine the nc and nt degrees for a location in the 
! wavefunction; i is location in the wavefunction
! x is the nt and y is the nc
! Should be in full, length m basis
SUBROUTINE get_loc(i,e,x,y)
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(OUT) :: e, x, y

  IF (i .GT. m/2) THEN
    e = 2
    y = (i - m/2 - 1)/dim_nt
    x = i - m/2 - y*dim_nt
    y = y + 1
  ELSE
    e = 1
    y = (i-1)/dim_nt 
    x = i - y*dim_nt
    y = y + 1
  END IF
END SUBROUTINE

! Make a plot that fills in the wavefunction probability
! density at points [low,high]*[low,high]
SUBROUTINE fill_wf_plot(sz,res,grid,low,high,wfn,plot)
  INTEGER, INTENT(IN) :: sz, res  ! Size of plot, resolution
  REAL*8, INTENT(IN) :: low, high  ! Bounds of plot
  COMPLEX*16, INTENT(IN), DIMENSION(m,1) :: wfn  ! The wavefunction in question (in unshifted HO basis)
  REAL*8, INTENT(OUT), DIMENSION(res,res,6) :: plot  ! The density plot
  REAL*8, INTENT(IN), DIMENSION(sz,res) :: grid  ! Filled in with 1d HO basis functions 
  COMPLEX*16, DIMENSION(res,res,2) :: plot_temp
  INTEGER :: i, j, k, e, nx, ny
  REAL*8 :: gap, x, y

  gap = (high - low) / (res - 1)
  plot_temp = 0.0d0
  ! Iterate over all wavefunction entries
  DO i = 1, m
    ! Determine the |e> X |Qc> X |Qt> composition
    CALL get_loc(i,e,nx,ny)
    ! Iterate across Qt coordinate
    DO j = 1, res
      x = gap*(j-1)
      ! Iterate across Qc coordinate
      DO k = 1, res
        y = gap*(k-1)
        ! Fetch values of HO wavefunctions of appropriate degree
        ! and that position on the grid. Fill in for proper diabatic electronic state.
        plot_temp(k,j,e) = plot_temp(k,j,e) + grid(nx,j)*grid(ny,k)*wfn(i,1)
      END DO
    END DO
  END DO

  DO j = 1, res
    DO k = 1, res
      DO i = 1, 2
        plot(j,k,i) = ABS(plot_temp(j,k,i))**2.0d0
      END DO
      plot(j,k,3) = REAL(REAL(plot_temp(j,k,1)))
      plot(j,k,4) = AIMAG(plot_temp(j,k,1))
      plot(j,k,5) = REAL(REAL(plot_temp(j,k,2)))
      plot(j,k,6) = AIMAG(plot_temp(j,k,2))
    END DO
  END DO
END SUBROUTINE

! Very simple 1-d grid integration to get the wavefunction
! density along the requested coordinate rather than as
! an overall 2d plot
SUBROUTINE integrate_out_dimension(coor,res,plot,low,high,out_grid)
  REAL*8, INTENT(IN), DIMENSION(res,res,6) :: plot  ! Wavefuction components in space 
  REAL*8, INTENT(OUT), DIMENSION(res,2) :: out_grid  ! Grid with integrated values
  INTEGER, INTENT(IN) :: res, coor  ! resolution of plot, coordinate 1 or 2 to integrate along
  REAL*8, INTENT(IN) :: low, high  ! Bounds of plot
  INTEGER :: i, j, k
  REAL*8 :: integrate_grid

  out_grid = 0.0d0
  DO i = 1, res
    DO k = 1, 2
      integrate_grid = 0.0d0
      DO j = 1, res
        IF (coor .EQ. 1) THEN
          integrate_grid = integrate_grid + plot(i,j,k)*(((high-low)/(1.0d0*res))**1.0d0)
        ELSE
          integrate_grid = integrate_grid + plot(j,i,k)*(((high-low)/(1.0d0*res))**1.0d0)
        END IF
      END DO
      out_grid(i,k) = integrate_grid
    END DO
  END DO

END SUBROUTINE

! Very simple grid integration to get an idea of how
! normalized/unnormalized a wavefunction plot is
! by summing along both dimensions. A good check
! of proper functionality.
REAL*8 FUNCTION integrate_grid(res,plot,low,high)
  REAL*8, INTENT(IN), DIMENSION(res,res,6) :: plot  ! Wavefunction components in space to integrate
  INTEGER, INTENT(IN) :: res   ! Resolution
  REAL*8, INTENT(IN) :: low, high  ! Low and high points of plot
  INTEGER :: i, j, k

  integrate_grid = 0.0d0
  DO i = 1, res
    DO j = 1, res
      DO k = 1, 2
        integrate_grid = integrate_grid + plot(i,j,k)*(((high-low)/(1.0d0*res))**2.0d0)
      END DO
    END DO
  END DO
END FUNCTION

! Evaluate the value of a harmonic oscillator wf at
! the given location; the variable input is assumed to
! be in dimensionless units SQRT(m*omega/hbar)
! One dimensional version
REAL*8 FUNCTION eval_ho_wf(deg,x,sz,hf)
  INTEGER, INTENT(IN) :: sz, deg  ! Size of grid, degree of polynomial
  REAL*8, INTENT(IN) :: x  ! Coordiante
  REAL*8, DIMENSION(sz,sz), INTENT(IN) :: hf  ! Hermite factors
  INTEGER :: i
  REAL*16 :: temp, temp2, tempx
  temp = 0.0d0
  tempx = x
  DO i = 1, deg
    temp = tempx**(i-1)*hf(deg,i) + temp
  END DO    
  temp2 = EXP(-(x**2.0d0)/2.0d0)*&
               (pi**(-0.25d0))/SQRT(fact_rn(deg - 1)*2.0d0**(deg-1))
  temp = temp*temp2
  eval_ho_wf = temp
END FUNCTION


! Plot grid prepared above; output the various
! values in "plot" which have wavefunction location
! information in a format that is convenient for
! gnuplotting as a heatmap
SUBROUTINE print_wf_plot(fnum,plot,low,high,res)
  INTEGER, INTENT(IN) :: fnum, res  ! File to write to, resolution of grid
  REAL*8, DIMENSION(res,res,6), INTENT(IN) :: plot  ! Wavefunction values
  REAL*8, INTENT(IN) :: low, high  ! Low and high extends of the grid
  REAL*8 :: gap  ! Gap between grid entries
  INTEGER :: i, j
  gap = (high-low)/(res - 1)
  DO i = 0, res - 1
    DO j = 0, res - 1
      WRITE(fnum,*) i*gap + low, j*gap + low, plot(i+1,j+1,1), plot(i+1,j+1,2), &
                    plot(i+1,j+1,3), plot(i+1,j+1,4), plot(i+1,j+1,5), plot(i+1,j+1,6)
    END DO
    WRITE(fnum,*) ""
  END DO
END SUBROUTINE

! Get average value of an operator on a wavefunction wfn; presumes same basis
! is used for operator X and wavefunction wfn
REAL*8 FUNCTION op_avg(wfn,X)
  COMPLEX*16, INTENT(IN), DIMENSION(n,n) :: X  ! Operator
  TYPE(wfunc), INTENT(IN) :: wfn  ! Wavefunction in
  COMPLEX*16, DIMENSION(n,1) :: wf
  REAL*8, DIMENSION(1,1) :: temp
  wf = wfn%fe
  temp = MATMUL(CONJG(TRANSPOSE(wf)),MATMUL(X,wf))
  op_avg = REAL(REAL(temp(1,1)))
END FUNCTION

! Takes an array of complex numbers and organizes
! them into an array of reals with twice the length
FUNCTION org_out(ent)
  COMPLEX*16, DIMENSION(n,1), INTENT(IN) :: ent
  REAL*8, DIMENSION(2*n) :: org_out
  INTEGER i
  org_out = 0.0d0
  DO i=1,n
    org_out(i*2 - 1) = REAL(REAL(ent(i,1)))
    org_out(i*2) = AIMAG(ent(i,1))
  END DO
END FUNCTION

! Will identify which eigenstate/noncollapsed
! state a wavefunction is in and set this
! information; 0 is not collapsed
INTEGER FUNCTION identify_estate(wfn)
  TYPE(wfunc), INTENT(IN) :: wfn
  INTEGER i, nz
  nz = 0
  DO i = 1, n
    IF (ABS(wfn%fe(i,1)) .GE. 1.0d-12) THEN
      IF (nz .EQ. 0) THEN
        nz = i 
      ELSE
        nz = 0
        EXIT
      END IF
    END IF 
  END DO
  identify_estate = nz
END FUNCTION

! Initialize an uncollapsed wavefunction
! This used to be important when dynamic 
! memory was involved
SUBROUTINE init_wfn(wfn)
  TYPE(wfunc), INTENT(INOUT) :: wfn
  wfn%collapsed = .FALSE. 
  wfn%sn = 0
END SUBROUTINE

! Destory a wavefunction (free memory)
! This used to be important when dynamic
! memory was involved
SUBROUTINE destroy_wfn(wfn)
  TYPE(wfunc), INTENT(INOUT) :: wfn
  wfn%collapsed = .FALSE.
END SUBROUTINE

! From a single wavefunction form a density matrix
SUBROUTINE wfn_to_rho_2(wfn,rho)
  COMPLEX*16, DIMENSION(n,n), INTENT(OUT) :: rho
  COMPLEX*16, DIMENSION(1,n), INTENT(IN) :: wfn
  INTEGER :: i, j
  DO i = 1, n
    DO j = 1, n
      rho(i,j) = wfn(1,j)*CONJG(wfn(1,i)) 
    END DO
  END DO
END SUBROUTINE

! Cauchy principal value integral one
! as used in the correlation function calculations
REAL*4 FUNCTION cpi_1(omega)
  REAL*4, INTENT(IN) :: omega
  cpi_1 = spectral_density(1.0d0*omega,cbi)*((EXP(beta*omega*hbar) - 1.0d0)**(-1.0d0))
END FUNCTION

! Cauchy principal value integral two
! as used in the correlation function calculations
REAL*4 FUNCTION cpi_2(omega)
  REAL*4, INTENT(IN) :: omega
  cpi_2 = spectral_density(1.0d0*omega,cbi)*((EXP(beta*omega*hbar)-1.0d0)**(-1.0d0) + 1.0d0)
END FUNCTION

! Tell the quadpack routine that integrates fourier transforms
! about what the extra required value that it needs to do the
! integration (omegaij) is by setting a global
SUBROUTINE set_omega_cbi(omegaij)
  REAL*8, INTENT(IN) :: omegaij
  omega_cbi = REAL(omegaij,4)
END SUBROUTINE

! Argument for the integral to infinity part of the CPV integral
! for the correlation function; first kind of weighting
REAL*4 FUNCTION cpi_1_weighted(omega)
  REAL*4, INTENT(IN) :: omega
  cpi_1_weighted = spectral_density(1.0d0*omega,cbi)*((EXP(beta*omega*hbar) - 1.0d0)**(-1.0d0))/&
          (omega - omega_cbi)
END FUNCTION

! Argument for the integral to infinity part of the CPV integral
! for the correlation function; second kind of weighting
REAL*4 FUNCTION cpi_2_weighted(omega)
  REAL*4, INTENT(IN) :: omega
  cpi_2_weighted = spectral_density(1.0d0*omega,cbi)*((EXP(beta*omega*hbar)-1.0d0)**(-1.0d0) + 1.0d0)/&
          (omega + omega_cbi)
END FUNCTION

! Read information from a trajectory file which details the progression
! of a wavefunction either by printing the full wavefunction or just
! its eigenstate and current phase.
! Read a double line record about a jump which has full wfn in it
! if necessary; returns 0 on success; -1 as estate is the flag for
! end of a trajectory
INTEGER FUNCTION read_dj_line(fnum,id,estate,time,wfn,jump_info)
  INTEGER, INTENT(IN) :: fnum  ! Read from
  INTEGER, INTENT(OUT) :: id, estate  ! unused value used to be important, wfn eigenstate
  REAL*8, INTENT(OUT) :: time  ! Jump time
  TYPE(wfunc), INTENT(OUT) :: wfn  ! Wavefunction read
  TYPE(jump_data), INTENT(OUT) :: jump_info  ! Informaiton about the jump
  INTEGER :: reason, i
  REAL*8, DIMENSION(2*n) :: temp
  REAL*8 :: re, im

  READ(fnum,*,IOSTAT=reason) id, time, estate
  jump_info%id = id
  jump_info%time = time
  jump_info%estate = estate
  
  IF (reason > 0) THEN
    WRITE(*,*) "Error ", reason
    read_dj_line = reason 
  ELSE IF (reason < 0) THEN
    read_dj_line = reason
  ELSE IF (estate .EQ. -1) THEN  ! Flag for no more jumps
    read_dj_line = 0
  ELSE
    IF (estate .EQ. 0) THEN  ! Read full wfn
      READ(fnum,*,IOSTAT=reason) (temp(i), i=1, 2*n)
      IF (reason .NE. 0) THEN  ! Read error
        WRITE(*,*) "Error ", reason
        read_dj_line = reason
      ELSE  ! Successful read
        DO i = 1, n
          wfn%fe(i,1) = CMPLX(temp(2*i -1),temp(2*i))
        END DO 
        jump_info%wfn = wfn
        read_dj_line = 0
      END IF
    ELSE  ! Read two pieces of eigendata rather than an array
      READ(fnum,*,IOSTAT=reason) re, im
      IF (reason .NE. 0) THEN  ! Read error
        read_dj_line = reason
        WRITE(*,*) "Error ", reason
      ELSE
        wfn%fe = 0.0d0
        wfn%fe(estate,1) = CMPLX(re,im)
        jump_info%wfn = wfn
        read_dj_line = 0
      END IF
    END IF
  END IF
END FUNCTION

! Give me an array of size x and I will compute the standard deviation
! of said array.
REAL*8 FUNCTION sdev(x,array)
  INTEGER, INTENT(IN) :: x
  REAL*8, DIMENSION(x), INTENT(IN) :: array
  REAL*8 :: average
  INTEGER :: i
  average = 0.0d0
  sdev = 0.0d0
  DO i = 1, x
    average = average + array(i)
  END DO
  average = average / x
  DO i = 1, x
    sdev = (array(i) - average)**2.0d0 
  END DO
  sdev = sdev / x
  sdev = SQRT(sdev)
END FUNCTION

END MODULE
