The following contains code and data to accompany 'Nonadiabatic transition paths
from quantum jump trajectories' Michelle C. Anderson, Addison J. Schile, and 
David T. Limmer.

LICENSE

This repository is distributed under the GPL license (see LICENSE).

The included TPT_CI is free software: you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation, either version 3 of the License, 
or (at your option) any later version.
TPT_CI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with TPT_CI. 
If not, see <https://www.gnu.org/licenses/>.

ZENODO 

This repository is archived permanently at: 10.5281/zenodo.6950370

PAPER DATA

Data and gnuplot scripts necessary to regenerate the figures in the paper are found
in folders titled 'figure_x'. Regeneration requires gnuplot. 
To generate the plots in file 'figure_x.ps' or both 'figure_x_right.ps' and
'figure_x_left.ps', run
'gnuplot plot_#.gs'

Note that in the code well 2 refers to the higher energy well, which is diabatic
state 1 in the paper writeup.

--------------------------------------------------------------------------------------------

PROGRAMS

All files are written in FORTRAN 90/95 and are found in the folder 'code.'

The programs necessary to carry out secular Redfield calculations, 
Lindblad jump calculations, Markov state model analysis, and plotting
functions are main, create_msm, markov_analysis, path_helper, sort, and
wfn_plot. The file params.f90 contains specifications to control the
programs, with a few controls being present in other program files, most
of these referring to very specific functions such as whether to impose
detailed balance when calculating committors in markov_analysis.
Note that detailed balance imposition is done after calculating both
forward and backward committors numerically and choosing, for each state in 
the Markov state model, the smaller of the two, be it forward or reverse, 
then using this smaller value to calculate its complement. This avoids 
catastrophic problems which result when important, small numbers
are overwritten with zero to enforce detailed balance.

The programs can be compiled from the given Makefile with

'make'

Most programs require Lapack and BLAS. Quadpack is also required
but the necessary file is included and compiled by make. Quadpack
is licensed as specified in 'quadpack.f90' under LGPL.

--------------------------------------------------------------------------------------------

Keep in mind that most programs are stochastic and or depend on solutions to 
systems of equations which may be poorly conditioned and that the path
selection algorithm for transition pathways can be noisy. In other words,
if the dominant pathways for a barrier crossing in one direction do not
match the dominant pathways for the same barrier crossing in the opposite
direction, this is likely a sign of a serious problem but if minor
pathways do not match exactly this is likely due to numerical differences
when solving the for the committors. Imposing detailed balance on a system,
even if that system already follows detailed balance rules by definition,
changes the rule by which the reverse committors are calculated and will have
an effect on the results. If this effect is not minor or confined to
less relevant reaction pathways, this may be indicative of a problem.

--------------------------------------------------------------------------------------------

To run a full analysis, including Markov state model and vertical 
excitation analysis, the following commands would need to be run
using the 'params.f90' file with the included specifications:

make  
./create_msm  # Generate a Markov state model

./main  # Generate Lindblad jump trajectories using default random seed

./sort  # Reads in data from 'main' and sorts according to where the trajectories ended into ending[12].txt

mkdir ending1

mkdir ending2

mv ending1.txt ending1/ph_jumps.txt  # Path helper reads from ph_jumps.txt

mv ending2.txt ending2/ph_jumps.txt

cp path_helper ending1/  # Organize according to all jumps witnessed in the trajectories

cp path_helper ending2/  

cp markov_mat.txt ending1/ # The equilibrium Markov state model is necessary for vertical analysis, too 

cp markov_mat.txt ending2/

cp evic.txt ending1/  # Contains information about the eigenvectors to avoid reforming the Hamiltonian

cp evic.txt ending2/

cd ending1

./path_helper  # Reads in data from 'sort' moved to ph_jumps.txt and catalogues what jumps occur

cd ../ending2

./path_helper  # Reads in data from 'sort' moved to ph_jumps.txt and catalogues what jumps occur

./markov_analysis  # Reads in data from evic.txt, markov_mat.txt and jumps_markov_mat.txt to find reactive paths

cd ../ending1

./markov_analysis  


The files of principal interest are then:

1) either ending's 'thermal_paths_1toB.txt' which contains the reaction pathways
between eigenstate 1 and eigenstate B as specified in 'markov_analysis_utils.f90' as well
as information about the average pre and post committor states and fluxes along reaction paths
2) 'ending1/graph_bt_to1.txt' which contains the dominant vertical relaxation
pathways bound for eigenstate 1 and average information about the post committor jumps and fluxes
3) 'ending2/graph_bt_toB.txt' which contains the dominant vertical relaxation
pathways bound for eigenstate B and average information about the post committor jumps and fluxes
4) either ending's 'thermal_comm_1toB.txt' which contains the thermal forwards and backwards 
committors when the source is eigenstate 1 and destination is eigenstate B


--------------------------------------------------------------------------------------------

Programs in General:

Programs were constructed and tested on Mac OS Mojave and compiled with gfortran 8.2.0. The
makefile settings assume that gfortran, Lapack and BLAS are available.
To change the compiler or library linking settings, edit 'Makefile'.

Note, all programs that call a setup function will write a file 'evic.txt'
which includes information about the average characteristics of energy eigenvectors of in the format:

'(coupling coordinate) (tuning coordinate) (energy) (diabatic population 1) (diabatic population 2)'.

A more human readable version is in 'eigenvector_info.txt' with columns labeled.
All programs also produce a file, 'plot_ev.txt' which merely contains:

'(eigenvecor number) (-100 or 100) (eigenvector energy)' 

which can be used to plot the energies of the eigenvectors.

It is recommended to use atomic units for all calculations. Although the value of reduced
Planck's constant can be changed in the parameters file, this option was not thoroughly 
tested.

--------------------------------------------------------------------------------------------

main [seed]

seed is an optional three digit number to seed the random number generator.
It must be precisely three digits and the default will be used if it is not
present or not precisely three digits.

Depending on parameters specified in 'params.f90', main  will run either
a secular Redfield calculation for a fixed period of time starting either
from a vertical excitation or an eigenstate, or a Lindblad jumps calculation
from either a vertical excitation or eigenstate which continues until 
the propagating wavefunction collapses to one of the two flag eigenstates
specifed in 'params.f90.' This can go on indefinitely if the flag eigenstates
are specified incorrectly.
Parameters to specify the number of baths as well as their cut off frequencies,
coupling strengths and temperatures are available. The basis size 'm' must
always be (coupling basis size)*(tuning basis size)*2 but 'n', the truncated
basis size, can be altered to be any even number less than or equal to m.

In the case of a secular run, diabatic populations will be output to
'evolution.txt' with the columns indicating '(time) (average diabatic electronic population
1) (average diabatic electronic population 2) (average coupling coordinate) (average
tuning coordinate). The file 'state_prop.txt' will include (iteration) (time) and
the populations of the first seventeen eigenstates if 'n' in params is seventeen or higher.
Note that in the accompanying paper, the labels of the diabatic states were
reserves relative to the code implementation.

In the case of a Lindblad jumps run, 'jump_data.txt' will include the complete
record of all wavefunction propogations, starting from the initial wavefunction
followed by a print out of the wavefunction at every time at which a jump occurs. The
wavefunction will be printed with the following format: 

'0 (time) 0

REAL(coefficients(1)) IMAGINARY(coefficients(1)) REAL(coefficients(2)) IMAGINARY(coefficients(2))...'

if the wavefunction has not collapsed to an eigenstate. The coefficients refer to the
contribution of each energy eigenstate to the overall wavefunction
If the wavefunction has collapsed to an eigenstate the output will read:

'0 (time) (collapsed eigenstate)
REAL(coefficients(collapsed eigenstate)), IMAGINARY(coefficients(collapsed eigenstate))'

meaning only data for the sole occupied entry will be printed to save space. The leading
zero in the specification line, 'id' in the code, was important when OMP functions were
employed but is no longer relevant and always set to zero.

--------------------------------------------------------------------------------------------

create_msm

Prepares a Markov State Model of a quantum system from Secular Redfield calculations.
Runs a short secular calculation starting from each eigenstate in turn for the 
specified system and assembles the populations for each run into a Markov state model
printed to 'markov_mat.txt' with the following format:

'(n=matrix dimension) (t=duration of propogation)

Equilibrium populations of all energy eigenstates in the system

Populations of eigenstates at t following initialization in energy eigenstate 1

Populations of eigenstates at t following initialization in energy eigenstate 2

...

Populations of eigenstates at t following initialization in energy eigenstate n'

--------------------------------------------------------------------------------------------

markov_analysis

Deals with analysis of a Markov State Model in which reactive pathways from eigenstate
1 to an eigenstate specified by B are of interest.  Eigestate B is specified by a parameter 
'B_estate' in 'markov_analysis_utils.f90'.
Reads in 'markov_mat.txt' and, if available, 'jumps_markov_mat.txt'. The first is used to 
determine equilibrium committors. Committors from reactant (energy eigenstate 1) to product
(energy eigenstate labeled by B) will be written in 'thermal_comm_Bto1.txt in the format:

'(eigenstate) (committor to 1) (committor to B)'.

The format for 'thermal_comm_1toB.txt' is the same with the destination and origin reversed. 
The reactive pathway information is in 'thermal_paths_1toB.txt' and 'thermal_paths_1toB_plot.txt' 
with the latter being designed for plotting and the former being human readable and 
containing extra information.


'thermal_paths_1toB.txt' prints the total reactive flux then the reactive paths in order of
largest flux. The format for this is:

'(pathway flux) (fraction total flux on this pathway) (energy eigenstates on pathway in order)'

In the next section, the same information is then printed again, but this time the eigenstates 
are replaced by the energy of the eigenstates for each entry on the path.

After pathway information, the total flux is again printed followed by specific information
about the jump at which each given pathway exceeds a committor probability of 50%, i.e. each
committor jump is printed along with the total flux that passes through that jump and the
percentage of total flux that passes through that jump.

The average over all paths of energy, coupling coordinate, tuning coordinate, diabatic state 1 character,
diabatic state 2 character, and forward committor probability for the jump just before
committment is then displayed along with the same data for the postcomittor jump and the difference
observed over the committor jump.

The last section of the file is a running sum of the total and total percent reactive flux accounted
for as a function of reactive pathways. The pathways are in order of decreasing flux.

The file 'flux.txt' shows the reactive flux and rate for thermal passage from reactants to products
and from products to reactants.

--------------------------------------------------------------------------------------------

path_helper

Prepares a graph where vertices are eigenstates and edges are the number of jumps observed
between vertices over an ensemble of Lindblad jump trajectories.
Reads through the data in 'ph_jumps.txt', which should be the renamed 'ending2.txt' or
'ending1.txt.' It records every observed jump to determine how many jumps occur between
any given energy eigenstate. The jump data will be written in the same format as a Markov
matrix in 'jumps_markov_mat.txt':

'(dimension of matrix)

(Jumps from 1-->1 are not recorded) (Jumps from 1-->2) (Jumps from 1-->3) .... (Jumps from 1-->uncollapsed are impossible)

(Jumps from 2-->1) (Jumps from 2-->2 are not recorded)... (Jumps from 2-->uncollapsed are impossible)

....

(Jumps from uncollapsed-->1) (Jumps from uncollaped-->2)... (Jumps from uncollapsed to uncollapsed not recorded)'


Note that because the uncollapsed wavefuntion has a state in this matrix, the matrix will be one 
entry larger than the equilibrium Markov state model matrix.

Other files printed are 'jumps_markov_mat_s#.txt' which contains a subset of the jump data for block 
averaging. Five of these are written, each containing a fifth of the data.

The file 'end_time.txt' and 'collapse_time.txt' include information about when trajectories finish and 
collapse in the format:
'(limit time) (number trajectories collapsed/completed between previous limit time and this one) (cummulative number collapsed/complete)'

The file 'num_jumps_histogram.txt' includes data to plot a histogram of the number of jumps between distinct
eigenstates. The format is
'(start of bin) (end of bin) (number of trajectories with total jump number between the bin bounds) 
(number of trajectories with total jumps lower than end of bin, inclusive)'

--------------------------------------------------------------------------------------------

sort

Splits the ensemble in two based on the final eigenstate.
Reads data in 'jump_data.txt' and sorts all of the wavefunction propogations recorded there
by whether they end in flag state 1 or flag state 2. Writes the results to 'ending1.txt' and
'ending2.txt' respectively in the exact same format as 'jump_data.txt'

--------------------------------------------------------------------------------------------

wfn_plot

This does not need to read any information. Plots all wavefunctions specified in 'params.f90' between bound1 and
 bound2 and between bound3 and bound4. The wavefunctions are written to files 'CI_data/wfn_plots/ef_2d_i.txt' and 
'CI_data/wfn_plots/c2_1d_i.txt' and 'CI_data/wfn_plots/c1_1d_i.txt'. The containing folders have to be made before the
program will run.

'CI_data/wfn_plots/ef_2d_i.txt' has the format:

'(x coordinate) (y coordinate) (wfn density in diabatic state 1) (wfn density in diabatic state 2) 
(REAL(wfn value in diabatic state 1)) (IMAGINARY(wavefunction value in diabatic state 2)) 
(REAL(wavefunction value in diabatic state 1)) (IMAGINARY(wavefunction value in diabatic state 2))'

'CI_data/wfn_plots/cx_1d_i.txt' has integrated out one dimension or the other to show a 1d wavefunction with format:
(index) (coordinate) (wavefunction density in diabatic state 1) (wavefunction density in diabatic sate 2)

Optionally, wfn_plot will attempt to open the file 'in_wfns.txt' and read in 2*n real numbers
from the file which are interpreted to be the real and imaginary parts of every complex
number required to specify a full truncated basis wavefunction. It will then plot this
wavefunction. Note that wfn_plot will not normalize this wavefunction. The basis must be
the trunacted energy eigenbasis.
