
INTRODUCTION
------------

MF-CG-TOID-MC is software to initialize, simulate, and analyze Monte Carlo
simulations of the Molecular Foundry Coarse-grained Model for Peptoids
(MF-CG-TOID).

LICENSE
-------

MF-CG-TOID-MC is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

MF-CG-TOID-MC is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
MF-CG-TOID-MC.  If not, see <http://www.gnu.org/licenses/>.

PUBLICATION
-----------

Details of the MF-CG-TOID model and simulation algorithms appear in our paper:

[1] T. K. Haxton, R. V. Mannige, R. N. Zuckermann, and S. Whitelam, "Modeling
sequence-specific polymers using anisotropic coarse-grained sites allows
quantitative comparison with experiment."

AUTHORS, CITATION, AND CONTACT INFORMATION
-------------------------------------------

MF-CG-TOID-MC was developed by Thomas Haxton at the Molecular Foundry, Lawrence 
Berkeley National Laboratory.

We ask that users please cite our publication [1] in any publication presenting 
results produced with MF-CG-TOID.

The latest version of MF-CG-TOID-MC can be found at 
nanotheory.lbl.gov/peptoid_code.

Contacts: 
Thomas Haxton (tomhaxton@gmail.com)
Stephen Whitelam (swhitelam@lbl.gov) 

PLATFORM
--------

MF-CG-TOID-MC has been tested on Linux and OS X.

REQUIREMENTS
------------

1. C compiler (e.g. gcc) to compile source code
2. VMD to view movies of molecular trajectories
3. Plotting software (e.g. gnuplot) to plot output data

INSTALLATION
------------

Run "make" in the directory "code".

CONTENTS/USAGE
--------------

monolayer

   Performs a Monte Carlo search over 44 degrees of freedom to find a low-energy
   monolayer configuration

   Selected options:
      -inputfile (optional): Monolayer configuration file (if no inputfile is
         given, the code first solves for a low-energy configuration of a single
         peptoid chain)
      -base: Directory for output files
      -chemistry: Series of integers denoting the sequence(s) of the peptoids.
         For a monodisperse collection of N block-4n peptoids, the syntax is:
            -chemistry 1 1 N 2 1 2 n 3 n
      -interchainspacing: Initial spacing between parallel chains in angstroms
      -temperature: Fictitious temperature in K used in the Monte Carlo search
      -cycles: Number of Monte Carlo cycles

   Inputs:
      <inputfile> (optional)

   Outputs:
      <base>/c<cycles>.cg.xyz: Time series of coarse-grained coordinates in XYZ
         format, recorded each time the Monte Carlo search finds a new lowest
         energy configuration.  Load in VMD Tk Console with 
         "mol new c<cycles>.cg.xyz".
      <base>/c<cycles>.cg.vmd.source: VMD Tcl script for loading bonds, colors,
         and radii of coarse-grained coordinates.  After loading XYZ file, run
         Tcl script with "source c<cycles>.cg.source".
      <base>/c<cycles>.allatom.xyz: Time series of all-atom coordinates
      <base>/c<cycles>.allatom.source: VMD Tcl script for all-atom coordinates
      <base>/c<cycles>.best.energy: Time series of total energy and components
         of the energy, recorded each time the search finds a new lowest energy
      <base>/c<cycles>.best.config: Monolayer configuration file (time series of
         the 44 degrees of freedom, recorded each time the search finds a new
         lowest energy)

   Example script: monolayer_search.sh (32 sec on 2.9 GHz processor)

bilayer_from_monolayer_enumerate

   Calculates energies of bilayers formed by sandwiching two monolayers
   together, iterating over a rectilinear grid of separations

   Selected options:
      -inputfile: Monolayer configuration file
      -outputfile: Output file for printing separations and energies
      -chemistry: Series of integers denoting the sequence(s) of the peptoids.
         For a monodisperse collection of N block-4n peptoids, the syntax is:
            -chemistry 1 1 N 2 1 2 n 3 n
         Note that N must account for both leaves of the monolayer, but the
         system size does not have to match the system size of <inputfile>
      -firstx -lastx -xspace -firsty -lasty -yspace -firstz -lastz -zspace:
         Parameters defining the rectilinear grid.  x parameters are in units of
         the residue spacing, y parameters are in units of the polymer spacing,
         and z parameters are in angstroms.

   Inputs:
      <inputfile>

   Outputs:
      <outputfile>

   Example script: bilayer_sandwich.sh (15 sec on 2.9 GHz processor)

bilayer_from_monolayer

   Performs a Monte Carlo search over 84 degrees of freedom to find a low-energy
   bilayer configuration, starting from a monolayer configuration (e.g. one
   minimized by "monolayer") and given a starting leaf separation (e.g. found
   with "bilayer_from_monolayer_enumerate")

   Selected options:
      -inputfile: Monolayer configuration file
      -base: Directory for output files
      -chemistry: Series of integers denoting the sequence(s) of the peptoids.
         For a monodisperse collection of N block-4n peptoids, the syntax is:
            -chemistry 1 1 N 2 1 2 n 3 n
      -leafxoffsetfrac -leafyoffsetfrac -leafspacing: Starting leaf separation,
         in units of the residuespacing (x), polymer spacing (7), and angstroms
      -temperature: Fictitious temperature in K used in the Monte Carlo search
      -cycles: Number of Monte Carlo cycles

   Inputs:
      <inputfile>

   Outputs:
      <base>/c<cycles>.cg.xyz: Time series of coarse-grained coordinates in XYZ
         format, recorded each time the Monte Carlo search finds a new lowest
         energy configuration.  Load in VMD Tk Console with 
         "mol new c<cycles>.cg.xyz".
      <base>/c<cycles>.cg.vmd.source: VMD Tcl script for loading bonds, colors,
         and radii of coarse-grained coordinates.  After loading XYZ file, run
         Tcl script with "source c<cycles>.cg.source".
      <base>/c<cycles>.allatom.xyz: Time series of all-atom coordinates
      <base>/c<cycles>.allatom.source: VMD Tcl script for all-atom coordinates
      <base>/c<cycles>.best.energy: Time series of total energy and components
         of the energy, recorded each time the search finds a new lowest energy
      <base>/c<cycles>.best.config: Bilayer configuration file (time series of
         the 84 degrees of freedom, recorded each time the search finds a new
         lowest energy)

   Example script: bilayer_search.sh (1 min 3 sec on 2.9 GHz processor)

peptoid

   Performs Monte Carlo simulations of the MF-CG-TOID model

   Selected options:
      -interface: Denotes whether (1) or not (0) the simulation has an air-water
         interface
      -ictype: Initial condition type
      -interface, -ictype:
         0, 0: Simulation coordinate file (continuing a simulation without an
            interface)
         0, 1: Dilute 2d lattice (or single chain) in low-energy configuration
         0, 3: Bilayer configuration from file
         0, 4: Dilute, random, non-overlapping solution (or single chain)
            replicated from single-chain conformation from file
         0, 5: Free-floating (non-periodic) bilayer configuration from file
         0, 6: Bilayer made from low-energy single-chain conformations
         0, 7: Stack of bilayers made from replicating a single bilayer
         1, 0: Simulation coordinate file (continuing a simulation with an
            interface)
         1, 2: Monolayer configuration from file
         1, 3: Bilayer configuration from file
         1, 4: Prescribed monolayer configuration
         1, 6: Dilute 2d lattice (or single chain) in low-energy configuration
         1, 7: Single chain from file
      -inputfile: Simulation coordinate file, monolayer configuration file,
         bilayer configuration file, or N/A, depending on <interface> and
         <ictype>
      -base: Directory for output files
      -chemistry: Series of integers denoting the sequence(s) of the peptoids.
         For a monodisperse collection of N block-4n peptoids, the syntax is:
            -chemistry 1 1 N 2 1 2 n 3 n
      -xchains: Number of chains in the x (chain backbone) direction (for 
         monolayers and bilayers initialized from configuration files or low-
         energy configurations)
      -temperature: Temperature in K (default 300)
      -surfacepressure: Surface pressure in kcal/mol/angstrom^2 (default 0)
      -density: Concentration in residues/angstrom^3 (used only for 2d lattice
         initial conditions)
      -molarity: Concentration in moles/L (used only for dilute, random solution
         and free-floating bilayers)
      -wholebilayertranslatefreq: Frequency per chain to translate entire 
         bilayer (used for stack of bilayers)
      -shiftbilayergapfreq: Frequency per chain to attempt "evaporation moves,"
         changing the gap between bilayers (used for stack of bilayers)
      -stackx, -stacky, -stackz: Separation between stacked bilayers, in 
         angstroms 
      -reset: Denotes whether (1) or not (0; default) to reset time stamp
      -cycles: Number of Monte Carlo cycles
      -maxhours: Maximum numer of hours before exiting the simulation

   Inputs:
      <inputfile>: Simulation coordinate file, monolayer configuration file,
         bilayer configuration file, or N/A, depending on <interface> and 
         <ictype>

   Outputs:
      <base>/vmd.cg.xyz: Time series of coarse-grained coordinates in XYZ
         format.  Load in VMD Tk Console with "mol new c<cycles>.cg.xyz".
      <base>/vmd.cg.source: VMD Tcl script for loading bonds, colors, and radii 
         of coarse-grained coordinates.  After loading XYZ file, run Tcl script 
         with "source c<cycles>.cg.source".
      <base>/vmd.allatom.xyz: Time series of all-atom coordinates
      <base>/vmd.allatom.source: VMD Tcl script for all-atom coordinates
      <base>/timeseries: Time series of energy components and simulation box
         dimensions
      <base>/trajectory: Time series of coarse-grained site coordinates
      <base>/com: Time series of center-of-mass position (when simulation
         contains only one chain)
      <base>/coord: Simulation coordinate file containing coordinates at the end
         of the simulation

   Example scripts:
      single_chain.sh (1 min 39 sec on 2.9 GHz processor)
      single_chain_interface.sh (2 min 6 sec on 2.9 GHz processor)
      monolayer.sh (2 hr 19 min on 2.9 GHz processor)
      bilayer.sh (4 hr 30 min on 2.9 GHz processor)
      bilayer_stack.sh (54 min on 2.9 GHz processor)
      bilayer_stack_evaporate.sh (4 hr 48 min on 2.9 GHz processor)

post_analysis

   Calculates X-ray and neutron scattering from estimated all-atom positions and
   pair distribution functions between coarse-grained sites

   Selected options:
      -chainlength: Number of residues in peptoid chains (assumes a monodisperse
         system)
      -Nchains: Number of chains
      -mincycle: Minimum Monte Carlo cycle used for averaging
      -nframes: Maximum number of frames used for averaging (Simulations output
         1 frame per 100000 cycles by default)
      -trajectoryfile: Input time series of coarse-grained site coordinates
      -xrdcode: Code denoting whether scattering calculation should use
         periodicity in x and y (-xrdcode 2; appropriate for monolayers and
         free-floating bilayers) or x, y, and z (-xrdcode 3; appropriate for
         stacks of bilayers)
      -directory: Directory for output files

   Inputs:
      <trajectoryfile>

   Outputs:
      <directory>/gr.<firstframe>-<lastframe>.oppositeleaf.2dnorm
      <directory>/gr.<firstframe>-<lastframe>.oppositeleaf.3dnorm
      <directory>/gr.<firstframe>-<lastframe>.sameleaf.2dnorm
      <directory>/gr.<firstframe>-<lastframe>.sameleaf.3dnorm
      <directory>/gr.<firstframe>-<lastframe>.samepoly.2dnorm
      <directory>/gr.<firstframe>-<lastframe>.samepoly.3dnorm
         Radial distribution functions between coarse-grained sites, separated
         by sites on the same polymer, different polymer but same leaf, and
         different leaves, and normalized to compare to a uniform distribution
         in two or three dimensions
      <directory>/histo.<firstframe>-<lastframe>.normr
      <directory>/histo.<firstframe>-<lastframe>.normr.rleftdotn
      <directory>/histo.<firstframe>-<lastframe>.normr.rrightdotn
      <directory>/histo.<firstframe>-<lastframe>.rleftdotn
      <directory>/histo.<firstframe>-<lastframe>.rrightdotn
      <directory>/histo.<firstframe>-<lastframe>.tripleproduct
      <directory>/histo.<firstframe>-<lastframe>.tripleproduct.rleftdotn
      <directory>/histo.<firstframe>-<lastframe>.tripleproduct.rrightdotn
         One- and two-dimensional histograms of the scalar coordinates appearing
         in the backbone bonded interaction
      <directory>/phenhisto.<firstframe>-<lastframe>
      <directory>/phenmixedhisto.<firstframe>-<lastframe>
      <directory>/phenmixedsubtracthisto.<firstframe>-<lastframe>
         Histograms for the scalar coordinates appearing in the Gay-Berne
         potential for the phenylethyl-phenylethyl nonbonded interaction
      <directory>/sidehisto.<firstframe>-<lastframe>.ndotn
      <directory>/sidehisto.<firstframe>-<lastframe>.rdotn
      <directory>/sidehisto.<firstframe>-<lastframe>.rparallel
      <directory>/sidehisto.<firstframe>-<lastframe>.rperp
         Histograms for the scalara coordinates appearing in the sidechain
         bonded interaction
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.raw
         Two-dimensional scattering spectrum in the x-y plane, broken down by 
         atom type and before multiplying by an atomic scattering factor
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.weighted
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.neutron.weighted
         Two-dimensional (x-y) scattering spectra by atom type after multiplying  
         by the appropriate scattering factors for X-ray and neutron scattering
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.total
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.neutron.total
         Total x-y X-ray and neutron scattering spectra
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.radav.raw
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.radav.weighted
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.radav.neutron.
         weighted
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.total
      <directory/xrd.allatom.<firstframe>-<lastframe>.faceon.neutron.total
         Radial averages (one dimensional spectra) in the x-y plane
      <directory/xrd.allatom.<firstframe>-<lastframe>.transverse.raw
      <directory/xrd.allatom.<firstframe>-<lastframe>.transverse.weighted
      <directory/xrd.allatom.<firstframe>-<lastframe>.transverse.neutron.
         weighted
      <directory/xrd.allatom.<firstframe>-<lastframe>.transverse.total
      <directory/xrd.allatom.<firstframe>-<lastframe>.transverse.neutron.total
         Transverse (z) scattering spectra

   Example scripts:
      post_analysis_monolayer.sh (3 min 6 sec on 2.9 GHz processor)
      post_analysis_bilayer.sh (6 min 18 sec on 2.9 GHz processor)
      post_analysis_stack.sh (1 min 53 sec on 2.9 GHz processor)

EXAMPLE SCRIPTS
---------------

Example scripts are in the directory "example_scripts."  Output generated from
these scripts are in the directory "output."

EXAMPLE GNUPLOT SCRIPTS
-----------------------

Example scripts for plotting data from the example output files are in the
directory "example_gnuplot_scripts."  To run them, use "load <script>" in
gnuplot.