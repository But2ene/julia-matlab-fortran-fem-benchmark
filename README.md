Julia, Matlab, Octave, and Fortran FEM Benchmark and Comparison
===============================================================

Benchmark and comparison of Julia, Matlab, Octave, and Fortran for a
2D Poisson problem solved on a unit square. The problem is discretized
with Q1 bilinear Lagrange finite elements.

The Octave and Matlab code is derived from the production
[FEATool Multiphysics](https://www.featool.com/) code. The Fortran 77
code is using the Feat2D FEM library found in the FEM CFD Solver
[FeatFlow](http://www.featflow.de) from which the Julia code is a
direct port. The problem setup is identical and the codes really are
equivalent and not _dummied down_ to perform better in the benchmark,
as they still feature all relevant code paths to allow different
finite element shape functions, variable coefficients, unstructured
grids, etc. Thus the results are relevant and comparable.


Installation and Running
------------------------

- Download and unzip or clone the repository

- Set up and install Matlab,
  [Octave](https://www.gnu.org/software/octave/),
  [Julia](http://julialang.org/), and a Fortran compiler (tested with
  gfortran and the Intel Fortran compiler)

- Edit the *testrun_param.txt* file, which contains three parameters

        N0    - the grid resolution of the coarsest test grid (grid level)
                (the number of cells in x and y-directions)
        NLEV  - The number of grid levels to run tests for
        NRUNS - Number of test runs for each grid level,
                the timings are averaged for all runs

- On Windows edit the _OCTAVE_, _MATLAB_, and _JULIA_ paths in the
  *run_tests.bat* script and execute the script to automatically
  run the benchmarks and generate the output files

- On other systems the *run_matlab.m*, *run_julia.jl*, and
  *run_fortran.sh* shell scripts can be run manually. Afterwards the
  postprocessing script *src_matlab/process_results.m* file can be run
  to generate the results and output files.

- The _NNWORK_ parameter in the main Fortran source file
  *src_fortran/src/featfem.f* controls the static memory allocation
  and might have to be increased and recompiled to run > 4 GB runs


Results and Output
------------------

After running the tests the _output_ directory will contain results
files, _table.txt_ with tabulated times, and _results.html_ which
shows [Plotly](https://plot.ly/) comparison graphs. Sample output is
shown here below

![Assembly Timings](https://raw.githubusercontent.com/precisesimulation/julia-matlab-fortran-fem-benchmark/tree/master/output/fig_assembly.jpg)

![Solver Timings](https://raw.githubusercontent.com/precisesimulation/julia-matlab-fortran-fem-benchmark/tree/master/output/fig_solve.jpg)

![Total Timings](https://raw.githubusercontent.com/precisesimulation/julia-matlab-fortran-fem-benchmark/tree/master/output/fig_total.jpg)


    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |  Octave |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |   t_solve |     t_tot |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 2.60e-003 | 2.60e-003 | 1.65e-002 | 1.30e-002 | 2.61e-003 | 0.00e+000 | 4.34e-003 | 4.17e-002 |
    |      64 | 8.67e-004 | 8.70e-004 | 1.91e-002 | 1.22e-002 | 5.21e-003 | 3.47e-003 | 2.20e-002 | 6.37e-002 |
    |     128 | 7.82e-003 | 5.21e-003 | 3.82e-002 | 2.60e-002 | 1.22e-002 | 6.08e-003 | 9.33e-002 | 1.89e-001 |
    |     256 | 1.65e-002 | 2.43e-002 | 1.35e-001 | 9.64e-002 | 3.99e-002 | 3.21e-002 | 5.14e-001 | 8.57e-001 |
    |     512 | 7.64e-002 | 8.94e-002 | 4.92e-001 | 7.28e-001 | 1.58e-001 | 1.38e-001 | 2.74e+000 | 4.42e+000 |
    |    1024 | 2.96e-001 | 3.43e-001 | 1.93e+000 | 7.33e+000 | 6.70e-001 | 5.39e-001 | 1.82e+001 | 2.93e+001 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|

    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |  Matlab |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |   t_solve |     t_tot |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 1.88e-003 | 1.87e-003 | 2.81e-003 | 3.26e-003 | 2.66e-003 | 2.12e-003 | 4.19e-003 | 1.88e-002 |
    |      64 | 7.05e-004 | 1.09e-003 | 6.61e-003 | 4.15e-003 | 5.95e-003 | 9.56e-003 | 2.17e-002 | 4.97e-002 |
    |     128 | 1.44e-003 | 6.26e-003 | 1.62e-002 | 1.24e-002 | 2.28e-002 | 4.07e-002 | 9.02e-002 | 1.90e-001 |
    |     256 | 5.21e-003 | 2.53e-002 | 5.54e-002 | 5.06e-002 | 8.99e-002 | 1.47e-001 | 4.46e-001 | 8.20e-001 |
    |     512 | 3.05e-002 | 9.79e-002 | 2.03e-001 | 3.39e-001 | 3.67e-001 | 6.10e-001 | 2.53e+000 | 4.18e+000 |
    |    1024 | 1.24e-001 | 3.95e-001 | 8.15e-001 | 3.30e+000 | 1.51e+000 | 2.52e+000 | 1.49e+001 | 2.36e+001 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|

    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |   Julia |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |   t_solve |     t_tot |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 1.67e-004 | 5.81e-004 | 8.92e-004 | 4.27e-004 | 4.90e-005 | 2.95e-004 | 7.65e-003 | 1.01e-002 |
    |      64 | 5.45e-004 | 2.16e-003 | 3.56e-003 | 1.86e-003 | 7.68e-005 | 1.48e-003 | 3.03e-002 | 4.00e-002 |
    |     128 | 1.80e-003 | 7.35e-003 | 1.43e-002 | 6.77e-003 | 1.29e-004 | 1.09e-002 | 1.23e-001 | 1.64e-001 |
    |     256 | 9.42e-003 | 9.52e-002 | 6.13e-002 | 2.69e-002 | 2.45e-004 | 2.55e-002 | 5.23e-001 | 7.41e-001 |
    |     512 | 1.60e-001 | 1.32e-001 | 3.00e-001 | 1.14e-001 | 4.58e-004 | 1.85e-001 | 2.66e+000 | 3.55e+000 |
    |    1024 | 4.82e-001 | 7.00e-001 | 9.57e-001 | 5.18e-001 | 9.76e-004 | 5.70e-001 | 1.37e+001 | 1.69e+001 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|

    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    | Fortran |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |   t_solve |     t_tot |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 0.00e+000 | 0.00e+000 | 1.74e-003 | 0.00e+000 | 0.00e+000 | 0.00e+000 | 8.68e-004 | 2.60e-003 |
    |      64 | 0.00e+000 | 0.00e+000 | 6.08e-003 | 8.68e-004 | 0.00e+000 | 0.00e+000 | 0.00e+000 | 6.94e-003 |
    |     128 | 3.47e-003 | 0.00e+000 | 2.26e-002 | 0.00e+000 | 8.68e-004 | 0.00e+000 | 1.39e-002 | 4.08e-002 |
    |     256 | 6.94e-003 | 1.91e-002 | 6.42e-002 | 8.68e-003 | 1.74e-003 | 0.00e+000 | 4.51e-002 | 1.46e-001 |
    |     512 | 3.04e-002 | 9.90e-002 | 2.20e-001 | 4.25e-002 | 0.00e+000 | 0.00e+000 | 1.80e-001 | 5.72e-001 |
    |    1024 | 1.40e-001 | 3.85e-001 | 8.98e-001 | 1.62e-001 | 0.00e+000 | 0.00e+000 | 7.71e-001 | 2.36e+000 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|



License
-------

Copyright (C) 2013-2016 J.S. Hysing.

Author: J.S. Hysing, info@featool.com

Keywords: Finite Element, FEM, Julia, Matlab, Octave, Fortran, Benchmark

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see
[http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).

Feat2D License
--------------

[FeatFlow](www.featflow.de) and the Feat library is Copyright (C) of
S. Turek et al. and may fall under a different license agreement.
