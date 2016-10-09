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


Results and Output
------------------

After running the tests the _output_ directory will contain results
files, _table.txt_ with tabulated times, and _results.html_ which
shows [Plotly](https://plot.ly/) comparison graphs. Sample output is
shown here below

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

<script src="https://cdn.plot.ly/plotly-latest.min.js"></script><div id="div1" style="width:800px;height:600px;"></div><div id="div2" style="width:800px;height:600px;"></div><div id="div3" style="width:800px;height:600px;"></div><div id="div4" style="width:800px;height:600px;"></div><div id="div5" style="width:800px;height:600px;"></div><div id="div6" style="width:800px;height:600px;"></div><div id="div7" style="width:800px;height:600px;"></div><div id="div8" style="width:800px;height:600px;"></div><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.002601,0.000867,0.007819,0.016495,0.076397,0.296113],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.001883,0.000705,0.001437,0.005214,0.030529,0.123996],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000167,0.000545,0.001796,0.009416,0.160459,0.481670],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000000,0.000000,0.003472,0.006944,0.030382,0.139757],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'Grid',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div1', data, layout);</script><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.002603,0.000870,0.005208,0.024311,0.089421,0.342604],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.001865,0.001090,0.006259,0.025341,0.097896,0.395367],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000581,0.002161,0.007345,0.095171,0.131756,0.699990],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000000,0.000000,0.000000,0.019097,0.098958,0.385417],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'Matrix Pointers',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div2', data, layout);</script><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.016496,0.019110,0.038199,0.134566,0.492242,1.932311],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.002810,0.006614,0.016168,0.055360,0.203359,0.815141],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000892,0.003558,0.014309,0.061276,0.299522,0.956528],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.001736,0.006076,0.022569,0.064236,0.220486,0.897569],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'Matrix Assembly',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div3', data, layout);</script><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.013026,0.012152,0.026046,0.096366,0.727509,7.326098],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.003261,0.004148,0.012350,0.050575,0.338916,3.304037],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000427,0.001860,0.006766,0.026866,0.113717,0.517585],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000000,0.000868,0.000000,0.008681,0.042535,0.162326],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'RHS Assembly',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div4', data, layout);</script><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.002606,0.005208,0.012151,0.039935,0.157999,0.669752],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.002664,0.005954,0.022820,0.089912,0.367051,1.508859],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000049,0.000077,0.000129,0.000245,0.000458,0.000976],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000000,0.000000,0.000868,0.001736,0.000000,0.000000],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'Boundary Conditions',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div5', data, layout);</script><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000000,0.003471,0.006082,0.032126,0.138037,0.538640],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.002115,0.009563,0.040710,0.147190,0.610436,2.516062],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000295,0.001477,0.010867,0.025498,0.184641,0.569815],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000000,0.000000,0.000000,0.000000,0.000000,0.000000],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'Sparse Allocation',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div6', data, layout);</script><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.004341,0.021985,0.093278,0.513568,2.742554,18.166663],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.004191,0.021653,0.090213,0.446080,2.529865,14.911531],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.007645,0.030280,0.122768,0.522664,2.660694,13.664485],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.000868,0.000000,0.013889,0.045139,0.179688,0.770833],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'Solver',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div7', data, layout);</script><script>var tr1={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.041674,0.063664,0.188784,0.857367,4.424159,29.272180],type:'scatter',name:'Octave',line:{color:'rgb(255,0,0)',width:2}};var tr2={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.018790,0.049727,0.189957,0.819671,4.178053,23.574994],type:'scatter',name:'Matlab',line:{color:'rgb(0,255,0)',width:2}};var tr3={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.010056,0.039956,0.163980,0.741137,3.551247,16.891049],type:'scatter',name:'Julia',line:{color:'rgb(0,0,255)',width:2}};var tr4={x:[32.000000,64.000000,128.000000,256.000000,512.000000,1024.000000],y:[0.002604,0.006944,0.040799,0.145833,0.572049,2.355903],type:'scatter',name:'Fortran',line:{color:'rgb(255,0,255)',width:2}};var layout={title:'Total Time',xaxis:{title:'1/h',tickvals:[32,64,128,256,512,1024]},yaxis:{title:'cpu [s]'}};data=[tr1,tr2,tr3,tr4];Plotly.newPlot('div8', data, layout);</script>


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
along with this program. If not, see http://www.gnu.org/licenses/.

Feat2D License
--------------

[FeatFlow](www.featflow.de) and the Feat library is Copyright (C) of
S. Turek et al. and may fall under a different license agreement.
