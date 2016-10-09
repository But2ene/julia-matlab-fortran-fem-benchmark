% FEM Benchmark Matlab and Octave test run script.

% Copyright 2013-2016 Precise Simulation, Ltd.
% J.S. Hysing 161008.


clear all


addpath( 'src_matlab' )


% n0    = 16
% nlev  = 3
% nruns = 10
params = load( 'testrun_param.txt' );
n0    = params(1);
nlev  = params(2);
nruns = params(3);


fem_poisson(n0);
timings = zeros(nlev,8);
for ilev = 1:nlev
  n = n0*2^(ilev-1);
  timings_ilev = zeros(1,7);

  for i = 1:max(3,nruns)
    if i ~= 1 && i ~= nruns
      [foo,timings_i] = fem_poisson( n );
      timings_ilev = timings_ilev + timings_i;
    end
  end
  timings_ilev = timings_ilev/(nruns - 2);

  timings(ilev,1) = n0*2^(ilev-1);
  timings(ilev,2:end) = timings_ilev;
end


save( 'output/output_matlab.txt', 'timings', '-ascii' )


exit
