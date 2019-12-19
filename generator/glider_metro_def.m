function opts_metro=glider_metro_def(ifask,opts,permmode)
% opts_metro=glider_metro_def(ifask,opts,permmode) gets options for the Metropolis algorithm
%
%  ifask: 0 just to fill in with defaults, 1 to ask for new values
%  opts: initial option values, can be omitted
%  permmode: 0->sets up options for glider_metro (check states are changed)
%            1->sets up options for cstats_metro (checks are permuted)
%            2->sets up options for donut_metro (checks are permuted inside of matching donuts)
%            3->sets up options for donutw_metro (checks are permuted inside of matching donuts, weighted)
%            0 if omitted
%
%  10 Feb 2013:  add adjustnf_[enable|factor|show] to implement adjustment of flip number
%   to attempt to salvage iterations in which a large number of checks are
%   flipped, so taht they are too close and conflict with each other
%
%  18 Feb 2013:  added fields to enable a "donut" that does not fully go to
%  the corners, such as
%  * 1 1 1 *
%  1 1 0 1 1
%  * 1 1 1 *
%  this requires specifying the "excluded" regions, and also the gliders that need to be checked 
%  to ensure that block counts are unchanged,
%  (as it is not obvious how to compute these on the fly)
%     in the above example, unchanged blocks are
%     1 1      1    1 1 1    
%     1 1    1 1 1    1
%  donut_exclude is an n x 2 array of the excluded locations, using the
%  same conventions as donut.inds, defaults to []
%  donut_checkvecs is a cell array of n x 2 arrays of gliders to check for
%  preservation of blockcounts, defaults to cell(0)
%
%  above changes demonstrated in tbln_demo.
%
%   See also:  GLIDER_METRO, GLIDER_METRO_DEMO, CSTATS, CSTATS_METRO_DEMO, DONUT_METRO, DONUTW_METRO,
%     GETCORRS_FARNESS, TBLN_DEMO.
%
if (nargin<=1) opts=[]; end
if (nargin<=2) permmode=0; end
opts_metro=opts;
opts_metro=filldefault(opts_metro,'numiters',10000);
opts_metro=filldefault(opts_metro,'sampfreq_map',1000);
opts_metro=filldefault(opts_metro,'showfreq_map',1000);
opts_metro=filldefault(opts_metro,'sampfreq_stat',100);
opts_metro=filldefault(opts_metro,'pr_meth',1); %method for calculating probabilities: 0->global, 1->just local
opts_metro=filldefault(opts_metro,'nf_meth',0); %method for calculating number to flip: 0->absolute, 1->fraction
opts_metro=filldefault(opts_metro,'nf_abs',1); %absolute number to flip
opts_metro=filldefault(opts_metro,'nf_frac',0.01); %fractional number to flip
opts_metro=filldefault(opts_metro,'nf_dist',0); %shape of distribution of number to flip
opts_metro=filldefault(opts_metro,'mapubi_bc',1); %0 for rectangular boundary conditions, 1 for periodic
opts_metro=filldefault(opts_metro,'permmode',permmode);
opts_metro=filldefault(opts_metro,'temperature',0.001);
opts_metro=filldefault(opts_metro,'userlabel',[]); %show block count check messages
opts_metro=filldefault(opts_metro,'donut_multwtpwr',0.5);
opts_metro=filldefault(opts_metro,'donut_show',0); %show some calculations in donut mode
opts_metro=filldefault(opts_metro,'donut_showbc',0); %show block count check messages
opts_metro=filldefault(opts_metro,'donutw_maxsrch',40);  %maximum number of swaps that are searched on each iteration to find valid swaps
opts_metro=filldefault(opts_metro,'donutw_maxcomp',5); %maximum number of swaps that are compared on each iteration
opts_metro=filldefault(opts_metro,'donutw_fracbest',1); %fraction of times that the best weighted perm is taken
opts_metro=filldefault(opts_metro,'donutw_distpwr',2); %power law used to determine distance
opts_metro=filldefault(opts_metro,'map_show_name',[]);%prefix name for map window
opts_metro=filldefault(opts_metro,'map_show_gen0',1); %set to show map on generation 0 (if showfreq_map>0)
%
opts_metro=filldefault(opts_metro,'adjustnf_show',0); %set to 1 to show adjustment of flip number
opts_metro=filldefault(opts_metro,'adjustnf_enable',1); %set to 1 to allow adjustment of flip number, possibly "saving" iterations in which too many checks are swapped
opts_metro=filldefault(opts_metro,'adjustnf_factor',2); %thinning factor for nf on each iteration -- higher values are faster but less effective, must be strictly > 1
%
opts_metro=filldefault(opts_metro,'donut_exclude',[]); % list of locations in donut to exclude
opts_metro=filldefault(opts_metro,'donut_checkvecs',cell(0)); % list of gliders to check (must be supplied if donut_exclude is non-empty)
%
if (ifask==0)
    return
end
if (permmode==0) disp('Enter options for Metropolis algorithm, to run by changing states of checks'); end
if (permmode==1) disp('Enter options for Metropolis algorithm, to run by permuting states of checks'); end
if (permmode==2) disp('Enter options for Metropolis algorithm, to run by permuting insides of donuts'); end
opts_metro.numiters=getinp('number of iterations','d',[1 Inf],opts_metro.numiters);
opts_metro.sampfreq_map=getinp('map sampling frequency (0 for never)','d',[0 Inf],opts_metro.sampfreq_map);
opts_metro.showfreq_map=getinp('map showing frequency (0 for never)','d',[0 Inf],opts_metro.showfreq_map);
opts_metro.sampfreq_stat=getinp('statistics sampling frequency (0 for never)','d',[0 Inf],opts_metro.sampfreq_stat);
if (permmode==0) %for permutation and donut modes, calculation is always global
    opts_metro.pr_meth=getinp('probability calculation method (0->global, 1->local, 2->both and check)','d',[0 2],opts_metro.pr_meth);
end
opts_metro.nf_meth=getinp('flip or permute number calculation method (0->absolute, 1->fraction)','d',[0 1],opts_metro.nf_meth);
opts_metro.nf_abs=getinp('absolute number of checks to flip or permute','d',[1+permmode Inf],max(1+permmode,opts_metro.nf_abs));
opts_metro.nf_frac=getinp('fractional number of checks to flip or permute','f',[0 1],opts_metro.nf_frac);
disp(' options for flip or permute number distribution (mean specified above)')
disp(' 0->fixed');
disp(' 1->uniform');
disp(' 2->Bernoulli (near-Gaussian)');
disp(' 3->Poisson');
disp(' 4->exponential');
opts_metro.nf_dist=getinp('choice','d',[0 4],opts_metro.nf_dist);
opts_metro.mapubi_bc=getinp('0 for rectangular boundary conditions, 1 for periodic','d',[0 1],opts_metro.mapubi_bc);
if (getinp('1 for a user label','d',[0 1],min(1,length(opts_metro.userlabel))))
    opts_metro.userlabel=getinp('user label','s',[0 1],opts_metro.userlabel);
else
    opts_metro.userlabel=[];
end
if (permmode==1)
    opts_metro.temperature=getinp('temperature','f',[0 Inf],opts_metro.temperature);
end
if (ismember(permmode,[2 3]))
    opts_metro.donut_multwtpwr=getinp('power law for weighting multiplicities','f',[0 Inf],opts_metro.donut_multwtpwr);
    opts_metro.donut_show=getinp('1 to show intermediate calculations','d',[0 1],opts_metro.donut_show);
    opts_metro.donut_showbc=getinp('1 to show block count check','d',[0 1],opts_metro.donut_showbc);
end
if (permmode==3)
    opts_metro.donutw_maxsrch=getinp('maximum number of swaps searched on each iteration','d',[1 Inf],opts_metro.donutw_maxsrch);
    opts_metro.donutw_maxcomp=getinp('maximum number of swaps compared on each iteration','d',[1 Inf],opts_metro.donutw_maxcomp);
    opts_metro.donutw_fractbest=getinp('fraction of time to used best weighted pemutation','f',[0 1],opts_metro.donutw_fracbest);
    opts_metro.donutw_distpwr=getinp('power law used to determine distances','f',[0 Inf],opts_metro.donutw_distpwr);
end
return


