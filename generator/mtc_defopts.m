function newopts=mtc_defopts(opts)
%newopts=mtc_defopts(opts) defines options for mtc modules
% opts (may be omitted): non-default options
% newopts: full set of options, omitted values of opts filled in
%
%   See also:  FILLDEFAULT, MTC_AUGCOORDS, MTC_GETCORRS_P2X2, MTC_CGSSTRUCT_MERGE.
%
if (nargin<1) opts=[]; end
opts=filldefault(opts,'iflog',0); %1 for logging (general)
opts=filldefault(opts,'iflog_maxaug',0); %1 for logging in mtc_maxaug
opts=filldefault(opts,'iflog_selectmeths',0); %1 for logging in mtc_selectmeths
%
opts=filldefault(opts,'tol_match',10^-6); %tolerance for matches
opts=filldefault(opts,'tol_nonneg',10^-6); % tolerance for values being non-negative (values > -tol_nonneg considered non-negative)
opts=filldefault(opts,'tol_Pickard',10^-6); % tolerance for Pickard
opts=filldefault(opts,'boundary_makeconds',10^-4); % boundary allowance in makeconds; 10^-5 fails for ng=7, AB_1_3, dir=1
%
opts=filldefault(opts,'probe',0.1); %how far to probe to determine the available methods
%
opts=filldefault(opts,'nowarn_getcorrs',0); %suppress warnings in mtc_getcorrs_p2x2
opts=filldefault(opts,'nowarn_merge',0); %suppress warnings in mtc_cgsstruct_merge
opts=filldefault(opts,'nowarn_exceednaive',1); % suppress warnings in mtc_maxaug if naive bouindaries are exceeded
opts=filldefault(opts,'nowarn_minpgtzero',0); % suppress warnings in mtc_maxaug if p>0 on boundary
%
opts=filldefault(opts,'mismatch_action','keep_orig');
opts=filldefault(opts,'mismatch_action_options',{'keep_orig','keep_addin','keep_average','error'}); %possible actions on mismatch
%
opts=filldefault(opts,'niters_maxaug',20); % bisection iteration limit for mtc_maxaug
opts=filldefault(opts,'boundary_maxaug',0.05); %values for checking boundary in maxaug; can be an increasing sequence; 0 to disable
%
opts=filldefault(opts,'augcoords_bc_zshort',1); % allow shortcut in augcoords in bc plane if b or c is <=tol_match, taking one coord to be zero
opts=filldefault(opts,'augcoords_bc_mshort',1); % allow shortcut in augcoords in bc plane if b or c is <=tol_match, bypassing maxent for Markov
opts=filldefault(opts,'augcoords_bd_zshort',1); % allow shortcut in augcoords in bd plane if b or d is <=tol_match, taking one coord to be zero
%
opts=filldefault(opts,'selectmeths_onlyPickard',0); %only Pickard options
opts=filldefault(opts,'selectmeths_noPickard_Markov',0); %set to 1 to remove Pickard Markov methods (bc and bd planes)
opts=filldefault(opts,'selectmeths_noPickard_consensus',0); %set to 1 to remove Pickard consensus methods (bc plane)
opts=filldefault(opts,'selectmeths_noPickard_MarkovTheta',1); %set to 1 to remove Markov theta methods (bd plane)
opts=filldefault(opts,'selectmeths_noPickard_ZeroTheta',0); %set to 1 to remove zero theta methods (bd plane)
%
opts=filldefault(opts,'maxng_3x3stabcheck',4); %maximum number of gray levels to do 3x3 stability check for NoPickTT methods
opts=filldefault(opts,'maxng_teestabcheck',4); %maximum number of gray levels to do tee stability check for NoPickBT methods
%
opts=filldefault(opts,'alloyplot_blackdir_deg',240); %standard direction of the bias-to--black direction (p=[1 0 0]) in an alloy plot
opts=filldefault(opts,'alloyplot_graydir_deg',0); %standard direction of the bias-to-gray direction (p=[0 1 0]) axis in an alloy plot
%
newopts=opts;
end

 