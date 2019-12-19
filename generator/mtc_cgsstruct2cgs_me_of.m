function v=mtc_cgsstruct2cgs_me_of(x,probs,modblocks,mtcs,opts)
% v=mtc_cgsstruct2cgs_me_of(x,probs,modblocks,mtcs,opts)  is the objective
% function for mtc_cgsstruct2cgs_me
%
% x: a list of modulator values, a column of length (ng-1)^ncheck 
% probs: a set of block probabilities, a column of length ng^ncheck
% modblocks:  the impact of each modulator on each probability, size
%    ng^ncheck x (ng-1)^ncheck
% mtcs: structure returned by mtc_define
% opts: options 
%   opts.mode: what is to be calculated: 'entropy','maxminprob' (default: entropy)
%   opts.entropy_negprob:  value returned for entropy if any of the probabilities are < 0
%      defaults to -Inf
%
%   See also:  MTC_DEFINE, MTC_CGSSTRUCT2CGS_ME, HISTINFO.
%
ng=mtcs.ng; %number of gray levels
nchecks=mtcs.nchecks; %number of checks
if (nargin<=4)
    opts=[];
end
opts=filldefault(opts,'mode','entropy');
opts=filldefault(opts,'entropy_negprob',-Inf);
newprobs=probs+modblocks*x;
switch opts.mode
    case 'maxminprob'
        v=-min(newprobs);
    case 'entropy'
        if (min(newprobs)>=0)
            v=-histinfo(newprobs);
        else
            v=-opts.entropy_negprob;
        end
end
%disp([v x'])
return
