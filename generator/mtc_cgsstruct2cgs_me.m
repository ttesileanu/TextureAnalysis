function [cgs_vals_me,cgs_vals_loworder,cgs_loworder_list,opts_used]=mtc_cgsstruct2cgs_me(cgs_struct,mtcs,opts)
% [cgs_vals_me,cgs_vals_loworder,cgs_loworder_list,opts_used]=mtc_cgsstruct2cgs_me(cgs_struct,mtcs,opts)
% finds the cgs values that agrees with the input cgs_struct up to order nchecks-1
% and adjusts the highest-order coordinate groups to achieve maxmium entropy
%
% cgs_struct:  a structure with fields from mtcs.coord_groups{*}.name, each of which 
%   is a vector of size [1 mtcs.ng]; all groups of order nchecks are ignored
% mtcs: structure returned by mtc_define
% opts: options related to minimization algorithm
%     opts.tol: a tolerance, defaults to 10^-6
%     opts.verbose: set to 1 to log activity, defaults to 0
%     opts.verify_mrms: set to 1 to verify root-mean-squared probabilities
%         are at minimum with highest-order coords set to 0; defaults to 0
%     opts.fminsearch_*:  options to control fminsearch (number of iterations, tolerance, etc.)
%     opts.cgshint_struct: a structure of highest-order coordinate groups,
%         as a hint for a place to start.  Any lower-order coordinate
%         groups are ignored and replaced by low-order values from cgs_struct
%
% cgs_vals_me: an array of size [length(mtcs.coord_groups),ng] that agree
% with the inputs on order nchecks-1 and below, and is otherwise max entropy
% cgs_vals_loworder: an array of size [length(mtcs.coord_groups),ng] that agree
%   with the inputs on order nchecks-1 and below, has other coords set to unbiased values
% cgs_loworder_list:  a list of the entries that are present in cgs_struct
% opts_used: options used and some intermediate outputs
%    d.loworder.sumprob: sum of probabilities before optimization, should be 1
%    d.loworder.minprob: the minimum probability of any configuration before optimization
%    d.loworder.entropy2x2: the initial entropy in a 2x2 patch.  This is NOT the entropy 
%      per unit area of the texture, since it is not assumed to be Markov and the
%      entropies of the marginals (1x2, 2x1, 1x1) are not taken into account.
%      Nevertheless for the Pickard case, maximizing entropy2x2 is the same as
%      maximizing entropy per area, since the marginals are not affected by the
%      biases of the highest-order coord groups
%
% Strategy to find a feasible starting point (all probs >=0):
%  loworder: cgs_struct, with all high-order coefs set to 0
%  maxmin:   search from loworder to maximize the minimum probability
%  hint: using optional cgshint_struct 
%  hintmaxmin: start search from hint to maximize the minimum probability
%
% This does not check for missing or extra fields in cgs_struct
%
%   See also:  MTC_DEFINE, MTC_CGSSTRUCT2CGS, MTC_AUGCOORDS, MTC_PROBS2CGS, MTC_CGSS2PROBS, MTC_CGST2PROBS, MTC_CGS2CGSSTRUCT,
%     MTC_CGSSTRUCT2CGS_ME_OF, HISTINFO.
%
ng=mtcs.ng; %number of gray levels
nchecks=mtcs.nchecks; %number of checks
reshape_matrix=repmat(ng,[1 nchecks]);
if (nargin<=2)
    opts=[];
end
opts=filldefault(opts,'tol',10^-6);
opts=filldefault(opts,'verbose',0);
opts=filldefault(opts,'verify_mrms',0);
opts=filldefault(opts,'cgshint_struct',[]);
%
opts_used=opts;
verbose=opts.verbose;
cgshint_struct=opts.cgshint_struct;
ifhint=double(~isempty(cgshint_struct));
%
opts=filldefault(opts,'fminsearch_MaxFunEvals_mult',1000); %multiplies the number of modulators
opts=filldefault(opts,'fminsearch_MaxFunEvals_lower',1000); %lower limit
opts=filldefault(opts,'fminsearch_MaxFunEvals_upper',100000); %absolute limit
opts=filldefault(opts,'fminsearch_MaxIter_mult',400); %multiplies the number of modulators
opts=filldefault(opts,'fminsearch_MaxIter_lower',400); %lower limit
opts=filldefault(opts,'fminsearch_MaxIter_upper',10000); %absolute limit
opts=filldefault(opts,'fminsearch_TolFun',10^-6); %function tolerance (standard: 10^-4)
opts=filldefault(opts,'fminsearch_TolX',10^-6); %probability tolerance (standard: 10^-4)
opts=filldefault(opts,'fminsearch_display','off'); %'off','notify','final','iter'
%
counts=mtc_countparams(mtcs);
nmodulators=counts.coord_byorder(nchecks); %number of expected modulators
%
opts_fminsearch=optimset('fminsearch');
opts_fminsearch_use=opts_fminsearch;
%
%maximum number of function evals is bounded from below, and upper limit depends on number of modulators, up to some absolute limit
opts_fminsearch_use=optimset(opts_fminsearch_use,'MaxFunEvals',...
    max(opts.fminsearch_MaxFunEvals_lower,min(opts.fminsearch_MaxFunEvals_upper,nmodulators*opts.fminsearch_MaxFunEvals_mult)));
%maximum number of iterations is bounded from below, and upper limit depends on number of modulators, up to some absolute limit
opts_fminsearch_use=optimset(opts_fminsearch_use,'MaxIter',...
    max(opts.fminsearch_MaxIter_lower,min(opts.fminsearch_MaxIter_upper,nmodulators*opts.fminsearch_MaxIter_mult)));
opts_fminsearch_use=optimset(opts_fminsearch_use,'TolFun',opts.fminsearch_TolFun);
opts_fminsearch_use=optimset(opts_fminsearch_use,'TolX',opts.fminsearch_TolX);
opts_fminsearch_use=optimset(opts_fminsearch_use,'display',opts.fminsearch_display);
opts_used.fminsearch=opts_fminsearch_use;
%
opts_used.errors=[];
%
% parse the low-order starting point and, if available, the hint
cgs_vals_loworder=ones(length(mtcs.coord_groups),ng)/ng;
cgs_loworder_list=[];
cgs_vals_hint=ones(length(mtcs.coord_groups),ng)/ng;
for icgs=1:length(mtcs.coord_groups)
    cgsname=mtcs.coord_groups{icgs}.name;
    %
    cg=mtcs.coord_groups{icgs};
    coord_nums=cg.coord_num;
    lincomb_indices=mtcs.coord_indices(coord_nums);
    corder=sum(mtcs.lincomb_list(lincomb_indices(1),:)>0);
    %
    if isfield(cgs_struct,cgsname)
        if (corder<nchecks)
            cgs_vals_loworder(icgs,:)=cgs_struct.(cgsname);
            cgs_loworder_list=[cgs_loworder_list,icgs];
            cgs_vals_hint(icgs,:)=cgs_struct.(cgsname);
        else
            opts_used.errors=strvcat(opts_used.errors,sprintf('Error: parameter of order ncheck included (%s), will ignore',cgsname));
        end
    end
    if isfield(cgshint_struct,cgsname)
        if (corder==nchecks)
            %high-order coords from hint
            cgs_vals_hint(icgs,:)=cgshint_struct.(cgsname);
        end
    end
end
%
cgs_vals_me=[];
%
if (verbose>0)
    disp('Initial probabilities calculated with unbiased highest-order parameters.');
end
startpoints=[];
nstarts=0;
%first we have to see if this is feasible
%with the high-order coefs unbiased, what are the probabilites?
d.loworder=mtc_me_util(mtc_cgst2probs(cgs_vals_loworder,mtcs),opts);
opts_used.d.loworder=d.loworder;
nstarts=nstarts+1;
startpoints.type{nstarts}='loworder';
startpoints.entropy2x2(nstarts)=d.loworder.entropy2x2;
startpoints.probs{nstarts}=d.loworder.probs;
%
if (abs(sum(d.loworder.probs(:))-1)>opts.tol)
    opts_used.errors=strvcat(opts_used.errors,sprintf('Error: probabilities not normalized, sum to %10.8f, will ignore',sum(d.loworder.probs(:))));
end
if (d.loworder.minprob>=0)
    if (verbose>0)
        disp('Initial probabilities are all non-negative.');
    end
else
    if (verbose>0)
        disp('Initial probabilities have some negative values.  Will seek other values of highest-order parameters.');
    end
end
%
% do a brief analysis of the hint
%
d.hint=mtc_me_util(mtc_cgst2probs(cgs_vals_hint,mtcs),opts);
opts_used.d.hint=d.hint;
if (abs(sum(d.hint.probs(:))-1)>opts.tol)
    opts_used.errors=strvcat(opts_used.errors,sprintf('Error: hint probabilities not normalized, sum to %10.8f, will ignore',sum(d.hint.probs(:))));
end
if (ifhint)
    nstarts=nstarts+1;
    startpoints.type{nstarts}='hint';
    startpoints.entropy2x2(nstarts)=d.hint.entropy2x2;
    startpoints.probs{nstarts}=d.hint.probs;
end
%
% set up the modulators (in probability values)
%
counts=mtc_countparams(mtcs);
nmodulators=counts.coord_byorder(nchecks); %number of expected modulators
parity_table=int2nary([0:2^nchecks-1]',2,nchecks);
parities=1-2*mod(sum(parity_table,2),2); %1 if even, -1 if odd
parity_block=reshape(parities,repmat(2,[1 nchecks]));
opts_used.parity_block=parity_block;
%
if (ng>=3)
    modulator_desc=2+int2nary([0:((ng-1)^nchecks-1)]',ng-1,nchecks);
else
    modulator_desc=repmat(ng,[1 nchecks]);
end
opts_used.modulator_desc=modulator_desc;
if (size(modulator_desc,1)~=nmodulators)
    opts_used.errors=strvcat(opts_used.errors,sprintf('Error: modulators expected: %6.0f, found: %6.0f, cannot proceed',...
        nmodulators,size(modulator_desc,1)));
    disp(opts_used.errors(end,:));
    return
end
%this will indicate how each free parameter influences the block probabilities
modblocks=zeros(ng^nchecks,nmodulators);
for imod=1:nmodulators
    %disp(sprintf('imod: %6.0f modulator',imod));
    modulator=modulator_desc(imod,:);
    %disp(modulator);
    modrep=repmat(modulator,size(parity_table,1),1);
    %coords is set up to form a parallelopiped of dimension nchecks, with one corne
    %tethered at (1,1,1,1..1), and the far corner given by the modulator
    coords=zeros(size(parity_table));
    coords(find(parity_table==0))=1;
    coords(find(parity_table>=1))=modrep(find(parity_table>=1));
    %
    coord_nums=1+nary2int(coords-1,ng,2);
    %disp([coords,coord_nums,parities])
    %
    modblock=zeros(1,ng^nchecks);
    modblock(coord_nums)=parities;
    modblocks(:,imod)=modblock;
end
if (verbose>0)
    disp('Modulators calculated.');
end
% verify: before reshaping, abs(modblocks(:,imod)) should sum to 2^nchecks
sumcheck=sum(abs(modblocks),1);
if ~all(sumcheck==2^nchecks)
   opts_used.errors=strvcat(opts_used.errors,sprintf('Error: modulator sum check fails, cannot proceed'));
    disp(opts_used.errors(end,:));
   return
end
% reshape so that last dimension of modblocks indexes the modulators
% and the first nchecks dimensions is the modulator itself
modblocks_reshaped=reshape(modblocks,[reshape_matrix nmodulators]);
% verify: after reshaping, all marginals should sum to 0
marginalcheck=0;
for idim=1:nchecks
    marginalsum=sum(modblocks_reshaped,idim);
    if max(abs(marginalsum(:)))>0
       marginalcheck=1;
       opts_used.errors=strvcat(opts_used.errors,sprintf('Error: modulator marginal check fails on dim %1.0f, cannot proceed',idim));
       disp(opts_used.errors(end,:));
    end
end
if (marginalcheck>0)
    return
end
%
opts_used.modblocks=modblocks; %this is [ng^nchecks (ng-1)^nchecks]
if (verbose>0)
    disp('Modulators verified.');
end
%
if (opts.verify_mrms>0)
    %
    %calculate probabilities with least mean-square -- but this should be the
    %same as the initial probabilities, since the initial cooordinates are
    %orthogonal to the highest-order coords
    %
    modblocks_proj=modblocks*(inv(modblocks'*modblocks))*modblocks';
    d.mrms=mtc_me_util(reshape(d.loworder.probs(:)'-(d.loworder.probs(:)')*modblocks_proj,reshape_matrix),opts);
    opts_used.d.mrms=d.mrms;
    nstarts=nstarts+1;
    startpoints.type{nstarts}='mrms';
    startpoints.entropy2x2(nstarts)=d.mrms.entropy2x2;
    startpoints.probs{nstarts}=d.mrms.probs;
    if max(abs(d.mrms.probs(:)-d.loworder.probs(:)))>opts.tol
        opts_used.errors=strvcat(opts_used.errors,sprintf('Error: minimum rms probs do not match initial probs, ignored'));
        if (verbose>0)
            disp(opts_used.errors(end,:));
        end
    else
        if (verbose>0)
            disp('Minimum rms probs verified.');
        end
    end
end
%
%find the maxmin starting point by starting at loworder and maximizing the minimum probability
%
maxmin.fminsearch=mtc_me_fminsearch(d.loworder.probs,modblocks,mtcs,setfield([],'mode','maxminprob'),zeros(nmodulators,1),opts_fminsearch_use);
d.maxmin=mtc_me_util(reshape(maxmin.fminsearch.probvecs,reshape_matrix),opts);
d.maxmin.fminsearch=maxmin.fminsearch;
opts_used.d.maxmin=d.maxmin;
nstarts=nstarts+1;
startpoints.type{nstarts}='maxmin';
startpoints.entropy2x2(nstarts)=d.maxmin.entropy2x2;
startpoints.probs{nstarts}=d.maxmin.probs;
opts_used.startpoints=startpoints;
%
%find the hintmaxmin starting point by starting at hint and maximizing the minimum probability
%
if (ifhint)
    hintmaxmin.fminsearch=mtc_me_fminsearch(d.hint.probs,modblocks,mtcs,setfield([],'mode','maxminprob'),zeros(nmodulators,1),opts_fminsearch_use);
    d.hintmaxmin=mtc_me_util(reshape(hintmaxmin.fminsearch.probvecs,reshape_matrix),opts);
    d.hintmaxmin.fminsearch=hintmaxmin.fminsearch;
    opts_used.d.hintmaxmin=d.hintmaxmin;
    nstarts=nstarts+1;
    startpoints.type{nstarts}='hintmaxmin';
    startpoints.entropy2x2(nstarts)=d.hintmaxmin.entropy2x2;
    startpoints.probs{nstarts}=d.hintmaxmin.probs;
    opts_used.startpoints=startpoints;
end
%
%summarize all starting points
%
maxstart=max(startpoints.entropy2x2);
whichstart=find(startpoints.entropy2x2==maxstart);
if ~isempty(whichstart)
    whichstart=min(whichstart);
end
if (verbose>0)
    for istart=1:nstarts
        if (istart==whichstart)
            symb='*';
        else
            symb=' ';
        end
        disp(sprintf('Start point %12s: %s entropy %9.6f minprob %9.6f maxprob %9.6f',...
            startpoints.type{istart},symb,startpoints.entropy2x2(istart),min(startpoints.probs{istart}(:)),max(startpoints.probs{istart}(:))));
    end
end
if ~isempty(whichstart)
    if (verbose>0)
        disp(sprintf('Using %s as starting point to maximize entropy.',startpoints.type{whichstart}));
    end
    maxentropy.fminsearch=mtc_me_fminsearch(startpoints.probs{whichstart},modblocks,mtcs,...
        setfield([],'mode','entropy'),zeros(nmodulators,1),opts_fminsearch_use);
    d.maxentropy=mtc_me_util(reshape(maxentropy.fminsearch.probvecs,reshape_matrix),opts);
    d.maxentropy.fminsearch=maxentropy.fminsearch;
    opts_used.d.maxentropy=d.maxentropy;
    probs_me=d.maxentropy.probs;
    cgs_vals_me=mtc_probs2cgs(probs_me,mtcs);
    if (verbose>0)
        disp(sprintf('Final point             :   entropy %9.6f minprob %9.6f maxprob %9.6f',...
            d.maxentropy.entropy2x2,min(d.maxentropy.probs(:)),max(d.maxentropy.probs(:))));
    end
else
    opts_used.d.maxentropy.entropy2x2=NaN;
    opts_used.errors=strvcat(opts_used.errors,sprintf('Error:  cannot find a feasible starting point.'));
    if (verbose>0)
        disp(opts_used.errors(end,:));
    end
end
if isempty(opts_used.errors)
    if (verbose>0) 
        disp('No errors in mtc_cgsstruct2cgs_me.');
    end
else
    disp('Errors in mgc_cgsstruct2cgs_me:');
    disp(opts_used.errors);
end
return

function results=mtc_me_fminsearch(probs_init,modblocks,mtcs,opts_of,x0,opts_fminsearch_use)
[mods_result,fval,exitflag,output]=...
    fminsearch(@(x) mtc_cgsstruct2cgs_me_of(x,probs_init(:),modblocks,mtcs,opts_of),x0,opts_fminsearch_use);
results.probvecs=probs_init(:)+modblocks*mods_result;
results.mods=mods_result;
results.fval=-fval; %fmin miminimzes the negative of the lowest probability or entropy, flip to get maximum value
results.exitflag=exitflag;
results.output=output;
results.entropy=-mtc_cgsstruct2cgs_me_of(mods_result,probs_init(:),modblocks,mtcs,setfield(opts_of,'mode','entropy'));
return

function dnew=mtc_me_util(probs,opts)
%clean up slightly negative probabilities and calculate entropies, normalizations
dnew.probs=probs;
dnew.sumprob=sum(probs(:));
dnew.minprob=min(probs(:));
if (dnew.minprob>=-opts.tol)
    probs=max(probs,0);
    probs=probs/sum(probs(:));
    dnew.probs=probs;
    dnew.entropy2x2=histinfo(probs);
    dnew.sumprob=sum(probs(:));
    dnew.minprob=min(probs(:));
else
    dnew.entropy2x2=NaN;
end
return
 