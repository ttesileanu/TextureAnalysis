function [map_synth,samp,optsused]=donut_metro(ng,donut,map_init,opts)
% [map_synth,samp,optsused]=donut_metro(ng,donut,map_init,opts) runs
% the Metropolis algorithm to find maximum-entropy textures with constrained
% 2x2 block counts, by swapping pixels inside of 3x3 donuts
%
%  2 Aug 2010:  branch added to catch possibility that there is nothing to swap
%  3 Aug 2010:  fixed logic related to cuts and choice, slightly increasing
%      probability that a high number of checks are swapped
%  4 Aug 2010:  added display of generation 0, and a prefix for the display
%      window name
%  5 Aug 2010:  began mods to that map_init can have a third dimension; donuts are mixed
%    across maps.  Calculation is done on a "concatenated" map, with
%    the stack of maps reorganized along dimension 2 Note that the intermediate saved maps
%    in "samp" keep this concatenation.  This requires opts.mapubi_bc=0.
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
%  donut_checkvecs is a cell array of n x 2 arrays of gliders to check for
%  preservation of blockcounts, defaults to cell(0)
%
%  above changes demonstrated in tbln_demo.
%
% 31 Uag 2015: begin debugging for ng>2
%
%  ng: number of gray levels, typically 2
%  donut: glider (as a structure created by glider_addcoords) specifying the donut
%  map_init:  map to start with, also determines size of final map
%  opts:  set of options, e.g., made by GLIDER_METRO_DEF with permmode=2
%
%  map_synth:  a synthesized map
%  samp:    intermediate results, including 
%     samp.map.map_synth:  snapshots of intermediate maps
%     samp.map.coords: coordinates changed at each snapshot
%     samp.stat.blockcounts: block counts within the glider
%     samp.bctable_setup:  setup of tables of block counts
%     samp.stat.bctable:  block counts on larger blocks
%     samp.sumperm:  summary of number of checks permuted (taking into account
%        that no perms are done if thre are block count failures)
%     samp.inside_ptrs, samp.barriers:  inside pointers and barriers, for checking that boundary
%        conditions and barriers between concatenated maps are in place
%
%  optsused:  option values used
%
%   See also:  DONUT_METRO_DEMO, GLIDER_METRO_DEF, GLIDER_MAPUBI, TEXTURELIB_MEE_DEMO, DONUT_METROW,
%     BTC_NOPICK, DILUMETRO_TEST, TBLN_DEMO.
%
nr=size(map_init,1);
nc_single=size(map_init,2);
nm=size(map_init,3);
nc=nc_single*nm;
if (nargin<=3)
    opts=[];
    ifask=1;
else
    ifask=0;
end
opts=glider_metro_def(ifask,opts,2);
optsused=opts;
samp=cell(0);
%concatenate maps to make the working array
map_synth_cat=zeros(nr,nc);
for im=1:nm
    map_synth_cat(:,[1:nc_single]+(im-1)*nc_single)=map_init(:,:,im);
end
%
% check arguments
%
%if (ng>2)
%    warning(sprintf('number of gray levels is %3.0f, algorithm may fail if >2',ng));
%end
if (ng<=max(map_synth_cat(:)))
    warning(sprintf('number of gray levels is %3.0f, but map appears to have more levels',ng));
end
%
m=donut.matrix;
inside=1-m;
if (sum(m(:)==1)~=(prod(size(m))-1)) | (sum(inside(:)==1)~=1)
    error('donut does not have the correct number of 1''s and 0''s');
    return
else
    invec=[find(sum(inside,2)==1),find(sum(inside,1)==1)];
    if any(invec<=1) | any(invec>=size(m))
        warning('donut inside is touching its edge.');
    end
end
% create vecs to check the block counts
%
% for standard donut, this is 2x2; for double donut, 3x3
%  **this approach does not work for "donut" templates that are missing a
%  corner.  For example, the donut template
%   1 1 1 
% 1 1 0 1 1
%   1 1 1
%
% would be expected to preserve a 2x2 block, AND a 1 x 3 block
% This is why missing-corner donuts are exclued by above test
%
if (isempty(opts.donut_checkvecs))
    if (~isempty(opts.donut_exclude))
        warning('checks are excluded from donut, but checkvecs is not supplied.');
    end
    checkvecs=cell(0);
    checkvec=[];
    for ix=0:invec(1)-1
        for iy=0:invec(2)-1
            checkvec=[checkvec;[ix iy]];
        end
    end
    checkvecs{1}=checkvec;
else
    checkvecs=opts.donut_checkvecs;
end
%
inside_ptrs=cell(0);
for c=1:2
    if (opts.mapubi_bc==0)
        inside_ptrs{c}=[1:(size(map_synth_cat,c)-size(m,c)+1)]+invec(c)-1;
    else
        inside_ptrs{c}=1+mod(invec(c)-1+[1:(size(map_synth_cat,c))],size(map_synth_cat,c));
    end
end
%add barriers between concatenated maps if more than one map in the stack
barriers=[];
if nm>1
    if (opts.mapubi_bc>0)
        warning('more than one map sent to donut_metro, but mapubi_bc is nonzero; proceeding as if it is 0.');
    end
    bsize=[min(inside_ptrs{2})-1,nc-max(inside_ptrs{2})];
    disp(sprintf('donut_metro is putting barriers (size %1.0f %1.0f) between %2.0f maps',bsize,nm));
    for im=1:nm
        barriers=[barriers, (im-1)*nc_single+[1:bsize(1) (nc_single-bsize(2)+1):nc_single]];
    end
end
samp.inside_ptrs=inside_ptrs;
samp.barriers=barriers;
%
% determine target mean number permuted on each iteration
%
nfmean=0;
if (opts.nf_meth==0)
    nfmean=min(nr*nc,max(1,opts.nf_abs));
end
if (opts.nf_meth==1)
    nfmean=min(nr*nc,max(1,round(nr*nc*opts.nf_frac)));
end
optsused.nf=nfmean; %mean number flipped on each iteration
if (nfmean<2)
    warning('less than two checks to be permuted on average per iteration.')
    map_synth=donut_metro_decat(map_synth_cat,nr,nc_single,nm);
    return;
end
hw=waitbar(0,sprintf('Doing %g Metropolis iterations with permuting',opts.numiters));
if ~isempty(opts.userlabel)
    set(hw,'Name',opts.userlabel);
end
if (opts.showfreq_map>0)
    hms=figure;
    set(gcf,'Position',[100 100 800 800]);
    set(gcf,'NumberTitle','off');
    %show the starting map
    if (opts.map_show_gen0>0)
        set(hms,'Name',deblank(cat(2,opts.userlabel,sprintf('%s iter %g',opts.map_show_name,0))));
        imagesc(map_synth_cat,[0 ng-1]);colormap gray;title(get(gcf,'Name'));axis equal;axis tight;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
    end
end
%
%set up some tables for block counts: extremal shapes of rectangular blocks
mct=floor(log(prod(size(map_synth_cat)))/log(ng)); %no more kinds of blocks than checks
bctable=[];
ibc=0;
for rmax=1:floor(sqrt(mct))
    cmax=floor(mct/rmax);
    ibc=ibc+1;
    bctable.rcmax(ibc,:)=[rmax,cmax];
    if (rmax<cmax)
        ibc=ibc+1;
        bctable.rcmax(ibc,:)=[cmax,rmax];
    end
end
bctable.rcmax=sortrows(bctable.rcmax);
for ibc=1:size(bctable.rcmax,1)
    bctable.label{ibc}=sprintf('%g x %g',bctable.rcmax);
    bctable.gstruct{ibc}=glider_addcoords(setfield([],'matrix',ones(bctable.rcmax(ibc,:))));
end
samp.bctable_setup=bctable;
samp.ng=ng;
samp.donut=donut;
%
nsamps_map=0;
nsamps_stat=0;
tic;
goodct=0; %count of number of successful iterations (block counts preserved) 
savect=0; %count (included in goodct) of iterations that are salvaged by thinning procedure 10_Feb_2013
badct=0; %counts of number of unsuccessful iterations (swapped checks are too close, and donuts overlap)
goodct_tot=0; %count of number of checks involved in permutations on successful iterations
savect_tot=0; %count (included in in goodct_tot) of checks involved in permutations on salvaged iterations
badct_tot=0; %count of number of checks involved in permutations on unsuccessful iterations
;
iters_noswap=0;
%
dtemplate=donut.matrix;
if ~isempty(opts.donut_exclude)
    for k=1:size(opts.donut_exclude,1)
        dtemplate(opts.donut_exclude(k,1),opts.donut_exclude(k,2))=0;
    end
    [d1,d2]=find(dtemplate>0);
    donut_mustmatch=[d1,d2];
else
    donut_mustmatch=donut.inds;
end
optsused.donut_mustmatch_matrix=dtemplate;
optsused.donut_mustmatch=donut_mustmatch;
%
for iter=1:opts.numiters
    %
    nf_target=0;
    while (nf_target<2)
        switch opts.nf_dist
            case 0 %constant
                nf_target=nfmean;
            case 1 %uniform
                nf_target=1+floor(2*nfmean*rand(1));
            case 2 % Bernoulli
                nf_target=sum(rand(1,2*nfmean)<0.5);
            case 3 % Poisson
                nf_target=sum(rand(1,1000*nfmean)<0.001);
            case 4 % exponential
                nf_target=round(-nfmean*log(rand(1)));
        end
    end
    insides=map_synth_cat(inside_ptrs{1},inside_ptrs{2});
    %insides assigns a unique number to the inside of the donut
    %ubi assigns a unique number to the outside of the donut
    [bci,ubi]=glider_mapubi(map_synth_cat,donut_mustmatch,ng,opts); %get the block counts and unique block indices
    %ubh assigns a unique number to entire donut, but, insides are set to Inf on barrier columns
    %to prevent swapping, since they don't result in counts
    map_synth_cat_bar=map_synth_cat;
    map_synth_cat_bar(inside_ptrs{1},barriers)=Inf;
    [bch,ubh]=glider_mapubi(map_synth_cat_bar,cat(1,donut_mustmatch,invec),ng,opts); %get the block counts and unique block indices
    bci_wrap=reshape(bch,[length(bci),ng])'; %each column corresponds to an entry in bci(1,:)
    %
    masked_allinds=find(ubh(:)==Inf);
    % we want to find the block indices that have large entries in ubi but also
    % more than one nonzero element in bci_wrap
    %
    %Aug 5, 2010: changed logic so that one looks for both elements of bci_wrap nonzero,
    % rather than max of bci_wrap less than bci (so that it also works if some pixels are masked out)
    %bci_eligible=find(max(bci_wrap)<bci); %changed Aug 5 this 
    %bci_eligible=find(min(bci_wrap)>0); %same logic as achanged Aug 5, so that it looks for nonzeros
    bci_eligible=find(sum(double(bci_wrap>0))>=2); %31 Aug 15, needed for ng>2
    %bci_eligible
    %[bci_wrap(:,bci_eligible);bci(bci_eligible)]
    %
    %if (isempty(barriers))
    %    bci_eligible_old=find(max(bci_wrap)<bci);
    %    if length(bci_eligible_old)==length(bci_eligible)
    %        if all(bci_eligible_old==bci_eligible)
    %            disp(' check is OK')
    %        else
    %            bci_eligible
    %            bci_eligible_old
    %        end
    %    else
    %        bci_eligible
    %        bci_eligible_old
    %    end
    %end
    %
    multmax=max(bci(bci_eligible));
    nf_try=nf_target;
    nf=min(nf_try,multmax);
    ubimult=bci_eligible(find(bci(bci_eligible)>=nf));
    %ubimult
    %
    if size(nf,2)==0 %branch added 2 Aug 2010
        iters_noswap=iters_noswap+1;
        if iters_noswap==1
            disp(sprintf(' iteration %7.0f: no blocks to swap (no other such exceptions individually logged)',iter));
        end
        permlist=[];
        newperm=[];
        allx=0;
        ally=0;
        nf=0;
        nfmod=0;
    else
        nbci_okmult=length(find(bci>=nf)); %number of unique block indices with sufficiently large multiplicities
        nbci_okdiff=length(bci_eligible); %number of unique block indices that have different insides
        nbci_okboth=length(ubimult); % number of unique block indices satisfying both criteria
        %disp(sprintf(' nf_try %5.0f nf %5.0f length(bci>=nf) %7.0f length(bci_eligible) %7.0f length(ubimult) %7.0f',...
        %    nf_try,nf,nbci_okmult,nbci_okdiff,nbci_okboth))
        %[ubimult',bci(ubimult)',bci_wrap(:,ubimult)']
        % ubimult lists the unique block indices which have at least nf copies, and not identical insides
        %
        %now choose a particular unique block, using multiplicities as weighted by donut_multwtpwr
        bcilist=bci(ubimult);
        bcilist_wt=bcilist.^opts.donut_multwtpwr;
        cuts=cumsum([0 bcilist_wt])/sum(bcilist_wt);
        choice=min(length(bcilist_wt),max(find(cuts<rand(1)))); %changed 3 Aug 2010
        %choice=min(length(bcilist_wt)-1,max(find(cuts<rand(1))));
        %
        %disp('cuts')
        %disp(cuts)
        %disp('bcilist_wt as a row');
        %disp((bcilist_wt(:))');
        %disp('length(bcilist_wt)')
        %disp(length(bcilist_wt))
        %disp('choice')
        %disp(choice)
        %
        ubi_choice=ubimult(choice); %this is the unique block index
        if (opts.adjustnf_show==1)
            disp(sprintf(' iteration %7.0f: nf_target %7.0f mult_max %7.0f nf %7.0f ubi_choice %10.0f',iter,nf_target,multmax,nf,ubi_choice));
        end
        %disp([bci(ubi_choice),bci_wrap(:,ubi_choice)']); %counts of blocks with index ubi_choice, with any (bci) or each (bci_wrap) check inside
        allinds=find(ubi(:)==ubi_choice);
        %disp('allinds')
        %disp(allinds')
        %disp('intersection with allinds')
        %disp(intersect(allinds',masked_allinds'))
        allinds=setdiff(allinds,masked_allinds);
        %disp('masked allinds')
        %disp(allinds')
        %
        allinsides=insides(allinds);
        [allx,ally]=ind2sub(size(insides),allinds);
        mallx=allx;
        mally=ally;
        if (opts.mapubi_bc==1) %offset the indices by 1 in periodic-boundary-condition mode
            allx=1+mod(allx,size(map_synth_cat,1));
            ally=1+mod(ally,size(map_synth_cat,2));
        end
        allxd=1+mod(allx+invec(1)-2,size(map_synth_cat,1)); %locatoin of center pixel
        allyd=1+mod(ally+invec(2)-2,size(map_synth_cat,2));
        if (opts.donut_show==1)
            xl=(size(donut.matrix,1)-1)/2;
            yl=(size(donut.matrix,2)-1)/2;
            for k=1:size(allx,1)
                disp([k allx(k) ally(k) ubi(mallx(k),mally(k))])
                if (allxd(k)>xl & allxd(k)<size(map_synth_cat,1)-xl  & allyd(k)>yl & allyd(k)<size(map_synth_cat,2)-yl)
                    disp(map_synth_cat(allxd(k)+[-xl:xl],allyd(k)+[-yl:yl]));
                end
            end
        end
        %[[1:length(allinds)]' allinds allx ally allinsides] %show all indices and coordinates
        %
        %choose a subset of nf gliders that have different insides:
        %determine counts of how many of each kind of donut to take
        %aim to divide the nf checks equally among ng possibilities,
        %but there are only bci_wrap(ig,ubi_choice) of each kind
        %have to do some kind of water-filling here
        ncounts=zeros(1,ng);
        navail=bci_wrap(:,ubi_choice)';
        %unassigned=nf;
        unassigned=sum(navail); %change 5 Aug 2010 in case any of the blocks to be switched were masked out
        nfmod=sum(navail); %as above
        %disp('bci_wrap')
        %disp(bci_wrap)
        %disp('ubi_choice')
        %disp(ubi_choice)
        %disp('navail')
        %navail
        while unassigned>0
            usable=find(navail>0);
            fill=min(navail(usable));
            ntake=min(fill,floor(unassigned/length(usable)));
            ntake1=min(fill,ceil(unassigned/length(usable)));
            ncounts(usable)=ncounts(usable)+ntake;
            navail(usable)=navail(usable)-ntake;
            unassigned=unassigned-ntake*length(usable);
            if (ntake1>ntake) %should lead to termination on next step
                extras=mod(unassigned,length(usable));
                q=randperm(length(usable));
                uneven=q([1:extras]);
                ncounts(usable(uneven))=ncounts(usable(uneven))+1;
                navail(usable(uneven))=navail(usable(uneven))-1;
                unassigned=unassigned-length(uneven);
            end
        end
        %select, for each ig, ncounts(ig) checks that have value ig
        if opts.donut_show==1
            disp('bci_wrap and ncounts');
            disp(bci_wrap(:,ubi_choice)');
            disp(ncounts)
        end
        nthin=1; %this gets progressively increased to adjust the effective value of ncounts
        %ncnz=min(ncounts(ncounts>0));
        ncnz=max(ncounts);
        if_swapok=0;
        %iterate until
        % nthin is too high
        % OR swapping is successful
        % OR adjustment of nf is disabled, and an adjustment has been attempted
        while (nthin<=ncnz) & (if_swapok==0) & ((opts.adjustnf_enable==1) | (nthin==1)) 
            ncounts_use=ceil(ncounts/nthin);
            permlist=[];
            for ig=1:ng
                if (ncounts_use(ig)>0)
                    avail_ptrs=find(allinsides==(ig-1));
                    select=randperm(length(avail_ptrs));
                    used_ptrs=avail_ptrs(select(1:ncounts_use(ig)));
                    permlist=[permlist;used_ptrs]; %list of indices to permute
                end
            end
            %[[1:length(allinds)]' allinds allx ally allinsides] %show all indices and coordinates
            %ncounts
            %permlist
            %
            newperm=[1:length(permlist)]; %find a random permutation that changes at least something
            while (all(newperm==[1:length(permlist)]))
                newperm=randperm(length(permlist));
            end
            %iteration step: global method
            %calculate new map
            map_new=map_synth_cat;
            for k=1:length(permlist)
                map_new(allxd(permlist(k)),allyd(permlist(k)))=...
                    map_synth_cat(allxd(permlist(newperm(k))),allyd(permlist(newperm(k))));
            end
            if (opts.adjustnf_show==1)
                nchanged=(sum(~(double(map_new(:)==map_synth_cat(:)))));
                disp(sprintf(' iteration %7.0f : npixels changed %7.0f, nfmod %7.0f, length(permlist) %7.0f ncnz %7.0f nthin %7.0f ncounts_use %7.0f %7.0f %7.0f',...
                    iter,nchanged,nfmod,length(permlist),ncnz,nthin,ncounts_use));
            end
            if (nthin==1) %only need to do this the first time through
                bc_old=[];
                for k=1:length(checkvecs)
                    bc_old{k}=glider_mapubi(map_synth_cat,checkvecs{k},ng,opts);
                end
            end
            bc_new=[];
            checks_ok=repmat(-1,[1 length(checkvecs)]); %-1 means not checked yet, 0 means mismatch, 1 means counts all match (and are OK)
            all_checks_ok=1;
            k=1;
            while (k<=length(checkvecs)) & all_checks_ok==1;
                bc_new{k}=glider_mapubi(map_new,checkvecs{k},ng,opts);
                if (all(bc_old{k}==bc_new{k}))
                    checks_ok(k)=1;
                else
                    all_checks_ok=0;
                    checks_ok(k)=0;
                %    k
                %    [bc_old{k};bc_new{k}]
                end
                k=k+1;
            end    
            %this could fail if some of the pixels to be swapped are adjacent,
            %e.g., swapping the interior two pixels in
            %[CCCC;BBAA;DDDD]
            if all_checks_ok==1
                map_synth_cat=map_new;
                goodct=goodct+1;
                goodct_tot=goodct_tot+nfmod;
                if_swapok=1; %we keep this iteration and exit the loop successfully
                if (nthin>1)
                    savect=savect+1;
                    savect_tot=savect_tot+round(nfmod/nthin);
                end
            end
            nthin=nthin*opts.adjustnf_factor;
        end %loop on nthin<=ncnz, if_swapok==0
        if (if_swapok==0)
            if (opts.donut_showbc) | (opts.adjustnf_show)
                disp(sprintf(' iteration %7.0f: block count check fails so swapping not done',iter));
                if (length(checkvecs)>1)
                    disp('checks_ok:');
                    disp(checks_ok)
                end
            end
            badct=badct+1;
            badct_tot=badct_tot+nfmod;
        end
    end %end of branch if there is nothing to swap
    %
    %sample the statistics
    if (mod(iter,opts.sampfreq_stat))==0
        nsamps_stat=nsamps_stat+1;
        samp.stat.iter(nsamps_stat)=iter;
        for ibc=1:size(bctable.rcmax,1)
            samp.stat.bctable{ibc}(:,nsamps_stat)=glider_mapubi(map_synth_cat,bctable.gstruct{ibc}.inds,ng,opts);
        end
        samp.stat.coords{nsamps_stat}=[allx(permlist),ally(permlist)];
        samp.stat.newperm{nsamps_stat}=newperm;
        samp.stat.nf(nsamps_stat)=nfmod;
    end
    %sample the map
    if (mod(iter,opts.sampfreq_map))==0
        nsamps_map=nsamps_map+1;
        samp.map.iter(nsamps_map)=iter;
        samp.map.map_synth(:,:,nsamps_map)=map_synth_cat;
    end
    %show the map
    if (mod(iter,opts.showfreq_map))==0
        figure(hms);
        set(hms,'Name',deblank(cat(2,opts.userlabel,sprintf('%s iter %g',opts.map_show_name,iter))));
        imagesc(map_synth_cat,[0 ng-1]);colormap gray;title(get(gcf,'Name'));axis equal;axis tight;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
    end
    waitbar(iter/opts.numiters,hw);
end
eltime=toc;
disp(sprintf('elapsed time %8.3f sec',eltime));
disp(sprintf('%g iterations, req target %10.4f perms per iteration, map is %g by %g',...
    opts.numiters,nfmean,size(map_synth_cat)));
disp(sprintf('yielding target of %12.3f perms per check',opts.numiters*nfmean/prod(size(map_synth_cat))));
disp(sprintf('%g iterations with no change in block counts (%g saved), %g with changes (and suppressed)',goodct,savect,badct));
disp(sprintf('%g iterations with no blocks to swap',iters_noswap));
samp.sumperm.iters_perm_done=goodct;
samp.sumperm.iters_perm_done_save=savect;
samp.sumperm.iters_perm_not_done=badct;
samp.sumperm.tot_perm_done=goodct_tot;
samp.sumperm.tot_perm_done_save=savect_tot;
samp.sumperm.tot_perm_not_done=badct_tot;
optsused.elapsed_time=eltime;
close(hw);
if opts.showfreq_map>0
    close(hms)
end
map_synth=donut_metro_decat(map_synth_cat,nr,nc_single,nm);
return
%
function map_synth=donut_metro_decat(map_synth_cat,nr,nc_single,nm)
%de-concatenate the maps
map_synth=zeros(nr,nc_single,nm);
for im=1:nm
    map_synth(:,:,im)=map_synth_cat(:,[1:nc_single]+(im-1)*nc_single);
end
return


    
