function [img,optsused,errs,metro]=mtc_flfs_makemaps(method,opts,area,dict,mtcopts)
% [img,optsused,errs,metro]=mtc_flfs_makemaps(method,opts,area,dict,mtcopts)
%  creates multilevel maps via the falling-leaves/falling sticks construction.
%  
%  Intended to be generalizable but at present (Oct 2015) it only does falling sticks
%
% method:  a structure returned by mtc_augcoords
%   method.name='FLFS'
%   method.flfs:  has information for the construction from mtc_augcoords; subfields are:
%         flfs: [1x1 struct], the flfs setup
%         flfs_configs: [1x14 struct], from flfs_enumerate
%         flfs_factors: [1x1 struct], from flfs_bc_factors
%         flfs_eivec: [14x1 double], the stable set of probabilities of configurations
%         flfs_weight: [1x1 struct], fields b, c, bb, and cc, indicating how much these correlations are diluted
%         flfs_p_hv:  the relative probability of horizontal stick (b); vertical has probability 1-flfs_p_hv
%         flfs_probs: {[ngxng double]  [ngxng double]}: 2-point probabilities for horiz (b) and vertical (c) sticks
% opts.show: set to 1 to draw the image, 2 to show the image and also intermediates
% opts.err_rept: 0 to abort on errors, 1 to report them into errs
% opts.nmaps: number of maps to make (defaults to 1)
% opts.burnin:  number of "burn-in" rows
% opts.minarea: minimum area to make, to allow for adequate Metropolis scrambling
% opts.metro_opts:  Metropolis options
% opts.metro_show: set to show Metropolis details
% opts.verbose: 1 for a log, 2 for a detailed log
% opts.onaxis: 1 if on axis (i.e., at most one param value is nonzero); this bypasses Metropolis
%
% area: size of img generated, defaults to [256 256]
%
% dict: dictionary of corrs; created with btc_define if not passed.
% mtcopts: options structure set by MTC_DEFOPTS, if passed
%
%  output
%
% img: the image; size(img) = area
% optsused: the options used
% errs:  string to report errors
% metro: results of call to donut_metro
%
%   See also MTC_AUGCOORDS, BTC_MAKEMAPS, FLFS_DEFINE, FLFS_BC_FACTORS, GENMRFMG_1D.
%

flfs_coverage_default=1; %number of times that a check needs to be covered by a leaf or stick
% can change this in method.flfs.flfs.coverage; default changed from 3 on 11 Aug 2016
flfs_reservoir_size_default=2; %size of the reservoir of leaves or sticks, as a multiple of
%  the area of the map being generated
%  can change this in method.flfs.flfs.reservoir_size; 
flfs_reservoir_multistart_default=5; %threshold for creating multiple reservoirs, if
%  the transition probability is low.  If minimum transition probability*reservoir length does not 
%  exceed this, then ng reservoirs are generated, one for each starting gray level value
%  can change this in method.flfs.flfs.reservior_multistart_default; 
flfs_minarea_multiplier_default=4; %multiplier to compute minimum size of map to make
%
metro=[];
ard=[256 256];
if (nargin >= 2)
   if (length(area) >= 2)
      ard=[area(1) area(2)];
   else
      ard=[area(1) area(1)];
   end
end
p_hv=method.flfs.flfs_p_hv; % relative probability of horizontal stick (b); vertical has probability 1-flfs_p_hv
probs=method.flfs.flfs_probs; % {[ngxng double]  [ngxng double]}: 2-point probabilities for horiz (b) and vertical (c) sticks
ng=size(probs{1},1);
%
opts=filldefault(opts,'show',0);
opts=filldefault(opts,'err_rept',0);
opts=filldefault(opts,'nmaps',1);
opts=filldefault(opts,'burnin',8);
opts=filldefault(opts,'minarea',flfs_minarea_multiplier_default*ng^8);
opts=filldefault(opts,'verbose',0);
opts=filldefault(opts,'metro',[]);
mtcopts=mtc_defopts(mtcopts);
opts.mtcopts=mtcopts;
%
optsused=opts;
%
img=zeros([ard opts.nmaps]);
stats=[];
errs=[];
%
if ~(strcmp(method.name,'FLFS'))
    estring='mtc_flfs_makemaps called, but method is not FLFS';
    if (opts.err_rept==0)
        error(estring);
    else
        errs=strvcat(errs,estring);
        return
    end
end
%
%determine number of gray levels and get flfs params
%
flfs=method.flfs.flfs; 
flfs=filldefault(flfs,'coverage',flfs_coverage_default);
flfs=filldefault(flfs,'reservoir_size',flfs_reservoir_size_default);
flfs=filldefault(flfs,'reservoir_multistart',flfs_reservoir_multistart_default);
%
coverage=flfs.coverage;
reservoir_size=flfs.reservoir_size;
ncats=length(flfs.cats);
nsubcats=length(flfs.subcats);
%generate probabilities of each category
p_cats=[p_hv,1-p_hv]; %probability of each category, specific to b,c FLFS model
%
% determine layout of multiple-map display if needed
%
[nr,nc]=nicesubp(opts.nmaps,ard(2)/ard(1)); 
layout=[nr*ard(1) nc*ard(2)];
ard_super_basic=layout+opts.burnin;
if prod(ard_super_basic)>=opts.minarea
    ard_super=ard_super_basic;
else
    fneeded=sqrt(opts.minarea/prod(ard_super_basic));
    ard_super=ceil(ard_super_basic*fneeded);
end
if (opts.verbose)
    disp(sprintf('map size  : %4.0f x %4.0f pixels',ard));
    disp(sprintf('layout in supermap: %4.0f x %4.0f maps   (%5.0f requested)',nr,nc,opts.nmaps));
    disp(sprintf('supermap to make  : %4.0f x %4.0f pixels (including burnin of %4.0f pixels but not minimum area constraint, area is %12.0f pixels)',...
        ard_super_basic,opts.burnin,prod(ard_super_basic)));
    disp(sprintf('supermap to make  : %4.0f x %4.0f pixels (including burnin of %4.0f pixels and also minimum area constraint of      %12.0f pixels)',...
        ard_super,opts.burnin,opts.minarea));
end
%
% build the reservoir of falling_sticks
%
reservoir_start_max=ceil(prod(ard_super)*reservoir_size);
reservoir_length=reservoir_start_max+max(ard_super)-1;
%
% if the block probabilities contain values that are very close to zero,
% then make several reservoirs with different start points
%
reservoir_count=[1 1];
for irc=1:2
    if (min(method.flfs.flfs_probs{irc}(:))*reservoir_length<flfs.reservoir_multistart)
        reservoir_count(irc)=ng;
    end
end
for irc=1:2
    reservoir{irc}=cell(0);
    for ires=1:reservoir_count(irc)
        if (reservoir_count(irc)>1)
            firstpxl=mod(ires-1,ng);
        else
            firstpxl=[];
        end
        reservoir{irc}{ires}=genmrfmg_1d(reservoir_length,ones(1,ng)/ng,method.flfs.flfs_probs{irc},firstpxl); %row or column process
        if (irc==1)
            reservoir{irc}{ires}=reservoir{irc}{ires}';
        end
        if (opts.verbose)
            disp(sprintf(' for dim %2.0f, reservoir size = %5.0f %5.0f, reservoir %2.0f of %2.0f made',...
                irc,size(reservoir{irc}{ires}),ires,reservoir_count(irc)));
        end
    end
end
optsused.flfs_stats.reservoir_length=reservoir_length;
optsused.flfs_stats.reservoir=reservoir;
%
%here do the falling-sticks construction in the supermap
%
imsuper=zeros(ard_super);
rcmax=max(size(imsuper)); %need a square array to avoid biasing row and column probabilities
dropped=0; %number of leaves/sticks dropped
dropped_cat=zeros(1,ncats); %number of leaves/sticks dropped in each category
misses=0; %counts number of times a leaf or stick misses imsuper entirely
cat_list=[]; %cat_list(ihit)=category of the ith hit
least_covered=0;
coverage_pointers=zeros(ard_super); %which leaf/stick is in which location
coverage_depth=zeros(ard_super);
coverage_coords=[]; %location of each leaf/stick
coverage_types=zeros([ard_super 2]); %one map for each coverage type
while least_covered<coverage
    %get a random number distributed as p_subcats
    icat=1+sum(rand(1)>cumsum(p_cats));
    cat_type=flfs.cats{icat}.type;
    hit=0;
    mask=zeros(size(imsuper)); %to accommodate more general leaves and sticks
    switch cat_type
        case 'row'
            ir=ceil(rcmax*rand(1));
            if (ir<=size(imsuper,1))
                hit=1;
                mask(ir,:)=1;
                coords=[ir,NaN];
            end
        case 'col'
            ic=ceil(rcmax*rand(1));
            if (ic<=size(imsuper,2))
                hit=1;
                mask(:,ic)=1;
                coords=[NaN,ic];
            end
    end
    if (hit==1)
        dropped=dropped+1;
        dropped_cat(icat)=dropped_cat(icat)+1;
        cat_list(dropped)=icat;
        coverage_pointers=coverage_pointers.*(1-mask)+dropped*mask;
        coverage_depth=coverage_depth+mask;
        least_covered=min(coverage_depth(:));
        coverage_coords(dropped,:)=coords;
    else
        misses=misses+1;
    end
end
%now find the contiguous regions of each leaf/stick that is uncovered
dropped_unique=unique(coverage_pointers(:));
for mord=1:length(dropped_unique)
    dropped=dropped_unique(mord);
    dcat=cat_list(dropped);
    dcoord=coverage_coords(dropped,:);
    cat_type=flfs.cats{dcat}.type;
    switch cat_type %not generic for falling leaves
        case 'row'
            stick=coverage_pointers(dcoord(dcat),:);
        case 'col'
            stick=coverage_pointers(:,dcoord(dcat));
    end
    stick=stick(:)'; %make it horizontal
    overlays=setdiff(unique(stick(:)),dropped);
    %disp('stick')
    %disp(stick);
    overlay_locs=[0 find(stick~=dropped) length(stick)+1];
    for iseg=1:length(overlay_locs)-1
        seg_ends=[overlay_locs(iseg)+1 overlay_locs(iseg+1)-1];
        seg_length=seg_ends(2)-seg_ends(1)+1;
        %disp(sprintf(' seg %4.0f: ends: %5.0f %5.0f    length %5.0f',iseg,seg_ends,seg_length));
        if (seg_length>0)
            mask=zeros(ard_super);
            switch cat_type %not generic for falling leaves
                case 'row'
                    mask(dcoord(dcat),seg_ends(1):seg_ends(2))=1;
                case 'col'
                    mask(seg_ends(1):seg_ends(2),dcoord(dcat))=1;
            end
            mlocs=find(mask(:)==1);
            ctype_1d=coverage_types(:,:,dcat);
            ctype_1d(mlocs)=1+rand(1); %asssign a constant random value to a leaf or stick
            coverage_types(:,:,dcat)=reshape(ctype_1d,ard_super);
            %
            %insert into img
            %
            istart=ceil(rand(1)*reservoir_start_max); %random starting point
            iresno=ceil(rand(1)*reservoir_count(dcat));
            leaf_stick=reservoir{dcat}{iresno}(istart+[0:seg_length-1]);% get a leaf or stick
            imsuper_1d=imsuper(:); %insert, using mask
            imsuper_1d(mlocs)=leaf_stick(:);
            imsuper=reshape(imsuper_1d,ard_super);
        end %set_length>0
    end
end
%
optsused.flfs_stats.misses=misses;
optsused.flfs_stats.dropped=dropped;
optsused.flfs_stats.dropped_cat=dropped_cat;
optsused.flfs_stats.cat_list=cat_list;
optsused.flfs_stats.coverage_pointers=coverage_pointers;
optsused.flfs_stats.coverage_depth=coverage_depth;
optsused.flfs_stats.coverage_coords=coverage_coords;
optsused.flfs_stats.imsuper=imsuper;
if (opts.show>1)
   figure;
   set(gcf,'Position',[100 100 1200 800]);
   set(gcf,'NumberTitle','off');
   set(gcf,'Name','FLFS internals');
   %
   subplot(2,3,1);
   imagesc(coverage_depth);
   axis equal;axis off;
   colormap('jet');
   colorbar;
   title('coverage depth');
   %
   subplot(2,3,4);
   hist(coverage_depth(:));
   title('coverage depth');
   %
   subplot(2,3,2);
   imagesc(coverage_pointers);
   axis equal;axis off;
   colormap('jet');
   colorbar;
   title('coverage pointers');
   %
   coverage_pointers_ord=zeros(size(coverage_pointers));
   for mord=1:length(dropped_unique)
       coverage_pointers_ord(find(coverage_pointers==dropped_unique(mord)))=mord;
   end
   subplot(2,3,5);
   imagesc(coverage_pointers_ord);
   axis equal;axis off;
   colormap('jet');
   colorbar;
   title('coverage pointers, ordinal');
   %
   subplot(2,3,3);
   imagesc(coverage_types(:,:,1),[0 2]);
   axis equal;axis off;
   colormap('jet');
   title('type 1 segments');
   subplot(2,3,6);
   imagesc(coverage_types(:,:,2),[0 2]);
   axis equal;axis off;
   colormap('jet');
   title('type 2 segments');
end
%
% apply Metropolis algorithm
%
if opts.metro_show>0
    if opts.metro_show>1
        disp(sprintf('in mtc_flfs_makemaps ready for Metropolis algorithm, %s, with metro_opts:','FLFS'));
        disp(opts.metro_opts)
        disp(sprintf('map size: %5.0f x %5.0f',size(imsuper)));
    else
        disp(sprintf('in mtc_flfs_makemaps ready for Metropolis algorithm, %s.','FLFS'));
    end
end
donut=[];
donut.name='donut';
donut.matrix=[1 1 1;1 0 1; 1 1 1];
donut=glider_addcoords(donut);
if opts.onaxis==1
    disp('Metropolis algorithm not required, onaxis=1');
else
    if ~opts.metro_opts.numiters==0
        disp(sprintf('Calling Metropolis algorithm donut_metro with numiters=%7.0f',opts.metro_opts.numiters));
        [imsuper_new,metro_samp,metro_optsused]=donut_metro(ng,donut,imsuper,setfield(opts.metro_opts,'map_show_name','FLFS'));
        metro=[];
        metro.samp=metro_samp;
        metro.optsused=metro_optsused;
        metro.frac_changed=mean(abs(imsuper_new(:)-imsuper(:)));
        metro.numiters=opts.metro_opts.numiters;
        if opts.metro_show>0
            disp(sprintf('iterations: %5.0f, fraction of checks changed: %8.5f',...
                opts.metro_opts.numiters,metro.frac_changed));
        end
        optsused.flfs_stats.imsuper_new=imsuper_new;
        imsuper=imsuper_new;
    else
        disp('Metropolis algorithm not called, numiters=0');
    end
end
%cut up and put results into img
for imap=1:opts.nmaps
    ic=mod(imap-1,nc)+1;
    ir=floor((imap-1)/nc)+1;
    img(:,:,imap)=imsuper(opts.burnin+ard(1)*(ir-1)+[1:ard(1)],opts.burnin+ard(2)*(ic-1)+[1:ard(2)]);
end
%
%plot if requested
if (opts.show>0)
   figure;
   imagesc(img(:,:,1),[0 ng-1]);
   axis equal;axis off;colormap('gray');
end
return
