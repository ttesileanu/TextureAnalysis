function [mtcs,dict_used,opts_used]=mtc_define(ng,dict,opts)
% [mtcs,dict_used,opts_used]=mtc_define(ng,dict,opts) is the first step in setting up
% coordinates for 2x2 textures with multiple gray levels
%
% ng: number of gray levels
% dict: (optional) output of btc_define or gtc_define.
% opts: (optional) options
%    opts.iflog=1 to log verifications (defaults to 0)
%
% mtcs: structure defining the local texture coordinates and grouping them by linear combinations
%   key variables:
%      mtcs.ng:      number of gray levels
%      mtcs.nchecks: number of checks 
%      mtcs.lincomb_count: number of linear combinations (ng^nchecks)
%      mtcs.lincomb_list: list of the linear combinations assigned to each check, by check number
%      mtcs.coord_count: number of coordinates (taking into account translation constraints)
%      mtcs.coord_indices:  the linear combination assigned to each coordinate, as index into lincomb_list
%      mtcs.lincomb_to_coord:  the coordinate assigned to each linear combination, as index into [1:mtcs.coord_count]
%      mtcs.coord_names: coordinate names, as strings
%      mtcs.coord_groups_count: number of coord groups
%      mtcs.coord_groups: groups the coordinates into disjoint subsets if their linear combinations are related by a scalar
% dict_used: dict used, or as supplied
% opts_used: options used
%
% in contrast with GTC_DEFINE, this generalizes to multiple gray levels, rather than other shapes
%
% See also:  BTC_DEFINE, TRICO_CVT, GTC_DEFINE.
%
if (nargin<=1)
    dict=[];
end
if isempty(dict)
    dict=btc_define;
end
dict_used=dict;
if (nargin<=2)
    opts=[];   
end
opts=filldefault(opts,'iflog',0);
opts_used=opts;
%
if isfield(dict,'boundingbox') %if dict comes from gtc_define
    region_extent=max(dict.boundingbox,[],1);
else
    region_extent=[2 2];
end
if ~isfield(dict,'nchecks')
    dict.nchecks=prod(region_extent);
end
mtcs=dict;
mtcs.ng=ng;
mtcs.region_extent=region_extent;
nchecks=dict.nchecks;
checkdef=dict.checkdef;
nbtc=length(dict.checks);
if ~isfield(dict,'trans_vecs') %if dict comes from btc_define
    for ibtc=1:nbtc
        mults=dict.mults(:,ibtc);
        mtcs.trans_vecs{ibtc}=zeros(2,0);
        for irow=0:mults(1)-1
            for icol=0:mults(2)-1
               mtcs.trans_vecs{ibtc}(:,end+1)=[irow;icol];
            end
        end
    end
end
%
letters=fieldnames(mtcs.checkdef);
if ~isfield(mtcs,'region_mask') | ~isfield(mtcs,'region_mask_letters') %these are defined in gtc_define
    mtcs.region_mask=zeros(mtcs.region_extent);
    mtcs.region_mask_letters=cell(mtcs.region_extent);
    for icheck=1:mtcs.nchecks
        checkcoords=mtcs.checkdef.(letters{icheck});
        mtcs.region_mask(1+checkcoords(1),1+checkcoords(2))=1;
        mtcs.region_mask_letters{1+checkcoords(1),1+checkcoords(2)}=letters{icheck};
    end  
end
%
mtcs.v=[];
%
mtcs.checkdef_list=zeros(nchecks,2);
checkdef_fields=fields(dict.checkdef);
for ifield=1:length(checkdef_fields)
    mtcs.checkdef_list(ifield,:)=dict.checkdef.(checkdef_fields{ifield});
end
%these are the different subsets of checks, without taking into account matches due to translation
mtcs.subset_list=int2nary([0:2^nchecks-1]',2,nchecks);
nsubsets=size(mtcs.subset_list,1);
mtcs.subset_index=zeros(nsubsets,1);
mtcs.subset_index_whichtrans=zeros(nsubsets,1);
%establish lists of which checks define the correlation, and of their translates
ntrans_tot=0;
ntrans_by_order=zeros(1,nchecks);
for ibtc=1:nbtc
    s.order=dict.order(ibtc);
    s.check_letters=dict.checks{ibtc};
    s.check_vecs=[];
    s.check_nums=[];
    s.check_nums_trans=[];
    for icheck=1:s.order
        s.check_vecs(icheck,:)=dict.checkdef.(s.check_letters(icheck)); %a list of the defining checks, as (row,column) offsets
        s.check_nums(1,icheck)=find(all(repmat(s.check_vecs(icheck,:),nchecks,1)==mtcs.checkdef_list,2)==1); %a list of the defining checks, as indexes into checkdef_list
    end
    ntrans=size(mtcs.trans_vecs{ibtc},2);
    for itransvec=1:ntrans
        translated=s.check_vecs+repmat(mtcs.trans_vecs{ibtc}(:,itransvec)',s.order,1);
        for icheck=1:s.order
            s.check_nums_trans(itransvec,icheck)=find(all(repmat(translated(icheck,:),nchecks,1)==mtcs.checkdef_list,2)==1);
        end
        checks_used=ismember([1:nchecks],s.check_nums_trans(itransvec,:));
        typeno=find(all(repmat(checks_used,nsubsets,1)==mtcs.subset_list,2)==1);
        mtcs.subset_index(typeno)=ibtc;
        mtcs.subset_index_whichtrans(typeno)=itransvec;
    end
    %
    ntrans_tot=ntrans_tot+ntrans;
    ntrans_by_order(s.order)=ntrans_by_order(s.order)+ntrans;
    %
    if isfield(dict,'codel') %if there is a short name
        s.codel=dict.codel(ibtc);
    end
    if isfield(dict,'name_order')
        s.name_order=dict.name_order{ibtc};
    end
    if isfield(dict,'name_order_aug')
        s.name_order_aug=dict.name_order_aug{ibtc};
    end
    s.lincomb=zeros(0,nchecks);
    mtcs.v.(dict.checks{ibtc})=s;
end
%
%check that each subset except the null subseet is used
%
mtcs.unused_subsets=length(find(mtcs.subset_index==0));
mtcs.ntrans_tot=ntrans_tot;
mtcs.ntrans_by_order=ntrans_by_order;
if ~(mtcs.unused_subsets==1) 
    warning(sprintf('number of unused subsets should be 1, but is %5.0f',mtcs.unused_subsets));
else
    if (opts.iflog==1)
        disp(sprintf('mtc_define check: number of unused subsets should be 1, is %5.0f',mtcs.unused_subsets));
    end
end
if ~(ntrans_tot==2^nchecks-1)
    warning(sprintf('number of total placements including translation should be %5.0f but is %5.0f',2^nchecks-1,ntrans_tot));
    disp(' by order:')
    disp(ntrans_by_order);
else
    if (opts.iflog==1)
        disp(sprintf('mtc_define check: number of total placements including translation should be %5.0f, is %5.0f',2^nchecks-1,ntrans_tot));
    end
end
%
% now generate all possible tuples of [0:ng-1] and sort them into the texture family
%
% instead of mtcs.lincomb_list=int2nary([0:ng^nchecks-1]',ng,nchecks);
% use first an ordering by order (number of checks involved) and inside that, a lexicographical order
%
lincomb_list_raw=int2nary([0:ng^nchecks-1]',ng,nchecks);
lincomb_orders=sum((lincomb_list_raw>0),2);
lincomb_orders_sorted=sortrows([lincomb_orders,(lincomb_list_raw>0),lincomb_list_raw]);
mtcs.lincomb_list=fliplr(lincomb_orders_sorted(:,(2+nchecks):end));
mtcs.lincomb_count=size(mtcs.lincomb_list,1);
%
lincomb_nz=double(mtcs.lincomb_list>0);
nlincomb=size(mtcs.lincomb_list,1);
mtcs.lincomb_subset=zeros(nlincomb,1); %which subset of checks defines this linear combination, after translation
mtcs.lincomb_trans=zeros(nlincomb,1); %which translation (index into mtcs.v.{checkletters}.check_nums_trans)
%
%what are the allowed multipliers for a linear combination?
%poss_factors=[2:floor(ng/2)]; %factors not including 1 or ng
mtcs.ng_factors=unique(factor(ng));
%
mtcs.disallowed_mults=zeros(nlincomb,length(mtcs.ng_factors));
%this only is executed if ng is composite
for ifactor=1:length(mtcs.ng_factors)
    %if multiplying by the factor increases the gcd of any entry in the linear combination, it can't be used
    %this will happen if at least one member of lincomb_list is not a multiple of the highest prime power that divides ng
    mtcs.disallowed_mults(:,ifactor)=any(gcd(mtcs.lincomb_list*mtcs.ng_factors(ifactor),ng)>gcd(mtcs.lincomb_list,ng),2);
end
%
% now determine the coordinates
%  coordinates are the linear combinations with mtcs_lincomb_trans=1
%  they may not be listed in the same order as mtcs.checks
%
mtcs.lincomb_mult=zeros(nlincomb,1); %the multiplier for the first nonzero check
mtcs.lincomb_baseharm=zeros(nlincomb,1); %the base harmonic, the gcd of all of the multipliers in the linear comb and ng 
for ilincomb=2:nlincomb % skip the first entry, [0 0 0 ... 0]
    sub_index=find(all(repmat(lincomb_nz(ilincomb,:),nsubsets,1)==mtcs.subset_list,2)==1); %which subset of checks is used
    ibtc=mtcs.subset_index(sub_index);
    mtcs.lincomb_subset(ilincomb)=ibtc;
    mtcs.lincomb_trans(ilincomb)=mtcs.subset_index_whichtrans(sub_index);
    mtcs.lincomb_mult(ilincomb)=mtcs.lincomb_list(ilincomb,min(find(mtcs.lincomb_list(ilincomb,:)>0))); %first nonzero value
    mtcs.v.(mtcs.checks{ibtc}).lincomb(end+1,:)=mtcs.lincomb_list(ilincomb,:);
    igcd=ng;
    for imult=find(mtcs.lincomb_list(ilincomb,:)>0);
        igcd=gcd(igcd,mtcs.lincomb_list(ilincomb,imult));
    end
    mtcs.lincomb_baseharm(ilincomb)=igcd;
end
mtcs.lincomb_baseharm(1)=ng; %by definition
mtcs.lincomb_period=ng./mtcs.lincomb_baseharm;
mtcs.allowed_divisors=zeros(nlincomb,length(mtcs.ng_factors));
%
%create the generalized totient function values, see graylevs_mtc.doc
% slow but clear
%
periods=unique(mtcs.lincomb_period);
mtcs.period_list=periods;
mtcs.totients=cell(1,max(periods));
mtcs.totients_spaced=zeros(ng,length(periods));
for iperiod=1:length(periods)
    d=periods(iperiod);
    for b=0:d-1
        mtcs.totients{d}(b+1)=totient(d,b);
    end
    mtcs.totients_spaced(1+(ng/d)*[0:d-1],iperiod)=mtcs.totients{d};
end
%this only is executed if ng is composite
for ifactor=1:length(mtcs.ng_factors)
    %see graylev_mtc.doc for this.  But it should be equivalent to 1-mtcs.disallowed_mults
    mtcs.allowed_divisors(2:end,ifactor)=double(mod(ng./mtcs.lincomb_baseharm(2:end),mtcs.ng_factors(ifactor))>0);
end
mtcs.allowed_divisors(1,:)=1; %all divisors allowed for 0
if ~(all(mtcs.disallowed_mults+mtcs.allowed_divisors==1))
    warning(sprintf('alternative computations of disallowed multipliers and allowed divisors are discrepant.'))   
else
    if (opts.iflog==1)
        disp(sprintf('mtc_define check: disallowed multipliers and allowed divisors agree.'));
    end
end
%
mtcs.unused_lincomb=length(find(mtcs.lincomb_subset==0));
if ~(mtcs.unused_lincomb==1) 
    warning(sprintf('number of unused linear combinations should be 1, but is %5.0f',mtcs.unused_lincomb));
else
    if (opts.iflog==1)
        disp(sprintf('mtc_define check: number of unused linear combinations should be 1, is %5.0f',mtcs.unused_lincomb));
    end
end
%
mtcs.coord_indices=find(mtcs.lincomb_trans==1);
mtcs.coord_count=length(mtcs.coord_indices);
for imult=1:ng-1
    mtcs.coord_multcount(imult)=sum(mtcs.lincomb_mult(mtcs.coord_indices)==imult);
end
if ~(min(mtcs.coord_multcount)==max(mtcs.coord_multcount))
    warning('number of coordinates with each multiplier must be equal.  Instead, they are')
    disp(mtcs.coord_multcount)
else
    if (opts.iflog==1)
        disp(sprintf('mtc_define check: number coordinates with each multiplier must be equal, and they are equal, to %5.0f',min(mtcs.coord_multcount)));
    end
end
%
%determine mapping of linear combinations to coords by translating the
%linear combination back to primary position
%
mtcs.lincomb_to_coord=zeros(nlincomb,1);
for ilincomb=2:nlincomb
    whichtrans=mtcs.lincomb_trans(ilincomb);
    ibtc=mtcs.lincomb_subset(ilincomb);
    checks=mtcs.checks{ibtc};
    lincomb_orig=mtcs.lincomb_list(ilincomb,:);
    checks_orig=mtcs.v.(checks).check_nums_trans(whichtrans,:);
    checks_primary=mtcs.v.(checks).check_nums_trans(1,:);
    lincomb_primary=zeros(1,nchecks);
    lincomb_primary(checks_primary)=lincomb_orig(checks_orig);
    mtcs.lincomb_to_coord(ilincomb)=find(all(repmat(lincomb_primary,mtcs.coord_count,1)==mtcs.lincomb_list(mtcs.coord_indices,:),2));
end
if ~all(mtcs.lincomb_to_coord(mtcs.coord_indices)==[1:mtcs.coord_count]')
    warning('discrepancy between mapping from linear combinatinos to coordinates, and coordinates to linear combinations.');
else
    if (opts.iflog==1)
        disp('mtc_define check: consistent mapping from linear combinatinos to coordinates, and coordinates to linear combinations.');
    end
end
%create coordinate names, one for each subset of linear combinations that differ by a translation
mtcs.coord_names=cell(1,mtcs.coord_count);
mtcs.coord_subset=mtcs.lincomb_subset(mtcs.coord_indices);
for icoord=1:mtcs.coord_count
    coord_index=mtcs.coord_indices(icoord);
    mtcs.lincomb_subset(coord_index);
    lincomb=mtcs.lincomb_list(coord_index,:);
    lincomb_nz=lincomb(find(lincomb>0));
    lincomb_string=[];
    %insert an underscore before each multiplier, to allow for multidigit multipliers
    for inz=1:length(lincomb_nz)
        lincomb_string=cat(2,lincomb_string,sprintf('_%1.0f',lincomb_nz(inz)));
    end
    mtcs.coord_names{icoord}=cat(2,dict.checks{mtcs.lincomb_subset(coord_index)},lincomb_string);
end
%now group the linear combinations that differ only by a multiplier
coord_group_index=zeros(mtcs.coord_count,1);
coord_groups=cell(0);
coord_groups_multlist=zeros(0,ng-1);
ngroups=0;
lincomb_coords=mtcs.lincomb_list(mtcs.coord_indices,:);
coord_group_ptrs=struct();
while any(coord_group_index==0)
    ngroups=ngroups+1;
    first=min(find(coord_group_index==0));
    lincomb=lincomb_coords(first,:);
    matches=[]; %keeps track of unique matches, and no entry if no match
    lincombs=[];
    mults=[1:ng-1];
    disallowed_mults=find(mtcs.disallowed_mults(mtcs.coord_indices(first),:)>0);
    for ifactor=1:length(disallowed_mults)
        prime_factor=mtcs.ng_factors(disallowed_mults(ifactor));
        mults=setdiff(mults,prime_factor*[1:ng/prime_factor]);
    end
    for imult_ptr=1:length(mults)
        imult=mults(imult_ptr);
        lincomb_m=mod(lincomb*imult,ng);
        lincombs=[lincombs;lincomb_m];
        newmatch=find(all(repmat(lincomb_m,mtcs.coord_count,1)==lincomb_coords,2));
        if ~isempty(newmatch) & coord_group_index(newmatch)==0
            matches=[matches newmatch];
        end
    end
    colsort=min(find(min(lincombs,[],1)>0)); %column to sort on
    if (isempty(colsort)) 
        ngroups
        disallowed_mults
        mults
        lincombs
        colsort
    else
        %find the unique linear combinations generated by multiplication and sort by colsort
        %if colsort values differ, then the rows must differ too
        %lincombs
        [colsorted,ci,cj]=unique(lincombs(:,colsort));
        %ci
        %cj
    end   
    coord_groups{ngroups}.name_firstcoord=mtcs.coord_names{first}; %first coordinate in lincomb
    coord_groups{ngroups}.name=mtcs.coord_names{matches(ci(1))}; %coordinate in lincomb with lowest first coefficient
    coord_group_index(matches)=ngroups;
    coord_groups{ngroups}.coord_num_allmults=matches; %each coordinate number used only once
    coord_groups{ngroups}.baseharm=mtcs.lincomb_baseharm(mtcs.coord_indices(first));
    coord_group_ptrs.(coord_groups{ngroups}.name)=ngroups;
    %create a list of which coordinates correspond to the multiples of this linear combination
    matches_bymult=zeros(1,ng-1); %keeps track of all matches, in order of the multipliers, and has a 0 entry if no match
    lincomb=lincomb_coords(matches(ci(1)),:);
    numlist=[];
    for imult_ptr=1:length(mults) %redo calculation similar to above but using the lincomb with lowest first coefficient as starting point
        imult=mults(imult_ptr);
        lincomb_m=mod(lincomb*imult,ng);
        newmatch=find(all(repmat(lincomb_m,mtcs.coord_count,1)==lincomb_coords,2));
        if ~isempty(newmatch)
            matches_bymult(imult)=newmatch;
            if ~ismember(newmatch,numlist)
                numlist=[numlist newmatch];
            end
        end
    end
    coord_groups_multlist(ngroups,:)=matches_bymult; %coordinate numbers may be repeated,and 0's may occur
    coord_groups{ngroups}.coord_num=numlist; %each coordinate number used only once, beginning with lowest linear comb
end
mtcs.coord_group_ptrs=coord_group_ptrs;
mtcs.coord_group_index=coord_group_index;
mtcs.coord_groups=coord_groups;
mtcs.coord_groups_count=length(coord_groups);
mtcs.coord_groups_multlist=coord_groups_multlist;
if (isprime(ng))
    if ~(length(mtcs.coord_groups)==size(mtcs.coord_groups_multlist,1))
        warning(sprintf('Since ng (=%2.0f) is prime, all coord groups should have same length, and this fails',ng));
    else
        if (opts.iflog==1)
            disp(sprintf('mtc_define check: ng (=%2.0f) is prime, all coord groups should have same length, %5.0f, and this holds',ng,size(mtcs.coord_groups_multlist,2)));
        end
    end
end
%verify that all coordinate numbers occur exactly once in the coordinate groups
coord_num_list=[];
for igroup=1:length(mtcs.coord_groups)
    coord_num_list=[coord_num_list,mtcs.coord_groups{igroup}.coord_num];
end
if ~(length(coord_num_list)==mtcs.coord_count)
    warning(sprintf('total number of coordinates in groups is %5.0f, should be %5.0f',length(coord_num_list),mtcs.coord_count))
else
    if (opts.iflog==1)
        disp(sprintf('mtc_define check: total number of coordinates in groups is %5.0f, and is %5.0f',length(coord_num_list),mtcs.coord_count))
    end
end
coord_num_list=unique(coord_num_list);
if  ~((coord_num_list(1)==1) & (coord_num_list(end)==mtcs.coord_count))
    warning('not every coordinate number appears in a group.  Omitted coordinate numbers:')
    disp(setdiff([1:mtcs.coord_count],coord_num_list));
else
    if (opts.iflog==1)
        disp('mtc_define_check: every coordinate number appears in a group')
    end
end
if (opts.iflog==1)
    disp(sprintf('mtc_define check: %2.0f gray levels, %2.0f checks, %5.0f check subsets, %5.0f linear combs, %5.0f  coordinates, %5.0f  coordinate groups',...
        ng,mtcs.nchecks,size(mtcs.subset_list,1),size(mtcs.lincomb_list,1),mtcs.coord_count,mtcs.coord_groups_count));
end
return

