function [dict,opts_used]=gtc_define(region,opts)
% [dict,opts_used]=gtc_define(region,opts) defines a structure "dict" that specifies the names and conventions
% for the correlation parameters on a general region
%
%  derived from btc_define
%
% to add:
%   convenient definition of region via strings
%   cross-referencing to btc (i.e., which coords are already in btc)
%   tag and sort subsets by symmetry
% done:
%   sort the order of subsets within each order by size of bounding box
%
% region:
%   defines the region shape.  region(ic,1), region(ic,2) are the x and y coordinates
%   of each check.  x is row, where 0 is top row. y is column, where 1 is right column
%   For a standard 2x2 region, region=[0 0;0 1;1 0;1 1];
% opts: 
%   ifshow:  show a table
%   ifshow_fig: show a figure
%   ordernames: cell array of order names, typically
%   {'gamma','beta','theta','alpha'}, for the coordinates found in btc_define
%   
% checks are (note anomaly in upper left for consistency with btc_define
%   A B E F G 
%   C D H I J
%   K L M N O
%   P Q R S T
%   U V W X Y
%   
%  dict.order(k) is the order (number of checks) in the kth parameterm
%  dict.checkcoords{k}.z is a  dict.order(k) x 2 array of coordinates of the checks that constitute a paramater
%
%    See also: BTC_DEFINE, INT2NARY.
%
regmax=[5 5]; %maximum region
%define check letters, with ABCD anomaly in upper left for consistency with btc_define
check_letters=cell(5,5);
for ix=1:2
    for iy=1:2
        check_letters{ix,iy}=char(double('A')+(iy-1)+2*(ix-1));
    end
    for iy=3:5
        check_letters{ix,iy}=char(double('A')+(iy+1)+3*(ix-1));
    end
end
for ix=3:5
    for iy=1:5
        check_letters{ix,iy}=char(double('A')+(iy-1)+5*(ix-1));
    end
end
xy={'x','y'};
%
if (nargin<2)
    opts=[];
end
opts=filldefault(opts,'ifshow',0);
opts=filldefault(opts,'ifshow_fig',0);
opts=filldefault(opts,'ordernames',{'gamma','beta','theta','alpha'});
dict=[];
opts_used=opts;
%
%check that region is properly defined
%
for ixy=1:2
    if min(region(:,ixy))~=0
        error(sprintf('region %s-values do not start at 0',xy{ixy}));
    end
    if max(region(:,ixy))>=regmax(ixy)
        error(sprintf('region %s-values equal or exceed  %1.0f',xy{ixy},regmax(ixy)));
    end
end
if size(unique(region,'rows'),1)<size(region,1)
    error('region definition contains duplicates.');
end
region=sortrows(region);
dict.nchecks=size(region,1);
dict.region=region;
%
for ic=1:dict.nchecks
    dict.checkdef.(check_letters{region(ic,1)+1,region(ic,2)+1})=region(ic,:);
    dict.region_mask(region(ic,1)+1,region(ic,2)+1)=1;
end
dict.region_mask_letters=cell(size(dict.region_mask));
for ic=1:dict.nchecks
    dict.region_mask_letters{region(ic,1)+1,region(ic,2)+1}=...
        check_letters{region(ic,1)+1,region(ic,2)+1};
end
%
%determine list of binary coordinates by seeing which subsets of region are
%equivalent up to translation, keep this in trans_vecs, and keep total
%multiplicity in mult  (dict.mults does not generalize)
%
% for each parameter:
% dict.order(k)
% dict.checkcoords{k}
% dict.checks{k}
% dict.mults(k): number of placements
% dict.trans_vecs{k}: mults(k) x 2 set of translation vectors
%
%one way to retrieve some region definitions
% gliderlib=glider_choose;
% gliderlib=glider_choose;
% gliderlib
% gliderlib = 
%        box: [2x2 double]
%        tee: [2x3 double]
%     tri_lr: [2x2 double]
% region=getfield(glider_addcoords(setfield([],'matrix',gliderlib.box)),'inds')-1
%
%[img,ga]=eo_glider('even',1);
% region=[0 0;ga]
%     can use glider='even','cross','zigzag','oblong',
%        'triangle','tee','wye','foot','el'

dict.order=[];
dict.checkcoords=cell(0);
dict.checks=cell(0);
dict.mults=[];
dict.trans_vecs=cell(0);
dict.boundingbox=zeros(0,2);
nary=int2nary([0:2^dict.nchecks-1]',2,dict.nchecks);
for iord=1:dict.nchecks %do this order by order
    nary_ord=nary(find(sum(nary,2)==iord),:);
    checksels=nchoosek([1:dict.nchecks],iord);
    %set up table of coords: cxy(:,:,1) has x-coord, cxy(:,:,2) has y-coord
    cxy=zeros([size(checksels,1) size(checksels,2),2]);
    for ixy=1:2
        cxy(:,:,ixy)=reshape(region(checksels,ixy),size(checksels));
    end
    cxyrel=cxy-repmat(cxy(:,1,:),[1 iord 1]);
    cxyrel=[cxyrel(:,:,1),cxyrel(:,:,2)];
    %code to ensure that unique_ptr points to first row that matches 
    %old versions of matlab don't allow a third argument for unique of 'first') 
    [unique_rows_flipped,unique_ptr_flipped,unique_which_flipped]=unique(flipud(cxyrel),'rows');
    [unique_rows,unique_ptr_old,unique_which]=unique(cxyrel,'rows');
    unique_ptr=size(cxyrel,1)+1-unique_ptr_flipped;
    for icoord=1:size(unique_rows,1)
        dict.order(end+1)=iord;
        urow=cxy(unique_ptr(icoord),:); %find the first set of absolute coords that led to this set of relative coords
        urow=sortrows(reshape(urow,[iord,2]));
        %translate it so that it is in the region            
        dict.checkcoords{end+1}=urow;
        dict.boundingbox(end+1,:)=max(urow,[],1)-min(urow,[],1)+1;
        %
        copies=find(unique_which==icoord);
        mults=sum(unique_which==icoord);
        dict.mults(end+1)=mults;
        first=min(find(unique_which==icoord)); %first example of this subset
        %
        which_letters=[];
        for icheck=1:iord
            which_letters(1,icheck)=check_letters{cxy(first,icheck,1)+1,cxy(first,icheck,2)+1};
        end
        dict.checks{end+1}=char(which_letters);
        %
        trans_vecs=[];
        for imult=1:mults
            trans_vecs(:,imult)=squeeze(cxy(copies(imult),1,:));
        end
        trans_vecs=trans_vecs-repmat(trans_vecs(:,1),[1 mults]); %translation is always relative to first example
        dict.trans_vecs{end+1}=trans_vecs;
    end
end
sortorder=[dict.order',prod(dict.boundingbox,2)];
[so,sptr]=sortrows(sortorder);
duns=dict;
for ic=1:length(sptr)
    dict.checkcoords{ic}=duns.checkcoords{sptr(ic)};
    dict.checks{ic}=duns.checks{sptr(ic)};
    dict.trans_vecs{ic}=duns.trans_vecs{sptr(ic)};
end
dict.order=duns.order(sptr); %should not change
dict.boundingbox=duns.boundingbox(sptr,:);
dict.mults=duns.mults(sptr);

%for each subregion, find the one that is maximally pushed to upper left
% and find the translation vectors to it
return
 
%fields below this point are from btc_define
dict.checks{1}=['A'];
dict.checks{2}=['A','B'];
dict.checks{3}=['A','C'];
dict.checks{4}=['A','D'];
dict.checks{5}=['B','C'];
dict.checks{6}=['B','C','D'];
dict.checks{7}=['A','C','D'];
dict.checks{8}=['A','B','C'];
dict.checks{9}=['A','B','D'];
dict.checks{10}=['A','B','C','D'];
%
dict.codel=['g','b','c','d','e','t','u','v','w','a'];%  single-letter code
%dict.order=[ 1   2   2   2   2   3   3   3   3   4]; % order of correlation
dict.posit=[ 1   1   2   3   4   4   3   1   2   1]; %subscript position in corrs
dict.rot=  [ 0   0   1   0   1   0   1   2   3   0]; %number of 90-deg rotations w.r.t. standard   
dict.mults=ones(2,length(dict.checks));
dict.inpickard=zeros(2,length(dict.checks));
for k=1:length(dict.checks)
    dict.checkcoords{k}=[];
    dict.order(1,k)=length(dict.checks{k});
    dict.name_order{1,k}=dict.ordernames{dict.order(k)};
    dict.name_full{1,k}=sprintf('%s(%1.0f)',dict.name_order{1,k},dict.posit(1,k));
    if (ismember(dict.order(k),[2 3]))
        dict.name{1,k}=dict.name_full{1,k};
    else
        dict.name{1,k}=dict.name_order{1,k};
    end
    for i=1:length(dict.checks{k})
%        dict.checkcoords{k}(i,:)=getfield(dict.checkdef,dict.checks{k}(i));
        dict.checkcoords{k}(i,:)=dict.checkdef.(dict.checks{k}(i));
    end
    % do multiplicities and augmented names (beta_horiz, beta_vert)
    dict.name_order_aug{1,k}=dict.name_order{1,k};
    if (dict.order(1,k)==2)
        if all(sum(dict.checkcoords{k},1)==1)
            dict.name_order_aug{1,k}=cat(2,dict.name_order_aug{1,k},'_diag');
        else
            dict.name_order_aug{1,k}=cat(2,dict.name_order_aug{1,k},'_hv');
        end
    end
    %multiplicity of possible placements along each axis
    dict.mults(:,k)=(2-(max(dict.checkcoords{k},[],1)-min(dict.checkcoords{k},[],1)))';
    dict.mult(1,k)=prod(dict.mults(:,k));
    %do involvements in Pickard
    if (dict.mult(1,k)>1)
        dict.inpickard(:,k)=1;
    else
        if (ismember(dict.order(1,k),[2 3]))
            if ismember('A',dict.checks{k}) & ismember('D',dict.checks{k})
                dict.inpickard(1,k)=1;
            end
            if ismember('B',dict.checks{k}) & ismember('C',dict.checks{k})
                dict.inpickard(2,k)=1;
            end
        end
    end
end
dict.name_order_aug_unique=cellstr(unique(dict.name_order_aug))'; %added 21 Dec 2010 since this is slow to recompute
%
if (opts.ifshow)
    disp(' uid letter rot name_full  name_order      name   checks used name_order_aug  mults    inpickard')
    for k=1:length(dict.checks)
        lstar=' ';
        if (dict.rot(k)==0)
            lstar='*';
        end
        disp(sprintf(' %2.0f    %1s%1s    %1.0f    %8s   %8s   %8s     %4s     %12s  %3.0f %3.0f %3.0f  %1.0f  %1.0f',...
            k,dict.codel(k),lstar,dict.rot(k),dict.name_full{k},dict.name_order{k},dict.name{k},dict.checks{k},dict.name_order_aug{k},...
            dict.mults(:,k),dict.mult(k),dict.inpickard(:,k)));
    end
end
if (opts.ifshow_fig)
    figure;
    set(gcf,'Position',[200 200 1200 500]);
    for k=1:length(dict.checks)
        subplot(2,5,k);
        map=repmat(0.5,4,4);
        map(2:3,2:3)=0;
        for i=1:size(dict.checkcoords{k},1);
            map(2+dict.checkcoords{k}(i,1),2+dict.checkcoords{k}(i,2))=1;
        end
        imagesc(map,[0 1]);
        axis square;
        title(sprintf('param %2.0f: %s=%s',k,dict.codel(k),dict.name{k}));
        axis off;
        colormap gray;
    end
end
opts_used=opts;
return

