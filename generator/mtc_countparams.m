function counts=mtc_countparams(mtcs,opts) 
% counts=mtc_countparams(mtcs) sumarizes and counts the parameters in mtcs
%
% mtcs:  a structure produced by mtc_define
% opts: iflog=1 to log, defaults to 0
%
% counts: a structure that counts the number of coordinate groups and number of coordinates, by order
%    counts.coord_group_list{corder}:  a list of coordinate groups for each order
%
%  See also:  MTC_DEFINE, MTC_FREEPARAMS.
%
if (nargin<=1) opts=[]; end
opts=filldefault(opts,'iflog',0);
ng=mtcs.ng;
nchecks=mtcs.nchecks;
coord_count=mtcs.coord_count;
coord_groups_count=mtcs.coord_groups_count;
coord_byorder=zeros(1,nchecks);
coord_groups_byorder=zeros(1,nchecks);
coord_byorder_verify=zeros(1,nchecks);
if (opts.iflog==1)
    disp(sprintf(' number of gray levels: %4.0f',ng));
    disp(sprintf(' number of      checks: %4.0f',nchecks));
end
%count coords by group
counts.coord_group_list=cell(1,nchecks);
counts.coord_list=cell(1,nchecks);
nfreeparams_byorder=zeros(1,nchecks);
for ig=1:coord_groups_count
    cg=mtcs.coord_groups{ig};
    coord_nums=cg.coord_num;
    ncoords=length(coord_nums); %number of coordinates that map into this coordinate group
    lincomb_indices=mtcs.coord_indices(coord_nums);
    corder=sum(mtcs.lincomb_list(lincomb_indices(1),:)>0);
    coord_groups_byorder(corder)=coord_groups_byorder(corder)+1;
    coord_byorder(corder)=coord_byorder(corder)+ncoords;
    counts.coord_group_list{corder}(end+1)=ig; %which groups have this order
    counts.coord_list{corder}(end+[1:ncoords])=coord_nums;
    %
    multlist=mtcs.coord_groups_multlist(ig,:);
    multlist=multlist(find(multlist>0));
    nfreeparams_byorder(corder)=nfreeparams_byorder(corder)+length(unique(multlist));
end
%count coords one by one
lincomb_indices_all=mtcs.coord_indices;
orders_all=sum(mtcs.lincomb_list(lincomb_indices_all,:)>0,2);
for iorder=1:nchecks
    coord_byorder_verify(iorder)=sum(orders_all==iorder);
end
counts.coord_byorder=coord_byorder;
counts.coord_byorder_verify=coord_byorder_verify;
counts.coord_groups_byorder=coord_groups_byorder;
if (opts.iflog==1)
    disp('         coord_groups   freeparams coords     verify')
    for iorder=1:nchecks
        disp(sprintf(' order %2.0f:  %7.0f    %7.0f    %7.0f    %7.0f',...
            iorder,coord_groups_byorder(iorder),nfreeparams_byorder(iorder),coord_byorder(iorder),coord_byorder_verify(iorder)))
    end
    disp(sprintf('    total:  %7.0f               %7.0f    %7.0f',sum(coord_groups_byorder),sum(coord_byorder),sum(coord_byorder_verify)));
    disp(sprintf('   verify:  %7.0f               %7.0f',coord_groups_count,coord_count));
end
return

