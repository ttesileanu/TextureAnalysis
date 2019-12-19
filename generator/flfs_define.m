function flfs=flfs_define(cats,mtcs,opts)
%flfs=flfs_define(cats,mtcs,opts) sets up a structure that defines a falling-leaves, falling-stick structure
%
% cats: a list of categories for falling-leaves or falling-sticks
%   each is a glider shape that can fall onto the glider defined by mtcs;
%      in principle, they can be duplicated, so that each contains different kinds of correlation structure 
%   cats{icat}.type={'row','col'} (only options allowed at present, but in future could be finite regions)
%   if empty, it is assumed that cats{1}='row',cats{2}='col'
% mtcs: a glider structure, typically set up by mtc_define
%   only the spatial arrangement matters; fields that depend on ng are not used
% opts (optional) options
%    opts.iflog=1 to log verifications (defaults to 0)
%
% flfs.cats: a structure for each category, flfs.cats{ic}.type: category type
% flfs.checkdef: (from mtcs) check definitions, i.e., coordinates of each check letter
% flfs.nchecks: (from mtcs) number of checks
% flfs.region_mask:(from mtcs) 1's indicating the glider region
% flfs.region_mask_letters:(from mtcs) letters indicating the glider region
% flfs.region_extent (from mtcs) height and width of region
% flfs.subcats: a structure for each subcategory (i.e., of a position of each cagegory within the glider)
%      subcats{is}.cat:      parent category number
%      subcats{is}.char:     a letter (a-z) to be used in enumerations
%      subcats{is}.type:     category type
%      subcats{is}.position: offset of category-glider inside of main glider (one entry is NaN if a row or col category)
%      subcats{is}.letters:  letter designation of positoin within the parent glider
%      subcats{is}.coords:   coordinates of the letters
%      subcats{is}.disjoint: subcategories that are disjoint
%      subcats{is}.precedes: subcategories that are disjoint and have a higher subcategory index
%
% flfs=flfs_define([],mtc_define(4,gtc_define([0 0;0 1;0 2;1 2])),setfield([],'iflog',1))
%
% See also:  MTC_DEFINE, FLFS_SEQREDUCE, FLFS_ENUMERATE, FLFS_ENUMERATE_DEMO.
%
from_mtcs={'checkdef','nchecks','region_mask','region_mask_letters','region_extent'};
if (nargin<=2)
    opts=[];
end
opts=filldefault(opts,'iflog',0);
if (isempty(cats))
    cats{1}.type='row';
    cats{2}.type='col';
end
flfs=[];
flfs.cats=cats;
for ifm=1:length(from_mtcs)
    flfs.(from_mtcs{ifm})=mtcs.(from_mtcs{ifm});
end
letters=fieldnames(flfs.checkdef);
%
% create subcategory list for each category
%
%pointers from cats to subcats in cats_to_subcats{ic}
%ponters from subcats to cats in subcats{is}.cats
flfs.subcats=[];
for ic=1:length(cats)
    flfs.cats{ic}.subcats=[];
    switch cats{ic}.type
        case {'row','col'}
            if strcmp(cats{ic}.type,'col')
                id=1;
            end
            if strcmp(cats{ic}.type,'row')
                id=2;
            end
            ids=3-id;
            nsubcats=flfs.region_extent(ids);
            for isubcat=1:nsubcats
                pos=isubcat-1;
                subcat.cat=ic; %set up pointer from subcategory to category
                subcat.type=cats{ic}.type;
                subcat.position=[NaN NaN]; %zero-index position of index cell of the category-glider within the overall glider
                subcat.position(ids)=pos;
                subcat.letters=[];
                subcat.coords=[];
                for icheck=1:flfs.nchecks
                    checkcoords=flfs.checkdef.(letters{icheck});
                    if (checkcoords(ids)==pos)
                        subcat.letters=cat(2,subcat.letters,letters{icheck});
                        subcat.coords=[subcat.coords;checkcoords];
                    end
                end
                flfs.subcats{end+1}=subcat; %
                flfs.cats{ic}.subcats(1,end+1)=length(flfs.subcats); %set up pointer from categories to subcategories
            end
        otherwise
            warning(sprintf('flfs_define: unknown category type %s, ignored.'));
    end %switch
end
for is=1:length(flfs.subcats)
    flfs.subcats{is}.char=char(double('a')-1+is);
    disjoint=[];
    for idj=1:length(flfs.subcats)
        if isempty(intersect(flfs.subcats{is}.letters,flfs.subcats{idj}.letters))
            disjoint=[disjoint,idj];
        end
    end
    flfs.subcats{is}.disjoint=disjoint;
    flfs.subcats{is}.precedes=disjoint(find(disjoint>is));
    if isempty(flfs.subcats{is}.precedes) %to avoid explicit 1x0 arrays
        flfs.subcats{is}.precedes=[];
    end
end
if (opts.iflog)
    disp(flfs.region_mask_letters);
    for ic=1:length(flfs.cats)
        disp(sprintf('category %2.0f:  type %6s',ic,flfs.cats{ic}.type));
        for is=flfs.cats{ic}.subcats
            disp(sprintf(' subcategory tag %2.0f (%2s): (parent category %2.0f, type %6s) position [%5.0f %5.0f] letters %s',...
                is,flfs.subcats{is}.char,flfs.subcats{is}.cat,flfs.subcats{is}.type,flfs.subcats{is}.position,flfs.subcats{is}.letters));
        end
    end
end
return
