function  [reduced,uncovered_cells,covered_cells,region_covered]=flfs_seqreduce(recent,flfs)
%[reduced,uncovered_cells,covered_cells,region_covered]=flfs_seqreduce(recent,flfs) reduces a sequence of fallen leaves to one in which each leaf is essential
%
% recent:  a row array of recently-fallen subcategories (1 to length(flfs.subcats)), recent(1) is the most recent
% flfs: structure defining falling-leaves/sticks categories and subcategories, returned by flfs_define
%
% reduced: a row array of the subcategories that are visible
% uncovered_cells:  a string of the region's cells that are not covered (uppercase)
% covered_cells:  a string with uppercase letters where cells are uncovered, lowercase indicating coverage by leaves;
%   tokens are in the order of the cells of the region.  For example, AdC
%   means that original cells A and C are uncovered, and original cell B is
%   replaced by the leaf of subcategory d
% region_covered:  a cell array with the information of covered_cells
%
% See also:  MTC_DEFINE, FLFS_ENUMERATE.
%
nsubcats=length(flfs.subcats);
x=recent(find((recent<=nsubcats)&(recent>=1))); %valid subcats only 
letters=fieldnames(flfs.checkdef);
%
letters_all=char(letters)';
uncovered_cells=letters_all;
covered_cells=letters_all;
region_covered=flfs.region_mask_letters;
reduced=[];
while ~isempty(uncovered_cells) & ~isempty(x)
    t=x(1); %consider the most recent subcategory that is not used yet
    x=x(find(x~=t)); %never need to look at this subcategory again
    newly_covered=intersect(uncovered_cells,flfs.subcats{t}.letters);
    %add it if it masks any letters
    if ~isempty(newly_covered)
        covered_pointers=find(ismember(letters_all,newly_covered));
        covered_cells(covered_pointers)=flfs.subcats{t}.char;
        uncovered_cells=setdiff(uncovered_cells,flfs.subcats{t}.letters);
        for ichar=1:length(newly_covered)
            coords=flfs.checkdef.(newly_covered(ichar));
            region_covered{1+coords(1),1+coords(2)}=flfs.subcats{t}.char;
        end
        reduced=[reduced,t];
        %now try further reduction according to precedence
        if length(reduced)>1
            %this is a pointer to potential positions in reduced that can be preceded by the current last token
            disjoint_with=find(ismember(reduced(1:end-1),flfs.subcats{t}.disjoint)); %leaves that are disjoint
            can_precede=find(ismember(reduced(1:end-1),flfs.subcats{t}.precedes)); %leaves that are disjoint and lower lexicographically
%            disp('reduced')
%            disp(reduced)
%            disp('can_precede')
%            disp(can_precede)
%            disp('disjoint_with')
%            disp(disjoint_with)
            if ~isempty(can_precede)
                % look for the longest consecutive sequence in disjoint_with that ends with length(reduced)-1
                need_match=[0:length(disjoint_with)-1]+[length(reduced)-length(disjoint_with)];
%               disp('need_match')
%               disp(need_match)
                matches=double(need_match==disjoint_with);
                last_mismatch=max(find(matches==0)); %pointer to last mismatch counting backwards, i.e., most recent
%               disp('matches')
%               disp(matches)
%               disp('last_mismatch')
%               disp(last_mismatch)
                if matches(end)==1
                    if isempty(last_mismatch)
                        last_mismatch_pos=need_match(1);
                    else
                        last_mismatch_pos=need_match(last_mismatch+1);
                    end
                    %last_mismatch_pos
                    can_swap=[last_mismatch_pos:length(reduced)-1]; %pointer to longest consecutive string of disjoint leaves
                    %find most recent in can_swap that is also in can_precede
                    first_swap=min(intersect(can_swap,can_precede));
                    if ~isempty(first_swap)
                        %the last token (in position length(reduced) floats to the
                        %position after first_swap
                        downshift=[first_swap:length(reduced)-1];
    %                   disp('downshift')
    %                   disp(downshift);
                        reduced=reduced([[1:first_swap-1],length(reduced),downshift]);
                    end %can_swap
                end
            end %possible swaps
        end %length>1
    end %a letter was added
%   reduced
%   disp('at cycle end: reduced, x, uncovered_cells')
%   reduced
%   x
%   uncovered_cells
end
return
