function [result_struct,counts,msgs_warn,msgs_all,opts_used]=mtc_cgsstruct_merge(orig_struct,addin_struct,opts)
% [result_struct,counts,msgs_warn,msgs_all,opts_used]=mtc_cgs_struct_merge(orig_struct,addin_struct,opts) merges one
% coordinate group structure into another.
%
% orig_struct:  a coordinate group structure, some fields may be absent
% addin_struct:  a coordinate group structure to be added in to orig_struct
% opts:  options for merging and reporting
%    opts.mismatch_action: options are
%        'keep_orig','keep_addin','keep_average','error'; default, set by mtc_defopts, is keep_orig
%    opts.nowarn_merge: suppress warnings if mismatch, default, set by mtc_defopts, is 0
%
% result_struct: the merged structure
% counts:  a summary of actions
%   counts.n_orig:   number of fields in orginal structure
%   counts.n_addin:  number of fields in structure to be added in
%   counts.n_result: number of fields in resulting structure
%   counts.n_new:    number of fields that are added in (i.e., absent in orig but present in addin)
%   counts.n_matches: number of fields present in original and addin, that match
%   counts.n_mismatches:  number of fields present in original and addin, that do not match
%
%   counts.n_result=counts.n_orig+counts.n_new;
%   counts.n_addin=counts.n_new+counts.n_matches+counts.n_mismatches
%
% msgs_all:  all actions
% msgs_warn: only warning actions
%
% opts_used:  options with missing entries filled in by mtc_defopts.
%
%  Note that validity of fields are not checked.  If there is a length mismatch, an error will occur.
%  Note also that an empty field and an explicit unbiased field are treated differently --
%    an empty field matches with anything, but an explict unbiased field must match within opts.tol_match
%
%   See also:  MTC_DEFOPTS, MTC_AUGCOORDS, MTC_CGSSTRUCT_CLEAN, MTC_CGSSTRUCT_CHECK.
%
if (nargin<=2)
    opts=[];
end
opts=mtc_defopts(opts);
opts_used=opts;
%
if isempty(orig_struct)
    orig_struct=struct();
end
if isempty(addin_struct)
    addin_struct=struct();
end
fnames_orig=fieldnames(orig_struct);
fnames_addin=fieldnames(addin_struct);
counts.n_orig=length(fnames_orig);
counts.n_addin=length(fnames_addin);
counts.n_new=0;
counts.n_matches=0;
counts.n_mismatches=0;
%
result_struct=orig_struct;
counts.n_result=counts.n_orig;
%
msgs_all=[];
msgs_warn=[];
%
for iadd=1:counts.n_addin
    fname=fnames_addin{iadd};
    if isfield(orig_struct,fname)
        v_orig=orig_struct.(fname);
        v_addin=addin_struct.(fname);
        match=0;
        if length(v_orig)==length(v_addin)
            dev=max(abs(v_orig-v_addin));
            if dev<opts.tol_match
                match=1;
            else
                mmsg=sprintf('dev: %f, s.b. < %f',dev,opts.tol_match);
            end
        else
            mmsg='unequal lengths';
        end
        if (match)
            msgs_all=strvcat(msgs_all,sprintf('field %15s: match in original and addin structures',fname));
            counts.n_matches=counts.n_matches+1;
        else
            counts.n_mismatches=counts.n_mismatches+1;
            mismatch_msg=sprintf('field %15s: mismatch (%s) in original and addin structures',fname,mmsg);
            switch opts.mismatch_action
                case 'keep_orig'
                    mismatch_msg=cat(2,mismatch_msg,'; original kept.');
                    msgs_all=strvcat(msgs_all,mismatch_msg);
                    msgs_warn=strvcat(msgs_warn,mismatch_msg);
                    if (opts.nowarn_merge==0)
                        disp(mismatch_msg);
                    end                   
                case 'keep_addin'
                    result_struct.(fname)=v_addin;
                    mismatch_msg=cat(2,mismatch_msg,'; addin kept.');
                    msgs_all=strvcat(msgs_all,mismatch_msg);
                    msgs_warn=strvcat(msgs_warn,mismatch_msg);
                    if (opts.nowarn_merge==0)
                        disp(mismatch_msg);
                    end                   
                case 'keep_average'
                    result_struct.(fname)=(v_orig+v_addin)/2;
                    mismatch_msg=cat(2,mismatch_msg,'; average kept.');
                    msgs_all=strvcat(msgs_all,mismatch_msg);
                    msgs_warn=strvcat(msgs_warn,mismatch_msg);
                    if (opts.nowarn_merge==0)
                        disp(mismatch_msg);
                    end                   
                case 'error'
                    msgs_all=strvcat(msgs_all,mismatch_msg);
                    msgs_warn=strvcat(msgs_warn,mismatch_msg);
                    if (opts.nowarn_merge==0)
                        disp(mismatch_msg);
                    end                   
                    error(mismatch_msg);
            end %mismatch_option
        end % match
            
    else %field to add in, no match
        result_struct.(fname)=addin_struct.(fname);
        counts.n_result=counts.n_result+1;
        counts.n_new=counts.n_new+1;
        msgs_all=strvcat(msgs_all,sprintf('field %15s: added in',fname));
    end
end
return
