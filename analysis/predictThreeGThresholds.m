function [thresholds, threshold_locs] = predictThreeGThresholds(noise_radius2, ...
    stats_mapping, directions, trafo, map_locs)
% predictThreeGThresholds Predict thresholds for G=3 axes given noise
% radius in transformed continuous-texture coordinates.
%   thresholds = predictThreeGThresholds(noise_radius2, stats_mapping, ...
%           directions, trafo, map_locs)
%   predicts the thresholds in the G=3 texture trajectories identified by
%   the `stats_mapping` cell array, given a squared noise radius
%   `noise_radius` in continuous-texture coordinates, after they are
%   transformed using the transformation matrix `trafo`. Each entry in
%   `stats_mapping` is a cell array of continuous statistics obtained for
%   patches at a fixed location in G=3 space. Each entry is thus a matrix
%   having 10 columns and an arbitrary number of columns, depending on how
%   many patches are generated at each G=3 location. The locations along the
%   G=3 axis are those identified by the vector `map_locs`. The cell array
%   `directions` matches the `stats_mapping` array, and contains the strings
%   'pos' or 'neg'. These show whether a threshold should be looked for in
%   the positive direction (towards higher values in `map_locs`), or in the
%   negative direction (towards lower values in `map_locs`) along the axis.
%
%   [thresholds, threshold_locs] = predictThreeGThresholds(...) also
%   returns the continuous-texture space locations of the thresholds that
%   have been found.

n = length(directions);

threshold_locs = cell(1, n);
thresholds = zeros(1, n);
for i = 1:n
    found_loc = [];
    found_thresh = [];
    if strcmp(directions{i}, 'pos')
        crt_dot_locs = map_locs(map_locs >= 0);
        crt_stats = stats_mapping{i}(map_locs >= 0);
    else
        crt_dot_locs = fliplr(map_locs(map_locs <= 0));
        crt_stats = fliplr(stats_mapping{i}(map_locs <= 0));
    end
    for k = 1:length(crt_dot_locs)-1
        loc1 = mean(crt_stats{k}, 1);
        loc2 = mean(crt_stats{k+1}, 1);
        dot_diff = crt_dot_locs(k+1) - crt_dot_locs(k);
        dloc = loc2 - loc1;
        
        p_a = dloc(:)'*trafo*dloc(:);
        p_b = dloc(:)'*trafo*loc1(:);
        p_c = loc1(:)'*trafo*loc1(:);
        
        crt_disc = p_b^2 + p_a*(noise_radius2 - p_c);
        if crt_disc >= 0
            crt_t = (-p_b + sqrt(crt_disc)) / p_a;
        else
            crt_t = [];
        end
        if ~isempty(crt_t) && crt_t >= 0 && (crt_t < 1 || k == length(crt_dot_locs)-1)
            found_loc = loc1 + crt_t*dloc;
            found_thresh = crt_dot_locs(k) + crt_t*dot_diff;
        end
    end
    if ~isempty(found_thresh)
        thresholds(i) = abs(found_thresh);
        threshold_locs{i} = found_loc;
    else
        thresholds(i) = nan;
    end
end

end