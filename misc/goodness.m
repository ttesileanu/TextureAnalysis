function stats = goodness(measurements, predictions)
% GOODNESS Calculate several measures of goodness of fit.
%   stats = GOODNESS(measurements, predictions) returns a set of
%   goodness-of-fit statistics to check how well the `predictions` match
%   the `measurements`. The following quantities are returned as fields in
%   the `stats` structure:
%     'relErrorMean'
%     'relErrorMedian'
%       Mean and median relative error, where the relative error for each
%       measurement is calculated as 2 * abs(y - x) / (y + x), with `x` the
%       measurement value, `y` the prediction.
%     'logErrorMean'
%     'logErrorMedian'
%       Mean and median error in log space, defined as abs(log(x) - log(y))
%       for measurement `x` and prediction `y`.
%     'logRms'
%       Root-mean-squared error in log space,
%           sqrt(mean((log(measurements) - log(predictions)).^2)) .
%       This is the D measure defined in the paper's SI.
%     'absRms'
%       Root-mean-squared error,
%           sqrt(mean((measurements - predictions).^2)) .
%     'absRmsOverVar'
%       Mean-squared error divided by the variance in the data,
%           absRms^2 / var(measurements)
%       This can be thought of as "unexplained variance", although note
%       that this quantity is not guaranteed to be smaller than 1, so
%       interpreting it as a fraction is problematic.

absErrors = abs(measurements(:) - predictions(:));

means = (measurements(:) + predictions(:)) / 2;
relErrors = absErrors ./ means;

x = log(measurements(:));
y = log(predictions(:));
logErrors = abs(x - y);

stats = struct;
stats.relErrorMean = mean(relErrors);
stats.relErrorMedian = median(relErrors);

stats.logErrorMean = mean(logErrors);
stats.logErrorMedian = median(logErrors);

stats.logRms = sqrt(mean(logErrors .^ 2));

stats.absRms = sqrt(mean(absErrors .^ 2));

stats.absRmsOverVar = stats.absRms^2 / var(measurements, 1);

end
