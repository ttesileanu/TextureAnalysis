function [gain, predictions, details] = getPredictionsFromTernaryStats(...
    niEv, ppData, varargin)
% getPredictionsFromTernaryStats Calculate predictions for the ternary
% texture case starting from natural image patches.
%   [gain, predictions] = getPredictionsFromTernaryStats(niEv, ppData)
%   returns the `gain` matrix obtained by applying the "variance = salience"
%   idea based on efficient coding. This uses the texture statistics
%   contained in the `niEv` matrix and assumes isotropic input and output
%   noise. The function also uses this gain matrix to generate threshold
%   predictions in the directions indicated by the `groups` and
%   `directions` fields of the structure `ppData`, and finally scales
%   these to optimally match the `thresholds` from `ppData` (see `fitScale`
%   for details).
%
%   Note that the first output argument, `gain`, is the gain matrix scaled
%   to include the best fit coefficient returned from `fitScale`. The
%   unscaled matrix is returned in the `details` structure (see below).
%
%   The "variance = salience" method is accurate when output noise is
%   negligible *and* the signal is small compared to the input noise. If
%   this isn't the case, the options passed to 'efficientCodingOptions' can
%   be used to perform a more accurate calculation.
%
%   [..., details] = getPredictionsFromTernaryStats(...) returns some
%   details of the calculation, including the `a` and `mse` outputs of
%   fitScale.
%
%   Options:
%    'fitScaleOptions'
%       Options to be passed to `fitScale`.
%    'efficientCodingOpts'
%       Options to be passed to `solveLinearEfficientCoding`.
%
%   See also: solveLinearEfficientCoding, gainsToThresholds, fitScale.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('fitScaleOptions', {}, @(c) iscell(c) && isvector(c));
parser.addParameter('efficientCodingOptions', {}, @(c) iscell(c) && isvector(c));

% show defaults
if nargin == 1 && strcmp(niEv, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% extend NI to full probability values, unless it's already been done
if size(niEv, 2) == 66
    niEv = expandTextureStats(niEv, 3);
end

% extend PP directions to the full dimensionality of texture space
directionsExt = cellfun(@(v) v - 1/3, ternaryextdir(ppData.groups, ppData.directions), ...
    'uniform', false);

% calculate the gain matrix
niCov = cov(niEv);
unscaledGain = solveLinearEfficientCoding(niCov, params.efficientCodingOptions{:});

% calculate thresholds from the gains
unscaledPredictions = gainsToThresholds(unscaledGain, directionsExt);

% now use fitscale to match predictions to measurements
[aCoeff, predictions, predMse] = fitScale(unscaledPredictions, ppData.thresholds, ...
    'stds', ppData.thresholdIntervals, params.fitScaleOptions{:});

% prepare the details structure
details.niCov = niCov;
details.unscaledPredictions = unscaledPredictions;
details.aCoeff = aCoeff;
details.predMse = predMse;
details.unscaledGain = unscaledGain;

% rescale the gain to give the optimally-scaled predictions when used with
% gainsToThresholds
gain = unscaledGain / aCoeff;

end