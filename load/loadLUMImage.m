function LUMImage = loadLUMImage(fname)
% loadLUMImage Load an image in luminosity format.
%   LUMImage = loadLUMImage(fname) loads an image in luminosity format. If
%   the file is a Matlab .mat file, then it must contain a 'LUMImage'
%   variable, and this is directly returned. Otherwise the file must be an
%   image file, in which case it will be processed using imageToLUM.m.
%
% See also: imageToLUM.

[~, ~, ext] = fileparts(fname);

LUMImage = [];
if strcmpi(ext, '.mat')
    contents = open(fname);    
    
    try
        LUMImage = contents.LUM_Image;
    catch
    end
end
    
if isempty(LUMImage)
    LUMImage = imageToLUM(fname);
end

end