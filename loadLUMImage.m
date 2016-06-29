function LUM_Image = loadLUMImage(fname)
% loadLUMImage Load an image in luminosity format.
%   LUM_Image = loadLUMImage(fname) loads an image in luminosity format. If
%   the file is a Matlab file, then it must contain a 'LUM_Image' variable,
%   and this is directly returned. Otherwise the file must be an image
%   file, in which case it will be processed using processImage.m.
%
% See also: processImage.

[~, ~, ext] = fileparts(fname);

LUM_Image = [];
if strcmpi(ext, '.mat')
    contents = open(fname);    
    
    try
        LUM_Image = contents.LUM_Image;
    catch
    end
end
    
if isempty(LUM_Image)
    LUM_Image = processImage(fname);
end

end