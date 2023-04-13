% Creates a list of all pictures and their exposure values in a certain directory.
%
% Note that the directory must only contain images which are named according to the
% naming conventions, otherwise this function fails.
% 
% Filename naming conventions:
% The filename of a picture must contain two numbers specifying the
% exposure time of that picture. The first number specifies the numerator,
% the second one the denominator. E.g. "image_1_15.jpg" specifies that this
% image has been taken with an exposure of 1/15 second.
function [filenames, exposures, numExposures] = readDir(dirName)

thumbs_file = fullfile(dirName,'thumbs.db');
if exist (thumbs_file, 'file')
    delete(thumbs_file)
end
filelist = dir(dirName);
for i = 3:size(filelist,1)
    filenames{i-2} = strcat(dirName,filelist(i).name);
end
i = 1;
for filename = filenames
    [s f] = regexp(filename, '(\d+)');
    string = filename{1};
    %        numerator = string(s(1):f(1)); %GZH (fen zi)
    %        denominator = string(s(2):f(2));    (fen mu)
    %         numerator = string(s{1}(1):f{1}(1));
    %         denominator = string(s{1}(2):f{1}(2));
    numerator = string(s{1}(end-1):f{1}(end-1));
    denominator = string(s{1}(end):f{1}(end));
    exposure = str2num(numerator) / str2num(denominator);
    exposures(i) = exposure;
    i = i + 1;
end
% sort ascending by exposure
[exposures indices] = sort(exposures);
filenames = filenames(indices);
% then inverse to get descending sort order
exposures = exposures(end:-1:1);
filenames = filenames(end:-1:1);
numExposures = size(filenames,2);




