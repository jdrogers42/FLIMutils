% Matlab script for converting Becker and Hickl SDT files from FLIM into 
% matlab files using the bioformats package for matlab:
% http://downloads.openmicroscopy.org/bio-formats/5.1.3/ 
% Bioformats will read the file into a sctructure. This script will find
% all files in the current working director, read them in, convert the data
% from structures to arrays, and save each to a *.mat file of the same 
% filename. Each mat file with then contain two variables, 'data' which is
% a data cube with dimensions timesteps X ypixels X xpixels, and 'info'
% which contains all of the meta data from that sdt file. 
% 
% Written by Jeremy D. Rogers <jdrogers5@wisc.edu> UW-Madison
% on July 21, 2015 [20150721]
% Last updated: 20150904
%
% Notes: use a matching pattern as an argument to select only certain files
% for example, sdt2mat('2015_*.sdt') will only convert filenames starting 
% with the '2015_'
%
% TODO: take path as argument so we don't need to change directory

function sdt2mat(matchpattern)

%% if no argument is given, assume you want to convert all sdt files in CWD
if ~exist('matchpattern')
    display('no matching pattern, converting all *.sdt files')
    matchpattern = '*.sdt';
end

%% read in files and save as mat files
% only need to do this once since it is slow. After the data is saved as a
% *.mat file, reading is much faster. 
files = dir(matchpattern); % get list of all files with sdt extension
for ii=1:size(files,1)
    tmp = bfopen(files(ii).name); % bioformats tools for matlab needs to be installed
    for i = 1:size(tmp{1,1},1)
        data(i,:,:) = tmp{1,1}{i,1}; % converts structure to 3D matrix
    end
    tmpname = tmp{2}.keys;
    tmpval = tmp{2}.elements;
    for i=1:tmp{2}.size
        info(i).name = tmpname.nextElement;
        info(i).val = tmpval.nextElement;
    end 
    save(files(ii).name(1:(end-4)), 'data','info')
    display(['saved ' files(ii).name(1:(end-4)) '.mat'])
end
%clear tmp tmpname tmpval files i ii data info % no need to clear if used as a function, butdoesn't hurt
end
