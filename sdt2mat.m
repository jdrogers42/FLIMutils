% Matlab script for converting Becker and Hickl SDT files from FLIM into 
% matlab files using the bioformats package for matlab:
% https://www.openmicroscopy.org/bio-formats/downloads/
% Bioformats will read the file into a sctructure. This script will find
% all files in the current working directory, read them in, convert the data
% from structures to arrays, and save each to a *.mat file of the same 
% filename. Each mat file will then contain two variables, 'data' which is
% a data cube with dimensions timesteps X ypixels X xpixels X channels, and 'info'
% which contains all of the meta data from that sdt file. 
%
% Note the docs for bioformats MATLAB: https://docs.openmicroscopy.org/bio-formats/5.3.4/developers/matlab-dev.html?highlight=matlab
% says:
%   This function returns an n-by-4 cell array, where n is the number of series in the dataset. If s is the series index between 1 and n:
%   The data{s, 1} element is an m-by-2 cell array, where m is the number of planes in the s-th series. If t is the plane index between 1 and m:
%   The data{s, 1}{t, 1} element contains the pixel data for the t-th plane in the s-th series.
%   The data{s, 1}{t, 2} element contains the label for the t-th plane in the s-th series.
%   The data{s, 2} element contains original metadata key/value pairs that apply to the s-th series.
%   The data{s, 3} element contains color lookup tables for each plane in the s-th series.
%   The data{s, 4} element contains a standardized OME metadata structure, which is the same regardless of the input file format, and contains common metadata values such as physical pixel sizes
% 
% 
% Written by Jeremy D. Rogers <jdrogers5@wisc.edu> UW-Madison
% on July 21, 2015 [20150721]
% Last updated: 20211120 (updated to bioformats 6.7 and fixed issue with reading multiple channels)
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
for ifile=1:size(files,1)
    tic;
    %tmp = bfopen(files(ifile).name); % bioformats tools for matlab needs to be installed
    % tmp is now a cell array with tmp{s,1}{t,1} such that s is the series
    % number and t-th plane (time axis). Channel are stacked on the t-axis.
    % Pre-allocate the array: time X imY X imX X channel X series
    
    % changed to use bfReader which is easier to access metadata and also
    % potentially improves reading efficiency for large files in the future:
    reader = bfGetReader(files(ifile).name)
    nSeries = reader.getSeriesCount; % number of series
    nChannels = reader.getSizeC; % Channels
    nTime = reader.getSizeT; % time bins
    nZ = reader.getSizeZ; 
    nY = reader.getSizeY;
    nX = reader.getSizeX;
    
    data = zeros(nTime,nY,nX,nZ,nChannels,nSeries);
    dimorder = reader.getDimensionOrder; % TODO: use this to unwrap in case it changes in the future, but for now, will hard code it
        
    for iSeries = 1:nSeries % loop over Series
      reader.setSeries(iSeries-1)
      for iChannels = 1:nChannels % loop over channels
        for iTime = 1:nTime % loop over time points (multiple channels are wrapped here)
          for iZ = 1:nZ
            % the getIndex method linearizes the wrapped t,c, and z planes
            data(iTime,:,:,iZ,iChannels,iSeries) = bfGetPlane(reader,reader.getIndex(iZ - 1, iChannels -1, iTime - 1) + 1);
            %tmp{iseries,1}{itime,1}; % converts structure to 3D matrix
          end
        end
      end
    end
    
    % while some meta data for preallocating the array is best accessed
    % with the reader, the full metadata is not easily converted from java
    % array to matlab type, and therefore cannot be saved in the mat file.
    % This is horribly inefficient and should be fixed, but for now, just
    % read in the neta data using bfopen and do the rest with the reader.
    tmp = bfopen(files(ifile).name);
    metadata = tmp{1, 2};
    
%     tmpname = tmp{2}.keys;
%     tmpval = tmp{2}.elements;
%     for itime=1:tmp{2}.size
%         info(itime).name = tmpname.nextElement;
%         info(itime).val = tmpval.nextElement;
%     end
    metadataKeys = metadata.keySet().iterator();
    for i=1:metadata.size()
        key{i} = metadataKeys.nextElement();
        value{i} = metadata.get(key{i});
        %fprintf('%s = %s\n', key{i}, value{i})
    end

    % save the meta data to info:
    info{1} = "data is a FLIM data cube with dimensions (T, Y, X, Z, Channel, Series); metadata keys are in info{2} and values are in info{3}";
    info{2} = key
    info{3} = value
    %metadata = reader.getCoreMetadataList;
    
    % save the mat file
    save(files(ifile).name(1:(end-4)), 'data','info', 'metadata')
    reader.close
    t=toc;
    display(['saved ' files(ifile).name(1:(end-4)) '.mat in ' num2str(t) ' sec'])
end
%clear tmp tmpname tmpval files i ii data info % no need to clear if used as a function, butdoesn't hurt
end
