% [s g] = phasorvals(datacube [, fractionaltimes])
%
% Function takes a FLIM data cube and returns two image arrays of s and g
% values. To ensure the correct period is used, the optional arguments
% include the time values and repitition rate of the laser. 
%
% fractionaltimes (optional) are the times for each time bin divided by the laser
% repetition period. If this argument is not included, the default is to
% assume that the data spans the full repition period
%
% TODO: 
% * Better results may be possible by 'filling' in missing data within the
% laser repitition cycle with estimated counts based on a fit of the whole
% data cube
% * 

function [s g] = phasorvals(datacube, fractionaltimes)

datacube = double(datacube); % ensure data is double
switch nargin
    case 2
      if length(fractionaltimes) ~= size(datacube,1), 'error: length of fractionaltimes is not consistent with datacube',end
    case 1
      fractionaltimes = ((1:size(datacube,1))-1)/size(datacube,1);
    otherwise
      'I am confused about the inputs'
end

%% make sin and cos arrays
 fcos=cos(2*pi*fractionaltimes)';
 fsin=sin(2*pi*fractionaltimes)';

%% construct phasor arrays
 % Note: this normalizes the existing data to 1 for each pixel, which is probably the best
 % that can be done without more information, but will be an error when
 % there is missing data as happens when the reprate is longer than the
 % data acquisition (due to the time required to reset the integrator). A
 % better estimate may be to fit an exponential (or double) to the average 
 % lifetime curve from the whole image and then normalize the data to the 
 % fraction that exists within the data set accounting for the lost counts
 % during integrator reset. This may be tried in the future.
 g = squeeze(sum(datacube.*repmat(fcos,[1,size(datacube,2),size(datacube,3)]))./sum(datacube,1));
 s = squeeze(sum(datacube.*repmat(fsin,[1,size(datacube,2),size(datacube,3)]))./sum(datacube,1));
 
end 