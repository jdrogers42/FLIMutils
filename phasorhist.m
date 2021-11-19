% phist = phasorhist(s, g, N)
%
% Function takes the g and s arrays from FLIM data and bins the data to
% make a histrogram array that can be plotted on the 'universal circle'.
% The N argument is the histogram resolution in the x (g) axis, the y (s)
% axis will be half this value.
%
% TODO: 
% * do better checking on bins to make sure sampling is reasonable, and
% possibly just warn if it is not (this may be better as part of
% phasorvals() since that function knows about the counts)
% * consider allowing non-sqaure bins
% * consider allowing range to be greater than universal circle extent

function phist = phasorhist(s,g,bins)

% check that bins makes sense (even number)
if bins/2 ~= round(bins/2), 'bins should be an even number', end

% check that bins makes sense (sampling)


%% bin data into grid of points (need cell array for hist3
 centers{1} = [1:(bins/2)]/(bins+1);
 centers{2} = [1:bins]/(bins+1);
 
%% use hist3 to calculate histogram, but plot seperately 
 phist = hist3([reshape(s,1,[]); reshape(g,1,[])]', centers);

end

