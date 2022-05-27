% Example script to plot phasor histograms

%% load data
load coumarin2.mat; 

%% average decay over the frame
datamean = squeeze(mean(mean(data,2),3));

%% shift so peak laser power is at zero time point
% start each at peak (actual laser pulse may be back a few ticks, but this
% is easy to make consistent)
% Use the average of all pixels to get robust stats on the peak
tshift=find(datamean==max(datamean))
datameanshifted=circshift(datamean,-tshift);
datashifted=circshift(data,-tshift,1);

%normalize to max
datameanshiftednorm=datameanshifted/max(datameanshifted);

%% phasor histogram
% need the fractional times for each data point to properly scale the
% phasor sine and cosine to the laser rep rate
fractimes = 12.5*(1:256)/256;
dropout = 196;
fractimes(dropout:end) = fractimes(dropout:end)+0*2.5;
fractimes = fractimes/12.5;

% compute phasor vals 
[s g] = phasorvals(datashifted,fractimes);

%create hist vals
phist = phasorhist(s,g,400);

% labels and universal circle
% attempt to make time scale meaningful
% time period of data is p ns (12.5ns for 80MHz reprate)
p=12.5;
% the points along the universal circle correspond to pure lifetimes with
% lifetime t0 such that g = p^2/(p^2+t0.^2) and s = (p*t0)./(p^2+t0.^2)
t0s = logspace(log10(0.25),log10(45),10); % a set of lifetimes to label u-plot
g_t0s = p^2./(p^2+(2*pi*t0s).^2); 
s_t0s = (2*pi*p*t0s)./((2*pi*t0s).^2+p^2);


figure(1)
subplot(321), semilogy(fractimes*p,datameanshifted),xlabel('t [ns]'),ylabel('normalized mean signal')
subplot(322), imagesc(squeeze(mean(data,1))), title('intensity image')%, caxis([.28 .39])
subplot(323), imagesc(s), title('s')%, caxis([.28 .39])
subplot(324), imagesc(g), title('g')%, caxis([.42 .58]),caxis([.4 .6])

subplot(313)
imagesc([0 1],[0 0.5],phist), axis xy%, caxis([0 100])
xlabel('g'),ylabel('s'), title('phasor histogram')
gcoords = 0:.01:1; hold on, uc = plot(gcoords,sqrt(.25-(gcoords-.5).^2),g_t0s,s_t0s,'+'); text(g_t0s,s_t0s,cellstr(num2str(t0s',2))); hold off
