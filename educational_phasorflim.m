
% questions for Matt
% 1. what about centroid of g as opposed to g-max?
% 2. can we segment and look at different centroids
% 3. what about getting tau2 or tau3 from data
% 4. mode locking and rep rate stability
% 5. 

% Matlab script for phasor analysis of FLIM data. This is an exploratory
% script and includes some code blocks for exploring several aspects of 
% phasors including the affect (artifacts) of shifts or frequency scaling 
% on exponentials. It also includes reading in SDT files using the
% bioformats package for matlab http://downloads.openmicroscopy.org/bio-formats/5.1.3/ 
% 
% Written by Jeremy D. Rogers <jdrogers5@wisc.edu> UW-Madison
% on July 21, 2015 [20150721]
% Last updated: 20150904
%
% Notes:
% Each section can be run seperately with ctrl+Enter
% If you run the entire script (i.e. using F5), some sections are disabled 
% by "if 0 ... end" blocks. You can run them by changing 0 to 1 or run
% the contents block by block with ctrl+enter.


%% play with some basic test functions 
N = 1000; % size of arrays 
tmax = N; % max time (arbitrary units)
t = (1:tmax)-1;
unicircg = (1:N)/N; unicircs = sqrt(.5^2-(unicircg-.5).^2);

t1 = .40*N; % lifetimes of free and bound
t2 = 0.05*N;

% determine concentration scaling coefficient
f1 = 0.5*exp(-t/t1);  % first exponential decay
f2 = 0.5*exp(-t/t2);  % second
bcoeff = 1-sum(f2)/sum(f1);

% Note on normalization: The exponential functions should be sample
% probability functions with integral equal to coefficient a. 
f1 = exp(-t/t1)/(t1/N); % first exponential decay
f2 = exp(-t/t2)/(t2/N); % second
f = f1/2+f2/2; %f = f+0.2*exp(-t/(0.05*N));
fcos = cos(2*pi*t/tmax);
fsin = sin(2*pi*t/tmax);

figure(1) % show a basic exponential
subplot(211)
plot(t,f,t,f1,t,f2,t,fcos,'k--',t,fsin,'k--'),ylim([-1 4])
subplot(212)
g = sum(f.*fcos)/sum(f); s = sum(f.*fsin)/sum(f);
g1 = sum(f1.*fcos)/sum(f1); s1 = sum(f1.*fsin)/sum(f1);
g2 = sum(f2.*fcos)/sum(f2); s2 = sum(f2.*fsin)/sum(f2);
plot(g,s,'o',g1,s1,'o',g2,s2,'o',unicircg,unicircs,'k')
xlabel('g [cos/real]'), ylabel('s [sin/imag]')


%% 
in = 0:.05:1-.05;
for iii=1:length(in)

a1 = .5;   % weight fraction and time constant
a1 = in(iii);
%a1 = (bcoeff-1)*(iii/length(in)-1)./(1+bcoeff*(iii/length(in)-1));
as(iii)=a1;
a2 = 1-a1;

% Note on normalization: The exponential functions should be sample
% probability functions with integral equal to coefficient a. 
f1 = a1*exp(-t/t1)/(t1/N); % first exponential decay
f2 = a2*exp(-t/t2)/(t2/N); % second
f = f1+f2; %f = f+0.2*exp(-t/(0.05*N));
fcos = cos(2*pi*t/tmax);
fsin = sin(2*pi*t/tmax);

figure(1),clf
subplot(311)
plot(t,f,t,f1,t,f2,t,fcos,t,fsin),ylim([-1 5])

% How the offset or shift affect the point on the phasor plot?
% As the function is shifted to the right, the phasor values circumscribe a
% circle centered at the origin. Try this for different t1 values and you 
% will see that the circle changes radius (large radius corresponds to fast
% decay). 
% The 'correct' shift value that corrsponds to the exponential  
% starting at time zero (blue circle) corresponds to the point on this blue
% circle that intersects the universal circle. 
% Intersections points for fast decay are on the right, long on the left.
for i=1:100;
g1(i) = sum(circshift(f1,[0,round(N*(i-1)/100)]).*fcos)/sum(f1);
s1(i) = sum(circshift(f1,[0,round(N*(i-1)/100)]).*fsin)/sum(f1);
g2(i) = sum(circshift(f2,[0,round(N*(i-1)/100)]).*fcos)/sum(f2);
s2(i) = sum(circshift(f2,[0,round(N*(i-1)/100)]).*fsin)/sum(f2);
g(i) = sum(circshift(f,[0,round(N*(i-1)/100)]).*fcos)/sum(f);
s(i) = sum(circshift(f,[0,round(N*(i-1)/100)]).*fsin)/sum(f);
end
unicircg = (1:N)/N; unicircs = sqrt(.5^2-(unicircg-.5).^2);
subplot(312)
sn=1;
plot(g1,s1,g1(sn),s1(sn),'bo',g2,s2,g2(sn),s2(sn),'ro',g,s,g(sn),s(sn),'yo',unicircg,unicircs,'k')
xlim([-1 1]),ylim([-1 1])

% What is the affect of using the wrong frequency to evaluate the phasor?
% FFT provides lots of freq points.
% Check what the affect of using different (incorrect) period of sinusoid
% is by evaluating fft and looking at each freq. The 'correct' freq to look
% at is one period per array which corresponds to the second point in the
% fft (highlighted by circle markers in the plot). 
% The value of pure exponentials lie on the universal circle.
% Superpositions of multiple exponentials produces a point that lies on a
% line connecting the two. The position along this line is not a linear
% function of the two amplitudes. This is probably more like the integral
% and an analytical mapping can probably be made to determine the scale.
ft = fft(f)/sum(f); ft1 = fft(f1)/sum(f1);ft2 = fft(f2)/sum(f2);
subplot(313)
plot(real(ft1),-imag(ft1),real(ft2),-imag(ft2),real(ft),-imag(ft),...
     real(ft1(2)),-imag(ft1(2)),'bo',real(ft2(2)),-imag(ft2(2)),'ro',real(ft(2)),-imag(ft(2)),'yo')
xlim([0 1]),ylim([-.5 .5])
% Interesting to note: It appears (need to confirm this mathematically)
% that a sum of 3 exponentials would result in a triangular space with
% vertices on the universal circle in which the phasor point of the total
% function can lie. In some of the data, the historgram appears somewhat
% traingular. I suspect that this says something about the number of pure
% states. Perhaps the pure states can be extrapolated with confidence
% intervals.
bs(iii)=sqrt((real(ft(2))-real(ft1(2)))^2+(imag(ft(2))-imag(ft1(2)))^2)...
        /sqrt((real(ft2(2))-real(ft1(2)))^2+(imag(ft2(2))-imag(ft1(2)))^2);
drawnow
pause(.1)
end



%bcorr = (2*bs(round(length(as)/2))-1)/(bs(round(length(as)/2))-1);
%plot(as,bs,'-+',as,1-as./(1+(as-1)*bcoeff))


%% try two different combinations that result in same point in phasor space
N = 1000; % size of arrays 
tmax = N; % max time (arbitrary units)
t = (1:tmax)-1;


t1 = 0.80*N; % lifetimes of free and bound
t2 = 0.1*N;
t3 = 0.3*N;
t4 = 0.05*N;


a1 = .35;   % weight fraction and time constant
a2 = 1-a1;
a3 = .75;   % weight fraction and time constant
a4 = 1-a1;


for i=1:100
% loop through shift
shift = (i-1)/100*2*pi;
% loop through freqs (comment out to use shift above
%shift=0; tmax = N*i/20;


% Note on normalization: The exponential functions should be sample
% probability functions with integral equal to coefficient a. 
f1 = a1*exp(-t/t1)/(t1/N); % first exponential decay
f2 = a2*exp(-t/t2)/(t2/N); % second
f3 = a3*exp(-t/t3)/(t3/N); % 
f4 = a4*exp(-t/t4)/(t4/N); % 
fa = f1+f2; 
fb = f3+f4;

fcos = cos(2*pi*t/tmax-shift);
fsin = sin(2*pi*t/tmax-shift);

figure(1),clf
subplot(311)
semilogy(t,fa,t,f1,t,f2,t,fb,':',t,f3,':',t,f4,':',t,fcos,t,fsin),ylim([-1 5])

g1(i) = sum(f1.*fcos)/sum(f1);
g2(i) = sum(f2.*fcos)/sum(f2);
ga(i) = sum(fa.*fcos)/sum(fa);
g3(i) = sum(f3.*fcos)/sum(f3);
g4(i) = sum(f4.*fcos)/sum(f4);
gb(i) = sum(fb.*fcos)/sum(fb);

s1(i) = sum(f1.*fsin)/sum(f1);
s2(i) = sum(f2.*fsin)/sum(f2);
sa(i) = sum(fa.*fsin)/sum(fa);
s3(i) = sum(f3.*fsin)/sum(f3);
s4(i) = sum(f4.*fsin)/sum(f4);
sb(i) = sum(fb.*fsin)/sum(fb);
drawnow
end

unicircg = (1:N)/N; unicircs = sqrt(.5^2-(unicircg-.5).^2);
subplot(312)
sn=1;
plot(ga,sa,ga(sn),sa(sn),'bo',gb,sb,gb(sn),sb(sn),'ro',unicircg,unicircs,'k',...
    g1(sn),s1(sn),'bo',g2(sn),s2(sn),'bo',g3(sn),s3(sn),'ro',g4(sn),s4(sn),'ro')
xlim([-1 1]),ylim([-1 1])

subplot(313)
%plot(real(ft1),-imag(ft1),real(ft2),-imag(ft2),real(ft),-imag(ft),...
%     real(ft1(2)),-imag(ft1(2)),'bo',real(ft2(2)),-imag(ft2(2)),'ro',real(ft(2)),-imag(ft(2)),'yo')
xlim([0 1]),ylim([-.5 .5])




%% process real data
if 1
% select which data to look at (comment out all but 1)
%load 2014-12-30_TG25_B6JOb(15wks)_10G-1.mat
%load 2014-12-30_TG25_B6JOb(15wks)_Rote-1.mat
load 2014-12-30_IRF_740_60X_zoom2_Coumarin-01.mat

% List all peices of header information saved in the files:
% for i=1:134;display([info(1).name ': ' num2str(info(1).val)]),end

% display each time point (psuedo-movie)
%figure(1),for i=1:size(data,1), imagesc(squeeze(data(i,:,:))),drawnow; end

figure(1),clf % average across time to show a 'brightfield' image
subplot(211),imagesc(squeeze(sum(data)))

% average across pixels to show the mean lifetime signal
subplot(212),semilogy(squeeze(mean(mean(data,2),3)))

%% try phasor plot with varying frequencies
for i=1:128
    shift = -59;
    scalefreq = .001+i/64;
    tmp = double(circshift(data(:,100:115,100:115),shift,1));
    fcos=cos(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1));
    fsin=sin(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1));
    meanflim = squeeze(mean(mean(tmp,2),3));
    figure(1),clf
    subplot(211)
    plot(1:size(tmp,1),meanflim/max(meanflim),1:size(tmp,1),fsin,1:size(tmp,1),fcos)
    % construct phasor arrays
    g = sum(tmp.*repmat(fcos,[1,size(tmp,2),size(tmp,3)]))./sum(tmp,1);
    s = sum(tmp.*repmat(fsin,[1,size(tmp,2),size(tmp,3)]))./sum(tmp,1);
    subplot(212)
    plot(reshape(g,1,[]),reshape(s,1,[]),'.'),xlim([0 1]),ylim([0 .5])
    gcoords = 0:.01:1; hold on, uc = plot(gcoords,sqrt(.25-(gcoords-.5).^2)); hold off
    drawnow
end

%% try phasor plot with varying shift
for i=1:100
    shift = -1*i;
    scalefreq = 1/1;
    tmp = double(circshift(data(:,80:115,80:115),shift,1));
    fcos=cos(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1));
    fsin=sin(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1));
    meanflim = squeeze(mean(mean(tmp,2),3));
    figure(1),clf
    subplot(211)
    plot(1:size(tmp,1),meanflim/max(meanflim),1:size(tmp,1),fsin,1:size(tmp,1),fcos)
    % construct phasor arrays
    g = sum(tmp.*repmat(fcos,[1,size(tmp,2),size(tmp,3)]))./sum(tmp,1);
    s = sum(tmp.*repmat(fsin,[1,size(tmp,2),size(tmp,3)]))./sum(tmp,1);
    subplot(212)
    plot(reshape(g,1,[]),reshape(s,1,[]),'.'),xlim([-1 1]),ylim([-1 1])%,xlim([0 1]),ylim([-1 1])
    gcoords = 0:.01:1; hold on, uc = plot(gcoords,sqrt(.25-(gcoords-.5).^2)); c = plot(gcoords*2-1,sqrt(.3-(gcoords*2-1).^2)); hold off
    drawnow
    rad(i,:,:)=sqrt(g.^2+s.^2);
end

%% try making 2D histogram phasor plot
%file = '2014-12-30_IRF_740_60X_zoom2_Coumarin-01.mat';fileID = fopen('Coumarinh1.int');globalsdata = fread(fileID,[256,256],'float32')';fclose(fileID);
file = '2014-12-30_TG25_B6JOb(15wks)_10G-1.mat';fileID = fopen('TG25-10G(nomedfilter)h1.int');globalsdata = fread(fileID,[256,256],'float32')';fclose(fileID);
%file = '2014-12-30_TG25_B6JOb(15wks)_Rote-1.mat';fileID = fopen('TG25-Rotenone(nomedfilter)h1.int');globalsdata = fread(fileID,[256,256],'float32')';fclose(fileID);
%file = '2014-12-30_IRF_740_60X_zoom2_Coumarin-01.mat';fileID = fopen('TG25Coumarin(nomedfilternoautocenter).int');globalsdata = fread(fileID,[256,256],'float32')';fclose(fileID);
globalsdata = globalsdata(1:end/2,1:end); % get rid of useless padding from globals

load(file) % using string allows for labeling plots with filename later
bins = 256; % number of bins in g and s space
timecrop = [26 231];%[26 231]; % start and end bins in time (should be 1 to 9 ns)
shift = (timecrop(1)-1)-58; % -58 puts coumerin right on the universal circle (manually), but looks like globals uses a fractional shift closer to 58.25 or so (interpolation or puts the shift in the sin and cos functions numerically)
shift=shift-14; % another 14 makes data fit to globals output for cells

% time scale params
laserrep = 1/80000000*1e9; datatimerange = 10; timeperpix = datatimerange/size(data,1);
scalefreq = 1;%timeperpix*(timecrop(2)-timecrop(1))/laserrep; %datatimerange/laserrep;

%tmp = double(circshift(data(:,100:115,100:115),shift,1));
%tmp = repmat([1:1000]',[1 128 128]); jitter = shiftdim(repmat(rand(128,128),[1 1 1000]),2); tmp = exp(-tmp./(100*jitter+100));
tmp = data(timecrop(1):timecrop(2),:,:);
tmp = double(circshift(tmp(:,1:end,1:end),shift,1)); 
tmp = tmp(1:end-0,:,:);
% fcos=cos(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1)-0.5*pi/256);
% fsin=sin(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1)-0.5*pi/256);
fcos=cos(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1));
fsin=sin(scalefreq*2*pi*(1:size(tmp,1))'/size(tmp,1));
meanflim = squeeze(mean(mean(tmp,2),3));

% plots
figure(2),clf
subplot(311)
plot(1:size(tmp,1),meanflim/max(meanflim),1:size(tmp,1),fsin,1:size(tmp,1),fcos)

% construct phasor arrays
g = sum(tmp.*repmat(fcos,[1,size(tmp,2),size(tmp,3)]))./sum(tmp,1);
s = sum(tmp.*repmat(fsin,[1,size(tmp,2),size(tmp,3)]))./sum(tmp,1);

% bin data into grid of points (need cell array for hist3
centers{1} = [1:(bins/2)]/(bins+1);
centers{2} = [1:bins]/(bins+1);
% use hist3 to calculate histogram, but plot seperately 
%g = 0*g+.995; s=0*s+.49; for confirming bin positions
a = hist3([reshape(s,1,[]); reshape(g,1,[])]', centers);
subplot(312)
    imagesc([0 1],[0 0.5],a), axis xy
    xlabel('g'),ylabel('s'),title(file(1:end-4),'Interpreter','none')
    gcoords = 0:.01:1; hold on, uc = plot(gcoords,sqrt(.25-(gcoords-.5).^2)); hold off
    %xlim([.2 .4]), ylim([.4 .5]) % use for coumarin
    %xlim([.3 .6]), ylim([.1 .25]) % use for cells
    %delete(uc)
subplot(313),
    imagesc([0 1],[0 0.5],(globalsdata)), axis xy
    hold on, ucg = plot(gcoords,sqrt(.25-(gcoords-.5).^2)); hold off
    xlabel('g'),ylabel('s'),title('Globals output')
    %xlim([.2 .4]), ylim([.4 .5]) % use for coumarin
    %xlim([.3 .6]), ylim([.1 .25]) % use for cells
    
%%
figure(3), hold on % Use surface instead
    [X Y] = meshgrid(centers{2},flipud(centers{1}));
%     mldata = surf(X, Y, a,'EdgeColor', 'None', 'FaceColor', 'Red','FaceAlpha', 0.5)
%     glbls = surf(X, Y, b,'EdgeColor', 'None', 'FaceColor', 'Blue','FaceAlpha', 0.5)
a1=.0001*a; b1=.0001*b; % scale down the z axis to make rotations look nicer
    mldata = surf(X(26:64,77:154), Y(26:64,77:154), a1(26:64,77:154),'EdgeColor', 'Red', 'FaceAlpha', 0.5)
    glbls = surf(X(26:64,77:154), Y(26:64,77:154), b1(26:64,77:154),'EdgeColor', 'Blue', 'FaceAlpha', 0.5)
    view(2)
    hold off

figure(4)
    image([0 1],[0 0.5],a, 'AlphaData', .5, 'CData',a/max(max(a))*64*.6),axis xy
    hold on
    colormap('Jet')
    image([0 1],[0 0.5],b, 'AlphaData', .5,'CData',b/max(max(b))*64*.4+.6*64)
    ucg = plot(gcoords,sqrt(.25-(gcoords-.5).^2)); hold off
    xlabel('g'),ylabel('s'),title('Mountain plot')
    %xlim([.2 .4]), ylim([.4 .5]) % use for coumarin
    %xlim([.3 .6]), ylim([.1 .25]) % use for cells
    
    
end


%% Notes about experiment parameters
% laser rep rate 80MHz
% data is 10 ns 
% crop is bin 26 (25 starting from 0) to 230 which is 1ns to 9ns
% coumarin is 2.5ns, think the componants bound and free .4 and 2.3-2.7
% for 10g tau1 = 435.51; tau2 2746.04, a1 70.84, a2 = 29.16, taumean =
% 107.46
% for rotenone, weird shift, differnent image gives tau1 414; tau2 2436, a1
% = 74; a2 =26

if 0
    %% Plot data as analyzed by the globals method
    %fileID = fopen('TG25-10G-1h1.int'); 
    %fileID = fopen('TG25-Rotenone-1h1.int'); 
    %fileID = fopen('Coumarinh1.int');
    fileID = fopen('TG25-10G(nomedfilter)h1.int');
    %fileID = fopen('TG25-Rotenone(nomedfilter)h1.int');
    A = fread(fileID,[256,256],'float32')';
    fclose(fileID)
    figure(5),imagesc([0 1],[0 0.5],A(1:128,:)), axis xy
end