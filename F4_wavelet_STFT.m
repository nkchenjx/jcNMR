% First coding: 2024/11/01, MATLAB version 2024a

% Copyright (c) 2024 Jixin Chen @ Ohio University
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clear;

dataPoint = 1024;

    

    NMRfreq = 500; % MHz
    Ts = 1E-4;% numcycle1ppm/dataPoint/NMRfreq; % unit s. Length of x to satisify ~cycle1ppm cycles for NMRfreq Hz.
    y = zeros(1, dataPoint);
    yf = y; %[y, zeros(1, length(y)*10)]; % if zero filling add a few times length of zeros after y.
    yfft = fft(yf);
    fs = 1/Ts;
    n = length(yf);
    fshift = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm




    %% set options
    NMRfreq = 500; % MHz
    % dataPoint = 1024; %number of data points
    xnoise = 1E-11;
    ynoise = 0.001;
    
    dfr = 0.001; %relative high resolution to get a smooth curve
    fr = -10.24:dfr:10.24; %frequency region shift (ppm)
    
    t = 0:dataPoint-1; % time points
    %numcycle1ppm = 10;
    Ts = 1E-4;% numcycle1ppm/dataPoint/NMRfreq; % unit s. Length of x to satisify ~cycle1ppm cycles for NMRfreq Hz.
    t = t*Ts;
    tn = t + randn(1,length(t))*xnoise; % add a little noise to time at ps level  
    
    cp1 = 1;
    cp2 = [1, 1]/2;
    cp3 = [1, 2, 1]/4;
    cp4 = [1, 3, 3, 1]/8;
    cp5 = [1, 4, 6, 4, 1]/16;
    cp6 = [1, 5, 10, 10, 5, 1]/32;
    cp7 = [1, 6, 15, 20, 15, 6, 1]/64;
    cp8 = [1, 7, 21, 35, 35, 21, 7, 1]/128;
    cp9 = [1, 8, 28, 56, 70, 56, 28, 8, 1]/256;
    
    
    meanHtype = 10; % average different H atom types 
    meanHt = 4; % average H atom in each type. Average H meanHtype*meanHt   from  https://doi.org/10.1016/j.jmr.2010.12.008
    maxJ = 20; % Hz, max J coupling value
    minJ = 0.5; % Hz, min J coupling value
    meanSigma = 4; % Hz, mean peak width
    minSigma = 1; % Hz, min peak width
    maxShift = 10; % ppm in max shift
    minShift = -5; % ppm min shift
    meanDamping = 0.5; % mean damping rate in s
    minDamping = 0.01; % min damping rate in second
    meanCP = 3; % average couping is 3.
    
    % equations:
    dampGauss = @(p, x, t) p(1)/sqrt(2*pi)/p(2)*dfr*exp(-((p(3)-x).^2/2/p(2)^2))*exp(-t/p(4));  % p1, amp; 2, sigma; 3, shift; 4, damping
    NMR = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))); % p1, amp; 2, shift; 3, pahse; %using +1i rotation only.
    
    %% initialize parameters
    numHtypes = ceil(abs(randn(1)*(meanHtype))); % number of types of H atoms
   
    H_num = ceil(abs(randn(1, numHtypes)*meanHt)); % number of H in each type
    H_shift = (maxShift-minShift)*abs(rand(1, numHtypes)) + minShift; % the chemical shift of each on in ppm 
    H_cp =  ceil(meanCP*abs(randn(1, numHtypes))); % coupling type of each H atom type
    H_J = maxJ*abs(rand(1, numHtypes)) + minJ; % the J value of the coupling
    H_sigma = meanSigma*abs(randn(1, numHtypes)) + minSigma;  %Hz
    H_phase0 = abs(randn(1, numHtypes)*pi/2); % set a initial phase as a Guassian centered at 0 and std pi/2
        H_phase0 = rem(H_phase0, 2*pi);
    H_damping = meanDamping*abs(randn(1, numHtypes)) + minDamping; % set the damping rate of each type
    %H_damping = meanDamping*(randn(1, numHtypes)) + minDamping; % set the damping rate of each type
    
    y = zeros(1, dataPoint);
    yf = y; %[y, zeros(1, length(y)*10)]; % if zero filling add a few times length of zeros after y.
    yfft = fft(yf);
    fs = 1/Ts;
    n = length(yf);
    fshift = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm

    %% calculate truth
    peakTrue = zeros(1, length(fr));
    peakTrueb = zeros(1, length(fshift));
    for i = 1:numHtypes
        amp = H_num(i);
        shift = H_shift(i); % ppm
        cp = H_cp(i); % type
        J = H_J(i); % Hz
        J = J/NMRfreq; % ppm
        if cp>9
            cp = meanCP;
        end
        switch cp
            case 1
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp1;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp1;
            case 2
                shift = shift + [-0.5, 0.5]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp2;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp2;
            case 3
                shift = shift + [-1, 0, 1]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp3;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp3;
            case 4
                shift = shift + [-1.5, -0.5, 0.5, 1.5]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp4;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp4;
            case 5
                shift = shift + [-2, -1, 0, 1, 2]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp5;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp5;
            case 6
                shift = shift + [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp6;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp6;
            case 7
                shift = shift + [-3, -2, -1, 0, 1, 2, 3]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp7;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp7;
            case 8
                shift = shift + [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp8;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp8;
            case 9
                shift = shift + [-4, -3, -2, -1, 0, 1, 2, 3, 4]*J;
                peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) = peakTrue(interp1(fr, 1:length(fr), shift, 'nearest', 'extrap')) + amp*cp9;
                peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) =  peakTrueb(interp1(fshift, 1:length(fshift), shift, 'nearest', 'extrap')) + amp*cp9;
            otherwise % if >9 then use average 
                printf('wrong coupling type');
        end
    end
    peakTrue = peakTrue/max(peakTrue); %normalize
    % figure; plot(fr, peakTrue); jcPlotStyle;
    
    
    %% calculate FID over time
    
    peakTime = zeros(length(fr), length(t)); % frequency domain signal over time
    FID = zeros(1,length(t));
    for j = 1:length(t)
        time = tn(j);
        for i = 1:numHtypes
            amp = H_num(i);
            shift = H_shift(i); % ppm
            cp = H_cp(i); % type
            J = H_J(i); % Hz
            J = J/NMRfreq; % ppm
            sigma = H_sigma(i)/NMRfreq; %ppm
            damping = H_damping(i);
            phase0 = H_phase0(i);
            if cp>9
                cp = meanCP;
            end
            switch cp
                case 1
                    amp = amp*cp1;
                case 2
                    shift = shift + [-0.5, 0.5]*J;
                    amp = amp*cp2;
                case 3
                    shift = shift + [-1, 0, 1]*J;
                    amp = amp*cp3;
                case 4
                    shift = shift + [-1.5, -0.5, 0.5, 1.5]*J;
                    amp = amp*cp4;
                case 5
                    shift = shift + [-2, -1, 0, 1, 2]*J;
                    amp = amp*cp5;
                case 6
                    shift = shift + [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]*J;
                    amp = amp*cp6;
                case 7
                    shift = shift + [-3, -2, -1, 0, 1, 2, 3]*J;
                    amp = amp*cp7;
                case 8
                    shift = shift + [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]*J;
                    amp = amp*cp8;
                case 9
                    shift = shift + [-4, -3, -2, -1, 0, 1, 2, 3, 4]*J;
                    amp = amp*cp9;
                otherwise % if >9 then use average 
                    printf('wrong coupling type');
            end
            % peak = peakCalc(fr, dfr, shift, amp, sigma, damping, time);
            
            p = zeros(1, 4);
            p(2) = sigma;
            p(4) = damping;
            peak = zeros(1, length(fr));
            % dampGauss = @(p, x, t) p(1)/sqrt(2*pi)/p(2)*dfr*exp(-((p(3)-x).^2/2/p(2)^2))*exp(-t/p(4));  % p1, amp; 2, sigma; 3, shift; 4, damping
            for s = 1:length(shift)
                p(1) = amp(s);
                p(3) = shift(s);
                peaki = dampGauss(p, fr, time);
    
                para = zeros(1,3);
                fs = interp1(fr, 1:length(fr), shift(s) - sigma*5, 'nearest', 'extrap'); % 5 sigma
                fe = interp1(fr, 1:length(fr), shift(s) + sigma*5, 'nearest', 'extrap'); % 5 sigma
                for f = fs:fe
                    para(1) = peaki(f);
                    para(2) = fr(f);
                    para(3) = phase0;
                    FID(j) = FID(j) + NMR(para, time);
                end
    
                peak = peak + peaki;
            end
            % NMR = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))); % p1, amp; 2, shift; 3, pahse;
    
    
           peakTime(:, j) = peakTime(:, j) + peak';
        end
    end
    
    
    % peakTime = peakTime/max(peakTime(:));
    FID = FID/max(abs(FID));
    noiseR = randn(1,length(FID), 'like', real(FID))*ynoise;
    noiseI = randn(1,length(FID), 'like', imag(FID))*ynoise;
    FID = FID + noiseR + 1i*noiseI; % add ynoise white noise with equal weight in real and imaginary
    % figure; hold on; plot(fr, peakTime(:,1)); plot(fr, peakTime(:,end)); jcPlotStyle;
    % figure; plot(t, real(FID)); jcPlotStyle;
    
    
    %% FFT
    y = FID;
    yf = y; %[y, zeros(1, length(y)*10)]; % if zero filling add a few times length of zeros after y.
    yfft = fft(yf);
    fs = 1/Ts;
    n = length(yf);
    fshift = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm
    yshift = fftshift(yfft);
    % figure; hold on;
    % plot(fshift,real(yshift)); 
    % plot(fshift, imag(yshift));
    % xlabel('Shift (ppm)')
    % ylabel('Magnitude')
    % jcPlotStyle;
    % figure; plot(fshift, abs(yshift)); jcPlotStyle;
    

    %% save results
    % frb = fshift;
    % peakTrueb = binData(peakTrue, 20);
    peakTrueb = peakTrueb/max(peakTrueb);
    yshift = yshift/max(abs(yshift));
    % figure; plot(frb, peakTrueb); jcPlotStyle;
    
    
    trainingData = FID';
    trainingDataFFT = yshift';
    trainingKey = peakTrue;
    trainingKeyb = peakTrueb;
    trainingKeyTime = peakTime;

figure; plot(fr,trainingKey); hold on; plot(fshift, peakTrueb); ax = gca; ax.XDir = 'reverse'; title('simulated true')

ts = 1:50:length(t);
figure; waterfall(fr, t(ts), trainingKeyTime(:,ts)'); ax = gca; ax.XDir = 'reverse'; 
view(30,45)
ylabel('Time (s)')
xlabel('shift (\delta)')
zlabel('Magnitude')
title('true over time')

figure; plot(fshift, abs(trainingDataFFT)); ax = gca; ax.XDir = 'reverse'; title('FFT of FID')

%-----------
F4.a.f = fr;
F4.a.t = t(ts);
F4.a.y = trainingKeyTime(:,ts)';
F4.a.title = 'true over time';
F4.b.t = t;
F4.b.FID = trainingData;
F4.b.f = fshift;
F4.b.FFT = trainingDataFFT;
F4.b.title = 'FID and FFT';

%% STFT

x = trainingData;

[y,f,t] = stft(x, fs, Window=kaiser(256,5),OverlapLength=220, FFTLength=1024);

y = abs(y);
f = flip(f); %reverse

figure; waterfall(f/NMRfreq,t,squeeze(abs(y))')

ax = gca;
ax.XDir = 'reverse';
view(30,45)

ylabel('Time (s)')
xlabel('shift (\delta)')
zlabel('Magnitude')
title('short-time FT')

%--------
F4.c.f = f/NMRfreq;
F4.c.t = t;
F4.c.y = squeeze(abs(y))';
F4.c.title = 'short-time FT';

%% 2D NMR from STFT  %from Noda's paper: Applied Spectroscopy, 1993, 47, 1329-1336
spec = y';
wave = f/NMRfreq;
x = wave;
time = t;
refspec = mean(spec,1);
figure; plot(x,refspec); title('reference spectrum (average spectrum)')
dynspec = [];
for i = length(time)
    dynspec(i,:) = abs(spec(i,:) - refspec);
end
FFTwave1 = fft(dynspec,[], 1);
FFTwave2 = ifft(dynspec, [], 1);

tdir = [];
for i = 1:length(x)
    for j = 1:length(x)
        tdir(i,j) = 1/pi/length(time)*sum(FFTwave1(:,i).*FFTwave2(:,j),1);
    end
end

figure; imagesc(x, x, real(tdir)); title('2d intensity');
figure; contourf(x, x, real(tdir), 100); title('2d intensity');
figure; surf(x, x, real(tdir), 'EdgeColor', 'none'); title('2d intensity');

tdir2 = real(tdir);
tdir3 = tdir2/max(max(tdir2))*100;

x = wave;
[X,Y] = meshgrid(x,x);
figure; surf(X, Y, tdir3, 'EdgeColor', 'none'); title('2d intensity'); colormap(jet); colormapeditor;
ax = gca; ax.XDir = 'reverse'; ax.YDir = 'reverse';




%% DLSTFT
x = real(trainingData');
dlx = dlarray(x);

[y,f,t] = dlstft(dlx, seconds(1/fs), 'Window', rectwin(256), 'OverlapLength',220, 'FFTLength',256, 'DataFormat','CTB');
y = extractdata(y);
f = extractdata(f);
t = seconds(t);

figure; waterfall(f/NMRfreq,t,squeeze(abs(y))')
ax = gca;
ax.XDir = 'reverse';
view(30,45)
ylabel('Time (s)')
xlabel('shift (\delta)')
zlabel('Magnitude')
title('Deep learning short-time FT')

%-----------
F4.d.f = f/NMRfreq;
F4.d.t = t;
F4.d.y = squeeze(abs(y))';
F4.d.title = 'Deep learning short-time FT';

%% wavelet transform
mtlb = trainingData;
[cfs,frq] = cwt(mtlb, fs, 'TimeBandwidth', 120);
tms = (0:numel(mtlb)-1)/fs;

figure
subplot(2,1,1)
plot(tms,real(mtlb))
axis tight
title("Signal and Scalogram")
xlabel("Time (s)")
ylabel("Amplitude")
subplot(2,1,2)
% surface(tms,frq/NMRfreq,abs(cfs(:,:,1)))
waterfall(frq/NMRfreq, tms, abs(cfs(:,:,1))' + abs(cfs(:,:,2))'); title('wavelet')
ax = gca;
ax.XDir = 'reverse';
view(30,45)
% axis tight
% shading flat
ylabel("Time (s)")
xlabel("shift (\delta)")
% set(gca,"yscale","log")

%-------
data = mtlb;
ts = (0:length(data)-1)*Ts;
fr = ((-10:0.01:10)*NMRfreq)';
tp = 20;
phi = 0;

y = [];
winw = length(data)*Ts/4;
movl = length(data)*Ts/30;
sigma = winw/6; % 2 times sigma each side of the window
[y, f, t] = jcwavelet(data, ts, fr, phi, winw, movl, sigma);

figure; waterfall(f/NMRfreq, t, abs(y)'); title('jcwavelet');
ax = gca;
ax.XDir = 'reverse';

%-------
F4.e.f =  f/NMRfreq;
F4.e.t = t;
F4.e.y = abs(y)';
F4.e.title = 'jcwavelet';



%% plot results
figure; 
subplot(3,2, 1);
waterfall(F4.a.f, F4.a.t, abs(F4.a.y)); ax = gca; ax.XDir = 'reverse'; 
view(30,45)
ylabel('Time (s)')
xlabel('\delta (ppm)')
zlabel('Magnitude')
title('Simulation parameters')

subplot(3,2,2);
plot(F4.b.t, F4.b.FID); title('Simulated FID');
xlabel('time (s)'); ylabel('Magnitude')

subplot(3,2,3);
plot(F4.b.f, abs(F4.b.FFT)); ax = gca; ax.XDir = 'reverse'; title('FFT of FID')
xlabel('Shift (\delta)'); ylabel('Magnitude')

subplot(3,2,4);
waterfall(F4.c.f,F4.c.t, abs(F4.c.y));
ax = gca;
ax.XDir = 'reverse';
view(30,45)
ylabel('Time (s)')
xlabel('\delta (ppm)')
zlabel('Magnitude')
title('Short-time FT')

subplot(3,2,5);
waterfall(F4.d.f,F4.d.t, abs(F4.d.y));
ax = gca;
ax.XDir = 'reverse';
view(30,45)
ylabel('Time (s)')
xlabel('\delta (ppm)')
zlabel('Magnitude')
title(F4.d.title)

subplot(3,2,6)
waterfall(F4.e.f, F4.e.t, abs(F4.e.y)); title('Wavelet');
ax = gca;
ax.XDir = 'reverse';
view(30,45)
ylabel('Time (s)')
xlabel('\delta (ppm)')
zlabel('Magnitude')






