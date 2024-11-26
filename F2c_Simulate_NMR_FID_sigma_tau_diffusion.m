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


%% an example for NMR decay
% this example takes 100 s on a Intel i7 cpu. No GPU or parallel funciton has been used yet.

clear

%  #3, diffusion Gaussion Peaks------------------------------------------------------------------------------------------

%% generate raw data of continuous peak
NMRfreq = 500; % MHz
dataPoint = 2048; %number of data points
xnoise = 1E-11;
ynoise = 0.01;


dfr = 0.001; %relative high resolution to get a smooth curve
fr = -10:dfr:10; %frequency region shift (ppm)


x = 0:dataPoint-1; % time points
%numcycle1ppm = 10;
Ts = 1E-4;% numcycle1ppm/dataPoint/NMRfreq; % unit s. Length of x to satisify ~cycle1ppm cycles for NMRfreq Hz.
x = x*Ts;
xn = x + randn(1,length(x))*xnoise; % add a little noise to time at ps level  

peaks = [1, 1,  0,    1E8,  0.005, 0.00001]; % Amp, freq(ppm), phase,    tau_damp, sigma, diffusion (Gaussian)
         % 1, -1,  0,   6,  0.1, 0.0001;
         % 1, 2,  0,    4,  0.01, 0.001;
         % %1, -2,  0,   4,  0.01, 0.001;
         % 1, 4,  0,    6,  0.01, 0.01;
         % 1, -4,  0,   6,  0.01, 0.01;
         % 1, 7,  0,    10, 0.01, 0.00001];

%NMR3 = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq+NMRfreq*1E6)*x-p(3)));
NMR3 = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3)));


y = x-x;  % time domain signal
peakf = zeros(length(fr), length(x)); % frequency domain signal over time
parfor (i = 1:length(x),10)
    t = xn(i); %time
    for p = 1:size(peaks,1)
        sigma = peaks(p,5) + sqrt(2*peaks(p,6)*t); % diffusive broaderning sigma = sigma(0) +sqrt(2Dt)
        peak = peaks(p,1)/sqrt(2*pi)/sigma*dfr*exp(-((fr-abs(peaks(p,2))).^2/2/sigma^2)).*exp(-t/peaks(p,4));
        para = zeros(1,3);
        fs = interp1(fr, 1:length(fr), peaks(p,2)-sigma*5, 'nearest'); % 5 sigma
        fe = interp1(fr, 1:length(fr), peaks(p,2)+sigma*5, 'nearest'); % 5 sigma
        for f = fs:fe
        % for f = 1:length(fr)
            para(1) = peak(f);
            para(2) = peaks(p,2)/abs(peaks(p,2))*fr(f);
            para(3) = peaks(p,3);
            y(i) = y(i) + NMR3(para, t);
        end
        peakf(:,i) = peakf(:,i) + peak';
    end
end
% figure; plot(fr, peakf(:, 1)); hold on; plot(fr, peakf(:, 100)); plot(fr, peakf(:, 1000));


noiseR = randn(1,length(x), 'like', real(y))*ynoise;
noiseI = randn(1,length(x), 'like', imag(y))*ynoise;
y = y + noiseR + 1i*noiseI; % add 10% noise white noise with equal weight 


%% load raw data

yf = y; %add zeros for zero padding or filling
yfft = fft(yf);
fs = 1/Ts; % maximum frquency (per data point) in FT
n = length(yf);
fshift = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm
yshift = fftshift(yfft);
% figure; hold on;
% plot(fshift,real(yshift)); 
% plot(fshift, imag(yshift));
% xlabel('Shift (ppm)')
% ylabel('Magnitude')
% jcPlotStyle;
% figure; plot(fshift, abs(yshift)); 



figure; plot(x, y); jcPlotStyle; xlim([0,0.12]); ylim([-1.1, 1.1]);

figure; hold on; scatter(fshift, real(yshift)/max(real(yshift)), 'o')
plot(fr, peakf(:,1)/max(peakf(:,1))); 
plot(fr, peakf(:,100)/max(peakf(:,100))); 
plot(fr, peakf(:,200)/max(peakf(:,200))); 
plot(fr, peakf(:,500)/max(peakf(:,500)));
plot(fr, peakf(:,1000)/max(peakf(:,1000)));
jcPlotStyle; xlim([0.94, 1.06]); ylim([-0.1, 1.1]);


