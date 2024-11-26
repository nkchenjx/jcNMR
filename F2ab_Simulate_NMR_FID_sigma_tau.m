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

clear
%  #1, 1 data point peaks --------------------------------------------------------------

%###### adjust amp, freq, sigma, phase and damp values to show their
%different effects on simulations and FFT results.

%% generate raw data of pure peak spike
NMRfreq = 500; % MHz; TMS reference frequency which is set to 0 ppm chemical shift.
dataPoint = 2^11; %number of data points
xnoise = 1E-11;
ynoise = 0.001;


dfr = 0.001; %relative high resolution to get a smooth curve
fr = -10:dfr:10; %frequency region shift (ppm)

x = 0:dataPoint-1; % time points
% numcycle1ppm = 10; Ts = numcycle1ppm/dataPoint/NMRfreq; % unit s. Length of x to satisify ~numcycle1ppm cycles for NMRfreq Hz.
Ts = 1E-4; % 0.5E-4 window 15 ppm, 1E-4 window 10 ppm
x = x*Ts;
xn = x + randn(1,length(x))*xnoise; % add a little noise to time at ps level  


peaks = [1, 1, 0,    0.1, 0.005;] % Amp, shift(ppm), phase,     tau_damp (s), sigma(ppm)  (Gaussian peak damp in exponential decay)
         % 1  1, 0,   0.1, 0.01;
         % 1, 3, 0,   0.1, 0.015;
         % 1, 7, 0,   0.1, 0.02];
% peaks(:,3) = peaks(:,3) + AcuisitionPhase;
%paraTrue = [3, 4, 1, 3,    2, 6, 0, 0.1,   1, 6.5, 0, 2,  0]; 
            % amp, freq, phase, damp

peakf = fr-fr;
for i = 1:size(peaks,1)
    peak = peaks(i,1)/sqrt(2*pi)/peaks(i,5)*dfr*exp(-((peaks(i,2)-fr).^2/2/peaks(i,5)^2));
    sumpeak = sum(peak)
    peakf = peakf + peak;
end
figure; plot(fr, peakf);


%NMR = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq)*x-p(3))).*exp(-x/p(4));
NMR = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))).*exp(-x/p(4));
% NMRfreq*1E6*ts = 2*n*pi  where n is integer such that its frequency

%% just one point data
y = x-x;  % time domain signal
parfor (i = 1:length(x), 12)
    t = xn(i); %time
    para = zeros(1, 4);
    for p = 1:size(peaks,1)
        peak = peaks(p,1)/sqrt(2*pi)/peaks(p,5)*dfr*exp(-((peaks(p,2)-fr).^2/2/peaks(p,5)^2));
        fs = interp1(fr, 1:length(fr), peaks(p,2)-peaks(p,5)*5, 'nearest'); % 5 sigma
        fe = interp1(fr, 1:length(fr), peaks(p,2)+peaks(p,5)*5, 'nearest'); % 5 sigma
        for f = fs:fe
        % for f = 1:length(fr)
            para(1) = peak(f);
            para(2) = fr(f);
            para(3) = peaks(p,3);
            para(4) = peaks(p,4);
            y(i) = y(i) + NMR(para, t);
        end
    end
end

noiseR = randn(1,length(x), 'like', real(y))*ynoise;
noiseI = randn(1,length(x), 'like', imag(y))*ynoise;
y = y + noiseR + 1i*noiseI; % add ynoise white noise with equal weight in real and imaginary


y1 = y;
figure; plot(x, real(y1)); jcPlotStyle; xlim([0,0.12]); ylim([-1.1, 1.1]);

x = x;
y = y1;

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

% figure; plot(fshift, abs(yshift)/max(abs(yshift))); hold on; plot(fr, peakf/max(peakf)); plot(fshift, real(yshift)/max(real(yshift)))
% jcPlotStyle; xlim([0.55, 1.45]); ylim([-0.1, 1.1]);

figure; hold on; plot(fr, peakf/max(peakf)); scatter(fshift, real(yshift)/max(real(yshift)));
jcPlotStyle; xlim([0.81, 1.19]); ylim([-0.1, 1.1]);

%-------------------------------------
    % End. 
    % by Jixin Chen @ Ohio University 2024

