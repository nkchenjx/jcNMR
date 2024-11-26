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
dataPoint = 2048; %number of data points
xnoise = 1E-11;
ynoise = 0.001;

dfr = 0.001; %relative high resolution to get a smooth curve
fr = -10:dfr:10; %frequency region shift (ppm)

x = 0:dataPoint-1; % time points
% numcycle1ppm = 10; Ts = numcycle1ppm/dataPoint/NMRfreq; % unit s. Length of x to satisify ~numcycle1ppm cycles for NMRfreq Hz.
Ts = 1E-4; % 0.5E-4 window 15 ppm, 1E-4 window 10 ppm
x = x*Ts;
xn = x + randn(1,length(x))*xnoise; % add a little noise to time at ps level  


peaks = [1, 0.5,   0,    0.1; % Amp, shift(ppm),  phase,     tau_damp (s)
         1  1,    0,    0.1;
         1,  5,   0,    0.1;
         1, 8,    0,    0.1];

%NMR = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))).*exp(-x/p(4));
NMR1 = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))).*exp(-x/p(4));
NMR1a = @(p, x) p(1)*exp(-1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))).*exp(-x/p(4));
% NMRfreq*1E6*ts = 2*n*pi  where n is integer such that its frequency
% signal cancels out over time.

%% 1, one point intensity
y1 = NMR1(peaks(1,1:4), xn) + NMR1(peaks(2, 1:4), xn) + NMR1(peaks(3, 1:4), xn) + NMR1(peaks(4, 1:4), xn);

NMRsamplingTrapezoidx = 1:5000; % numerically calculation the acquisition window of each data point
NMRsamplingTrapezoidx = NMRsamplingTrapezoidx/length(NMRsamplingTrapezoidx)*100E-9; % s; acquisition time 100 ns
peaksa = peaks;
i = 1;
xc = xn(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);



figure; plot(xc, yc); jcPlotStyle; title('point #1');

i = 10;
xc = xn(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
figure; plot(xc, yc); jcPlotStyle; title('point #10');



figure; plot(x,y1); title('1-point');jcPlotStyle;

%FFT
x = x;
y = y1;
noiseR = randn(1,length(x), 'like', real(y))*ynoise;
noiseI = randn(1,length(x), 'like', imag(y))*ynoise;
y = y + noiseR + 1i*noiseI; % add ynoise white noise with equal weight in real and imaginary

yf = y; %[y, zeros(1, length(y)*10)]; % if zero filling add a few times length of zeros after y.
yfft = fft(yf);
fs = 1/Ts;
n = length(yf);
fshift = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm
yshift = fftshift(yfft);
figure; hold on;
plot(fshift,real(yshift)); 
plot(fshift, imag(yshift));
xlabel('Shift (ppm)')
ylabel('Magnitude')
jcPlotStyle;
figure; plot(fshift, abs(yshift)); jcPlotStyle

fshift1 = fshift;
yshift1 = yshift;

%

%% 2, Trapezoidx with Gaussian tails

NMRsamplingTrapezoidx = 1:6000; % numerically calculation the acquisition window of each data point
NMRsamplingTrapezoidx = NMRsamplingTrapezoidx/length(NMRsamplingTrapezoidx)*120E-9; % s; acquisition time 100 ns
NMRsamplingTrapezoidy = ones(1, 6000); % set the numerical step to be 100 points per two pi of NMRfreqE6 Hz
xt = 1:500; % make the acquisition window about 90 ns
TM1 = exp(-(xt-500).^2/2/100^2); NMRsamplingTrapezoidy(xt) = NMRsamplingTrapezoidy(xt).*TM1;
xt = 1:1000;
TM2 = exp(-(xt-1000).^2/2/200^2); NMRsamplingTrapezoidy(end+1-xt) = NMRsamplingTrapezoidy(end+1-xt).*TM2;
NMRsamplingTrapezoidy = NMRsamplingTrapezoidy/0.3351; % correct the acquization amplifier. Obtained from setting shift 1ppm and acquire to match the amplitude expected near time 0.
AcuisitionPhase = -1.5049; %to correct the phase shift of the acquisition at shift = 1 pmm. Obtained form setting shift to 0 and check at a given acq trapezoid.


% phase - 1.885 to correct the correction phase shift e.g. input 1.885 is
% phase shift  = 0 for shift =1 ppm;
peaksa = peaks;
peaksa(:,3) = peaksa(:,3) + AcuisitionPhase;
ya1 = zeros(1, length(x));
for i = 1:length(x)
    xc = x(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
    yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
    yc = yc.*NMRsamplingTrapezoidy; % true signal times acquisition trapezoid.
    ya1(i) = sum(yc);
end

i = 1;
xc = x(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
yc2 = yc.*NMRsamplingTrapezoidy; % true signal times acquisition trapezoid.
ya1(i) = sum(yc2);

figure; plot(xc, yc); hold on; plot(xc, NMRsamplingTrapezoidy); hold on; plot(xc, yc2);jcPlotStyle; title('point #1')


i = 10;
xc = x(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
yc2 = yc.*NMRsamplingTrapezoidy; % true signal times acquisition trapezoid.
ya1(i) = sum(yc2);

figure; plot(xc, yc); hold on; plot(xc, NMRsamplingTrapezoidy); hold on; plot(xc, yc2);jcPlotStyle; title('point #10')




figure; plot(x,ya1);jcPlotStyle


% FFT
x = x;
y = ya1;
noiseR = randn(1,length(x), 'like', real(y))*ynoise;
noiseI = randn(1,length(x), 'like', imag(y))*ynoise;
y = y + noiseR + 1i*noiseI; % add ynoise white noise with equal weight in real and imaginary

yf = y; %[y, zeros(1, length(y)*10)]; % if zero filling add a few times length of zeros after y.
yfft = fft(yf);
fs = 1/Ts;
n = length(yf);
fshift = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm
yshift = fftshift(yfft);
figure; hold on;
plot(fshift,real(yshift)); 
plot(fshift, imag(yshift));
xlabel('Shift (ppm)')
ylabel('Magnitude')
jcPlotStyle;
figure; plot(fshift, abs(yshift)); jcPlotStyle

fshift2 = fshift;
yshift2 = yshift;

%

%% 3, Trapezoidx 2
NMRsamplingTrapezoidx = 1:5000; % numerically calculation the acquisition window of each data point
NMRsamplingTrapezoidx = NMRsamplingTrapezoidx/length(NMRsamplingTrapezoidx)*100E-9; % s; acquisition time 100 ns
NMRsamplingTrapezoidy = ones(1, 5000); % set the numerical step to be 100 points per two pi of NMRfreqE6 Hz
xt = 1:50; % make the acquisition window about 90 ns
TM1 = xt/50; NMRsamplingTrapezoidy(xt) = NMRsamplingTrapezoidy(xt).*TM1;
xt = 1:100;
TM2 = xt/100; NMRsamplingTrapezoidy(end+1-xt) = NMRsamplingTrapezoidy(end+1-xt).*TM2;
NMRsamplingTrapezoidy = NMRsamplingTrapezoidy/10.135; % correct the acquization amplifier. Obtained from setting shift 1ppm and acquire to match the amplitude expected near time 0.
AcuisitionPhase = -3.14; %to correct the phase shift of the acquisition at shift = 1 pmm. Obtained form setting shift to 0 and check at a given acq trapezoid.


% phase - 1.885 to correct the correction phase shift e.g. input 1.885 is
% phase shift  = 0 for shift =1 ppm;
peaksa = peaks;
peaksa(:,3) = peaksa(:,3) + AcuisitionPhase;
ya2 = zeros(1, length(x));
for i = 1:length(x)
    xc = x(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
    yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
    yc = yc.*NMRsamplingTrapezoidy; % true signal times acquisition trapezoid.
    ya2(i) = sum(yc);
end

i = 1;
xc = x(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
yc2 = yc.*NMRsamplingTrapezoidy; % true signal times acquisition trapezoid.
ya2(i) = sum(yc2);

figure; plot(xc, yc); hold on; plot(xc, NMRsamplingTrapezoidy); hold on; plot(xc, yc2);jcPlotStyle; title('point #1')


i = 10;
xc = x(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
yc2 = yc.*NMRsamplingTrapezoidy; % true signal times acquisition trapezoid.
ya2(i) = sum(yc2);

figure; plot(xc, yc); hold on; plot(xc, NMRsamplingTrapezoidy); hold on; plot(xc, yc2);jcPlotStyle; title('point #10')



% FFT
x = x;
y = ya2;
 
yf = y; %[y, zeros(1, length(y)*10)]; % if zero filling add a few times length of zeros after y.
yfft = fft(yf);
fs = 1/Ts;
n = length(yf);
fshift = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm
yshift = fftshift(yfft);
figure; hold on;
plot(fshift,real(yshift)); 
plot(fshift, imag(yshift));
xlabel('Shift (ppm)')
ylabel('Magnitude')
jcPlotStyle;
figure; plot(fshift, abs(yshift)); jcPlotStyle

fshift3 = fshift;
yshift3 = yshift;


%%

figure; plot(x,y1, LineWidth=2); hold on; plot(x,ya1, LineWidth=1.5); plot(x, ya2, LineWidth=1); jcPlotStyle

figure; hold on; 
plot(fshift1, abs(yshift1)); 
plot(fshift2, abs(yshift2)); 
plot(fshift3, abs(yshift3)); 
jcPlotStyle



    %-------------------------------------
    % End. 
    % by Jixin Chen @ Ohio University 2024
%}