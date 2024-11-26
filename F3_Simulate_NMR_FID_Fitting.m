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
dataPoint = 20; %number of data points
xnoise = 1E-11;
ynoise = 0.01;


dfr = 0.001; % chemical shift step size, (ppm), dfr*NMRfreq (Hz).
fr = -10:dfr:10; %frequency region shift (ppm)


x = 0:dataPoint-1; % time points
%numcycle1ppm = 10;
Ts = 1E-4;% numcycle1ppm/dataPoint/NMRfreq; % unit s. Length of x to satisify ~cycle1ppm cycles for NMRfreq Hz.
x = x*Ts;
xn = x + randn(1,length(x))*xnoise; % add a little noise to time at ps level  

peaks = [1, -1,  1,      0.05, 0.005,  0.01;
         1, 1,   0,      1E8,  0.01, 0; % Amp, freq(ppm), phase,    tau_damp, sigma, diffusion (Gaussian)
         1, 2,   0,      0.15, 0.005,  0.1;
         1, 4,   0,      0.02, 0.01, 1E-1;
         1, 4.1, 2,      0.1,  0.01, 1E-2;
         1, 7,   1,      0.5,  0.005,  1];

%NMR3 = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq+NMRfreq*1E6)*x-p(3)));
NMR3 = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3)));


y = x-x;  % time domain signal
peakf = zeros(length(fr), length(x)); % frequency domain signal over time
parfor (i = 1:length(x),10)
    t = xn(i); %time
    for p = 1:size(peaks,1)
        sigma = peaks(p,5) + sqrt(2*peaks(p,6)*t); % diffusive broaderning sigma = sigma(0) +sqrt(2Dt)
        % peak = peaks(p,1)/sqrt(2*pi)/sigma*dfr*exp(-((fr-abs(peaks(p,2))).^2/2/sigma^2)).*exp(-t/peaks(p,4));
        peak = peaks(p,1)/sqrt(2*pi)/sigma*dfr*exp(-((fr-(peaks(p,2))).^2/2/sigma^2)).*exp(-t/peaks(p,4));
        para = zeros(1,3);

        fs = interp1(fr, 1:length(fr), peaks(p,2)-sigma*5, 'nearest'); % 5 sigma
        fe = interp1(fr, 1:length(fr), peaks(p,2)+sigma*5, 'nearest'); % 5 sigma
        for f = fs:fe
        %for f = 1:length(fr)
            para(1) = peak(f);
            %para(2) = peaks(p,2)/abs(peaks(p,2))*fr(f);
            para(2) = fr(f);
            para(3) = peaks(p,3);
            y(i) = y(i) + NMR3(para, t);
        end
        peakf(:,i) = peakf(:,i) + peak';
    end
end
figure; plot(fr, peakf(:, 1)); hold on; plot(fr, peakf(:, round(length(x)/4))); plot(fr, peakf(:, round(length(x)/2))); plot(fr, peakf(:, end));


noiseR = randn(1,length(x), 'like', real(y))*ynoise;
noiseI = randn(1,length(x), 'like', imag(y))*ynoise;
y = y + noiseR + 1i*noiseI; % add 10% noise white noise with equal weight 


 %% load raw data
 
x = x;
y = y;

yf = y;
%yf = [y, zeros(1,2048-length(y))]; %add zeros for zero padding or filling
yfft = fft(yf);
fs = 1/Ts; % maximum frquency (per data point) in FT
n = length(yf);
fshift1 = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm
yshift1 = fftshift(yfft);
figure; hold on;
plot(fshift1,real(yshift1)); 
plot(fshift1, imag(yshift1));
xlabel('Shift (ppm)')
ylabel('Magnitude')
jcPlotStyle;
figure; plot(fshift1, abs(yshift1)); jcPlotStyle;


yf = [y, zeros(1,2048-length(y))]; %add zeros for zero padding or filling
yfft = fft(yf);
fs = 1/Ts; % maximum frquency (per data point) in FT
f = (0:length(yfft)-1)*fs/length(yfft);
n = length(yf);
fshift2 = (-n/2:n/2-1)*(fs/n)/NMRfreq; % convert frequency Hz to shift in ppm
yshift2 = fftshift(yfft);
figure; hold on;
plot(fshift2,real(yshift2)); 
plot(fshift2, imag(yshift2));
xlabel('Shift (ppm)')
ylabel('Magnitude')
jcPlotStyle;
figure; hold on;  plot(fshift1, abs(yshift1)); plot(fshift2, abs(yshift2));  jcPlotStyle;
figure; hold on;  plot(fshift1, real(yshift1)); plot(fshift2, real(yshift2));  jcPlotStyle;


dr = abs(yshift2)-3; %background correction
dra = [sum(dr(880:960)), sum(dr(1080:1160)), sum(dr(1180:1280)),sum(dr(1400:1440)),sum(dr(1441:1480)),sum(dr(1700:1800))];
dra = dra/dra(1)

% x=[0.25 0.5 1 1.5 2 3 4 6 8];
% y=[19.21 18.15 15.36 14.10 12.89 9.32 7.45 5.24 3.01];
%------------END of loading data: x and y in row vectors--------


%

%% ----Setting up fitting model and parameters-------------
%           the rest can be a double expnential function, any custom function 
%           or a group of functions in a separated Matlab
%           function. Just pass the function handle to the fitting
%           funciton, e.g.
%           function [yfit1, yfit2, yfit3,...] = yournamedoubleExp(para, x1, x2, x3,...)
%                 
%           All functions are numerical solved with no efforts to fit
%           analytically with x and y data.
%-----------------------------------------

% set fitting options
option.maxiteration = 1000;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
option.precision = 1E-5;  % best searching precision, recommend 1 decimal better than desired. e.g want 0.01, set to 0.001.
%option.accuracy = 1E-5;
option.exps = 0.2;
option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.

% ----------------Attn: change below for different fitting equations-----------------
% set the fitting equation to double exponential decay with a base line%
NMR = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))).*exp(-x/p(4));
mdl = @(p, x) NMR(p(1:4), x) + NMR(p(5:8), x) + NMR(p(9:12), x) + NMR(p(13:16), x) + NMR(p(17:20), x) + NMR(p(21:24), x)+ p(end);
% make mdl into a separated function file for earlier versions of Matlab that does not
% support handles of combinded functions.


%  mdl = @(para, x) para(1)*exp(-(x/para(2))); 
% equation grammar: modle name 'mdl' use as function y = mdl(para, x), the rest is the equation.
% you can also generate a function to replace the above equation: 
% function newy = mdl(para, x)
% which will allow combination of equations and global fitting using
% different equations for different pieces of data that share some
% common parameters.

% initial guess
paraGuess = [1, -0.5, 0, 1,   1, 1.5, 2, 1,   1, 2.5, 0, 1,    1, 3.5, 0, 1,  1, 4.5, 0, 1,    1, 7.9, 0, 1,  0];  
%     paraGuess = [1, 33];  % amp, shift,  phase, damping.
% boundarys
bounds = [0, -2, -0.1, 1E-5,    0, 0, -0.1, 1E-5,     0, 1, -0.1, 1E-5,    0, 3, -0.1, 1E-5,       0, 3, -0.1, 1E-5,     0, 6, -0.1, 1E-5,     -0.5; % lower boundary
          2, 0, 6.3, 1E8,       2, 2, 6.3, 1E8,       2, 3, 6.3, 1E8,     2, 5, 6.3, 1E8,          2, 5, 6.3, 1E8,         2, 8, 6.3, 1E8,     0.5]; % upper boundary

%     bounds = [0, 0,;   % lower boundary
%               100, 100]; % upper boundary

%-------------------------------
d1 = paraGuess-bounds(1,:);
d2 = bounds(2,:)-paraGuess;
if prod(d1)*prod(d2)<=0
    display('WARNING: initial guess out of boundary');
end
%--------------END of fitting option setting, equation, initial guess,
%              and guessed parameter boundaries.


%------------------and start fitting:------------------------
% paraGuess = parafinal3;
 
tic
    % paraGuess = parafinal3;
    [paraHist3, parafinal3, paraBounds_95_3, chisq3, rsq3] = jcfit(mdl, x, y, paraGuess, bounds, option);
    % [parafinal3, yfit] = jcfit_LM(mdl, x, y, paraGuess, bounds, option);
% warning: the parameter 95% confidence lower and upper bounds are based on estimation of the local minimum,
% not considering global minimum and other local minima.
toc
%     fprintf(['\n rsq = ', num2str(rsq), '\n']);
% parafinal is the fitted results; yfit is the fittered curve; 
% use residual = y-yfit; to get the residual
% rsq: root mean sqare value best will be close to 1

%--------- plot results -----------------
yfit = mdl(parafinal3, x);
residual = y - yfit;
figure; hold on;  plot(x,real(y),'linewidth',1.5); plot(x,real(yfit),'linewidth',1.5); plot(x, real(residual),'linewidth',1.5);
title(['rsq = ', num2str(rsq3)]); jcPlotStyle;

draf = parafinal3([1,5,9,13,17,21])
% draf = draf/draf(1)

%-------------------------------------
% End. 
% by Jixin Chen @ Ohio University 2024

%}