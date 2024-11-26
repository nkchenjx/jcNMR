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



% %% simulate data
% NMRfreq = 500; % MHz; TMS reference frequency which is set to 0 ppm chemical shift.
% dataPoint = 1024; %number of data points
% xnoise = 1E-11;
% ynoise = 0.001;
% 
% NMRsamplingTrapezoidx = 1:6000; % numerically calculation the acquisition window of each data point
% NMRsamplingTrapezoidx = NMRsamplingTrapezoidx/length(NMRsamplingTrapezoidx)*120E-9; % s; acquisition time 100 ns
% NMRsamplingTrapezoidy = ones(1, 6000); % set the numerical step to be 100 points per two pi of NMRfreqE6 Hz
% xt = 1:500; % make the acquisition window about 90 ns
% TM1 = exp(-(xt-500).^2/2/100^2); NMRsamplingTrapezoidy(xt) = NMRsamplingTrapezoidy(xt).*TM1;
% xt = 1:1000;
% TM2 = exp(-(xt-1000).^2/2/200^2); NMRsamplingTrapezoidy(end+1-xt) = NMRsamplingTrapezoidy(end+1-xt).*TM2;
% NMRsamplingTrapezoidy = NMRsamplingTrapezoidy/0.3351; % correct the acquization amplifier. Obtained from setting shift 1ppm and acquire to match the amplitude expected near time 0.
% AcuisitionPhase = -1.5049; %to correct the phase shift of the acquisition at shift = 1 pmm. Obtained form setting shift to 0 and check at a given acq trapezoid.
% 
% dfr = 0.001; %relative high resolution to get a smooth curve
% fr = -10:dfr:10; %frequency region shift (ppm)
% 
% x = 0:dataPoint-1; % time points
% % numcycle1ppm = 10; Ts = numcycle1ppm/dataPoint/NMRfreq; % unit s. Length of x to satisify ~numcycle1ppm cycles for NMRfreq Hz.
% Ts = 1E-4; % 0.5E-4 window 15 ppm, 1E-4 window 10 ppm
% x = x*Ts;
% xn = x + randn(1,length(x))*xnoise; % add a little noise to time at ps level  
% 
% 
% peaks = [1,  0.5,   0,    0.1; % Amp, shift(ppm),  phase,     tau_damp (s)
%          1  1,    0,    0.1;
%          1,  5,   0,    0.1;
%          1, 8,    0,    0.1];
% 
% 
% %NMR = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))).*exp(-x/p(4));
% NMR1 = @(p, x) p(1)*exp(1i*(2*pi*(p(2)*NMRfreq + NMRfreq*1E6)*x-p(3))).*exp(-x/p(4));
% % NMRfreq*1E6*ts = 2*n*pi  where n is integer such that its frequency
% % signal cancels out over time.
% 
% 
% 
% % phase - 1.885 to correct the correction phase shift e.g. input 1.885 is
% % phase shift  = 0 for shift =1 ppm;
% peaksa = peaks;
% peaksa(:,3) = peaksa(:,3) + AcuisitionPhase;
% ya = zeros(1, length(x));
% for i = 1:length(x)
%     xc = xn(i)+ NMRsamplingTrapezoidx; % acquisition window of each data point
%     yc = NMR1(peaksa(1,1:4), xc) + NMR1(peaksa(2, 1:4), xc) + NMR1(peaksa(3, 1:4), xc) + NMR1(peaksa(4, 1:4), xc);
%     yc = yc.*NMRsamplingTrapezoidy; % true signal times acquisition trapezoid.
%     ya(i) = sum(yc);
% end
% 
% 
% y = NMR1(peaks(1,1:4), xn) + NMR1(peaks(2, 1:4), xn) + NMR1(peaks(3, 1:4), xn) + NMR1(peaks(4, 1:4), xn);
% 
% 
% 
% %%
% data = y';
% ts = (0:length(data)-1)*Ts;
% fr = ((-10:0.01:10)*NMRfreq)';
% tp = 20;
% phi = 0;
% 
% y = [];
% winw = length(data)*Ts/10;  % gaussian wavelet window width
% movl = length(data)*Ts/20; % window moving step size
% sigma = winw/4; % 2 times sigma each side of the window
% [y, f, t] = jcwavelet(data, ts, fr, phi, winw, movl, sigma);
% 
% figure; waterfall(f/NMRfreq, t, abs(y)'); title('jcwavelet');
% ax = gca;
% ax.XDir = 'reverse';


%function [y, f, t] = jcwavelet(data, ts, fr, tp, phi, sigma) % FID, time, frequency series, number of time piece, pahse, gaussian sigma
function [y, f, t] = jcwavelet(data, ts, fr, phi, winw, movl, sigma) %FID, time series (s), freq series (Hz), phase (rad), Gaussian window width (s), moving step (s), Gaussian sigma (s)   
    Ts = ts(2)-ts(1);
    tp = floor((ts(end) - ts(1)-winw)/movl);
    % pl = floor(length(data)/tp); %ind
    hwi= floor(winw/Ts/2); %half windows ind
    mi = floor(movl/Ts);

    s = sigma;
    wl = @(s, f, t, tc) 1/sqrt(2*pi)/s*Ts*exp(-(t-tc).^2/2/s^2).*exp(1i*(2*pi*f*t-phi));
    t = zeros(1,tp);
    for i = 1:tp
        tc = ts(1+ hwi+ mi*(i-1));
        ti = ts(1+mi*(i-1):1+2*hwi+mi*(i-1));
        d = data(1+mi*(i-1):1+2*hwi+mi*(i-1));
        for j = 1:length(fr)
              wlp = wl(s, fr(j), ti, tc);
              y(j,i) = sum(wlp'.*d);
        end
        t(1, i) = tc;
    end
    
    f = flip(fr);
    %t = ((1:tp)*pl-round(pl/2))*Ts;
end

