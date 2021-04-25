close all
clear 
clc

%% Dati
Fs = 250; % [Hz]

[idx, annotations] = textread ('sel33.q1c','%d %c');

signal = load('sel33.dat','-ascii');

ECG = signal(:,1);

%% 
% Find instants marked 't'

T = find(annotations == 't');

figure('Name','T waves')
plot(ECG)
hold on
plot(idx(T),ECG(idx(T)),'*')
title 'T waves highlighted on the entire signal'

figure
plot(ECG)
hold on
plot(idx(T),ECG(idx(T)),'*')
xlim([idx(1)-Fs, idx(end)+Fs])
title 'T-waves'

%% Pre-processing
% High-pass filtering (2 Hz)
HP = load('HP.mat');
y = filter(HP.Num,1,ECG);
delay1 = floor(length(HP.Num)/2);  % FIR filters has group delay
y(1:delay1) = [];                  % Correct the delay         
    
% Low-pass filtering (30 Hz)
LP = load('LP.mat');
y = filter(LP.Num,1,y);
delay2 = floor(length(LP.Num)/2);
y(1:delay2) = [];
    
% Be 'xm' the minimun of the derivative of the signal, and after the maximum
% of a positive T-wave, searched in a window of 200 ms after T-peak
 derivative = filter([1,-1],2,y);

 figure('Name','Derivative of the signal')
 plot(y)
 hold on
 plot(derivative)
 xlim([idx(1)-Fs, idx(end)+Fs])
 legend('Filtered signal','Derivative of signal')

window = (200e-3)*Fs;
minimum = zeros(1, length(T));
xm = zeros(1, length(T));
for i=1:length(T)
    [minimum(i),xm(i)] = min(derivative(idx(T(i)):idx(T(i))+window));
    xm(i) = xm(i) + idx(T(i));
end

figure('Name','xm')
plot(y)
hold on
plot(xm,y(xm),'g*')
xlim([idx(1)-Fs, idx(end)+Fs])
title 'x_m'

% Be 'xr' a point in the range between 200 and 400 ms after T peak, in the
% isoelectric segment, where the derivative is under a certein thresholde
% (near zero). If there is no a point that satisfy this requirement, choose
% the central point in this range.

threshold = 0.1;

xr = zeros(1, length(T));
for i =1:length(T)
    v = find(derivative(idx(T(i))+window : idx(T(i))+2*window) < threshold);
       % v contains all element under the threshold in the considered range.

    if sum(size(v)) == 0   % If v = [], then size(v) = [0 0] 
                           % (no element under the threshold)
       xr(i) = idx(T(i)) + 1.5*window;
    else
       xr(i) = v(length(v)); % If more than one element are under threshold, 
                             % I consider the last one.
       xr(i) = xr(i) + idx(T(i))+ window;  % Note: if xr(i) = [] ----> xr(i) + something = [] anyway.
    end
end

yr = y(xr);

figure('Name','xr')
plot(y)
hold on
plot(xr,yr,'r*')
xlim([idx(1)-Fs, idx(end)+Fs])
title 'x_r'

figure('Name','Range in which I''m looking for xi')
plot(y)
hold on
plot(xm,y(xm),'g*')
hold on
plot(xr,yr,'r*')
xlim([idx(1)-Fs, idx(end)+Fs])
title 'Range in which I''m looking for xi'
     
% Compute the trapezoidal area between each 'xm' and 'xr'   

%        xm,y(xm)       xend,y(xm)                xr,y(xm)
%          .________________.________________________. 
%           \                                        |
%            \                                       |
%              \                                     |
%                \                                   |
%                   \                                |
%                      \                             |
%                         \.________________________.|
%                      xend,y(xend)                    xr,y(xend)

% To find [xend,y(xend)] I consider a generic 'xi' in the range xm->xr to
% find the associated areas, and then I consider the maximum area to find
% the point xi=xend

B = xr - xm;
area = zeros(length(T), max(B));
for j = 1: length(T)
    for xi = xm(j) : xr(j)
        
        b = xr(j) - xi;
        h = y(xm(j))-y(xi);
         
        area(j,xi-xm(j)+1) = 0.5*(B(j)+b)*h;
    
    end
end
    
figure('Name','Area')
plot (area')
title 'Area'
   
%%
% I can finally find the final point of T-wave (xend) as the maximum of
% each area:

[maximum,xend] = max(area,[],2); % Maximum value of each row
for i = 1:length(T)
    xend(i) = xend(i) + xm(i);
end
xend = xend';

figure('Name','Identified points')
plot(y)
hold on
plot(idx(T),y(idx(T)),'k*')
hold on
plot(xm,y(xm),'g*')
hold on
plot(xend,y(xend),'c*')
hold on
plot(xr,y(xr),'r*')
xlim([idx(1)-Fs, idx(end)+Fs])
title 'Identified points'
legend('','T peak','x_m','x_e_n_d','x_r')