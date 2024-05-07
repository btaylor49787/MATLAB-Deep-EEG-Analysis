
%**************************************************************************
%	EEG_Analysis.m, Created on 01/23/2024
%	Project – Deep EEG Analysis
%   Marshall University
%   Ben Taylor
%**************************************************************************

clc; clear; clf;

fs = 5000;                  % sampling frequency
Ts = 1/fs;                  % sampling interval
fm = 50;                    % frequency of interest
Tp = 1;
T_pulse = 0.005;            % pulse duration
f_pulse = 1/T_pulse;        % pulse frequency   
R0 = 0.1;                   % model's radius
sigma = 0.1;                % conductivity

g1_x = R0;
g1_y = R0;

Npt = fs*Tp;
N_pulse = T_pulse/Ts;
N_seg = floor(Npt/N_pulse);
Tm = (0:1:(Npt-1))/fs;

x_pulse = zeros(1,N_pulse);

x_mes = zeros(1,Npt);
y_mes = zeros(1,Npt);

% *************************************************************************
% 	Generate pulse
% *************************************************************************

for i=1:N_pulse
    dummy1 = sin(2*pi*f_pulse*Tm(i));
    dummy2 = cos(2*pi*f_pulse*Tm(i));
    if dummy1 <= 0.0
        x_pulse(i) = 0.0;
    elseif dummy2 < 0.0
        x_pulse(i) = -1.0;
    else
        x_pulse(i) = 1.0;
    end
end

% *************************************************************************
% 	Generate x_mes data
% *************************************************************************

flag = round(round(10*rand(1,N_seg))/10);

for i=1:N_seg
    for j=1:N_pulse
        x_mes((i-1)*N_pulse+j) = x_pulse(j)*flag(i);
    end
end

% *************************************************************************
% 	Filter Data
% *************************************************************************

for i=1:Npt
    y_mes(i) = (1/(4*pi*0.1))*(x_mes(i)/R0);
end
 
% *************************************************************************
% 	Mapping Generation
% *************************************************************************

Mpt =100;
v_map = zeros(Mpt,Mpt);
r = zeros(Mpt,Mpt);

x_map = (1:1:Mpt)*2*R0/Mpt;
y_map = (1:1:Mpt)*2*R0/Mpt;

for i=1:Mpt
    for j=1:Mpt
        r(i,j) = sqrt((x_map(j)-g1_x)^2+(y_map(i)-g1_y)^2);
        if r(i,j) <= R0
            v_map(i,j) = (1/(4*pi*sigma))/r(i,j);
        end
    end
end

% *************************************************************************
% 	Random Points Generation
% *************************************************************************

% Radius of 1cm circle
R1 = .01; 

%// running variable
t = linspace(0,2*pi,fs);

% generation of 1cm circle
x_1 = 0 + R1*sin(t);
y_1 = 0 + R1*cos(t);

% generation of 10cm circle
x_10 = 0 + R0*sin(t);
y_10 = 0 + R0*cos(t);


% generates random coordinates on edge of 10cm circle
edge = [];
n_edge = 10;  % number of random points on edge of circle
theta = rand(1,n_edge)*2*pi;
xp = R0*cos(theta)+0;
yp = R0*sin(theta)+0;
edge = [edge; xp ;yp];
edge = transpose(edge);

i2=0;
pts=[];
% generates 10 random coordinates within a 1cm circle
while i2<10
% Create a random point vector
    pt = [R1*2*(rand(fs,1)-.5) R1*2*(rand(fs,1)-.5)];
% Select the first pair that lies inside circle
    pt = pt(find(sqrt( pt(:,1).^2 + pt(:,2).^2 )<R1,1),:);
    pts = [pts; pt]
    i2 = i2+1;
end

% Distance between the random point and edge of circle
diff=pts-edge;
sq=diff.^2;
s=sum(sq,2);
dist=sqrt(s)

% *************************************************************************
% 	Plot data
% *************************************************************************

figure(1);
subplot(211);
plot(Tm,x_mes);
title('x__mes: Inner Signal');
subplot(212);
plot(Tm,y_mes);
title('y__mes: Outer Signal');

figure(2);
plot(x_mes);
ylim([-1.1 1.1])

figure(3);
mesh(v_map);

disp('program finished')