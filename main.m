%%Wireless Communication Networks Hackathon 1
close all
clear variables
clc

%% Phase 1: Designing of LED transmitter
% 1.	Calculate the wave length and optical power of LED with energy bandgap of 1eV,
%       current consumption of 1mA, and quantum efficiency of 0.1

E = 1.60218e-19;        %Energy bandgap 1eV = 1.60218e-19 Joule
i_LED = 1e-3;           %Current consumption 1mA
qe = 0.1;               %Quantum efficiency 0.1
%Wavelength calculation:
h = 6.62607004e-34;     %Planck Constant
c = 3e8;                %Lightspeed
nu = E/h;               %Wave frequency
lambda = c/nu;          %Wavelength
q = 1.60217662e-19;     %electron charge

% 2.	What will be the approximate power density after 1 km if diameter of the LED of 1cm
L = 1e3;                %Distance 1km
r = 0.5e-2;             %Diameter = 1cm => radius = 0.5cm
theta = -atan(r/L):0.1e-6:atan(r/L);
x = L*cos(theta);       %Hypothenus of projection

P0 = (i_LED/q) * qe * h * (c/lambda);
P = P0./(pi*(theta*L).^2);
figure()
plot(theta,P)
title('Power density at distance of 1[km]')
xlabel('\theta [rad]')
ylabel('Power Density')
% ylim([0 10000])
% 3.	Draw the spectra bandwidth of LED as a function of temperature (250-330K),
%       what is the ratio between the maximum and the minimum values

T = 250:330;            %Temp scale 250K-330K
kb = 1.38064852e-23;
%Spectra calculation
delta_E_ph = 3.3*kb*T;
delta_nu = delta_E_ph/h;
delta_lambda = lambda^2 * delta_nu/c;

figure()
plot(T,delta_lambda)
xlabel('Temp [K]')
ylabel('Spectral Bandwidth')
title('Spectra Bandwidth of LED as a function of temprature')
%Ratio calculation
ratio = max(delta_lambda)/min(delta_lambda);

%% Phase 2: The performance of MWSK communication system
% The system includes 4 LED with Gaussian spectral distribution at peak wave length of
% 450nm, 520nm, 530nm, 670nm, and standard deviation of 20nm, 30nm, 40nm, 50nm respectively. M=2
%
lambda_vec = [450e-9 520e-9 530e-9 670e-9]';
sigma_vec = [20e-9 30e-9 40e-9 50e-9]';
M = 2;
% 1.	Draw the system
wavelength_scale = 200e-9:1e-9:1000e-9;
LEDs = zeros(length(lambda_vec),length(wavelength_scale));
for i = 1:length(lambda_vec)
    LEDs(i,:) = (1/(sqrt(2*pi*sigma_vec(i)^2)))*exp(-(wavelength_scale-lambda_vec(i)).^2/(2*sigma_vec(i)^2));
    
end
% LEDs = (1/(sqrt(2*pi*sigma_vec.^2))).*exp(-(wavelength_scale-lambda_vec).^2./(2*sigma_vec.^2));         %TODO: calculate dpectral dist of each gaussian and plot
figure()
title('LED Spectral Distribution')
hold on
xlabel('\lambda [m]')
ylabel('Spectral Density')
for i = 1:4
    plot(wavelength_scale, LEDs(i,:))
end
legend('450nm', '520nm', '530nm', '670nm')

% 2.	Allocate the symbols to specific wavelengths
symbols = ['00' '01' '11' '10'];
% 3.	Find by numerical calculation the bandwidth of ideal filters that minimize the leakage
%       of energy between the symbols

for i = 1:4
    pd(i,:) = makedist('Normal','mu',lambda_vec(i),'sigma',sigma_vec(i));
end
leakage1 = 0;
for i = 2:4
    p = cdf(pd(i), 2e-7:1e-9:4.81e-7);
    leakage1 = leakage1 + p(length(p));
    
end
leakage2 = 0;
for i = [1 3 4]
    p = cdf(pd(i), 2e-7:1e-9:559e-9) - [cdf(pd(i), 2e-7:1e-9:481e-9) zeros(1,559-481)];
    leakage2 = leakage2 + p(length(p));
end
leakage3 = 0;
for i = [1 2 4]
    p = cdf(pd(i), 2e-7:1e-9:591e-9) - [cdf(pd(i), 2e-7:1e-9:469e-9) zeros(1,591-469)];
    leakage3 = leakage3 + p(length(p));    
end
leakage4 = 0;
for i = 1:3
    p = 1 - cdf(pd(i), 2e-7:1e-9:5.955e-7);
    leakage4 = leakage4 + p(length(p));
end


%% Phase 3: Tone response of laser
alpha = [1 2 5] ;
w_o = [5 2 1];
freq = 1*10^9;

t= 0:0.01:500;


%Underdamped solution
syms V(t)
Vt = diff(V);
ode = diff(V,t,2) + diff(0.01*V,t,1) + diff(0.05*V,t,0) == sin(2*pi*(10^9)*t);
cond1 = V(0) == 0;
cond2 = Vt(0) == 0;


conds = [cond1 cond2];
VSol(t) = dsolve(ode,conds);
VSol = simplify(VSol);
VSolFun = matlabFunction(VSol);
figure()
hold on
fplot(VSolFun(t),[0,500]) 

%Critically damped solution
syms V(t)
Vt = diff(V);
ode = diff(V,t,2)+diff(0.025*V,t,1)+ diff(double(1/1600)*V,t,0) == sin(2*pi*(10^9)*t);
cond1 = V(0) == 0;
cond2 = Vt(0) == 0;


conds = [cond1 cond2];
VSol(t) = dsolve(ode,conds);
VSol = simplify(VSol);
VSolFun = matlabFunction(VSol);

fplot(VSolFun(t),[0,500])


%Overdamped solution
syms V(t)
Vt = diff(V);
ode = diff(V,t,2)+diff(0.05*V,t,1)+ diff(0.001*V,t,0) == sin(2*pi*(10^9)*t);
cond1 = V(0) == 0;
cond2 = Vt(0) == 0;


conds = [cond1 cond2];
VSol(t) = dsolve(ode,conds);

VSol = simplify(VSol);
VSolFun = matlabFunction(VSol);

fplot(VSolFun(t),[0,500]) 
title('Solutions of second order ODE for different dampenings')
xlabel('Time [sec]')
ylabel('Amplitude')
legend('Underdamped', 'Critically Damped', 'Overdamped')

%% Phase 4: Gaussian beam propagation
z_vec = 0:600;
lambda4 = 1e6;
n = 1.6;
w0=1*10^-3;
z0 = pi*w0^2/lambda;
I0 =1;

legend_arg = [];
figure();
for constant = 100:400:1700

    R = @(z) ((0.5*(((w0^2).*(1+(z./z0).^2)).^0.5).^2).*(log(constant*((((w0^2).*(1+(z./z0).^2)).^0.5).^2).*(w0^-2)))).^(0.5);
    plot(z_vec, R(z_vec));
    hold on;
end
title('R(z) for I(R,Z) = Constant');
legend('100','500','900','1300','1700');
xlabel('z [m]')
ylabel('Constant ratio')


%% Phase 5: Lamebrain light source
theta = -pi/2:(pi/2)/100:pi/2;
C1 = 1;
C2 = 0;
figure()
hold on
for C3 = 1:4
    I = C1 * (cos(abs(theta)-C2)).^C3;
    plot(theta,I)
end
legend('C3 = 1', 'C3 = 2', 'C3 = 3', 'C3 = 4')
xlabel('\theta [rad]')
ylabel('Power Density')

%% Phase 6: VLC system

x = -4:4;
y = -4:4;
height = 3;

[X, Y] = meshgrid(x,y);

distances = sqrt(X.^2 + Y.^2);
angles = atan(distances/height);

Noise = 0;

for i = 1:9
   for j = 1:9
      Noise = Noise + cos(abs(angles(i,j))) ;
   end
end
Noise = Noise - cos(abs(angles(5,5)));

SNIR = 1/Noise;

%% Phase 7: Modulation
t = (0:0.1:100);
w = 2*pi*8e12;
phi = (0:(2*pi)/100:2*pi)';


MachZender = cos(w*t) + cos(w*t + phi);
figure()
mesh(t,phi,MachZender)
title('Mach-Zehnder modulator output as a function of phase')
xlabel('Time [sec]')
ylabel('Phase [rad]')
zlabel('Modulator Amplitude')



%% Phase 8: Optical amplifier

Gain_dB = [5 10 20 -10];
GainAmp = 10.^(Gain_dB./10);
n_sp = 5;
AmpNum = 1:10;
Fn_Total = zeros(length(AmpNum),1);
Fn_Total_Atten = zeros(length(AmpNum),1);
figure()
hold on
for i = 1:3
    Fn = 2*n_sp;
    Fn_Atten = 2*n_sp;
    Fn_Total(1,1) = Fn;
    Fn_Total_Atten(1,1) =  Fn_Total(1,1) + (Fn_Atten-1)/((GainAmp(4)));
    for n = 2:length(AmpNum)
        Fn_Total(n,1) = Fn_Total(n-1,1) + (Fn-1)/((GainAmp(i))^(n-1));
        Fn_Total_Atten(n,1) =  Fn_Total(n,1) + (Fn_Atten-1)/((GainAmp(4)));
    end
    plot(AmpNum,Fn_Total_Atten)
    
end

title('Noise Figure as a function of number of units')
xlabel('Number of units (Amplifiers + Attenuator)')
ylabel('Noise Figure')
legend('Gain = 5dB', 'Gain = 10dB', 'Gain = 20dB')



