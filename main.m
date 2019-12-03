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
xlabel('\lambda [nm]')
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
    %     leakage1 = leakage1 +trapz(LEDs(i,1:282));
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
% % t = -5:0.01:5;          %time vector
% f = 1e9;                  %sin frequency
% w = 2*pi*10^9;            %sqrt of resonance freq.
% 
% % Over Damped, Critically Damped, Under Damped
% alpha = [0.1*w w 10*w];              %dampening rate
% syms y(t)
% Dy = diff(y);
% figure()
% for n = 1:3
% ode = diff(y,t,2) == -(2*alpha(n))*Dy - (w^2)*y + sin(2*pi*f*t);
% cond1 = y(0) == 0;
% cond2 = Dy(0) == 1;
% 
% conds = [cond1 cond2];
% ySol = dsolve(ode,conds);
% ySol = simplify(ySol);
% 
% 
% yFunc(t) = ySol(1);
% dt = 0:0.01:10-0.01;
% y_vect = zeros(1,length(dt));
% 
% i = 1;
% for dt_ = dt
%     y_vect(1,i) = double(yFunc(dt_));
%     i = i+1;
% end
% 
% plot(y_vect)
% hold on
% legend ('\alpha_1', '\alpha_2', '\alpha_3')
% xlabel('Time [sec]')
% ylabel('Spectral Bandwidth')
% title('Tone response of laser')
% end

close all
%% Phase 4: Gaussian beam propagation
% lambda4 = 1e-6;
% n = 1.6;        
% w0 = 1e-3;                      % Waist size of Laser beam
% z0 = (pi*(w0^2)*n)/lambda4;     % Rayleigh Range
% z = 0.01:0.01:10;               % Distance
% r =                             % Elevation Distance
% w_z = (w0^2)*(1+(z/z0)^2);      % Radis ar which the field amp. and intensity drop to 1/e and 1/e^2 of axial values
% I0 = 1;                         % Scaled Power Intensity
% 
% I = I0*((w0/w_z)^2)*exp(-2*((r/w_z)^2));


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




%% Phase 7: Modulation
t = 0:0.1:100;
w = 2*pi*8e12;
phi = (0:(2*pi)/100:2*pi)';


MachZender = cos(w*t) + cos(w*t + phi);
figure()

phi2 = 0:(2*pi)/1000:2*pi;
mesh(MachZender)


%% Phase 8: Optical amplifier

Gain_dB = [5 10 20 -10];
GainAmp = 10.^(Gain_dB./10);
n_sp = 5;
AmpNum = 1:100;
Fn_Total = zeros(length(AmpNum),1);
Fn_Total_Atten = zeros(length(AmpNum),1);
figure()
hold on
for i = 1:3
    Fn = 2*((GainAmp(i)-1)*n_sp) / GainAmp(i);
    Fn_Atten = 2*((GainAmp(4)-1)*n_sp) / GainAmp(4);
    Fn_Total(1,1) = Fn;
    for n = 1:length(AmpNum)
        if (n == 1)
            Fn_Total_Atten(n,1) =  Fn_Total(1,1) + (Fn_Atten-1)/((GainAmp(4)));
            continue;
        end
        Fn_Total(n,1) = Fn_Total(n-1,1) + (Fn-1)/((GainAmp(i))^(n-1));
        Fn_Total_Atten(n,1) =  Fn_Total(n,1) + (Fn_Atten-1)/((GainAmp(4)));
    end
    plot(AmpNum,Fn_Total_Atten)
end




