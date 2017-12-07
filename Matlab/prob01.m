function []=prob01()
% Problem 1
% MOLAR VOLUME AND COMPRESSIBILITY FACTOR FROM VAN DER WAALS EQUATION

% Mohammad Rahmani
% Chemical Engineering Department
% Amirkabir University of Technology
% m.rahmani@aut.ac.ir
% Thursday, Dec 7th, 2017

%% 0. Short description
% Calculate the molar volume and compressibility factor for gaseous ammonia at a pressure
% (a) P = 56 atm and a temperature T = 450 K using the van der Waals equation of state.
% (b) Repeat the calculations for the following reduced pressures: Pr = 1, 2, 4, 10, and 20.
% (c) How does the compressibility factor vary as a function of Pr.?
%
%
% Ref:
% [1]. Cutlip, M.B., et al., A collection of 10 numerical problems in chemical
%      engineering solved by various mathematical software packages.
%      Computer Applications in Engineering Education, 1998. 6(3): p. 169-180.
%      DOI: 10.1002/(SICI)1099-0542(1998)6:3<169::AID-CAE6>3.0.CO;2-B

%% 1. Problem data
%set the constnats
R  = 0.08206;		% general gas constant, atm.liter/g-mol.K

% Ammonia critical pressure and temperature
Pc = 111.3;         % Ammonia critical pressure, atm
Tc = 405.5;         % Ammonia critical temperature, K

% Calculate the van der Waals EOS constants for ammonia
a = 27/64*(R^2*Tc.^2/Pc); 
b = R*Tc/(8*Pc);

% Operating data
T  = 450;                   % Temp, K
Pr = [56/Pc 1 2 4 10 20];   % Reduced pressure, unitless

%% 2. Calculate the molar volume using van der Waals EOS
% The molar volume for ammonia is calculated using fzero
% for all pressures
n = length(Pr);
Vm = zeros(1,n);    % Molar volume, liter/g-mol
Z  = zeros(1,n);    % Gas compressibility, unitless
for j=1: n
    P = Pr(j)*Pc;
    Vguess = R*T/P;     % Initial guess, using the ideal gas law
    % Solve van der Waals equation for molar volume
    Vm(j)=fzero(@vander_Waals,Vguess);
    Z(j) = P*Vm(j)/(R*T);   % Compressibility, unitless
end

%% 3. Visualization
% 3.1. Display results 
disp(' Reduced pressure,  Molar volume,  Compressibility factor ')
for j=1:n
    fprintf('%7.4f  %18.5f %14.5f \n',Pr(j), Vm(j), Z(j))
end

% 3.2. Plot compressibility factor and molar volume versus reduced pressure
subplot(2,1,1)
plot(Pr,Z,'ro-')
title('Compressibility factor vs Reduced pressure')
ylabel('Compressibility factor')
xlabel('Reduced pressure')

subplot(2,1,2)
plot(Pr,Vm,'b*-')
title('Molar volume vs Reduced pressure')
ylabel('Molar volume, liter/g-mol')
xlabel('Reduced pressure')



%% Nested functions

    function res=vander_Waals(V)
        % implements the van der Waals equation of state
        % in form of f(v)=0
        res=P*V.^3 - P*b*V.^2 - R*T*V.^2 + a*V-a*b;
    end

% End of prob01 function
end