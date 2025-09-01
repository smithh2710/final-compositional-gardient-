comp = [0.0042,0.0069,0.5035,0.07849,0.06769,0.0104,0.032,0.0116,0.0155,0.0188,0.035,0.0375,0.0228,0.03264,0.02593,0.02926,0.01459,0.01646,0.01166,0.01237,0.00958,0.0063]';
M_wgt = [28.01,44.01 ,16.04 ,30.07 ,44.10 ,58.12 ,58.12 ,72.15 ,72.15 ,86.18 ,96 ,107 ,121 ,140.09 ,167.57 ,204.75 ,243.6 ,275.27 ,317.02 ,370.39 ,456.83 ,640.76]';
acentric = [0.0400 ,0.2250 ,0.0080 ,0.0980 ,0.1520 ,0.1760 ,0.1930 ,0.2270 ,0.2510 ,0.2960 ,0.3375 ,0.3743 ,0.4204 ,0.4833 ,0.5703 ,0.6847 ,0.7947 ,0.8805 ,0.9836 ,1.0983 ,1.2292 ,1.1596]' ;
Tc = [-146.95 ,31.05 ,-82.55 ,32.25 ,96.65 ,134.95 ,152.05 ,187.25 ,196.45 234.25 ,301.63 ,323.97 ,349.62 ,381.84 ,422.96 ,473.01 ,519.43 ,555.64 ,600.26 ,654.86 ,737.85 ,912.38]' + 273.15 ;
Pc = [33.94 ,73.76 ,46.00 ,48.84 ,42.46 ,36.48 ,38 ,33.84 ,33.74 ,29.69 ,29.59 ,27.37 ,25.01 ,22.59 ,20.12 ,17.91 ,16.4 ,15.55 ,14.71 ,13.96 ,13.16 ,12.27 ]'  * 1e5 ; 
comp = comp ./ sum(comp) ; 
M = M_wgt/1000;

% comp = [0.42, 0.69, 50.04, 7.85, 6.77, 1.04, 3.2, 1.16, 1.55, 1.88, 3.50 , 3.75, 2.28, 3.26, 2.59, 2.93 ,1.46 ,1.65 ,1.17 ,1.24 ,0.96 ,0.63]';
% comp = comp ./ sum(comp) ; 
% M_wgt = [28.0 ,44.0 ,16.0 ,30.1 ,44.1 ,58.1 ,58.1 ,72.2 ,72.2 ,86.2 ,96.0 ,107 ,121 ,140.1 ,167.6, 204.7 ,243.6 ,275.3 ,317.0 ,370.4 ,456.8 ,640.9]';
% 
% M = M_wgt/1000;
% Tc = [-147 ,31.1 ,-82.5 ,32.2 ,96.6 ,134.9 ,152.1 ,187.2 ,196.4 234.2 ,269.9 ,290.9 ,315.3 ,345.9 ,385 ,432.7 ,476.8 ,511.1 ,553.3 ,604.9 ,683.3 ,847.8]' + 273.15 ;
% Pc = [33.94 ,73.76 ,46.00 ,48.84 ,42.46 ,36.48 ,38 ,33.84 ,33.74 ,29.69 ,29.60 ,27.86 ,25.83 ,23.76 ,21.58 ,19.59 ,18.24 ,17.54 ,16.84 ,16.23 ,15.62 ,14.72 ]'  * 1e5 ; 
% acentric = [0.0400 ,0.2250 ,0.0080 ,0.0980 ,0.1520 ,0.1760 ,0.1930 ,0.2270 ,0.2510 ,0.2960 ,0.338 ,0.374 ,0.420 ,0.483 ,0.570 ,0.685 ,0.795 ,0.881 ,0.984 ,1.098 ,1.229 ,1.159]' ;

n = length(comp);
BIP = zeros(n,n); 

BIP(1,2) = -0.017    ; BIP(2,1) = -0.017 ; 
BIP(1,3) = 0.0311 ; BIP(3,1)= 0.0311;
BIP(1,4) = 0.0515 ; BIP(4,1) = 0.0515 ; 
BIP(1,5) = 0.0852 ; BIP(5,1)= 0.0852 ; 
BIP(1,6) = 0.1033 ; BIP(6,1)=0.1033; 
BIP(1,7) = 0.0800 ; BIP(7,1) = 0.0800;
BIP(1,8)= 0.0922 ; BIP(8,1) = 0.0922 ; 
BIP(1,9)= 0.1000 ; BIP(9,1)=0.1000;
BIP(1,10:n)= 0.0800 ; BIP(10:n,1)= 0.0800 ; 

BIP(2,3:10)= 0.120 ;BIP(3:10,2 )= 0.120 ; 
BIP(2,11:n)= 0.10 ; BIP(11:n , 2)= 0.10;

ncomp = length(comp);

molecular_weights = [28.01,44.01,16.04,30.07,44.10,58.12,58.12,72.15,72.15,86.18,96,107,121,140.09,167.57,204.75,243.6,275.27,317.02,370.39,456.83,640.76]'; 
Cp_coeffs = [
    % Component    C1        C2        C3          C4
    31.15,     -0.014,    2.68e-5,   -1.17e-8;    % N2 (from PVTsim)
    19.79,      0.073,   -5.60e-5,    1.72e-8;    % CO2
    19.25,      0.052,    1.20e-5,   -1.13e-8;    % C1
     5.41,      0.178,   -6.94e-5,    8.71e-9;    % C2
    -4.22,      0.306,   -1.59e-4,    3.21e-8;    % C3
    -1.39,      0.385,   -1.83e-4,    2.90e-8;    % iC4
     9.49,      0.331,   -1.11e-4,   -2.82e-9;    % nC4
    -9.52,      0.507,   -2.73e-4,    5.72e-8;    % iC5
    -3.63,      0.487,   -2.58e-4,    5.30e-8;    % nC5
    -4.41,      0.582,   -3.12e-4,    6.49e-8;    % C6
    -5.15,      0.676,   -3.65e-4,    7.66e-8;    % nC7
    -6.10,      0.771,   -4.20e-4,    8.85e-8;    % nC8
     3.14,      0.677,   -1.93e-4,   -2.98e-8;    % nC9
    25.20,      0.830,   -3.23e-4,    4.06e-8;    % C10-C11
    30.14,      0.993,   -3.87e-4,    4.86e-8;    % C12-C13
    36.83,      1.213,   -4.72e-4,    5.93e-8;    % C14-C16
    43.82,      1.443,   -5.62e-4,    7.06e-8;    % C17-C18
    49.52,      1.630,   -6.35e-4,    7.98e-8;    % C19-C21
    57.03,      1.878,   -7.31e-4,    9.19e-8;    % C22-C24
    66.63,      2.194,   -8.55e-4,    1.07e-7;    % C25-C29
    82.18,      2.706,   -1.05e-3,    1.32e-7;    % C30-C37
   115.26,      3.795,   -1.48e-3,    1.86e-7     % C38-C80
];

% Reference enthalpies
H_ig_ref_per_mass = [-20; 20; 0; 7.5; 15; 17; 17; 25; 25; 33; 2; 2; 2; 2; 93; 93; 93; 31; 31; 31; 8; 8];


temp = 93 + 273.15 ;
press = 284e5 ;
tol = 10^(-30);
maxiter = 1000;
% psat = 272 ; 
 R = 8.3144598;
% h_ref = 175;
% h = [167, 175, 204, 228, 327];
temp_gradient = - 0.025 ;
%pressb_ini = 272e5;*[comp_h, press_h, press_b_h] = main(h(3), comp, press, temp, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini);

 % [comp_h, press_h, press_b_h] = main_nonisothermal(h(5), comp, press, temp, temp_gradient, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini)



 %%
temp_ref = 93 + 273.15;
press_ref = 284e5;
h_ref = 175;
h = 125:5:350;

pressb_ini = 280e5;
p_iso = zeros(size(h));
pb_iso = zeros(size(h));
p_noniso = zeros(size(h));
pb_noniso = zeros(size(h));

fprintf('Calculating both models...\n');

for i = 1:length(h)
    % Isothermal
    [~, p_iso(i), pb_iso(i)] = main(h(i), comp, press_ref, temp_ref, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini);
    
    % Non-isothermal  
    [~, p_noniso(i), pb_noniso(i)] = main_nonisothermal(h(i), comp, press_ref, temp_ref, temp_gradient, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini, molecular_weights, Cp_coeffs, H_ig_ref_per_mass);
    
    if ~isnan(pb_iso(i))
        pressb_ini = pb_iso(i);
    end
end

% Plot
figure;
plot(p_iso/1e5, h, 'b-', 'LineWidth', 2);
hold on;
plot(pb_iso/1e5, h, 'b--', 'LineWidth', 2);
% plot(p_noniso/1e5, h, 'r-', 'LineWidth', 2);
plot(pb_noniso/1e5, h, 'r--', 'LineWidth', 2);

xlabel('Pressure (bar)');
ylabel('Depth (m)');
legend('P_{res} isothermal', 'P_{bub} isothermal', 'P_{res} non-isothermal', 'P_{bub} non-isothermal');
set(gca, 'YDir', 'reverse');
grid on;
title('Pressure vs Depth: Isothermal vs Non-isothermal');
% Experimental field data from Table 1
field_depth = [ 175, 204, 228, 327];
field_pressure = [284, 286, 287, 293] * 1e5; % Convert to Pa
field_psat = [272, 267, 265, 242] * 1e5; % Convert to Pa

% Add to your existing plot
plot(field_pressure/1e5, field_depth, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'black', 'DisplayName', 'Field P_{res}');
plot(field_psat/1e5, field_depth, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'black', 'DisplayName', 'Field P_{sat}');




%%
[K, comp_vap, comp_liq, phasefrac] = vaporliquideq(press, temp, comp, Pc, Tc, acentric, BIP, tol, maxiter);
pressb_ini = pressbubest_multicomp(comp, temp, Pc, Tc, acentric);
[pressb_h, ~] = pressbub_multicomp_newton(comp,pressb_ini, temp, Pc, Tc, acentric, BIP, tol, maxiter);

% R = 8.3144598;   
% c = [-4.23, -1.91 ,-5.20 ,-5.79 ,-6.35 ,-7.18 ,-6.49 ,-6.20 ,-5.12 ,1.42 ,8.45 ,9.20 ,10.32 ,11.12 ,10.48 ,6.57 ,-1.30 ,-9.89 ,-23.43 ,-43.30 ,-79.96 ,-151.44]'* 1e-6;
% c_mix = sum(comp .* c)* press_h / (R * temp );






