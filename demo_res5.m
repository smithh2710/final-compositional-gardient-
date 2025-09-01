comp = [0.809, 0.095, 36.146, 10.620, 12.374, 2.637, 4.759, 2.048, 2.080, 2.647, 1.496, 1.410, 1.329, 5.589, 4.161, 2.550, 2.447, 1.822, 1.583, 1.403, 1.110, 0.884]' ;
comp = comp ./ sum(comp) ; 
M_wgt = [28.014, 44.010, 16.043, 30.070, 44.097, 58.124, 58.124, 72.151, 72.151, 86.178, 96.000, 107.000, 121.000, 159.752, 234.113, 296.195, 357.716, 428.351, 504.597, 600.681, 729.800, 936.663]';

Tc = [-146.950, 31.050, -82.550, 32.250, 96.650, 134.950, 152.050, 187.250, 196.450, 234.250, 267.521, 287.936, 311.539, 370.742, 459.450, 523.321, 582.327, 645.445, 710.459, 789.006, 890.300, 1048.191]' + 273.15  ;

Pc = [33.94, 73.76, 46.00, 48.84, 42.46, 36.48, 38.00, 33.84, 33.74, 29.69, 30.59, 28.17, 25.65, 21.11, 16.89, 15.17, 14.13, 13.35, 12.78, 12.29, 11.86, 11.46]' * 1e5;

acentric = [0.0400, 0.2250, 0.0080, 0.0980, 0.1520, 0.1760, 0.1930, 0.2270, 0.2510, 0.2960, 0.3379, 0.3746, 0.4208, 0.5458, 0.7466, 0.8875, 1.0015, 1.0963, 1.1530, 1.1536, 1.0309, 0.5485]';

c = [-4.23, -1.64, -5.20, -5.79, -6.35, -7.18, -6.49, -6.20, -5.12, 1.39, 11.65, 14.54, 17.99, 26.20, 30.11, 24.59, 13.83, -3.45, -25.73, -57.40, -103.65, -178.08]' * 1e-6;

n = length(comp);
BIP = zeros(n,n);
BIP(1,2) = -0.017; BIP(2,1) = -0.017;
BIP(1,3) = 0.0311 ; BIP(3,1)= 0.0311;
BIP(1,4) = 0.0515 ; BIP(4,1) = 0.0515 ; 
BIP(1,5) = 0.0852 ; BIP(5,1)= 0.0852 ; 
BIP(1,6) = 0.1033 ; BIP(6,1)=0.1033; 
BIP(1,7) = 0.0800 ; BIP(7,1) = 0.0800;
BIP(1,8)= 0.0922 ; BIP(8,1) = 0.0922 ; 
BIP(1,9)= 0.1000 ; BIP(9,1)=0.1000;
BIP(1,10:n)= 0.0800 ; BIP(10:n,1)= 0.0800 ;


BIP(2,3:10) = 0.120; BIP(3:10,2) = 0.120;
BIP(2,11:n) = 0.10; BIP(11:n,2) = 0.10;

M = M_wgt / 1000;

temp_ref = 79.8 + 273.15;  
press_ref = 194e5;       
h_ref = 8 ;               
temp_gradient = - 0.025;    
tol = 10^(-30);
maxiter = 1000 ; 

% Heat capacity coefficients for Reservoir 5 (J/mol/K)
% Cp = C1 + C2*T + C3*T^2 + C4*T^3
Cp_coeffs = [
    31.15,    -0.014,    2.68e-5,   -1.17e-8;     % N2
    19.79,     0.073,   -5.60e-5,    1.72e-8;     % CO2
    19.25,     0.052,    1.20e-5,   -1.13e-8;     % C1
     5.41,     0.178,   -6.94e-5,    8.71e-9;     % C2
    -4.22,     0.306,   -1.59e-4,    3.21e-8;     % C3
    -1.39,     0.385,   -1.83e-4,    2.90e-8;     % iC4
     9.49,     0.331,   -1.11e-4,   -2.82e-9;     % nC4
    -9.52,     0.507,   -2.73e-4,    5.72e-8;     % iC5
    -3.63,     0.487,   -2.58e-4,    5.30e-8;     % nC5
    -4.41,     0.582,   -3.12e-4,    6.49e-8;     % C6
   -26.75,     0.581,   -1.99e-4,    0.000;       % C7
   -16.84,     0.611,   -2.32e-4,    0.000;       % C8
   -12.71,     0.677,   -2.66e-4,    0.000;       % C9
    -5.18,     0.891,   -3.64e-4,    0.000;       % C10-C14
     1.07,     1.302,   -5.33e-4,    0.000;       % C15-C19
     4.43,     1.648,   -6.73e-4,    0.000;       % C20-C23
     7.97,     2.001,   -8.13e-4,    0.000;       % C24-C28
    11.71,     2.403,   -9.72e-4,    0.000;       % C29-C33
    15.37,     2.838,   -1.14e-3,    0.000;       % C34-C39
    20.25,     3.390,   -1.36e-3,    0.000;       % C40-C47
    27.47,     4.137,   -1.66e-3,    0.000;       % C48-C58
    39.98,     5.360,   -2.14e-3,    0.000;       % C59-C80
];

% Reference ideal gas specific enthalpies at 273.15 K
% H^ig/(M*R) in K/g - Based on Table 8 for Reservoir 5
H_ig_ref_per_mass = [
    -20;      % N2        (Table 6)
     20;      % CO2       (Table 6)
      0;      % C1        (Table 6)
    7.5;      % C2        (Table 6)
     15;      % C3        (Table 6)
     17;      % iC4       (Table 6)
     17;      % nC4       (Table 6)
     25;      % iC5       (Table 6)
     25;      % nC5       (Table 6)
     33;      % C6        (Table 6)
      2;      % C7        (M < 150)
      2;      % C8        (M < 150)
      2;      % C9        (M < 150)
    877;      % C10-C14   (150 ≤ M ≤ 250)
    877;      % C15-C19   (150 ≤ M ≤ 250)
    880;      % C20-C23   (250 ≤ M ≤ 400)
    880;      % C24-C28   (250 ≤ M ≤ 400)
    639;      % C29-C33   (M ≥ 400)
    639;      % C34-C39   (M ≥ 400)
    639;      % C40-C47   (M ≥ 400)
    639;      % C48-C58   (M ≥ 400)
    639;      % C59-C80   (M ≥ 400)
];

 h = [0, 8, 21,28]; 
 pressb_ini = pressbubest_multicomp(comp, temp_ref, Pc, Tc, acentric);
[~, p_noniso, pb_noniso] = main_nonisothermal(h(3), comp, press_ref, temp_ref, temp_gradient, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini, M, Cp_coeffs, H_ig_ref_per_mass)
  
%% 

% pressb_ini = pressbubest_multicomp(comp, temp_ref, Pc, Tc, acentric);
% [pressb, comp_vap] = pressbub_multicomp_newton(comp, pressb_ini, temp_ref, Pc, Tc, acentric, BIP, tol, maxiter)
% % [~, p_iso(i), pb_iso(i)] = main(h(2), comp, press_ref, temp_ref, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini);


h = 0:0.5:28;
p_iso = zeros(size(h));
pb_iso = zeros(size(h));
p_noniso = zeros(size(h));
pb_noniso = zeros(size(h));

fprintf('Calculating both models...\n');

for i = 1:length(h)
    % Isothermal
    [~, p_iso(i), pb_iso(i)] = main(h(i), comp, press_ref, temp_ref, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini);
    
    % Non-isothermal  
    [~, p_noniso(i), pb_noniso(i)] = main_nonisothermal(h(i), comp, press_ref, temp_ref, temp_gradient, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini, M, Cp_coeffs, H_ig_ref_per_mass);
    
    if ~isnan(pb_iso(i))
        pressb_ini = pb_iso(i);
    end
end

figure;
plot(p_iso/1e5, h, 'b-', 'LineWidth', 2);
hold on;
plot(pb_iso/1e5, h, 'b--', 'LineWidth', 2);
% plot(p_noniso/1e5, h, 'r-', 'LineWidth', 2);
plot(pb_noniso/1e5, h, 'r--', 'LineWidth', 2);

% xlabel('Pressure (bar)');
% ylabel('Depth (m)');
% legend('P_{res} isothermal', 'P_{bub} isothermal', 'P_{bub} non-isothermal');
set(gca, 'YDir', 'reverse');
grid on;
% title('Pressure vs Depth: Isothermal vs Non-isothermal');

% Set x-axis with 10 bar intervals
xlim([150 200]);
xticks(150:10:200);

field_depth = [0, 8, 21, 28];
field_pressure = [NaN, 194.0, 194.1, 195.1] * 1e5;
field_psat = [NaN, 185.7, 175.4, 152] * 1e5;

% Add experimental data points
plot(field_pressure/1e5, field_depth, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'black', 'DisplayName', 'Field P_{res}');
plot(field_psat/1e5, field_depth, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'black', 'DisplayName', 'Field P_{sat}');

% Update legend to include experimental data
% legend('P_{res} isothermal', 'P_{bub} isothermal', 'P_{bub} non-isothermal', 'Field P_{res}', 'Field P_{sat}');

