comp = [1.227, 0.824, 47.221, 7.560, 4.414, 0.844, 2.141, 0.925, 1.126, 1.588, 2.765, 2.528, 2.310, 5.806, 4.434, 3.386, 2.586, 2.523, 1.761, 1.697, 1.317, 1.015]';

comp = comp ./ sum(comp);
M_wgt = [28.014, 44.010, 16.043, 30.070, 44.097, 58.124, 58.124, 72.151, 72.151, 86.178, 96.000, 107.000, 121.000, 146.526, 189.406, 235.798, 275.597, 323.260, 379.146, 447.348, 547.702, 754.228]';
molecular_weights = [28.014, 44.010, 16.043, 30.070, 44.097, 58.124, 58.124, 72.151, 72.151, 86.178, 96.000, 107.000, 121.000, 146.526, 189.406, 235.798, 275.597, 323.260, 379.146, 447.348, 547.702, 754.228]' ;
Tc = [126.200, 304.200, 190.600, 305.400, 369.800, 408.100, 425.200, 460.400, 469.600, 507.400, 537.031, 557.574, 581.291, 620.907, 677.320, 730.657, 772.695, 820.224, 872.393, 933.128, 1017.950, 1190.769]' ;

Pc = [3394.39, 7376.46, 4600.15, 4883.87, 4245.52, 3647.70, 3799.69, 3384.26, 3374.12, 2968.82, 2960.14, 2729.81, 2488.77, 2174.44, 1848.95, 1642.73, 1530.54, 1436.19, 1359.69, 1295.57, 1234.02, 1165.75]' * 1000;

acentric = [0.0400, 0.2250, 0.0080, 0.0980, 0.1520, 0.1760, 0.1930, 0.2270, 0.2510, 0.2960, 0.3375, 0.3742, 0.4204, 0.5020, 0.6263, 0.7484, 0.8421, 0.9406, 1.0341, 1.1144, 1.1606, 0.9203]';
c = [-4.23e-6, -1.64e-6, -5.20e-6, -5.79e-6, -6.35e-6, -7.18e-6, -6.49e-6, -6.20e-6, -5.12e-6, 1.39e-6, 7.40e-6, 0.000010, 0.000013, 0.000019, 0.000024, 0.000023, 0.000020, 0.000013, 9.23e-7, -0.000017, -0.000049, -0.000117];

n = length(comp);
BIP = zeros(n,n);


BIP(1,2) = -0.017; BIP(2,1) = -0.017;  % N2-CO2
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

temp_ref = 90.0 + 273.15;  
press_ref = 243.0e5;       
h_ref = 6.5;               
temp_gradient = -0.025;    
% depths = [0, 6.5, 11.5, 57.5]; 
pressb_ini = 243e5;
 h = 0:0.5:60;
% [comp_h, press_h, press_b_h] = main(h(2), comp, press, temp, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini);

% [C1, C2, C3, C4] where Cp = C1 + C2*T + C3*T^2 + C4*T^3
Cp_coeffs = [
    31.15,    -0.014,    2.68e-5,   -1.17e-8;    % N2
    19.79,     0.073,   -5.60e-5,    1.72e-8;    % CO2
    19.25,     0.052,    1.20e-5,   -1.13e-8;    % C1
    5.41,      0.178,   -6.94e-5,    8.71e-9;    % C2
    -4.22,     0.306,   -1.59e-4,    3.21e-8;    % C3
    -1.39,     0.385,   -1.83e-4,    2.90e-8;    % iC4
    9.49,      0.331,   -1.11e-4,   -2.82e-9;    % nC4
    -9.52,     0.507,   -2.73e-4,    5.72e-8;    % iC5
    -3.63,     0.487,   -2.58e-4,    5.30e-8;    % nC5
    -4.41,     0.582,   -3.12e-4,    6.49e-8;    % C6
    -17.97,    0.597,   -1.99e-4,    0.000;      % C7
    -6.89,     0.625,   -2.32e-4,    0.000;      % C8
    -2.30,     0.693,   -2.67e-4,    0.000;      % C9
    3.03,      0.833,   -3.31e-4,    0.000;      % C10-C12
    7.19,      1.077,   -4.29e-4,    0.000;      % C13-C15
    10.62,     1.345,   -5.34e-4,    0.000;      % C16-C18
    12.99,     1.573,   -6.24e-4,    0.000;      % C19-C21
    16.00,     1.852,   -7.33e-4,    0.000;      % C22-C25
    19.74,     2.177,   -8.59e-4,    0.000;      % C26-C29
    23.78,     2.574,   -1.01e-3,    0.000;      % C30-C35
    29.56,     3.158,   -1.24e-3,    0.000;      % C36-C44
    42.95,     4.441,   -1.75e-3,    0.000;      % C45-C80
];


% For Reservoir 2: H^ig/(M*R) in K/g
H_ig_ref_per_mass = [
    -20;   % N2
    20;    % CO2  
    0;     % C1
    7.5;   % C2
    15;    % C3
    17;    % iC4
    17;    % nC4
    25;    % iC5
    25;    % nC5
    33;    % C6
    2;     % C7 (M < 150)
    2;     % C8 (M < 150)
    2;     % C9 (M < 150)
    93;    % C10-C12 (150 ≤ M ≤ 250)
    93;    % C13-C15 (150 ≤ M ≤ 250)
    93;    % C16-C18 (150 ≤ M ≤ 250)
    31;    % C19-C21 (250 ≤ M ≤ 400)
    31;    % C22-C25 (250 ≤ M ≤ 400)
    31;    % C26-C29 (250 ≤ M ≤ 400)
    8;     % C30-C35 (M ≥ 400)
    8;     % C36-C44 (M ≥ 400)
    8;     % C45-C80 (M ≥ 400)
];



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

field_depth = [ 0, 6.5, 11.5, 57.5];
field_pressure = [243, 243, 243.3, 246.4] * 1e5; % Convert to Pa
field_psat = [243, 243,242.12,235.33] * 1e5; % Convert to Pa

% Add to your existing plot
plot(field_pressure/1e5, field_depth, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'black', 'DisplayName', 'Field P_{res}');
plot(field_psat/1e5, field_depth, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'black', 'DisplayName', 'Field P_{sat}');

