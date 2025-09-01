
% Demonstrate the absolute enthalpy calculation

% Example conditions (Reservoir 1 from Pedersen paper)
temp = 93 + 273.15;  % 366.15 K
press = 284e5;       % 284 bar converted to Pa (adjust units as needed for your fugacity function)

% Example composition (22 components) - replace with actual values
comp = [0.42, 0.69, 50.04, 7.85, 6.77, 1.04, 3.2, 1.16, 1.55, 1.88, 3.50 , 3.75, 2.28, 3.26, 2.59, 2.93 ,1.46 ,1.65 ,1.17 ,1.24 ,0.96 ,0.63]';

% Ensure composition sums to 1
comp = comp / sum(comp);

Pc = [33.94 ,73.76 ,46.00 ,48.84 ,42.46 ,36.48 ,38 ,33.84 ,33.74 ,29.69 ,29.60 ,27.86 ,25.83 ,23.76 ,21.58 ,19.59 ,18.24 ,17.54 ,16.84 ,16.23 ,15.62 ,14.72 ]' * 1e5 ;

Tc = [-147 ,31.1 ,-82.5 ,32.2 ,96.6 ,134.9 ,152.1 ,187.2 ,196.4 234.2 ,269.9 ,290.9 ,315.3 ,345.9 ,385 ,432.7 ,476.8 ,511.1 ,553.3 ,604.9 ,683.3 ,847.8]' + 273.15; % Convert to K

acentric = [0.0400 ,0.2250 ,0.0080 ,0.0980 ,0.1520 ,0.1760 ,0.1930 ,0.2270 ,0.2510 ,0.2960 ,0.338 ,0.374 ,0.420 ,0.483 ,0.570 ,0.685 ,0.795 ,0.881 ,0.984 ,1.098 ,1.229 ,1.159]';

molecular_weights = [28.0 ,44.0 ,16.0 ,30.1 ,44.1 ,58.1 ,58.1 ,72.2 ,72.2 ,86.2 ,96.0 ,107 ,121 ,140.1 ,167.6, 204.7 ,243.6 ,275.3 ,317.0 ,370.4 ,456.8 ,640.9]';

n = length(comp);
BIP = zeros(n, n); % Simplified - use actual BIP values

BIP(1,2) = -0.017    ; BIP(2,1) = -0.017 ; 
BIP(1,3) = 0.031 ; BIP(3,1)= 0.031;
BIP(1,4) = 0.052 ; BIP(4,1) = 0.052 ; 
BIP(1,5) = 0.085 ; BIP(5,1)= 0.085 ; 
BIP(1,6) = 0.103 ; BIP(6,1)=0.103; 
BIP(1,7) = 0.0800 ; BIP(7,1) = 0.0800;
BIP(1,8)= 0.092 ; BIP(8,1) = 0.092 ; 
BIP(1,9)= 0.1000 ; BIP(9,1)=0.1000;
BIP(1,10:n)= 0.0800 ; BIP(10:n,1)= 0.0800 ; 

BIP(2,3:10)= 0.120 ;BIP(3:10,2 )= 0.120 ; 
% BIP(2,11:n)= 0.10 ; BIP(11:n , 2)= 0.10;


% Calculate absolute enthalpies
[H_abs_components, H_abs_mixture, H_ig_specific, H_res_specific, H_abs_specific] = ...
    calculate_absolute_enthalpy(temp, press, comp, Pc, Tc, acentric, BIP);

% Display results
fprintf('Absolute Enthalpy Calculation Results:\n');
fprintf('Temperature: %.1f K\n', temp);
fprintf('Pressure: %.1f bar\n', press/1e5);
fprintf('Mixture Absolute Enthalpy: %.1f J/mol\n', H_abs_mixture);

% Display detailed component results
fprintf('\nComponent Results:\n');
fprintf('%-8s %12s %12s %12s %12s %12s\n', 'Comp', 'H_abs(J/mol)', 'H_ig/(M*R)', 'H_res/(M*R)', 'H_abs/(M*R)', 'Units');
fprintf('%-8s %12s %12s %12s %12s %12s\n', '', '', '(K/g)', '(K/g)', '(K/g)', '');
fprintf(repmat('-', 1, 80));
fprintf('\n');

component_names = {'N2', 'CO2', 'C1', 'C2', 'C3', 'iC4', 'nC4', 'iC5', 'nC5', 'C6', ...
                   'nC7', 'nC8', 'nC9', 'C10-C11', 'C12-C13', 'C14-C16', 'C17-C18', ...
                   'C19-C21', 'C22-C24', 'C25-C29', 'C30-C37', 'C38-C80'};

for i = 1:length(H_abs_components)
    fprintf('%-8s %12.1f %12.1f %12.1f %12.1f\n', component_names{i}, ...
            H_abs_components(i), H_ig_specific(i), H_res_specific(i), H_abs_specific(i));
end

% Summary for Haase model usage
fprintf('\nFor Haase Model Compositional Grading:\n');
fprintf('Use H_abs/(M*R) values (K/g units) - these correspond to Table 8 in Pedersen paper\n');


plot( molecular_weights , H_ig_specific  ) ; hold on ; 
plot(molecular_weights , H_res_specific ) ; hold on ; 
plot(molecular_weights , H_abs_specific  ); hold off 
ylim([-50, 950]);
yticks(-50:100:950);

function [H_abs_components, H_abs_mixture, H_ig_specific, H_res_specific, H_abs_specific] = calculate_absolute_enthalpy(temp, press, comp, Pc, Tc, acentric, BIP)
% Calculate absolute enthalpy for each component and mixture
% Inputs:
%   temp - Temperature in Kelvin
%   press - Pressure in Pa (or bar, ensure consistency with fugacity function)
%   comp - Mole fraction vector for all components
%   Pc, Tc, acentric - Critical properties vectors
%   BIP - Binary interaction parameters matrix
% Outputs:
%   H_abs_components - Absolute enthalpy for each component (J/mol)
%   H_abs_mixture - Mixture absolute enthalpy (J/mol)
%   H_ig_specific - Ideal gas specific enthalpy H_ig/(M*R) (K/g)
%   H_res_specific - Residual specific enthalpy H_res/(M*R) (K/g)
%   H_abs_specific - Absolute specific enthalpy H_abs/(M*R) (K/g)

% Validate inputs
if length(comp) ~= length(Pc) || length(comp) ~= length(Tc) || length(comp) ~= length(acentric)
    error('Component vectors must have same length');
end

function molecular_weights = get_molecular_weights()
% Molecular weights (g/mol) for all components
molecular_weights = [28.0; 44.0; 16.0; 30.1; 44.1; 58.1; 58.1; 72.2; 72.2; 86.2; ...
                     96.0; 107.0; 121.0; 140.1; 167.6; 204.7; 243.6; 275.3; 317.0; 370.4; 456.8; 640.9];
end

% Get molecular weights
molecular_weights = get_molecular_weights();
R = 8.3144598; % J/mol/K

% Calculate ideal gas enthalpy for each component
H_ig_components = calculate_ideal_gas_enthalpy_components(temp);

% Calculate residual enthalpy for each component
[H_res_mixture, H_res_components] = calculate_residual_enthalpy_corrected(temp, press, comp, Pc, Tc, acentric, BIP);

% Calculate absolute enthalpy: H_abs = H_ig + H_res
H_abs_components = H_ig_components + H_res_components;

% Calculate mixture absolute enthalpy
H_abs_mixture = sum(comp .* H_abs_components);

% Calculate specific enthalpies: H/(M*R) in units K/g
H_ig_specific = H_ig_components ./ (molecular_weights * R);
H_res_specific = H_res_components ./ (molecular_weights * R);
H_abs_specific = H_abs_components ./ (molecular_weights * R);

end

function H_ig_components = calculate_ideal_gas_enthalpy_components(T)
% Calculate ideal gas enthalpy for each component at temperature T
% Output: H_ig_components - Ideal gas enthalpy in J/mol for each component

T_ref = 273.15; % Reference temperature (K)
R = 8.3144598;  % Gas constant (J/mol/K)

% Reference ideal gas specific enthalpies at 273.15 K (H^ig/(M*R) in K/g)
defined_components = [-20; 20; 0; 7.5; 15; 17; 17; 25; 25; 33]; % N2, CO2, C1-C6

% ALL C7+ should use Table 8 step function values (these are correct)  
pseudocomponent_values = [
    2;  % C7  (MW=96.0, ≤150)
    2;  % C8  (MW=107.0, ≤150) 
    2;  % C9  (MW=121.0, ≤150)
    2;  % C10-C11 (MW=140.1, ≤150)
    93; % C12-C13 (MW=167.6, 150-250) - HIGH aromatic value
    93; % C14-C16 (MW=204.7, 150-250) - HIGH aromatic value  
    93; % C17-C18 (MW=243.6, 150-250) - HIGH aromatic value
    93; % C19-C21 (MW=275.3, 250-400)
    31; % C22-C24 (MW=317.0, 250-400)
    31; % C25-C29 (MW=370.4, 250-400)
    8;  % C30-C37 (MW=456.8, >400)
    8   % C38-C80 (MW=640.9, >400)
];

H_ig_ref_per_mass = [defined_components; pseudocomponent_values];

% Molecular weights (g/mol)
molecular_weights = [28.0 ,44.0 ,16.0 ,30.1 ,44.1 ,58.1 ,58.1 ,72.2 ,72.2 ,86.2 ,96.0 ,107 ,121 ,140.1 ,167.6, 204.7 ,243.6 ,275.3 ,317.0 ,370.4 ,456.8 ,640.9]';

% Get heat capacity coefficients (CORRECTED)
Cp_coeffs = get_cp_coefficients_corrected();

n_components = length(H_ig_ref_per_mass);
H_ig_components = zeros(n_components, 1);

for i = 1:n_components
    % Get coefficients for this component
    C1 = Cp_coeffs(i, 1);
    C2 = Cp_coeffs(i, 2);
    C3 = Cp_coeffs(i, 3);
    C4 = Cp_coeffs(i, 4);
    
    % Calculate enthalpy change from integration: ∫Cp dT from T_ref to T
    dH_ig = C1 * (T - T_ref) + ...
            C2/2 * (T^2 - T_ref^2) + ...
            C3/3 * (T^3 - T_ref^3) + ...
            C4/4 * (T^4 - T_ref^4);
    
    % Calculate absolute ideal gas enthalpy
    H_ig_ref_absolute = H_ig_ref_per_mass(i) * molecular_weights(i) * R; % J/mol
    H_ig_components(i) = H_ig_ref_absolute + dH_ig;
end

end

function [H_res_mixture, H_res_components] = calculate_residual_enthalpy_corrected(temp, press, comp, Pc, Tc, acentric, BIP)
% CORRECTED residual enthalpy calculation
% Uses thermodynamic relation: H_res = -R*T^2 * (d ln φ / dT)_P

R = 8.3144598;
dT = 0.1; % K - temperature perturbation

% Calculate mixture residual enthalpy
[phi_minus, ~] = fugacitycoef_multicomp(comp, press, temp-dT/2, Pc, Tc, acentric, BIP);
[phi_plus, ~] = fugacitycoef_multicomp(comp, press, temp+dT/2, Pc, Tc, acentric, BIP);
dln_phi_dT = (log(phi_plus) - log(phi_minus)) / dT;
H_res_mixture = -R * temp^2 * sum(comp .* dln_phi_dT);

% Calculate partial molar residual enthalpies
n_comp = length(comp);
H_res_components = zeros(n_comp, 1);

% CORRECTED METHOD: Use partial molar property relation
% H̄i_res = H_res + (∂H_res/∂ni)_T,P,nj - this requires more complex derivatives
% Simplified approach: Use component fugacity coefficients directly

for i = 1:n_comp
    % Calculate partial molar residual enthalpy for component i
    % H̄i_res = -R*T^2 * (∂ ln φi / ∂T)_P,x
    
    dln_phi_i_dT = dln_phi_dT(i); % This is already (∂ ln φi / ∂T)
    H_res_components(i) = -R * temp^2 * dln_phi_i_dT;
end

end

function Cp_coeffs = get_cp_coefficients_corrected()
% CORRECTED heat capacity coefficients from PVTsim
% Units: J/mol/K

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

end

