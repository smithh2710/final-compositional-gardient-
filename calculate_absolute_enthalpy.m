function [H_abs_components, H_abs_mixture, H_ig_specific, H_res_specific, H_abs_specific] = calculate_absolute_enthalpy(temp, press, comp, Pc, Tc, acentric, BIP, molecular_weights, Cp_coeffs, H_ig_ref_per_mass)
% Calculate absolute enthalpy for each component and mixture
% Inputs:
% temp - Temperature in Kelvin
% press - Pressure in Pa (or bar, ensure consistency with fugacity function)
% comp - Mole fraction vector for all components
% Pc, Tc, acentric - Critical properties vectors
% BIP - Binary interaction parameters matrix
% molecular_weights - Molecular weights vector (g/mol)
% Cp_coeffs - Heat capacity coefficients matrix (n_comp x 4) [C1, C2, C3, C4] in J/mol/K
% H_ig_ref_per_mass - Reference ideal gas specific enthalpies at 273.15 K (H^ig/(M*R) in K/g)
% Outputs:
% H_abs_components - Absolute enthalpy for each component (J/mol)
% H_abs_mixture - Mixture absolute enthalpy (J/mol)
% H_ig_specific - Ideal gas specific enthalpy H_ig/(M*R) (K/g)
% H_res_specific - Residual specific enthalpy H_res/(M*R) (K/g)
% H_abs_specific - Absolute specific enthalpy H_abs/(M*R) (K/g)

% Validate inputs
if length(comp) ~= length(Pc) || length(comp) ~= length(Tc) || length(comp) ~= length(acentric)
    error('Component vectors must have same length');
end

if length(comp) ~= length(molecular_weights) || length(comp) ~= length(H_ig_ref_per_mass)
    error('Molecular weights and reference enthalpies must match component count');
end

if size(Cp_coeffs, 1) ~= length(comp) || size(Cp_coeffs, 2) ~= 4
    error('Cp_coeffs must be n_components x 4 matrix');
end

R = 8.3144598; % J/mol/K

% Calculate ideal gas enthalpy for each component
H_ig_components = calculate_ideal_gas_enthalpy_components(temp, molecular_weights, Cp_coeffs, H_ig_ref_per_mass);

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

function H_ig_components = calculate_ideal_gas_enthalpy_components(T, molecular_weights, Cp_coeffs, H_ig_ref_per_mass)
% Calculate ideal gas enthalpy for each component at temperature T
% Output: H_ig_components - Ideal gas enthalpy in J/mol for each component

T_ref = 273.15; % Reference temperature (K)
R = 8.3144598; % Gas constant (J/mol/K)

n_components = length(molecular_weights);
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

for i = 1:n_comp
    % Calculate partial molar residual enthalpy for component i
    % H̄i_res = -R*T^2 * (∂ ln φi / ∂T)_P,x
    dln_phi_i_dT = dln_phi_dT(i);
    H_res_components(i) = -R * temp^2 * dln_phi_i_dT;
end
end