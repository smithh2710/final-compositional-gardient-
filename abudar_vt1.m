 function [v_translated,v_pr, C , zfactor , Z_roots_VT ,c] = abudar_vt1(comp, press, temp, Pc, Tc, acentric,BIP, Vc, c1)
% ABUDAR_VT1 Apply Abudour et al. (2012) volume translation to PR EOS
% Inputs:
%   comp     - Mole fractions (will be normalized)
%   press    - Pressure (Pa)
%   temp     - Temperature (K)
%   Pc       - Critical pressures (Pa)
%   Tc       - Critical temperatures (K)
%   acentric - Acentric factors
%   M    - Molecular weights (kg/mol)
%   BIP      - Binary interaction parameters
%   Vc       - Critical volumes (m³/mol)
%   c1       - Volume translation parameter
% Outputs:
%   v_translated        - Translated molar volume (m³/mol)
%   volume_change       - Dimensionless volume change
%   rho_pr             - Original PR density (kg/m³)
%   rho_translated     - Corrected density (kg/m³)
%   density_improvement - Percentage density improvement

% comp = [0.8635 , 0.1365]' ; 
% press = 423.5 * 0.00689475729 ; %MPa 
% temp =288.75 ; % k 
% Pc = [48.84 , 42.46]' * 10^(-1) ; % bar to MPa 
% Tc = [305.322 , 369.89]' ; % K 
% acentric = [0.0995 , 0.1521]' ; 
% M_gmol =[30.069 , 44.0956 ]' ; 
% M = M_gmol/1000 ; 
% 
% ncomp = length(comp); 
% BIP = zeros(ncomp);
% Vc = [145.8388 , 200 ]'; 
% c1 = [0.00993 , 0.00778]' ;

R = 8.3144598;  
zc_pr = 0.3074; 

% Input validation
if any(comp < 0) || abs(sum(comp) - 1) > 1e-6
    error('Invalid composition: must be positive and sum to 1');
end
if temp <= 0 || press <= 0
    error('Invalid temperature or pressure: must be positive');
end
if length(comp) ~= length(Pc) || length(comp) ~= length(Tc)
    error('Inconsistent array dimensions');
end

% Ensure column vectors and normalize composition
comp = comp(:); 
comp = comp/sum(comp);
Pc = Pc(:); 
Tc = Tc(:); 
acentric = acentric(:);

Vc = Vc(:);
c1 = c1(:);
ncomp = numel(comp);

Tr = temp ./ Tc;
alpha_gasem = exp((2 + 0.836.*Tr) .* (1 - Tr.^(0.134 + 0.508.*acentric - 0.0467.*acentric.^2)));
a_pr_i = 0.45724  .* (R .* Tc).^2 .* alpha_gasem ./ Pc;
b_pr_i = 0.07780 .* R .* Tc ./ Pc;

a = zeros(ncomp, ncomp);
b = zeros(ncomp, ncomp);
% am = zeros(ncomp, ncomp);
% bm = zeros(ncomp, ncomp); 

for i = 1:ncomp
    for j = 1:ncomp
        a(i,j) = sqrt(a_pr_i(i) * a_pr_i(j)) * (1 - BIP(i,j)); 
        b(i,j) = 0.5 * (b_pr_i(i) + b_pr_i(j)); % Dij = 0
    end
end


% for i = 1:ncomp
%     for j = 1:ncomp
%    am(i,j)= comp(i)*comp(j) * a(i,j);
%    bm(i,j)= comp(i)*comp(j) * b(i,j);
%     end
% end
% a_mix = sum(am);
% b_mix = sum(bm);

a_mix = comp.' * a * comp;
b_mix = comp.' * b * comp;


[~, zfactor] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP);
v_pr = zfactor * R * temp / press;


if v_pr <= b_mix
    warning('PR volume less than covolume - possible phase misidentification');
end

% Calculate mixture critical properties using Chueh-Prausnitz correlations
theta = comp .* Vc.^(2/3); 
theta = theta / sum(theta);  % Surface fractions
Tc_m = sum(theta .* Tc);
Vc_m = sum(theta .* Vc);
omega_m = sum(comp .* acentric);
Pc_m = (0.2905 - 0.085*omega_m) * R * Tc_m / Vc_m;

% Calculate distance function d^PR (Equation 2.96 from Abudour et al.)
% ∂P/∂v at constant T for PR EOS
DPDV = R * temp / (v_pr - b_mix)^2 - ...
       a_mix * (2 * v_pr + 2 * b_mix) / ((v_pr^2 + 2 * b_mix * v_pr - b_mix^2)^2);
d_pr = (v_pr^2 / (R * Tc_m)) * DPDV;


C1_m = sum(comp .* c1);

cm = (R * Tc_m / Pc_m) * (C1_m - (0.004 + C1_m) * exp(-2*d_pr));

Vcm_pr = (R * Tc_m / Pc_m) * zc_pr;
delta_pr = Vcm_pr - Vc_m;

v_translated = v_pr + cm - delta_pr * (0.35 / (0.35 + d_pr));

if v_translated <= 0
    warning('Non-physical volume after translation, reverting to PR volume');
    v_translated = v_pr;
end


% M_mix = sum(comp .* M); % kg/mol
% rho_pr = M_mix / v_pr;                % Original PR density (kg/m³)
% rho_translated = M_mix / v_translated; % Corrected density (kg/m³)
 

c = (v_translated-v_pr );
C = c * press / (R * temp); % Dimensionless volume change
% zfactor_vt = press * v_translated / ( R * temp ) ; 


bvt = b_mix - c ;  
    
% Calculate dimensionless parameters
A = a_mix * press / (R * temp^2);
B = b_mix * press / (R * temp);    % Original B
BVT = bvt * press / (R * temp);    % BVT = B for standard translation
    
% a1 =1 ; 
% a2 = BVT  + 4*C - 1;
% a3 = 2*C^2 - 4*C*BVT + 2*BVT - 4*C + A - 1;
% a4 = -A*BVT - (2*C^2 - BVT^2)*(1 + BVT);


a1 =1 ; 
a2 = B  + 4*C - 1;
a3 = 2*C^2 - 4*C*B + 2*B - 4*C + A - 1;
a4 = -A*B - (2*C^2 - B^2)*(1 + B);

% a1 = 1 ; 
% a2 = 3*C -1 ; 
% a3 = 2 * C^2 - 2*C*B - B - 3*C +A -1 ;
% a4 = -A*B - (C*B + 2*C^2)* (1+ B) ; 

coeffs_VT = [a1, a2, a3, a4];
Z_roots_VT = real ( roots(coeffs_VT) ) ; 


  end