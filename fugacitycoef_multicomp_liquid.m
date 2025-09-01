%% CALCULATE THE FUGACITY COEFFICIENT AND Z-FACTOR OF MULTI-COMPONENT SYSTEMS
% -------------------------------------------------------------------------
% The Definition of Variables.
% comp    : composition
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
% BIP     : binary interaction parameter
% fugcoef : fugacity coefficient
% zfactor : compressibility factor
% -------------------------------------------------------------------------
% In this function, the minimum z-factor is automatically chosen.
function [fugcoef, zfactor] = fugacitycoef_multicomp_liquid(comp, press, temp, pressc, tempc, acentric, BIP)

[A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric);

[Amix, Bmix, Amix2] = calcabmix(comp, A, B, BIP);

zfactor = calczfactor(Amix, Bmix);

if (size(zfactor,1) > 1)
    zfactor = min(zfactor);
end

fugcoef = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2);

% c = [-4.23, -1.64, -5.20, -5.79, -6.35, -7.18, -6.49, -6.20, -5.12, 1.39, 11.65, 14.54, 17.99, 26.20, 30.11, 24.59, 13.83, -3.45, -25.73, -57.40, -103.65, -178.08]' * 1e-6;
%  R = 8.3144598;
%     fugcoef = zeros(size(comp));
% 
%     for i = 1:length(comp)
%         fugcoef(i) = fugcoef_pr(i) / exp((c(i) * press) / (R * temp));
%     end



end
