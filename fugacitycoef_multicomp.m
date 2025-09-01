
function [fugcoef, zfactor] = fugacitycoef_multicomp(comp, press, temp, pressc, tempc, acentric, BIP)

[A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric);

[Amix, Bmix, Amix2] = calcabmix(comp, A, B, BIP);

zfactor = calczfactor(Amix, Bmix);

if (size(zfactor,1) > 1)
    zfactor = choosezfactor(zfactor, comp, A, B, Amix, Bmix, Amix2);
end

fugcoef = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2);

% c = [-4.23, -1.64, -5.20, -5.79, -6.35, -7.18, -6.49, -6.20, -5.12, 1.39, 11.65, 14.54, 17.99, 26.20, 30.11, 24.59, 13.83, -3.45, -25.73, -57.40, -103.65, -178.08]' * 1e-6;
%  R = 8.3144598 ;
%     fugcoef = zeros(size(comp));
% 
%     for i = 1:length(comp)
%         fugcoef(i) = fugcoef_pr(i) / exp(( (c(i) ) * press) / (R * temp));
%     end


end

%% SEARTCH AND RETURN AN APPROPRIATE Z-FACTOR
% Calculate dimensionless excess gibbs free energy, and return the z
% factor which minimizes the gibbs free energy.
function minzfactor = choosezfactor(zfactor, comp, A, B, Amix, Bmix, Amix2)

gibbsenergy = [];

for i = 1:size(zfactor,1)
    fugcoef = calcfugcoef_multicomp(zfactor(i), A, B, Amix, Bmix, Amix2);
    g = calcgibbsenergy(comp, fugcoef);
    gibbsenergy = cat(1,gibbsenergy,g);
end

[~, index] = sort(gibbsenergy);
minzfactor = zfactor(index(1));

end

function g = calcgibbsenergy(comp, fugcoef)

ncomp = size(comp,1);
g = 0;

for i = 1:ncomp
    if comp(i) ~= 0
        g = g + comp(i)*log(comp(i)*fugcoef(i));
    end
end

end