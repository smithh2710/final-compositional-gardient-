%% CALCULATE DIMENSIONLESS ATTRACTION AND COVOLUME, A & B
% The Definition of Variables.
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
% ncomp   : the number of components

function [A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric)

ncomp = size(pressc,1);
m = zeros(ncomp,1);
alpha = zeros(ncomp,1);
A = zeros(ncomp,1);
B = zeros(ncomp,1);
% L = zeros(ncomp,1);
% M = zeros(ncomp,1);
omegaa = 0.45724;
omegab = 0.07780;

for i = 1:ncomp
    
    % Calculate m.
    if acentric(i) > 0.491
        m(i) = 0.37964 + 1.48503*acentric(i) - 0.164423*acentric(i)^2 + 0.016666*acentric(i)^3;
    else
        m(i) = 0.37464 + 1.54226*acentric(i) - 0.26992*acentric(i)^2;
    end
    
    % Calculate reduced pressure and temperature.
    pressr = press/pressc(i);
    tempr = temp/tempc(i);
    
    % Calculate alpha.
     alpha(i) = ( 1 + m(i)*(1 - sqrt(tempr)) )^2;
% 
% alpha(i) = exp((2 + 0.836 * tempr) * (1 - tempr^(0.134 + 0.508 * acentric(i) - 0.0467 * acentric(i)^2)));

%  L(i) = 0.0728 + 0.66693* acentric(i) + 0.0925* acentric(i)^2 ; %%% two alpha
% M(i) = 0.8788 - 0.2258 * acentric(i) + 0.1695*acentric(i)^2 ; 
%  alpha(i) = (tempr)^(2*( M(i) - 1 ) ) * exp(L(i)*(1 - tempr^(2*M(i)))); 



    A(i) = omegaa*alpha(i)*pressr/tempr^2;
    B(i) = omegab*pressr/tempr;
    
end

end
