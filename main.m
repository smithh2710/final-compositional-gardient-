%% isothermal reservoirs 

function [comp_h, press_h, pressb_h] = main(h, comp_ref, press_ref, temp, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini)

    if nargin < 11
        try
            pressb_ini = pressbubest_multicomp(comp_ref, temp, Pc, Tc, acentric);
        catch
            pressb_ini = 250e5; % Default fallback
        end
    end

    R = 8.3144598;
    g = 9.80665;
    n = length(comp_ref);
    tol = 1e-12;
    maxiter = 1000;

    % Calculate target fugacities using gravity segregation
    [fugcoef_ref, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp, Pc, Tc, acentric, BIP);
    f_ref = fugcoef_ref .* comp_ref * press_ref;
    f_h = f_ref .* exp((M .* g * (h - h_ref)) ./ (R * temp));

    initial_guess = [comp_ref; press_ref];

    % Solve system
    fun = @(x) residual_fugacity(x(1:n), x(end), f_h, temp, Pc, Tc, acentric, BIP);
    options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

    try
        solution = fsolve(fun, initial_guess, options);
        comp_h = solution(1:n);
        press_h = solution(end);
        
        % Ensure positive compositions and normalize
        comp_h = max(comp_h, 1e-10);
        comp_h = comp_h / sum(comp_h);
    catch
        warning('fsolve failed at depth %.1f m', h);
        comp_h = comp_ref;
        press_h = press_ref;
    end

    % Calculate bubble point pressure
    try
        [pressb_h, ~] = pressbub_multicomp_newton(comp_h, pressb_ini, temp, Pc, Tc, acentric, BIP, tol, maxiter);
        
        % Validate bubble point result
        if ~isreal(pressb_h) || pressb_h <= 0 || ~isfinite(pressb_h)
            pressb_h = NaN;
        end
    catch
        pressb_h = NaN;
    end
end

function F = residual_fugacity(comp_h, press_h, f_h, temp, pressc, tempc, acentric, BIP)
    n = length(f_h);
    [fugcoef_h, ~] = fugacitycoef_multicomp(comp_h, press_h, temp, pressc, tempc, acentric, BIP);
    F = zeros(n+1, 1);
    for i = 1:n
        F(i) = comp_h(i) * press_h * fugcoef_h(i) - f_h(i);
    end
    F(n+1) = sum(comp_h) - 1;
end