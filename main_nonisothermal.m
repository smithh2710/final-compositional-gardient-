function [comp_h, press_h, pressb_h] = main_nonisothermal(h, comp_ref, press_ref, temp_ref, temp_gradient, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini, molecular_weights, Cp_coeffs, H_ig_ref_per_mass)
    % Non-isothermal compositional grading using Haase model
    % Solves the full system directly with fsolve
    
    if nargin < 12
        try
            pressb_ini = pressbubest_multicomp(comp_ref, temp_ref, Pc, Tc, acentric);
        catch
            pressb_ini = 250e5; % Default fallback
        end
    end
    
    % Handle missing enthalpy parameters with default values (for backward compatibility)
    if nargin < 15
        error('Missing enthalpy parameters: molecular_weights, Cp_coeffs, H_ig_ref_per_mass required');
    end
    
    % R = 8.3144598;
    % g = 9.80665;
    n = length(comp_ref);
    tol = 1e-12;
    maxiter = 1000;
    
    % Calculate temperature at new depth
    temp_h = temp_ref + temp_gradient * (h - h_ref);
    
    % Calculate reference fugacities (same for both isothermal and non-isothermal)
    [fugcoef_ref, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP);
    f_ref = fugcoef_ref .* comp_ref * press_ref;
    
    % Initial guess: use isothermal solution
    [comp_guess, press_guess, ~] = main(h, comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP, M, h_ref, pressb_ini);
    initial_guess = [comp_guess; press_guess];
    
    % Solve non-isothermal system
    fun = @(x) residual_fugacity_nonisothermal(x(1:n), x(end), f_ref, comp_ref, press_ref, temp_ref, temp_h, h, h_ref, Pc, Tc, acentric, BIP, M, molecular_weights, Cp_coeffs, H_ig_ref_per_mass);
    options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
    
    try
        solution = fsolve(fun, initial_guess, options);
        comp_h = solution(1:n);
        press_h = solution(end);
        
        % Ensure positive compositions and normalize
        comp_h = max(comp_h, 1e-10);
        comp_h = comp_h / sum(comp_h);
    catch
        warning('Non-isothermal fsolve failed at depth %.1f m', h);
        comp_h = comp_ref;
        press_h = press_ref;
    end
    
    % Calculate bubble point pressure at final temperature
    try
        [pressb_h, ~] = pressbub_multicomp_newton(comp_h, pressb_ini, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        
        % Validate bubble point result
        if ~isreal(pressb_h) || pressb_h <= 0 || ~isfinite(pressb_h)
            pressb_h = NaN;
        end
    catch
        pressb_h = NaN;
    end
end

function F = residual_fugacity_nonisothermal(comp_h, press_h, f_ref, comp_ref, press_ref, temp_ref, temp_h, h, h_ref, Pc, Tc, acentric, BIP, M, molecular_weights, Cp_coeffs, H_ig_ref_per_mass)
    % Non-isothermal residual function - calculates thermal terms within fsolve iterations
    
    R = 8.3144598;
    g = 9.80665;
    n = length(f_ref);
    
    % Calculate current fugacity coefficients
    [fugcoef_h, ~] = fugacitycoef_multicomp(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP);
    
    % Calculate absolute enthalpies at reference state (UPDATED CALL)
    [H_abs_ref, H_abs_mix_ref, ~, ~, ~] = calculate_absolute_enthalpy(temp_ref, press_ref, comp_ref, Pc, Tc, acentric, BIP, molecular_weights, Cp_coeffs, H_ig_ref_per_mass);
    
    % Calculate absolute enthalpies at current state (UPDATED CALL)
    [H_abs_h, H_abs_mix_h, ~, ~, ~] = calculate_absolute_enthalpy(temp_h, press_h, comp_h, Pc, Tc, acentric, BIP, molecular_weights, Cp_coeffs, H_ig_ref_per_mass);
    
    % Calculate average molecular weights
    M_avg_ref = sum(comp_ref .* M);
    M_avg_h = sum(comp_h .* M);
    
    % Calculate thermal terms for each component
    delta_T = temp_h - temp_ref;
    F = zeros(n+1, 1);
    
    for i = 1:n
        % Gravity term
        gravity_term = M(i) * g * (h - h_ref) / (R * temp_h);
        
        % Thermal term: M_i * [(H_i_abs/M_i) - (H_mix_abs/M_avg)] * ΔT / (R*T^2)
        % PHYSICS CORRECTION: Use correct sign for thermal diffusion
        ref_specific_enthalpy = H_abs_mix_ref / M_avg_ref;
        component_specific_enthalpy = H_abs_h(i) / M(i);
        thermal_term = M(i) * (component_specific_enthalpy - ref_specific_enthalpy) * delta_T / (R * temp_h^2);
        

        % thermal_term = M(i) * ( ref_specific_enthalpy -  component_specific_enthalpy  ) * delta_T / (R * (temp_ref)^2 );
        
    
        % Non-isothermal fugacity balance:
        % ln(φᵢʰ zᵢʰ Pʰ) - ln(φᵢʰ° zᵢʰ° Pʰ°) = gravity_term - thermal_term
        ln_f_current = log(comp_h(i) * press_h * fugcoef_h(i));
        ln_f_ref = log(f_ref(i));
        
        F(i) = ln_f_current - ln_f_ref - gravity_term + thermal_term;
    end
    
    % Mole fraction constraint
    F(n+1) = sum(comp_h) - 1;
end


