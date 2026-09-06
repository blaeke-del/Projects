function runProgram()
    % -------------------- CONSTANTS --------------------    
    % General Constants
    P_amb = 101325; % Ambient Pressure (Pa)
    R = 8314; % Universal Gas Constant (J/mol/K)
    g = 9.81; % Acceleration of Gravity (m/s^2)
    SF = 1.5; % Safety Factor
    
    % Water Constants
    y = 1.3;
    M = 18.015; % Molecular Mass of Propellant (g/mol)
    p_prop = 1000; % Density of Propellant (kg/m^3)
    
    % 316 Steel Properties
    S = 150e+6; % Yield Strength of Steel (MPa)
    p_stl = 8000; % Density of Steel (kg/m^3)

    % -------------------- INPUTS --------------------
    variables = input(['Input Array: ' ...
                    'Initial Temperature (C), ' ...
                    'Tank Inner Diameter (m), ' ...
                    'Tank Inner Height (m), ' ...
                    'Nozzle Throat Diameter (m), ' ...
                    'Discharge Coefficient, ' ...
                    'Added Weight (N)... ']);
    
    while true
        plot_type = input(['Characteristic to plot: ' ...
                            'Temperature, ' ...
                            'Pressure, ' ...
                            'Thrust, ' ...
                            'Exhaust Velocity, ' ...
                            'Mass Flow Rate, ' ...
                            'Mass...  '], "s");
        
        if ismember(plot_type, ["Temperature", "Pressure", "Thrust", ...
                        "Exhaust Velocity", "Mass Flow Rate", "Mass"])
            break;
        else
            disp("Invalid Input");
        end
    end
                           
    % Initial Temperature
    T_max = variables(1); % Maximum Temperature of Propellant (C)

    % Tank Dimensions
    D_i = variables(2); % Inner Diameter of Pressure Vessel (m)
    H = variables(3); % Inner Height of Pressure Vessel (m)
    
    % Nozzle Characteristics
    D_t = variables(4); % Nozzle Throat Diameter (m)
    C_d = variables(5); % Discharge Coefficient

    % Extra Weight
    W_extra = variables(6); % Extra Added Weight (N)
    
    % -------------------- CALCULATED VALUES --------------------
    % General Values
    P_max = XSteam('psat_T', T_max) * 10^5; % Maximum Pressure of Propellant (Pa)
    R_sp = R / M; % Specific Gas Constant of Water (J/kg/K)

    % Water Characteristics
    V = pi/4 * D_i^2 * H; % Volume of Water (m^3)
    m_max = p_prop * V; % Mass of Water (kg)
    W_prop = m_max * g; % Weight of Water (N)

    % Nozzle Characteristics
    A_t = pi/4 * D_t^2; % Throat Area (m^2)
    PR_chk = (2 / (y+1))^(y / (y-1)); % Nozzle Choke Pressure Ratio
    PR = P_max / P_amb; % Nozzle Design Pressure Ratio
    M_e = sqrt(2 / (y-1) * (PR^((y-1)/y)-1)); % Nozzle Mach Number
    ER = 1/M_e * (2/(y+1) * (1 + ((y-1)/2) * M_e^2))^((y+1)/(2*(y-1))); % Nozzle Expansion Ratio
    
    % Required Wall Thickness
    t_pv = SF * P_max * D_i / (2 * S); % Pressure Vessel Casing Min. Thickness (m)
    std_values = [1/16, 1/8, 3/16, 1/4, 5/16, 3/8, 7/16, 1/2, 9/16, 5/8] / 39.37;
    t_pv_corr = interp1(std_values, std_values, t_pv, 'next'); % Pressure Vessel Casing Corrected Thickness (m)

    % Recalculate if Thin-Walled Approximation Doesn't Apply
    if t_pv_corr / (D_i / 2) > 0.1
        t_pv = (D_i/2)*sqrt((S/SF + P_max) / (S/SF - P_max)) + (D_i/2);
        t_pv_corr = interp1(std_values, std_values, t_pv, 'next'); % Pressure Vessel Casing Corrected Thickness (m)
    end
    
    % Weight Calculations
    t_blkhead = SF * D_i * sqrt(0.3 * P_max / S); % Bulkhead Min. Thickness (m)
    
    V_casing = ((D_i/2 + t_pv_corr)^2 - (D_i/2)^2) * pi * (H + 2*t_blkhead); % Pressure Vessel Casing Volume (m^3)
    V_blkhead = D_i * pi/4 * t_blkhead; % Approximate Bulkhead Volume (m^3)    
    V_pv = V_casing + 2 * V_blkhead; % Total Pressure Vessel Volume (m^3)

    W_pv = V_pv * p_stl; % Pressure Vessel Weight (N)
    W_total = W_pv + W_prop + W_extra; % Total System Weight (N)

    % -------------------- TRANSIENT CALCULATOR --------------------
    disp("Starting transient calculation...")

    % Initializing Variables
    t = 0; % Current Time (s)
    delta_t = 0.1; % Time Step (s)
    T = T_max; % Current Temperature (C)
    m = m_max; % Current Water Mass (kg)
    P_sat = P_max;
    U = m * XSteam('uL_T', T); % Water Internal Energy (kJ)
    
    % Initializing Arrays for Transient Information
    T_values = []; % Temperature
    P_values = []; % Pressure
    F_values = []; % Thrust
    v_values = []; % Exhaust Velocity
    mdot_values = []; % Mass Flow Rate
    m_values = []; % Mass
    t_values = []; % Time
    
    % Residual Temperature-Finding Algorithm
    function delta = residual(T, u, v)
        u_lqd = XSteam('uL_T', T); % Liquid Specific Internal Energy (kJ/kg)
        u_vpr = XSteam('uV_T', T); % Vapor Specific Internal Energy (kJ/kg)
        
        v_lqd = XSteam('vL_T', T); % Liquid Specific Volume (m^3/kg)
        v_vpr = XSteam('vV_T', T); % Vapor Specific Volume (m^3/kg)
        
        x = (u - u_lqd)/(u_vpr - u_lqd); % Vapor Quality
        v_mix = v_lqd + x*(v_vpr - v_lqd); % Specific Volume of Water Mixture (m^3/kg)
        
        delta = v_mix - v; % Minimizing Function
    end
    
    % Transient Thermodynamics Calculator
    while true
        % Calculate Relevant Thermodynamic Values
        P_sat = XSteam('psat_T', T) * 10^5;
        P_e = P_sat * (1 + (y-1)/2 * M_e^2) ^ (-y/(y-1));
        u_lqd = XSteam('uL_T', T);
        u_vpr = XSteam('uV_T', T);
    
        % Calculate Flow Values
        v_e = sqrt(2*y/(y-1)*R_sp*(T+273.15) * (1 - (P_e / P_sat)^((y-1)/y))); % Exhaust Velocity (m/s)
        mdot = C_d * A_t * P_sat * sqrt(y * (2/(y+1))^((y+1)/(y-1)) / (R_sp * T)); % Mass Flow Rate (kg/s)
        F = mdot * v_e + (P_e - P_amb) * A_t;
        
        % Calculate Specific Internal Energy and Volume
        h_vpr = XSteam('hV_T', T); % Water Vapor Enthalpy (kJ/kg)  
        m = m - mdot * delta_t; % Water Mass (kg)
        U = U - mdot*h_vpr*delta_t; % Water Internal Energy (kJ)
        u = U / m; % Water Specific Internal Energy (kJ/kg)
        v = V / m; % Water Specific Volume (m^3/kg)
        
        % Find Temperature by Minimizing Residual
        T_new = fzero(@(T) residual(T,u,v), T);
        
        % Store Relevant Values in Arrays
        T_values(end+1) = T_new;
        P_values(end+1) = P_sat;
        F_values(end+1) = F;
        v_values(end+1) = v_e;
        mdot_values(end+1) = mdot;
        m_values(end+1) = m;
        t_values(end+1) = t;

        % Store Values for Next Iteration
        T = T_new;
        t = t + delta_t;

        % Check for Choked Condition, End if False
        choked = P_sat > 1.5 * P_amb / PR_chk;
        if ~choked
            break;
        end
    end

    % -------------------- PLOTTING --------------------
    disp("Plotting...")

    % Plot Desired Values
    if plot_type == "Temperature"
        values = T_values;
    elseif plot_type == "Pressure"
        values = P_values;
    elseif plot_type == "Thrust"
        values = F_values;
    elseif plot_type == "Exhaust Velocity"
        values = v_values;
    elseif plot_type == "Mass Flow Rate"
        values = mdot_values;
    elseif plot_type == "Mass"
        values = m_values;
    end

    plot(t_values, values, 'g-')
    hold on
    yline(values(end), 'r-', 'Min. Choked Condition', 'LabelHorizontalAlignment', 'left')
    
    % Calculate Impulse & Specific Impulse
    time = linspace(0, t_values(end), 100);
    F_fit = polyfit(t_values, F_values, 3);
    F_approx = @(time) polyval(F_fit, time);
    imp = integral(F_approx, 0, t_values(end));
    isp = imp / W_prop;
    
    % Calculate Hover Time
    m_fit = polyfit(t_values, m_values, 3);
    m_approx = @(time) polyval(m_fit, time);
    W_prop_avg = integral(m_approx, 0, t_values(end)) / t_values(end);
    t_hover = imp / (W_prop_avg + W_pv + 100);

    % Label Plot
    xlabel('Time (s)')
    if plot_type == "Temperature"
        ylabel(plot_type + " (C)")
    elseif plot_type == "Pressure"
        ylabel(plot_type + " (Pa)")
    elseif plot_type == "Thrust"
        ylabel(plot_type + " (N)")
    elseif plot_type == "Exhaust Velocity"
        ylabel(plot_type + " (m/s)")
    elseif plot_type == "Mass Flow Rate"
        ylabel(plot_type + " (kg/s)")
    elseif plot_type == "Mass"
        ylabel(plot_type + " (kg)")
    end
    
    % Find Plot Borders
    ax = gca;
    x_max = ax.XLim(2);
    y_max = ax.YLim(2);

    x_offset = diff(ax.XLim); 
    y_offset = diff(ax.YLim);
    
    % Display Configuration Variables
    text(x_max - x_offset*0.02, y_max - y_offset*0.02, ...
        {['Initial Temperature: ' num2str(T_max) ' C'], ...
        ['Tank Inner Diameter: ' num2str(D_i) ' m'], ...
        ['Tank Inner Height: ' num2str(H) ' m'], ...
        ['Nozzle Throat Diameter: ' num2str(D_t) ' m'], ...
        ['Discharge Coefficient: ' num2str(C_d)], ...
        ['Additional Weight: ' num2str(W_extra) ' N']}, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    
    % Display Impulse & Time Data
    text(x_max - x_offset*0.02, y_max - y_offset*0.4, ...
        {['Total Impulse: ' num2str(imp) ' N-s'], ...
        ['Specific Impulse: ' num2str(isp) ' s'], ...
        ['Hover Duration: ' num2str(t_hover) ' s'], ...
        ['Max Thrust Duration: ' num2str(t_values(end)) ' s']}, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    disp("Finished!")
end