classdef FourChamberModelUniversal < handle%Make it like Python classes
    % FourChamberModelUniversal: State-space representation of the 4-Chamber transport model.
    %
    %   This class defines the physics-based grey-box model for particle transport, 
    %   dividing the tokamak into four distinct reservoirs:
    %     1. Core Ions (N_c)
    %     2. Divertor Ions (N_d)
    %     3. Divertor Neutrals (N_{n,d})
    %     4. Main Chamber Neutrals (N_{n,c})
    %
    %   Depending on how configured, the model estimates crucial parameters such as ionization times, confinementtimes, and diffusion coefficients based on experimental input-output data.
    properties
        name = '4 Chamber';

        %% --- Architecture settings ---
        reactor = 'MASTU'; % Options: 'MASTU' or 'TCV'
        N_chain_coreion_divion = 4; % Default number of diffusion steps between the core ions and divertor ions.
        N_chain_divion_divneu = 0;  % Default number of diffusion steps between divertor ions and divertor neutrals.

        midplane_injection = true; % Set to true if you inject through the midplane, if divertor set to false.
        direct_input = true; % Determines if the B matrix directly inserts into the ions (true) or neutrals (false)
        ionisation_splitting = true; % Determines if ionisation splits like the B matrix

        diagnostics_used = 'DAFIG';% Which diagnostics are used. Either DAFIG, FIGFIG, FIGDA or DADA. (TCV only has FIGFIG)

        input_name % Set in constructor
        output_names = {'Core Neutrals', 'Core Ions', 'Divertor Ions', 'Divertor Neutrals'};

        output_weight = diag([1, 1, 5, 1]);% How much weight each output has. Can be "noise", can also be diag([1, 1, 5, 1]); (third channel has 20 times more weight). Or an empty array [] if you want it all equally.
        
        %% --- Parameter settings ---
        % If you don't want to have a certain parameter have effect make it very large here (and then disable it's fitting in the margin settings part.
        %default_parameters = {'Ionization time', 0.003; 'Pumping time 1', 3; 'Pumping time 4', 0.005; 'Leaking time', 0.05*10^10; 'Confinement time', 0.01; 'Recycling time', 0.02; 'gamma', 0.5; 'frac', 0.5 ; 'C1', 300; 'C2', 60; 'C3', 100; 'C4', 500};%No Leackage
        
        %%{
        % Default V2
        default_parameters = { ...
            'Ionization time', 0.003; ... %Make 0.1 if using FIGs 0.003 when using DA/IF
            'Confinement time', 0.05; ...
            'Pumping time 1', 1; ...           
            'Div Ionization time', 0.01*1e9; ... % Can disable/make infinite (maybe see same comment as ion)
            'Recombination time', 0.01; ...
            'Pumping time 4', 0.05; ...
            'Leaking time', 0.05*1e9; ...% Can disable/make infinite
            'k_diff_core_div', 100/16; ... %Can/should be disabled via other settings.
            'k_diff_div_neu', 100/16; ... %Can/should be disabled via other settings.
            'Divertor fraction', 0.5; ...
            'Ionization fraction', 0.5; ... 
            'C1', 300; ...
            'C2', 60; ...
            'C3', 100; ...
            'C4', 500 ...
         };
        %}

        %{
        % Thelosen result
        % Midplane injection MASTU, no diffusion chain, direct_input = true, ionisation_splitting=true
        default_parameters = { ...
            'Ionization time', 0.001; ...      % tau_ion
            'Confinement time', 0.03; ...      % tau_con
            'Pumping time 1', 2.905; ...       % tau_pump1 (Core/Vessel)
            'Div Ionization time', 1e9; ...    % Not present in Thelosen (Infinite)
            'Recycling time', 0.01; ...        % tau_rec
            'Pumping time 4', 0.001; ...       % tau_pump4 (Divertor)
            'Leaking time', 0.0222; ...        % tau_leak (See note above)
            'k_diff_core_div', 100/16; ...        % Unused (N_chain=0)
            'k_diff_div_neu', 100/16; ...         % Unused (N_chain=0)
            'Divertor fraction', 0.25; ...     % 1 - gamma (1 - 0.75)
            'Ionization fraction', 0.7; ...    % 1 - frac (1 - 0.3)
            'C1', 3000; ...                    % Output Scaling 1
            'C2', 60; ...                      % Output Scaling 2
            'C3', 100; ...                     % Output Scaling 3
            'C4', 500 ...                      % Output Scaling 4
         };
        %}

        %{
        % Demonstration of the chain model core->div
        default_parameters = { ...
            'Ionization time', 0.003; ... 
            'Confinement time', 0.01; ...
            'Pumping time 1', 3; ...           
            'Div Ionization time', 0.01; ... 
            'Recombination time', 0.0005; ...% Set ver low to show effect of chain diffusion core ions divertor ions (use with midplane_injection, direct_input, ionisation_splitting alll true)
            'Pumping time 4', 0.005; ...
            'Leaking time', 1e9; ...% infinite/no leak
            'k_diff_core_div', 100/16; ... 
            'k_diff_div_neu', 100/16; ...
            'Divertor fraction', 0.5; ...
            'Ionization fraction', 0.5; ... 
            'C1', 300; ...% Calibration
            'C2', 60; ...
            'C3', 100; ...
            'C4', 500 ...
         };
        %}

        %{
        % Default V1
        default_parameters = { ...
            'Ionization time', 0.003; ... 
            'Confinement time', 0.01; ...
            'Pumping time 1', 3; ...           
            'Div Ionization time', 0.01; ... % Can disable/make infinite
            'Recombination time', 0.02; ...
            'Pumping time 4', 0.03; ...
            'Leaking time', 0.05*1e9; ...% Can disable/make infinite
            'k_diff_core_div', 100/16; ... %Can/should be disabled via other settings.
            'k_diff_div_neu', 100/16; ... %Can/should be disabled via other settings.
            'Divertor fraction', 0.5; ...
            'Ionization fraction', 0.5; ... 
            'C1', 300; ...% Calibration
            'C2', 60; ...
            'C3', 100; ...
            'C4', 500 ...
         };
        %}

        %% --- Margin settings ---
        % If you don't want to fit a certain parameter you can set the margins here to 1. If you don't want it to have any effect you can consider making the value of it very large/small at the parameter settings.

        %%{
        % No leak, no re-ionisation
        default_margin_factors = [ ...
            100, 100; ...   % Ionization time
            100, 100; ...   % Confinement time
            100, 100; ...   % Pumping time 1
            100, 100; ...   % Div Ionization time
            100, 100; ...   % Recycling time
            100, 100; ...   % Pumping time 4
            100, 100; ...   % Leaking time
            100, 100; ...   % k_diff_core_div
            100, 100; ...   % k_diff_div_neu
            inf, 2; ...   % Divertor fraction
            inf, 2; ...   % Ionization fraction
            100, 100; ... % C1
            100, 100; ... % C2
            100, 100; ... % C3
            100, 100 ...  % C4
            ];
        %}

        %{
        % Only do the C matrix.
        default_margin_factors = [ ...
            1, 1; ...   % Ionization time
            1, 1; ...   % Confinement time
            1, 1; ...   % Pumping time 1
            1, 1; ...   % Div Ionization time
            1, 1; ...   % Recycling time
            1, 1; ...   % Pumping time 4
            1, 1; ...   % Leaking time
            1, 1; ...   % k_diff_core_div
            1, 1; ...   % k_diff_div_neu
            1, 1; ...   % Divertor fraction
            1, 1; ...   % Ionization fraction
            100, 100; ... % C1
            100, 100; ... % C2
            100, 100; ... % C3
            100, 100 ...  % C4
            ];
        %}
    end

    methods
        %% Constructor
        function obj = FourChamberModelUniversal()
            if obj.midplane_injection
                obj.input_name = {'Valve midplane'};
            else
                obj.input_name = {'Valve divertor'};
            end
        end
        
        %% Define the matrices of the chamber model
        function A = getAMatrix(obj, t_ion, t_confinement, t_pump1, t_ion_div, t_rec, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4)
            V1 = 1;
            V2 = 1;
            V3 = 1;
            V4 = 1;

            total_states = 4 + obj.N_chain_coreion_divion + obj.N_chain_divion_divneu;
            A = zeros(total_states, total_states);

            % --- Standard Blocks ---
            if obj.direct_input || obj.ionisation_splitting
                A(1:4, 1:4) = ...
                   [-1/t_ion-1/t_pump1-1/(t_leak*V1), 0, 0, 1/(t_leak*V4);%Common Flux core
                    (1-div_frac)/t_ion, -1/t_confinement, 0, 0;% Core plasma
                    div_frac/t_ion, 0, -1/t_rec, 1/t_ion_div;% Scrape off layer / Div Ions 
                    1/(t_leak*V1), 0, 1/t_rec, -1/t_pump4-1/(t_leak*V4)-1/t_ion_div];% Common flux divertor / Div Neutrals
            elseif obj.midplane_injection
                A(1:4, 1:4) = ...
                   [-1/t_ion-1/t_pump1-1/(t_leak*V1), 0, 0, 1/(t_leak*V4);%Common Flux core
                    1/t_ion, -1/t_confinement, 0, 0;% Core plasma
                    0, 0, -1/t_rec, 1/t_ion_div;% Scrape off layer / Div Ions 
                    1/(t_leak*V1), 0, 1/t_rec, -1/t_pump4-1/(t_leak*V4)-1/t_ion_div];% Common flux divertor / Div Neutrals
            else
                A(1:4, 1:4) = ...
                   [-1/t_ion-1/t_pump1-1/(t_leak*V1), 0, 0, 1/(t_leak*V4);%Common Flux core
                    1/t_ion, -1/t_confinement, 0, 1/t_ion_div*(1-div_frac);% Core plasma
                    0, 0, -1/t_rec, 1/t_ion_div*div_frac;% Scrape off layer / Div Ions 
                    1/(t_leak*V1), 0, 1/t_rec, -1/t_pump4-1/(t_leak*V4)-1/t_ion_div];% Common flux divertor / Div Neutrals
            end
            
            % --- Chain 1: Core Ions (2) -> Divertor Ions (3) ---
            if obj.N_chain_coreion_divion == 0
                A(3, 2) = 1/t_confinement;
            else
                A(5, 2) = 1/t_confinement; 
                chain_indices = (1:obj.N_chain_coreion_divion) + 4;
                node_sequence = [chain_indices, 3]; 

                for k = 1:(length(node_sequence) - 1)
                    u = node_sequence(k);     
                    d = node_sequence(k+1);   
                    A(d, u) = A(d, u) + k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                    A(u, u) = A(u, u) - k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                    A(u, d) = A(u, d) + k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                    A(d, d) = A(d, d) - k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                end
            end

            % --- Chain 2: Divertor Ions (3) -> Divertor Neutrals (4) ---
            if obj.N_chain_divion_divneu > 0
                A(3, 3) = A(3, 3) + 1/t_rec; 
                A(4, 3) = A(4, 3) - 1/t_rec; 
                start_idx_c2 = 4 + obj.N_chain_coreion_divion + 1;
                A(3, 3) = A(3, 3) - 1/t_rec;      
                A(start_idx_c2, 3) = 1/t_rec;     

                chain2_indices = (1:obj.N_chain_divion_divneu) + (4 + obj.N_chain_coreion_divion);
                node_sequence_c2 = [chain2_indices, 4]; 

                for k = 1:(length(node_sequence_c2) - 1)
                    u = node_sequence_c2(k);
                    d = node_sequence_c2(k+1);
                    A(d, u) = A(d, u) + k_diff_div_neu*(obj.N_chain_divion_divneu^2);
                    A(u, u) = A(u, u) - k_diff_div_neu*(obj.N_chain_divion_divneu^2);
                    A(u, d) = A(u, d) + k_diff_div_neu*(obj.N_chain_divion_divneu^2);
                    A(d, d) = A(d, d) - k_diff_div_neu*(obj.N_chain_divion_divneu^2);
                end
            end
        end

        function B = getBMatrix(obj, t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4)
            if obj.midplane_injection
                if obj.direct_input
                    B = [1 - ion_frac;              % Core Neutrals
                        (1 - div_frac) * ion_frac;  % Core Ions
                        div_frac * ion_frac;        % Divertor Ions
                        0];                         % Divertor Neutrals
                else
                    B = [1; 0; 0; 0];
                end
            else % Divertor injection
                if obj.direct_input
                    B = [0;                          
                        (1 - div_frac) * ion_frac;   
                        div_frac * ion_frac;         
                        1 - ion_frac];               
                else
                    B = [0; 0; 0; 1];
                end
            end
            
            B = [B; zeros(obj.N_chain_coreion_divion + obj.N_chain_divion_divneu, 1)];
        end

        function C_mat = getCMatrix(obj, t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4)
            C_mat = [C1, 0, 0, 0;
                0, C2, 0, 0;
                0, 0, C3, 0;
                0, 0, 0, C4];
            C_mat = [C_mat, zeros(4, obj.N_chain_coreion_divion + obj.N_chain_divion_divneu)];
        end

        function D = getDMatrix(obj, t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4)
            D = [0; 0; 0; 0];
        end

        function [A, B, C, D] = getMatrices(obj, t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4, Ts)
            % The most important function of this class used by idgrey.
            % Arguments must match the order in default_parameters.
            A = obj.getAMatrix(t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4);
            B = obj.getBMatrix(t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4);
            C = obj.getCMatrix(t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4);
            D = obj.getDMatrix(t_ion, t_confinement, t_pump1, t_ion_div, t_recycling, t_pump4, t_leak, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, C4);
        end

        function [data, settings] = preProccesData(obj, raw_data, settings)
            % preProccesData Cleans, aligns, and normalizes experimental data.
            % ---- Extract data from raw_data ----
            if strcmpi(obj.reactor, 'MASTU')
                try
                    valve_vessel = raw_data.valve.LFSV_BOT_L09;
                catch
                    valve_vessel = [];
                end
                try
                    valve_divertor = raw_data.valve.LFSD_BOT_L0506;
                catch
                    valve_divertor = [];
                end
            
                fig_midplane = raw_data.FIG.HM12;
                IF = raw_data.IF;
                DA_midplane = raw_data.DA.HM10ET;
                L_pol = raw_data.fd;
                DA_divertor = raw_data.DA.HL02SXDT;
                fig_divertor = raw_data.FIG.HL11;
            elseif strcmpi(obj.reactor, 'TCV')
                valve_vessel = raw_data.GV.v1;
                valve_vessel.data = valve_vessel.u_m;
                
                valve_divertor = valve_vessel; 

                fig_midplane.data = raw_data.APG.p_mid;
                fig_midplane.time = raw_data.APG.t;
                
                fig_divertor.data = raw_data.APG.p_div;
                fig_divertor.time = raw_data.APG.t;

                IF.data = raw_data.FIR.ne_ctr;
                IF.time = raw_data.FIR.time_ctr;

                L_pol.data = raw_data.CIII_50.Lf;
                L_pol.time = raw_data.CIII_50.time;
            else
                error("No reactor set in chamber model class")
            end
        
            % ---- Align and prepare data ---
            if obj.midplane_injection 
                if strcmpi(obj.diagnostics_used,'FIGFIG')
                    allignedData = allignData({valve_vessel, fig_midplane, IF, L_pol, fig_divertor}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                elseif strcmpi(obj.diagnostics_used,'DAFIG')
                    allignedData = allignData({valve_vessel, {DA_midplane, IF}, IF, L_pol, fig_divertor}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                elseif strcmpi(obj.diagnostics_used,'FIGDA')
                    allignedData = allignData({valve_vessel, fig_midplane, IF, L_pol, {DA_divertor, IF}}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                elseif strcmpi(obj.diagnostics_used,'DADA')
                    allignedData = allignData({valve_vessel, {DA_midplane, IF}, IF, L_pol, {DA_divertor, IF}}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                else
                    error("Diagnostics configured wrong")
                end
            else 
                if strcmpi(obj.diagnostics_used,'FIGFIG')
                    allignedData = allignData({valve_divertor, fig_midplane, IF, L_pol, fig_divertor}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                elseif strcmpi(obj.diagnostics_used,'DAFIG')
                    allignedData = allignData({valve_divertor, {DA_midplane, IF}, IF, L_pol, fig_divertor}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                elseif strcmpi(obj.diagnostics_used,'FIGDA')
                    allignedData = allignData({valve_divertor, fig_midplane, IF, L_pol, {DA_divertor, IF}}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                elseif strcmpi(obj.diagnostics_used,'DADA')
                    allignedData = allignData({valve_divertor, {DA_midplane, IF}, IF, L_pol, {DA_divertor, IF}}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                else
                    error("Diagnostics configured wrong")
                end
                
                %allignedData = allignData({valve_divertor, fig_midplane, IF, L_pol, fig_divertor}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                %%allignedData = allignData({valve_divertor, fig_midplane, IF, L_pol, {DA_divertor, IF}}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                %allignedData = allignData({valve_divertor, {DA_midplane, IF}, IF, L_pol, {DA_divertor, IF}}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
            end
            
            % ---- (Re-) construct and save normalisation
            num_signals = numel(allignedData);

            normalization_matrix = zeros(num_signals, 2);
            for i = 1:num_signals
                offset = allignedData{i}.normalisation_offset_used;
                scaling = allignedData{i}.normalisation_scaling_used;

                normalization_matrix(i, 1) = offset;
                normalization_matrix(i, 2) = scaling;
            end
            
            settings.normalize_data_used = normalization_matrix;
            settings.normalize_data = settings.normalize_data_used;
            
            % ---- Make idddata object ----
            U_meas = [];
            for i = 1:length(obj.input_name)
                U_meas = [U_meas, allignedData{i}.data];
            end
            U_meas = U_meas';

            Y_meas = [];
            for i = (length(obj.input_name)+1):num_signals
                Y_meas = [Y_meas, allignedData{i}.data];
            end
        
            dt_array = diff(allignedData{1}.time);
            Ts = median(dt_array);
        
            data = iddata(Y_meas, U_meas.', Ts, 'Name', settings.chamber_model.name);
            data.InputName = settings.chamber_model.input_name;
            data.OutputName = settings.chamber_model.output_names;
            data.Tstart = allignedData{1}.time(1);
            data.TimeUnit = 's';

            data = iddata_remove_nan(data);
        end
    end
end