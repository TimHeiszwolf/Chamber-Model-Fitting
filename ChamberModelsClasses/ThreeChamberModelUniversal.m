classdef ThreeChamberModelUniversal < handle%Make it like Python classes
    % ThreeChamberModelUniversal: State-space representation of a 3-Chamber model.
    %
    %   This class defines a physics-based grey-box model for particle transport, 
    %   simplifying the system into three distinct reservoirs by lumping the neutral gas:
    %     1. Core Ions (N_c)
    %     2. Divertor Ions (N_d)
    %     3. Neutrals (N_n) - Lumped main chamber and divertor neutral gas
    %
    %   The model estimates parameters such as core confinement time, recombination, 
    %   and lumped ionization times based on experimental input-output data.
    properties
        name = '3 Chamber';

        %% --- Architecture settings ---
        reactor = 'TCV'; % Options: 'MASTU' or 'TCV'
        N_chain_coreion_divion = 0; % Default number of diffusion steps between the core ions and divertor ions.
        N_chain_divion_neu = 0;  % Default number of diffusion steps between divertor ions and divertor neutrals.

        midplane_injection = false; % Set to true if you inject through the midplane, if divertor set to false.
        direct_input = false; % Determines if the B matrix directly inserts into the ions (true) or neutrals (false)
        ionisation_splitting = false; % Determines if ionisation splits like the B matrix
        
        recycle_divertor_ions_to_neutrals = false; % If true, ions go back to Neutrals. If false, they leave the system.

        input_name % Set in constructor
        output_names = {'Neutrals', 'Core Ions', 'Divertor Ions'};

        output_weight = diag([1, 1, 5]);% How much weight each output has. Can be "noise", can also be diag([1, 1, 5, 1]); (third channel has 20 times more weight). Or an empty array [] if you want it all equally.
        
        %% --- Parameter settings ---
        % Adapted to match 4-Chamber Universal names. 
        % Removed 'Leaking time'.
        
        %{
        % Default V2
        default_parameters = { ...
            'Ionization time', 0.1; ... %Make 0.1 if using FIGs 0.003 when using DA/IF
            'Confinement time', 0.05; ...
            'Pumping time 1', 0.05; ...% Either 0.05 if recycle to neutral, otherwise maybe 1?
            'Recombination time', 0.01; ...      
            'k_diff_core_div', 100/16; ...%Can/should be disabled via other settings.
            'k_diff_div_neu', 100/16; ... %Can/should be disabled via other settings. % Only used if recycling to neutrals is ON
            'Divertor fraction', 0.5; ...
            'Ionization fraction', 0.5; ...
            'C1', 300; ...
            'C2', 60; ...
            'C3', 100 ...
         };
        %}

        %%{
        % Derks et al. 2025 (TCV) settings
        % Note: Requires direct_input = false, ionisation_splitting = true, recycle_divertor_ions_to_neutrals = true and no diffusion chains.
        default_parameters = { ...
            'Ionization time', 0.00625; ...    % tau_i (derived from EIRENE)
            'Confinement time', 0.025; ...     % tau_c
            'Pumping time 1', 0.5; ...         % tau_p
            'Recombination time', 0.001; ...   % tau_r (Fast divertor dynamics)
            'k_diff_core_div', 100/16; ...        % Unused (Set N_chain=0)
            'k_diff_div_neu', 100/16; ...         % Unused (Set N_chain=0)
            'Divertor fraction', 0.77; ...     % gamma (SOL transparency)
            'Ionization fraction', 0.5; ...    % Unused (input goes 100% to neutrals)
            'C1', 300.0; ...                   % Neutrals Scaling (Matches Derks C3)
            'C2', 60.0; ...                    % Core Scaling (Matches Derks C1)
            'C3', 100.0 ...                    % Divertor Scaling (Matches Derks C2)
         };
        %}

        %{
        % Default V1
        default_parameters = { ...
            'Ionization time', 0.003; ...
            'Confinement time', 0.01; ...
            'Pumping time 1', 0.03; ...
            'Recombination time', 0.01; ...      
            'k_diff_core_div', 100/16; ...%Can/should be disabled via other settings.
            'k_diff_div_neu', 100/16; ... %Can/should be disabled via other settings. % Only used if recycling to neutrals is ON
            'Divertor fraction', 0.5; ...
            'Ionization fraction', 0.5; ...
            'C1', 300; ...
            'C2', 60; ...
            'C3', 100 ...
         };
        %}

        %% --- Margin settings ---
        %%{
        default_margin_factors = [ ...
            100, 100; ...   % Ionization time
            100, 100; ...   % Confinement time
            100, 100; ...   % Pumping time 1
            100, 100; ...   % Recombination time
            100, 100; ...   % k_diff_core_div
            1, 1; ...       % k_diff_div_neu
            inf, 2; ...     % Divertor fraction
            inf, 2; ...     % Ionization fraction
            100, 100; ...   % C1
            100, 100; ...   % C2
            100, 100 ...    % C3
            ];
        %}

        %{
        % Locked version
        default_margin_factors = [ ...
            1, 1; ...   % Ionization time
            1, 1; ...   % Confinement time
            1, 1; ...   % Pumping time 1
            1, 1; ...   % Recombination time
            1, 1; ...   % k_diff_core_div
            1, 1; ...       % k_diff_div_neu
            1, 1; ...     % Divertor fraction
            1, 1; ...     % Ionization fraction
            100, 100; ...   % C1
            100, 100; ...   % C2
            100, 100 ...    % C3
            ];
        %}
    end

    methods
        %% Constructor
        function obj = ThreeChamberModelUniversal()
            if obj.midplane_injection
                obj.input_name = {'Valve midplane'};
            else
                obj.input_name = {'Valve divertor'};
            end
        end
        
        %% Define the matrices of the chamber model
        % Input arguments matched to default_parameters list (t_leak removed)
        function A = getAMatrix(obj, t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3)
            V1 = 1;
            
            % Base states: 1=Neutrals, 2=Core Ions, 3=Divertor Ions
            total_states = 3 + obj.N_chain_coreion_divion + obj.N_chain_divion_neu;
            A = zeros(total_states, total_states);

            % --- Standard Blocks (3x3) ---
            % x1: Neutrals
            % x2: Core Plasma (Ions)
            % x3: Divertor Ions
            
            % Neutrals Loss: Ionization + Pumping (Leakage removed)
            A(1,1) = -1/t_ion - 1/t_pump1;
            
            % Core Ions Loss: Confinement
            A(2,2) = -1/t_confinement;
            
            % Divertor Ions Loss: Recombination
            A(3,3) = -1/t_recomb;

            % --- Couplings ---
            if obj.direct_input || obj.ionisation_splitting
                % Split Ionization source based on div_frac
                A(2,1) = (1-div_frac)/t_ion; % To Core Ions
                A(3,1) = div_frac/t_ion;     % To Divertor Ions
            elseif obj.midplane_injection
                % Standard Midplane: All ionization goes to Core, then Diffuses
                A(2,1) = 1/t_ion;
            else
                % Standard Divertor: 
                A(2,1) = 1/t_ion; 
            end
            
            % --- Chain 1: Core Ions (2) -> Divertor Ions (3) ---
            if obj.N_chain_coreion_divion == 0
                % Direct Flow
                A(3, 2) = 1/t_confinement;
            else
                % Diffusion Chain
                % 1. Core (2) -> Start of Chain (4)
                A(4, 2) = 1/t_confinement; 
                
                chain_indices = (1:obj.N_chain_coreion_divion) + 3; % Starts at 4
                node_sequence = [chain_indices, 3]; % Ends at Divertor Ions (3)

                for k = 1:(length(node_sequence) - 1)
                    u = node_sequence(k);     
                    d = node_sequence(k+1);   
                    A(d, u) = A(d, u) + k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                    A(u, u) = A(u, u) - k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                    A(u, d) = A(u, d) + k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                    A(d, d) = A(d, d) - k_diff_core_div*(obj.N_chain_coreion_divion^2);   
                end
            end

            % --- Chain 2: Divertor Ions (3) -> Neutrals (1) ---
            if obj.recycle_divertor_ions_to_neutrals
                if obj.N_chain_divion_neu == 0
                    % Direct recycling: Div Ions (3) -> Neutrals (1)
                    A(1, 3) = A(1, 3) + 1/t_recomb; 
                    % Note: A(3,3) already has -1/t_recomb
                else
                    % Chain recycling
                    start_idx_c2 = 3 + obj.N_chain_coreion_divion + 1;
                    
                    A(start_idx_c2, 3) = 1/t_recomb;     % Gain to start of chain

                    chain2_indices = (1:obj.N_chain_divion_neu) + (3 + obj.N_chain_coreion_divion);
                    node_sequence_c2 = [chain2_indices, 1]; % Ends at Neutrals (1)

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
        end

        function B = getBMatrix(obj, t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3)
            
            if obj.direct_input
                B = [1 - ion_frac;              % Neutrals
                    (1 - div_frac) * ion_frac;  % Core Ions
                    div_frac * ion_frac];       % Divertor Ions
            else
                B = [1; 0; 0];
            end
            
            B = [B; zeros(obj.N_chain_coreion_divion + obj.N_chain_divion_neu, 1)];
        end

        function C_mat = getCMatrix(obj, t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3)
            C_mat = [C1, 0, 0;
                     0, C2, 0;
                     0, 0, C3];
            
            C_mat = [C_mat, zeros(3, obj.N_chain_coreion_divion + obj.N_chain_divion_neu)];
        end

        function D = getDMatrix(obj, t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3)
            D = [0; 0; 0];
        end

        function [A, B, C, D] = getMatrices(obj, t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3, Ts)
            A = obj.getAMatrix(t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3);
            B = obj.getBMatrix(t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3);
            C = obj.getCMatrix(t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3);
            D = obj.getDMatrix(t_ion, t_confinement, t_pump1, t_recomb, k_diff_core_div, k_diff_div_neu, div_frac, ion_frac, C1, C2, C3);
        end

        function [data, settings] = preProccesData(obj, raw_data, settings)
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
                L_pol = raw_data.fd; % Used for Divertor Ions
                DA_midplane = raw_data.DA.HM10ET; 
                DA_divertor = raw_data.DA.HL02SXDT;
                fig_divertor = raw_data.FIG.HL11; 
                
            elseif strcmpi(obj.reactor, 'TCV')
                % TCV Mapping
                valve_vessel = raw_data.GV.v1;
                valve_vessel.data = valve_vessel.u_m;
                
                valve_divertor = valve_vessel; 
                DA_midplane = [];
                DA_divertor = [];

                fig_midplane.data = raw_data.APG.p_mid;
                fig_midplane.time = raw_data.APG.t;
                
                % Use L_pol for Div Ions
                L_pol.data = raw_data.CIII_50.Lf;
                L_pol.time = raw_data.CIII_50.time;
                
                IF.data = raw_data.FIR.ne_ctr;
                IF.time = raw_data.FIR.time_ctr;
                fig_divertor.data = raw_data.APG.p_div;
                fig_divertor.time = raw_data.APG.t;
            else
                error("No reactor set in chamber model class")
            end
            
            % ---- Align and prepare data ---
            
            if obj.midplane_injection
                allignedData = allignData({valve_vessel, fig_midplane, IF, L_pol}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
                %allignedData = allignData({valve_vessel, {{DA_midplane, IF}, IF, L_pol}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
            else 
                allignedData = allignData({valve_divertor, fig_divertor, IF, L_pol}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
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