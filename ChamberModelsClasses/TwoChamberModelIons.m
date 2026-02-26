classdef TwoChamberModelIons < handle
    % TwoChamberModelIons: State-space representation of a 2-Chamber ion model.
    %
    %   This class defines a simplified physics-based grey-box model for particle 
    %   transport focusing exclusively on the plasma state, divided into:
    %     1. Core Ions (N_c)
    %     2. Divertor Ions (N_d)
    %
    %   It estimates parameters such as core confinement time and divertor recycling.
    properties
        name = '2 Chamber Ions';

        %% --- Architecture settings ---
        reactor = 'MASTU'; % Options: 'MASTU' or 'TCV'
        N_chain_coreion_divion = 4; % Default number of diffusion steps between the core ions and divertor ions.

        midplane_injection = true; % Set to true if you inject through the midplane, if divertor set to false.
        ionisation_splitting = true; % Determines if the input splits between Core and Divertor based on 'Divertor fraction'
        
        input_name % Set in constructor
        output_names = {'Core Ions', 'Divertor Ions'};

        output_weight = diag([1, 5]); % Higher weight for Divertor Ions typically
        
        %% --- Parameter settings ---
        % Reduced set of parameters (No Ionization or Pumping times)
        % C2 = Core Ions Scaling, C3 = Divertor Ions Scaling (Kept naming for consistency)
        
        default_parameters = { ...
            'Confinement time', 0.05; ...
            'Recombination time', 0.01; ...      
            'k_diff_core_div', 100/16; ... % Only used if N_chain > 0
            'Divertor fraction', 0.5; ...  % Used if ionisation_splitting = true
            'C2', 60; ...                  % Core Scaling
            'C3', 100 ...                  % Divertor Scaling
         };

        %% --- Margin settings ---
        default_margin_factors = [ ...
            100, 100; ...   % Confinement time
            100, 100; ...   % Recombination time
            100, 100; ...   % k_diff_core_div
            inf, 2; ...     % Divertor fraction
            100, 100; ...   % C2
            100, 100 ...    % C3
            ];
    end

    methods
        %% Constructor
        function obj = TwoChamberModelIons()
            if obj.midplane_injection
                obj.input_name = {'Valve midplane'};
            else
                obj.input_name = {'Valve divertor'};
            end
        end
        
        %% Define the matrices of the chamber model
        function A = getAMatrix(obj, t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3)
            % Base states: 1=Core Ions, 2=Divertor Ions
            total_states = 2 + obj.N_chain_coreion_divion;
            A = zeros(total_states, total_states);

            % --- Standard Blocks ---
            % x1: Core Ions
            % x2: Divertor Ions
            
            % Core Ions Loss: Confinement
            A(1,1) = -1/t_confinement;
            
            % Divertor Ions Loss: Recombination
            A(2,2) = -1/t_recomb;

            % --- Couplings: Core -> Divertor ---
            if obj.N_chain_coreion_divion == 0
                % Direct Flow: Core(1) -> Divertor(2)
                A(2, 1) = 1/t_confinement;
            else
                % Diffusion Chain
                % 1. Core (1) -> Start of Chain (3)
                A(3, 1) = 1/t_confinement; 
                
                % Chain indices start at 3
                chain_indices = (1:obj.N_chain_coreion_divion) + 2; 
                node_sequence = [chain_indices, 2]; % Ends at Divertor Ions (2)

                for k = 1:(length(node_sequence) - 1)
                    u = node_sequence(k);     
                    d = node_sequence(k+1);   
                    diff_rate = k_diff_core_div*(obj.N_chain_coreion_divion^2);
                    
                    A(d, u) = A(d, u) + diff_rate;   
                    A(u, u) = A(u, u) - diff_rate;   
                    A(u, d) = A(u, d) + diff_rate;   
                    A(d, d) = A(d, d) - diff_rate;   
                end
            end
        end

        function B = getBMatrix(obj, t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3)
            
            % B-Matrix handles injection location and splitting
            if obj.ionisation_splitting
                % If splitting is ON, we distribute based on div_frac regardless of location
                % (This assumes the user wants to model 'Direct Ionization' into divertor)
                B = [1 - div_frac;   % Core Ions
                     div_frac];      % Divertor Ions
            else
                % If splitting is OFF, stick to standard injection points
                if obj.midplane_injection
                    B = [1;  % Core Ions
                         0]; % Divertor Ions
                else
                    % Divertor Injection
                    B = [0;  % Core Ions
                         1]; % Divertor Ions
                end
            end
            
            % Append zeros for any diffusion chain states
            B = [B; zeros(obj.N_chain_coreion_divion, 1)];
        end

        function C_mat = getCMatrix(obj, t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3)
            % Output 1: Core Ions (Scaled by C2 to match naming convention)
            % Output 2: Divertor Ions (Scaled by C3)
            C_mat = [C2, 0;
                     0, C3];
            
            % Append zeros for diffusion chain states
            C_mat = [C_mat, zeros(2, obj.N_chain_coreion_divion)];
        end

        function D = getDMatrix(obj, t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3)
            D = [0; 0];
        end

        function [A, B, C, D] = getMatrices(obj, t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3, Ts)
            A = obj.getAMatrix(t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3);
            B = obj.getBMatrix(t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3);
            C = obj.getCMatrix(t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3);
            D = obj.getDMatrix(t_confinement, t_recomb, k_diff_core_div, div_frac, C2, C3);
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
            
                IF = raw_data.IF;        % Core Ions Proxy
                L_pol = raw_data.fd;     % Divertor Ions Proxy
                
            elseif strcmpi(obj.reactor, 'TCV')
                valve_vessel = raw_data.GV.v1;
                valve_vessel.data = valve_vessel.u_m;
                
                valve_divertor = valve_vessel; 

                IF.data = raw_data.FIR.ne_ctr;
                IF.time = raw_data.FIR.time_ctr;

                L_pol.data = raw_data.CIII_50.Lf;
                L_pol.time = raw_data.CIII_50.time;
            else
                error("No reactor set in chamber model class")
            end
            
            % ---- Align and prepare data ---
            % We only align Input, Core Ions (IF), and Div Ions (L_pol)
            % We ignore Neutrals data (FIG/DA)
            
            if obj.midplane_injection
                allignedData = allignData({valve_vessel, IF, L_pol}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
            else 
                allignedData = allignData({valve_divertor, IF, L_pol}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
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