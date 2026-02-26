classdef TwoChamberModelDivertor < handle
    % TwoChamberModelDivertor: State-space representation of a 2-Chamber divertor model.
    %
    %   This class defines a physics-based grey-box model isolating the divertor 
    %   region's particle transport, divided into:
    %     1. Divertor Ions (N_d)
    %     2. Divertor Neutrals (N_{n,d})
    %
    %   It estimates localized processes such as divertor recombination and ionization.
    properties
        name = '2 Chamber Divertor';

        %% --- Architecture settings ---
        reactor = 'MASTU'; % Options: 'MASTU' or 'TCV'
        
        % Diffusion chains
        N_chain_neu_ion = 4; % Chain from Neutrals -> Ions (Ionization delay)
        N_chain_ion_neu = 4; % Chain from Ions -> Neutrals (Recycling delay)

        direct_input = true; % If true, enables 'Ionization fraction' to split input between Neutrals and Ions
        
        input_name % Set in constructor
        output_names = {'Divertor Neutrals', 'Divertor Ions'};

        output_weight = diag([1, 5]); % Typically ions (L_pol) need more weight than neutrals (FIG)
        
        %% --- Parameter settings ---
        default_parameters = { ...
            'Div Ionization time', 0.01; ...
            'Recombination time', 0.005; ...
            'Pumping time 4', 0.05; ... % Pumping of Divertor Neutrals
            'Recycling fraction', 0.5; ... % Fraction of recombined ions that become neutrals again
            'k_diff_neu_ion', 100/16; ... % Used if N_chain_neu_ion > 0
            'k_diff_ion_neu', 100/16; ... % Used if N_chain_ion_neu > 0
            'Ionization fraction', 0.1; ... % Used if direct_input = true
            'C3', 100; ... % Divertor Ions Scaling
            'C4', 500 ...  % Divertor Neutrals Scaling
         };

        %% --- Margin settings ---
        default_margin_factors = [ ...
            100, 100; ...   % Div Ionization time
            100, 100; ...   % Recombination time
            100, 100; ...   % Pumping time
            inf, 2; ...     % Recycling fraction
            100, 100; ...   % k_diff_neu_ion
            100, 100; ...   % k_diff_ion_neu
            inf, 2; ...     % Ionization fraction
            100, 100; ...   % C3
            100, 100 ...    % C4
            ];
    end

    methods
        %% Constructor
        function obj = TwoChamberModelDivertor()
            obj.input_name = {'Valve divertor'};
        end
        
        %% Define the matrices of the chamber model
        function A = getAMatrix(obj, t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4)
            % Base states: 1=Divertor Neutrals, 2=Divertor Ions
            total_states = 2 + obj.N_chain_neu_ion + obj.N_chain_ion_neu;
            A = zeros(total_states, total_states);

            % --- Standard Losses ---
            % Neutrals (1): Loss via Ionization and Pumping
            A(1,1) = -1/t_ion_d - 1/t_pump;
            
            % Ions (2): Loss via Recombination
            A(2,2) = -1/t_rec;

            % --- Path 1: Neutrals (1) -> Ions (2) (Ionization) ---
            if obj.N_chain_neu_ion == 0
                % Direct: Neutrals -> Ions
                A(2, 1) = 1/t_ion_d;
            else
                % Diffusion Chain: Neutrals(1) -> Chain -> Ions(2)
                % Chain 1 indices start at 3
                A(3, 1) = 1/t_ion_d; 
                
                chain1_indices = (1:obj.N_chain_neu_ion) + 2; 
                node_sequence1 = [chain1_indices, 2]; 

                for k = 1:(length(node_sequence1) - 1)
                    u = node_sequence1(k);     
                    d = node_sequence1(k+1);   
                    diff_rate = k_neu_ion*(obj.N_chain_neu_ion^2);
                    
                    A(d, u) = A(d, u) + diff_rate;   
                    A(u, u) = A(u, u) - diff_rate;   
                    A(u, d) = A(u, d) + diff_rate;   
                    A(d, d) = A(d, d) - diff_rate;   
                end
            end
            
            % --- Path 2: Ions (2) -> Neutrals (1) (Recycling) ---
            % Flux = f_rec * (1/t_rec) * N_ions
            if obj.N_chain_ion_neu == 0
                % Direct: Ions -> Neutrals
                A(1, 2) = A(1, 2) + f_rec/t_rec;
            else
                % Diffusion Chain: Ions(2) -> Chain -> Neutrals(1)
                % Chain 2 indices start after Chain 1
                start_idx_c2 = 2 + obj.N_chain_neu_ion + 1;
                
                A(start_idx_c2, 2) = f_rec/t_rec; 

                chain2_indices = (1:obj.N_chain_ion_neu) + (2 + obj.N_chain_neu_ion);
                node_sequence2 = [chain2_indices, 1];

                for k = 1:(length(node_sequence2) - 1)
                    u = node_sequence2(k);
                    d = node_sequence2(k+1);
                    diff_rate = k_ion_neu*(obj.N_chain_ion_neu^2);
                    
                    A(d, u) = A(d, u) + diff_rate;   
                    A(u, u) = A(u, u) - diff_rate;   
                    A(u, d) = A(u, d) + diff_rate;   
                    A(d, d) = A(d, d) - diff_rate;   
                end
            end
        end

        function B = getBMatrix(obj, t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4)
            if obj.direct_input
                % Split input between Neutrals (1) and Ions (2)
                B = [1 - f_ion;   % Neutrals
                     f_ion];      % Ions
            else
                % Standard: Input goes to Neutrals
                B = [1; 
                     0];
            end
            
            % Append zeros for diffusion chain states
            B = [B; zeros(obj.N_chain_neu_ion + obj.N_chain_ion_neu, 1)];
        end

        function C_mat = getCMatrix(obj, t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4)
            % Output 1: Neutrals (Scaled by C4)
            % Output 2: Ions (Scaled by C3)
            C_mat = [C4, 0;
                     0, C3];
            
            % Append zeros for diffusion chain states
            C_mat = [C_mat, zeros(2, obj.N_chain_neu_ion + obj.N_chain_ion_neu)];
        end

        function D = getDMatrix(obj, t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4)
            D = [0; 0];
        end

        function [A, B, C, D] = getMatrices(obj, t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4, Ts)
            A = obj.getAMatrix(t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4);
            B = obj.getBMatrix(t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4);
            C = obj.getCMatrix(t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4);
            D = obj.getDMatrix(t_ion_d, t_rec, t_pump, f_rec, k_neu_ion, k_ion_neu, f_ion, C3, C4);
        end

        function [data, settings] = preProccesData(obj, raw_data, settings)
            % ---- Extract data from raw_data ----
            if strcmpi(obj.reactor, 'MASTU')
                try
                    valve_divertor = raw_data.valve.LFSD_BOT_L0506;
                catch
                    valve_divertor = [];
                end
            
                L_pol = raw_data.fd;       % Ions Proxy
                fig_divertor = raw_data.FIG.HL11; % Neutrals Proxy
                
            elseif strcmpi(obj.reactor, 'TCV')
                valve_divertor = raw_data.GV.v1;
                valve_divertor.data = valve_divertor.u_m;

                L_pol.data = raw_data.CIII_50.Lf;
                L_pol.time = raw_data.CIII_50.time;
                
                fig_divertor.data = raw_data.APG.p_div;
                fig_divertor.time = raw_data.APG.t;
            else
                error("No reactor set in chamber model class")
            end
            
            % ---- Align and prepare data ---
            % Order: Input, Neutrals, Ions
            allignedData = allignData({valve_divertor, fig_divertor, L_pol}, settings.normalize_data, settings.detrend_data, settings.filter_data, settings.smooth_data, settings.sample_time, settings.start_time, settings.end_time);
            
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