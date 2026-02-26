%% GetData.m - Experimental Data Ingestion and Preparation
% This script fetches raw experimental data (e.g., from MAST-U) for a specific 
% shot, calculates the fundamental excitation frequencies and periods required 
% for frequency-domain analysis (LPM), and saves the packaged data.
% It is designed to work on the DIFFER Rekenserver.

%% Start clean
clc; close all; clearvars;
%% Add functions and data path
mydir  = which('GetData');                                         %present folder path
idcs   = strfind(mydir,filesep);                                            %find filesep indices
maindir = mydir(1:idcs(end-1));%"/home/emc/heiszwolf/Downloads/mastu_rstbx-OWN/";% MASTU-RSTBX location
addpath(genpath(maindir));                                                  %add full main folder
addpath(genpath('/Differ/Data/MAST-U'));                                    %add datafolder


%% User Settings
% Choose Shot
% Known shots. Not all might work, but good good start.
%MU02 https://wiki.differ.nl/mediawiki/index.php/MASTU_-_MU02_-_EXH04/RT22-05
%shot = 47080;settings.valveplenumpressure=750;

%shot = 47116;settings.valveplenumpressure=750;
shot = 47118;settings.valveplenumpressure=750;
%shot = 47119;settings.valveplenumpressure=750;

%MU03 https://wiki.differ.nl/mediawiki/index.php/MASTU_-_MU03_-_RT22-05_%2B_EXH08_%2B_SUP
%shot = 49297;settings.valveplenumpressure=750;
%shot = 49298;settings.valveplenumpressure=750;
%shot = 49299;settings.valveplenumpressure=750;

%shot = 49300;settings.valveplenumpressure=750;

%shot = 49301;settings.valveplenumpressure=750;
%shot = 49302;settings.valveplenumpressure=750;

%MU04 https://wiki.differ.nl/mediawiki/index.php/MASTU_-_MU04_-_CTRL-02/DIV-07/WPTE
%shot = 50643;settings.valveplenumpressure=500;
%shot = 50648;settings.valveplenumpressure=500;
%shot = 50654;settings.valveplenumpressure=500;

%shot = 50862;settings.valveplenumpressure=500;
%shot = 50863;settings.valveplenumpressure=500;
%shot = 50864;settings.valveplenumpressure=500;

% Data settings/creation
cut_signal_to_size = true;
data.name = 'Test';
data.version = 1.1;
data.shot = shot;

% Valve settings
settings.valvespec = 'flowrate';
%settings.valvespec = 'measured';
%settings.valvespec = 'requested';
%settings.valvespec = 'sysid';
settings.take_requested_from_measured = 2;
settings.calculate_flowrate = 1;
settings.assume_confirm_fdcut = true;%???


%% Get valve and IF data
valves = ["LFSV_BOT_L09", "HFS_MID_U02", "LFSD_TOP_U0102", "LFSD_BOT_L0506"];
succesfull_valves = '';
last_succesfull_valve = "";% Needed for the other data retrival.

for i = 1:length(valves)
    current_valve = valves{i};
    
    try
        settings.valve = string(current_valve);% String needed to convert to string scaler
        settings.take_requested_from_measured = true;
    
        input_type = 'valve';
        output_type = 'IF';% Use this as output because it seems most reliable? Another advantage is that we only needs this data once.
        settings.filter_elms = false;
    
        % Get FRF parameters for this shot
        [f0,f_exc,t1,P,settings] = get_frfpars_mastu(shot, input_type, output_type, settings);
        % Get input and output
        [u,y] = get_inputandoutput(input_type, output_type,shot,settings);

        data.valve.(current_valve) = u;
        data.IF = y;

        data.valve.(current_valve).fs = get_fs(input_type,settings.valvespec,settings.take_requested_from_measured);
        data.IF.fs = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);

        if cut_signal_to_size
            try
                fliptop_output=false;
                fs_u = get_fs(input_type,settings.valvespec,settings.take_requested_from_measured);
                fs_y = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);
                 if fs_u<=fs_y
                     [~,id1]=min(abs(u.time-t1));
                     t1 = u.time(id1);
                 else
                      [~,id1]=min(abs(y.time-t1));
                     t1 = y.time(id1);
                 end
                uc  = frf_cutsignal(u.data,u.time,t1,f_exc(1),fs_u,P);
                uc.time = uc.time';
                uc.data = uc.data';
                uc.data_detrend = uc.data_detrend';
                yc  = frf_cutsignal(y.data,y.time,t1,f_exc(1),fs_y,P,fliptop_output);
                yc.time = yc.time';
                yc.data = yc.data';
                yc.data_detrend = yc.data_detrend';
    
                data.valve.(current_valve).cut = uc;
                data.IF.cut = yc;
    
                %[FRF,U,Up,Y,Ydet,Yp] = frf_lpm_main(uc,yc,f0,f_exc,P,false);
            catch ME
                disp(current_valve + " detrend error: " + ME.message)
            end
        end
    
        succesfull_valves = fieldnames(data.valve);
        last_succesfull_valve = current_valve;
    catch ME
        disp(current_valve + " error: " + ME.message)
    end
    
end
if length(succesfull_valves)<1
    error('No valve data retrieved succesfully')
else
    disp(" ")
    disp("Succesfull valves:")
    disp(succesfull_valves)
    disp(" ")
end


%% Get FIG data
figs = ["HM12", "HL11", "HU08"];
succesfull_figs = '';

for i = 1:length(figs)
    current_fig = figs{i};
    
    try
        settings.valve = last_succesfull_valve;
        settings.take_requested_from_measured = true;
        
        settings.FIGlocation = string(current_fig);% String needed to convert to string scaler
        disp(settings.FIGlocation)

        input_type = 'valve';
        output_type = 'FIG';
        settings.filter_elms = false;
    
        % Get FRF parameters for this shot
        [f0,f_exc,t1,P,settings] = get_frfpars_mastu(shot, input_type, output_type, settings);
        % Get input and output
        [u,y] = get_inputandoutput(input_type, output_type, shot, settings);

        % data.valve.(current_valve) = u;
        data.FIG.(current_fig) = y;

        data.FIG.(current_fig).fs = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);

        if cut_signal_to_size
            try
                fliptop_output=false;
                fs_u = get_fs(input_type,settings.valvespec,settings.take_requested_from_measured);
                fs_y = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);
                 if fs_u<=fs_y
                     [~,id1]=min(abs(u.time-t1));
                     t1 = u.time(id1);
                 else
                      [~,id1]=min(abs(y.time-t1));
                     t1 = y.time(id1);
                 end
                uc  = frf_cutsignal(u.data,u.time,t1,f_exc(1),fs_u,P);
                uc.time = uc.time';
                uc.data = uc.data';
                uc.data_detrend = uc.data_detrend';
                yc = frf_cutsignal(y.data,y.time,t1,f_exc(1),fs_y,P,fliptop_output);
                yc.time = yc.time';
                yc.data = yc.data';
                yc.data_detrend = yc.data_detrend';
    
                %data.valve.(current_valve).cut = uc;
                data.FIG.(current_fig).cut = yc;
            catch ME
                disp(current_valve + " detrend error: " + ME.message)
            end
        end

        succesfull_figs = fieldnames(data.FIG);
    catch ME
        disp(current_fig + " error: " + ME.message)
    end
    
end
if length(succesfull_figs)<1
    warning('No fig data retrieved succesfully')
    data.name = append(data.name, "NoFIG");
else
    disp(" ")
    disp("Succesfull figs:")
    disp(succesfull_figs)
    disp(" ")
end



%% Get DA data
DAs = ["HM10ER", "HM10ET",...
    "HU10OSPR", "HU10OSXDT",...
    "HL01SXDR", "HL02OSPR", "HL02SXDT",...
    "HE05ISPR"];
succesfull_DAs = '';

for i = 1:length(DAs)
    current_DA = DAs{i};
    
    try
        settings.valve = last_succesfull_valve;

        %settings.valvespec = 'flowrate';
        %settings.valvespec='measured';
        %settings.valvespec='requested';
        %settings.valvespec='sysid';
        settings.take_requested_from_measured = true;
        
        settings.scope = string(current_DA);% String needed to convert to string scaler

        input_type = 'valve';
        output_type = 'DA';
        settings.filter_elms = false;
    
        % Get FRF parameters for this shot
        [f0,f_exc,t1,P,settings] = get_frfpars_mastu(shot, input_type, output_type, settings);
        % Get input and output
        [u,y] = get_inputandoutput(input_type, output_type,shot,settings);
    
        data.DA.(current_DA) = y;
        data.FIG.(current_DA).fs = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);

        if cut_signal_to_size
            try
                fliptop_output=false;
                fs_u = get_fs(input_type,settings.valvespec,settings.take_requested_from_measured);
                fs_y = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);
                 if fs_u<=fs_y
                     [~,id1]=min(abs(u.time-t1));
                     t1 = u.time(id1);
                 else
                      [~,id1]=min(abs(y.time-t1));
                     t1 = y.time(id1);
                 end
                uc  = frf_cutsignal(u.data,u.time,t1,f_exc(1),fs_u,P);
                uc.time = uc.time';
                uc.data = uc.data';
                uc.data_detrend = uc.data_detrend';
                yc  = frf_cutsignal(y.data,y.time,t1,f_exc(1),fs_y,P,fliptop_output);
                yc.time = yc.time';
                yc.data = yc.data';
                yc.data_detrend = yc.data_detrend';
    
                %data.valve.(current_valve).cut = uc;
                data.DA.(current_DA).cut = yc;
            catch ME
                disp(current_valve + " detrend error: " + ME.message)
            end
        end

        succesfull_DAs = fieldnames(data.DA);
    catch ME
        disp(current_DA + " error: " + ME.message)
    end
    
end
if length(succesfull_DAs)<1
    warning('No DA data retrieved succesfully')
    data.name = append(data.name, "NoDA");
else
    disp(" ")
    disp("Succesfull DAs :")
    disp(succesfull_DAs)
    disp(" ")
end

%% Get Lpol data (using fd (fast diagnostics))
succesfull_fd = '';
%try
    settings.valve = last_succesfull_valve;
    settings.take_requested_from_measured = true;
    
    input_type ='valve';
    output_type = 'fd';
    settings.filter_elms = false;
    
    settings.fdspec = 'FB_50'; %for MWI Lpol
    % settings.fdspec = 'FB_50_fix'; %for MWI Lpol
    settings.specline=settings.fdspec;
    settings.referenceshot=0; %speficy 0 to disable
    settings.inversion=0;
    settings.fixequil=0;

    [f0,f_exc,t1,P,settings] = get_frfpars_mastu(shot, input_type, output_type, settings);
    % Get input and output
    [u,y] = get_inputandoutput(input_type, output_type,shot,settings);

    % data.valve.(current_valve) = u;
    data.fd = y;
    data.fd.data = data.fd.data.';%Fix oerientation of data.
    data.fd.fs = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);

    if cut_signal_to_size
        try
            fliptop_output=false;
            fs_u = get_fs(input_type,settings.valvespec,settings.take_requested_from_measured);
            fs_y = get_fs(output_type,settings.valvespec,settings.take_requested_from_measured);
             if fs_u<=fs_y
                 [~,id1]=min(abs(u.time-t1));
                 t1 = u.time(id1);
             else
                  [~,id1]=min(abs(y.time-t1));
                 t1 = y.time(id1);
             end
            uc  = frf_cutsignal(u.data,u.time,t1,f_exc(1),fs_u,P);
            uc.time = uc.time';
            uc.data = uc.data';
            uc.data_detrend = uc.data_detrend';
            yc  = frf_cutsignal(y.data,y.time,t1,f_exc(1),fs_y,P,fliptop_output);
            yc.time = yc.time';
            yc.data = yc.data';
            yc.data_detrend = yc.data_detrend';
    
            %data.valve.(current_valve).cut = uc;
            data.fd.cut = yc;
        catch ME
            disp(current_valve + " detrend error: " + ME.message)
        end
    end
    succesfull_fd = 'Yes fd succesfull';
%catch ME
%    disp("FD L_pol error: " + ME.message)
%    warning('No fd/L_pol data retrieved succesfully')
%    data.name = append(data.name, "Nofd");
%end

%% Thomson scattering
% NOTE: this part required a sperate file for the data. 
succesfull_TS = '';
try
    %load("Data/DTS/DTSdata/DTS_47118.mat")
    load(append("Data/DTS/DTSdata/DTS_", string(shot), ".mat"))

    data.TS.n_e.time = n_e.time;
    data.TS.n_e.data = n_e.signal.';%data.TS.n_e.data(1,:)
    data.TS.n_e.dn_e = dn_e.signal.';
    data.TS.n_e.fs = 1/median(diff(n_e.time));

    data.TS.T_e.time = T_e.time;
    data.TS.T_e.data = T_e.signal.';%data.TS.T_e.data(1,:)
    data.TS.T_e.dn_e = dT_e.signal.';
    data.TS.T_e.fs = 1/median(diff(T_e.time));

    succesfull_TS = 'Yes TS succesfull';
catch ME
    disp("TS error: " + ME.message)
    warning('No TS data retrieved succesfully')
    data.name = append(data.name, "NoTS");
end


%% Saving data
disp(" ")
disp("---------------------------------------")
disp(" ")
check="Y";
if (length(succesfull_figs)<1)&&(length(succesfull_DAs)<1)&&(length(succesfull_fd)<1)
    check = input("No FIGs, DAs and fd retrieved, still want to save? (Y/N): ", 's');
end

if strcmpi(check,"Y")
    disp("Succesfull valves:")
    disp(succesfull_valves)
    disp("Succesfull figs:")
    disp(succesfull_figs)
    disp("Succesfull DAs :")
    disp(succesfull_DAs)
    disp("Succesfull fd :")
    disp(succesfull_fd)
    disp("Succesfull TS :")
    disp(succesfull_TS)
    disp(" ")
    disp(fieldnames(data))

    disp('Saving data')
    save(append(append(append(string(shot), string(settings.valvespec)), data.name),'.mat'))
    %save(append(string(shot), data.name, '.mat'))
    disp('Finished succesfully!')
else
    disp("Not saving data")
end