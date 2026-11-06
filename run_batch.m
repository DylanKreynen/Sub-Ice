clear all
close all
clc

%% Sub-Ice: run map_multi_channel.m for a batch of configuration files
%  Runs map_multi_channel.m once for every configuration file listed below
%  (or found in config_dir), e.g. to map channels on different ice shelves
%  or DEMs, or with different settings. If a run fails, the error is
%  reported and the batch continues with the next configuration file.
%
%  Configuration files have the same format as config.m. Relative paths in
%  them (input data, output directories) are relative to the Sub-Ice folder
%  (the folder this script is run from), wherever the configuration files
%  are stored. Use a different proj_subdir and/or file_prefix in each
%  configuration file to avoid overwriting output.
%
%  (c) Dylan Kreynen
%  University of Oslo
%  2024-2026


%% user specifiable variables
%  (update as required)

% configuration files to run: list them explicitly, e.g.
% config_fns = ["./config/config_shelf_A.m", "./config/config_shelf_B.m"];
% or leave empty ([]) to run all files in config_dir matching config_pattern
config_fns = [];
config_dir = './config/';
config_pattern = 'config*.m';


%% list configuration files

if isempty(config_fns)
    listing = dir(fullfile(config_dir, config_pattern));
    config_fns = string(fullfile({listing.folder}, {listing.name}));
end
no_config = numel(config_fns);
if no_config == 0
    error(append("No configuration files found (", fullfile(config_dir, config_pattern), "). "))
end
disp(append("Found ", string(no_config), " configuration file(s) to run. "))


%% run map_multi_channel.m for each configuration file

run_ok = false(no_config, 1);
run_msg = strings(no_config, 1);
for i = 1:no_config
    disp(append("Configuration file ", string(i), "/", string(no_config), ": ", config_fns(i), ". "))
    [run_ok(i), run_msg(i)] = run_one(config_fns(i));
    close all
end

disp("Batch finished. ")
disp(table(config_fns(:), run_ok, run_msg, 'VariableNames', {'config_file', 'success', 'error_message'}))


%% local functions

function [ok, msg] = run_one(config_fn)
% run map_multi_channel.m with the given configuration file, in this
% function's own workspace (map_multi_channel.m starts with "clear all",
% which would otherwise clear the batch list above)
    ok = false;
    msg = "";
    if ~isfile(config_fn)
        msg = "configuration file not found";
        disp(append("Error: ", msg, ". Skipping. "))
        return
    end
    try
        setappdata(groot, 'subice_config_fn', char(config_fn))
        map_multi_channel
        % (variables of this function were cleared by map_multi_channel.m)
        ok = true;
        msg = "";
    catch err
        ok = false;
        if isappdata(groot, 'subice_config_fn')
            rmappdata(groot, 'subice_config_fn')
        end
        msg = string(err.message);
        disp(append("Error: ", msg, " Continuing with next configuration file. "))
    end
end
