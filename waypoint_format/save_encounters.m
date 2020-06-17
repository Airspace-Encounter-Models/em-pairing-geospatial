function save_encounters(filename, encounters, varargin)
% function save_encounters(filename, encounters, varargin)
%
% This function is intended as a helper to save_scripts and save_waypoints
%
% Inputs:
%       filename: name of output file
%
%       encounters: structure of encounters to save to file specified by
%       "filename".  For a guide to the format of this structure, see
%       the help documentation in load_scripts.m.
%
%       varargin: parameter/value pairs
%
%           +   parameter: numupdatetype
%
%               value (character string): data type of "num_update"
%               variable: the number of updates scripted for some aircraft
%               in some encounter (default is 'uint8')
%
%           +   parameter: append
%
%               value (logical): set to true to append encounters to the
%               end of the file specified by "filename".  Set to false
%               (default) to write a new file with the name specified by
%               "filename".
%
%           +   parameter: floattype
%
%               value (character string): data type of encounter initial
%               conditions and updates (default is 'double')

%% Handle varargin input
opts = inputParser;
opts.addParamValue('numupdatetype', 'uint8', @ischar);
opts.addParamValue('append', false, @islogical);
opts.addParamValue('floattype', 'double', @ischar);
opts.parse(varargin{:});

floattype = opts.Results.floattype;
numupdatetype = opts.Results.numupdatetype;

%% Append encounters to existing scripts file
if opts.Results.append
    fid = fopen(filename, 'r+');
    num_encounters = size(encounters, 2);
    num_ac = size(encounters, 1);
    
    % Error in opening file
    if fid == -1 
        fid = fopen(filename, 'w');
        fwrite(fid, num_encounters, 'uint32');
        fwrite(fid, num_ac, 'uint32');
    else %no error in opening file
        % Go to beginning of file
        frewind(fid); 
        file_num_encounters = fread(fid, 1, 'uint32');
        if isempty(file_num_encounters)
            fid = fopen(filename, 'w');
            fwrite(fid, num_encounters, 'uint32');
            fwrite(fid, num_ac, 'uint32');
        else
            file_num_ac = fread(fid, 1, 'uint32');
            if file_num_ac ~= num_ac
                error('Must have the same number of aircraft when appending encounter files')
            end
            fseek(fid, 0, 'bof');
            totalencounters = num_encounters + file_num_encounters;
            fwrite(fid, totalencounters, 'uint32');
            fseek(fid, 0, 'eof');
        end
    end
else % Overwrite existing file or create new one
    fid = fopen(filename, 'w');
    num_encounters = size(encounters, 2);
    num_ac = size(encounters, 1);
    fwrite(fid, num_encounters, 'uint32');
    fwrite(fid, num_ac, 'uint32');
end

%% Save encounters
for i = 1:num_encounters
    for j = 1:num_ac
        fwrite(fid, encounters(j, i).initial, floattype);
    end
    for j = 1:num_ac
        num_update = size(encounters(j, i).update, 2);
        if ~isequal(num_update, cast(num_update, numupdatetype))
            error('Number of updates is greater than numupdatetype allows: encounter %0.0f, aircraft %0.0f (see JIRA ticket ACASXCSIM-196)', i, j);
        else
            fwrite(fid, num_update, numupdatetype);
        end
        fwrite(fid, encounters(j, i).update, floattype);
    end
end
fclose(fid);
