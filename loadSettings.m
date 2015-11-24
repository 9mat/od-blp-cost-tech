function [ settings ] = loadSettings
%LOADSETTINGS Summary of this function goes here
%   Detailed explanation goes here

% default values
settings.tolercafe = 1e-4;
settings.tolertech = 1e-3;
settings.tolerbert = 1e-5;
settings.maxitercafe = 400;
settings.maxitertech = 400;
settings.maxbertrand = 5000;

fileid = fopen('settings.txt', 'r');
while ~feof(fileid)
    varname = fscanf(fileid, '%s', 1);
    if isempty(varname)
        continue;
    end
    value = fscanf(fileid, '%f', 1);
    settings.(varname) = value;
end
fclose(fileid);
end

