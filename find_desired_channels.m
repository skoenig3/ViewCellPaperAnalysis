function [desired_channels] = find_desired_channels(cfg,channel_type)
% written by Seth Konig August 19, 2014
% determine which channel corresponds to the indeces in the data stucture
% array.
%
% Inputs:
%   1) cfg: configuration and data file data
%   3) channel_type: which type of channel are you looking for i.e. units,
%   LFPs, or eye [data]
%
% updated 10/29/15; input was hdr, data,channel_type
% Code rechecked for bugs October 17, 2016 SDK


% find_desired_channels
mychannels = cfg.channel;
% mychannels = cell(1,length(data));
% for chan = 1:length(data)
%     mychannels{chan} = hdr.label{data(chan).whichchannel};
% end

channel_type = lower(channel_type);

switch channel_type
    case {'sig','units','unit'}
        temp = strfind(mychannels,'sig');
        temp(cellfun(@isempty,temp)) = {0};
        desired_channels = find(cell2mat(temp));
    case {'lfp','lfps'}
        temp = strfind(mychannels,'AD');
        temp(cellfun(@isempty,temp)) = {0};
        desired_channels = find(cell2mat(temp));
    case 'eye'
        temp = strfind(mychannels,'X');
        temp(cellfun(@isempty,temp)) = {0};
        desired_channels = find(cell2mat(temp));
        temp = strfind(mychannels,'Y');
        temp(cellfun(@isempty,temp)) = {0};
        desired_channels(2) = find(cell2mat(temp));
    case 'pupil'
        temp = strfind(mychannels,'pupil');
        temp(cellfun(@isempty,temp)) = {0};
        desired_channels = find(cell2mat(temp));
end
end
