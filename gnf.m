function h = gnf(mode)
%% h = gnf(mode)
% Get New Figure: get a handle number to a new (unused) figure without
% opening it.
%
% SYNTAX
%  h = hfigure;
%  h = hfigure(mode);
%
% EXMAPLE
%  h = hfigure;
%  figure(h);
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

if nargin < 1
    mode = 'next';
end

%% Find all the currently open figures (including the hidden ones)
% https://stackoverflow.com/a/4540637
if verLessThan('matlab','8.4') % MATLAB R2014a and earlier
    O = findall(0, 'Type', 'figure');
    H = O;
else % MATLAB R2014b and later
    O = findall(groot, 'Type', 'figure');
    H = cell(1,length(O));
    if isempty(O)
        H = [];
    else
        [H{:}] = O.Number; H = cell2mat(H);
    end
end

%% Generate a new handle
switch mode
    case 'next' % Next smallest handle
        % https://stackoverflow.com/a/1586940/2531987
        h = find(sort(H) > 1:length(H),1);
        if isempty(h)
            h = length(H) + 1;
        end
    case 'random'
        % This approach can go inot an infinite loop when you have more
        % than N figures oopen. Given that N is large, I think one will
        % have more serious performance problems with N figures open than
        % this function. Anyway though:
        % WIP: Figure out a better way of doing this
        N = 1e3; % Max number of figures
        h = randi([1,N]);
        while any(h == H)
            h = randi(1,N);
        end
end