function [varargout] = dbkg(X,o,dim,varargin)
%% [Y,dY,P,S,Mu] = dbkg(X,o,dim,...)
% This function fits (& subtracts) the n'th order polynomial background
% function from the data along the 'dim' dimention.
%
% SYNTAX
% ...
% 
% EXAMPLES
% ...
%  for j = 1:16;
%      A(:,:,j) = conv2(magic(64),magic(j),'same');
%  end
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

%% Parse input
narginchk(1,Inf);
nargoutchk(0,5);

msgID = 'dbkg:InvalidInput'; % InvalidInput message ID

assert(isnumeric(X),msgID,'`X` must be a numeric array');
sizeX = size(X);

% Default of 'dim': Find first non-singleton dimension of X
if nargin < 3 || isempty(dim)
    dim = find(sizeX ~= 1,1);
else
    assert(isnumeric(dim),msgID,'`dim` must be numeric, not a %s',class(dim));
    assert(isequal(size(dim),[1 1]),msgID,'`dim` must be a single number, not %s',mat2str(size(dim)));
    assert(dim >= 1,msgID,'`dim` supplied (%d) must be >= 1',dim);
    assert(dim == fix(dim),'`dim` supplied(%d) mist be an integer',dim);
    assert(dim <= length(sizeX),msgID,'`dim` supplied (%d) cannot be larger than number of dimentions of X (%d)',dim,length(sizeX));
end

% Default of 'o': 0'th order polynomial - mean of the slice
if nargin < 2 || isempty(0)
    o = 0;
else
    assert(isnumeric(o),msgID,'`o` must be numeric, not a %s',class(o));
    assert(isequal(size(o),[1 1]),msgID,'`o` must be a single number, not %s',mat2str(size(o)));
    assert(o >= 0,msgID,'`o` must be >= 0. A value of %d was supplied.',o);
end

valid.output = @(x) any(cellfun(@(y) strcmp(x,y),{'subtract', 'fit', 'divide'}));
valid.num = @(x) isnumeric(x) && isequal(size(x),[1 1]);

p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'output','subtract',valid.output);
addOptional(p,'plot',false);
addOptional(p,'plotpause',0,valid.num);
parse(p,varargin{:});
opt = p.Results;

%% Stuff
cleaners = {};
warningState = warning('off','backtrace');
cleaners{end+1} = onCleanup(@() warning(warningState)); %#ok<NASGU>

teapot = MException('dbkg:Error418','I''m a teapot'); % Idiot error. This should never happen.

%% Select @fun
% The form of these is [y,dy,p,s,mu] = SubBkgN(s,o,output);
switch nargout % Switch between functions with different N of outputs
    case 1
        fun = @SubBkg1;
    case 2
        fun = @SubBkg2;
    case 3
        fun = @SubBkg3;
    case 4
        fun = @SubBkg4;
    case 5
        fun = @SubBkg5;
    otherwise
        throw(teapot);
end
    

if opt.plot, opt.figure = figure; end

%% Apply dfun and process the outputs
F = dfun(X,fun,dim,{o,opt},'advinput',true);

if nargout >= 2, dY = cell2mat(F(:,2));   end
if nargout >= 3,  P = cell2mat(F(:,3));   end
if nargout >= 4,  S = F(:,4);             end
if nargout >= 5, Mu = cell2mat(F(:,5)')'; end
if nargout >  1,  Y = cell2mat(F(:,1));   end % If nargout == 1, F=F;

switch nargout
    case 1, varargout = {F};
    case 2, varargout = {Y,dY};
    case 3, varargout = {Y,dY,P};
    case 4, varargout = {Y,dY,P,S};
    case 5, varargout = {Y,dY,P,S,Mu};
    otherwise
        throw(teapot);
end
end

function [y,dy,p,s,mu] = SubBkg5(x,S,o,opt)
%% [y,dy,p,s,mu] = SubBkgN(s,o,output)
% Subtract polynomial background (All the outputs)

n = length(x); % Number of opits in the slice
sI = 1:n; % Slice indices

if o >= n, o = n-1; end % Cant fit a polynomial of order >= n of points

[p,s,mu] = polyfit(sI,x,o); % Fit n'th order
[y,dy] = polyval(p,1:n,s,mu);

% WIP: Move plotfun to dfun
if opt.plot
    % Convert subs to a string
    subs = S.subs;
    for i = 1:length(subs)
        if isnumeric(subs{i}), subs{i} = num2str(subs{i}); end
    end
    subs = ['[' strjoin(subs,' ') ']'];
    
    figure(opt.figure);
    plot(sI,x,sI,y,sI,y+dy,'r--',sI,y-dy,'r--');
    grid on;
    xlim([1 n]); set(gca,'XTickLabel','');
    title(sprintf('%s slice though X',subs));
    drawnow;
    pause(opt.plotpause);
end

switch opt.output
    case 'subtract'
        y = x - y;
    case 'divide'
        y = x./y - 1;
        dy = dy./y; % Output fractional error
    case 'fit'
        % y = y;
    otherwise
        teapot = MException('dbkg:Error418','I''m a teapot'); % Idiot error. This should never happen.
        throwAsCaller(teapot);
end

end


function y = SubBkg1(x,S,o,opt)
% Subtract polynomial background (1 output)
[y,~,~,~,~] = SubBkg5(x,S,o,opt);
end

function [y,dy] = SubBkg2(x,S,o,opt)
% Subtract polynomial background (2 outputs)
[y,dy,~,~,~] = SubBkg5(x,S,o,opt);
end

function [y,dy,p] = SubBkg3(x,S,o,opt)
% Subtract polynomial background (3 outputs)
[y,dy,p,~,~] = SubBkg5(x,S,o,opt);
end

function [y,dy,p,s] = SubBkg4(x,S,o,opt)
% Subtract polynomial background (4 outputs)
[y,dy,p,s,~] = SubBkg5(x,S,o,opt);
end