function [varargout] = dbkg(X,o,dim,varargin)
%% [Y,dY,P,S,Mu] = dbkg(X,o,dim,...)
% This function fits (& subtracts) the n'th order polynomial background function
% from the data along the 'dim' dimention
%
% SYNTAX
%
% EXAMPLES
%
%  for j = 1:16;
%      A(:,:,j) = conv2(magic(64),magic(j),'same');
%  end

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

valid.output = @(x) any(cellfun(@(y) strcmp(x,y),{'subtract', 'fit'}));

p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'output','subtract',valid.output);
addOptional(p,'plot',false);
parse(p,varargin{:});
opt = p.Results;

%% Stuff
cleaners = {};
warningState = warning('off','backtrace');
cleaners{end+1} = onCleanup(@() warning(warningState)); %#ok<NASGU>

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
        fun = @SubBkgN;
    otherwise
        throw(teapot);
end

% Process plot options
if opt.plot
    opt.plot = 'slice+plotfun';
    plotfun = @(x,o,~) SubBkg1(x,o,'fit'); % Use the `fit` output as the @plotfun
else
    opt.plot = 'none';
    plotfun = [];
end

%% Apply dfun and process the outputs
Y = dfun(X,fun,dim,{o,opt.output},'plot',opt.plot,'plotfun',plotfun);

keyboard % <- WIP

switch nargout
    case 1
    case 2
    case 3
    case 4
    case 5
    otherwise
        throw(teapot);
end

varargout = {Y};
if nargout > 1
    order = [dim 1:dim-1 dim+1:length(size(X))]; % Permute dimention of interest to beginning
    iOrder = [2:dim 1 dim+1:order(end)]; % Inverse permute
    sOrder = [1:dim-1 length(size(X)) dim:length(size(X))-1];
    P = permute(P,iOrder);
    varargout{end+1} = P;
end
if nargout > 2
    S = permute(S,sOrder);
    varargout{end+1} = S;
end
if nargout > 3
    Mu = permute(Mu,iOrder);
    varargout{end+1} = Mu;
end

end

function [y,dy,p,s,mu] = SubBkgN(x,o,output)
%% [y,dy,p,s,mu] = SubBkgN(s,o,output)
% Subtract polynomial background (All the outputs)

n = length(x); % Number of opits in the slice
sI = 1:n; % Slice indices

if o >= n, o = n-1; end % Cant fit a polynomial of order >= n of points

[p,s,mu] = polyfit(sI,x,o); % Fit n'th order
[fit,dfit] = polyval(p,1:n,s,mu);

switch output
    case 'subtract'
        y = x - fit;
        dy = dfit;
    case 'divide'
        y = x./fit - 1;
        dy = dfit; % WIP
    case 'fit'
        y = fit;
        dy = dfit;
    otherwise
        % Idiot error. This should never happen.
        teapot = MException('dfun:Error418','I''m a teapot');
        throwAsCaller(teapot);
end
end

function y = SubBkg1(s,o,output)
% Subtract polynomial background (1 output)
[y,~,~,~,~] = SubBkgN(s,o,output);
end

function [y,dy] = SubBkg2(s,o,output)
% Subtract polynomial background (2 outputs)
[y,dy,~,~,~] = SubBkgN(s,o,output);
end

function [y,dy,p] = SubBkg3(s,o,output)
% Subtract polynomial background (3 outputs)
[y,dy,p,~,~] = SubBkgN(s,o,output);
end

function [y,dy,p,s] = SubBkg4(s,o,output)
% Subtract polynomial background (4 outputs)
[y,dy,p,s,~] = SubBkgN(s,o,output);
end