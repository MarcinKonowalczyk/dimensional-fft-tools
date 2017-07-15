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
%{
WIP
sizeX = size(X);
notDim = 1:length(sizeX); notDim(dim) = [];
notDimSizeX = sizeX; notDimSizeX(dim) = [];

if length(notDimSizeX) > 1
    sliceDone = false(notDimSizeX);
else
    sliceDone = false(1,notDimSizeX);
end

if nargin < 3 || isempty(dim)
    % Find first nonsingleton dimention of X
    dim = find(sizeX ~= 1,1);
end

if nargin < 2 || isempty(o)
    % Default to a 0'th order polynomial
    o = 0;
end
%}

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

%% Create @fun

if nargout > 1
    %
    P = zeros([o+1 notDimSizeX]);
end
if nargout > 2
    S = cell(notDimSizeX);
end
if nargout > 3
    Mu = zeros([2 notDimSizeX]);
end

iiDisplayStep = fix(numel(X)./10); % Step to display ~ 20 notifications

Y = fdim(X,fun,dim);

for ii = 1:numel(X) % For each element of X
    [subs{:}] = ind2sub(sizeX,ii); % Convert linear index to subscripts
    sNotDim.subs = subs(notDim);
    
    if ~mod(ii,iiDisplayStep), fprintf('%3.3f %% done\n',ii./numel(X).*100), end
    if subsref(sliceDone,sNotDim), continue, end
    
    subs{dim} = ':';
    sDim.subs = subs;
    
    slice = subsref(X,sDim);
    slice = slice(:)';
    
    [p,s,mu] = polyfit(1:n,slice,o);
    
    if nargout > 1 % Polynomial fit coeffs
        sForP.type = '()';
        sForP.subs = {':' , subs{notDim}};
        P = subsasgn(P,sForP,p);
    end
    if nargout > 2 % Error estimates structure
        sForS.type = '{}';
        sForS.subs = subs(notDim);
        S = subsasgn(S,sForS,s);
    end
    if nargout > 3 %
        sForMu.type = '()';
        sForMu.subs = {':' , subs{notDim}};
        Mu = subsasgn(Mu,sForMu,mu);
    end
    
    % Plot of the fit
    if opt.plot > 1
        figure(1);
        plot(1:n,slice,'.',1:n,fit,'-');
        grid on; drawnow; pause(0.01);
    end
    
    % Assign into Y
    Y = subsasgn(Y,sDim,slice);
    
    % Mark slice as done
    sliceDone = subsasgn(sliceDone,sNotDim,true);
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

function [y,dy,p,s,mu] = SubBkgN(s,o,output)
%% [y,dy,p,s,mu] = SubBkgN(s,o,output)
% Subtract polynomial background
% All the outputs

n = length(s); % Number of opits in the slice
sI = 1:n; % Slice indices

if o >= n, o = n-1; end % Cant fit a polynomial of order >= n of points

[p,s,mu] = polyfit(sI,s,o); % Fir n'th order
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