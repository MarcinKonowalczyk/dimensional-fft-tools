function [varargout] = dbkg(X,o,dim,varargin)
%% [Y,P,S,Mu] = dbkg(X,o,dim,...)
% This function fits (subtracts) the n'th order polynomial background function
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
nargoutchk(0,4);

sizeX = size(X);
notDim = 1:length(sizeX); notDim(dim) = [];
notDimSizeX = sizeX; notDimSizeX(dim) = [];

if length(notDimSizeX) > 1
    sliceDone = false(notDimSizeX);
else
    sliceDone = false(1,notDimSizeX);
end

%{
if isequal(X_size,[1 1]);
    if nargout < 2
        varargout = {0};
    else % nargout == 2
        varargout = {X,NaN};
    end
    return;
end
%}

if nargin < 3 || isempty(dim);
    % Find first nonsingleton dimention of X
    dim = find(sizeX ~= 1,1);
end

if nargin < 2 || isempty(o)
    % Default to a 0'th order polynomial
    o = 0;
end

%
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'output','subtract');
addOptional(p,'plot',0);
parse(p,varargin{:});
opt = p.Results;
%

%%

cleaners = {};
warningState = warning('off','backtrace');
cleaners{end+1} = onCleanup(@() warning(warningState));

%%
Y = zeros(sizeX);

n = size(X,dim);

sDim.type = '()';
sNotDim.type = '()';
subs = cell(1,ndims(X)); % Initialise subs

if nargout > 1
    P = zeros([o+1 notDimSizeX]);
end
if nargout > 2
    S = cell(notDimSizeX);
end
if nargout > 3
    Mu = zeros([2 notDimSizeX]);
end

iiDisplayStep = fix(numel(X)./10); % Step to display ~ 20 notifications

for ii = 1:numel(X) % For each element of X
    [subs{:}] = ind2sub(sizeX,ii); % Convert linear index to subscripts
    sNotDim.subs = subs(notDim);
    
    if ~mod(ii,iiDisplayStep), fprintf('%3.3f %% done\n',ii./numel(X).*100), end;
    if subsref(sliceDone,sNotDim), continue, end;
    
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
    
    [fit,~] = polyval(p,1:n,s,mu); % <- WIP: delta
    
    % Plot of the fit
    if opt.plot > 1
        figure(1);
        plot(1:n,slice,'.',1:n,fit,'-');
        grid on; drawnow; pause(0.01);
    end
    
    % Subtract(/divide/fit) background
    switch opt.output
        case 'subtract'
            slice = slice - fit;
        case 'divide'
            slice = slice./fit - 1;
        case 'fit'
            slice = fit;
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



