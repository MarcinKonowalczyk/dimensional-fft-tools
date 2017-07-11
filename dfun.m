function Y = dfun(X,fun,dim,varargin)
%% Y = dfun(X,fun,dim,...)
% This function fits (subtracts) the n'th order polynomial background function
% from the data along the 'dim' dimention
%
% SYNTAX
%
% EXAMPLES
%
%  for j = 1:16;
%      A(:,:,j) = conv2(peaks(64),magic(j),'same');
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

%{
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'output','subtract');
addOptional(p,'plot',0);
parse(p,varargin{:});
opt = p.Results;
%}

%%

cleaners = {};
warningState = warning('off','backtrace');
cleaners{end+1} = onCleanup(@() warning(warningState));

%% Determine the nature of the output of @fun
if nargout(fun) > 1
    opt.mode = 'cell';
else
    % Apply @fun to a dummy slice
    [subs{:}] = ind2sub(sizeX,randi([1 numel(X)])); % Take a random point in X
    subs{dim} = ':'; % Make it a slice in the correct dimention
    sDummy.type = '()'; sDummy.subs = subs; % Create the subscript struct
    sDO = feval(fun,subsref(X,sDummy)); % Dummy output of @fun; sampleDummyOutput

    % Determine the output mode
    % WIP: Add support for the output to be larger than a vector in the matrix mode (has to add dimentions to Y)
    if isnumeric(sDO) % Numeric output
        if isequal(size(sDO),[1 1])
            opt.mode = 'numeric-one';
        elseif isvector(sDO);
            opt.mode = 'numeric-vector';
        else
            opt.mode = 'cell';
        end
    elseif islogical(sDO) % Logical output
        if isequal(size(sDO),[1 1])
            opt.mode = 'logical-one';
        elseif isvector(sDO);
            opt.mode = 'logical-vector';
        else
            opt.mode = 'cell';
        end
    else
        opt.mode = 'cell';
    end
end

%% Initialise the output variable according to the opt.mode
switch opt.mode
    case 'numeric-one'
        sizeY = sizeX; sizeY(dim) = 1;
        Y = zeros(sizeY,'like',sDO);
    case 'numeric-vector'
        sizeY = sizeX; sizeY(dim) = length(sDO);
        Y = zeros(sizeY,'like',sDO);
    case 'logical-one'
        sizeY = sizeX; sizeY(dim) = 1;
        Y = false(sizeY);
    case 'logical-vector'
        sizeY = sizeX; sizeY(dim) = length(sDO);
        Y = false(sizeY);
    case 'cell'
        %
        sizeY = sizeX; sizeY(dim) = 1;
        Y = cell(sizeY);
    otherwise
        error('dfun:E418','I''m a teapot');
end
n = size(X,dim);

sDim.type = '()';
sNotDim.type = '()';
subs = cell(1,ndims(X)); % Initialise subs

iiDisplayStep = fix(numel(X)./10); % Step to display ~ 20 notifications

flagWarning = false;
for ii = 1:numel(X) % For each element of X
    [subs{:}] = ind2sub(sizeX,ii); % Convert linear index to subscripts
    sNotDim.subs = subs(notDim);
    
    if ~mod(ii,iiDisplayStep), fprintf('%3.3f %% done\n',ii./numel(X).*100), end;
    if subsref(sliceDone,sNotDim), continue, end;
    
    subs{dim} = ':';
    sDim.subs = subs;
    
    slice = subsref(X,sDim);
    slice = slice(:)';
    
    % Collect outputs of feval
    % (https://uk.mathworks.com/matlabcentral/answers/96038-how-can-i-capture-an-unknown-number-of-output-arguments-in-a-cell-array-in-matlab-7-5-r2007b#answer_216954)
    sDO = cell(1,nargout(fun));
    [sDO{:}] = feval(fun,slice);
    
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



