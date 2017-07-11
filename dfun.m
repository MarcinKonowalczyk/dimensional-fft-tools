function Y = dfun(X,fun,dim,funargs,varargin)
%% Y = dfun(X,fun,dim,funsrgs,...)
% Applies a function specified by the function handle '@fun' to each 1D
% slice of 'X' along the dimention 'dim'.
%
% SYNTAX
%  Y = dfun(X,fun);
%  Y = dfun(X,fun,dim);
%  Y = dfun(X,fun,dim,funargs);
%  Y = dfun(X,fun,dim,funargs,...);
%
%  X       - Input matrix. This can be anything really. Y = X if x is an
%            empty object
%  fun     - Funtion handle to the function to apply to each slice though X
%  dim     - Dimention to slice X by. If this input is left empty then it
%            defaults to the first non-singleton dimention of X it can find
%  funargs - Additional arguments to pass to @fun. For each slice though X,
%            @fun will be called with: `fun(sliceThoughX,funargs{:})`.
%            Regardless of the dimentinoality of X, `sliceThoughX` will
%            allways be a collumn vector.
%
% OPTIONS
%  'verbosity' - Flag which controlls the output of the function to the
%  (true)        console. This can be usefull when the @fun is expected to
%                take a long time, or X is large.
%  'plot'      - Flag which controlls plotting of the slices though X in
%  (flase)       the folr loop. Use only when neccessary - e.g. for
%                debugging your @fun.
%                
% EXAMPLES
%  Work In Progress
%   for j = 1:16;
%       A(:,:,j) = conv2(peaks(64),magic(j),'same');
%   end

%% Parse and validate input
narginchk(2,Inf); % X and @fun are required
nargoutchk(0,1);

if isempty(X) || isempty(fun); Y = X; return; end; % Return Y = X if X or @fun are empty
sizeX = size(X);

assert(isa(fun,'function_handle'),'dfun:invalid-input','''fun'' must be a function handle, not %s',class(fun));

% Default of 'dim': Find first nonsingleton dimention of X
if nargin < 3 || isempty(dim);
    dim = find(sizeX ~= 1,1);
else
    assert(isnumeric(dim),'dfun:invalid-input','''dim'' must be numeric, not a %s');
    assert(isequal(size(dim),[1 1]),'dfun:invalid-input','''dim'' must be a single number, not %s',mat2str(size(dim)));
end

% Default of 'funargs': No additional @fun arguments
if nargin < 4 || isempty(funargs)
    funargs = {};
elseif ~isa(funargs,'cell');
    % Put funargs into a cell for input into feval
    % Nothing to assert about funargs. It can be anything.
    funargs = {funargs};
end

% Parse the options
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'verbosity',true);
addOptional(p,'plot',false);
% WIP: option to time the execution of the @funs ?? (<- not sure if necessary)
parse(p,varargin{:});
opt = p.Results;


%% Stuff
cleaners = {};
warningState = warning('off','backtrace');
cleaners{end+1} = onCleanup(@() warning(warningState)); %#ok<NASGU>

teapot = MException('dfun:E418','I''m a teapot'); % Idiot error. This should never happen.

%% Determine the nature of the output of @fun
if nargout(fun) > 1
    opt.mode = 'cell';
else
    % Apply @fun to a dummy slice
    subs = cell(1,ndims(X)); % Initialise subs
    [subs{:}] = ind2sub(sizeX,randi([1 numel(X)])); % Take a random point in X
    subs{dim} = ':'; % Make it a slice in the correct dimention
    sDummy.type = '()'; sDummy.subs = subs; % Create the subscript struct
    sDO = feval(fun,subsref(X,sDummy),funargs{:}); % Dummy output of @fun; sampleDummyOutput

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
    case {'numeric-one', 'numeric-vector'}
        % Y is a matrix where the dimention specified by 'dim' has the outputs of @fun
        sizeY = sizeX; sizeY(dim) = length(sDO);
        Y = zeros(sizeY,'like',sDO);
    case {'logical-one', 'logical-vector'}
        % Y is a logical array of T/F vector outputs of @fun
        sizeY = sizeX; sizeY(dim) = length(sDO);
        Y = false(sizeY);
    case 'cell'
        % Y is a cell array of the outputs of @fun
        sizeY = sizeX; sizeY(dim) = nargout(fun);
        Y = cell(sizeY);
    otherwise
        throw(teapot);
end

%% Loop though each element of X
sliceDoneDim = 1:length(sizeX); sliceDoneDim(dim) = []; % Dimentions to index sliceDone
sliceDoneSize = sizeX; sliceDoneSize(dim) = []; % Dimentions of sliceDone

% Create sliceDone matrix
if length(sliceDoneSize) > 1
    sliceDone = false(sliceDoneSize);
else
    sliceDone = false(1,sliceDoneSize); % Special case for 1D sliceDone
end

sSliceX.type = '()';
sSliceDone.type = '()';

subs = cell(1,ndims(X)); % Initialise subs

iiDisplayStep = fix(numel(X)./10); % Step to display ~ 20 notifications

for xi = 1:numel(X)
    % Display % done message
    if ~mod(xi,iiDisplayStep) && opt.verbosity, fprintf('%3.3f %% done\n',xi./numel(X).*100), end;
    
    % Get the index of the xi'th element of X
    [subs{:}] = ind2sub(sizeX,xi); % Convert linear index to subscripts
    
    % Check if slice already done
    sSliceDone.subs = subs(sliceDoneDim);
    if subsref(sliceDone,sSliceDone), continue, end;
    
    % Slice X
    sSliceX.subs = subs; sSliceX.subs{dim} = ':';
    slice = subsref(X,sSliceX);
    slice = slice(:); % Ensures that slice is a column vector
    
    % Collect outputs of feval
    switch opt.mode
        % WIP: test test test
        % WIP: Add try, catch mechanics
        case {'numeric-one', 'logical-one'}
            sliceOutput = feval(fun,slice,funargs{:}); % Evaluate the @fun on slice
            sSliceY.type = '()'; sSliceY.subs = subs; sSliceY.subs{dim} = 1; % Create the subscript struct
        case {'numeric-vector', 'logical-vector'}
            sliceOutput = feval(fun,slice,funargs{:}); % Evaluate the @fun on slice
            sSliceY.type = '()'; sSliceY.subs = subs; sSliceY.subs{dim} = ':'; % Create the subscript struct
        case 'cell'
            % https://uk.mathworks.com/matlabcentral/answers/96038-how-can-i-capture-an-unknown-number-of-output-arguments-in-a-cell-array-in-matlab-7-5-r2007b#answer_216954
            sliceOutput = cell(1,nargout(fun)); % Prepare the slice output cell
            [sliceOutput{:}] = feval(fun,slice,funargs{:}); % Evaluate the @fun on slice
            sSliceY.type = '{}'; sSliceY.subs = subs; sSliceY.subs{dim} = nargout(fun); % Create the subscript struct
        otherwise
            throw(teapot);
    end
    
    % Put the outputs into Y
    Y = subsasgn(Y,sSliceY,sliceOutput);

    % Mark slice as done
    sliceDone = subsasgn(sliceDone,sSliceDone,true);
end
end