function [Y,D] = partcell(X,dim,varargin)
%% [varargout] = partcell(X,dim)
% This function goes though the cell X in the specified dimention, and
% attempts to construct a numeric/logical output form each slice specified 
% by I along it. If I is [] then it goes though every slice.
% WIP: This explanation is not very good looking...
%
% X   : {100,300,5,20};
% dim : 3;
% Y   : {5} where each elem is a size:[100 300 squeeze(n m ...) 20] where 
%       n,m,... are the dimentions of the element of X along dim
% WIP: option `nosqueeze`

%% Parse the input
narginchk(1,2);

msgID = 'dbkg:InvalidInput'; % InvalidInput message ID

assert(iscell(X),msgID,'X must be a cell');
sizeX = size(X);

% Default of 'dim': Find first non-singleton dimension of X
if nargin < 2 || isempty(dim)
    dim = find(sizeX ~= 1,1);
else
    assert(isnumeric(dim),msgID,'`dim` must be numeric, not a %s',class(dim));
    assert(isequal(size(dim),[1 1]),msgID,'`dim` must be a single number, not %s',mat2str(size(dim)));
    assert(dim >= 1,msgID,'`dim` must be > 1. A value of %d  was supplied.',dim);
    assert(dim <= length(sizeX),msgID,'`dim` supplied (%d) cannot be larger than number of dimentions of X (%d)',dim,length(sizeX));
end

if nargin < 3 || isempty(I)
    I = 1:sizeX(dim);
else
    % ...
end

valid.bool = @(x) islogical(x) && isequal(size(x),[1 1]); % Is a valid T/F flag

% Input parser
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'verbosity',true,valid.bool);
% WIP: option to time the execution of the @funs ?? (<- not sure if necessary)
parse(p,varargin{:});
opt = p.Results; % Container for user-supplied options
clear varargin

%% Stuff
cleaners = {};
warningState = warning('off','backtrace'); % Set warning backtrace to `off`
cleaners{end+1} = onCleanup(@() warning(warningState)); % Return warnings to previous state when out of scope

teapot = MException('dfun:Error418','I''m a teapot'); % Idiot error. This should never happen.

%% Figure out the slices
sizeXI = 1:length(sizeX);
sizeY = sizeX(sizeXI ~= dim);

Y = cell(1,length(I));

for Ii = 1:length(I) % For each of the slices along the specified dimention
    i = I(Ii);
    sX.type = '()';
    sX.subs = {};
    for is = 1:length(sizeX)
        sX.subs{end+1} = ':';
    end
    sX.subs{dim} = i;
    
    %keyboard
    slice = subsref(X,sX);
    sizeSlice = size(slice);
    
    %% Check subset of elements for type and size consistency
    % Choose log10(<n of elements>) + <n of dimentions> indices without repetition
    ns = numel(slice);
    randomI = randperm(ns,ceil(log10(ns))+length(sizeSlice));
    
    classes = cell(1,numel(randomI));
    sizes = cell(1,numel(randomI));
    for is = 1:numel(randomI) % Loop though trial slice indices
        subs = cell(1,length(sizeSlice)); % Initialise subs
        [subs{:}] = ind2sub(sizeSlice,randomI(is)); % Convert linerar index to subscript
        ris.type = '{}'; ris.subs = subs;
        trialSliceElement = subsref(slice,ris);
        classes{is} = class(trialSliceElement);
        sizes{is} = size(trialSliceElement);
    end
    classInconsistency = any(cellfun(@(x) ~isequal(x,classes{1}),classes));
    sizeInconsistency = any(cellfun(@(x) ~isequal(x,sizes{1}),sizes));
    
    keyboard
    
    % Skip to the next interation if inconsistent
    % WIP on adding 'force' options
    if classInconsistency
        Y{Ii} = [];
        if opt.verbosity, warning(msgID,'Unexpected class inconsistency. Skipping this slice.'); end;
        continue % Skip slice
    else
        % Take the first element of the trial classes as the class of the whole slice
        if isempty(trialSliceElement)
            Y{Ii} = trialSliceElement;
            continue
        elseif isnumeric(trialSliceElement)
            if isequal(size(trialSliceElement),[1 1])
                mode = 'numeric-one'; % single number element
            elseif isvector(trialSliceElement)
                mode = 'numeric-vector'; % 1D element
            else
                mode = 'numeric-matrix'; % Higher dimentional beeing
            end
        elseif islogical(trialSliceElement)
            if isequal(size(trialSliceElement),[1 1])
                mode = 'logical-one'; % single number element
            elseif isvector(trialSliceElement)
                mode = 'logical-vector'; % 1D element
            else
                mode = 'logical-matrix'; % Higher dimentional beeing
            end
        else
            % Unrecognised mode. Don't convert.
            Y{Ii} = slice;
            continue;
        end
    end
    
    if sizeInconsistency
        Y{Ii} = [];
        if opt.verbosity, warning(msgID,'Unexpected size inconsistency. Skipping this slice.'); end
        continue % Skip slice
    else
        % Take the size of the first element of the trial as the size of each element in the slice
        sizeSliceElement = sizes{1};
    end
    
    % Assumig that's handled...
    
    %% Initialise slice matrix
    
    keyboard
    switch mode
        case 'numeric'
        case 'logical'
        otherwise
            throw(teapot);
    end
    
    %% Convert slice to matrix
    for is = 1:numel(slice)
        [subs{:}] = ind2sub(sizeSlice,ri);
    end
    keyboard
end

