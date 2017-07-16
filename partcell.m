function [Y] = partcell(X,dim,I)
%% [varargout] = partcell(X,dim,I)
% This function goes though the cell X in the specified dimention, and
% attempts to construct a numeric/logical output form each slice specified 
% by I along it. If I is [] then it goes though every slice.
% WIP: This explanation is not very good looking...
%
% X : {100,300,5};
% Y : {5} where each elem is a size:[100 300]

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
%% Figure out the slices
sizeXI = 1:length(sizeX);
sizeY = sizeX(sizeXI ~= dim);

Y = cell(1,length(I));

for Ii = 1:length(I) % For each of the slices along the specified dimention
    i = I(Ii);
    sX.type = '()';
    sX.subs = {};
    for j = 1:length(sizeX)
        sX.subs{end+1} = ':';
    end
    sX.subs{dim} = i;
    
    %keyboard
    slice = subsref(X,sX);
    sizeSlice = size(slice);
    
    %% Check subset of elements for type coherance
    
    % Choose log10(n) + 1 indices without repetition
    ns = numel(slice);
    randomI = randperm(ns,fix(log10(ns)+1)); % Randomly permute indices of the slice
    
    for j = 1:numel(randomI)
        ri = randomI(j);
        [subs{:}] = ind2sub(sizeSlice,ri);
        ris.type = '{}';
        ris.subs = subs;
        keyboard
    end
end

