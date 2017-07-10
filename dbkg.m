function [Y,P,S,Mu] = dwin(X,o,dim)
%% [Y] = dwin(X,n,dim)
% This function subtracts the n'th order polynomial background function
% from the data along the 'dim' dimention
%
% SYNTAX


%% Parse input
narginchk(1,3);
nargoutchk(0,4);

sizeX = size(X);
notDim = 1:length(sizeX); notDim(dim) = [];
notDimSizeX = sizeX; notDimSizeX(dim) = [];

if length(notDimSizeX) > 1
    sliceDone = false(notDimSizeX);
else
    sliceDone = false(1,notDimSizeX);
end

keyboard
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

%{
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'abs',true);
addOptional(p,'trim',true);
parse(p,varargin{:});
opt = p.Results;
%}

Y = zeros(size(X));

n = size(X,dim);

%order = [dim 1:dim-1 dim+1:length(size(X))]; % Permute dimention of interest to beginning
%iorder = [2:dim-1 1 dim:order(1)]; % Inverse permute

S.type = '()';
subs = cell(1,ndims(X)); % Initialise subs
for ii = 1:numel(X) % For each element of X
    [subs{:}] = ind2sub(sizeX,ii); % Convert linear index to subscripts
    subs(dim); % Subscript of the dimention of interest
    S.subs = subs(notDim);
    
    if subsref(sliceDone,S), continue, end;
    
    Y = subsasgn(Y,S,true);
    
    sliceDone = subsasgn(sliceDone,S,true);
end


p = zeros(o+1,n);
for j = 1:numel(X)
    
end





%% Subscript along the specified dimention
S.type = '()';
S.subs = num2cell(repmat(':',1,length(size(X))));

iseven = ~mod(o,2); % Is number of fft points even
if iseven
    len = o/2+1;
    f = linspace(0,0.5,len)./dt;
else
    len = ceil(o/2);
    f = linspace(0,0.5 - 1/o, len)./dt;
end

S.subs{dim} = 1:len;
Y = subsref(Y,S);

if opt.trim
    if iseven
        S.subs{dim} = 2:len-1;
    else
        S.subs{dim} = 2:len;
    end
    Y = subsasgn(Y,S,2.*subsref(Y,S));
end

% Output
if nargout <2
    varargout = {Y};
else % nargout == 2
    varargout = {Y,f};
end

end



