function [Y] = dwin(X,win,dim)
%% [Y] = dwin(X,dim)
% This function appl
%
% SYNTAX


%% Parse input
narginchk(2,Inf);
nargoutchk(0,2);

X_size = size(X);
if isequal(X_size,[1 1]);
    if nargout <2
        varargout = {X};
    else % nargout == 2
        varargout = {X,NaN};
    end
    return;
end

if nargin < 4 || isempty(dim);
    % Find first nonsingleton dimention of X
    dim = find(X_size ~= 1,1);
end

if nargin < 3 || isempty(n)
    n = size(X,dim);
end

t = t(:);
dt = mean(diff(t));
assert(dt>0,'afft:Invalid t','''t'' must be monotonically increasing vector');

p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'abs',true);
addOptional(p,'trim',true);
parse(p,varargin{:});
opt = p.Results;

%% FFT core
Y = fft(X,n,dim);

%% Take the absolute value
if opt.abs
    Y = abs(Y);
end

%% Subscript along the specified dimention
S.type = '()';
S.subs = num2cell(repmat(':',1,length(size(X))));

iseven = ~mod(n,2); % Is number of fft points even
if iseven
    len = n/2+1;
    f = linspace(0,0.5,len)./dt;
else
    len = ceil(n/2);
    f = linspace(0,0.5 - 1/n, len)./dt;
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



