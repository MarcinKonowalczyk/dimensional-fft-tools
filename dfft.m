function [varargout] = dfft(t,X,n,dim,varargin)
%% [f,Y] = dfft(t,X,n,dim,...)
% This function takes the array X and applies a fourier transform along
% each 1D slice along the dimention specified by `dim`. It also calculates 
%
% SYNTAX
% 

%% Parse input
narginchk(1,Inf);
nargoutchk(0,5);

msgID = 'dfft:InvalidInput'; % InvalidInput message ID

% Process required arguments
assert(isnumeric(X),msgID,'`X` must be a numeric array');
sizeX = size(X);

assert(isnumeric(t),msgID,'t must be numeric');
assert(isvector(t),msgID,'t must be a vector');
t = t(:); % Make sure t is a column vector

dt = mean(diff(t));
assert(dt>0,'dimensional-fft-tools:Invalid t','''t'' must be monotonically increasing vector');

% Default of 'dim': Find first non-singleton dimension of X
if nargin < 4 || isempty(dim)
    dim = find(sizeX ~= 1,1);
else
    assert(isnumeric(dim),msgID,'`dim` must be numeric, not a %s',class(dim));
    assert(isequal(size(dim),[1 1]),msgID,'`dim` must be a single number, not %s',mat2str(size(dim)));
    assert(dim >= 1,msgID,'`dim` supplied (%d) must be >= 1',dim);
    assert(dim == fix(dim),'`dim` supplied(%d) mist be an integer',dim);
    assert(dim <= length(sizeX),msgID,'`dim` supplied (%d) cannot be larger than number of dimentions of X (%d)',dim,length(sizeX));
end

% Default of 'n': Length of X along dim
if nargin < 3 || isempty(n)
    n = sizeX(dim);
else
    assert(isnumeric(n),msgID,'`n` must be numeric, not a %s',class(n));
    assert(isequal(size(n),[1 1]),msgID,'`n` must be a single number, not %s',mat2str(size(n)));
    assert(n >= 1,msgID,'`n` supplied (%d) must be >= 1',n);
    assert(n == fix(n),'`n` supplied(%d) mist be an integer',n);
end

p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'abs',true);
addOptional(p,'trim',true);
parse(p,varargin{:});
opt = p.Results;

%% Stuff
cleaners = {};
warningState = warning('off','backtrace');
cleaners{end+1} = onCleanup(@() warning(warningState)); %#ok<NASGU>

teapot = MException('dbkg:Error418','I''m a teapot'); % Idiot error. This should never happen.

%% Seleft @fun


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
    f = linspace(0,0.5 - 1/(2*n), len)./dt;
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



