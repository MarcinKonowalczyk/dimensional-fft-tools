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

valid.bool = @(x) islogical(x) && isequal(size(x),[1 1]); % Is a valid T/F flag
valid.outputNames = {'trim', 'shift', 'fast'};
valid.output = @(x) any(cellfun(@(y) strcmp(x,y),valid.outputNames));

p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'abs',true,valid.bool);
addOptional(p,'output','trim',valid.output);
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

iseven = ~mod(n,2); % If number of fft points even
% WIP: add more modes to opt.trim: 'shift', 'trim', 'fast'
switch opt.output
    case 'trim' % Trim to the positive frequency componets only
        if iseven
            len = n/2+1;
            f = linspace(0,0.5,len)./dt;
        else
            len = ceil(n/2);
            db = 1/n; % delta-bin - bin spacing
            f = linspace(0,0.5 - db/2, len)./dt;
        end
        S.subs{dim} = 1:len;
        Y = 2*subsref(Y,S);
    case 'shift' % Shift the spectrum with fftshift
        len = n;
        if iseven
            % create 1-too-long freqyency vector and trim the end to avoid
            % doubling up on nyquist
            f = linspace(-0.5,0.5,len+1)./dt;
            f(end) = [];
        else
            db = 1/n; % delta-bin - bin spacing
            f = linspace(-0.5 + db/2,0.5 - db/2, len)./dt;
        end
        %f = fftshift(f,2);
        %f(end) = [];
        Y = fftshift(Y,dim);
    case 'fast'
        if nargout > 1 % Make frequency axis only if needed
            if iseven
                % WIP: do this better
                f = linspace(0,1,n+1)./dt;
                f(end) = [];
            else
                db = 1/n; % delta-bin - bin spacing
                f = linspace(0 + db/2,1 - db/2, n)./dt;
            end
        end
    otherwise
        throw(teapot);
end

% Output
if nargout <2
    varargout = {Y};
else % nargout == 2
    varargout = {Y,f};
end
end

function y = fftshift(x,dim)
%% x = fftshift(x,dim)
% In-built fftshift implementation. It is analogous to Matlab's `fftshift`
% function. It shifts zero-frequency component of the spectrum to the
% center.

k = floor(size(x,dim)./2); % Half the length of the spectrum along `dim`
y = circshift(x,k,dim); % Circular shift
end


