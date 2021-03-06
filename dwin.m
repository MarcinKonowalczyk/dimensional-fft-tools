function [Y] = dwin(X,type,dim,varargin)
%% [Y] = dwin(X,o,dim,...)
% This function applies a specified window to each 1D slice of the data
% along the 'dim' dimention
%
% SYNTAX
%
% EXAMPLES
%
%  for j = 1:16;
%      A(:,:,j) = conv2(magic(64),magic(j),'same');
%  end
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

%% Parse input
narginchk(1,Inf);
nargoutchk(0,5);

msgID = 'dwin:InvalidInput'; % InvalidInput message ID

assert(isnumeric(X),msgID,'`X` must be a numeric array');
sizeX = size(X);

% Default of 'dim': Find first non-singleton dimension of X
if nargin < 3 || isempty(dim)
    dim = find(sizeX ~= 1,1);
else
    assert(isnumeric(dim),msgID,'`dim` must be numeric, not a %s',class(dim));
    assert(isequal(size(dim),[1 1]),msgID,'`dim` must be a single number, not %s',mat2str(size(dim)));
    assert(dim >= 1,msgID,'`dim` supplied (%d) must be >= 1',dim);
    assert(dim == fix(dim),'`dim` supplied(%d) mist be an integer',dim);
    assert(dim <= length(sizeX),msgID,'`dim` supplied (%d) cannot be larger than number of dimentions of X (%d)',dim,length(sizeX));
end

valid.output = @(x) any(cellfun(@(y) strcmp(x,y),{'subtract', 'fit', 'divide'}));
valid.num = @(x) isnumeric(x) && isequal(size(x),[1 1]);
valid.windowNames = {'', 'none', 'rect', 'trigle','parzen','welch','sine','hann', 'hamm','flat'};
valid.bool = @(x) islogical(x) && isequal(size(x),[1 1]); % Is a valid T/F flag
valid.windowType = @(x) any(cellfun(@(y) strcmp(x,y),valid.windowNames)) || valid.bool(x);

% Default of 'type': 'rect'
if nargin < 2 || isempty(0)
    type = 'rect';
else
    assert(ischar(type) || islogical(type),msgID,'`type` must be string, not a %s',class(type));
    assert(feval(valid.windowType,type),msgID,'Invalid `type`: %s',type);
end

% Handle null cases of `type` i.e. if `type` is a null window or `flase`
nullNames = {'', 'none', 'rect'};
if any(cellfun(@(y) strcmp(type,y),nullNames)) || ( feval(valid.bool,type) && ~type )
    Y = X;
    return
elseif ( feval(valid.bool,type) && type )
    type = 'hann'; % If `type` is true, default to hann window.
end

p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'plot',false);
addOptional(p,'plotpause',0,valid.num);
parse(p,varargin{:});
opt = p.Results;

%% Stuff
cleaners = {};
warningState = warning('off','backtrace');
cleaners{end+1} = onCleanup(@() warning(warningState)); %#ok<NASGU>

teapot = MException('dwin:Error418','I''m a teapot'); % Idiot error. This should never happen.

%% Select @fun
% WIP: add more windows
switch type
    case 'trigle', fun = @trigle;
    case 'parzen', fun = @persistent_parzen;
    case 'welch', fun = @welch;
    case 'sine', fun = @sine;
    case 'hann', fun = @hann;
    case 'hamm', fun = @hamm;
    case 'flat', fun = @flat;
    otherwise
        throw(teapot);
end

if opt.plot
    opt.plot = 'fun';
else
    opt.plot = 'none';
end

%% Apply dfun
Y = dfun(X,fun,dim,[],'plot',opt.plot,'plotpause',opt.plotpause);
end

function xw = trigle(x)
% Triangular window
N = length(x);
n = 0:(N-1);
alpha = (N-1)./2;
w = 1-abs(n./alpha - 1);
cpg = trapz(w)/N;
w = w./cpg;
xw = x.*w;
end

function xw = welch(x)
% Welch window
N = length(x);
n = 0:(N-1);
alpha = (N-1)./2;
w = 1-(n./alpha - 1).^2;
cpg = trapz(w)/N;
w = w./cpg;
xw = x.*w;
end

function xw = sine(x)
% Sine window
N = length(x);
n = 0:(N-1);
alpha = (N-1)./pi;
w = sin(n./alpha);
cpg = trapz(w)/N;
w = w./cpg;
xw = x.*w;
end

function xw = hann(x)
% Hanning window
N = length(x);
n = 0:(N-1);
alpha = (N-1)./(2*pi);
w = 0.5*(1-cos(n./alpha));
cpg = trapz(w)/N;
w = w./cpg;
xw = x.*w;
end

function xw = hamm(x)
% Hamming window
N = length(x);
n = 0:(N-1);
alpha = (N-1)./(2*pi);
w = 0.54 - 0.46*(cos(n./alpha));
cpg = trapz(w)/N;
w = w./cpg;
xw = x.*w;
end

function xw = flat(x)
% Flattop window
% uk.mathworks.com/help/signal/ref/flattopwin.html
N = length(x);
n = 0:(N-1);
alpha = (N-1)./(2*pi);
a = [0.21557895 0.41663158 0.277263158 0.083578947 0.006947368];
w = a(1) - a(2).*cos(n./alpha) + a(3).*cos(2.*n./alpha) - a(4).*cos(3.*n./alpha) + a(5).*cos(4.*n./alpha);
cpg = trapz(w)/N;
xw = x.*w./cpg;
end

function xw = parzen(x)
% Parzen window
N = length(x);
[xw] = core_parzen(x,N);
end

function xw = persistent_parzen(x)
% Parzen window with persistance
% This save
persistent N;
persistent w;

if isempty(N) || N ~= length(x);
    N = length(x);
else
    xw = x.*w;
    return;
end
[xw,w] = core_parzen(x,N);
end

function [xw,w] = core_parzen(x,N)
% Core of the Parzen window functions
% https://uk.mathworks.com/help/signal/ref/parzenwin.html
n = abs(linspace(-(N-1)/2,(N-1)/2,N));
w = zeros(1,N);

isCenter = n <= (N-1)/4;
w(isCenter) = 1 - 6*(n(isCenter)/N*2).^2 + 6*(n(isCenter)/N*2).^3;
w(~isCenter) = 2*(1-n(~isCenter)/N*2).^3;

cpg = trapz(w)/N;
w = w./cpg;
xw = x.*w;
end
