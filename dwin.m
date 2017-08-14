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
valid.windowNames = {'', 'none', 'rect', 'hann', 'hamm'};
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

teapot = MException('dbkg:Error418','I''m a teapot'); % Idiot error. This should never happen.

keyboard

%% Select @fun
% The form of these is [y,dy,p,s,mu] = SubBkgN(s,o,output);
switch nargout % Switch between functions with different N of outputs
    case 1
        fun = @SubBkg1;
    case 2
        fun = @SubBkg2;
    case 3
        fun = @SubBkg3;
    case 4
        fun = @SubBkg4;
    case 5
        fun = @SubBkg5;
    otherwise
        throw(teapot);
end
    

if opt.plot, opt.figure = figure; end

%% Apply dfun and process the outputs
F = dfun(X,fun,dim,{type,opt},'advinput',true);

if nargout >= 2, dY = cell2mat(F(:,2));   end
if nargout >= 3,  P = cell2mat(F(:,3));   end
if nargout >= 4,  S = F(:,4);             end
if nargout >= 5, Mu = cell2mat(F(:,5)')'; end
if nargout >  1,  Y = cell2mat(F(:,1));   end % If nargout == 1, F=F;

switch nargout
    case 1, varargout = {F};
    case 2, varargout = {Y,dY};
    case 3, varargout = {Y,dY,P};
    case 4, varargout = {Y,dY,P,S};
    case 5, varargout = {Y,dY,P,S,Mu};
    otherwise
        throw(teapot);
end
end