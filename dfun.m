function Y = dfun(X,fun,dim,funargs,varargin)
%% Y = dfun(X,fun,dim,funargs,...)
% Applies a function specified by the function handle '@fun' to each 1D
% slice of 'X' along the dimention 'dim'.
%
% SYNTAX
%  Y = dfun(X,fun);
%  Y = dfun(X,fun,dim);
%  Y = dfun(X,fun,dim,funargs);
%  Y = dfun(X,fun,dim,funargs,...);
%
% INPUT
%  X       - Input matrix. This can be anything really. Y = X if x is an
%            empty object
%  fun     - Function handle to the function to apply to each slice though X
%  dim     - Dimension to slice X by. If this input is left empty then it
%            defaults to the first non-singleton dimension of X it can find
%  funargs - Additional arguments to pass to @fun. For each slice though X,
%            @fun will be called with: `fun(sliceThoughX,funargs{:})`.
%            Regardless of the dimensionality of X, `sliceThoughX` will
%            always be a row vector.
% OUTPUT
%  Y       - Depending on @fun, this is either a numeric matrix - if @fun's
%            output is a single number or a vector (preserves class), or a
%            cell array otherwise - if @fun's output is a multidimensional
%            matrix, or if it has more than one output.
%
% OPTIONS
%  Options are given as name-value pairs into varargin.
%   'verbosity' - Flag which controls the output of the function to the
%   (false)       console. This can be useful when the @fun is expected to
%                 take a long time, or X is large.
%   'plot'      - Controls plotting of the slices though X in the for loop.
%   (false)       Use only when necessary - e.g. for debugging your @fun.
%                 The allowed values can be found below in the section
%                 PLOTTING.
%   'plotfun'   - An additional funciton handle to use with some plotting
%   ([])          modes (see PLOTTING below)
%   'plotpause' - Delay (in secconds) between individual slice plots when
%   (0)           using the `plot` mode.
%   'plotfargs' - Equivalent of the `funargs` input, but for use with the
%   (funargs)     @plotfun. The default value is the same as funargs. You
%                 can set it to `{}` fo the @plotfun to recieve no
%                 additonal arguments.
%   'advinput'  - Controls the nature of input into @fun. If this option is
%   (false)       set to 'true', the input into @fun changes to
%                 `fun(sliceThoughX,S,funargs{:})`, where `S` is a structure
%                 describing the access method to the particular slice
%                 though X. See documentation for `subsref` for details.
%                 This should be used with caution as it changes the format
%                 to the @fun required by the function.
%
% PLOTTING
%  dfun allows one to plot the slices though X as slicing takes place. This
%  is a tool for debugging your @fun and getting a better idea of the
%  nature of your data, rather than for `everyday` use. The avalibe
%  plotting modes are:
%   'none'          - Dont plot anything.
%   'slice'         - Plot each slice.
%   true / false    - Translated to 'slice' and 'none' respectivelly.
%   ''              - Translated to 'none'.
%   'fun'           - Plot the @fun only.
%   'plotfun'       - Plot the additional @plotfun only.
%   'slice+fun'     - Plot each slice and @fun overlayed on the same graph
%   'slice+plotfun' - Plot each slice and @plotfun overlayed on the same graph
% 
%  The delay between the plots can be controlled with the 'plotpause'
%  option. The 'plotfun' is not evaluated if it is not needed for plotting.
%  If 'advinput' is true, @plotfun also recieves an additional input
%  argument.
%
% EXAMPLE 1 - Near IR spectra (2D matrix) presentation using function handles
%  load spectra % loads matlab sample data
%  xBkg = [36 97 199 329]; % x for background correction
%  NIR = dfun(NIR,@(y) y - polyval(polyfit(xBkg,y(xBkg),2),1:length(y)),2);
%  dNIR = dfun(NIR,@(x) x - NIR(1,:),2); % Difference spectrum
%  rmsNIR = dfun(dNIR,@(x) sqrt(mean(x.^2)),1);
%  % Plot
%  figure(1); subplot(2,2,1); contourf(NIR,'edgecolor','none'); grid on; title('NIR');
%  set(gca,'XTickLabel','','YTickLabel','');
%  subplot(2,2,2); contourf(dNIR,'edgecolor','none'); grid on; title('\deltaNIR');
%  set(gca,'XTickLabel','','YTickLabel','');
%  subplot(2,2,[3,4]); plot(rmsNIR); grid on; axis tight; title('NIR spectral rms');
%  set(gca,'XTickLabel','');
%
% EXAMPLE 2 - When possible, use build-in functions
%   load spectra % loads matlab sample data
%   tic; NIR_std_1 = std(NIR,1,2); t_base = toc;
%   tic; NIR_std_2 = dfun(NIR,@std,2,{1},'verbosity',false); t_dfun = toc;
%   fprintf('The result of dfun is correct: %i\n',isequal(NIR_std_1,NIR_std_2));
%   fprintf('dfun is ~%.1f times slower than in-build std() function\n',t_dfun./t_base);
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

%% Parse and validate input
narginchk(2,Inf); % X and @fun are required
nargoutchk(0,1);

msgID = 'dfun:InvalidInput'; % InvalidInput message ID

assert(~isa(X,'cell'),msgID,'`X` musn''t be a cell array');
if isempty(X) || isempty(fun); Y = X; return; end % Return Y = X if X or @fun are empty
sizeX = size(X);

assert(isa(fun,'function_handle'),msgID,'@fun must be a function handle, not %s',class(fun));
nargoutFun = nargout(fun); % Number of outputs of the original @fun
assert(nargoutFun ~= 0,msgID,'@fun must be a function which provides outputs');

% Default of 'dim': Find first non-singleton dimension of X
if nargin < 3 || isempty(dim)
    dim = find(sizeX ~= 1,1);
else
    assert(isnumeric(dim),msgID,'`dim` must be numeric, not a %s',class(dim));
    assert(isequal(size(dim),[1 1]),msgID,'`dim` must be a single number, not %s',mat2str(size(dim)));
    assert(dim >= 1,msgID,'`dim` must be > 1. A value of %d  was supplied.',dim);
end

% Default of 'funargs': No additional @fun arguments
if nargin < 4 || isempty(funargs)
    funargs = {};
elseif ~isa(funargs,'cell')
    % Put funargs into a cell for input into feval
    % Nothing to assert about funargs. It can be anything.
    funargs = {funargs};
end

%% Parse the options
% Validators
valid.bool = @(x) islogical(x) && isequal(size(x),[1 1]); % Is a valid T/F flag
valid.num = @(x) isnumeric(x) && isequal(size(x),[1 1]);
valid.plotNames = {'', 'none', 'slice', 'fun', 'plotfun', 'slice+fun', 'slice+plotfun'};
valid.plot = @(x) any(cellfun(@(y) strcmp(x,y),valid.plotNames)) || valid.bool(x);
valid.fun = @(x) isa(x,'function_handle') || isempty(x);

% Input parser
p = inputParser;
p.KeepUnmatched = true;
addOptional(p,'verbosity',false,valid.bool);
addOptional(p,'plot','none',valid.plot);
addOptional(p,'plotfun',[],valid.fun);
addOptional(p,'plotfargs',funargs,@iscell);
addOptional(p,'plotpause',0,valid.num);
addOptional(p,'advinput',false,valid.bool);
% WIP: option to time the execution of the @funs ?? (<- not sure if necessary)
parse(p,varargin{:});
opt = p.Results; % Container for user-supplied options
flag = struct; % Container for inferred options

clear varargin

%% Process plot options
if isempty(opt.plot), opt.plot = 'none'; end % Coerce empty opt.plot to 'none'
% Coersce T/F plot options to `slice` and `none` respectively
if valid.bool(opt.plot)
    if opt.plot
        opt.plot = 'slice';
    else
        opt.plot = 'none';
    end
end
flag.plot = ~strcmp(opt.plot,'none'); % Whether to plot

% Figure out whether meaningful @plotfun exists
% WIP: I think the code below can be written in a neater way
flag.plotfun = flag.plot && ~isempty(opt.plotfun); % Plot and @plotfun supplied
if flag.plotfun % Check whether plotfun is needed by opt.plot
    flag.plotfun = any(cellfun(@(y) strcmp(opt.plot,y),{'plotfun','slice+plotfun'}));
    if ~flag.plotfun && opt.verbosity, warning(msgID,'@plotfun is supplied but unused'); end
end

clear valid p

%% Stuff
cleaners = {};
warningState = warning('off','backtrace'); % Set warning backtrace to `off`
cleaners{end+1} = onCleanup(@() warning(warningState)); % Return warnings to previous state when out of scope

teapot = MException('dfun:Error418','I''m a teapot'); % Idiot error. This should never happen.

%% Create internal @fun (and @plotfun if needed)
% Coerce nargoutFun
if nargoutFun == -1
    % @fun is probably a function handle. These will always give only one output.
    % It's also possible @fun's only output is 'varargout'. Therefore issue a warning about that.
    nargoutFun = 1;
    if opt.verbosity; warning(msgID,'Assuming that @fun has only one output argument'); end
elseif nargoutFun < -1
    % @fun has 'varargout'
    % Number of arguments out minus varargout
    nargoutFun = abs(nargoutFun + 1);
    if opt.verbosity; warning(msgID,'Ignoring `varargout` in the output of the @fun. This may cause errors.'); end
elseif nargoutFun == 0
    throw(teapot);
end

% Redefine @fun for internal use
if opt.advinput
    fun = @(x,S) feval(fun,x,S,funargs{:});
else
    fun = @(x) feval(fun,x,funargs{:});
end

if flag.plotfun
    % Coersce nargoutPlotFun (analogous to above)
    switch nargout(opt.plotfun)
        case -1
            %nargoutPlotFun = 1;
            if opt.verbosity; warning(msgID,'Assuming that @plotfun has only one output argument'); end
        case -2
            %nargoutPlotFun = 1; % Number of arguments out minus varargout
            if opt.verbosity; warning(msgID,'Ignoring `varargout` in the output of the @plotfun. This may cause errors.'); end
        otherwise
            % Abort plot since plotting with invalid @plotfun requested
            flag.plotfun = false; flag.plot = false;
            if ~flag.plotfun && opt.verbosity, warning(msgID,'@plotfun supplied must have exactly one oputput. No plot will be shown.'); end
    end
    
    % Define @plotfun for internal use
    if opt.advinput
        plotfun = @(x,S) feval(opt.plotfun,x,S,opt.plotfargs{:});
    else
        plotfun = @(x) feval(opt.plotfun,x,opt.plotfargs{:});
    end
end

%% Determine the nature of the output of @fun (and @plotfun if needed)
if nargoutFun > 1
    flag.mode = 'cell';
else
    % Apply @fun to a dummy slice
    subs = cell(1,ndims(X)); % Initialise subs
    [subs{:}] = ind2sub(sizeX,randi([1 numel(X)])); % Take a random point in X
    subs{dim} = ':'; % Make it a slice in the correct dimension
    sDummy.type = '()'; sDummy.subs = subs; % Create the subscript struct
    % Dummy output of @fun; sampleDummyOutput
    sliceDummy = subsref(X,sDummy); sliceDummy = sliceDummy(:)';
    if opt.advinput
        subs{dim} = []; subs = cell2mat(subs);
        sDO = fun(sliceDummy,subs);
    else
        sDO = fun(sliceDummy);
    end
    
    % Determine the output mode
    flag.mode = slice2mode(sDO);
end

if flag.plotfun % @plotfun has only one output and is needed
    % Apply @plotfun to a dummy slice
    if opt.advinput
        sDO2 = plotfun(sliceDummy,subs);
    else
        sDO2 = plotfun(sliceDummy);
    end
    
    % Determine the plot mode
    flag.plotmode = slice2mode(sDO2);
    npc = strcmp(flag.plotmode,'cell'); % No Plot Condition
    npc = npc || isempty(sDO2);
    % Check for cases where the plot would contain just one point
    npc = npc || strcmp(opt.plot,'fun') && any(cellfun(@(y) strcmp(flag.mode,y),{'number-one', 'logical-one'}));
    npc = npc || strcmp(opt.plot,'plotfun') && any(cellfun(@(y) strcmp(flag.plotmode,y),{'number-one', 'logical-one'}));
    
    if npc
        if opt.verbosity, warning(msgID,'@plotfun supplied has invalid output. No plot will be shown.'); end
        flag.plot = false; flag.plotfun = flase;
    end
end

%% Initialise the output variable according to the flag.mode
switch flag.mode
    case {'numeric-one', 'numeric-vector'}
        % Y is a matrix where the dimension specified by `dim` has the outputs of @fun
        sizeY = sizeX; sizeY(dim) = length(sDO);
        Y = zeros(sizeY,'like',sDO);
    case {'logical-one', 'logical-vector'}
        % Y is a logical array of T/F vector outputs of @fun
        sizeY = sizeX; sizeY(dim) = length(sDO);
        Y = false(sizeY);
    case 'cell'
        % Y is a cell array of the outputs of @fun
        sizeY = sizeX; sizeY(dim) = nargoutFun;
        Y = cell(sizeY);
    otherwise
        throw(teapot);
end

clear sDO sDO2 sizeY

%% Loop though each element of X
sliceDoneDim = 1:length(sizeX); sliceDoneDim(dim) = []; % Dimensions to index sliceDone
sliceDoneSize = sizeX; sliceDoneSize(dim) = []; % Dimensions of sliceDone

% Create sliceDone matrix
if length(sliceDoneSize) > 1
    sliceDone = false(sliceDoneSize);
else
    sliceDone = false(1,sliceDoneSize); % Special case for 1D sliceDone
end

% Prepare slicers
sSliceX.type = '()';
sSliceDone.type = '()';

% Initialise subs
subs = cell(1,ndims(X));

% Open a new figure if needed
if flag.plot
    fh = figure;
    % Name the figure if running MATLAB R2014b or later
    if ~verLessThan('matlab','8.4'), fh.Name = 'dfun plot'; end
end

% Prepare for verbose output
if opt.verbosity
    iiDisplayStep = fix(numel(X)./100); % Step to display ~ 20 notifications
    stringFun = func2str(fun);
    if ~strcmp(stringFun(1),'@'), stringFun = ['@' stringFun]; end % Make sure stringfun starts with `@`
    fprintf('fdim running with @fun = %.30s: 000%%',stringFun); % Limit the function name to be only 30 characters long
    cleaners{end+1} = onCleanup(@() fprintf('\n')); %#ok<NASGU> % Print a new line character when out of scope
end

for xi = 1:numel(X)
    % Display % done message
    if opt.verbosity && ~mod(xi,iiDisplayStep), fprintf('\b\b\b\b%03.f%%',xi./numel(X).*100), end
    
    % Get the index of the xi'th element of X
    [subs{:}] = ind2sub(sizeX,xi); % Convert linear index to subscripts
    
    % Check if slice already done
    sSliceDone.subs = subs(sliceDoneDim);
    if subsref(sliceDone,sSliceDone), continue, end
    
    % Slice X
    sSliceX.subs = subs; sSliceX.subs{dim} = ':';
    slice = subsref(X,sSliceX);
    slice = slice(:)'; % Ensures that slice is a row vector
    
    % Collect outputs of feval
    % https://uk.mathworks.com/matlabcentral/answers/96038-how-can-i-capture-an-unknown-number-of-output-arguments-in-a-cell-array-in-matlab-7-5-r2007b#answer_216954
    % WIP: Add try, catch mechanics
    % WIP: Catch warnings and throw them only once. Print the progress bar correctly.
    sliceOutput = cell(1,nargoutFun); % Prepare the sliceOutput cell to catch the outputs fo feval
    if opt.advinput
        [sliceOutput{:}] = fun(slice,sSliceX); % Evaluate the @fun on slice
    else
        [sliceOutput{:}] = fun(slice); % Evaluate the @fun on slice
    end
    
    slicePlotOutput = [];
    if flag.plotfun
        if opt.advinput
            slicePlotOutput = plotfun(slice,sSliceX); % Evaluate the @fun on slice
        else
            slicePlotOutput = plotfun(slice); % Evaluate the @fun on slice
        end
    end
    
    % Prepare to slice into Y
    switch flag.mode
        % WIP: test test test
        case {'numeric-one', 'logical-one'}
            sliceOutput = sliceOutput{1}; % Only one output
            sSliceY.type = '()'; sSliceY.subs = subs; sSliceY.subs{dim} = 1; % Create the subscript struct
        case {'numeric-vector', 'logical-vector'}
            sliceOutput = sliceOutput{1}; % Only one output
            sSliceY.type = '()'; sSliceY.subs = subs; sSliceY.subs{dim} = ':'; % Create the subscript struct
        case 'cell'
            sSliceY.type = '{}'; sSliceY.subs = subs; sSliceY.subs{dim} = nargoutFun; % Create the subscript struct
        otherwise
            throw(teapot);
    end
    
    % Slice the outputs into Y
    Y = subsasgn(Y,sSliceY,sliceOutput);
    
    % Mark slice as done
    sliceDone = subsasgn(sliceDone,sSliceDone,true);
    
    % Plot
    if flag.plot, plotSlice(fh,slice,sliceOutput,slicePlotOutput,opt,flag,sSliceX); end
end
end

function mode = slice2mode(y)
%% mode = slice2mode(y)
% Converts a (dummy) slice to the appropriate `mode` string
% WIP: Add support for the output to be larger than a vector in the matrix mode (has to add dimensions to Y)

if isnumeric(y) % Numeric output
    if isequal(size(y),[1 1])
        mode = 'numeric-one';
    elseif isvector(y)
        mode = 'numeric-vector';
    else
        mode = 'cell';
    end
elseif islogical(y) % Logical output
    if isequal(size(y),[1 1])
        mode = 'logical-one';
    elseif isvector(y)
        mode = 'logical-vector';
    else
        mode = 'cell';
    end
else
    mode = 'cell';
end
end

function plotSlice(fh,slice,sliceOutput,slicePlotOutput,opt,flag,sSliceX)
%% dfun_plot(fh,slice,sliceOutput,opt,flag,sSliceX)
% Plot the slice to the figure specified by fh

n = length(slice);
teapot = MException('dfun:Error418','I''m a teapot');

%% Colors
colors.slice = [.1 .1 .1];
colors.true  = [.1 .9 .1];
colors.false = [.9 .1 .1];
colors.fun   = [.1 .1 .9];

%% Convert subs to a string
subs = sSliceX.subs;
for i = 1:length(subs)
    if isnumeric(subs{i}), subs{i} = num2str(subs{i}); end
end
subs = ['[' strjoin(subs,' ') ']'];

%% Switch for different plot modes
figure(fh);

switch opt.plot
    case 'slice'
        X = 1:length(slice);
        plot(X,slice,'color',colors.slice);
        title(sprintf('%s slice through X',subs));
    case {'fun', 'plotfun'} % Case for both @fun and @plotfun only
        if strcmp(opt.plot,'fun')
            Y = sliceOutput; X = 1:length(sliceOutput);
            mode = flag.mode;
            name = '@fun';
        else
            Y = slicePlotOutput; X = 1:length(slicePlotOutput);
            mode = flag.plotmode;
            name = '@plotfun';
        end
        switch mode
            case 'logical-vector'
                tempY = double(plot(1:n,Y));
                % Plot the `true` and 'false' part separatelly
                plot(X(Y),tempY(Y),'.','color',colors.true); hold on;
                plot(X(~Y),tempY(~Y),'.','color',colors.false); hold off;
                ylim([-0.5 1.5]); % Set Y limit to see the T/F vector well
                set(gca,'YTick',[0 1],'YTickLabel',{'T' 'F'});
            case 'numeric-vector'
                plot(X,Y,'color',colors.fun);
            otherwise
                throw(teapot);
        end
        title(sprintf('%s of the %s''th slice through X',name,subs));
    case {'slice+fun', 'slice+plotfun'}
        X = 1:length(slice);
        p1 = plot(X,slice,'color',colors.slice); hold on;
        if strcmp(opt.plot,'slice+fun')
            Y = sliceOutput;
            mode = flag.mode;
            name = '@fun';
        else
            Y = slicePlotOutput;
            mode = flag.plotmode;
            name = '@plotfun';
        end
        switch mode
            case 'numeric-one'
                plot([min(X) max(X)],[Y Y],'--','color',colors.fun);
            case 'numeric-vector'
                X = 1:length(Y);
                plot(X,Y,'color',colors.fun);
            case 'logical-one'
                if Y
                    set(p1,'color',colors.true);
                else
                    set(p1,'color',colors.false);
                end
            case 'logical-vector'
                if length(slice) == length(Y)
                    % The function to plot is a logical of the same length as the slice
                    plot(X(Y),slice(Y),'o','color',colors.true);
                    plot(X(~Y),slice(~Y),'o','color',colors.false);
                else
                    % Scale Y to be visible
                    tempY = double(Y); tempY = tempY .* std(slice) + mean(slice);
                    X = 1:length(Y);
                    plot(X(Y),tempY(Y),'.','color',colors.true);
                    plot(X(~Y),tempY(~Y),'.','color',colors.false);
                end                
            otherwise
                throw(teapot);
        end
        hold off;
        title(sprintf('%s slice through X + %s',subs,name));
        legend('slice',name);
    otherwise
        throw(teapot);
end
grid on; xlim([1 n]); set(gca,'XTickLabel','');
drawnow; pause(opt.plotpause);
end