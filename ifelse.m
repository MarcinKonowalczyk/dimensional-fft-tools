function output = ifelse(condition,iftrue,iffalse)
%% output = ifelse(condition,iftrue,iffalse)
% One-line if/else statement in a function form. This can be useful for
% more versatile one-line implementations using dfun. No input validation
% is done on the input.
%
% Written by Marcin Konowalczyk (although, admitedly, there wasn't much to write here)
% Timmel Group @ Oxford University

narginchk(3,3); % Must have 3 inputs

if condition
    output = iftrue;
else
    output = iffalse;
end
end

