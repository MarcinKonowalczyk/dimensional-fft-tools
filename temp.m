function [varargout] = temp(x,y,n)

[p,S,mu] = polyfit(x,y,n);
out = {p,S,mu};
end

