function [p,S,mu] = temp(x,y,n)

[p,S,mu] = polyfit(x,y,n);
out = {p,S,mu};
end

