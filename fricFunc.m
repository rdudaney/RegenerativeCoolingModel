function [F] = fricFunc(f,Re)

F = 4.0*log10(Re*sqrt(f))-0.4-1/sqrt(f);

end
