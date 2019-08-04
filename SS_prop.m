function [k] = SS_prop(T)

T = T - 273.15; % K to C

if (T >= 21.1) && (T <= 815.6)
    k = 0.015*T+13.99;
end

end

