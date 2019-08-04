function [Rho,Cp,k,Visc] = LH2_prop(T)

T = T *9/5; % K to R

%% Rho
if (T >= 50) && (T <= 500)
    Rho = 0.2416*T^(-1.1462);
end

%% Cp
if (T >= 50) && (T < 100)
    p1 = -0.0007536;
    p2 = 0.14716;
    p3 = -2.912;
    
    Cp = p1*T^2 + p2*T + p3;
  
elseif (T >= 100) && (T <= 500)
    p1 = -5.1014e-012;
    p2 = 8.5258e-009;
    p3 = -5.4440e-006;
    p4 = 0.0016345;
    p5 = -0.22622;
    p6 = 15.192;
    
    Cp = p1*T^5 + p2*T^4 + p3*T^3 + p4*T^2 + p5*T + p6 ;
end

%% k
if (T >= 50) && (T <= 500)
    p1 = 5.5393e-021;
    p2 = -1.0669e-017;
    p3 = 8.3866e-015;
    p4 = -3.4167e-012;
    p5 = 7.4316e-010;
    p6 = -7.5117e-008;
    p7 = 3.88e-006;
    
    k = p1*T^6 + p2*T^5 + p3*T^4 + p4*T^3 + p5*T^2 + p6*T + p7 ;
end


%% Visc
if (T >= 50) && (T < 125)
    p1 = -1.568e-012;
    p2 = 5.536e-010;
    p3 = -6.426e-008;
    p4 = 2.725e-006;
    
    Visc = p1*T^3 + p2*T^2 +  p3*T + p4 ;
  
elseif (T >= 125) && (T <= 450)
    p1 = 4.4e-010;
    p2 = 2.25e-007;

    Visc = p1*T + p2;
end


%% Unit Conversion
Rho = Rho/(2.20462*0.0254^3); % lb/in^3 to kg/m^3
Cp = Cp * 1055.06 * 2.20462 *9/5; % Btu/lb-R to J/kg-K
k = k*1055.06/0.0254*9/5; % Btu/in-s-R to W/m-K
Visc = Visc/(2.20462*0.0254); % lb/in-s to kg/(m-s)

end

