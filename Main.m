clear
close all
clc

%% Phase I
% Given
F = 16500*4.44822; % N
OF = 5.0;
Pc = 475*101325/14.7; % Pa
N = 180; % Number of tubes
g = 9.81; % m/s^2

%% i
coor =load('NozzleCoordinates2.txt');
z = coor(:,1)*0.0254; %m 
r = coor(:,2)*0.0254; %m 
Per_tube = coor(:,3)*0.0254; %m 
A_tube = coor(:,4)*0.0254^2; %m^2 
d_tube = 4.*A_tube./Per_tube;
t_tube = 0.0115*0.0254; %m 

A = pi*r.^2; % m^2

[A_t,I] = min(A); % Throat Area (m^2)
z_t = z(I); % Location of throat (m^2)

out = 't p visc rho gam mach prand ivac';
[DATA] = func(out,[],A(end)/A_t,OF,Pc);

T0 = DATA(1,1); % K
P0 = DATA(1,2)*10^5; % Pa
gam = DATA(1,5);
Pr = DATA(1,7);
Ivac = DATA(3,8);


for i = 1:length(z)
    if i <= I
        [M(i), T_ratio, P_ratio, RHO_ratio, A_ratio] = flowisentropic(gam, A(i)/A_t, 'sub');
    else
        [M(i), T_ratio, P_ratio, RHO_ratio, A_ratio] = flowisentropic(gam, A(i)/A_t, 'sup');
    end
    
    P(i) = P_ratio * P0; % Pa
    T(i) = T_ratio * T0; % K
end

% figure(1) 
% plot(z,r)
% 
% figure(2) 
% plot(z,M)
% 
% figure(3) 
% plot(z,T)

%% ii

for i = 1:length(z)

    if i == 1
        z1 = z(1);
        z2 = (z(i)+z(i+1))/2;
        r1 = r(1);
        r2 = (r(i)+r(i+1))/2;
    elseif i == length(z)
        z1 = (z(i)+z(i-1))/2;
        z2 = z(end);
        r1 = (r(i)+r(i-1))/2;
        r2 = r(end);
    else
        z1 = (z(i)+z(i-1))/2;
        z2 = (z(i)+z(i+1))/2;
        r1 = (r(i)+r(i-1))/2;
        r2 = (r(i)+r(i+1))/2;
    end
    S(i) = sqrt((z2-z1)^2+(r2-r1)^2);
end

%% Phase 2

P0_fuel = 1000*101325/14.7; % Pa
T0_fuel = 100*5/9; % K
W = F/Ivac; % kg/s
W_fuel_tot = (OF+1)^-1*W;
W_fuel = W_fuel_tot/N;

error=1e-3;
incr = 0.1; % Temperature increment (K)
stop_count = 50000;

out2 = 't mach gam prand rho visc son cond';

for i = length(z):-1:1
    i
    % Previous Temeperature and Pressure
    if i == length(z)
        P_l_prev = P0_fuel;
        T_l_prev = T0_fuel;
        dz = 0;
    else
        P_l_prev = P_tube(i+1);
        T_l_prev = T_tube(i+1);
        dz = z(i+1)-z(i);
    end
    
    % Chamber properties
    if i < I
        [DATA] = func(out2,A(i)/A_t,[],OF,Pc);
        T_inf = DATA(2,1);
        M_inf = DATA(2,2);
        gam_inf = DATA(2,3);
        Pr_inf = DATA(2,4);
        Rho_inf = DATA(2,5);
        Visc_inf = DATA(2,6)*10^-5;
        son_inf = DATA(2,7);
        k_inf = DATA(2,8)/10;
        
        v_inf = M_inf*son_inf;
    else
        [DATA] = func(out2,[],A(i)/A_t,OF,Pc);
        T_inf = DATA(3,1);
        M_inf = DATA(3,2);
        gam_inf = DATA(3,3);
        Pr_inf = DATA(3,4);
        Rho_inf = DATA(3,5);
        Visc_inf = DATA(3,6)*10^-5;
        son_inf = DATA(3,7);
        k_inf = DATA(3,8)/10;
        
        v_inf = M_inf*son_inf;
    end
    
    % LH2 properties (based on previous Temp)
    [rho_l,cp_l,k_l,visc_l] = LH2_prop(T_l_prev);
    v_l = W_fuel/(rho_l*A_tube(i));
    Pr_l = cp_l*visc_l/k_l;
    
    % Calculatre Tr, Rel, Reg, hg, hl
    r_fact = Pr_inf^(1/3);
    Tr = T_inf*(1+(gam_inf-1)/2*M_inf^2*r_fact);
    
    Re_g = Rho_inf*v_inf*(2*r(i))/Visc_inf;
    Re_l = rho_l*v_l*d_tube(i)/visc_l;
    
    hg = k_inf/(2*r(i))*0.026*Re_g^0.8*Pr_inf^0.4;
    
    % Iterate through Tw guess
    delta_q = error+1;
    Twg_guess = 27; % K
    j = 0;
    while (delta_q >=error) && (j <= stop_count)
       j = j+1;
       Twg_guess = Twg_guess + incr;
       
       ql = hg*(Tr - Twg_guess);
       [k_ss] = SS_prop(Twg_guess);
       Twl = Twg_guess - ql*t_tube/k_ss;
       
%        [~,~,~,visc_l_w] = LH2_prop(Twl);
%        hl = k_l/d_tube(i)*0.0214*Re_l^0.8*Pr_l^0.4*(visc_l/visc_l_w)^0.14;
       hl = k_l/d_tube(i)*0.0214*Re_l^0.8*Pr_l^0.4;
        
       qr = hl*(Twl - T_l_prev);
       
       delta_q = abs(qr/ql-1);
    end
    
    % save vectors
    Tr_vect(i) = Tr; 
    Re_l_vect(i) = Re_l;
    Pr_l_vect(i) = Pr_l;
    ql_vect(i) = ql;
    qr_vect(i) = qr;
    
    hl_vect(i) = hl;
    hg_vect(i) = hg;
    Twl_vect(i) = Twl;
    q_dot(i) = qr;
    delta_q_vect(i) = delta_q;
    Twg_vect(i) = Twg_guess;
    % error report
    if (delta_q > error)
        er_rep(i) = 1;
    else
        er_rep(i) = 0;
    end
    
    % Update Pressure and Temperature
    f(i)= fsolve(@(f) fricFunc(f,Re_l),0.1);
    P_tube(i) = P_l_prev - rho_l*g*dz-(0.5*rho_l*v_l^2)*4*f(i)*S(i)/d_tube(i);
    T_tube(i) = T_l_prev + qr*S(i)*Per_tube(i)/(2*cp_l*W_fuel);
    
end
