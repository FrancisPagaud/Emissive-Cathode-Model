
%% Code for the article Pagaud, Dolique, Plihon 2023, PSST

% Author : Francis Pagaud, 20/09/2022

%% Parameters
global interr Ih r sigma Ieth_only epsilon_0 x_max V_b t_test dx dt alpha beta Te n n0 Vp debut_step fin_step width gamma power V

global u_t_past_2 t_test_2
u_t_past_2 = 0;
t_test_2 = 0;
Ieth_only = 0;

% Global parameters
r = 2.53e-4;
sigma = 5.67e-8;
epsilon_0 = 0.26;
alpha = 0.15;           %Fraction of radiation to the neighbours. Le calcul donne environ 9.8 %
beta = 1.013;           %Ion bombardment correction factor
I = 16.2;               % EXPERIMENTAL PARAMETERS
Vb = -62;
n0 = 1e18;              %Initial density
n = n0;
Te = 4.5;               %Electron temperature
t_test = 0;
Vp = -1;                %Initial plasma potential
V = Vp;
gamma = 6.6;         %Vp fit = -1 - 6.6*Ib_sim^0.335
                     %Vp(t) = Vp - gamma*Ieth^(1/power)
power = 2.985;
debut = 5;  %time to reach stationary T profile before plasma (s)
fin = 12;   %End of plasma
t_max = fin+1; % End of simulation (time for temperature relaxation)

% Domaines spatiaux et temporels
x_max = 198e-3 ; % total length of filament

dx = 1e-3;          %Spatial step
dt = 1e-3;          %Time step
width = 0.01/dx; %Transition length for self heating via dissipative radiations (10 mm)
debut_step = round(debut/dt);       %Step for plasma start
fin_step = round(fin/dt);           %Step for plasma end

x = (0:dx:x_max) ; %Spatial domain
t = (0:dt:t_max) ; %Temporal domain

plot_results = 1;
save_results = 0;


e = 1.602e-19;
M = 40*1.67e-27;        %Argon mass
eps_T = 75;             %Energy for ion-electron recombination
R = 5e-2;               %Radius of the plasma column


for Ih = I          %%% SCANNING ALL EXPERIMENTAL PARAMETERS
    for V_b = Vb
%         close all
        interr = 1;         %Boolean. interr == 1 : plasma stopped. 
                            %It is the case when runaway behaviour or
                            %before debut_step / after fin_step
        saving_name = strcat("C:\Users\Francis\Documents\These\Figures\matlab\Pyrometre\Simu_res\Ih_", num2str(Ih), "_Vb_", num2str(V_b), "_Vp_", num2str(Vp), "_beta_", num2str(beta), "_gamma_", num2str(gamma));

        %%% EQUATION : 
        %%% c(x,t,u,du/dx) * du/dt = x^−m * d/dx [x^m f(x,t,u,du/dx)] + s(x,t,u,du/dx)

        % m parameter
        m = 0 ;

        % This line sets the time step for the calculus, it can fasten the
        % numerical solving.
        options = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'Stats', 'on', 'NonNegative',1, 'MaxStep', 5e-3) ;
        sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t, options);   %%% SOLVING THE SYSTEM DEFINED BELOW

        % Extract the temperature solution as u(time, space).
        u = sol(:,:,1);

        pas = 10;   %Lighter solution for manipulating results: dt*10
        len = (fin_step-debut_step)/pas+1;
        U_h = zeros(len, length(x));        %Bias along the cathode
        Ieth_cum = zeros(len, length(x));   %Sum of I_b currents from cathode to plasma between s = 0 and x
        Ieth_tot = zeros(1,len);            %Current I_b from cathode to the plasma over time
        Vp_step = zeros(1,len);             %Plasma potential over time
        N = zeros(1,len);                   %Plasma density over time
        
        %%% Initial state
        Ieth_tot(1) = trapz(Ieth(u(debut_step,:),dx)) + Iis(n0,Te,x_max);
        Ieth_cum(1,:) = cumtrapz((0:x_max/dx),Ieth(u(debut_step,:),dx)) + Iis(n0,Te,0:dx:x_max);
        U_h(1,:) = V_b + cumtrapz((0:x_max/dx),rho(u(debut_step,:))*dx/(pi*r^2).*Ih);
        N(1) = n0;
        Vp_step(1) = Vp;
        
        %%% Over time
        for i = 2:len
            U_h(i,:) = V_b + cumtrapz((0:x_max/dx),rho(u(debut_step+pas*(i-1),:))*dx/(pi*r^2).*(Ih+Ieth_tot(i-1)-Ieth_cum(i-1,:)));
            %%% Ieth_s is the space-charged limited solution according to
            %%% Ye and Takamura 2000, PoP
            Ieth_tot(i) = trapz( min(Ieth_s(Vp_step(i-1), U_h(i-1,:), Te, N(i-1), dx), Ieth(u(debut_step+pas*(i-1),:),dx)) ) + Iis(N(i-1),Te,x_max) ;
            Ieth_x = min(Ieth_s(Vp_step(i-1), U_h(i-1,:), Te, N(i-1), dx), Ieth(u(debut_step+pas*(i-1),:),dx));     %Thermionic emission along cathode at step i. Necessary for power balance
            Ieth_cum(i,:) = cumtrapz((0:x_max/dx),min(Ieth_s(Vp_step(i-1), U_h(i-1,:), Te, N(i-1), dx), Ieth(u(debut_step+pas*(i-1),:),dx))) + Iis(N(i-1),Te,0:dx:x_max);
            Vp_step(i) = Vp - gamma*Ieth_tot(i)^(1/power);
            N(i) = n0 + 3*sum(Ieth_x.*(Vp_step(i)-U_h(i,:)))/(pi*R^2*e*eps_T*sqrt(e*Te/M));    % Solving power balance
        end
        Ieth_max = max(Ieth_tot);

        if save_results == 1
            t_save = t(debut_step:10:end);
            u_save = u(debut_step:10:end,:);
            save(strcat(saving_name, '.mat'), 'I_h', 'V_b', 'Vp', 'alpha', 'beta', 'gamma', 'power', 'Te', 'x', 't_save', 'debut', 'fin', 'u_save', 'N', 'Vp_step', 'Ieth_tot', 'Ieth_cum', 'U_h', 'pas', 'interr', 'Ieth_max')
        end

        %% Plot figure: T(s,t), U_h(s,t), I_b(t)

        if plot_results == 1

%             close all
            k_max = 6;
            color = flipud(hot(2*k_max)) ;

            fh = figure(40);
            subplot("Position", [0.08, 0.15, 0.25, 0.78])
            hold on

            for kk = 1:k_max
                plot(x,u(debut_step+1+round(length(t(debut_step+2:fin_step))*kk/k_max) ,:), 'LineWidth', 4, 'Color', color(k_max + kk,:), 'DisplayName', strcat('t =', num2str(round(length(t(debut_step+1:fin_step))*dt*kk/k_max,1)), ' s'))
            end
            plot(x,u(debut_step,:), '--', 'LineWidth', 4, 'Color', color(k_max,:), 'DisplayName', 'Initial value')

            legend('Location', 'northeast', 'box', 'off', 'NumColumns', 2)

            xlabel('s (m)', 'FontSize', 28)
            ylabel('T (K)', 'FontSize', 28)
            ax = gca();
            ax.FontSize = 20;
            ax.LabelFontSizeMultiplier = 1.55;
            box('on')


            subplot("Position", [0.4, 0.15, 0.25, 0.35])
            patch([t(debut_step:pas:fin_step)-debut+0.001 nan],[Ieth_tot(1:len) nan],[t(debut_step:pas:fin_step)-debut+0.001 nan],[t(debut_step:pas:fin_step)-debut+0.001 nan], 'edgecolor', 'interp', 'Linewidth', 5); 
            colormap(color(k_max:k_max*2,:))
            % plot(t,Ieth_tot, '-', 'LineWidth', 5, 'Color', color(k_max,:))
            hold on
            scatter(t(1),Ieth_tot(1), 250, color(k_max,:), 'filled')
            for kk = 1:k_max
                scatter(t(debut_step+1+round(length(t(debut_step+1:fin_step))*kk/k_max))-debut,Ieth_tot(1+round(length(t(debut_step+1:pas:fin_step))*kk/k_max)), 250, color(k_max + kk,:), 'filled')
            end

            % legend('Location', 'best', 'Interpreter', 'latex')

            xlabel('t (s)', 'FontSize', 28)
            ylabel('I_{em} (A)', 'FontSize', 28)
            ax = gca();
            ax.FontSize = 20;
            ax.LabelFontSizeMultiplier = 1.55;
            box('on')


            subplot("Position", [0.4, 0.58, 0.25, 0.35])
            patch([t(debut_step:pas:fin_step)-debut+0.001 nan],[U_h(1:len,end)'-V_b nan],[t(debut_step:pas:fin_step)-debut+0.001 nan],[t(debut_step:pas:fin_step)-debut+0.001 nan], 'edgecolor', 'interp', 'Linewidth', 5); 
            colormap(color(k_max:k_max*2,:))
            % plot(t,Ieth_tot, '-', 'LineWidth', 5, 'Color', color(k_max,:))
            hold on
            scatter(t(1),U_h(1,end)-V_b, 250, color(k_max,:), 'filled')
            for kk = 1:k_max
                scatter(t(debut_step+1+round(length(t(debut_step+1:fin_step))*kk/k_max))-debut,U_h(1+round(length(t(debut_step+1:pas:fin_step))*kk/k_max),end)-V_b, 250, color(k_max + kk,:), 'filled')
            end

            % legend('Location', 'best', 'Interpreter', 'latex')

            ylabel('U_h (V)', 'FontSize', 28)
            ax = gca();
            ax.XTickLabel = {};
            ax.FontSize = 20;
            ax.LabelFontSizeMultiplier = 1.55;
            box('on')



            subplot("Position", [0.72, 0.15, 0.25, 0.78])
            plot([x(1) x(end)], [Vp Vp], '-.', 'linewidth', 3, 'Color', color(k_max,:))
            hold on
            plot(x,U_h(1,:), '--', 'LineWidth', 4, 'Color', color(k_max,:))


            for kk = 1:k_max
                plot(x,U_h(1+round(length(t(debut_step+1:10:fin_step))*kk/k_max-1),:), 'LineWidth', 4, 'Color', color(k_max + kk,:))
                plot([x(1) x(end)], [Vp_step(1+round(length(t(debut_step+1:10:fin_step))*kk/k_max-1)) Vp_step(1+round(length(t(debut_step+1:10:fin_step))*kk/k_max-1))], '-.', 'linewidth', 3, 'Color', color(k_max + kk,:))
            end


            xlabel('s (m)', 'FontSize', 28)
            ylabel('U_h (V)', 'FontSize', 28)
            ax = gca();
            ax.FontSize = 20;
            ax.LabelFontSizeMultiplier = 1.55;
            box('on')
            pause(0.5)
          
            fh.WindowState = 'maximized';

            if save_results == 1
                saveas(gcf, strcat(saving_name, '.png'))
            end
            

        end
    end
end
%% Fonctions en parametres de pdepe

% --------------------------------------------------------------
% EQUATION DIFFERENTIELLE


function [c,f,s] = pdex1pde(x, t,u,DuDx)
global Ih r sigma u_inter interr u_t_past u_t_past_2 Ieth_tot Ieth_only Ieth_cum V_b U_h t_test t_test_2 dx dt alpha beta x_max Te n n0 Vp V debut_step fin_step width gamma power
e = 1.602e-19;
E_i = 15.78;        %Ionisation energy Ar+
W = 4.54;           %Work fct of tungsten for ion bombardment
k_b = 1.38e-23;
eps_T = 75;
M = 40*1.67e-27;
R = 5e-2;

%%% x is a scalar between 0 and x_max. ind is the index it is referred to
%%% in the spatial domain
len = length(x-dx/2:dx:x_max);
ind = round((x_max+dx)/dx)+1-len;
      
calcul_n = 0;   %Parameter to compute plasma density. Is conducted only when 
                %we get to a new time step
if (ind == 1) && (t_test ~= round(t/dt))   %Going back to the filament start
    %%% Here is a heavy technical loop: the way Matlab solve it is
    %%% optimized to save time. But with a step function such as plasma
    %%% ignition, discrepancies prevents the convergence of the algorithm,
    %%% so the solver has to take a step back in time. If it is the case,
    %%% the integral computed at time step i-1 becomes wrong. The following
    %%% "if" loops enable to correct that effect in case of 1 or 2 steps
    %%% back in a row of the PDE solver. 
    
    %%% It could easily be improved.
    
    disp('New loop')
    disp(num2str(round(t/dt)))  %Time index
    calcul_n = 1;   %boolean. calcul_n = 1 : increase in plasma density, cf line 232
    u_t_past_3 = u_t_past_2;
    u_t_past_2 = u_t_past;
    u_t_past = u_inter;
    
    
    
    if t_test < round(t/dt)     %We advanced in time
        Ieth_only = trapz( min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_inter,dx)) );  %Calculation of integrals using temperature at timestep i-1
        Ieth_x = min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_inter,dx));
        Ieth_tot = Ieth_only + Iis(n, Te, x_max);
        Ieth_cum = cumtrapz((0:x_max/dx), min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_inter,dx)) ) + Iis(n, Te, (0:dx:x_max));
        U_h = V_b + cumtrapz((0:x_max/dx),rho(u_inter)*dx/(pi*r^2).*(Ih+Ieth_tot-Ieth_cum));
        t_test_2 = 0;
    elseif t_test_2 == 0        %A step back occurred
        disp('Step back #1')
        Ieth_only = trapz( min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_t_past_2,dx)) ); %Calculation of integrals using temperature at timestep i-2
        Ieth_x = min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_t_past_2,dx));
        Ieth_tot = Ieth_only + Iis(n, Te, x_max);
        Ieth_cum = cumtrapz((0:x_max/dx), min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_t_past_2,dx)) ) + Iis(n, Te, (0:dx:x_max));
        U_h = V_b + cumtrapz((0:x_max/dx),rho(u_t_past_2)*dx/(pi*r^2).*(Ih+Ieth_tot-Ieth_cum));
        t_test_2 = 1;
    else                        %Two steps back in a row
        disp('STEP BACK #2')
        Ieth_only = trapz( min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_t_past_3,dx)) ); %Calculation of integrals using temperature at timestep i-3
        Ieth_x = min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_t_past_3,dx));
        Ieth_tot = Ieth_only + Iis(n, Te, x_max);
        Ieth_cum = cumtrapz((0:x_max/dx), min(Ieth_s(V, U_h, Te, n, dx), Ieth(u_t_past_3,dx)) ) + Iis(n, Te, (0:dx:x_max));
        U_h = V_b + cumtrapz((0:x_max/dx),rho(u_t_past_3)*dx/(pi*r^2).*(Ih+Ieth_tot-Ieth_cum));
        t_test_2 = 0;
    end
    t_test = round(t/dt);
end



% EQUATION
% c(x,t,u,du/dx) * du/dt = x^−m * d/dx [x^m f(x,t,u,du/dx)] + s(x,t,u,du/dx)

%%% Plasma ignition
plasma_on = 1;
if round(t/dt) < debut_step || round(t/dt) > fin_step || interr ~= 1
    Ieth_tot = 0;
    Ieth_cum = zeros(1,x_max/dx+1);
    plasma_on = 0;
end


%%% In case of runaway behaviour: disruption to avoid numerical issues
if round(t/dt) > debut_step && round(t/dt) < fin_step && interr == 1 && Ieth_tot > 25
    plasma_on = 0;
    interr = t/dt-20;
end


%%% Plasma density increase
if plasma_on == 1 && calcul_n == 1
    V = Vp - plasma_on*gamma*Ieth_tot^(1/power);
    n =n0 + 3*sum(Ieth_x.*(V-U_h))/(pi*R^2*e*eps_T*sqrt(e*Te/M));      
end
    

%%% Saving temperature at time step i for the integral at time step i+1
u_inter(ind) = u;


%%% TIME-DERIVATIVE TERM: thermal inertia
c = Cpm(u)*pi*r^2;


%%% SELF-HEATING duE TO INwARDSZ DISSIPATIVE RADIATIONS
%%% Boundaries
start_h = 0.01/dx;
end_h = 0.188/dx;

%%% Taking into account the increase in length of the loops over s
reach_start = 0.015/dx;     %First loop : 15 mm
reach_fin = 0.04/dx;        %Last loop : 40 mm
reach_l = reach_interp(ind, start_h+reach_start, end_h, reach_start, reach_fin);    %Define the position of the previous turn
reach_r = reach_interp(ind, start_h, end_h-reach_fin, reach_start, reach_fin);      %Define the position of the next turn

if ind > start_h-width+reach_l && ind < end_h+width
    alpha_l = 1;                                        %In the middle of the spiral: presence of a neighbouring filament before the actual position "ind"
    if ind < start_h+reach_l
        alpha_l = (ind+width-start_h-reach_l)/width;    %Increase progressively between 0 and 1 if in the "width" region
    elseif ind > end_h
        alpha_l = (end_h+width-ind)/width;              %Decrease progressively between 1 and 0 if in the "width" region
    end
    heat_left = sigma*alpha_l*alpha*epsilon(u)*epsilon(u_t_past(ind-reach_l))*2*pi*r*u_t_past(ind-reach_l)^4;   %Absorption of heat from previous loop
else
    heat_left = 0;      %No neighbour before
end

if ind > start_h-width && ind < end_h+width-reach_r 
    alpha_r = 1;                                         %In the middle of the spiral: presence of a neighbouring filament after the actual position "ind"
    if ind < start_h
        alpha_r = (ind+width-start_h)/width;            %Increase progressively between 0 and 1 if in the "width" region
    elseif ind > end_h-reach_r
        alpha_r = (end_h-reach_r+width-ind)/width;      %Decrease progressively between 1 and 0 if in the "width" region
    end
    heat_right = sigma*alpha_r*alpha*epsilon(u)*epsilon(u_t_past(ind+reach_l))*2*pi*r*u_t_past(ind+reach_r)^4;  %Absorption of heat from next loop
else
    heat_right = 0;     %No neighbour after
end

%%% Self-heating
self_heat = heat_right + heat_left;

%%% Thermionic cooling
Ieth_cool = plasma_on*min(Ieth_s(V, U_h(ind), Te, n, dx), Ieth(u,dx))*(W+2*k_b/e*u)/dx;

%%% Dissipative radiation
rad_dissip = sigma*epsilon(u)*2*pi*r*(u^4);

%%% Ohmic heating
DC_heat = rho(u)/(pi*r^2).*(Ih+Ieth_tot-Ieth_cum(ind)).^2;

%%% Ion bombardment
ion_bombardment = plasma_on*Iis(n,Te,1)*(beta*max((V-U_h(ind)),0) + (E_i-W)*max((V-U_h(ind)),0)/(V-U_h(ind)));


%%% SOURCE TERMS
s = self_heat - rad_dissip + ion_bombardment + DC_heat - Ieth_cool;


%%% CONDUCTION TERM
f = lambd(u)*DuDx*pi*r^2;
    
end


% --------------------------------------------------------------
% INITIAL CONDITION
% Raw power balance between homogeneous dissipative radiations and ohmic 
% heating before plasma

function u0 = pdex1ic(~)
global Ih r sigma epsilon_0 x_max u_inter u_t_past Ieth_tot Ieth_cum U_h V_b dx
T0 = (2000:1:3200);

y0_ext = rho(T0)*x_max/(pi*r^2)*Ih^2-sigma*epsilon_0*x_max*2*pi*r*T0.^4;    
index = find( abs(y0_ext) == min(abs(y0_ext)) );    %Powers balance at 0
T_0_ext = T0(index);
u0 = T_0_ext;   %Homogeneous initial condition

%%% Raw initial state
u_inter = u0*ones(1,x_max/dx+1);
u_t_past = u_inter;
Ieth_tot = trapz(Ieth(u_inter,dx));
Ieth_cum = cumtrapz((0:x_max/dx),Ieth(u_inter,dx));
U_h = V_b + cumtrapz((0:x_max/dx),rho(u_inter)*dx/(pi*r^2).*(Ih+Ieth_tot-Ieth_cum));
end

% --------------------------------------------------------------
% BOUNDARY CONDITIONS

% BC are of the form p + q * f(du/dx) = 0 
% l and r correspond to BCs at x_min for l (left) and x_max for r (right)
function [pl,ql,pr,qr] = pdex1bc(~,ul,~,ur,~)

lambda_Cu = 40;     %Effective conductivity of copper considering the great thermal contact resistance
r_Cu = 2e-3;        %Copper rods radius
L_Cu = 0.35;        %Copper rods length
T_Cu = 300;         %Ambient temperature at the other end of copper rods
pl = -(ul-T_Cu)/L_Cu*pi*r_Cu^2*lambda_Cu;
ql = 1;
pr = (ur-T_Cu)/L_Cu*pi*r_Cu^2*lambda_Cu;
qr = 1;
end