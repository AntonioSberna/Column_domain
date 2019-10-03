

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                        %%%   
%%%    Dominio di resistenza per sezioni presso-inflesse   %%%
%%%         (valido solo per sezioni rettangolari)         %%%
%%%                                                        %%%
%%%          antoniopio.sberna@studenti.polito.it          %%%
%%%                                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear 
clear all

%
% Geometrie

% Altezza della sezione [mm]
h = 3700;
% Larghezza trave [mm]
b = 360;
%Ricoprimento copriferro [mm]
d_prime = 40;
%Area armatura inferiore [mm^2]
A_s = 3308.1;
%Area armatura superiore [mm^2]
A_sprime = 3308.1;


%
% Materiali

% Resistenza a compressione del calcestruzzo [MPa]
f_ck = 25;
% Tensione di snervamento acciaio [GPa]
f_yk = 450;


%
% Azioni

% Forza normale agente [kN]
N_Ed = 1355.3;
% Momento agente [kN*m]
M_Ed = 3938.6;



% Passo di analisi
passo = 0.1;



%%%
%%% Fine input
%%%

%%% Parametri acciaio
% Modulo di Young acciaio [GPa]
E_s = 200;
% Deformazione allo snervamento acciaio [%]
epsilon_yd = 1.96;
% Deformazione a rottura acciaio 
epsilon_ud = 67.5;

% Resistenza a compressione del cls di calcolo
f_cd = 0.85 * f_ck / 1.5;
% Tensione di snervamento acciaio di progetto [GPa]
f_yd = f_yk /1.15;


%altezza utile
d = h - d_prime;
%deformazione fibra estrema cls
epsilon_c = 3.5;
% Deformazione acciaio per epsilon_cls nullo
epsilon_s_cls0 = - epsilon_c * d_prime / h;

% Definizione variabile 
temp = 1;

% Dimensione matrice
dim = ceil((-epsilon_s_cls0 + epsilon_ud)/passo + 1);

% Definizione matrici output
%sigma = zeros(700,1);
%sigma2 = zeros(700,1);
%eps = zeros(700,1);
Rd = zeros(dim,4);


% Campo 5
    %epsilon_c = epsilon_s = epsilon_sprime = 2 0/00    

    % Force of respectively upper and lower steel
    S_prime = f_yd * A_sprime / 1000;
    S = f_yd * A_s / 1000;
    
    % Force of concrete
    C =  f_cd * b*h /1000;
    
    % Forza resistente
    N_rd = C + S_prime + S;    
    
    % Momento resistente
    M_rd = ((h /2 - d_prime) * (S_prime - S))/1000;
    
    
    Rd(temp,1) = N_rd;
    Rd(temp,2) = M_rd; % C'è qualcosa che non va
    Rd(temp,3) = -2;

    
    temp = temp+1;



% Campo 3, 4 e 4a
for epsilon_s = epsilon_s_cls0:passo:epsilon_ud
    x_u = (epsilon_c /(epsilon_s + epsilon_c)) * d;
    epsilon_sprime = (epsilon_c /x_u) * (d_prime + x_u);
    
    % Tension on lower steel
    if epsilon_s < epsilon_yd
        if epsilon_s > - epsilon_yd
            sigma_s = E_s * epsilon_s;
        else
            sigma_s = -f_yd;
        end
    else
        sigma_s = f_yd;
    end
    
    % Tension on upper steel
    if epsilon_sprime < epsilon_yd
        if epsilon_sprime > - epsilon_yd
            sigma_sprime = E_s * epsilon_sprime;
        else
            sigma_sprime = -f_yd;
        end
    else
        sigma_sprime = f_yd;
    end
    
    % Position of upper steel force
    if d_prime < x_u
        i = -1;
    else
        i = 1;
    end
   
    % Force of respectively upper and lower steel
    S_prime = i * sigma_sprime * A_sprime / 1000;
    S = sigma_s * A_s / 1000;
    
    % Calcolo coefficienti di riempimento
    [beta_1,beta_2] = coeff_riemp(epsilon_c);
    
    % Position of concrete force
    y_c = beta_2 * x_u;
  
    
    % Force of concrete
    C = -1 * beta_1 * f_cd * x_u * b / 1000;
    
    % Forza resistente
    N_rd = C + S_prime + S;    
    
    % Momento resistente
    M_rd = (-(h/2 - y_c) * C + (h /2 - d_prime) * (S - S_prime))/1000;
    

    %x_uc(temp,1) = x_u;
    %eps(temp,1) = epsilon_sprime;
    %sigma(temp,1) = sigma_s;
    %sigma2(temp,1) = sigma_sprime;
    
    Rd(temp,1) = -N_rd;
    Rd(temp,2) = M_rd;
    
    Rd(temp,3) = epsilon_s;
    Rd(temp,4) = x_u/d;
    
    temp = temp+1;
    
end



      %Limite dominio per domini con armatura non simmetrica
      lim(1,1) = Rd(1,1);
      lim(1,2) = Rd(1,2);
      lim(2,1) = Rd(dim,1);
      lim(2,2) = Rd(dim,2);
      
      
        
% Plot del dominio di resistenza
%plot(N_Ed,M_Ed,'+')
plot(Rd(:,1),Rd(:,2))
xlabel('Sforzo normale N_{rd} [kN]')
ylabel('Momento  M_{rd} [kN\cdot m]')
hold on
plot(N_Ed,M_Ed,'+')
if A_s ~= A_sprime
plot(lim(:,1),lim(:,2), '--')
end
hold off

%Determinazione coefficienti di riempimento
function [beta1, beta2] = coeff_riemp(epsilon)
    beta = [
        0.1 0.0492 0.3347;
        0.2 0.0967 0.3362;
        0.3 0.1425 0.3377;
        0.4 0.1867 0.3393;
        0.5 0.2292 0.3409;
        0.6 0.2700 0.3426;
        0.7 0.3092 0.3443;
        0.8 0.3467 0.3462;
        0.9 0.3825 0.3480;
        1.0 0.4167 0.3500;
        1.1 0.4492 0.3520;
        1.2 0.4800 0.3542;
        1.3 0.5092 0.3564;
        1.4 0.5367 0.3587;
        1.5 0.5625 0.3611;
        1.6 0.5867 0.3636;
        1.7 0.6092 0.3663;
        1.8 0.6300 0.3690;
        1.9 0.6492 0.3720;
        2.0 0.6667 0.3750;
        2.1 0.6825 0.3782;
        2.2 0.6970 0.3814;
        2.3 0.7101 0.3846;
        2.4 0.7222 0.3878;
        2.5 0.7333 0.3909;
        2.6 0.7436 0.3939;
        2.7 0.7531 0.3968;
        2.8 0.7619 0.3996;
        2.9 0.7701 0.4022;
        3.0 0.7778 0.4048;
        3.1 0.7849 0.4072;
        3.2 0.7917 0.4095;
        3.3 0.7980 0.4118;
        3.4 0.8039 0.4139;
        3.5 0.8095 0.4160];
    
    if epsilon == 3.5
        beta1 = 0.8095;
        beta2 = 0.4160;
    elseif epsilon < 3.5 && epsilon > 0
        i = 1;
        
        while epsilon > beta(i,1)
            i = i + 1;
        end
    
        beta1 = ((epsilon - beta(i+1,1))/(beta(i,1) - beta(i+1,1))) * beta(i,2) - ((epsilon - beta(i,1))/(beta(i,1) - beta(i+1,1))) * beta(i+1,2);
        beta2 = ((epsilon - beta(i+1,1))/(beta(i,1) - beta(i+1,1))) * beta(i,3) - ((epsilon - beta(i,1))/(beta(i,1) - beta(i+1,1))) * beta(i+1,3);
    end 
    
end

