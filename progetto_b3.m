% Progetto di regolatore per sistema satellite in orbita
% 
% Controlli Automatici T
% Umberto Laghi
% Marco Missiroli
%
% SPECIFICHE
% 
% 1) Errore a regime nullo a un gradino w(t) = 1(t) e d(t) = 1(t) pari a 0.1
%
% 2) Margine di fase >= 55°
%
% 3) S% <= 1%
%
% 4) Ta,5 <= 0.05s
%
% 5) Attenuazione di almeno 35dB per d(t)
%       con [omega_d_min, omega_d_MAX] = [0,0.5]
%
% 6) Attenuazione di almeno 70dB per n(t)
%       con [omega_n_min, omega_n_MAX] = [2.5*10^4,10^5]


clear all; close all; clc

%estremi del grafico
omega_plot_min = 1e-5;
omega_plot_max = 1e7;

%% Parametri

beta1 = 0.25; % [N/m^2] Coefficiente di attrito
beta2 = 2.5; % [N/m^2] Coefficiente di attrito
m = 250; % [kg] Massa
k = 3.5; % Parametro adimensionale
rho_e = 3e8; % [m] Distanza tra terra e satellite
k_g = 6.67e-11; % [N*m^2/kg] Costante di gravitazione universale
M = 5.98e24; % [kg] Massa della Terra

%% Specifiche

% ampiezze gradini
WW = -3e-7;
DD = 2e-8;
NN = 3e-8;

% errore a regime
e_star = 0;

% margine di fase
m_fase = 55;

% Sovraelongazione massima e tempo d'assestamento all'1%
S_100_spec = 0.01;
T_a5_spec = 0.05;

% attenuazione disturbo sull'uscita
A_d = 35;
%omega_d_min = 0; %non lo esplicito
omega_d_MAX = 0.5;

% attenuazione disturbo di misura
A_n = 70;
omega_n_min = 2.5e4;
omega_n_MAX = 1e5;

%Matrici al punto di equilibrio
%x_e = [rho_e, 0, sqrt((k_g*M)/(rho_e^3))];
x1_e = rho_e;
x2_e = 0;
x3_e = sqrt((k_g*M)/(rho_e^3));
u_e =  beta2*sqrt((k_g*M)/rho_e);

A_e = [0,                                                          1,                      0;
     (k-1)*((-3*k_g*M)/(rho_e^3)),                           -beta1/m,        -2*(k-1)*sqrt((k_g*M)/rho_e);
     -(beta2/(m*(rho_e^2)))*sqrt((k_g*M)/rho_e), -(2/rho_e)*sqrt((k_g*M)/(rho_e^3)), -beta2/m];

B_e = [0; 0; 1/(m*rho_e)];

C_e = [0, 0, 1];

D_e = [0];

%% Creazione sistema

% funzione di trasferimento
s = tf('s');
[N,D]=ss2tf(A_e,B_e,C_e,D_e);
GG=tf(N,D);
zpk(GG)
%{
N = (1/(m*rho_e))*(s^2 + s*(beta1/m) + (3*(k-1)*k_g*M)/(rho_e^3));
D = s^3 + (s^2)*((beta1+beta2)/m) + s*((beta1*beta2)/(m^2)-((k-1)*k_g*M)/(rho_e^3)) + ((k-1)*k_g*M*beta2)/(m*(rho_e^3));
GG = N/D
%}
%% Diagramma di Bode

figure(1);
bode(GG,{omega_plot_min,omega_plot_max});
grid on, zoom on;

% return;
%% Regolatore statico - proporzionale con polo nell'origine

% mu statico scelto in maniera aleatoria
mu_s =1;
RR_s = mu_s / s;

% Sistema esteso
GG_e = RR_s*GG;


%% Diagrammi di Bode di Ge con specifiche

figure(2);
hold on;

% Calcolo specifiche S% => Margine di fase
xi = 0.826;
S_100 = 100*exp(-pi*xi/sqrt(1-xi^2))
Mf_spec = xi*100; % Mf_spec = 82.6°

% Disegna diagrammi di Bode di G_e con i margini di stabilità ed esce
% Impostare a 0 per i diagrammi con le specifiche
if 0
    margin(GG_e,{omega_plot_min,omega_plot_max});
    grid on, zoom on;
    return;
end

% Specifiche su d
omega_d_min = omega_plot_min; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -600; -600];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 400; 400];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = omega_plot_min; % lower bound per il plot
omega_Ta_MAX = 300/(Mf_spec*T_a5_spec); % omega_c >= 300/(Mf*T^*)
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -600; -600];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;


% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

% STOP qui per le specifiche
if 0
    return;
end


%% Design del regolatore dinamico

%Ci troviamo nello scenario B, quindi è necessaria una rete anticipatrice
% Rete anticipatrice

Mf_star = Mf_spec;
omega_c_star = 150;
[mag_omega_c_star, arg_omega_c_star, ~] = bode(GG_e, omega_c_star);

mag_omega_c_star_db = 20*log10(mag_omega_c_star);

M_star = 10^(-mag_omega_c_star_db/20);
phi_star = Mf_star - 180 - arg_omega_c_star;

% Formule di inversione
tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha = (cos(phi_star*pi/180) - 1/M_star)/(M_star - cos(phi_star*pi/180));

check_flag = cos(phi_star*pi/180) - 1/M_star
if check_flag < 0
    disp('Errore: alpha negativo');
    return;
end

R_a = (1 + tau*s)/(1 + alpha*tau*s); % rete anticipatrice

% return;

%% Diagrammi di Bode con specifiche includendo regolatore dinamico

%metto un mu_d per far rispettare l'attenuazione di A_n
mu_d = 10^(-0.1); %guadagno -2dB
%mu_d = 1;
R_d = mu_d*R_a;
LL = R_d*GG_e; % funzione di anello

figure(3);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
legend(Legend_arg);

% STOP qui per sistema con controllore dinamico + specifiche
if 0
    return;
end

% Funzione di sensibilità complementare
FF=LL/(1+LL);

margin(FF,{omega_plot_min,omega_plot_max});
margin(LL,{omega_plot_min,omega_plot_max});
%patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+M_f,-180+M_f,-270,-270],'g','FaceAlpha',0.2,'EdgeAlpha',0); 
legend('G_e(s)','Margine di fase', 'F(s)', 'L(s)');
hold off;
%Il fatto che |F(s)| sia = 0       per omega << omega_c 
%                        = |L(s)|  per omega >> omega_c 
%conferma che sia giusto

%% Check prestazioni in anello chiuso

% FUNZIONI DI SENSITIVITA'
% Funzione di sensitività
SS = 1/(1+LL);

% DIAGRAMMI DI BODE
% Rappresento le funzioni di sensitività e L(s) in uno stesso grafico
% in modo da mostrare che rispettano le specifiche
figure(4)
hold on;
bode(LL,{omega_plot_min,omega_plot_max})
bode(FF,{omega_plot_min,omega_plot_max})
bode(SS,{omega_plot_min,omega_plot_max})
legend('L(s)','F(s)', 'S(s)');
grid on, zoom on, hold off;


%% RISPOSTE ANELLO CHIUSO 

%  Disturbo in uscita
omega_d = 0.125;
omega_n = 2.5e4;

tt = (0:1e-2:200)';% 200 secondi con passo 0.01

% segnali in ingresso: riferimento e disturbo di misura
ww = WW*heaviside(tt-100);
dd = DD*(sin(omega_d*tt)+sin(2*omega_d*tt)+sin(3*omega_d*tt)+sin(4*omega_d*tt));
nn = NN*(sin(omega_n*tt)+sin(2*omega_n*tt)+sin(3*omega_n*tt)+sin(4*omega_n*tt));

% Risposte in anello CHIUSO
y_w = lsim(FF,ww,tt);
y_d = lsim(SS,dd,tt);
y_n = lsim(-FF,nn,tt);
y_tot = y_w + y_n + y_d;

%valuto anche l'errore, che deve essere nullo nella situazione statica
e_w = lsim(SS, ww, tt);
e_d = lsim(-SS, dd, tt);
e_n = lsim(FF, nn, tt);
e_tot = e_w + e_n + e_d;
%e_tot = ww - y_tot;
%% RISPOSTE IN ANELLO CHIUSO

figure(5)
hold on, grid on, zoom on
plot(tt,ww,'m')
plot(tt,y_w,'b')
grid on
legend('ww','y_w')

figure(6)
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('d(t)','y_d')

figure(7)
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('n(t)','y_n')

figure(8)
hold on, grid on, zoom on
plot(tt,ww,'m')
plot(tt,y_tot,'b')
grid on
legend('w(t)','y_{tot}(t)')

figure(9)
hold on, grid on, zoom on
plot(tt,e_tot,'b')
grid on
legend('e(t)')

% Risposta al gradino
figure(10);
WW = 3e-7;
T_simulation = 0.15;
[y_step,t_step] = step(WW*FF, 0.15);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[WW*(1+S_100_spec),WW*(1+S_100_spec),WW+(0.5e-7),WW+(0.5e-7)],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
ylim([0,WW+0.5e-7]);

% vincolo tempo di assestamento al 5%
LV = abs(evalfr(WW*FF,1e-7)) % valore limite gradino: W*F(0)
patch([T_a5_spec,T_simulation,T_simulation,T_a5_spec],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5_spec,T_simulation,T_simulation,T_a5_spec],[LV*(1+0.05),LV*(1+0.05),LV+1e-7,LV+1e-7],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%step(FF);
%grid on, zoom on;
%I requisiti di sovraelongazione e tempo di assestamento al 5% sono soddisfatti
%Il margine di fase è rispettato poichè nella nuova frequenza di taglio
%omega_c non si è nella zona proibita

%Ampiezza L in omega_n <= -70 -> requisito disturbo n rispettato
mag2db(abs(evalfr(LL,omega_n)))

%Ampiezza L in omega_d >= 35 -> requisito disturbo n rispettato
mag2db(abs(evalfr(LL,omega_d)))

%% PUNTO 5
%{
devo plottare la traiettoria del sistema non linearizzato, considerando
anche d(t) e n(t)
uso ode45 
%}

%{

% Risposta al gradino
figure(4);


[y_step,t_step] = step(WW*FF, 1);
plot(t_step,y_step,'b');



Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

% Diagramma di bode sistema in anello chiuso
figure(5);
bode(FF, {omega_plot_min,omega_plot_max})
grid on, zoom on

if 0
    return;
end

%}

%{
%% ROBA VECCHIA

% Funzione di sensitività
%FF ce l'ho già
SS = 1/(1+LL);
%figure(6);



%Faccio le trasformate di Laplace dei segnali e trovo la trasformata di y(t)
W = laplace(ww);
D = laplale(dd);
N = laplace(nn);

Y = FF*(W - N) + SS*D;
%Antitrasformo e la confronto con omega(t)
y = ilaplace(Y);

figure(9);
plot(tt, y);
hold on, grid on, zoom on;
plot(tt, ww);
hold off;
%}