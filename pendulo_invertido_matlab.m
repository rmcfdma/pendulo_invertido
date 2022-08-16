%% 1 - Preparação do Script.
close all;
clear all;
clc;

%% 2 - Parâmetros de simulação.
tamanho_legenda = 14;   % Tamanho da fonte da legendo dos gráficos 
tamanho_titulo  = 14;   % Tamanho dos titulos dos gráficos 
espessura_linha = 1;    % Espessura das linhas dos gráficos
limites_grafico = [-1.5 1.5 -0.8 0.8];  % Limites eixo x
s = 1;                  % escala do pêndulo

%% 3 - Definição da Variáveis.
syms x1 x2 x3 x4 u l M m g J dt
x = [x1 x2 x3 x4];      % Variáveis de estados           

%% 4 - Especificações Quanser IP02 (Long).
l = 0.641;          % Comprimento da Haste
m = 0.230;          % Massa do Pêndulo
M = 0.57;           % Massa do Carrinho
g = 9.81;           % Aceleração da Gravidade
J =  m*(l/2)^2;     % Momento de Inércia de uma Barra Delgada do
l_carrinho = 0.10;  % Largura do Carrinho
h_carrinho = 0.15;  % Altura do Carrinho
w_carrinho = 0.061; % comprimento
% As especificações a seguir são fornecidas em resposta a um ponto
% de ajuste de posição do carrinho de onda quadrada de ± 30 mm.
% Tr <= 1.5s % Tempo de subida menor ou igual a 1.5 segundos
% mod(alpha) <= 1 degree % alpha menor ou igual a 1 grau
% V < 10v % Tensão menor que 10 volts
% Carrinho +- 30mm (0.03m ou 3cm)
% Rack dimensions (L x W x H) 102 x 15 x 6.1 cm

%% 5 - Parâmetros construtivos do motor do Pêndulo.
etag = 0.90;         % Planetary geabox efficiency
Kg = 3.71;           % Planetary gearbox gear ratio
Kt = 7.68*10^(-3);   % Motor current-torque constant (N*(m/A)) 
Rm =  2.6;           % Motor armature resistance (ohms)
rmp = 6.35*10^(-3);  % Motor pinion radius (m) 
Km =  7.68*10^(-3);  % Motor back-emf constant (V/(rad/s))
Vm = 6;              % Tensão Nominal do Motor (V)
etam = 0.69;         % Motor efficiency
omega_max = 628.3;   % Maximum motor speed (rad/s)
Jeq = 3.90*10^(-7)   % Momento de Inercia Armadura kg*m2

% v_carrinho_max = rmp*omega_max; % Velocidade Máxima do Carrinho
% u_max_1 = ((etag*Kg*Kt)/(Rm*rmp))*(-((Kg*Km*v_carrinho_max)/rmp) + etam*Vm); % Força máxima entregue pelo motor (Quanser)
% u_max_2 = ((Kg*Km*Vm)/(Rm*rmp))-((Kg^2*Km^2*-v_carrinho_max)/(Rm*rmp^2));    % Força máxima entregue pelo motor

%% 6 - Condições Iniciais.
ci  = [0 0 pi/180 0];    % Condições iniciais do sistema
ci2 = [0 0 pi/180 0];    % Condições iniciais do estimador de estados ekf
ci3 = [0 0 pi/180 0];    % Condições iniciais do estimador de kf

%% 7 - Equações de Estados Não-Lineares.
f = [x(2);
    (u+m*l*x(4)^2*sin(x(3))-m*g*sin(x(3))*cos(x(3)))/(M+m-m*cos(x(3))^2);
     x(4);  
    (u*cos(x(3))-(M+m)*g*sin(x(3))+m*l*x(4)^2*sin(x(3))*cos(x(3)))/(m*l*cos(x(3))^2-(M+m)*l)]; % Vetor de estados não-lineares
h = [x(1);x(3)];  % Vetor de saída

%% 8 - Handles Não-Lineares.
f_h = matlabFunction(f); % Transforma f simbólico em um vetor de funções handles.
h_h = matlabFunction(h); % Transforma h simbólico em um vetor de funções handles.

f = f_h;
h = h_h;
%% 9 - Derivada (Jacobiano - Linearização) das Equações de Estado.
% Jfx = simplify(jacobian(f,x));       % Jacobiano de f em relação a x
% Jfu = simplify(jacobian(f,u));       % Jacobiano de f em relação a u
% Jhx = simplify(jacobian(h,x));       % Jacobiano de h em relação a x
% Jhu = simplify(jacobian(h,u));       % Jacobiano de f em relação a u

%% 10 - Jacobianos para as Matrizes de transição do EKF.
f_j = simplify(jacobian(x' + f*dt,x)); % Derivada da espressão do Método de Euler em realação aos estados 
g_j = simplify(jacobian(x' + f*dt,u)); % Derivada da espressão do Método de Euler em realação ao controle
h_j = simplify(jacobian(h,x));         % Derivada da saída em realação ao controle

%% 11 - Linearização.
x0 = [0 0 0 0]; % Ponto de equilíbrio (estados)
u0 = 0;         % Ponto de equilíbrio (controle)

%% 12 - Linearização no ponto x0 e u0.
% Jfx_lin = subs(Jfx,{x1,x2,x3,x4},x0);  % x = x0 em Jfx
% Jfu_lin = subs(Jfu,{u},u0);            % u = u0 em Jfu
% Jhx_lin = subs(Jhx,{x1,x2,x3,x4},x0);  % x = x0 em Jhx
% Jhu_lin = subs(Jhu,{u},u0);            % u = u0 em Jhu

%% 13 - Matrizes numéricas A = Jfx_lin, B = Jfu_lin, C = Jgx_lin e D = Jgu_lin 
% A = Jfx_lin; % Matriz de estados numérica
% B = Jfu_lin; % Matriz de Entrada numérica
% C = Jhx_lin; % Matriz de Saída numérica
% D = Jhu_lin; % Matriz de transmissão direta numérica

%% 14 - Matrizes simbólicas das Equações de Estado Linearizadas. (Basta calcular os jacobianos para M,m,g e l simbólicos)
A = [ 0 1 0 0; 0 0 -(g*m)/M 0; 0 0 0 1; 0 0 -(g*(M + m))/(-l*M) 0]; % Matriz A simbólica
B = [0; 1/M; 0; -(1/(l*M))]; % Matriz B simbólica
C = [1 0 0 0; 0 0 1 0];      % Matriz C simbólica
D = [0;0];                   % Matriz D simbólica

%% 15 - Matrizes de Transição - Jacobianos de (x + f(x,u)*dt) em relação a x e u, e de h em relação a x.
% F_j = @(x,u,dt)[1 dt 0 0;
%     0 1 -(dt*(14743*x(4)^2*cos(x(3))-451260*cos(x(3))^2+225630))/(23000*cos(x(3))^2-80000)-(4600*dt*cos(x(3))*sin(x(3))*((14743*sin(x(3))*x(4)^2)/100000+u-(22563*cos(x(3))*sin(x(3)))/10000))/(23*cos(x(3))^2-80)^2 (14743*dt*x(4)*sin(x(3)))/(500*(23*sin(x(3))^2 + 57));
%     0 0 1 dt;
%     0 0 (4600000*dt*cos(x(3))*sin(x(3))*((14743*cos(x(3))*sin(x(3))*x(4)^2)/100000-(981*sin(x(3)))/125+u*cos(x(3))))/(641*(23*cos(x(3))^2-80)^2)-(dt*(784800*cos(x(3))-29486*x(4)^2*cos(x(3))^2+100000*u*sin(x(3))+14743*x(4)^2))/(14743*cos(x(3))^2-51280) (14743*dt*x(4)*sin(2*x(3)))/(100000*((14743*cos(2*x(3)))/200000 - 87817/200000)) + 1]; % Matriz de transição de estados não-linear na forma de função handle reduzida      
% G_j = @(x,dt)[0;-(100*dt)/(23*cos(x(3))^2 -
% 80);0;(100000*dt*cos(x(3)))/(14743*cos(x(3))^2 - 51280)]; % Matriz de
% transição de entrada não-linear na forma de função handle reduzida
% H_j = [1 0 0 0; 0 0 1 0]; Matriz de transição de saída linear na forma de
% função handle reduzida

F_j = matlabFunction(f_j); % Matriz de transição de estados não-linear na forma de função handle criada pelo Matlab
G_j = matlabFunction(g_j); % Matriz de transição de entrada não-linear na forma de função handle criada pelo Matlab
H_j = matlabFunction(h_j); % Matriz de transição de saída linear na forma de função handle criada pelo Matlab

F = F_j;                   % Trocando o nome da função
H = H_j;                   % Trocando o nome da função
%% 16 - Variáveis Temporais.
tempo_simulacao = 10;      % Tempo de simulação

%% 17 - Matriz de transição para o sistema linearizado.
% Pelo ZOH. Mapeia corretamente do plano s para o plano z.
% F = expm(A*dt);          % Matriz de transição de estados linearizada
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B;  % Matriz de entrada de estados linearizada

% Pelo método de Euler. Pode mapear do plano s para o plano z incorretamente.
F2 = eye(4,4)+dt*A;        % Matriz de transição de estados linearizada
% G = dt*B                 % Matriz de transição de entrada linearizada

%% 18 - Teste de Controlabilidade.
cc = ctrb(A,B);               % Obtém a matriz de controlabilidade
%cc = [B A*B (A^2)*B (A^3)*B] % n-1 = 3, n = 4
[linhas,~] = size(cc);        % Obtém as dimensões de cc
posto = rank(cc);             % Obtém o posto de ob. Como o posto = 4 = n, o sistema é controlável

%% 19 - Teste de Observabilidade.
ob = obsv(A,C);                   % Obtém a matriz de observabilidade
%ob = [C; C*A; C*(A^2); C*(A^3)]  % n-1 = 3, n = 4
[linhas,colunas] = size(ob);      % Obtém as dimensões de ob
posto = rank(ob);                 % Obtém o posto de ob. Como o posto = 4 = n, o sistema é observável

%% 20 - Projeto do LQR.
Q = [  100 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];       % Matriz de ponderação dos estados
R = 1;                 % Matriz de poderação do controle
K = lqr( A, B, Q, R);  % Matriz de Ganho de realimentação
Q_lqr = Q;             % Troca de variável
R_lqr = R;             % Troca de variável

%% 21 - Tentativa frustrada de encontrar o período de amostragem pelo controle digital.
% syms s z T
% Sistema em malha aberta
% sys_ss = ss(A,B,C,D);       % Sistema em malha aberta representado no espaço de estados
% [num,den] = ss2tf(A,B,C,D); % Sistema mimo 1x2 em malha aberta representado em função de transferência
% sys1 = tf(num(1,:),den);    % Função de transferência para a posição do carrinho
% sys2 = tf(num(2,:),den);    % Função de transferência para o ângulo da haste
% sys = [sys1;sys2];          % Vetor coluna com os dois sistemas     
% r = roots(den);             % Obtém as raizes do denominador
%  
% % Espansão em frações parciais no domínio de Laplace
% fp =  vpa(-0.0544/(s+4.6346) + 0.0544/(s-4.6346) + 1.0880e-16/s + 1.25/(s^2),4);
% 
% % Tranformada Z de fp
% fpz = vpa((-0.0544*z)/(z-exp(-4.6346*T)) +(0.0544*z)/(z+exp(4.6346*T)) + (1.0880e-16*z)/(z-1) + (1.25*T*z)/((z-1)^2));
% 
% % Sistema em malha aberta no domínio Z
% [numz,denz] = numden(fpz);
% numz = collect(numz);
% denz = collect(denz);
% Gz = vpa(numz/denz,4);
% % Numerador de (1 + KG(z))
% den_tz = (1 + K(1)*Gz);
% [num,den] = numden(den_tz);
% numTz = vpa(collect(num),4);    % Coloca em odem decrescente de potência
% a = 4.6346;
% % T = 0.336181;
% den_tz = vpa(z^4*exp(a*T) - z^3*(2*exp(a*T) - 1.544*exp(2*a*T) + 12.5*T*exp(a*T) + 0.456) + z^2*(12.5*T - 3.088*exp(2*a*T) - 12.5*T*exp(2*a*T) + 0.912) + z*(2*exp(a*T) + 1.544*exp(2*a*T) + 12.5*T*exp(a*T) - 0.456) - 1.0*exp(a*T),4);
% % den_tz = sym2poly(den_tz);
% raizes = roots(den_tz);
% vetor = [exp(a*T)  -(2*exp(a*T) - 1.544*exp(2*a*T) + 12.5*T*exp(a*T) + 0.456) (12.5*T - 3.088*exp(2*a*T) - 12.5*T*exp(2*a*T) + 0.912) (2*exp(a*T) + 1.544*exp(2*a*T) + 12.5*T*exp(a*T) - 0.456) -exp(a*T)]
% dd = subs(vetor)
% abs(vpa(roots(dd),4))
% r = solve(den_tz,z, 'MaxDegree', 4)

%% 22 - Sistema em malha aberta.
sysmass = ss(A,B,C,D);      % Sistema mimo 1x2 em malha aberta representadono espaço de estados
[num,den] = ss2tf(A,B,C,D); % Sistema mimo 1x2 em malha aberta representado em função de transferência
sys1 = tf(num(1,:),den);    % Função de transferência para a posição do carrinho
sys2 = tf(num(2,:),den);    % Função de transferência para o ângulo da haste
sys = [sys1;sys2];          % Vetor coluna com os dois sistemas

%% 23 - Sistema em malha fechada.
sysmfss = ss(A-B*K,B,C,D);            % Espaço de estados da função de transferência dos sistema em malha fechada
[nummf,denmf]= ss2tf(A-B*K,B,C,D);    % Passando do espaço de estados para função de transferência dos sistema em malha fechada
sysmftf1 = tf(nummf(1,:),denmf);      % Função de transferência para a posição do carrinho
sysmftf2 = tf(nummf(2,:),denmf);      % Função de transferência para o ângulo da haste
sysmftf = [sysmftf1;sysmftf2];        % Vetor coluna com os dois sistemas

%% 24 - Análise das margens de estabilidade no domínio da frequência
[mag1,phase1,wout1] = bode(sysmftf);  % Diagrama de Bode para o sistema 1

%% 25 - Polos em Malha Aberta
% p_ma = pole(sysmass);            % Polos do sistema em malha aberta
%diagrama_pz1 = pzplot(sysmass);   % Diagrama de polos e zeros
% bode(sysmass);                   % Diagrama de Bode 
% a = findobj(gca,'type','line');  % Encontrando os objetos do tipo linha
% for i = 1:length(a)
%     set(a(i),'markersize',18);   % Muda o tamanho da marcação do polo no gráfico
%     set(a(i), 'linewidth',3);    % Muda a espessura da linha
% end
% ylabel('Eixo imaginário','FontSize',18);           % Label do eixo y
% xlabel('Eixo Real','FontSize',18);                 % Label do eixo x
% title('Diagrama de Polos e Zeros','FontSize',18);  % Titulo do diagrama de polos e zeros
% grid on;                                           % Grade

%% 26 - Polos em malha fechada
% diagrama_pz1 = pzplot(sysmftf);  % Diagrama de polos e zeros para Q11 = 100
% polos = eig(A-B*K);              % Polos do sistema em malha fechada
% frequences = abs(polos);         % Frequência dos polos do sistema em MF
% a = findobj(gca,'type','line');  % Encontrando os objetos do tipo linha
% for i = 1:length(a)
%     set(a(i),'markersize',18);   % Muda o tamanho da marcação do polo no gráfico
%     set(a(i), 'linewidth',3);    % Muda a espessura da linha
% end
% ylabel('Eixo imaginário','FontSize',18);           % Label do eixo y
% xlabel('Eixo Real','FontSize',18);                 % Label do eixo x
% title('Diagrama de Polos e Zeros','FontSize',18);  % Titulo do diagrama de polos e zeros
% grid on;                                           % Grade

% bode(sysmfss);                                     % Diagrama de Bode (somente módulo) dos sistemas 1 e 2
% ylabel('Módulo','FontSize',18);                    % Label do eixo y
% xlabel('Frequência','FontSize',18);                % Label do eixo x
% title('Diagrama de Bode','FontSize',18);           % Titulo do diagrama Bode (somente módulo)
% grid on;                                           % Grade

%% 27 - Cáculo da Largura de Banda.
omega_bw = max(abs(eig(sysmfss)));   % Valor da largura de banda pelo comando do matlab

% opt1 = bodeoptions;                % Muda as opções do diagrama de Bode para o sistema 1        
% opt1.Grid = 'on';                  % Habilita a grade
% opt1.FreqScale = 'log';            % Apresenta a grade com escalalogarítmica
% opt1.Title.String = '';            % Apaga o título do diagrama de Bode
% opt1.Title.FontSize = 14;          % Fote do título caso houvesse
% opt1.XLabel.String = 'Frequência'; % Texto do eixo x
% opt1.XLabel.FontSize = 14;         % Fonte do eixo x
% opt1.YLabel.String = {'Magnitude','Fase'}; % Texto do eixo y
% opt1.YLabel.FontSize = 14;         % Fonte do eixo y
% bodeplot(sysmftf(1),opt1);         % Valor da largura de banda por inspeção do diagrama de bode do sistema em malha fechada: 3.67 rad/s

% opt2 = bodeoptions;                % Muda as opções do diagrama de Bode para o sistema 2
% opt2.Grid = 'on';                  % Habilita a grade
% opt2.FreqScale = 'log';            % Apresenta a grade com escalalogarítmica
% opt2.Title.String = '';            % Apaga o título do diagrama de Bode
% opt2.Title.FontSize = 14;          % Fote do título caso houvesse
% opt2.XLabel.String = 'Frequência'; % Texto do eixo x
% opt2.XLabel.FontSize = 14;         % Fonte do eixo x
% opt2.YLabel.String = {'Magnitude','Fase'}; % Texto do eixo y
% opt2.YLabel.FontSize = 14;         % Fonte do eixo y
% bodeplot(sysmftf(2),opt2);         % Valor da largura de banda por inspeção do diagrama de bode do sistema 2 que de infinito


%% 28 - Cálculo do Período de Amostragem pelo Cálculo Numerico, Teorema de Nyquist e valores usuais.
% Todos os exemplos encontrados de estabilidade numérica não envolviam 
% matrizes, logo os resultados para a estabilidade marginal aqui 
% apresentados foram deduzidos intuitivamente, utilizando o módulo dos
% auto-valores de (A-BK), que forneceram o resultado esperado quando 
% alterado qualquer valor da matriz Q do LQR, com exceção do elemento 
% Q(3,3), onde dt_marginal não tornou o algorítmo marginalmente estável,
% conforme esperado, por um motivo ainda desconhecido.

% Valores para a estabilidade marginal
I = eye(4);                          % Cria uma matriz identidade 4x4
M = diag(eig(A-B*K));                % Cria uma matriz diagonal com os autovalores de A-BK
S = 2*I*vpa(inv(abs(M)),6);          % Faz 2*I dividido por uma matriz diagonal que contém o módulo dos autovalores 
S_norma_2 = 2/(norm(M,2));           % Faz 2*I dividido por uma matriz diagonal que contém o módulo dos autovalores 
dt_marginal_1 = S(1,1);              % Pega o elemento (1,1) de S que forneceu o valor de dt para a estabilidade marginal
dt_marginal_2 = S(2,2);              % Pega o elemento (2,2) de S que forneceu o valor de dt para a estabilidade marginal
dt_marginal_3 = S(3,3);              % Pega o elemento (3,3) de S que forneceu o valor de dt para a estabilidade marginal
dt_marginal_4 = S(4,4);              % Pega o elemento (4,4) de S que forneceu o valor de dt para a estabilidade marginal
dt_marginal_norma = S_norma_2;       % dt pela norma 2 de M

% Limite segundo o Teorema de Amostragem de Nyquist
omega_nyquist = 2*omega_bw;          % Frequência de Nyquist
dt_nyquist = (2*pi)/(omega_nyquist); % Máximo período de amostragem segundo o Teorema de Nyquist

% Valor ideal de 20 a 40 omega_bw (Franklin. pg. 504)
omega_s_min = 20*omega_bw;           % Valor mínimo usual da frequência de amostragem
omega_s_max = 40*omega_bw;           % Valor máximo usual da frequência de amostragem
omega_s_escolhido = 35*omega_bw;     % Valor escolhido para a frequancia de amostragem dentro da faixa usual
dt_max = (2*pi)/(omega_s_min);       % Período de amostragem máximo usual
dt_min = (2*pi)/(omega_s_max);       % Período de amostragem mínimo usual
dt_escolhido = round((2*pi)/(omega_s_escolhido),2); % Príodo de amostragem escolhido
dt = 0.02;                           % Período de amostragem escolhido para as simulações

%% 29 - Definindo os parâmetros probabilísticos.
         
% Variância
sigma_quadrado_w = 0.01;             % Covariância do ruído do sistema
sigma_quadrado_v = 0.01;             % Covariância do ruído dos sensores

% Desvio Padrão
sigma_w = sqrt(sigma_quadrado_w);    % Desvio padrão do ruído do sistema
sigma_v = sqrt(sigma_quadrado_v);    % Desvio padrão do ruído do sistema

bias = 0;                            % Acurácia
media = 0;                           % Média

% Matriz de covariância        
P =  eye(4,4);                       % Criando a matriz de covariância entre o sistema e as estimativas  
P2 = eye(4,4);                       % Criando a matriz de covariância entre o sistema e as estimativas  
Px = diag(P);                        % Criando o vetor acumulador de diagonal principal das matrizes covariância
Px2 = diag(P2);                      % Criando o vetor acumulador de diagonal principal das matrizes covariância

Q_w = sigma_quadrado_w * eye(4,4)*0.01; % Matriz de incertezas do processo
Q_v = sigma_quadrado_v * eye(2,2)/0.01; % Matriz de incertezas da medição

%% 30 - Variáveis utilizadas no Método Iterativo.
x = ci';             % Condições iniciais
X = x;               % Inicializando o acumulador de estados
y = [x(1);x(3)];     % Inicializando o a saída
Y = y;               % Inicializando o acumulador de saídas
x2 = ci';            % Condições iniciais do estimador
X2 = x2;             % Inicializando o acumulador de condições iniciais estimadas (Não utilizado)
y2 = [x2(1);x2(3)];  % Inicializando a medição
Y2 = y2;             % Inicializando acumulador de saídas estimadas (Não utilizado)
x3 = ci';            % Condições iniciais
X3 = x3;             % Inicializando o acumulador de estados
y3 = [x3(1);x3(3)];  % Inicializando o acumulador de medições
Y3 = y3;             % Inicializando o acumulador de saída
x4 = ci';            % Condições iniciais
X4 = x4;             % Inicializando o acumulador de estados
y4 = [x4(1);x4(3)];  % Inicializando a medição
Y4 = y4;             % Inicializando o acumulador de medições
u3 = -K*x3;          % Inicializando o sinal de controle
U3 = u3;             % Inicializando o acumulador de sinal de controle
u4 = -K*x4;          % Inicializando o sinal de controle
U4 = u4;             % Inicializando o acumulador de sinal de controle
xhat = ci2';         % Estado estimado
Xhat = xhat;         % Inicializando o acumulador de estados estimados
yhat = [xhat(1);xhat(3)];     % Saída estimada 
Yhat = yhat;         % Inicializando o acumulador de saídas estimadas
xhat2 = ci3';        % Estado estimado
Xhat2 = xhat2;       % Inicializando o acumulador de estados estimados
yhat2 = [xhat2(1);xhat2(3)];  % Saída estimada 
Yhat2 = yhat2;       % Inicializando o acumulador de saídas estimadas
yt = y - yhat;       % Erro de estimativa inicial
Yt = yt;             % Acumulando os erro de estimativa
yt2 = y2 - yhat2;    % Erro de estimativa inicial
Yt2 = yt2;           % Acumulando os erro de estimativa
u = -K*x;            % Sinal de controle inicial
U = u;               % Acumulando o sinal de controle 
u2 = -K*x2;          % Sisnal de controle inicial do sistema 2
U2 = u2;             % Acumulando o sinal de controle do sistema 2
Q = Q_w;             % Matriz de incertezas do processo
R = Q_v;             % Matriz de incertezas das medições
l_k = P*C'*inv(C*P*C'+R);    % Covariância do erro inicial
L_k = l_k;                   % Inicializando o acumulador de covariância do erro
l_k2 = P2*C'*inv(C*P2*C'+R); % Covariância do erro inicial
L_k2 = l_k2;         % Inicializando o acumulador de covariância do erro
gravar = 0;          % Variável que indica se a função animar_pendulo deve ou não gravar a animação
t = 0;               % Tempo inicial
H = [];              % Covariância predita
w =  media + sigma_w * randn(4,1); % Criando o ruído no sistema inicial
v =  media + sigma_v * randn(2,1); % Criando o ruído no no sensor inicial

F = F_j;            % Matriz de Transição de estados para o sistema não-linear
f = f_h;            % Vetor de estados não-lineares
F2 = eye(4) + A*dt; % Matriz de Transição de estados para o sistema linearizado

%% 31 - Simulação do sistema não-linear sem LQR (Utilizando comandos do Matlab).
% t0 = 0;                        % Tempo inicial
% tf = tempo_simulacao;          % Tempo final
% % K = [0 0 0 0];               % Controle cancelado
% opts = odeset('RelTol',1e-5);  % Erro relativo
% [t,X] = ode45(@equacoes_nao_lineares,[t0 tf],ci,opts,K); % Ou -> ode45(@(t,x)equacoes_nao_lineares_lqr(t,x,K),[to tf],ci,opts);
% U = -K*X';                     % Sinal de controle
% Y = X(:,[1 3]);                % Seleciona os estados de saída
% 
% plotar_sistema('um_com_lqr',t',X',Y',U',[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y,[],[],[],[],[],s,l,l_carrinho,h_carrinho,'ode45','','','','','','ode45_lqr','Animação do modelo não-linear.',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,0); % Animação

%% 32 - Simulação do sistema linearizado sem LQR (Utilizando comandos do Matlab).
% tf = 10;                         % Tempo final
% % K = [0 0 0 0];                 % Controle cancelado
% sys = ss(A-B*K,B,C,D);           % Sistema em malha fechada linearizado
% [Y,t,X] = initial(sys,ci,tf);    % Simula o sistema linearizado
% U = -K*X';                       % Sinal de controle

% plotar_sistema('um_com_lqr',t',X',Y',U',[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y,[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Initial','','','','','','initial_lqr','Animação do modelo linearizado.',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,0); % Animação

%% 33 - Simulação dos sistemas linear e não linear com LQR (Utilizando comandos do Matlab).
dt = 0;
t0 = 0;                            % Tempo inicial
tempo_simulacao = 10;              % Tempo de simulação
tf = tempo_simulacao;              % Tempo final
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);                % Tolerância relativa e absoluta
[t1,x1] = ode45(@equacoes_nao_lineares,[t0 tf],ci,opts,K); % Simulando o sistema não linear
Y = x1(:,[1 3]);                   % Estados medidos não-linear
u1 = -K*x1';                       % Sinal de controle do sistem não-linear

sys = ss(A-B*K,B,C,D);             % Sistema em malha fechada linearizado
% [y2,t2,x2] = initial(sys,ci,tf); % Simula o sistema linearizado
Y2 = [];T2 = [];X2 = [];           % Acumuladores
for i = 1:2:length(t1)             % Loop para pegar os elementos de t1 de 3 em 3
    if i ~= length(t1)             % se i ainda não chegou no final
      [y2,t2,x2] = initial(sys,ci,[t1(i,1),t1(i+1,1)]);  % simula o sistema e recolhe de dois em dois pontos de t1
      Y2 = [Y2 y2'];               % Acumula y2
      X2 = [X2 x2'];               % Acumula x2
      T2 = [T2 t2'];               % Acumula t2
    end
end

u2 = -K*X2;;                                  % Sinal de controle do sistem linearizado    
limites_grafico = [-0.03 0.03 -0.025 0.025];  % Limites para o gráfico de 1º
s = 0.03;                                     % Escala para 1º 
% plotar_sistema('dois_sem_filtragem_com_lqr',t1',x1',x1(:,[1 3])',u1',T2',X2',Y2',u2',[],[],[],[],[],[],[],[],[],[],espessura_linha,'Não-Linear','Linearizado','',tamanho_legenda,tamanho_titulo)
animar_pendulo(2,Y,Y2',[],[],[],[],s,l,l_carrinho,h_carrinho,'ode45','Initial','','','','','initial_ode45_initial_lqr','Animação dos modelos não-linear e linearizado.',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,1); % Animação

%% 34 - Simulação do sistema para os diferentes valores calculados para o período de amostragem T.
% dt = dt_nyquist;                % Período de amostragem calculado pelo Teorema de Nyquist
% dt = dt_max;                    % Período de amostragem máximo usual sugerido pela bibliografia
% dt = dt_min;                    % Período de amostragem mínimo usual sugerido pela bibliografia
% dt = dt_marginal_1;             % T(1,1) calculado pelo método do auto-valor para a estabilidade numérica do Método de Euler
% dt = dt_marginal_2;             % T(2,2) calculado pelo método do auto-valor para a estabilidade numérica do Método de Euler
% dt = dt_marginal_3;             % T(3,3) calculado pelo método do auto-valor para a estabilidade numérica do Método de Euler
% dt = dt_marginal_4;             % T(4,4) calculado pelo método do auto-valor para a estabilidade numérica do Método de Euler
% dt = dt_marginal_norma;         % T calculado pela norma 2
% dt = dt_animacao;               % Período de amostragem escolhido para as animações
% dt = dt_simulacao;              % Período de amostragem escolhido para as simulações
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                   % Lei de controle 
%     x = x + (A*x + B*u)*dt;     % Discretização diferencial pelo método de Euler Direto
%     y = C*x;                    % Medição
%     X = [X x];                  % Acumulando o estado
%     Y = [Y y];                  % Acumulando a saída
%     U = [U u];                  % Acumulando o sinal de controle
%     t = [t i];                  % Acumulando o tempo
% end

% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Medido','','','','','','linearizado_recursivo_com_lqr_dt','Animação do modelo não-linear (recursivo)',Q_lqr,R_lqr,0,0,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 35 - Sistema não-linear (medido estocástico e esperança) sem controle e (com e sem) ruído.
% for i = 0:dt:tempo_simulacao
%       u = 0;                                       % Lei de controle da esperança do sistema 
%       u2 =0;                                       % Lei de controle do sistema real 
%       x = x + (f(u,x(2),x(3),x(4)) + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w
%       x2 = x2 + f(u2,x2(2),x2(3),x2(4))*dt;        % Discretização diferencial pelo método de Euler Direto
%       y = h(x(1),x(3));                            % Adicionando o ruído na medição
%       y2 = h(x2(1),x2(3));                         % Adicionando o ruído na medição
%       X = [X x];                                   % Acumulando o estado
%       X2 = [X2 x2];                                % Acumulando o estado
%       Y = [Y y];                                   % Acumulando a saída
%       Y2 = [Y2 y2];                                % Acumulando a saída
%       U = [U u];                                   % Sinal de controle da esperança do sistema 
%       U2 = [U2 u2];                                % Sinal de controle do sistema real
%       t = [t i];                                   % Acumulando o tempo
%       w =  media + sigma_w * randn(4,1);           % Criando o ruído no sistema
%       v =  media + sigma_v * randn(2,1);           % Criando o ruído no no sensor
% end

% plotar_sistema('dois_sem_filtragem_sem_lqr',t,X,Y,U,t,X2,Y2,U2,[],[],[],[],[],[],[],[],[],[],espessura_linha,'Medição','Esperança','',tamanho_legenda,tamanho_titulo)
% animar_pendulo(2,Y',Y2',[],[],[],[],s,l,l_carrinho,h_carrinho,'Estocástico','Esperança','','','','','nao_linear_recursivo_medido_esperanca_com_lqr','Animação do modelo não-linear (recursivo)',Q_lqr,R_lqr,0,0,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 36 - Sistema não-linear sem controle e sem ruído.
% for i = 0:dt:tempo_simulacao
%       u = 0;                             % Lei de controle 
%       x = x + f(u,x(2),x(3),x(4))*dt;    % Discretização diferencial pelo método de Euler Direto
%       y = h(x(1),x(3));                  % Medição
%       X = [X x];                         % Acumulando o estado
%       Y = [Y y];                         % Acumulando a saída
%       U = [U u];                         % Sinal de controle
%       t = [t i];                         % Acumulando o tempo
% end

% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Medido','','','','','','nao_linear_recursivo_sem_lqr','Animação do modelo não-linear (recursivo).',zeros(4),0,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 37 - Sistema não-linear com controle LQR e sem ruídos.
% for i = 0:dt:100 
%     u = -K*x
%     x = x + (f(u,x(2),x(3),x(4)))*dt;    % Discretização diferencial pelo método de Euler Direto das Equações de estados não-lineares   
%     y = C*x;   
%     X = [X x];                           % Acumulando o estado
%     Y = [Y y];                           % Acumulando a saída
%     U = [U u];                           % Acumulando o sinal de controle 
%     t = [t i];                           % Acumulando o Tempo
% end

% % plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Não-Linear','','','','','','nao_linear_recursivo_com_lqr','Animação do modelo não-linear (recursivo)',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 38 - Sistema não-linear com controle LQR e com ruido apenas no sistema.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                               % Lei de controle
%     x = x + (f(u,x(2),x(3),x(4)) + w)*dt;   % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = h(x(1),x(3));                       % Medição
%     X = [X x];                              % Acumulando o estado
%     Y = [Y y];                              % Acumulando a saída
%     U = [U u];                              % Acumulando o sinal de controle
%     t = [t i];                              % Acumulando o Tempo
%     w =  media + sigma_w * randn(4,1);      % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1);      % Criando um novo ruído no no sensor     
% end
% 
% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Medido','','','','','','nao_linear_w_recursivo_com_lqr','Animação do modelo não-linear (Recursivo).',Q_lqr,R_lqr,0,0,sigma_quadrado_w,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 39 - Sistema não-linear com controle LQR e com ruido no sistema e no sensor.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                                  % Lei de controle
%     x = x + (f(u,x(2),x(3),x(4)) + w)*dt;      % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = h(x(1),x(3)) + v;                      % Medição com ruído v
%     X = [X x];                                 % Acumulando o estado
%     Y = [Y y];                                 % Acumulando a saída
%     U = [U u];                                 % Acumulando o sinal de controle
%     t = [t i];                                 % Acumulando o tempo
%     w =  media + sigma_w * randn(1,1);         % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1);         % Criando um novo ruído no no sensor
% end
 
% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Medido','','','','','','nao_linear_w_v_recursivo_com_lqr','Animação do modelo não-linear (Recursivo).',Q_lqr,R_lqr,0,0,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação


%% 40 - Sistema linearizado sem controle e sem ruido.  
% for i = 0:dt:tempo_simulacao
%     u = 0;                      % Lei de controle  
%     x = x + (A*x + B*u)*dt;     % Discretização diferencial pelo método de Euler Direto
%     y = C*x;                    % Medição
%     X = [X x];                  % Acumulando o estado
%     Y = [Y y];                  % Acumulando a saída
%     U = [U u];                  % Acumulando o sinal de controle
%     t = [t i];                  % Acumulando o tempo
% end

% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Medido','','','','','','linearizado_recursivo_sem_lqr','Animação do modelo linearizado (Recursivo).',zeros(4),0,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 41 - Sistema linearizado com controle LQR e sem ruido. 
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                   % Lei de controle 
%     x = x + (A*x + B*u)*dt;     % Discretização diferencial pelo método de Euler Direto
%     y = C*x;                    % Medição
%     X = [X x];                  % Acumulando o estado
%     Y = [Y y];                  % Acumulando a saída
%     U = [U u];                  % Acumulando o sinal de controle
%     t = [t i];                  % Acumulando o tempo
% end

% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Medido','','','','','','linearizado_recursivo_lqr','Animação do modelo linearizado (Recursivo).',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 42 - Simulação dos sistemas linearizados pelo algoritmo recursivo e via rotinas do Matlab em conjunto.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                    % Lei de controle 
%     x = x + (A*x + B*u)*dt;      % Discretização diferencial pelo método de Euler Direto
%     y = C*x;                     % Medição
%     X = [X x];                   % Acumulando o estado
%     Y = [Y y];                   % Acumulando a saída
%     U = [U u];                   % Acumulando o sinal de controle
%     t = [t i];                   % Acumulando o tempo
% end

% tf = tempo_simulacao;            % Tempo final
% sys = ss(A-B*K,B,C,D);           % Sistema em malha fechada linearizado
% [Y2,t2,X2] = initial(sys,ci,tf); % Simula o sistema linearizado
% U2 = -K*X2';                     % Sinal de controle

% plotar_sistema('dois_sem_filtragem_com_lqr',t,X,Y,U,t2',X2',Y2',U2',[],[],[],[],[],[],[],[],[],[],espessura_linha,'Recursivo','Initial','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(2,Y',Y2,[],[],[],[],s,l,l_carrinho,h_carrinho,'Linearizado (Recursivo)','Linearizado (initial)','','','','','linearizado_recursivo_initial_lqr','Animação do modelo linearizado.',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 43 - Simulação dos sistemas não-lineares pelo algoritmo recursivo e via rotinas do Matlab em conjunto.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                            % Lei de controle
%     x = x + f(u,x(2),x(3),x(4))*dt;      % Discretização diferencial pelo método de Euler Direto
%     y = h(x(1),x(3));                    % Medição do sistema linearizado - recursivo
%     X = [X x];                           % Acumulando o estado
%     Y = [Y y];                           % Acumulando a saída
%     U = [U u];                           % Acumulando o sinal de controle
%     t = [t i];                           % Acumulando o tempo
% end

% t0 = 0;                        % Tempo inicial
% tf = tempo_simulacao;          % Tempo final
% opts = odeset('RelTol',1e-5);  % Erro relativo
% [t2,X2] = ode45(@equacoes_nao_lineares,[t0 tf],ci,opts,K); % Ou -> ode45(@(t,x)equacoes_nao_lineares_lqr(t,x,K),[to tf],ci,opts);
% U2 = -K*X2';                   % Sinal de controle
% Y2 = X2(:,[1 3]);              % Medição do sistema não-linear - ode45
 
% plotar_sistema('dois_sem_filtragem_com_lqr',t,X,Y,U,t2',X2',Y2',U2',[],[],[],[],[],[],[],[],[],[],espessura_linha,'Recursivo','ode45','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(2,Y',Y2,[],[],[],[],s,l,l_carrinho,h_carrinho,'Não-Linear (Recursivo)','Não-Linear (ode45)','','','','','nao_linear_recursivo_ode45_lqr','Animação do modelo não-linear.',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 44 - Sistema linearizado com controle LQR e com ruido aditivo no sistema.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                          % Lei de controle 
%     x = x + (A*x + B*u + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w 
%     y = C*x;                           % Medição
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída
%     U = [U u];                         % Acumulando o sinal de controle
%     t = [t i];                         % Acumulando o tempo
%     w =  media + sigma_w * randn(1,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end

% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Modelo Linearizado (Recursivo)','','','','','','linearizado_recursivo_w_lqr','Animação do modelo linearizado (recursivo)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 45 - Sistema linearizado com controle LQR e com ruido no sistema e no sensor.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                          % Lei de controle
%     x = x + (A*x + B*u + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = C*x + v;                       % Medição com ruído v
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída
%     U = [U u];                         % Acumulando o sinal de controle
%     t = [t i];                         % Acumulando o tempo
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end

% plotar_sistema('um_com_lqr',t,X,Y,U,[],[],[],[],[],[],[],[],[],[],[],[],[],[],espessura_linha,'','','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% (1,Y',[],[],[],[],[],s,l,l_carrinho,h_carrinho,'Modelo Linearizado (Recursivo)','','','','','','linearizado_recursivo_w_v_lqr','Animação do modelo linearizado (recursivo)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 46 - Sistemas Não-Linear e Linearizado com LQR e sem ruídos simulados em conjunto.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                               % Lei de controle - Linearizado
%     u2 = -K*x2;                             % Lei de controle - Não-Linear
%     x = x + (A*x + B*u)*dt;                 % Discretização diferencial pelo método de Euler Direto - Linearizado
%     x2 = x2 + f(u2,x2(2),x2(3),x2(4))*dt;   % % Discretização diferencial pelo método de Euler Direto - Não-Linear
%     y = C*x;                                % Medição com ruído v - Linearizado
%     y2 = h(x2(1),x2(3));                    % Medição com ruído v - Não-Linear
%     X = [X x];                              % Acumulando o estado - Linearizado
%     X2 = [X2 x2];                           % Acumulando o estado - Não-Linear
%     Y = [Y y];                              % Acumulando a saída - Linearizado
%     Y2 = [Y2 y2];                           % Acumulando a saída - Não-Linear
%     U = [U u];                              % Acumulando o sinal de controle - Linearizado
%     U2 = [U2 u2];                           % Acumulando o sinal de controle - Não-Linear
%     t = [t i];                              % Acumulando o tempo     
% end

% limites_grafico = [-0.03 0.03 -0.025 0.025];
% s = 0.03;

% plotar_sistema('dois_sem_filtragem_com_lqr',t,X,Y,U,t,X2,Y2,U2,[],[],[],[],[],[],[],[],[],[],espessura_linha,'Linearizado','Não-Linear','',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(2,Y',Y2',[],[],[],[],s,l,l_carrinho,h_carrinho,'Linearizado (Recursivo)','Não-Linear (Recursivo)','','','','','linearizado_nao_linear_recursivo_lqr','Animação dos modelos linearizado e não-linear (Recursivo).',Q_lqr,R_lqr,0,0,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 47 - KF - Filtro de Kalman seguindo o sistema linearizado sem ruido e a esperança.
% for i = 0:dt:tempo_simulacao
%     u = 0;                             % Lei de controle do sistema estocástico - KF
%     u2 = 0 ;                           % Lei de controle da esperança
%     x = x + (A*x + B*u + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w  - KF
%     x2 = x2 + (A*x2 + B*u2)*dt;        % Discretização diferencial pelo método de Euler Direto - Esperança 
%     y = C*x;                           % Medição
%     y2 = C*x2;                         % Medição da esperança
%     % Predição
%     xhat  = xhat + (A*xhat + B*u)*dt;  % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F2*P*F2'+ Q;                  % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     X2 = [X2 x2];                      % Acumulando o estado
%     Y2 = [Y2 y2];                      % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     U2 = [U2 u2];                      % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Linearizado','Esperança Linearizada','KF','','','','kf_sem_ruido','Animação do modelo linearizado (Filtro de Kalman)',zeros(4),0,Q_w,Q_v,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 48 - KF - Filtro de Kalman seguindo o sistema linearizado com ruido no sistema.
% for i = 0:dt:tempo_simulacao
%     u = 0;                             % Lei de controle do sistema estocástico - KF
%     u2 = 0 ;                           % Lei de controle da esperança
%     x = x + (A*x + B*u + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w  - KF
%     x2 = x2 + (A*x2 + B*u2)*dt;        % Discretização diferencial pelo método de Euler Direto - Esperança 
%     y = C*x;                           % Medição
%     y2 = C*x2;                         % Medição da esperança
%     % Predição
%     xhat  = xhat + (A*xhat + B*u)*dt;  % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F2*P*F2'+ Q;                  % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     X2 = [X2 x2];                      % Acumulando o estado
%     Y2 = [Y2 y2];                      % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     U2 = [U2 u2];                      % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Linearizado','Esperança Linearizada','KF','','','','kf_w','Animação do modelo linearizado (Filtro de Kalman)',zeros(4),0,Q_w,Q_v,sigma_quadrado_w,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 49 - KF - Filtro de Kalman seguindo o sistema linearizado com ruido no sistema e no sensor.
% for i = 0:dt:tempo_simulacao
%     u = 0;                             % Lei de controle do sistema estocástico - KF
%     u2 = 0 ;                           % Lei de controle da esperança
%     x = x + (A*x + B*u + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w  - KF
%     x2 = x2 + (A*x2 + B*u2)*dt;        % Discretização diferencial pelo método de Euler Direto - Esperança 
%     y = C*x + v;                       % Medição com ruído v
%     y2 = C*x2;                         % Medição da esperança
%     % Predição
%     xhat  = xhat + (A*xhat + B*u)*dt;  % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F2*P*F2'+ Q;                  % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     X2 = [X2 x2];                      % Acumulando o estado
%     Y2 = [Y2 y2];                      % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     U2 = [U2 u2];                      % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Linearizado','Esperança Linearizada','KF','','','','kf_w_v','Animação do modelo linearizado (Filtro de Kalman)',zeros(4),0,Q_w,Q_v,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 50 - LQG (Filtro de Kalman + LQR) sem ruído.
% for i = 0:dt:tempo_simulacao
%     u = -K*xhat;                       % Lei de controle do sistema estocástico - KF
%     u2 = -K*x2;                        % Lei de controle da esperança
%     x = x + (A*x + B*u)*dt;            % Discretização diferencial pelo método de Euler Direto com ruído w  - KF
%     x2 = x2 + (A*x2 + B*u2)*dt;        % Discretização diferencial pelo método de Euler Direto - Esperança 
%     y = C*x;                           % Medição com ruído v
%     y2 = C*x2;                         % Medição da esperança
%     % Predição
%     xhat  = xhat + (A*xhat + B*u)*dt;  % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F2*P*F2'+ Q;                  % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     X2 = [X2 x2];                      % Acumulando o estado
%     Y2 = [Y2 y2];                      % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     U2 = [U2 u2];                      % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Linearizado','Esperança Linearizada','KF','','','','lqg_sem_ruido','Animação do modelo linearizado (Filtro de Kalman)',Q_lqr,R_lqr,Q_w,Q_v,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 51 - LQG (Filtro de Kalman + LQR) com ruído no sistema.
% for i = 0:dt:tempo_simulacao
%     u = -K*xhat;                       % Lei de controle do sistema estocástico - KF
%     u2 = -K*x2;                        % Lei de controle da esperança
%     x = x + (A*x + B*u + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w  - KF
%     x2 = x2 + (A*x2 + B*u2)*dt;        % Discretização diferencial pelo método de Euler Direto - Esperança 
%     y = C*x;                           % Medição com ruído v
%     y2 = C*x2;                         % Medição da esperança
%     % Predição
%     xhat  = xhat + (A*xhat + B*u)*dt;  % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F2*P*F2'+ Q;                    % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     X2 = [X2 x2];                      % Acumulando o estado
%     Y2 = [Y2 y2];                      % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     U2 = [U2 u2];                      % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Linearizado','Esperança Linearizada','KF','','','','lqg_w','Animação do modelo linearizado (Filtro de Kalman)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 52 - LQG (Filtro de Kalman + LQR) com ruído no sistema e nos sensores.
% for i = 0:dt:tempo_simulacao
%     u = -K*xhat;                       % Lei de controle do sistema estocástico - KF
%     u2 = -K*x2;                        % Lei de controle da esperança
%     x = x + (A*x + B*u + w)*dt;        % Discretização diferencial pelo método de Euler Direto com ruído w  - KF
%     x2 = x2 + (A*x2 + B*u2)*dt;            % Discretização diferencial pelo método de Euler Direto - Esperança 
%     y = C*x + v;                       % Medição com ruído v
%     y2 = C*x2;                         % Medição da esperança
%     % Predição
%     xhat  = xhat + (A*xhat + B*u)*dt;  % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F2*P*F2'+ Q;                    % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     X2 = [X2 x2];                      % Acumulando o estado
%     Y2 = [Y2 y2];                      % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     U2 = [U2 u2];                      % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end
 
% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Linearizado','Esperança Linearizada','KF','','','','lqg_w_v','Animação do modelo linearizado (Filtro de Kalman)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 53 - EKF (Filtro de Kalman Estendido) sem ruídos.
% for i = 0:dt:tempo_simulacao   
%     u = 0;                                 % Lei de controle do sistema estocástico 
%     u2 = 0;                                % Lei de controle da esperança 
%     x =  x + f(u,x(2),x(3),x(4))*dt;       % Discretização diferencial pelo método de Euler Direto com ruído w
%     x2 =  x2 + f(u2,x2(2),x2(3),x2(4))*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = C*x;                               % Medição do sistema estocástico Não-Linear
%     y2 = C*x2;                             % Medição da esperança Não-Linear
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;              % Estado predito
%     yhat = C*xhat;                                              % Medição do estado  predito
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Cáculo da covariância 
%     % Resíduos
%     S = C*P*C'+R;                          % Resíduo da covariância
%     yt = y - yhat;                         % Resíduo de medição 
%     % Atualização
%     L = P*C'*inv(S);                       % Ganho de Kalman
%     xhat = xhat + L*(yt);                  % Estado Atualizado
%     yhat = C*xhat;                         % Medição do estado atualizado
%     P = P - L*C*P;                         % Atualização da covariância   
%     % Acumuladores
%     X = [X x];                             % Acumulando o sistema estocástico
%     Y = [Y y];                             % Acumulando a medição do sistema estocástico
%     X2 = [X2 x2];                          % Acumulando a esperança do sistema
%     Y2 = [Y2 y2];                          % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle do sistema estocástico
%     U2 = [U2 u2];                          % Acumulando a lei de controle da esperança
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','Esperança Não-Linear','EKF','','','','ekf_w_sem_lqr','Animação do modelo não-linear (Filtro de Kalman Estendido)',zeros(4),0,Q_w,Q_v,sigma_quadrado_w,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 54 - EKF (Filtro de Kalman Estendido) com ruído aditivo no sistema.
% for i = 0:dt:tempo_simulacao   
%     u = 0;                                 % Lei de controle do sistema estocástico 
%     u2 = 0;                                % Lei de controle da esperança 
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     x2 =  x2 + f(u2,x2(2),x2(3),x2(4))*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = C*x;                               % Medição do sistema estocástico Não-Linear
%     y2 = C*x2;                             % Medição da esperança Não-Linear
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;              % Estado predito
%     yhat = C*xhat;                                              % Medição do estado  predito
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Cáculo da covariância 
%     % Resíduos
%     S = C*P*C'+R;                          % Resíduo da covariância
%     yt = y - yhat;                         % Resíduo de medição 
%     % Atualização
%     L = P*C'*inv(S);                       % Ganho de Kalman
%     xhat = xhat + L*(yt);                  % Estado Atualizado
%     yhat = C*xhat;                         % Medição do estado atualizado
%     P = P - L*C*P;                         % Atualização da covariância   
%     % Acumuladores
%     X = [X x];                             % Acumulando o sistema estocástico
%     Y = [Y y];                             % Acumulando a medição do sistema estocástico
%     X2 = [X2 x2];                          % Acumulando a esperança do sistema
%     Y2 = [Y2 y2];                          % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle do sistema estocástico
%     U2 = [U2 u2];                          % Acumulando a lei de controle da esperança
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','Esperança Não-Linear','EKF','','','','ekf_w_sem_lqr','Animação do modelo não-linear (Filtro de Kalman Estendido)',zeros(4),0,Q_w,Q_v,sigma_quadrado_w,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 55 - EKF (Filtro de Kalman Estendido) com ruidos aditivos no sistema e nos sensores.
% for i = 0:dt:tempo_simulacao   
%     u = 0;                                 % Lei de controle do sistema estocástico 
%     u2 = 0;                                % Lei de controle da esperança 
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     x2 =  x2 + f(u2,x2(2),x2(3),x2(4))*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = C*x + v;                           % Medição do sistema estocástico Não-Linear
%     y2 = C*x2;                             % Medição da esperança Não-Linear
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;              % Estado predito
%     yhat = C*xhat;                                              % Medição do estado  predito
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Cáculo da covariância 
%     % Resíduos
%     S = C*P*C'+R;                          % Resíduo da covariância
%     yt = y - yhat;                         % Resíduo de medição 
%     % Atualização
%     L = P*C'*inv(S);                       % Ganho de Kalman
%     xhat = xhat + L*(yt);                  % Estado Atualizado
%     yhat = C*xhat;                         % Medição do estado atualizado
%     P = P - L*C*P;                         % Atualização da covariância   
%     % Acumuladores
%     X = [X x];                             % Acumulando o sistema estocástico
%     Y = [Y y];                             % Acumulando a medição do sistema estocástico
%     X2 = [X2 x2];                          % Acumulando a esperança do sistema
%     Y2 = [Y2 y2];                          % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle do sistema estocástico
%     U2 = [U2 u2];                          % Acumulando a lei de controle da esperança
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Não-Linear','Estimativa com EKF','Esperança Não-Linear',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','Esperança Não-Linear','EKF','','','','ekf_w_v_sem_lqr','Animação do modelo não-linear (Filtro de Kalman Estendido)',zeros(4),0,Q_w,Q_v,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 56 - EKF (Filtro de Kalman Estendido) com LQR sem ruídos. 
% for i = 0:dt:tempo_simulacao   
%     u = -K*xhat;                           % Lei de controle do sistema estocástico 
%     u2 = -K*x2;                            % Lei de controle da esperança 
%     x =  x + (f(u,x(2),x(3),x(4)))*dt;     % Discretização diferencial pelo método de Euler Direto com ruído w
%     x2 =  x2 + f(u2,x2(2),x2(3),x2(4))*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = C*x;                               % Medição do sistema estocástico Não-Linear
%     y2 = C*x2;                             % Medição da esperança Não-Linear
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;              % Estado predito
%     yhat = C*xhat;                                              % Medição do estado  predito
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Cáculo da covariância 
%     % Resíduos
%     S = C*P*C'+R;                          % Resíduo da covariância
%     yt = y - yhat;                         % Resíduo de medição 
%     % Atualização
%     L = P*C'*inv(S);                       % Ganho de Kalman
%     xhat = xhat + L*(yt);                  % Estado Atualizado
%     yhat = C*xhat;                         % Medição do estado atualizado
%     P = P - L*C*P;                         % Atualização da covariância   
%     % Acumuladores
%     X = [X x];                             % Acumulando o sistema estocástico
%     Y = [Y y];                             % Acumulando a medição do sistema estocástico
%     X2 = [X2 x2];                          % Acumulando a esperança do sistema
%     Y2 = [Y2 y2];                          % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle do sistema estocástico
%     U2 = [U2 u2];                          % Acumulando a lei de controle da esperança
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','Esperança Não-Linear','EKF','','','','ekf_com_lqr','Animação do modelo não-linear (Filtro de Kalman Estendido)',Q_lqr,R_lqr,Q_w,Q_v,0,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 57 - EKF (Filtro de Kalman Estendido) com LQR e com ruído no sistema.
% for i = 0:dt:tempo_simulacao   
%     u = -K*xhat;                           % Lei de controle do sistema estocástico 
%     u2 = -K*x2;                            % Lei de controle da esperança 
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     x2 =  x2 + f(u2,x2(2),x2(3),x2(4))*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = C*x;                               % Medição do sistema estocástico Não-Linear
%     y2 = C*x2;                             % Medição da esperança Não-Linear
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;              % Estado predito
%     yhat = C*xhat;                                              % Medição do estado  predito
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Cáculo da covariância 
%     % Resíduos
%     S = C*P*C'+R;                          % Resíduo da covariância
%     yt = y - yhat;                         % Resíduo de medição 
%     % Atualização
%     L = P*C'*inv(S);                       % Ganho de Kalman
%     xhat = xhat + L*(yt);                  % Estado Atualizado
%     yhat = C*xhat;                         % Medição do estado atualizado
%     P = P - L*C*P;                         % Atualização da covariância   
%     % Acumuladores
%     X = [X x];                             % Acumulando o sistema estocástico
%     Y = [Y y];                             % Acumulando a medição do sistema estocástico
%     X2 = [X2 x2];                          % Acumulando a esperança do sistema
%     Y2 = [Y2 y2];                          % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle do sistema estocástico
%     U2 = [U2 u2];                          % Acumulando a lei de controle da esperança
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end

% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','Esperança Não-Linear','EKF','','','','ekf_w_com_lqr','Animação do modelo não-linear (Filtro de Kalman Estendido)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,0,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 58 - EKF (Filtro de Kalman Estendido) com LQR e ruídos no sistema e no sensor.
% for i = 0:dt:tempo_simulacao   
% u = -K*xhat;                               % Lei de controle do sistema estocástico 
%     u2 = -K*x2;                            % Lei de controle da esperança 
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     x2 =  x2 + f(u2,x2(2),x2(3),x2(4))*dt; % Discretização diferencial pelo método de Euler Direto com ruído w
%     y = C*x + v;                           % Medição do sistema estocástico Não-Linear
%     y2 = C*x2;                             % Medição da esperança Não-Linear
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;              % Estado predito
%     yhat = C*xhat;                                              % Medição do estado  predito
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Cáculo da covariância 
%     % Resíduos
%     S = C*P*C'+R;                          % Resíduo da covariância
%     yt = y - yhat;                         % Resíduo de medição 
%     % Atualização
%     L = P*C'*inv(S);                       % Ganho de Kalman
%     xhat = xhat + L*(yt);                  % Estado Atualizado
%     yhat = C*xhat;                         % Medição do estado atualizado
%     P = P - L*C*P;                         % Atualização da covariância   
%     % Acumuladores
%     X = [X x];                             % Acumulando o sistema estocástico
%     Y = [Y y];                             % Acumulando a medição do sistema estocástico
%     X2 = [X2 x2];                          % Acumulando a esperança do sistema
%     Y2 = [Y2 y2];                          % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle do sistema estocástico
%     U2 = [U2 u2];                          % Acumulando a lei de controle da esperança
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end
 
% plotar_sistema('dois_com_filtragem',t,X,Y,U,t,Xhat,Yhat,[],t,X2,Y2,U2,Yt,[],Px,[],L_k,[],espessura_linha,'Medição Linear','Estimativa com KF','Esperança Linearizada',tamanho_legenda,tamanho_titulo) % Plot do sistema
% animar_pendulo(3,Y',Y2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','Esperança Não-Linear','EKF','','','','ekf_w_v_com_lqr','Animação do modelo não-linear (Filtro de Kalman Estendido)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str(rad2deg(ci),4),mat2str(rad2deg(ci2),4),mat2str(rad2deg(ci3),4),limites_grafico,tempo_simulacao,0); % Animação

%% 59 - EKF e LQG (Filtro de Kalman Estendido e Regulador Linear Quadrático Gaussiano) sem ruídos. 
% dt = 0.02;               % Redefinição período de amostragem para testes
% F = F_j;                 % Matriz de Transição de estados para o sistema não-linear
% f = f_h;                 % Vetor de estados não-lineares
% F2 = eye(4) + A*dt;      % Matriz de Transição de estados para o sistema linearizado
% tempo_simulacao = 5;     % Redefinição do tempo de simulação para testes
% for i = 0:dt:tempo_simulacao 
%     u = -K*xhat;                             % Lei de controle - EKF
%     u2 = -K*xhat2;                           % Lei de controle - KF 
%     u3 = -K*x3;                              % Lei de controle - Esperança Não-Linear
%     u4 = -K*x4;                              % Lei de controle - Esperança Linearizada 
%     x =  x + f(u,x(2),x(3),x(4))*dt;         % Discretização da equação diferencial pelo método de Euler Direto - Não-Linear
%     x2 = x2 + (A*x2 + B*u2)*dt;              % Discretização da equação diferencial pelo método de Euler Direto - Linearizado
%     x3 = x3 + f(u3,x3(2),x3(3),x3(4))*dt;    % Discretização da equação diferencial pelo método de Euler Direto - Esperança não-linear
%     x4 = x4 + (A*x4 + B*u4)*dt;             % Discretização da equação diferencial pelo método de Euler Direto - Esperança linearizada 
%     y =  C*x;                                % Medição - EKF
%     y2 = C*x2;                               % Medição - KF
%     y3 = C*x3;                               % Medição - Esperança Não-Linear
%     y4 = C*x4;                               % Medição - Esperança Linearizada 
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;   % Estado predito - EKF
%     xhat2  = xhat2 + (A*xhat2 + B*u2)*dt;                        % Estado predito - KF 
%     yhat = C*xhat;                                   % Medição do estado predito - EKF     
%     yhat2 = C*xhat2;                                 % Medição do estado predito - KF
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Predição da covariância  - EKF
%     P2  = F2*P2*F2'+ Q;                              % Predição da covariância do erro - KF
%     % Resíduos
%     S = C*P*C'+R;                         % Resíduo da covariância - EKF
%     S2 = C*P2*C'+R;                       % Resíduo da covariância - KF
%     yt = y - yhat;                        % Resíduo de medição - EKF 
%     yt2 = y2 - yhat2;                     % Resíduo de medição - KF 
%     % Atualização
%     L = P*C'*inv(S);                      % Ganho de Kalman - EKF
%     L2 = P2*C'*inv(S2);                   % Ganho de Kalman - KF
%     xhat = xhat + L*(yt);                 % Estado Atualizado - EKF
%     xhat2 = xhat2 + L2*(yt2);             % Estado Atualizado - KF
%     yhat = C*xhat;                        % Medição do estado atualizado - EKF
%     yhat2 = C*xhat2;                      % Medição do estado atualizado - KF
%     P = (I - L*C)*P;                      % Atualização da covariância - EKF 
%     P2 = (I - L2*C)*P2;                   % Atualização da covariância - KF 
%     %Acumuladores
%     X = [X x];                            % Acumulando o sistema estocástico - Não-linear
%     X2 = [X2 x2];                         % Acumulando o sistema estocástico - Linearizado
%     X3 = [X3 x3];                         % Acumulando a esperança - Não-linear
%     X4 = [X4 x4];                         % Acumulando a esperança - Linearizada
%     Y = [Y y];                            % Acumulando a medição do sistema estocástico Não-Linear
%     Y2 = [Y2 y2];                         % Acumulando a medição do sistema estocástico Linearizada
%     Y3 = [Y3 y3];                         % Acumulando a medição da esperança Não-Linear
%     Y4 = [Y4 y4];                         % Acumulando a medição da esperança Linearizada
%     Xhat = [Xhat xhat];                   % Acumulando a estimativa atualizada - EKF
%     Xhat2 = [Xhat2 xhat2];                % Acumulando a estimativa atualizada - KF
%     Yhat = [Yhat yhat];                   % Medição da estimativa atualizada - EKF
%     Yhat2 = [Yhat2 yhat2];                % Medição da estimativa atualizada - KF
%     Px = [Px diag(P)];                    % Acumulando a covariância atualizada - EKF
%     Px2 = [Px2 diag(P2)];                 % Acumulando a covariância atualizada - KF
%     L_k = [L_k L];                        % Acumulando o Ganho de Kalman - EKF
%     L_k2 = [L_k2 L2];                     % Acumulando o Ganho de Kalman - KF
%     Yt = [Yt yt];                         % Acumulando o erro de medição - EKF
%     Yt2 = [Yt2 yt2];                      % Acumulando o erro de medição - KF
%     U = [U u];                            % Acumulando a lei de controle - EKF
%     U2 = [U2 u2];                         % Acumulando a lei de controle - KF
%     t = [t i];                            % Acumulando o tempo 
% end

% Medição do sistema estocástico não-linear 
% plotar_sistema('kf_ekf_sem_ruido',t,X,Y,U,t,Xhat,Yhat,U2,t,Xhat2,Yhat2,[],Yt,Yt2,Px,Px2,L_k,L_k2,espessura_linha,'Medição (Sistema Não-Linear)','EKF','KF',tamanho_legenda,tamanho_titulo) % Plot do sistema

% Medição do sistema estocástico linearizado 
% plotar_sistema('kf_ekf_sem_ruido',t,X2,Y2,U2,t,Xhat2,Yhat2,U,t,Xhat,Yhat,[],Yt,Yt2,Px,Px2,L_k,L_k2,espessura_linha,'Medição (Sistema Linearizado)','KF','EKF',tamanho_legenda,tamanho_titulo) % Plot do sistema

%  limites_grafico = [-0.03 0.03 -0.025 0.025];
%  s = 0.03;

% animar_pendulo(3,Y',Yhat2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','KF','EKF','','','','ekf_lqr_lqg_seguidor','Animação dos modelos linearizado com \bf LQG \rm e não-linear com \bf LQR \rm e \bf Filtro de Kalman Estendido \rm.',Q_lqr,R_lqr,Q_w,Q_v,0,0,strcat(num2str(dt),'s'),mat2str([round(ci(1,1),3) ci(1,2) round(rad2deg(ci(1,3))) ci(1,4)] ,4),mat2str([round(ci3(1,1),3) ci3(1,2) round(rad2deg(ci3(1,3))) ci3(1,4)] ,4),mat2str([round(ci2(1,1),3) ci2(1,2) round(rad2deg(ci2(1,3))) ci2(1,4)] ,4),limites_grafico,tempo_simulacao,1); % Animação

%% 60 - EKF e LQG (Filtro de Kalman Estendido e Regulador Linear Quadrático Gaussiano) com ruídos no sistema e na medição. 
% dt = 0.001;            % Redefinição período de amostragem para testes
% F = F_j;               % Matriz de Transição de estados para o sistema não-linear
% f = f_h;               % Vetor de estados não-lineares
% F2 = eye(4) + A*dt;    % Matriz de Transição de estados para o sistema linearizado
% tempo_simulacao = 5;   % Redefinição do tempo de simulação para testes
% for i = 0:dt:tempo_simulacao       
%     u = -K*xhat;                             % Lei de controle - EKF
%     u2 = -K*xhat2;                           % Lei de controle - KF 
%     u3 = -K*x3;                              % Lei de controle - Esperança Não-Linear
%     u4 = -K*x4;                              % Lei de controle - Esperança Linearizada 
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt;   % Discretização da equação diferencial pelo método de Euler Direto - Não-Linear
%     x2 = x2 + ((A*x2 + B*u2) + w)*dt;        % Discretização da equação diferencial pelo método de Euler Direto - Linearizado
%     x3 =  x3 + f(u3,x3(2),x3(3),x3(4))*dt;   % Discretização da equação diferencial pelo método de Euler Direto - Esperança não-linear
%     x4 = x4 + (A*x4 + B*u4)*dt;              % Discretização da equação diferencial pelo método de Euler Direto - Esperança linearizada 
%     y = C*x + v;                             % Medição - Estocástico Não-Linear
%     y2 = C*x2 + v;                           % Medição - Estocástico Linearizado
%     y3 = C*x3;                               % Medição - Esperança Não-Linear
%     y4 = C*x4;                               % Medição - Esperança Linearizada 
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt;   % Estado predito - EKF
%     xhat2  = xhat2 + (A*xhat2 + B*u2)*dt;            % Estado predito - KF 
%     yhat = C*xhat;                                   % Medição do estado predito - EKF     
%     yhat2 = C*xhat2;                                 % Medição do estado predito - KF
%     P  = F(dt,u,xhat(3),xhat(4))*P*F(dt,u,xhat(3),xhat(4))'+ Q; % Predição da covariância  - EKF
%     h = P;
%     P2  = F2*P2*F2'+ Q;                              % Predição da covariância do erro - KF
%     % Resíduos
%     S = C*P*C'+R;                         % Resíduo da covariância - EKF
%     S2 = C*P2*C'+R;                       % Resíduo da covariância - KF
%     yt = y - yhat;                        % Resíduo de medição - EKF 
%     yt2 = y2 - yhat2;                     % Resíduo de medição - KF 
%     % Atualização
%     L = P*C'*inv(S);                      % Ganho de Kalman - EKF
%     L2 = P2*C'*inv(S2);                   % Ganho de Kalman - KF
%     xhat = xhat + L*(yt);                 % Estado Atualizado - EKF
%     xhat2 = xhat2 + L2*(yt2);             % Estado Atualizado - KF
%     yhat = C*xhat;                        % Medição do estado atualizado - EKF
%     yhat2 = C*xhat2;                      % Medição do estado atualizado - KF
%     P = P - L*C*P;                        % Atualização da covariância - EKF 
%     P2 = P2 - L2*C*P2;                    % Atualização da covariância - KF 
%     %Acumuladores
%     X = [X x];                            % Acumulando o sistema estocástico - Não-linear
%     X2 = [X2 x2];                         % Acumulando o sistema estocástico - Linearizado
%     X3 = [X3 x3];                         % Acumulando a esperança - Não-linear
%     X4 = [X4 x4];                         % Acumulando a esperança - Linearizada
%     Y = [Y y];                            % Acumulando a medição do sistema estocástico Não-Linear
%     Y2 = [Y2 y2];                         % Acumulando a medição do sistema estocástico Linearizada
%     Y3 = [Y3 y3];                         % Acumulando a medição da esperança Não-Linear
%     Y4 = [Y4 y4];                         % Acumulando a medição da esperança Linearizada
%     Xhat = [Xhat xhat];                   % Acumulando a estimativa atualizada - EKF
%     Xhat2 = [Xhat2 xhat2];                % Acumulando a estimativa atualizada - KF
%     Yhat = [Yhat yhat];                   % Medição da estimativa atualizada - EKF
%     Yhat2 = [Yhat2 yhat2];                % Medição da estimativa atualizada - KF
%     Px = [Px diag(P)];                    % Acumulando a covariância atualizada - EKF
%     Px2 = [Px2 diag(P2)];                 % Acumulando a covariância atualizada - KF
%     L_k = [L_k L];                        % Acumulando o Ganho de Kalman - EKF
%     L_k2 = [L_k2 L2];                     % Acumulando o Ganho de Kalman - KF
%     H = [H h];                            % Acumulando a covariância  predita
%     Yt = [Yt yt];                         % Acumulando o erro de medição - EKF
%     Yt2 = [Yt2 yt2];                      % Acumulando o erro de medição - KF
%     U = [U u];                            % Acumulando a lei de controle - EKF
%     U2 = [U2 u2];                         % Acumulando a lei de controle - KF
%     U3 = [U3 u3];                         % Acumulando a lei de controle - KF
%     U4 = [U4 u4];                         % Acumulando a lei de controle - KF
%     t = [t i];                            % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);    % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);    % Criando um novo o ruído no no sensor
% end

% Medição do sistema estocástico não-linear 
% plotar_sistema('kf_ekf_com_ruido',t,X,Y,U,t,Xhat,Yhat,U2,t,Xhat2,Yhat2,[],Yt,Yt2,Px,Px2,L_k,L_k2,espessura_linha,'Medição do Sistema Não-Linear','EKF','KF',tamanho_legenda,tamanho_titulo) % Plot do sistema

% Medição do sistema estocástico linearizado 
% plotar_sistema('kf_ekf_com_ruido',t,X2,Y2,U2,t,Xhat2,Yhat2,U,t,Xhat,Yhat,[],Yt,Yt2,Px,Px2,L_k,L_k2,espessura_linha,'Medição do Sistema Linearizado','KF','EKF',tamanho_legenda,tamanho_titulo) % Plot do sistema

% animar_pendulo(6,Y',Y3',Yhat',Y2',Y4',Yhat2',s,l,l_carrinho,h_carrinho,'Medido Não-Linear','Esperança Não-Linear','EKF','Medido Linearizado','Esperança Linearizada','KF','ekf_lqr_lqg_w_v','Animação dos modelos linearizado e não-linear (Filtros de Kalman e Kalman Estendido)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str([round(ci(1,1),1) ci(1,2) rad2deg(ci(1,3)) ci(1,4)] ,4),mat2str([round(ci3(1,1),1) ci3(1,2) rad2deg(ci3(1,3)) ci3(1,4)] ,4),mat2str([round(ci2(1,1),1) ci2(1,2) rad2deg(ci2(1,3)) ci2(1,4)] ,4),limites_grafico,1); % Animação
% animar_pendulo(3,Y',Yhat2',Yhat',[],[],[],s,l,l_carrinho,h_carrinho,'Medido Não-Linear','KF','EKF','','','','ekf_lqr_lqg_filtragem','Animação dos modelos linearizado e não-linear (Filtros de Kalman e Kalman Estendido)',Q_lqr,R_lqr,Q_w,Q_v,sigma_quadrado_w,sigma_quadrado_v,strcat(num2str(dt),'s'),mat2str([round(ci(1,1),1) ci(1,2) rad2deg(ci(1,3)) ci(1,4)] ,4),mat2str([round(ci3(1,1),1) ci3(1,2) rad2deg(ci3(1,3)) ci3(1,4)] ,4),mat2str([round(ci2(1,1),1) ci2(1,2) rad2deg(ci2(1,3)) ci2(1,4)] ,4),limites_grafico,tempo_simulacao,1); % Animação


% Plot da covariância predita
% plot([H(1,1:4:end)',H(2,1:4:end)',H(3,1:4:end)',H(4,1:4:end)'],'lineWidth',2);
% legend('Posição do Carrinho','Velocidade do Carrinho','Ângulo da Haste','Velocidade da Haste')
% ylabel('Covariância Predita','FontSize',18);           % Label do eixo y
% xlabel('Tempo(s)','FontSize',18);                 % Label do eixo x
% title('Gráfico da Covariância Predita vs. Tempo','FontSize',18);  % Titulo do diagrama de polos e zeros
% grid on; 

%% 61 - Função Densidade de Probabilidades.
% % Código para a criação dó grafico de uma PDF  com média mu e desvio padrão sigma
% 
% % Código dos caracteres especiais utilizados no gráfico
% % char(160) - Espaço
% % char(40) - (
% % char(41) - )
% % char(963) - sigma
% % char(178) - ^2
% % char(58) - :
%             
% variancia = 1e-3;                 % Arbitrando a variância
% mu = media;                       % Arbitrando a média
% sigma = round(sqrt(variancia),4); % Cálculo do desvio padrão                                
% 
% x_length = 1002;                                     % Tamanho do vetor
% x = linspace(mu - 4*sigma, mu + 4*sigma, x_length);  % Cria um vetor de x_length posições iniciando em mu - 4*sigma até mu + 4*sigma
% j = pdf('Normal',x,mu,sigma);                        % Distribuição normal de x com média mu e desvio padrão sigma
% y = pdf('Normal',mu,mu,sigma);                       % Valor da pdf em mu
% y1 = pdf('Normal',mu-sigma,mu,sigma);                % Valor da pdf em mu-sigma que é igual ao valor da pdf em mu+sigma  
% y2 = pdf('Normal',mu-2*sigma,mu,sigma);              % Valor da pdf em mu-2*sigma que é igual ao valor da pdf em mu+2*sigma  
% y3 = pdf('Normal',mu-3*sigma,mu,sigma);              % Valor da pdf em mu-3*sigma que é igual ao valor da pdf em mu+3*sigma  
% y4 = pdf('Normal',mu-4*sigma,mu,sigma);              % Valor da pdf em mu-4*sigma que é igual ao valor da pdf em mu+4*sigma              
%           
% tamanho = x_length/8;  % Tamanho do desvio padrão a ser mostrado no gráfico
%             
% % Cria os pontos no eixo x para os desvios padrões 
% m1 = linspace(mu-sigma,mu+sigma,2*tamanho);     % Cria espaço para 2 desvios padrão entre mu-sigma e mu+sigma
% m2_e = linspace(mu-2*sigma,mu-sigma,tamanho);   % Cria espaço para 1 desvio padrão entre mu-2*sigma e mu-sigma
% m2_d = linspace(mu+sigma,mu+2*sigma,tamanho);   % Cria espaço para 1 desvio padrão entre mu+2*sigma e mu+sigma
% m3_e = linspace(mu-3*sigma,mu-2*sigma,tamanho); % Cria espaço para 1 desvio padrão entre mu-3*sigma e mu-2*sigma
% m3_d = linspace(mu+2*sigma,mu+3*sigma,tamanho); % Cria espaço para 1 desvio padrão entre mu+3*sigma e mu+2*sigma
% m4_e = linspace(mu-4*sigma,mu-3*sigma,tamanho); % Cria espaço para 1 desvio padrão entre mu-4*sigma e mu-3*sigma
% m4_d = linspace(mu+3*sigma,mu+4*sigma,tamanho); % Cria espaço para 1 desvio padrão entre mu+4*sigma e mu+3*sigma
%             
% % Cria os pontos em y relacionadas aos eixo x através da pdf
% j1 = pdf('Normal',m1,mu,sigma);        % Cria os pontos da curva baseados em m1
% j2_e = pdf('Normal',m2_e,mu,sigma);    % Cria os pontos da curva baseados em m2_e (esquerda)
% j2_d = pdf('Normal',m2_d,mu,sigma);    % Cria os pontos da curva baseados em m2_d (direita)
% j3_e = pdf('Normal',m3_e,mu,sigma);    % Cria os pontos da curva baseados em m3_e (esquerda)
% j3_d = pdf('Normal',m3_d,mu,sigma);    % Cria os pontos da curva baseados em m3_d (direita)
% j4_e = pdf('Normal',m4_e,mu,sigma);    % Cria os pontos da curva baseados em m4_e (esquerda)
% j4_d = pdf('Normal',m4_d,mu,sigma);    % Cria os pontos da curva baseados em m4_d (direita)
%           
% subplot(2,1,1);                        % Põe o gráfico da pdf na parte de cima
% % Cria as áreas a serem pintadas no gráfico da pdf
% a = zeros(4,1);                        % Cria um vetor para colocar as áreas
% a(4) = area([m4_e m4_d],[j4_e j4_d]);  % Pinta a área acima de m4_e m4_d até os pontos em j4_e j4_d da curva
% hold on;                               % Permite colocar as outras áreas no mesmo gráfico
% a(3) = area([m3_e m3_d],[j3_e j3_d]);  % Pinta a área acima de m3_e m3_d até os pontos em j3_e j3_d da curva
% a(2) = area([m2_e m2_d],[j2_e j2_d]);  % Pinta a área acima de m2_e m2_d até os pontos em j2_e j2_d da curva
% a(1) = area(m1,j1);                    % Pinta a área acima de m1_e m1_d até os pontos em j1_e j1_d da curva
% hold off;
%            
% lgd = legend(a, '68.26%','27.7%','4.003%','0.005%'); % Cria a legenda 
% lgd.FontSize = 20;
% 
% % RGB copiado do site https://www.rapidtables.com/web/color/RGB_Color.html
% RGB = validatecolor({'#99CCFF','#66B2FF','#3399FF','#0080FF'},'multiple'); % Cores escolhidas tons de azul
%  
% % Muda a cor de cada área para as cores escolhidas
% set(a(1),'FaceColor',RGB(1,:)); % Pinta a área a(1)
% set(a(2),'FaceColor',RGB(2,:)); % Pinta a área a(2)
% set(a(3),'FaceColor',RGB(3,:)); % Pinta a área a(3)
% set(a(4),'FaceColor',RGB(4,:)); % Pinta a área a(4)     
%  
% ax = gca; % Pega o objeto eixo do gráfico 
% % ax.XTick = [mu-4*sigma mu-3*sigma  mu-2*sigma mu-sigma mu mu+sigma mu+2*sigma  mu+3*sigma mu+4*sigma];      % Marcações do eixo x 
% % ax.XTickLabel = {'-4\sigma','-3\sigma','-2\sigma','-\sigma','\mu','\sigma','2\sigma','3\sigma','4\sigma'};  % Labels do eixo x 
% % ax.YTick = [];         % Apaga as marcações do eixo y
% % ax.YTickLabel = [];    % Apaga os labels do eixo y
  
% ax.XTick = [mu-4*sigma mu-3*sigma  mu-2*sigma mu-sigma mu mu+sigma mu+2*sigma  mu+3*sigma mu+4*sigma];        % Muda os marcadores do eixo x
% ax.XTickLabel = {num2str(-4*sigma), num2str(-3*sigma),  num2str(-2*sigma), num2str(-sigma), num2str(mu), num2str(sigma), num2str(2*sigma),  num2str(3*sigma), num2str(4*sigma)}; % Coloca os valores em cada marcador no eixo x
% ax.FontSize = 20;        % Fonte do texto do gráfico
 
%% 62 - Ruido Branco Gaussiano.
% % Código que cria um gráfico que contém um ruído de média mu e desvio
% % padrão sigma abaixo do gráfico da FDP
% 
% % variancia = 1e-4;         % Variância do ruído
% % mu = media;               % Média do ruído
% % sigma = sqrt(variancia);  % Desvio Padrão do ruído
% 
% x_length = 1002;                                       % Tamanho do vetor x na horizontal do gráfico
% x = 1:x_length;                                        % Vetor x de tamanho x_length
% % x = linspace(mu - 4*sigma, mu + 4*sigma, x_length);  % Vetor de 4000 posições
% vetor_ruido =  mu + sigma * randn(length(x),1);        % Cria o ruído propriamente dito
% subplot(2,1,2);                                        % Coloca o gráfico do ruído abaixo do gráfico da FDP
% p = plot(x,vetor_ruido);                               % Plota o ruído e aloca os pontos em p
% grid on;                                               % Habilita a grade no gráfico
% ylabel('Desvio Padrão - \sigma ');                     % Coloca um texo vertical no eixo y
% xlabel('Amostras');                                    % Coloca um texo horizontal no eixo x 
% ax = gca;                                              % Extrai os eixos do gráfico
% ax.YTick = [-4*sigma -3*sigma  -2*sigma -sigma mu sigma 2*sigma  3*sigma  4*sigma ]; % Reordena a marcação do eixo vertical do gráfico dividindo em sigmas
% ax.YTickLabel = {num2str(-4*sigma), num2str(-3*sigma),num2str(-2*sigma), num2str(-sigma), num2str(mu), num2str(sigma),num2str(2*sigma),  num2str(3*sigma), num2str(4*sigma)};% Coloca o texto nas marcações do eixo horizontal
% ax.YLim = [-4*sigma 4*sigma];                          % Altera os valores limites do eixo vertical do gráfico
% ax.XLim = [0 x_length];                                % Altera os valores limites do eixo horizontal do gráfico 
% ax.FontSize = 20;
