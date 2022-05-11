%% 1 - Preparação do Script.
close all;
clear all;
clc;

%% 2 - Definição da Variáveis.
syms x1 x2 x3 x4 u l M m g J dt
x = [x1 x2 x3 x4];           

%% 3 - Especificações Quanser IP02 (Long).
l = 0.641; % Comprimento da Haste
m = 0.230; % Massa do Pêndulo
M = 0.57;  % Massa do Carrinho
g = 9.81;  % Aceleração da Gravidade
J =  m*(l/2)^2; % Momento de Inércia de uma Barra Delgada do
l_carrinho = 0.45; % Largura do Carrinho
h_carrinho = 0.15; % Altura do Carrinho
w_carrinho = 0.061; % comprimento
% As especificações a seguir são fornecidas em resposta a um ponto
% de ajuste de posição do carrinho de onda quadrada de ± 30 mm.
% Tr <= 1.5s % Tempo de subida menor ou igual a 1.5 segundos
% mod(alpha) <= 1 degree % alpha menor ou igual a 1 grau
% V < 10v % Tensão menor que 10 volts
% Carrinho +- 30mm (0.03m ou 3cm)
% Rack dimensions (L x W x H) 102 x 15 x 6.1 cm


%% 4 - Condições Iniciais.
ci  = [0 0 pi/6 0]; % Condições iniciais do sistema
ci2 = [0 0 pi/180 0]; % Condições iniciais do estimador de estados

%% 5 - Equações de Estados Não-Lineares.
f = [x(2);
    (u+m*l*x(4)^2*sin(x(3))-m*g*sin(x(3))*cos(x(3)))/(M+m-m*cos(x(3))^2);
     x(4);
    (u*cos(x(3))-(M+m)*g*sin(x(3))+m*l*x(4)^2*sin(x(3))*cos(x(3)))/(m*l*cos(x(3))^2-(M+m)*l)];
h = [x(1);x(3)];

%% 6 - Handles Não-Lineares.
f_h = matlabFunction(f); % Transforma f simbólico em um vetor de funções handles.
h_h = matlabFunction(h); % Transforma h simbólico em um vetor de funções handles.


%% 7 - Jacobiano das Equações de Estado.
% Jfx = simplify(jacobian(f,x)); 
% Jfu = simplify(jacobian(f,u));
% Jgx = simplify(jacobian(h,x));
% Jgu = simplify(jacobian(h,u));

%% 8 - Jacobianos para as Matrizes de transição do EKF.
f_j = simplify(jacobian(x' + f*dt,x));
g_j = simplify(jacobian(x' + f*dt,u));
h_j = simplify(jacobian(h,x));

%% 9 - Linearização.
% x0 = [1 0 0 0]; % Ponto de equilíbrio
% u0 = 0;

%% 10 - Matrizes A = Jfx_lin, B = Jfu_lin, C = Jgx_lin e D = Jgu_lin.
% Jfx_lin = subs(Jfx,{x1,x2,x3,x4},x0);
% Jfu_lin = subs(Jfu,{x1,x2,x3,x4},x0);
% Jhx_lin = subs(Jhx,{x1,x2,x3,x4},x0);
% Jhu_lin = subs(Jhu,{x1,x2,x3,x4},x0);

%% 11 - Matrizes das Equações de Estado Linearizadas.
A = [ 0 1 0 0; 0 0 -(g*m)/M 0; 0 0 0 1; 0 0 -(g*(M + m))/(-l*M) 0];
B = [0; 1/M; 0; -(1/(l*M))];
C = [1 0 0 0; 0 0 1 0];
D = [0;0];

%% 12 - Matrizes de Transição - Jacobianos de (x + f(x,u)*dt) em relação a x e u, e de h em relação a x.
% F_j = @(x,u,dt)[1 dt 0 0;
%     0 1 -(dt*(14743*x(4)^2*cos(x(3))-451260*cos(x(3))^2+225630))/(23000*cos(x(3))^2-80000)-(4600*dt*cos(x(3))*sin(x(3))*((14743*sin(x(3))*x(4)^2)/100000+u-(22563*cos(x(3))*sin(x(3)))/10000))/(23*cos(x(3))^2-80)^2 (14743*dt*x(4)*sin(x(3)))/(500*(23*sin(x(3))^2 + 57));
%     0 0 1 dt;
%     0 0 (4600000*dt*cos(x(3))*sin(x(3))*((14743*cos(x(3))*sin(x(3))*x(4)^2)/100000-(981*sin(x(3)))/125+u*cos(x(3))))/(641*(23*cos(x(3))^2-80)^2)-(dt*(784800*cos(x(3))-29486*x(4)^2*cos(x(3))^2+100000*u*sin(x(3))+14743*x(4)^2))/(14743*cos(x(3))^2-51280) (14743*dt*x(4)*sin(2*x(3)))/(100000*((14743*cos(2*x(3)))/200000 - 87817/200000)) + 1];       
% G_j = @(x,dt)[0;-(100*dt)/(23*cos(x(3))^2 - 80);0;(100000*dt*cos(x(3)))/(14743*cos(x(3))^2 - 51280)];
% H_j = [1 0 0 0; 0 0 1 0];

F_j = matlabFunction(f_j);
G_j = matlabFunction(g_j);
H_j = matlabFunction(h_j);

%% 13 - Variáveis Temporais.
tempo_simulacao = 10;
dt_graficos = 0.001;
dt_animacao = 0.06;

%% 14 - Matriz de transição para o sistema linearizado.
% Pelo ZOH. Mapeia corretamente do plano s para o plano z.
% F = expm(A*dt);
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B;

% Pelo método de Euler. Pode mapear do plano s para o plano z incorretamente.
% F = eye(4,4)+dt*A;
% G = dt*B

%% 15 - Teste de Controlabilidade.
cc = ctrb(A,B);
%cc = [B A*B (A^2)*B (A^3)*B] % n-1 = 3, n = 4
[linhas,~] = size(cc);
rank(cc);

%% 16 - Teste de Observabilidade.
ob = obsv(A,C);
%ob = [C; C*A; C*(A^2); C*(A^3)] % n-1 = 3, n = 4
[linhas,colunas] = size(ob);
rank(ob); % Como o posto = 4 = n, o sistema é observável

%% 17 - Projeto do LQR.
Q = [  100 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
R = 1;
K = lqr( A, B, Q, R);  % Matriz de Ganho de realimentação
%% 18 - Tentativa frustrada de encontrar o período de amostragem pelo controle digital.
% syms s z T
% 
% % Sistema em malha aberta
% sys_ss = ss(A,B,C,D);       % Sistema em malha aberta representado no espaço de estados
% [num,den] = ss2tf(A,B,C,D); % Sistema mimo 1x2 em malha aberta representado em função de transferência
% sys1 = tf(num(1,:),den);    % Função de transferência para a posição do carrinho
% sys2 = tf(num(2,:),den);    % Função de transferência para o ângulo da haste
% sys = [sys1;sys2];          % Vetor coluna com os dois sistemas

% r = roots(den);
% 
% % Espansão em frações parciais
% fp =  vpa(-0.0544/(s+4.6346) + 0.0544/(s-4.6346) + 1.0880e-16/s + 1.25/(s^2),4);
% % Tranformada Z do sistema em malha aberta
% fpz = vpa((-0.0544*z)/(z-exp(-4.6346*T)) +(0.0544*z)/(z+exp(4.6346*T)) + 1.0880e-16 + (1.25*T*z)/((z-1)^2));
% 
% % Sistema em malha aberta no domínio Z
% [numz,denz] =numden(fpz);
% numz = collect(numz);
% denz = collect(denz);
% Gz = vpa(numz/denz,4);
% % Numerador de (1 + KG(z))
% den_tz = (1 + K(1)*Gz);
% [num,den] = numden(den_tz);
% num = vpa(collect(num),4);    % Coloca em odem decrescente de potência
% a = 4.6346;
% T = 0.336181;
% den_tz = z^4*exp(a*T) - z^3*(2*exp(a*T) - 1.544*exp(2*a*T) + 12.5*T*exp(a*T) + 0.456) + z^2*(12.5*T - 3.088*exp(2*a*T) - 12.5*T*exp(2*a*T) + 0.912) + z*(2*exp(a*T) + 1.544*exp(2*a*T) + 12.5*T*exp(a*T) - 0.456) - 1.0*exp(a*T);
% den_tz = sym2poly(den_tz);
% roots(den_tz);
% vetor = [exp(a*T)  -(2*exp(a*T) - 1.544*exp(2*a*T) + 12.5*T*exp(a*T) + 0.456) (12.5*T - 3.088*exp(2*a*T) - 12.5*T*exp(2*a*T) + 0.912) (2*exp(a*T) + 1.544*exp(2*a*T) + 12.5*T*exp(a*T) - 0.456) -exp(a*T)]
% dd = subs(vetor)
% abs(vpa(roots(dd),4))
% r = solve(den_tz,z, 'MaxDegree', 4)

%% 19 - Sistema em malha aberta.
[num,den] = ss2tf(A,B,C,D); % Sistema mimo 1x2 em malha aberta representado em função de transferência
sys1 = tf(num(1,:),den);    % Função de transferência para a posição do carrinho
sys2 = tf(num(2,:),den);    % Função de transferência para o ângulo da haste
sys = [sys1;sys2];          % Vetor coluna com os dois sistemas

%% 20 - Sistema em malha fechada.
sysmf = feedback(sys,[K(1) K(3)]); % Sistema em malha fechada

%% 21 - Cáculo da Largura de Banda.
omega_bw1 = bandwidth(sysmf(1)); % Valor da largura de banda pelo comando do matlab
% omega_bw2 = bandwidth(sysmf(2)); % Deu infinito
% bode(sysmf(1));   % Valor da largura de banda por inspeção do diagrama de bode do sistema em malha fechada: 3.67 rad/s
% bode(sysmf(2));   % Não foi possivel determinar a largura de banda

%% 22 - Cálculo do Período de Amostragem pelo Cálculo Numerico, Teorema de Nyquist e valores usuais.
% Todos os exemplos encontrados de estabilidade numérica não envolviam 
% matrizes, logo os resultados para a estabilidade marginal aqui 
% apresentados foram deduzidos intuitivamente, utilizando o módulo dos
% auto-valores de (A-BK), que forneceram o resultado esperado quando 
% alterado qualquer valor da matriz Q do LQR, com exceção do elemento 
% Q(3,3), onde dt_marginal não tornou o algorítmo marginalmente estável,
% conforme esperado, por um motivo ainda desconhecido.

% Valores para a estabilidade marginal
I = eye(4);                        % Cria uma matriz identidade 4x4
M = diag(eig(A-B*K));
% S = 2*I*vpa(inv(abs(M)),6); % Faz 2*I dividido por uma matriz diagonal que contém o módulo dos autovalores 
S = 2/(norm(M,2)); % Faz 2*I dividido por uma matriz diagonal que contém o módulo dos autovalores 
% dt_marginal = S(1,1);              % Pega o elemento (1,1) de S que forneceu o valor de dt para a estabilidade marginal

% Limite segundo o Teorema de Amostragem de Nyquist
omega_s = 2*omega_bw1;             % Frequência de Nyquist
dt_nyquist = (2*pi)/(omega_s);     % Máximo período de amostragem segundo o Teorema de Nyquist

% Valor ideal de 20 a 40 omega_bw (Franklin. pg. 504)
omega_s_min = 20*omega_bw1;        % Valor mínimo usual da frequência de amostragem
omega_s_max = 40*omega_bw1;        % Valor máximo usual da frequência de amostragem
omega_s_escolhido = 40*omega_bw1;% Valor escolhido para a frequancia de amostragem dentro da faixa usual
dt_max = (2*pi)/(omega_s_min);     % Período de amostragem máximo usual
dt_min = (2*pi)/(omega_s_max);     % Período de amostragem mínimo usual
dt_escolhido = (2*pi)/(omega_s_escolhido);   % Príodo de amostragem escolhido
% dt = dt_simulacao;                 % Período de amostragem escolhido para as simulações
dt = dt_animacao;                  % Período de amostragem escolhido para as animações
%% 23 - Definindo os parâmetros probabilísticos.

% Variância
sigma_quadrado_w = 0.005;          % Covariância do ruído do sistema
sigma_quadrado_v = 0.005;          % Covariância do ruído dos sensores

% Desvio Padrão
sigma_w = sqrt(sigma_quadrado_w);  % Desvio padrão do ruído do sistema
sigma_v = sqrt(sigma_quadrado_v);  % Desvio padrão do ruído do sistema

bias = 0;                          % Acurácia
media = 0;                         % Média

% Matriz de covariância        
P = eye(4,4);                      % Criando a matriz de covariância entre o sistema e as estimativas               
Px = diag(P);                      % Criando o vetor acumulador de diagonal principal das matrizes covariância

% Matrizes de Incertezas
Q_w = dt*sigma_w^2*eye(4,4);       % Matriz de incertezas do processo
Q_v = sigma_v^2*eye(2,2)/dt;       % Matriz de incertezas da medição

%% 24 - Variáveis utilizadas no Método Iterativo.
x = ci';             % Condições iniciais
X = x;               % Inicializando o acumulador de estados
y = [x(1);x(3)];     % Inicializando o a saída
Y = y;               % Inicializando o acumulador de saídas
x2 = ci';            % Condições iniciais do estimador
X2 = x2;             % Inicializando o acumulador de condições iniciais estimadas (Não utilizado)
y2 = [x2(1);x2(3)];  % Saída do estimador (Não utilizado)
Y2 = y2;             % Inicializando acumulador de saídas estimadas (Não utilizado)
xhat = ci2';         % Estado estimado
Xhat = xhat;         % Inicializando o acumulador de estados estimados
yhat = [xhat(1);xhat(3)];  % Saída estimada 
Yhat = yhat;               % Inicializando o acumulador de saídas estimadas
yt = y - yhat;       % Erro de estimativa inicial
Yt = yt;             % Acumulando os erro de estimativa
u = -K*x;            % Sinal de controle inicial
U = u;               % Acumulando o sinal de controle 
u2 = -K*x2;          % Sisnal de controle inicial do sistema 2
U2 = u2;             % Acumulando o sinal de controle do sistema 2
Q = Q_w;             % Matriz de incertezas do processo
R = Q_v;             % Matriz de incertezas das medições
l_k = P*C'*inv(C*P*C'+R); % Covariância do erro inicial
L_k = l_k;                % Inicializando o acumulador de covariância do erro
t = 0;                    % Tempo inicial
w =  media + sigma_w * randn(4,1); % Criando o ruído no sistema inicial
v =  media + sigma_v * randn(2,1); % Criando o ruído no no sensor inicial

%% 25 - Simulação do sistema não-linear sem LQR (Utilizando comandos do Matlab).
% t0 = 0;                        % Tempo inicial
% tf = tempo_simulacao;          % Tempo final
% K = [0 0 0 0];                 % Controle cancelado
% opts = odeset('RelTol',1e-5);  % Erro relativo
% [t1,x1] = ode45(@equacoes_nao_lineares,[t0 tf],ci,opts,K); % Ou -> ode45(@(t,x)equacoes_nao_lineares_lqr(t,x,K),[to tf],ci,opts);
% u1 = -K*x1';                   % Sinal de controle
% plotar_sistema(t1',x1',x1(:,[1 3])',u1',[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(x1(:,[1 3]),[],2,l,l_carrinho,h_carrinho,'Simulação do sistema não-linear sem LQR (Utilizando comandos do Matlab)',''); % Animação

%% 26 - Simulação do sistema linearizado sem LQR (Utilizando comandos do Matlab).
% tf = 10;                         % Tempo final
% ci = [0 0 pi/180 0];             % Condições iniciais
% K = [0 0 0 0];                   % Controle cancelado
% sys = ss(A-B*K,B,C,D);           % Sistema em malha fechada linearizado
% [y1,t1,x1] = initial(sys,ci,tf); % Simula o sistema linearizado
% u1 = -K*x1';                     % Sinal de controle
% plotar_sistema(t1',x1',x1(:,[1 3])',u1',[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(x1(:,[1 3]),[],2,l,l_carrinho,h_carrinho,'Simulação do sistema não-linear sem LQR (Utilizando comandos do Matlab)',''); % Animação

%% 27 - Simulação dos sistemas linear e não linear com LQR (Utilizando comandos do Matlab).
% t0 = 0;                          % Tempo inicial
% tf = tempo_simulacao;            % Tempo final
% [t1,x1] = ode45(@equacoes_nao_lineares,[t0 tf],ci,[],K); % Simulando o sistema não linear
% u1 = -K*x1';                     % Sinal de controle do sistem não-linear
% sys = ss(A-B*K,B,C,D);           % Sistema em malha fechada linearizado
% [y2,t2,x2] = initial(sys,ci,tf); % Simula o sistema linearizado
% u2 = -K*x2';                     % Sinal de controle do sistem linearizado
% plotar_sistema(t1',x1',x1(:,[1 3])',u1',t2',x2',y2',u2',[],[],2,'Não-Linear','Linearizado'); % Plot do sistema
% animar_pendulo(x1(:,[1 3]),[],2,l,l_carrinho,h_carrinho,'Simulação dos sistemas linear e não linear com LQR.',''); % Animação

%% 28 - Sistema não-linear (real e sua esperança) sem controle e (com e sem) ruído.
for i = 0:dt:tempo_simulacao
      u = -K*x;                                    % Lei de controle da esperança do sistema 
      u2 = -K*x2;                                  % Lei de controle do sistema real 
      x = x + (f_h(u,x(2),x(3),x(4)) + w)*dt;      % Estados do sistema (Método de Euler)
      x2 = x2 + f_h(u2,x2(2),x2(3),x2(4))*dt;      % Estados do sistema (Método de Euler)
      y = h_h(x(1),x(3)) + v;                      % Adicionando o ruído na medição
      y2 = h_h(x2(1),x2(3));                       % Adicionando o ruído na medição
      X = [X x];                                   % Acumulando o estado
      X2 = [X2 x2];                                % Acumulando o estado
      Y = [Y y];                                   % Acumulando a saída
      Y2 = [Y2 y2];                                % Acumulando a saída
      U = [U u];                                   % Sinal de controle da esperança do sistema 
      U2 = [U2 u2];                                % Sinal de controle do sistema real
      t = [t i];                                   % Acumulando o tempo
      w =  media + sigma_w * randn(4,1);           % Criando o ruído no sistema
      v =  media + sigma_v * randn(2,1);           % Criando o ruído no no sensor
end
% plotar_sistema(t,X,Y,U,t,X2,Y2,U2,[],[],2,'Sistema Medido','Esperança do Sistema');     % Plot do sistema.
animar_pendulo(Y',Y2',2,l,l_carrinho,h_carrinho,'Sistema Real','Esperança do Sistema'); % Animação

%% 29 - Sistema não-linear sem controle e sem ruído.
% for i = 0:dt:tempo_simulacao
%       u = 0;                             % Lei de controle 
%       x = x + f_h(u,x(2),x(3),x(4))*dt;  % x[n+1] = x[n]+f(x[n],u[n])*dt
%       y = h_h(x(1),x(3));                % Saída do sistema
%       X = [X x];                         % Acumulando o estado
%       Y = [Y y];                         % Acumulando a saída
%       U = [U u];                         % Sinal de controle
%       t = [t i];                         % Acumulando o tempo
% end
% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema Não-Linear sem controle e sem ruído',''); % Animação

%% 30 - Sistema não-linear com controle LQR e sem ruídos.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                            % Lei de controle
%     x = x + (f_h(u,x(2),x(3),x(4)))*dt;  % x[n+1] = x[n]+f(x[n],u[n])*dt
%     y = h_h(x(1),x(3));                  % Medição
%     X = [X x];                           % Acumulando o estado
%     Y = [Y y];                           % Acumulando a saída
%     U = [U u];                           % Acumulando o sinal de controle
%     t = [t i];                           % Acumulando o Tempo
% end
% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema não-linear com controle LQR e sem ruídos.',''); % Animação

%% 31 - Sistema não-linear com controle LQR e com ruido apenas no sistema.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                               % Lei de controle
%     x = x + (f_h(u,x(2),x(3),x(4)) + w)*dt; % x[n+1] = x[n]+f(x[n],u[n])*dt com ruído w
%     y = h_h(x(1),x(3));                     % Medição
%     X = [X x];                              % Acumulando o estado
%     Y = [Y y];                              % Acumulando a saída
%     U = [U u];                              % Acumulando o sinal de controle
%     t = [t i];                              % Acumulando o Tempo
%     w =  media + sigma_w * randn(4,1);      % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1);      % Criando um novo ruído no no sensor     
% end

% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema não-linear com controle LQR e com ruido apenas no sistema.',''); % Animação

%% 32 - Sistema não-linear com controle LQR e com ruido no sistema e no sensor.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                                  % Lei de controle
%     x = x + (f_h(u,x(2),x(3),x(4)) + w)*dt;    % x[n+1] =  x[n]+f(x[n],u[n])*dt com ruído w
%     y = h_h(x(1),x(3)) + v;                    % Medição com ruído v
%     X = [X x];                                 % Acumulando o estado
%     Y = [Y y];                                 % Acumulando a saída
%     U = [U u];                                 % Acumulando o sinal de controle
%     t = [t i];                                 % Acumulando o tempo
%     w =  media + sigma_w * randn(1,1);         % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1);         % Criando um novo ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema não-linear com controle LQR e com ruido no sistema e no sensor.',''); % Animação
 
%% 33 - Sistema linearizado sem controle e sem ruido.  
% for i = 0:dt:tempo_simulacao
%     u = 0;                      % Lei de controle  
%     x = x + (A*x + B*u)*dt;     % x[n+1] = (Ax[n] + Bu[n])*dt
%     y = C*x;                    % Medição
%     X = [X x];                  % Acumulando o estado
%     Y = [Y y];                  % Acumulando a saída
%     U = [U u];                  % Acumulando o sinal de controle
%     t = [t i];                  % Acumulando o tempo
% end
% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema linearizado sem controle e sem ruido .',''); % Animação

%% 34 - Sistema linearizado com controle LQR e sem ruido. 
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                   % Lei de controle 
%     x = x + (A*x + B*u)*dt;     % x[n+1] = (Ax[n] + Bu[n])*dt
%     y = C*x;                    % Medição
%     X = [X x];                  % Acumulando o estado
%     Y = [Y y];                  % Acumulando a saída
%     U = [U u];                  % Acumulando o sinal de controle
%     t = [t i];                  % Acumulando o tempo
% end
% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema linearizado com controle LQR e sem ruido.',''); % Animação

%% 35 - Sistema linearizado com controle LQR e com ruido aditivo no sistema.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                          % Lei de controle 
%     x = x + (A*x + B*u + w)*dt;        % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = C*x;                           % Medição
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída
%     U = [U u];                         % Acumulando o sinal de controle
%     t = [t i];                         % Acumulando o tempo
%     w =  media + sigma_w * randn(1,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema linearizado com controle LQR e com ruido aditivo no sistema.',''); % Animação

%% 36 - Sistema linearizado com controle LQR e com ruido no sistema e no sensor.
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                          % Lei de controle
%     x = x + (A*x + B*u + w)*dt;        % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = h_h(x(1),x(3)) + v;            % Medição com ruído v
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída
%     U = [U u];                         % Acumulando o sinal de controle
%     t = [t i];                         % Acumulando o tempo
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,[],[],[],[],[],[],2,'',''); % Plot do sistema
% animar_pendulo(Y',[],2,l,l_carrinho,h_carrinho,'Sistema linearizado com controle LQR e com ruido no sistema e no sensor.',''); % Animação

%% 37 - Sistemas Não-Linear e Linearizado com LQR e sem ruídos simulados juntos. 
% for i = 0:dt:tempo_simulacao
%     u = -K*x;                               % Lei de controle
%     u2 = -K*x2;                             % Lei de controle
%     x = x + (A*x + B*u)*dt;                 % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     x2 = x2 + f_h(u2,x2(2),x2(3),x2(4))*dt; % x[n+1] =  x[n]+f(x[n],u[n])*dt com ruído w
%     y = C*x;                                % Medição com ruído v
%     y2 = h_h(x2(1),x2(3));                  % Medição com ruído v
%     X = [X x];                              % Acumulando o estado
%     X2 = [X2 x2];                           % Acumulando o estado
%     Y = [Y y];                              % Acumulando a saída
%     Y2 = [Y2 y2];                           % Acumulando a saída
%     U = [U u];                              % Acumulando o sinal de controle
%     U2 = [U2 u2];                           % Acumulando o sinal de controle
%     t = [t i];                              % Acumulando o tempo
% end
% % plotar_sistema(t,X,Y,U,t,X2,Y2,U2,[],[],2,'Linearizado','Não-Linear'); % Plot do sistema
% animar_pendulo(Y',Y2',2,l,l_carrinho,h_carrinho,'Linearizado','Não-Linear'); % Animação

%% 38 - KF - Filtro de Kalman seguindo o sistema linearizado sem ruido.
% F = expm(A*dt);                    % Matriz de transição estados
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B; % Matriz de transição de saída
% 
% for i = 0:dt:tempo_simulacao
%     u = 0;                             % Lei de controle
%     x = x+(A*x + B*u)*dt;              % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = C*x;                           % Medição
%     % Predição
%     xhat  = F*x + G*u;                 % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F*P*F'+ Q;                    % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Equação de atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação
% 

%% 39 - KF - Filtro de Kalman seguindo o sistema linearizado com ruido no sistema.
% F = expm(A*dt);                    % Matriz de transição estados
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B; % Matriz de transição de saída
% 
% for i = 0:dt:tempo_simulacao
%     u = 0;                             % Lei de controle
%     x = x+(A*x + B*u + w)*dt;          % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = C*x;                           % Medição
%     % Predição
%     xhat  = F*x + G*u;                 % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F*P*F'+ Q;                    % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Equação de atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 40 - KF - Filtro de Kalman seguindo o sistema linearizado com ruido no sistema e no sensor.
% F = expm(A*dt);                        % Matriz de transição estados
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B; % Matriz de transição de saída
% 
% for i = 0:dt:tempo_simulacao
%     u = 0;                             % Lei de controle
%     x = x+(A*x + B*u + w)*dt;          % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = C*x + v;                       % Medição com ruído v
%     % Predição
%     xhat  = F*x + G*u;                 % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F*P*F'+ Q;                    % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Equação de atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 41 - LQG (Filtro de Kalman + LQR) sem ruído.
% F = expm(A*dt);                        % Matriz de transição estados
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B; % Matriz de transição de saída
% 
% for i = 0:dt:tempo_simulacao
%     u = -K*xhat;                       % Lei de controle
%     x = x+(A*x + B*u)*dt;              % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = C*x;                           % Medição
%     % Predição
%     xhat  = F*x + G*u;                 % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F*P*F'+ Q;                    % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Equação de atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 42 - LQG (Filtro de Kalman + LQR) com ruído no sistema.
% F = expm(A*dt);                        % Matriz de transição estados
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B; % Matriz de transição de saída
% 
% for i = 0:dt:tempo_simulacao
%     u = -K*xhat;                       % Lei de controle
%     x = x+(A*x + B*u + w)*dt;          % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = C*x;                           % Medição
%     % Predição
%     xhat  = F*x + G*u;                 % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F*P*F'+ Q;                    % Predição da covariância do erro
%     % Resíduos
%     S = C*P*C'+ R;                     % Resíduo da covariância
%     yt = y - yhat;                     % Resíduo da medição
%     % Atualização   
%     L = P*C'*inv(S);                   % Ganho de Kalman
%     xhat = xhat + L*(yt);              % Equação da atualização do estado estimado
%     yhat = C*xhat;                     % Medição do estado atualizado
%     P = P - L*C*P;                     % Equação de atualização da covariância
%     % Acumuladores
%     X = [X x];                         % Acumulando o estado
%     Y = [Y y];                         % Acumulando a saída não corrompida
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 43 - LQG (Filtro de Kalman + LQR) com ruído no sistema e nos sensores.
% F = expm(A*dt);                                         % Matriz de transição estados
% G = integral(@(t) expm(A*t),0,dt,'ArrayValued',true)*B; % Matriz de transição de saída
% 
% for i = 0:dt:tempo_simulacao
%     u = -K*xhat;                       % Lei de controle
%     x = x+(A*x + B*u + w)*dt;          % x[n+1] = (Ax[n] + Bu[n] + w)*dt 
%     y = C*x + v;                       % Medição com ruído v
%     % Predição
%     xhat  = F*x + G*u;                 % Predição do estado estimado 
%     yhat = C*xhat;                     % Medição do estado predito
%     P  = F*P*F'+ Q;                    % Predição da covariância do erro
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
%     Xhat = [Xhat xhat];                % Acumulando o estado estimado
%     Yhat = [Yhat yhat];                % Acumulando a saída estimada
%     Px = [Px diag(P)];                 % Acumulando a covariância atualizada
%     L_k = [L_k L];                     % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                      % Acumulando o erro de medição
%     U = [U u];                         % Acumulando a lei de controle
%     t = [t i];                         % Acumulando o tempo diferencial
%     w =  media + sigma_w * randn(4,1); % Criando um novo ruído no sistema
%     v =  media + sigma_v * randn(2,1); % Criando um novo ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado');
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 44 - EKF (Filtro de Kalman Estendido) sem ruídos.
% F = F_j;  % Matriz de Transiçãp de estados
% f = f_h;  % Equações de estados não-lineares

% for i = 0:dt:tempo_simulacao   
%     u = 0;                                 % Lei de controle  
%     x =  x + f(u,x(2),x(3),x(4))*dt;       % x[n+1] = f(u[n],x[n])*dt
%     y = C*x;                               % Medição
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
%     X = [X x];                             % Acumulando a esperança do sistema
%     Y = [Y y];                             % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle
%     t = [t i];                             % Acumulando o tempo 
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 45 - EKF (Filtro de Kalman Estendido) com ruído aditivo no sistema.
% F = F_j;  % Matriz de Transição de estados
% f = f_h;  % Equações de estados não-lineares

% for i = 0:dt:tempo_simulacao   
%     u = 0;                                 % Lei de controle  
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % x[n+1] = f(u[n],x[n])*dt
%     y = C*x;                               % Medição
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
%     X = [X x];                             % Acumulando a esperança do sistema
%     Y = [Y y];                             % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 46 - EKF (Filtro de Kalman Estendido) com ruidos aditivos no sistema e nos sensores.
% F = F_j;  % Matriz de Transição de estados
% f = f_h;  % Equações de estados não-lineares

% for i = 0:dt:tempo_simulacao   
%     u = 0;                                 % Lei de controle  
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % x[n+1] = f(u[n],x[n])*dt
%     y = C*x + v;                           % Medição com ruído v
%     % Predição
%     xhat = xhat + (f(u,xhat(2),xhat(3),xhat(4))+w)*dt;          % Estado predito
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
%     X = [X x];                             % Acumulando a esperança do sistema
%     Y = [Y y];                             % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 47 - EKF (Filtro de Kalman Estendido) com LQR sem ruídos. 
% F = F_j;  % Matriz de Transição de estados
% f = f_h;  % Equações de estados não-lineares

% for i = 0:dt:tempo_simulacao   
%     u = -K*xhat;                           % Lei de controle  
%     x =  x + f(u,x(2),x(3),x(4))*dt;       % x[n+1] = f(u[n],x[n])*dt
%     y = C*x;                               % Medição
%     % Predição
%     xhat = xhat + f(u,xhat(2),xhat(3),xhat(4))*dt; % Estado predito
%     yhat = C*xhat;                                 % Medição do estado  predito
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
%     X = [X x];                             % Acumulando a esperança do sistema
%     Y = [Y y];                             % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle
%     t = [t i];                             % Acumulando o tempo 
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 48 - EKF (Filtro de Kalman Estendido) com LQR e com ruído no sistema.
% F = F_j;  % Matriz de Transição de estados
% f = f_h;  % Equações de estados não-lineares
% 
% for i = 0:dt:tempo_simulacao   
%     u = -K*xhat;                           % Lei de controle  
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % x[n+1] = f(u[n],x[n])*dt
%     y = C*x;                               % Medição
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
%     X = [X x];                             % Acumulando a esperança do sistema
%     Y = [Y y];                             % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 49 - EKF (Filtro de Kalman Estendido) com LQR e ruídos no sistema e no sensor.
F = F_j;  % Matriz de Transição de estados
f = f_h;  % Equações de estados não-lineares
% for i = 0:dt:tempo_simulacao   
%     u = -K*xhat;                           % Lei de controle  
%     x =  x + (f(u,x(2),x(3),x(4)) + w)*dt; % x[n+1] = f(u[n],x[n])*dt
%     y = C*x + v;                           % Medição com ruído v
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
%     X = [X x];                             % Acumulando a esperança do sistema
%     Y = [Y y];                             % Acumulando a medição da esperança
%     Xhat = [Xhat xhat];                    % Acumulando a estimativa atualizada
%     Yhat = [Yhat yhat];                    % Medição da estimativa atualizada
%     Px = [Px diag(P)];                     % Acumulando a covariância atualizada
%     L_k = [L_k L];                         % Acumulando o Ganho de Kalman
%     Yt = [Yt yt];                          % Acumulando o erro de medição
%     U = [U u];                             % Acumulando a lei de controle
%     t = [t i];                             % Acumulando o tempo 
%     w =  media + sigma_w * randn(4,1);     % Criando um novo o ruído no sistema
%     v =  media + sigma_v * randn(2,1);     % Criando um novo o ruído no no sensor
% end
% plotar_sistema(t,X,Y,U,t,Xhat,Yhat,[],Yt,Px,2,'Medido','Estimado'); % Plot do sistema
% animar_pendulo(Y',Yhat',2,l,l_carrinho,h_carrinho,'Sistema Real','Sistema Estimado'); % Animação

%% 50 - Função Densidade de Probabilidades.
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
% variancia = 0.5;         % Arbitrando a variância
% mu = media;              % Arbitrando a média
% sigma = sqrt(variancia); % Cálculo do desvio padrão                                
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
% legend(a, '68.26%','27.7%','4.003%','0.005%'); % Cria a legenda              
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
% % ax.XTick = [mu-4*sigma mu-3*sigma  mu-2*sigma mu-sigma mu mu+sigma mu+2*sigma  mu+3*sigma mu+4*sigma];
% % ax.XTickLabel = {'-4\sigma','-3\sigma','-2\sigma','-\sigma','\mu','\sigma','2\sigma','3\sigma','4\sigma'};         
% 
% ax.XTick = [mu-4*sigma mu-3*sigma  mu-2*sigma mu-sigma mu mu+sigma mu+2*sigma  mu+3*sigma mu+4*sigma]; % Muda os marcadores do eixo x
% ax.XTickLabel = {num2str(-4*sigma), num2str(-3*sigma),  num2str(-2*sigma), num2str(-sigma), num2str(mu), num2str(sigma), num2str(2*sigma),  num2str(3*sigma), num2str(4*sigma)}; % Coloca os valores em cada marcador no eixo x
        
%% 51 - Ruido Branco Gaussiano.
% Código que cria um gráfico que contém um ruído de média mu e desvio
% padrão sigma abaixo do gráfico da FDP

% variancia = 0.5;         % Variância do ruído
% mu = media;              % Média do ruído
% sigma = sqrt(variancia); % Desvio Padrão do ruído
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

%% Plot do Ganho de Kalman

% y1 = L_k(:,1:2:end); 
% y2 = L_k(:,2:2:end);
%  
% subplot(2,2,1)
% p1_a = plot(t,y1(1,:)) 
% hold on;
% p1_b = plot(t,y2(1,:));
% grid on;
% title('Parte referente à posição do carrinho.')
% legend('Saida 1 - (Posição do Carrinho)','Saida 2 - (Ângulo da Haste)')
% hold off;
% 
% subplot(2,2,2)
% p2_a = plot(t,y1(3,:)) 
% hold on;
% p2_b = plot(t,y2(3,:));
% grid on;
% title('Parte referente ao ângulo da haste.')
% legend('Saida 1 - (Posição do Carrinho)','Saida 2 - (Ângulo da Haste)')
% hold off;
% 
% subplot(2,2,3)
% p3_a = plot(t,y1(2,:)) 
% hold on;
% p3_b = plot(t,y2(2,:));
% grid on;
% title('Parte referente á velocidade do carrinho.')
% legend('Saida 1 - (Posição do Carrinho)','Saida 2 - (Ângulo da Haste)')
% hold off;
% 
% subplot(2,2,4)
% p4_a = plot(t,y1(4,:)) 
% hold on;
% p4_b = plot(t,y2(4,:));
% grid on;
% title('Parte referente á velocidade angular da haste.')
% legend('Saida 1 - (Posição do Carrinho)','Saida 2 - (Ângulo da Haste)')
% hold off;
    
% L_1 = norm(L_k(:,1),2);
% L_2 = norm(L_k(:,2),2);
% for i = 3:2:length(L_k)
%  L_1 = [L_1 norm(L_k(:,i),2)];
%  L_2 = [L_2 norm(L_k(:,i+1),2)];
% end
% plot(t,L_1')
% hold on
% plot(t,L_2')
% hold off

