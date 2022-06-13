function plotar_sistema(tipo,t1,x1,y1,u1,t2,x2,y2,u2,t3,x3,y3,u3,e,e2,p,p2,L_k,L_k2,w,legenda1,legenda2,legenda3,v,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que plota os gráficos para os sistemas com ou sem filtro.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Glossário de variáveis
% tipo -> Escolha do tipo de simulação
% t1 -> Tempo de simulação para o sistema 1
% x1 -> Estados do sistema 1
% y1 -> Saida do sistema 1
% t2 -> Tempo de simulação para o sistema 2
% x2 -> Estados do sistema 2
% y2 -> Saida do sistema 2
% t3 -> Tempo de simulação para o sistema 3
% x3 -> Estados do sistema 3
% y3 -> Saida do sistema 3
% u1  -> Sinal de controle do sistema 1
% u2  -> Sinal de controle do sistema 1
% u3  -> Sinal de controle do sistema 1
% e1  -> Erro de estimativa entre a medição e o filtro 1 
% e1  -> Erro de estimativa entre a medição e o filtro 2 
% p   -> Covariância do erro para o filtro 1
% p2  -> Covariância do erro para o filtro 2
% legenda1 -> Legenda do sistema 1
% legenda2 -> Legenda do sistema 2
% legenda3 -> Legenda do sistema 3
% w  -> Espessura da linha do gráfico
% v  -> Tamanho da fonte da legenda e eixo y 
% m  -> Tamanho da fonte do título

fig1 =figure(1);             % cria uma nova figura
tl1 = tiledlayout(2,2);       % Layout da figura 2 linhas e 2 colunas
tl1.TileSpacing = 'compact';  % Diminui o espaço entre os graficos
tl1.Padding = 'compact';      % Diminui o espaço lateral dos gráficos
    

switch tipo
    case 'um_sem_lqr'      
        posicao_carrinho = nexttile;                                    % Gráfico da posição do carrinho
        p1 = plot(posicao_carrinho,t1,y1(1,:)','LineWidth',w);          % Plot da posição do carrinho do sistema 1 
        xlabel(posicao_carrinho,'Tempo [s]','FontSize',v)               % Texto do eixo x
        ylabel(posicao_carrinho,'Posição do Carrinho [m]','FontSize',v) % Texto do eixo y
        grid(posicao_carrinho,'on');                                    % Habilita a grade

        angulo_haste = nexttile;                                        % Gráfico da posição angular da haste
        p2 = plot(angulo_haste,t1,(180/pi)*y1(2,:)','LineWidth',w);     % Plot do ângulo da haste do sistema 1
        xlabel(angulo_haste,'Tempo [s]','FontSize',v)                   % Texto do eixo x
        ylabel(angulo_haste,strcat('Posição da Haste - [',char(176),']'),'FontSize',v) % Texto do eixo y
        grid(angulo_haste,'on');                                        % Habilita a grade  

        velocidade_carrinho = nexttile;                                 % Gráfico das velocidades
        p3 = plot(velocidade_carrinho,t1,x1(2,:),'LineWidth',w);        % Plot da velocidade linear do sistema 1   
        xlabel(velocidade_carrinho,'Tempo [s]','FontSize',v);           % Texto do eixo x
        ylabel(velocidade_carrinho,'Velocidade do Carrinho [m/s]','FontSize',v) % Texto do eixo y
        grid(velocidade_carrinho,'on');                                 % Habilita a grade

        velocidade_haste = nexttile;                                    % Gráfico do sinal de controle para um sistema
        p4 = plot(velocidade_haste,t1,(180/pi)*x1(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 1 
        xlabel(velocidade_haste,'Tempo [s]','FontSize',v);              % Texto do eixo x
        ylabel(velocidade_haste,'Velocidade da Haste - [graus/s]','FontSize',v); % Texto do eixo y
        grid(velocidade_haste,'on');                                    % Habilita a grade
    
    case 'um_com_lqr'  
        posicao_carrinho = nexttile;                                    % Gráfico da posição do carrinho
        p1 = plot(posicao_carrinho,t1,y1(1,:)','LineWidth',w);          % Plot da posição do carrinho do sistema 1 
        xlabel(posicao_carrinho,'Tempo [s]','FontSize',v)               % Texto do eixo x
        ylabel(posicao_carrinho,'Posição do Carrinho [m]','FontSize',v) % Texto do eixo y
        grid(posicao_carrinho,'on');                                    % Habilita a grade

        angulo_haste = nexttile;                                        % Gráfico da posição angular da haste
        p2 = plot(angulo_haste,t1,(180/pi)*y1(2,:)','LineWidth',w);     % Plot do ângulo da haste do sistema 1
        xlabel(angulo_haste,'Tempo [s]','FontSize',v)                   % Texto do eixo x
        ylabel(angulo_haste,strcat('Posição da Haste - [',char(176),']'),'FontSize',v) % Texto do eixo y
        grid(angulo_haste,'on');                                        % Habilita a grade  

        velocidades = nexttile;                                         % Gráfico das velocidades
        p3 = plot(velocidades,t1,x1(2,:),'LineWidth',w);                % Plot da velocidade linear do sistema 1
        hold(velocidades,'on');                                         % Retém o gráfico
        p4 = plot(velocidades,t1,(180/pi)*x1(4,:),'LineWidth',w);       % Plot da velocidade angular do sistema 1 
        legend(velocidades,'Velocidade Linear do Carrinho','Velocidade Angular da Haste','FontSize',v); % Legenda do gráfico
        xlabel(velocidades,'Tempo [s]','FontSize',v);                   % Texto do eixo x
        ylabel(velocidades,'Vel. Carrinho [m/s] - Haste [graus/s]','FontSize',v) % Texto do eixo y
        grid(velocidades,'on');                                         % Habilita a grade
        hold(velocidades,'off');  

        sinal_controle = nexttile;                                      % Gráfico do sinal de controle para um sistema
        p5 = plot(sinal_controle,t1,u1','LineWidth',w);                 % Plot do sinal de controle 
        xlabel(sinal_controle,'Tempo [s]','FontSize',v);                % Texto do eixo x
        ylabel(sinal_controle,'Sinal de Controle','FontSize',v);        % Texto do eixo y
        grid(sinal_controle,'on');                                      % Habilita a grade

    case 'dois_sem_filtragem_sem_lqr'  
        posicao_carrinho = nexttile;                                    % Gráfico da posição do carrinho
        p1_a = plot(posicao_carrinho,t1,y1(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 1
        hold(posicao_carrinho,'on');                                    % Retém o gráfico
        p1_b = plot(posicao_carrinho,t2,y2(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 2 
        legend(posicao_carrinho,'Velocidade Linear do Carrinho','Velocidade Angular da Haste','FontSize',v); % Legenda do gráfico
        xlabel(posicao_carrinho,'Tempo [s]','FontSize',v)               % Texto do eixo x
        ylabel(posicao_carrinho,'Posição do Carrinho [m]','FontSize',v) % Texto do eixo y
        grid(posicao_carrinho,'on');                                    % Habilita a grade
        hold(posicao_carrinho,'off');                                   % Libera o gráfico

        angulo_haste = nexttile;                                        % Gráfico da posição angular da haste
        p2_a = plot(angulo_haste,t1,(180/pi)*y1(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 1
        hold(angulo_haste,'on');                                        % Retém o gráfico
        p2_b = plot(angulo_haste,t2,(180/pi)*y2(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 2
        legend(angulo_haste,'Velocidade Linear do Carrinho','Velocidade Angular da Haste','FontSize',v); % Legenda do gráfico
        xlabel(angulo_haste,'Tempo [s]','FontSize',v)                   % Texto do eixo x
        ylabel(angulo_haste,strcat('Posição da Haste - [',char(176),']'),'FontSize',v) % Texto do eixo y
        grid(angulo_haste,'on');                                        % Habilita a grade 
        hold(angulo_haste,'off');                                       % Libera o gráfico

        velocidade_carrinho = nexttile;                                 % Gráfico das velocidades
        p3_a = plot(velocidade_carrinho,t1,x1(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 1 
        hold(velocidade_carrinho,'on');                                 % Retém o gráfico
        p3_b = plot(velocidade_carrinho,t2,x2(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 2   
        legend(velocidade_carrinho,'Velocidade Linear do Carrinho','FontSize',v); % Legenda do gráfico
        xlabel(velocidade_carrinho,'Tempo [s]','FontSize',v);           % Texto do eixo x
        ylabel(velocidade_carrinho,'Velocidade do Carrinho [m/s]','FontSize',v) % Texto do eixo y
        grid(velocidade_carrinho,'on');                                 % Habilita a grade
        hold(velocidade_carrinho,'off');                                % Libera o gráfico

        sinal_controle = nexttile;                                      % Gráfico do sinal de controle para um sistema
        p4_a = plot(sinal_controle,t1,(180/pi)*x1(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 1 
        hold(sinal_controle,'on');                                      % Retém o gráfico
        p4_b = plot(sinal_controle,t2,(180/pi)*x2(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 2 
        legend(sinal_controle,'Velocidade Linear do Carrinho','Velocidade Angular da Haste','FontSize',v); % Legenda do gráfico
        xlabel(sinal_controle,'Tempo [s]','FontSize',v);                % Texto do eixo x
        ylabel(sinal_controle,'Velocidade da Haste - [graus/s]','FontSize',v); % Texto do eixo y
        grid(sinal_controle,'on');                                      % Habilita a grade
        hold(sinal_controle,'off');                                     % Libera o gráfico

    case 'dois_sem_filtragem_com_lqr'  
        posicao_carrinho = nexttile;                                    % Gráfico da posição do carrinho
        p1_a = plot(posicao_carrinho,t1,y1(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 1
        hold(posicao_carrinho,'on');                                    % Retém o gráfico
        p1_b = plot(posicao_carrinho,t2,y2(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 2
        legend(posicao_carrinho,legenda1,legenda2,'FontSize',v);        % Legenda do gráfico
        xlabel(posicao_carrinho,'Tempo [s]','FontSize',v)               % Texto do eixo x
        ylabel(posicao_carrinho,'Posição do Carrinho [m]','FontSize',v) % Texto do eixo y
        grid(posicao_carrinho,'on');                                    % Habilita a grade
        hold(posicao_carrinho,'off');                                   % Libera o gráfico

        angulo_haste = nexttile;                                        % Gráfico da posição angular da haste
        p2_a = plot(angulo_haste,t1,(180/pi)*y1(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 1
        hold(angulo_haste,'on');                                        % Retém o gráfico
        p2_b = plot(angulo_haste,t2,(180/pi)*y2(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 2
        legend(angulo_haste,legenda1,legenda2,'FontSize',v);            % Legenda do gráfico
        xlabel(angulo_haste,'Tempo [s]','FontSize',v)                   % Texto do eixo x
        ylabel(angulo_haste,strcat('Posição da Haste - [',char(176),']'),'FontSize',v) % Texto do eixo y
        grid(angulo_haste,'on');                                        % Habilita a grade  
        hold(angulo_haste,'off');                                       % Libera o gráfico

        velocidades = nexttile;                                         % Gráfico das velocidades
        p3_a = plot(velocidades,t1,x1(2,:),'LineWidth',w);              % Plot da velocidade linear do sistema 1
        hold(velocidades,'on');                                         % Retém o gráfico
        p3_b = plot(velocidades,t1,(180/pi)*x1(4,:),'LineWidth',w);     % Plot da velocidade angular do sistema 1 
        p3_c = plot(velocidades,t2,x2(2,:),'LineWidth',w);              % Plot da velocidade linear do sistema 2 
        p3_d = plot(velocidades,t2,(180/pi)*x2(4,:),'LineWidth',w);     % Plot da velocidade angular do sistema 2 
        legend(velocidades,strcat(legenda1,char(160),'- Carrinho'),strcat(legenda1,char(160),'- Haste'),strcat(legenda2,char(160),'- Carrinho'),strcat(legenda2,char(160),'- Haste'),'FontSize',v); % Legenda do gráfico
        xlabel(velocidades,'Tempo [s]','FontSize',v);                   % Texto do eixo x
        ylabel(velocidades,'Vel. Carrinho [m/s] - Haste [graus/s]','FontSize',v) % Texto do eixo y
        grid(velocidades,'on');                                         % Habilita a grade  
        hold(velocidades,'off');                                        % Libera o gráfico

        sinal_controle = nexttile;                                      % Gráfico do sinal de controle 
        p4_a = plot(sinal_controle,t1,u1','LineWidth',w);               % Plot do sinal de controle para o sistema 1
        hold(sinal_controle,'on');                                      % Retém o gráfico
        p4_b = plot(sinal_controle,t2,u2','LineWidth',w);               % Plot do sinal de controle para o sistema 2
        legend(sinal_controle,legenda1,legenda2,'FontSize',v);          % Legenda do gráfico
        xlabel(sinal_controle,'Tempo [s]','FontSize',v);                % Texto do eixo x
        ylabel(sinal_controle,'Sinal de Controle','FontSize',v);        % Texto do eixo y
        grid(sinal_controle,'on');                                      % Habilita a grade  
        hold(sinal_controle,'off');                                     % Libera o gráfico
        
    case 'dois_com_filtragem'    
        posicao_carrinho = nexttile;                                    % Gráfico da posição do carrinho
        p1_a = plot(posicao_carrinho,t1,y1(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 1
        hold(posicao_carrinho,'on');                                    % Retém o gráfico
        p1_b = plot(posicao_carrinho,t2,y2(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 2 
        p1_c = plot(posicao_carrinho,t3,y3(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 3 
        legend(posicao_carrinho,legenda1,legenda2,legenda3,'FontSize',v); % Legenda do gráfico
        xlabel(posicao_carrinho,'Tempo [s]','FontSize',v)               % Texto do eixo x
        ylabel(posicao_carrinho,'Posição do Carrinho [m]','FontSize',v) % Texto do eixo y
        grid(posicao_carrinho,'on');                                    % Habilita a grade
        hold(posicao_carrinho,'off');                                   % Libera o gráfico

        angulo_haste = nexttile;                                        % Gráfico da posição angular da haste
        p2_a = plot(angulo_haste,t1,(180/pi)*y1(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 1
        hold(angulo_haste,'on');                                        % Retém o gráfico
        p2_b = plot(angulo_haste,t2,(180/pi)*y2(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 2
        p2_c = plot(angulo_haste,t3,(180/pi)*y3(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 3
        legend(angulo_haste,legenda1,legenda2,legenda3,'FontSize',v);   % Legenda do gráfico
        xlabel(angulo_haste,'Tempo [s]','FontSize',v)                   % Texto do eixo x
        ylabel(angulo_haste,strcat('Posição da Haste - [',char(176),']'),'FontSize',v) % Texto do eixo y
        grid(angulo_haste,'on');                                        % Habilita a grade 
        hold(angulo_haste,'off');                                       % Libera o gráfico

        velocidade_carrinho = nexttile;                                 % Gráfico das velocidades lineares
        p3_a = plot(velocidade_carrinho,t1,x1(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 1 
        hold(velocidade_carrinho,'on');                                 % Retém o gráfico 
        p3_b = plot(velocidade_carrinho,t2,x2(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 2   
        p3_c = plot(velocidade_carrinho,t3,x3(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 3   
        legend(velocidade_carrinho,legenda1,legenda2,legenda3,'FontSize',v); % Legenda do gráfico
        xlabel(velocidade_carrinho,'Tempo [s]','FontSize',v);           % Texto do eixo x
        ylabel(velocidade_carrinho,'Velocidade do Carrinho [m/s]','FontSize',v) % Texto do eixo y
        grid(velocidade_carrinho,'on');                                 % Habilita a grade
        hold(velocidade_carrinho,'off');                                % Libera o gráfico

        velocidade_angular = nexttile;                                  % Gráfico da velocidade angular da haste
        p4_a = plot(velocidade_angular,t1,(180/pi)*x1(4,:),'LineWidth',w);  % Plot da velocidade linear do sistema 1 
        hold(velocidade_angular,'on');                                      % Retém o gráfico
        p4_b = plot(velocidade_angular,t2,(180/pi)*x2(4,:),'LineWidth',w);  % Plot da velocidade linear do sistema 2   
        p4_c = plot(velocidade_angular,t3,(180/pi)*x3(4,:),'LineWidth',w);  % Plot da velocidade linear do sistema 3   
        legend(velocidade_angular,legenda1,legenda2,legenda3,'FontSize',v); % Legenda do gráfico
        xlabel(velocidade_angular,'Tempo [s]','FontSize',v);                % Texto do eixo x
        ylabel(velocidade_angular,'Velocidade do Carrinho [m/s]','FontSize',v) % Texto do eixo y
        grid(velocidade_angular,'on');                                      % Habilita a grade
        hold(velocidade_angular,'off');                                     % Libera o gráfico

        fig2 = figure(2);                                               % cria uma nova figura
        tl2 = tiledlayout(2,2);                                         % Layout da figura 2 linhas e 2 colunas
        tl2.TileSpacing = 'compact';                                    % Diminui o espaço entre os graficos
        tl2.Padding = 'compact';                                        % Diminui o espaço lateral dos gráficos

        sinal_controle = nexttile;                                      % Gráfico do sinal de controle para um sistema
        p5 = plot(sinal_controle,t1,u1','LineWidth',w);                 % Plot do sinal de controle  
        xlabel(sinal_controle,'Tempo [s]','FontSize',v);                % Texto do eixo x
        ylabel(sinal_controle,'Sinal de Controle','FontSize',v);        % Texto do eixo y
        grid(sinal_controle,'on');                                      % Habilita a grade

        erro_estimativa = nexttile;                                     % Gráfico do erro de estimativa 
        p6 = plot(erro_estimativa,t1,e','LineWidth',w);                 % Plot do erro de estimativa da posição do carrinho entre o sistema 1 e seus estimador  ,'-r'
        legend(erro_estimativa,'Posição do Carrinho','Ãngulo da Haste','FontSize',v); % Legenda do gráfico       
        xlabel(erro_estimativa,'Tempo [s]','FontSize',v);               % Texto do eixo x
        ylabel(erro_estimativa,'Erro de Estimativa','FontSize',v);      % Texto do eixo y
        grid(erro_estimativa,'on');                                     % Habilita a grade

        covariancia = nexttile;                                         % Gráfico da covariância
        p7 = plot(covariancia,t1,p','LineWidth',w);                     % Plot da covariância 
        legend(covariancia,'Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
        xlabel(covariancia,'Tempo [s]','FontSize',v);                   % Texto do eixo x
        ylabel(covariancia,'Covarância do Erro - EKF','FontSize',v);    % Texto do eixo y
        grid(covariancia,'on');                                         % Habilita a grade

        ganho_kalman = nexttile;                                        % Gráfico da norma 2 das duas colunas do Ganho de Kalman
        L_1 = norm(L_k(:,1),2);                                         % Norma 2 da primeira coluna de L_k
        L_2 = norm(L_k(:,2),2);                                         % Norma 2 da segunda coluna de L_k
         for i = 3:2:length(L_k)                                        % Começa o loop da terceira coluna até o fim
           L_1 = [L_1 norm(L_k(:,i),2)];                                % Acumula a norma das colunas ímpares de L_k
           L_2 = [L_2 norm(L_k(:,i+1),2)];                              % Acumula a norma das colunas pares de L_k
         end
         p8_a = plot(ganho_kalman,t1,L_1','LineWidth',w)                % Mostra L1 
         hold on                                                        % Retém o gráfico
         p8_b = plot(ganho_kalman,t1,L_2','LineWidth',w)                % Mostra L2  
         legend(ganho_kalman,'Posição do Carrinho','Ângulo da Haste','FontSize',v) % Legenda
         ylabel(ganho_kalman,'Ganho de Kalman - Norma 2','FontSize',v); % Texto do eixo y
         xlabel(ganho_kalman,'Tempo - s');                              % Texto do eixo x
         grid(ganho_kalman,'on');                                       % Habilita a grade 
         hold(ganho_kalman,'off')                                       % Libera o gráfico
           
    case 'kf_ekf'       
         posicao_carrinho = nexttile;                                    % Gráfico da posição do carrinho
         p1_a = plot(posicao_carrinho,t1,y1(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 1
         hold(posicao_carrinho,'on');                                    % Retém o gráfico
         p1_b = plot(posicao_carrinho,t2,y2(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 2
         p1_c = plot(posicao_carrinho,t3,y3(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 3
         legend(posicao_carrinho,legenda1,legenda2,legenda3,'FontSize',v); % Legenda do gráfico
         xlabel(posicao_carrinho,'Tempo [s]','FontSize',v)               % Texto do eixo x
         ylabel(posicao_carrinho,'Posição do Carrinho [m]','FontSize',v) % Texto do eixo y
         grid(posicao_carrinho,'on');                                    % Habilita a grade
         hold(posicao_carrinho,'off');                                   % Libera o gráfico
        
         angulo_haste = nexttile;                                        % Gráfico da posição angular da haste
         p2_a = plot(angulo_haste,t1,(180/pi)*y1(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 1
         hold(angulo_haste,'on');                                        % Retém o gráfico
         p2_b = plot(angulo_haste,t2,(180/pi)*y2(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 2
         p2_c = plot(angulo_haste,t3,(180/pi)*y3(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 3
         legend(angulo_haste,legenda1,legenda2,legenda3,'FontSize',v);   % Legenda do gráfico
         xlabel(angulo_haste,'Tempo [s]','FontSize',v)                   % Texto do eixo x
         ylabel(angulo_haste,strcat('Posição da Haste - [',char(176),']'),'FontSize',v) % Texto do eixo y
         grid(angulo_haste,'on');                                        % Habilita a grade
         hold(angulo_haste,'off');                                       % Libera o gráfico
         
         velocidade_carrinho = nexttile;                                 % Gráfico das velocidades lineares
         p3_a = plot(velocidade_carrinho,t1,x1(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 1
         hold(velocidade_carrinho,'on');                                 % Retém o gráfico
         p3_b = plot(velocidade_carrinho,t2,x2(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 2
         p3_c = plot(velocidade_carrinho,t3,x3(2,:),'LineWidth',w);      % Plot da velocidade linear do sistema 3
         legend(velocidade_carrinho,'Estocástico',legenda2,legenda3,'FontSize',v); % Legenda do gráfico
         xlabel(velocidade_carrinho,'Tempo [s]','FontSize',v);           % Texto do eixo x
         ylabel(velocidade_carrinho,'Velocidade do Carrinho [m/s]','FontSize',v) % Texto do eixo y
         grid(velocidade_carrinho,'on');                                 % Habilita a grade
         hold(velocidade_carrinho,'off'); 
         
         velocidade_angular = nexttile;                                  % Gráfico das velocidades
         p4_a = plot(velocidade_angular,t1,(180/pi)*x1(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 1
         hold(velocidade_angular,'on');                                      % Retém o gráfico
         p4_b = plot(velocidade_angular,t2,(180/pi)*x2(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 2
         p4_c = plot(velocidade_angular,t3,(180/pi)*x3(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 3
         legend(velocidade_angular,'Estocástico',legenda2,legenda3,'FontSize',v); % Legenda do gráfico
         xlabel(velocidade_angular,'Tempo [s]','FontSize',v);                % Texto do eixo x
         ylabel(velocidade_angular,'Velocidade da Haste [graus/s]','FontSize',v) % Texto do eixo y
         grid(velocidade_angular,'on');                                      % Habilita a grade
         hold(velocidade_angular,'off');                                     % Libera o gráfico
         
         fig2 = figure(2);                                               % Cria uma nova figura
         tl2 = tiledlayout(2,3);                                         % Layout da figura 2 linhas e 3 colunas
         tl2.TileSpacing = 'compact';                                    % Diminui o espaço entre os graficos
         tl2.Padding = 'compact';                                        % Diminui o espaço lateral dos gráficos
         
         sinal_controle = nexttile;                                      % Gráfico do sinal de controle
         p5_a = plot(sinal_controle,t1,u1','LineWidth',w);               % Plot do sinal de controle do sistema 1
         hold(sinal_controle,'on');                                      % Retém o gráfico
         p5_b = plot(sinal_controle,t1,u2','LineWidth',w);               % Plot do sinal de controle do sistema 2
         legend(sinal_controle,legenda2,legenda3,'FontSize',v); % Legenda do gráfico
         xlabel(sinal_controle,'Tempo [s]','FontSize',v);                % Texto do eixo x
         ylabel(sinal_controle,'Sinal de Controle','FontSize',v);        % Texto do eixo y
         grid(sinal_controle,'on');                                      % Habilita a grade
         hold(sinal_controle,'on');                                      % Libera o gráfico
         
         erro_estimativa1 = nexttile;                                    % Gráfico dos sinais de erro para a posição do carrinho
         p6_a = plot(erro_estimativa1,t1,e(1,:)','-r','LineWidth',w);    % Plot do erro de estimativa da posição do carrinho entre o sistema 1 e seus estimador  ,'-r'
         hold(erro_estimativa1,'on');                                    % Retém o gráfico % Libera o gráfico
         p6_b = plot(erro_estimativa1,t1,e2(1,:)',':b','LineWidth',w);   % Plot do erro de estimativa do ângulo da haste entre o sistema 1 e seus estimador  ,'-r' 1':b',
         p6_b.MarkerIndices = 1:500:length(e2)                           % Marcadores de 500 em 500
         legend(erro_estimativa1,'Posição do Carrinho','Ãngulo da Haste','FontSize',v); % Legenda do gráfico
         xlabel(erro_estimativa1,'Tempo [s]','FontSize',v);              % Texto do eixo x
         ylabel(erro_estimativa1,'Erro de Estimativa [saida 1]','FontSize',v); % Texto do eixo y
         grid(erro_estimativa1,'on');                                    % Libera o gráfico
         
         erro_estimativa2 = nexttile;                                    % Gráfico dos sinais de erro para o ângulo da haste
         p7_a = plot(erro_estimativa2,t2,e(2,:)','-r','LineWidth',w);    % Plot do erro de estimativa da posição do carrinho entre o sistema 1 e seu estimador  ,'-r'
         hold(erro_estimativa2,'on');                                    % Retém o gráfico
         p7_b = plot(erro_estimativa2,t2,e2(2,:)',':b','LineWidth',w);   % Plot do erro de estimativa do ângulo da haste entre o sistema 1 e seu estimador  ,'-r' 1':b',
         p7_b.MarkerIndices = 1:500:length(e2)                           % Marcadores de 500 em 500
         legend(erro_estimativa2,'Posição do Carrinho','Ãngulo da Haste','FontSize',v); % Legenda do gráfico
         xlabel(erro_estimativa2,'Tempo [s]','FontSize',v);              % Texto do eixo x
         ylabel(erro_estimativa2,'Erro de Estimativa [saida 2]','FontSize',v); % Texto do eixo y
         grid(erro_estimativa2,'on');                                    % Libera o gráfico
        
         covariancia1 = nexttile;                                        % Gráfico da covariância do filtro 1
         p8 = plot(covariancia1,t1,p','LineWidth',w);                    % Plot da covariância do filtro 1
         legend(covariancia1,'Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         xlabel(covariancia1,'Tempo [s]','FontSize',v);                  % Texto do eixo x
         ylabel(covariancia1,'Covarância do Erro - EKF','FontSize',v);   % Texto do eixo y
         grid(covariancia1,'on');                                        % Habilita a grade   
         
         covariancia2 = nexttile;                                        % Gráfico da covariância do filtro 2
         p9 = plot(covariancia2,t1,p2','LineWidth',w);                   % Plot da covariância do filtro 2
         legend(covariancia2,'Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         xlabel(covariancia2,'Tempo [s]','FontSize',v);                  % Texto do eixo x
         ylabel(covariancia2,'Covarância do Erro - KF','FontSize',v);    % Texto do eixo y
         grid(covariancia2,'on');                                        % Habilita a grade   
               
         ganho_kalman = nexttile;                                        % Gráfico da norma 2 das duas colunas do Ganho de Kalman
         L_1 = norm(L_k(:,1),2);                                         % Norma 2 da primeira coluna de L_k
         L_2 = norm(L_k(:,2),2);                                         % Norma 2 da segunda coluna de L_k
         L_3 = norm(L_k2(:,1),2);                                        % Norma 2 da primeira coluna de L_k2
         L_4 = norm(L_k2(:,2),2);                                        % Norma 2 da segunda coluna de L_k2
         for i = 3:2:length(L_k)                                         % Começa o loop da terceira coluna até o fim
             L_1 = [L_1 norm(L_k(:,i),2)];                               % Acumula a norma das colunas ímpares de L_k
             L_2 = [L_2 norm(L_k(:,i+1),2)];                             % Acumula a norma das colunas pares de L_k
             L_3 = [L_3 norm(L_k2(:,i),2)];                              % Acumula a norma das colunas ímpares de L_k2
             L_4 = [L_4 norm(L_k2(:,i+1),2)];                            % Acumula a norma das colunas pares de L_k2
         end
         p10_a = plot(ganho_kalman,t1,L_1','-p','LineWidth',w)           % Mostra L1 
         p10_a.MarkerIndices = 1:800:length(L_1)                         % Marcadores de 800 em 800
         hold(ganho_kalman,'on');                                        % Retém o gráfico
         p10_b = plot(ganho_kalman,t1,L_2','-x','LineWidth',w)           % Mostra L2
         p10_b.MarkerIndices = 1:600:length(L_1)                         % Marcadores de 600 em 600
         p10_c = plot(ganho_kalman,t3,L_3','-s','LineWidth',w)           % Mostra L3
         p10_c.MarkerIndices = 1:500:length(L_1)                         % Marcadores de 500 em 500
         p10_d = plot(ganho_kalman,t3,L_4','-d','LineWidth',w)           % Mostra L4
         p10_d.MarkerIndices = 1:700:length(L_1)                         % Marcadores de 700 em 700
         legend('Posição do Carrinho - EKF','Ângulo da Haste - EKF','Posição do Carrinho - KF','Ângulo da Haste - KF','FontSize',v) % Legenda
         ylabel('Ganho de Kalman - Norma 2','FontSize',v);               % Texto do eixo y
         xlabel('Tempo - s');                                            % Texto do eixo x
         grid(ganho_kalman,'on');                                        % Habilita a grade
         hold(ganho_kalman,'off');                                       % Libera o gráfico
     otherwise
        disp('Não escolheeu opção alguma')
end      
end

