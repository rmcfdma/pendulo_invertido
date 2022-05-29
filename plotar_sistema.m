function plotar_sistema(t1,x1,y1,u1,t2,x2,y2,u2,t3,x3,y3,u3,e,e2,p,p2,L_k,L_k2,w,legenda1,legenda2,legenda3,v,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que plota os gráficos para os sistemas com ou sem filtro.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Dicionário de variáveis
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

um = 0;               % Inicialização (não é para 1 gráfico)
dois = 0;             % Inicialização (não é para 2 gráfico)
tres = 0;             % Inicialização (não é para 3 gráfico)
filtragem = 0;        % Inicializa o indicador de filtragem (false)
controle1 = 0;        % Inicializa o indicador de controle do sistema 1


if isempty(x2) & isempty(x3) & isempty(u2)       % Se existir apenas um sistema
    um = 1;
elseif ~isempty(x2) & isempty(x3) & ~isempty(u2) % Se existir um segundo sistema
    dois = 1;
elseif ~isempty(x3) & ~isempty(u2)               % Se existir um terceiro sistema  
    tres = 1;
end
if length(L_k) ~= 0 & length(p) ~= 0             % Se existir um filtro
    filtragem = 1;
end    
if ~isempty(u1) & all(u1 ~= 0)                   % Se o controle estiver presente
    controle1 = 1;
end


    fig1 = figure(1);             % cria uma nova figura
    tl1 = tiledlayout(2,2);       % Layout da figura 2 linhas e 2 colunas
    tl1.TileSpacing = 'compact';  % Diminui o espaço entre os graficos
    tl1.Padding = 'compact';      % Diminui o espaço lateral dos gráficos
 
    nexttile;                                      % Gráfico da posição do carrinho
    p1 = plot(t1,y1(1,:)','LineWidth',w);          % Plot da posição do carrinho do sistema 1 
    xlabel('Tempo [s]','FontSize',v)               % Texto do eixo x
    ylabel('Posição do Carrinho [m]','FontSize',v) % Texto do eixo y
    grid on;                                       % Habilita a grade
    if dois & ~tres                                % Plotar dois sistemas
      hold on;                                     % Retém o gráfico
      p2 = plot(t2,y2(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 2
      legend(legenda1,legenda2,'FontSize',v)       % Legenda para dois sistemas
    elseif tres                                    % Plotar três sistema
      hold on;                                     % Retém o gráfico
      p3 = plot(t2,y2(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 2,':*g'
      %  p3.MarkerIndices = 1:250:length(y2);      % 
      p4 = plot(t3,y3(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 2,'-r',
      legend(legenda1,legenda2,legenda3,'FontSize',v)  % Legenda para dois sistemas    
    end
      
    nexttile;                                      % Gráfico da posição angular da haste
    p5 = plot(t1,(180/pi)*y1(2,:)','LineWidth',w); % Plot do ângulo da haste do sistema 1
    xlabel('Tempo [s]','FontSize',v)               % Texto do eixo x
    ylabel(strcat('Posição da Haste - [',char(176),']'),'FontSize',v) % Texto do eixo y
    grid on;                                       % Habilita a grade  
    if dois & ~tres                                % Plotar dois sistemas
      hold on;                                     % Retém o gráfico
      p6 = plot(t2,(180/pi)*y2(2,:)','LineWidth',w);  % Plot do ângulo da haste do sistema 2
      legend(legenda1,legenda2,'FontSize',v);      % Legenda do gráfico
    elseif tres                                    % Se forem 3 sistemas
      hold on;                                     % Retém o gráfico
      p7 = plot(t2,(180/pi)*y2(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 2,':*g'
      %  p7.MarkerIndices = 1:250:length(y2)           % Quantidade de marcadores no gráfico
      p8 = plot(t3,(180/pi)*y3(2,:)','LineWidth',w);   % Plot do ângulo da haste do sistema 2,'-r'
      legend(legenda1,legenda2,legenda3,'FontSize',v); % Legenda do gráfico
    end
     
    nexttile;                                      % Gráfico do sinal de controle
    plot(t1,x1(2,:),'LineWidth',w);                % Plot da velocidade linear do sistema 1
    xlabel('Tempo [s]','FontSize',v);              % Texto do eixo x
    grid on;                                       % Habilita a grade
    if ~dois & controle1 & ~tres                   % Se for somente um sistema e hover controle
        hold on;                                   % Retém o gráfico
        plot(t1,(180/pi)*x1(4,:),'LineWidth',w);   % Plot da velocidade angular do sistema 1 
        legend('Velocidade Linear do Carrinho','Velocidade Angular da Haste','FontSize',v); % Legenda do gráfico
        ylabel('Vel. Carrinho [m/s] - Haste [graus/s]','FontSize',v) % Texto do eixo y
        hold off;                                            % Libera o gráfico    
    elseif ~dois & ~controle1 & ~tres                        % Se for somente um sistema e não houver controle       
        ylabel('Velocidade do Carrinho [m/s]','FontSize',v)  % Texto do eixo y               
        nexttile;                                            % Gráfico da velocidade do carrinho
        plot(t1,(180/pi)*x1(4,:),'LineWidth',w);             % Plot da velocidade angular do sistema 1 
        ylabel('Velocidade da Haste [graus/s]','FontSize',v) % Texto do eixo y
        xlabel('Tempo [s]','FontSize',v);                    % Texto do eixo x
        grid on;       
    elseif dois & ~filtragem & ~tres               % Se for dois sistemas e não tratar-se de filtragem
        hold on;                                   % Retém o gráfico
        plot(t1,(180/pi)*x1(4,:),'LineWidth',w);   % Plot da velocidade angular do sistema 1 
        plot(t2,x2(2,:),'LineWidth',w)             % Plot da velocidade linear do sistema 2
        plot(t2,(180/pi)*x2(4,:),'LineWidth',w);   % Plot da velocidade angular do sistema 2
        vcarrinho = 'Velocidade do Carrinho -';    % Texto
        vhaste = 'Velocidade da Haste -';          % Texto
        legend(strcat(vcarrinho,char(160),legenda1),strcat(vhaste,char(160),legenda1),strcat(vcarrinho,char(160),legenda2),strcat(vhaste,char(160),legenda2),'FontSize',v); % Legenda do gráfico
        ylabel('Vel. Carrinho [m/s] - Haste [graus/s]','FontSize',v) % Texto do eixo y
        hold off;                                  % Libera o gráfico
     elseif dois & filtragem & ~tres               % Se forem dois sistemas e houver filtragem
        hold on;                                   % Retém o grpafico
        plot(t2,x2(2,:),'LineWidth',w)             % Plot da velocidade linear do sistema 2 
        legend('Esperança','Estimada','FontSize',v)% Legenda
        ylabel('Velocidade do Carrinho [m/s]','FontSize',v);  % Label do eixo y
        hold off;                                  % Libera o gráfico
        
        nexttile;                                  % Gráfico das velocidade linear e angular do carrinho e da haste
        plot(t1,(180/pi)*x1(4,:),'LineWidth',w);   % Plot da velocidade linear do sistema 1
        hold on                                    % retém o gráfico
        plot(t2,(180/pi)*x2(4,:),'LineWidth',w);   % Plot da velocidade angular do sistema 2
        legend('Esperança','Estimada','FontSize',v)% Legenda do gráfico
        xlabel('Tempo [s]','FontSize',v);          % Texto do eixo x
        grid on;                                   % Libera o gráfico   
        ylabel('Velocidade da Haste [graus/s]','FontSize',v); % Texto do eixo y
    elseif tres                                    % Se forem três sistema
        hold on;                                   % Retém o grpafico
        plot(t2,x2(2,:),'LineWidth',w)             % Plot da velocidade linear do sistema 2 
        plot(t3,x3(2,:),'LineWidth',w)             % Plot da velocidade linear do sistema 2 
        legend('Esperança','EKF','KF','FontSize',v)           % Legenda
        ylabel('Velocidade do Carrinho [m/s]','FontSize',v);  % Label do eixo y
        hold off;                                  % Libera o gráfico
        
        nexttile;                                     % Grafico das velocidade 
        plot(t1,(180/pi)*x1(4,:),'LineWidth',w);      % Plot da velocidade angular da haste do sistema 1
        hold on                                       % retém o gráfico
        plot(t2,(180/pi)*x2(4,:),'LineWidth',w);      % Plot da velocidade angular da haste do sistema 2
        plot(t3,(180/pi)*x3(4,:),'LineWidth',w);      % Plot da velocidade angular da haste do sistema 3
        legend('Esperança','EKF','KF','FontSize',v);  % Legenda do gráfico
        xlabel('Tempo [s]','FontSize',v);             % Texto do eixo x
        grid on;                                      % Libera o gráfico   
        ylabel('Velocidade da Haste [graus/s]','FontSize',v); % Texto do eixo y    
    end
        
     if um                                        % Se for apenas um sistema e houver controle
        nexttile;                                 % Gráfico do sinal de controle para um sistema
        plot(t1,u1','LineWidth',w);               % Plot do sinal de controle 
        xlabel('Tempo [s]','FontSize',v);         % Texto do eixo x
        ylabel('Sinal de Controle','FontSize',v); % Texto do eixo y
        grid on;                                  % Habilita a grade
        hold on;                                  % Retém o gráfico
     elseif dois & ~filtragem & ~tres             % Se forem apenas dois sistemas e não houver filtragem
        nexttile;                                 % Gráfico do sinal de controle para dois sistemas
        plot(t1,u1','LineWidth',w);               % Plot do sinal de controle 
        hold on;                                  % Retém o gráficoe
        plot(t2,u2','LineWidth',w);               % Plot do sinal de controle 
        legend(legenda1,legenda2,'FontSize',v);   % Legenda do gráfico    
        xlabel('Tempo [s]','FontSize',v);         % Texto do eixo x
        ylabel('Sinal de Controle','FontSize',v); % Texto do eixo y
        grid on;                                  % Habilita a grade
        hold on;                                  % Retém o gráfico
     elseif dois & filtragem & ~tres              % Se forem dois sistemas e houver filtragem
         fig2 = figure(2);                        % Cria uma segunda figura
         tl2 = tiledlayout(2,3);                  % Layout da figura 2 linhas e 2 colunas
         tl2.TileSpacing = 'compact';             % Diminui o espaço entre os graficos
         tl2.Padding = 'compact';                 % Diminui o espaço lateral dos graficos
         
         nexttile;                                % Gráfico para um sistema com filtro
         plot(t1,u1','LineWidth',w);              % Plot do sinal de controle com filtro
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Sinal de Controle','FontSize',v);% Texto do eixo y
         grid on;                                 % Habilita a grade
                         
         nexttile;                                % Gráfico do erro de estimativa sistema e filtro 1
         plot(t1,e(1,:)','LineWidth',w);          % Plot do erro de estimativa entre o sistema e o filtro 1  
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Erro de Estimativa [Saída 1]','FontSize',v); % Texto do eixo y
         grid on;                                 % Habilita a grade
         
         nexttile;                                % Gráfico do erro de estimativa sistema e filtro 2
         plot(t1,e(2,:)','LineWidth',w);          % Plot do erro de estimativa entre o sistema e o filtro 2  
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Erro de Estimativa [Saída 2]','FontSize',v); % Texto do eixo y
         grid on;                                 % Habilita a grade
         
         nexttile;                                % Gráfico da covariência
         plot(t1,p','LineWidth',w);               % Plot do ângulo da haste do sistema 1  
         legend('Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Covarância do Erro','FontSize',v); % Texto do eixo y
         grid on;                                 % Habilita a grade
         
         L_1 = L_k(:,(1:2:end));                  % Seleciona as colunas ímpares referente a saída 1
         L_2 = L_k(:,(2:2:end));                  % Seleciona as colunas pares referente a saída 2
         
         nexttile;                                % Gráfico do Ganho de Kalman referente à posição do carrinho
         plot(t1,L_1','LineWidth',w)              % Plota o Ganho de Kalman referente a saída 1
         legend('Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         ylabel('Ganho de Kalman [Saída 1]','FontSize',v); % Label do eixo vertical
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         grid on;                                 % Habilita a grade
         
         nexttile;                                % Gráfico do Ganho de Kalman referente ao ângulo da haste
         plot(t1,L_2','LineWidth',w)              % Plota o do Ganho de Kalman referente a saída 2
         legend('Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         ylabel('Ganho de Kalman [Saída 2]','FontSize',v); % Texto do eixo y
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         grid on;                                 % Habilita a grade
         
    elseif tres                                   % Se forem três sistemas
         fig2 = figure(2);                        % Cria uma segunda figura
         tl2 = tiledlayout(2,3);                  % Layout da figura 2 linhas e 3 colunas
         tl2.TileSpacing = 'compact';             % Diminui o espaço entre os graficos
         tl2.Padding = 'compact';                 % Diminui o espaço lateral dos graficos
         
         nexttile;                                % Gráfico dois sinais de controle referentes aos sistema e seus filtros
         plot(t2,u2','LineWidth',w);              % Plot do sinal de controle do sistema 1 ,'-r'
         hold on;                                 % Retém o gráfico
         pu3 = plot(t3,u3','LineWidth',w);        % Plot do sinal de controle do sistema 2 ,':b'
         % pu3.MarkerIndices = 1:50:length(u3)
         legend('EKF','KF','FontSize',v)          % Legenda
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Sinal de Controle','FontSize',v);% Texto do eixo y
         grid on;                                 % Habilita a grade
         
         nexttile;                                % Gráfico dos sinais de erro 
         plot(t1,e(1,:)','-r','LineWidth',w);     % Plot do erro de estimativa da posição do carrinho entre o sistema 1 e seus estimador  ,'-r'
         hold on                                  % Retém o Gráfico
         pe21 = plot(t1,e2(1,:)',':b','LineWidth',w);  % Plot do erro de estimativa da posição do carrinho entre o sistema 2 e seus estimador,':b' 
         pe21.MarkerIndices = 1:500:length(e2)
         legend('EKF','KF','FontSize',v)          % Legenda
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Erro de Estimativa [Saída 1]','FontSize',v); % Texto do eixo y
         grid on;
         
         nexttile;                                % Grafico dos erros de estimativa 
         plot(t1,e(2,:)','-r','LineWidth',w);          % Plot do erro de estimativa do ângulo da haste entre o sistema 1 e seus estimador  ,'-r' ,'-r'
         hold on;                                 % Retém o gráfico
         pe22 = plot(t1,e2(2,:)',':b','LineWidth',w);  % Plot do erro de estimativa do ângulo da haste entre o sistema 1 e seus estimador  ,'-r' 1':b',
         pe22.MarkerIndices = 1:500:length(e2)
         legend('EKF','KF','FontSize',v);         % Legenda   
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Erro de Estimativa [Saída 2]','FontSize',v); % Texto do eixo y
         grid on;                                 % Habilita a grade
         
         nexttile;                                % Gráfico da covariância do filtro 1
         plot(t1,p','LineWidth',w);               % Plot da covariância do filtro 1  
         legend('Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Covarância do Erro - EKF','FontSize',v); % Texto do eixo y
         grid on;                                 % Habilita a grade
         
         nexttile;                                % Gráfico da covariância do filtro 2 
         plot(t3,p2','LineWidth',w);              % Plot da covariância do filtro 2 
         legend('Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         xlabel('Tempo [s]','FontSize',v);        % Texto do eixo x
         ylabel('Covarância do Erro - KF','FontSize',v); % Texto do eixo y
         grid on;
                  
         nexttile;                                % Gráfico da norma 2 das duas colunas do Ganho de Kalman
         L_1 = norm(L_k(:,1),2);                  % Norma 2 da primeira coluna de L_k
         L_2 = norm(L_k(:,2),2);                  % Norma 2 da segunda coluna de L_k
         L_3 = norm(L_k2(:,1),2);                 % Norma 2 da primeira coluna de L_k2
         L_4 = norm(L_k2(:,2),2);                 % Norma 2 da segunda coluna de L_k2
         for i = 3:2:length(L_k)                  % Começa o loop da terceira coluna até o fim
           L_1 = [L_1 norm(L_k(:,i),2)];          % Acumula a norma das colunas ímpares de L_k
           L_2 = [L_2 norm(L_k(:,i+1),2)];        % Acumula a norma das colunas pares de L_k
           L_3 = [L_3 norm(L_k2(:,i),2)];         % Acumula a norma das colunas ímpares de L_k2
           L_4 = [L_4 norm(L_k2(:,i+1),2)];       % Acumula a norma das colunas pares de L_k2
         end
         p15 = plot(t1,L_1','-p','LineWidth',w)   % Mostra L1 
         p15.MarkerIndices = 1:800:length(L_1)
         hold on                                  % Retém o gráfico
         p16 = plot(t1,L_2','-x','LineWidth',w)   % Mostra L2
         p16.MarkerIndices = 1:600:length(L_1)
         p17 = plot(t3,L_3','-s','LineWidth',w)   % Mostra L3
         p17.MarkerIndices = 1:500:length(L_1)
         p18 = plot(t3,L_4','-d','LineWidth',w)   % Mostra L4
         p18.MarkerIndices = 1:700:length(L_1)
         hold off                                 % Libera o gráfico
         legend('Posição do Carrinho - EKF','Ângulo da Haste - EKF','Posição do Carrinho - KF','Ângulo da Haste - KF','FontSize',v) % Legenda
         ylabel('Ganho de Kalman - Norma 2','FontSize',v); % Texto do eixo y
         xlabel('Tempo - s');                     % Texto do eixo x
         grid on;                                 % Habilita a grade         
     end          

end

