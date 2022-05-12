function plotar_sistema(t1,x1,y1,u1,t2,x2,y2,u2,e,p,w,legenda1,legenda2)

% t1 -> Tempo de simulação para o primeiro sistema
% x1 -> Estados do primeiro sistema
% t2 -> Tempo de simulação para o segundo sistema
% x2 -> Estados do segundo sistema
% u  -> Sinal de controle
% e  -> Erro de estimativa
% p  -> Covariância do erro
% filtragem -> Se a plotagem será para o sistema filtrado
% dois -> Se a plotagem será para um ou dois sistemas
% legenda -> Legenda do(s) sistema(s)
% w  -> Espessura da linha do gráfico

dois = 0;             % Inicializa a quantidade de sistemas a serem plotados
filtragem = 0;        % Inicializa o indicador de filtragem

if length(x2) ~= 0    % Se existir um segundo sistema
    dois = 1;
end
if length(e) ~= 0 & length(p) ~= 0 % Se existir um filtro
    filtragem = 1;
end    

if filtragem          % Se o plot envolve os filtros KF ou EKF
    tiledlayout(3,2); % Layout da figura 3 linhas e 2 colunas
else                  % Se o plot não envolve os filtros KF ou EKF
    tiledlayout(2,2); % Layout da figura 2 linhas e 2 colunas
end
 
 nexttile;
 plot(t1,y1(1,:)','LineWidth',w);        % Plot da posição do carrinho do sistema 1 
 if dois
     hold on;                            % Retém o gráfico
     plot(t2,y2(1,:)','LineWidth',w);    % Plot da posição do carrinho do sistema 2
     legend(legenda1,legenda2)           % Legenda para dois sistemas
 end
 title('Posição do Carrinho vs. Tempo')  % Título do gráfico   
 xlabel('Tempo - s')                     % Texto do eixo x
 ylabel('Posição do carrinho - m')       % Texto do eixo y
 grid on;                                % Habilita a grade

 
 nexttile;
 plot(t1,(180/pi)*y1(2,:)','LineWidth',w);        % Plot do ângulo da haste do sistema 1
 if dois
     hold on;                                     % Retém o gráfico
     plot(t2,(180/pi)*y2(2,:)','LineWidth',w);    % Plot do ângulo da haste do sistema 2
     legend(legenda1,legenda2);                   % Legenda do gráfico
 end
 title('Ângulo da Haste vs. Tempo')  % Título do gráfico
 xlabel('Tempo - s')               % Texto do eixo x
 ylabel('Posição da Haste - º')    % Texto do eixo y
 grid on;                          % Habilita a grade
 if dois
     hold off;
 end
 
 nexttile;
 plot(t1,x1(2,:),'LineWidth',w);  % Plot da velocidade linear do sistema 1
 hold on;                         % Retém o gráfico
 plot(t1,x1(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 1
 vcarrinho = 'Velocidade do Carrinho - ';
 vhaste = 'Velocidade da Haste - ';
 if dois 
     plot(t2,x2(2,:),'LineWidth',w)  % Plot da velocidade linear do sistema 2
     plot(t2,x2(4,:),'LineWidth',w); % Plot da velocidade angular do sistema 2
     if ~filtragem
         legend(strcat(vcarrinho,char(160),legenda1),strcat(vhaste,char(160),legenda1),strcat(vcarrinho,char(160),legenda2),strcat(vhaste,char(160),legenda2)); % Legenda do gráfico
     else
         legend(strcat(vcarrinho,char(160),'Real'),strcat(vhaste,char(160),'Real'),strcat(vcarrinho,char(160),'Estimado'),strcat(vhaste,char(160),'Estimado')); % Legenda do gráfico
     end
 else
     legend('Velocidade Linear do Carrinho','Velocidade Angular da Haste'); % Legenda do gráfico
 end
 title('Velocidades Linear e Angular  vs. Tempo')  % Título do gráfico   
 xlabel('Tempo - s');                           % Texto do eixo x
 ylabel('Vel. do Carrinho m/s - Vel. da Haste rad/s') % Texto do eixo y
 grid on;                                       % Habilita a grade
 hold off;                                      % Libera o gráfico
 
 nexttile;
 plot(t1,u1','LineWidth',w);                    % Plot do sinal de controle 
 title('Sinal de Controle vs. Tempo')           % Título do gráfico
 if dois & length(u2) ~= 0
     hold on;
     plot(t2,u2','LineWidth',w); % Plot do sinal de controle 
     legend(legenda1,legenda2);  % Legenda do gráfico
 end
 xlabel('Tempo - s');                    % Texto do eixo x
 ylabel('Sinal de Controle');            % Texto do eixo y
 grid on;                                % Habilita a grade
 
 if filtragem
     nexttile;
     plot(t1,p','LineWidth',w);   % Plot do ângulo da haste do sistema 1  
     title('Covariância do Erro vs. Tempo')  % Título do gráfico
     legend('Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular'); % Legenda do gráfico
     xlabel('Tempo - s')          % Texto do eixo x
     ylabel('Covarância do Erro') % Texto do eixo y
     grid on;  
     
     nexttile;
     plot(t1,e','LineWidth',w);             % Plot do ângulo da haste do sistema 1  
     title('Erro de Estimativa entre as Saídas vs. Tempo') % Título do gráfico
     legend('Posição do Carrinho','Ângulo da Haste')       % Legenda do gráfico
     xlabel('Tempo - s')                    % Texto do eixo x
     ylabel('Erro de Estimativa vs. Tempo') % Texto do eixo y
     grid on;      
 end
end

