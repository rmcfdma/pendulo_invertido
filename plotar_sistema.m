function plotar_sistema(t1,x1,y1,u1,t2,x2,y2,u2,e,p,L_k,w,legenda1,legenda2,v,m)

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
% v  -> Tamanho da fonte da legenda e eixo y 
% m  -> Tamanho da fonte do título 

dois = 0;             % Inicializa a quantidade de sistemas a serem plotados
filtragem = 0;        % Inicializa o indicador de filtragem

if length(x2) ~= 0    % Se existir um segundo sistema
    dois = 1;
end
if length(L_k) ~= 0 & length(p) ~= 0 % Se existir um filtro
    filtragem = 1;
end    

    fig1 = figure(1);
    tl1 = tiledlayout(2,2);       % Layout da figura 2 linhas e 2 colunas
    tl1.TileSpacing = 'compact';  % Diminui o espaço entre os graficos
    tl1.Padding = 'compact';      % Diminui o espaço lateral dos gráficos
 
    nexttile;
    plot(t1,y1(1,:)','LineWidth',w);               % Plot da posição do carrinho do sistema 1 
    title('(a)','FontSize',m)                      % Título do gráfico   
    xlabel('Tempo - s')                            % Texto do eixo x
    ylabel('Posição do carrinho - m','FontSize',v) % Texto do eixo y
    grid on;                                       % Habilita a grade
    if dois                                        % Plotar dois sistemas
      hold on;                                     % Retém o gráfico
      plot(t2,y2(1,:)','LineWidth',w);             % Plot da posição do carrinho do sistema 2
      legend(legenda1,legenda2,'FontSize',v)       % Legenda para dois sistemas

    end
    
  
    nexttile;
    plot(t1,(180/pi)*y1(2,:)','LineWidth',w);      % Plot do ângulo da haste do sistema 1
    if dois                                        % Plotar dois sistemas
      hold on;                                     % Retém o gráfico
      plot(t2,(180/pi)*y2(2,:)','LineWidth',w);    % Plot do ângulo da haste do sistema 2
      legend(legenda1,legenda2,'FontSize',v);      % Legenda do gráfico
    end
    title('(b)','FontSize',m)                      % Título do gráfico
    xlabel('Tempo - s')                            % Texto do eixo x
    ylabel('Posição da Haste - grau','FontSize',v) % Texto do eixo y
    grid on;                                       % Habilita a grade  

 
    nexttile;
    plot(t1,x1(2,:),'LineWidth',w);                % Plot da velocidade linear do sistema 1
    xlabel('Tempo - s');                           % Texto do eixo x
    grid on;                                       % Habilita a grade
    title('(c)','FontSize',m)                      % Título do gráfico       
    if ~dois                                       % Se for somente um sistema
        hold on;                                   % Retém o gráfico
        plot(t1,x1(4,:),'LineWidth',w);            % Plot da velocidade angular do sistema 1 
        legend('Velocidade Linear do Carrinho','Velocidade Angular da Haste','FontSize',v); % Legenda do gráfico
        ylabel('Vel. do Carrinho m/s - Vel. da Haste rad/s','FontSize',v) % Texto do eixo y
        hold off;                                  % Libera o gráfico
    elseif dois & ~filtragem                       % Se for dois sistemas e não tratar-se de filtragem
        hold on;                                   % Retém o gráfico
        plot(t1,x1(4,:),'LineWidth',w);            % Plot da velocidade angular do sistema 1 
        plot(t2,x2(2,:),'LineWidth',w)             % Plot da velocidade linear do sistema 2
        plot(t2,x2(4,:),'LineWidth',w);            % Plot da velocidade angular do sistema 2
        vcarrinho = 'Velocidade do Carrinho - ';   
        vhaste = 'Velocidade da Haste - ';   
        legend(strcat(vcarrinho,char(160),legenda1),strcat(vhaste,char(160),legenda1),strcat(vcarrinho,char(160),legenda2),strcat(vhaste,char(160),legenda2),'FontSize',v); % Legenda do gráfico
        ylabel('Vel. do Carrinho m/s - Vel. da Haste rad/s','FontSize',v) % Texto do eixo y
        hold off;                                  % Libera o gráfico
     elseif dois & filtragem                       % Se for dois sistemas e tratar-se de filtragem
        hold on;                                   % Retém o grpafico
        plot(t2,x2(2,:),'LineWidth',w)             % Plot da velocidade linear do sistema 2 
        legend('Medida','Estimada','FontSize',v)            % Legenda
        ylabel('Velocidade do Carrinho m/s','FontSize',v);  % Label do eixo y
        hold off;                                  % Libera o gráfico
        
        nexttile;
        plot(t1,x1(4,:),'LineWidth',w);  % Plot da velocidade linear do sistema 1
        hold on                          % retém o gráfico
        plot(t2,x2(4,:),'LineWidth',w);  % Plot da velocidade angular do sistema 2
        legend('Medida','Estimada','FontSize',v) % Legenda do gráfico
        xlabel('Tempo - s');             % Texto do eixo x
        grid on;                         % Libera o gráfico
        title('(d)','FontSize',m)        % Título do gráfico         
        ylabel('Velocidade da Haste rad/s','FontSize',v); % Texto do eixo y
    end
        
     if ~dois          
        nexttile;
        plot(t1,u1','LineWidth',w);               % Plot do sinal de controle 
        title('(d)','FontSize',m)                 % Título do gráfico
        xlabel('Tempo - s');                      % Texto do eixo x
        ylabel('Sinal de Controle','FontSize',v); % Texto do eixo y
        grid on;                                  % Habilita a grade
        hold on;                                  % Retém o gráfico
     elseif dois & ~filtragem                     % Se for dois sistemas e não tratar-se de filtragem
        nexttile;
        plot(t1,u1','LineWidth',w);               % Plot do sinal de controle 
        hold on;                                  % Rtém o gráficoe
        plot(t2,u2','LineWidth',w);               % Plot do sinal de controle 
        legend(legenda1,legenda2,'FontSize',v);   % Legenda do gráfico    
        title('(d)','FontSize',m)                 % Título do gráfico
        xlabel('Tempo - s');                      % Texto do eixo x
        ylabel('Sinal de Controle','FontSize',v); % Texto do eixo y
        grid on;                                  % Habilita a grade
        hold on;                                  % Retém o gráfico
     elseif dois & filtragem                      % Se for dois sistemas e tratar-se de filtragem
         fig2 = figure(2);                        % Cria uma segunda figura
         tl2 = tiledlayout(2,2);                  % Layout da figura 2 linhas e 2 colunas
         tl2.TileSpacing = 'compact';             % Diminui o espaço entre os graficos
         tl2.Padding = 'compact';                 % Diminui o espaço lateral dos graficos
         
         nexttile;
         plot(t1,u1','LineWidth',w);               % Plot do sinal de controle 
         title('(e)','FontSize',m)                 % Título do gráfico
         xlabel('Tempo - s');                      % Texto do eixo x
         ylabel('Sinal de Controle','FontSize',v); % Texto do eixo y
         grid on;                                  % Habilita a grade
         
         nexttile;
         plot(t1,p','LineWidth',w);                % Plot do ângulo da haste do sistema 1  
         title('(f)','FontSize',m)                 % Título do gráfico
         legend('Posição do Carrinho','Velocidade do Carrinho','Posição Angular','Velocidade Angular','FontSize',v); % Legenda do gráfico
         xlabel('Tempo - s')                       % Texto do eixo x
         ylabel('Covarância do Erro','FontSize',v) % Texto do eixo y
         grid on;                                  % Habilita a grade
         
         nexttile;
         plot(t1,e','LineWidth',w);                % Plot do ângulo da haste do sistema 1  
         title('(g)','FontSize',m)                 % Título do gráfico
         legend('Posição do Carrinho','Ângulo da Haste','FontSize',v) % Legenda do gráfico
         xlabel('Tempo - s')                       % Texto do eixo x
         ylabel('Erro de Estimativa','FontSize',v) % Texto do eixo y
         grid on;
         
         nexttile;
         L_1 = norm(L_k(:,1),2);
         L_2 = norm(L_k(:,2),2);
         for i = 3:2:length(L_k)
           L_1 = [L_1 norm(L_k(:,i),2)];
           L_2 = [L_2 norm(L_k(:,i+1),2)];
         end
         plot(t1,L_1','LineWidth',w)
         hold on
         plot(t1,L_2','LineWidth',w)
         hold off
         title('(h)','FontSize',m);
         ylabel('Ganho de Kalman','FontSize',v);
         xlabel('Tempo - s');
         grid on;
         legend('Referente à Saída 1 - Posição do Carrinho','Referente à Saída 2 - Ângular da haste','FontSize',v) 
     end          
end
