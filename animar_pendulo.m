function  animar_pendulo(qtd,y1,y2,y3,y4,y5,y6,s,l_haste,l_carrinho,h_carrinho,titulo1,titulo2,titulo3,titulo4,titulo5,titulo6,nome_arquivo,titulo,Q_lqr,R_lqr,sigma_quadrado_w_Q,sigma_quadrado_v_R,sigma_quadrado_w,sigma_quadrado_v,dt,ci1,ci2,ci3,limites_grafico,tempo_simulacao,gravar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que realiza a animação do pêndulo organizando o layout conforme o % 
% número de sistema a ser simulado, inclusive com vetores de tamanhos      %
% diferentes e grava a animação no formato .mp4                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Glossário de variáveis
    % qtd -> Quantidade de gráficos por figura
    % y1 -> Saída do sistema 1
    % y2 -> Saída do sistema 2
    % y3 -> Saída do sistema 3
    % y4 -> Saída do sistema 4
    % x -> Estado atual do sistema
    % s -> Escala da animação para aumento ou diminuição das dimensão do desenho
    % l_haste -> Comprimento da haste
    % l_carrinho -> Comprimento do carrinho
    % h_carrinho -> Altura do carrinho
    % titulo1 -> Título do sistema 1
    % titulo2 -> Título do sistema 2
    % titulo3 -> Título do sistema 3
    % titulo4 -> Título do sistema 4
    % nome_arquivo -> Nome do arquivo de vídeo que irá ser gravado
    % titulo -> Título mestre no topo dos gráficos
    % Q -> Matriz de poderação dos estados do LQR
    % R -> Matriz de ponderação do controle do LQR
    % var_Q -> Variância da incerteza do processo
    % var_R -> Variância da incerteza da medição
    % var_w -> Variância do ruído aditivo no processo
    % var_v -> Variância do ruído aditivo na medição
    % dt -> Período de amostragem
    % ci1 -> Condições iniciais para o sistema 1
    % ci2 -> Condições iniciais para o sistema 2
    % ci3 -> Condições iniciais para o sistema 3
    % gravar -> Variável lógica que indica se a animação será gravada (0 - não gravar e 1 - gravar)   
    
%% Obtenção das dimensões do modelo físico e aplicação da escala 
    lh = l_haste*s;            % Comprimento da haste em escala
    lc = l_carrinho*s;         % Comprimento do carrinho em escala
    hc = h_carrinho*s;         % Altura do carrinho em escala
    posicao_trilho = [limites_grafico(1,1),-0.6*hc,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.2*hc]; % Trilho  
    metade_trilho = 0.814/2;   % Medida da metado do trilho para configurar os limites

%% Configuração do título
    variancias = strcat('\bf',char(963),char(178),'_Q \rm=',char(160),num2str(sigma_quadrado_w_Q(1,1)),',',char(160),char(160),'\bf',char(963),char(178),'_R \rm=',char(160),num2str(sigma_quadrado_v_R(1,1)),',',char(160),char(160),'\bf',char(963),char(178),'_w\rm =',char(160),num2str(sigma_quadrado_w),',',char(160),char(160),'\bf',char(963),char(178),'_v \rm=',char(160),num2str(sigma_quadrado_v));
    tempo = strcat('\bf Tempo \rm=',char(160),num2str(tempo_simulacao,2),'s'); 
    periodo_amostragem = strcat('\bf T \rm=',char(160),num2str(dt));      % Período de Amostragem
    q_lqr = strcat('\bf Q \rm= diag',char(160),mat2str(diag(Q_lqr)',4));  % Diagonal principal da matriz Q do LQR
    r_lqr = strcat('\bfR \rm=',char(160),num2str(R_lqr));                 % Matriz R do LQR
    subtitulo = strcat(tempo,',',char(160),periodo_amostragem,',',char(160),char(160),q_lqr,',',char(160),char(160),r_lqr,',',char(160),char(160),variancias); % Subtítulo    

%% Criação e abertura do arquivo de vídeo   
    if gravar                  % Se for para gravar o vídeo
        video = VideoWriter(strcat('animacao\wb_0_02\',nome_arquivo,'.mp4'),'MPEG-4'); % Cria o arquivo de video
        video.FrameRate = 50;  % Quantidade de quadros po segundo (FPS)
        open(video);           % Abre o arquivo de video
    end
    
%% Quantidade de gráficos
f = figure;                               % Cria a conteiner para o gráfico

if qtd == 1
    f.Position = [250 150 800 500];       % Posição do gráfico para 1 animação
    tl1 = tiledlayout(1,1);               % Layout para um sistemas
    tl1.TileSpacing = 'compact';          % Diminui o espaço entre os gráficos
    tl1.Padding = 'compact';              % Diminui o espaço ao redor dos gráficos
    title(tl1,{titulo;subtitulo});        % Título do gráfico 1
    xlabel(tl1,'Posição Horizontal - m'); % Texto do eixo x
    ylabel(tl1,'Posição Vertical - m');   % Texto do eixo y
    ax1 = nexttile;       
     for i = 1:length(y1)             % Loop até o tamanho de y1       
      if isgraphics(ax1)              % Se o gráfico for válido
        cla(ax1);                     % Apaga o gráfico corrente
        hold(ax1,'on');               % Retém o gráfico corrente
        axis(ax1,limites_grafico);    % Define os eixos x e y  
        grid(ax1,'on');               % Habilita a grade
        plot(ax1,[y1(i,1), y1(i,1)+lh*sin(y1(i,2))], [0, lh*cos(y1(i,2))],'blue','LineWidth',2); % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
        rectangle(ax1,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
        rectangle(ax1,'Position',[y1(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
        hold(ax1,'off');              % Libera o gráfico corrente        
        drawnow();                    % Atualiza a figura com os dados anteriores
        if ~isempty(gcf) & gravar     % Verifica se a figura corrente é válida e se gravar = 1 (true)
          F(i) = getframe(f);         % Obtém a imagem da figura corrente     
          writeVideo(video,F(i));     % Armazena a imagem    
        end
      end
     end
elseif qtd == 2
    f.Position = [200 200 900 400];       % Posição do gráfico para 2 animações  
    tl1 = tiledlayout(1,2);               % Layout para um sistemas
    tl1.TileSpacing = 'compact';          % Diminui o espaço entre os gráficos
    tl1.Padding = 'compact';              % Diminui o espaço ao redor dos gráficos
    title(tl1,{titulo;subtitulo});        % Título do gráfico 1
    xlabel(tl1,'Posição Horizontal - m'); % Texto do eixo x
    ylabel(tl1,'Posição Vertical - m');   % Texto do eixo y
    ax1 = nexttile;                       % Gráfico 1
    title(strcat(titulo1,char(160),char(8658),char(160),'ci =',char(160),ci1)); % Título do gráfico 1
    ax2 = nexttile;                       % Gráfico 2
    title(strcat(titulo2,char(160),char(8658),char(160),'ci =',char(160),ci2)); % Título do gráfico 2
    t1 = linspace(0,tempo_simulacao,length(y1));     % Vetor para o relógio baseado em y2
    t2 = linspace(0,tempo_simulacao,length(y2));     % Vetor para o relógio baseado em y2
    
    if length(y1) >= length(y2)                      % No caso de 2 vetores de tamanhos diferentes especifica o maior e o menor para o de maior tamanho continuar executando enquanto o menor para 
        maior = y1;
        tmaior = t1;
        menor = y2;
        tmenor = t2;
    else
        maior = y2;
        tmaior = t2;
        menor = y1;
        tmenor = t1;
    end      
     for i = 1:length(maior)                          % Realiza o loop até o tamanho de maior (maior continua e menor para caso tenham tamanhos diferentes)
            if isgraphics(ax1)                        % Verifica se o handle ax1 é valido
                cla(ax1);                             % Limpa a figura anterior
                hold(ax1,'on');                       % Retém o gráfico corrente
                axis(ax1,limites_grafico);            % Define os eixos x e y
                grid(ax1,'on');                       % Habilita a grade
                plot(ax1,[maior(i,1), maior(i,1)+lh*sin(maior(i,2))], [0, lh*cos(maior(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax1,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980])                           % Carrinho
                rectangle(ax1,'Position',[maior(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos                
                plot(ax1, maior(i,1)+lh*sin(maior(i,2)), lh*cos(maior(i,2)), '.b', 'MarkerSize',40) % Manopla
                hold(ax1,'off');                      % Retém o gráfico corrente
                text(ax1,0.002,limites_grafico(1,3)+0.002,strcat('Tempo :',char(160),num2str(tmaior(1,i),2),'s')); % Relógio para gráficos com pequenos limites
                %text(ax1,metade_trilho+0.1,limites_grafico(1,3)+0.2,strcat('Tempo :',char(160),num2str(t1(1,i),2),'s')); % Relógio

            end
            if i <= length(menor) & isgraphics(ax2)   % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax2);                             % Limpa a figura anterior
                hold(ax2,'on');                       % Retém o gráfico corrente
                axis(ax2,limites_grafico);            % Define os eixos x e y             
                grid(ax2,'on');                       % Habilita a grade
                plot(ax2,[menor(i,1), menor(i,1)+lh*sin(menor(i,2))], [0, lh*cos(menor(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax2,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax2,'Position',[menor(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos                
                plot(ax2, menor(i,1)+lh*sin(menor(i,2)),  lh*cos(menor(i,2)), '.b', 'MarkerSize',40)
                text(ax2,0.002,limites_grafico(1,3)+0.002,strcat('Tempo :',char(160),num2str(tmenor(1,i),2),'s')); % Relógio para gráficos com pequenos limites
                %text(ax2,metade_trilho+0.1,limites_grafico(1,3)+0.2,strcat('Tempo :',char(160),num2str(t1(1,i),2),'s')); % Relógio
                hold(ax2,'off');                      % Libera o gráfico corrente
            end
            drawnow();                                % Atualiza a figura com os dados anteriores 
            if ~isempty(gcf) & gravar                 % Verifica se a figura é válida e se gravar = 1 (true)
                F(i) = getframe(gcf);                 % Obtém a imagem da figura corrente
                writeVideo(video,F(i));               % Armazena a imagem
            end
        end     
elseif qtd == 3
    f.Position = [100 200 1100 400];             % Posição do gráfico para 3 animações
    tl1 = tiledlayout(1,3);                      % Layout para um sistemas
    tl1.TileSpacing = 'compact';                 % Diminui o espaço entre os gráficos
    tl1.Padding = 'compact';                     % Diminui o espaço ao redor dos gráficos
    title(tl1,{titulo;subtitulo});               % Título do gráfico 1
    xlabel(tl1,'Posição Horizontal -\bf m \rm'); % Texto do eixo x
    ylabel(tl1,'Posição Vertical -\bf m \rm');   % Texto do eixo y
    ax1 = nexttile;                                  % Gráfico da medição do sistema não-linear
    title(strcat(titulo1,char(160),char(8658),char(160),'ci =',char(160),ci1));   % Título do grafico 1
    ax2 = nexttile;                                  % Gráfico da esperança não-linear
    title(strcat(titulo2,char(160),char(8658),char(160),'ci =',char(160),ci2));   % Título do gráfico 2
    ax3 = nexttile;                                  % Gráfico do EKF
    title(strcat(titulo3,char(160),char(8658),char(160),'ci =',char(160),ci3));   % Título do gráfico 2
    t1 = linspace(0,tempo_simulacao,length(y1));     % Vetor para o relógio baseado em y2
    t2 = linspace(0,tempo_simulacao,length(y2));     % Vetor para o relógio baseado em y2
    t3 = linspace(0,tempo_simulacao,length(y3));     % Vetor para o relógio baseado em y2
         for i = 1:length(y1)                        % Loop até o tamanho de maior
            if isgraphics(ax1)                       % Verifica de ax1 é válido
                cla(ax1);                            % Limpa a figura anterior
                hold(ax1,'on');                      % Retém o gráfico corrente
                axis(ax1,limites_grafico);           % Define os eixos x e y
                grid(ax1,'on');                      % Habilita a grade
                plot(ax1,[y1(i,1), y1(i,1)+lh*sin(y1(i,2))], [0, lh*cos(y1(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax1,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980])               % Carrinho
                rectangle(ax1,'Position',[y1(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta')                  % Trilhos
                plot(ax1, y1(i,1)+lh*sin(y1(i,2)),  lh*cos(y1(i,2)), '.b', 'MarkerSize',40)               % Manopla do pêndulo
                if (y1(i,1) == -metade_trilho + lc/2 | y1(i,1) < -metade_trilho + lc/2)   % Se o carrinho esbarra ou passa da barra lateral esquerda
                    xline(ax1,-metade_trilho,'--r',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal'); % Barra lateral esquerda vermelha                       
                    xline(ax1,metade_trilho,'--k',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal');             % Barra lateral direita azul              
                elseif (y1(i,1) == metade_trilho - lc/2 | y1(i,1) > metade_trilho - lc/2) % Se o carrinho  esbarra ou passa da barra lateral direita
                    xline(ax1,-metade_trilho,'--k',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal'); % Barra lateral esquerda azul                        
                    xline(ax1,metade_trilho,'--r',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal');             % Barra lateral direita vermelha               
                else                                                                      % Caso o carrinho esteja entre as duas barras tudo azul
                    xline(ax1,-metade_trilho,'--k',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal');                        
                    xline(ax1,metade_trilho,'--k',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal');
                end
                %text(ax1,0.002,limites_grafico(1,3)+0.002,strcat('Tempo :',char(160),num2str(t1(1,i),2),'s'));          % Relógio para gráficos com pequenos limites
                text(ax1,metade_trilho+0.1,limites_grafico(1,3)+0.2,strcat('Tempo :',char(160),num2str(t2(1,i),2),'s')); % Relógio
                hold(ax1,'off');                      % Retém o gráfico corrente
            end
            if  isgraphics(ax2)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax2);                             % Limpa a figura anterior
                hold(ax2,'on');                       % Retém o gráfico corrente
                axis(ax2,limites_grafico);            % Define os eixos x e y             
                grid(ax2,'on');                       % Habilita a grade
                plot(ax2,[y2(i,1), y2(i,1)+lh*sin(y2(i,2))], [0, lh*cos(y2(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax2,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980])               % Carrinho
                rectangle(ax2,'Position',[y2(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta')                  % Trilhos
                plot(ax2, y2(i,1)+lh*sin(y2(i,2)),  lh*cos(y2(i,2)), '.b', 'MarkerSize',40)               % Manopla do pêndulo
                if (y2(i,1) == -metade_trilho + lc/2 | y2(i,1) < -metade_trilho + lc/2)   % Se o carrinho esbarra ou passa da barra lateral esquerda
                    xline(ax2,-metade_trilho,'--r',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal'); % Barra lateral esquerda vermelha                       
                    xline(ax2,metade_trilho,'--k',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal'); % Barra lateral direita azul              
                elseif (y2(i,1) == metade_trilho - lc/2 | y2(i,1) > metade_trilho - lc/2) % Se o carrinho  esbarra ou passa da barra lateral direita
                    xline(ax2,-metade_trilho,'--k',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal'); % Barra lateral esquerda azul                        
                    xline(ax2,metade_trilho,'--r',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal'); % Barra lateral direita vermelha               
                else                                                                      % Caso o carrinho esteja entre as duas barras tudo azul
                    xline(ax2,-metade_trilho,'--k',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal');                        
                    xline(ax2,metade_trilho,'--k',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal');
                end
                %text(ax2,0.002,limites_grafico(1,3)+0.002,strcat('Tempo :',char(160),num2str(t1(1,i),2),'s'));          % Relógio
                text(ax2,metade_trilho+0.1,limites_grafico(1,3)+0.2,strcat('Tempo :',char(160),num2str(t2(1,i),2),'s')); % Relógio
                hold(ax2,'off');                      % Libera o gráfico corrente
            end
            if  isgraphics(ax3)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax3);                             % Limpa a figura anterior
                hold(ax3,'on');                       % Retém o gráfico corrente
                axis(ax3,limites_grafico);            % Define os eixos x e y             
                grid(ax3,'on');                       % Habilita a grade
                plot(ax3,[y3(i,1), y3(i,1)+lh*sin(y3(i,2))], [0, lh*cos(y3(i,2))],'blue','LineWidth',2);      % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax3,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980])                   % Carrinho
                rectangle(ax3,'Position',[y3(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta')                      % Trilhos
                plot(ax3, y3(i,1)+lh*sin(y3(i,2)),  lh*cos(y3(i,2)), '.b', 'MarkerSize',40)                   % Manopla do pêndulo
                if (y3(i,1) == -metade_trilho + lc/2 || y3(i,1) < -metade_trilho + lc/2)                      % Se o carrinho esbarra ou passa da barra lateral esquerda
                    xline(ax3,-metade_trilho,'--r',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal'); % Barra lateral esquerda vermelha                       
                    xline(ax3,metade_trilho,'--k',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal'); % Barra lateral direita azul              
                elseif (y3(i,1) == metade_trilho - lc/2 || y3(i,1) > metade_trilho - lc/2)                    % Se o carrinho  esbarra ou passa da barra lateral direita
                    xline(ax3,-metade_trilho,'--k',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal'); % Barra lateral esquerda azul                        
                    xline(ax3,metade_trilho,'--r',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal'); % Barra lateral direita vermelha               
                else                                                                                          % Caso o carrinho esteja entre as duas barras tudo azul
                    xline(ax3,-metade_trilho,'--k',strcat('-',num2str(metade_trilho),char(160),'Limite'),'LineWidth',1,'LabelHorizontalAlignment','left','LabelOrientation','horizontal');                        
                    xline(ax3,metade_trilho,'--k',strcat('Limite +',num2str(metade_trilho)),'LineWidth',1,'LabelHorizontalAlignment','right','LabelOrientation','horizontal');
                end
                %text(ax3,0.002,limites_grafico(1,3)+0.002,strcat('Tempo :',char(160),num2str(t1(1,i),2),'s'));          % Relógio
                text(ax3,metade_trilho+0.1,limites_grafico(1,3)+0.2,strcat('Tempo :',char(160),num2str(t3(1,i),2),'s')); % Relógio
                hold(ax3,'off');                      % Libera o gráfico corrente
            end            
            drawnow();                                % Atualiza a figura com os dados anteriores                               
            if ~isempty(gcf) & gravar                 % Verifica se a figura é válida e se gravar = 1 (true)
                F(i) = getframe(gcf);                 % Obtém a imagem da figura corrente
                writeVideo(video,F(i));               % Armazena a imagem
            end
        end

    
elseif qtd == 4
    f.Position = [250 50 850 620];        % Posição do gráfico para 4 animações
    tl1 = tiledlayout(2,2);               % Layout para um sistemas
    tl1.TileSpacing = 'compact';          % Diminui o espaço entre os gráficos
    tl1.Padding = 'compact';              % Diminui o espaço ao redor dos gráficos
    %title(tl1,{titulo;subtitulo});       % Título do gráfico 1
    xlabel(tl1,'Posição Horizontal - m'); % Texto do eixo x
    ylabel(tl1,'Posição Vertical - m');   % Texto do eixo y
    ax1 = nexttile;                                  % Gráfico da medição do sistema não-linear
    title(strcat(titulo1,char(160),char(8658),char(160),'ci =',char(160),ci1));   % Título do grafico 1
    ax2 = nexttile;                                  % Gráfico da esperança não-linear
    title(strcat(titulo2,char(160),char(8658),char(160),'ci =',char(160),ci2));   % Título do gráfico 2
    ax3 = nexttile;                                  % Gráfico do EKF
    title(strcat(titulo3,char(160),char(8658),char(160),'ci =',char(160),ci3));   % Título do gráfico 3
    ax4 = nexttile;                                  % Gráfico do EKF
    title(strcat(titulo4,char(160),char(8658),char(160),'ci =',char(160),ci3));   % Título do gráfico 4

    
     for i = 1:length(y1)                             % Loop até o tamanho de maior
            if isgraphics(ax1)                        % Verifica de ax1 é válido
                cla(ax1);                             % Limpa a figura anterior
                hold(ax1,'on');                       % Retém o gráfico corrente
                axis(ax1,limites_grafico);            % Define os eixos x e y
                grid(ax1,'on');                       % Habilita a grade
                plot(ax1,[y1(i,1), y1(i,1)+lh*sin(y1(i,2))], [0, lh*cos(y1(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax1,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax1,'Position',[y1(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax1,'off');                      % Retém o gráfico corrente
            end
            if  isgraphics(ax2)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax2);                             % Limpa a figura anterior
                hold(ax2,'on');                       % Retém o gráfico corrente
                axis(ax2,limites_grafico);            % Define os eixos x e y             
                grid(ax2,'on');                       % Habilita a grade
                plot(ax2,[y2(i,1), y2(i,1)+lh*sin(y2(i,2))], [0, lh*cos(y2(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax2,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax2,'Position',[y2(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax2,'off');                      % Libera o gráfico corrente
            end
            if  isgraphics(ax3)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax3);                             % Limpa a figura anterior
                hold(ax3,'on');                       % Retém o gráfico corrente
                axis(ax3,limites_grafico);            % Define os eixos x e y             
                grid(ax3,'on');                       % Habilita a grade
                plot(ax3,[y3(i,1), y3(i,1)+lh*sin(y3(i,2))], [0, lh*cos(y3(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax3,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax3,'Position',[y3(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax3,'off');                      % Libera o gráfico corrente
            end
            if  isgraphics(ax4)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax4);                             % Limpa a figura anterior
                hold(ax4,'on');                       % Retém o gráfico corrente
                axis(ax4,limites_grafico);            % Define os eixos x e y             
                grid(ax4,'on');                       % Habilita a grade
                plot(ax4,[y4(i,1), y4(i,1)+lh*sin(y4(i,2))], [0, lh*cos(y4(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax4,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax4,'Position',[y4(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax4,'off');                      % Libera o gráfico corrente
            end
            drawnow();                                % Atualiza a figura com os dados anteriores                               
            if ~isempty(gcf) & gravar                 % Verifica se a figura é válida e se gravar = 1 (true)
                F(i) = getframe(gcf);                 % Obtém a imagem da figura corrente
                writeVideo(video,F(i));               % Armazena a imagem
            end
        end

elseif qtd == 6
    f.Position = [100 50 1100 620];       % Posição do gráfico para 6 figuras
    tl1 = tiledlayout(2,3);               % Layout para um sistemas
    tl1.TileSpacing = 'compact';          % Diminui o espaço entre os gráficos
    tl1.Padding = 'compact';              % Diminui o espaço ao redor dos gráficos
    title(tl1,{titulo;subtitulo});        % Título do gráfico 1
    xlabel(tl1,'Posição Horizontal - m'); % Texto do eixo x
    ylabel(tl1,'Posição Vertical - m');   % Texto do eixo y
    ax1 = nexttile;                                  % Gráfico da medição do sistema não-linear
    title(strcat(titulo1,char(160),char(8658),char(160),'ci =',char(160),ci1));   % Título do grafico 1
    ax2 = nexttile;                                  % Gráfico da esperança não-linear
    title(strcat(titulo2,char(160),char(8658),char(160),'ci =',char(160),ci2));   % Título do gráfico 2
    ax3 = nexttile;                                  % Gráfico do EKF
    title(strcat(titulo3,char(160),char(8658),char(160),'ci =',char(160),ci3));   % Título do gráfico 3
    ax4 = nexttile;                                  % Gráfico do EKF
    title(strcat(titulo4,char(160),char(8658),char(160),'ci =',char(160),ci1));   % Título do gráfico 4
    ax5 = nexttile;                                  % Gráfico do EKF
    title(strcat(titulo5,char(160),char(8658),char(160),'ci =',char(160),ci2));   % Título do gráfico 3
    ax6 = nexttile;                                  % Gráfico do EKF
    title(strcat(titulo6,char(160),char(8658),char(160),'ci =',char(160),ci3));   % Título do gráfico 4

    
     for i = 1:length(y1)                             % Loop até o tamanho de maior
            if isgraphics(ax1)                        % Verifica de ax1 é válido
                cla(ax1);                             % Limpa a figura anterior
                hold(ax1,'on');                       % Retém o gráfico corrente
                axis(ax1,limites_grafico);            % Define os eixos x e y
                grid(ax1,'on');                       % Habilita a grade
                plot(ax1,[y1(i,1), y1(i,1)+lh*sin(y1(i,2))], [0, lh*cos(y1(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax1,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax1,'Position',[y1(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax1,'off');                      % Retém o gráfico corrente
            end
            if  isgraphics(ax2)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax2);                             % Limpa a figura anterior
                hold(ax2,'on');                       % Retém o gráfico corrente
                axis(ax2,limites_grafico);            % Define os eixos x e y             
                grid(ax2,'on');                       % Habilita a grade
                plot(ax2,[y2(i,1), y2(i,1)+lh*sin(y2(i,2))], [0, lh*cos(y2(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax2,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax2,'Position',[y2(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax2,'off');                      % Libera o gráfico corrente
            end
            if  isgraphics(ax3)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax3);                             % Limpa a figura anterior
                hold(ax3,'on');                       % Retém o gráfico corrente
                axis(ax3,limites_grafico);            % Define os eixos x e y             
                grid(ax3,'on');                       % Habilita a grade
                plot(ax3,[y3(i,1), y3(i,1)+lh*sin(y3(i,2))], [0, lh*cos(y3(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax3,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax3,'Position',[y3(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax3,'off');                      % Libera o gráfico corrente
            end
            if  isgraphics(ax4)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax4);                             % Limpa a figura anterior
                hold(ax4,'on');                       % Retém o gráfico corrente
                axis(ax4,limites_grafico);            % Define os eixos x e y             
                grid(ax4,'on');                       % Habilita a grade
                plot(ax4,[y4(i,1), y4(i,1)+lh*sin(y4(i,2))], [0, lh*cos(y4(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax4,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax4,'Position',[y4(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax4,'off');                      % Libera o gráfico corrente
            end
             if  isgraphics(ax5)                      % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax5);                             % Limpa a figura anterior
                hold(ax5,'on');                       % Retém o gráfico corrente
                axis(ax5,limites_grafico);            % Define os eixos x e y             
                grid(ax5,'on');                       % Habilita a grade
                plot(ax5,[y5(i,1), y5(i,1)+lh*sin(y5(i,2))], [0, lh*cos(y5(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax5,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax5,'Position',[y5(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax5,'off');                      % Libera o gráfico corrente
            end
            if  isgraphics(ax6)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax6);                             % Limpa a figura anterior
                hold(ax6,'on');                       % Retém o gráfico corrente
                axis(ax6,limites_grafico);            % Define os eixos x e y             
                grid(ax6,'on');                       % Habilita a grade
                plot(ax6,[y6(i,1), y6(i,1)+lh*sin(y6(i,2))], [0, lh*cos(y6(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax6,'Position',posicao_trilho,'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax6,'Position',[y6(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax6,'off');                      % Libera o gráfico corrente
            end
            
            drawnow();                                % Atualiza a figura com os dados anteriores                               
            if ~isempty(gcf) & gravar                 % Verifica se a figura é válida e se gravar = 1 (true)
                F(i) = getframe(gcf);                 % Obtém a imagem da figura corrente
                writeVideo(video,F(i));               % Armazena a imagem
            end
        end

end
    
%% Fecha o arquivo de vídeo
if gravar           % Se gravar = 1 (true)
    close(video);   % Grava e fecha o arquivo  
end
end
