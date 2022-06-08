ffunction animar_pendulo(y1,y2,y3,y4,s,l_haste,l_carrinho,h_carrinho,titulo1,titulo2,titulo3,titulo4,nome_arquivo,titulo,Q,R,var_Q,var_R,var_w,var_v,dt,ci1,ci2,ci3,limites_grafico,gravar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que realiza a animação do pêndulo organizando o layout conforme o % 
% número de sistema a ser simulado, inclusive com vetores de tamanhos      %
% diferentes e grava a animação no formato .mp4                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Dicionário de variáveis
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
    
%% Configuração do título
    variancias = strcat(char(963),char(178),'_Q =',char(160),num2str(var_Q(1,1)),',',char(160),char(160),char(963),char(178),'_R =',char(160),num2str(var_R(1,1)),',',char(160),char(160),char(963),char(178),'_w =',char(160),num2str(var_w),',',char(160),char(160),char(963),char(178),'_v =',char(160),num2str(var_v));
    periodo_amostragem = strcat('T =',char(160),dt);           % Período de Amostragem
    Q_lqr = strcat('Q = diag',char(160),mat2str(diag(Q)',4));  % Diagonal principal da matriz Q do LQR
    R_lqr = strcat('R =',char(160),num2str(R));                % Matriz R do LQR
    subtitulo = strcat(periodo_amostragem,',',char(160),char(160),Q_lqr,',',char(160),char(160),R_lqr,',',char(160),char(160),variancias); % Subtítulo    

%% Criação e abertura do arquivo de vídeo   
    if gravar                  % Se for para gravar o vídeo
        video = VideoWriter(strcat('animacao\',nome_arquivo,'.mp4'),'MPEG-4'); % Cria o arquivo de video
        open(video);           % Abre o arquivo de video
    end
    
%% Trilho do pêndulo

 
    
%% Animação para 1,2 ou 3 sistemas    
if isempty(y2) & isempty(y3)          % Se a animação for para 1 sistema
    f = figure;                       % Cria a conteiner para o gráfico
    f.Position = [300 100 700 500];   % Tamanho e posição configurados para 2 gráficos [x1,y1,x2,y2]
    tl1 = tiledlayout(1,1);           % Layout para um sistemas
    tl1.TileSpacing = 'compact';      % Diminui o espaço entre os gráficos
    tl1.Padding = 'compact';          % Diminui o espaço ao redor dos gráficos
    title(tl1,{titulo;subtitulo});        % Título do gráfico 1
    xlabel(tl1,'Posição Horizontal - m'); % Texto do eixo x
    ylabel(tl1,'Posição Vertical - m');   % Texto do eixo y
    ax = nexttile;                        % Gráfico 1
    title(ax,strcat('ci =',char(160),ci1)); % Título do gráfico 1
   
    for i = 1:length(y1)             % Loop até o tamanho de y1       
      if isgraphics(ax)              % Se o gráfico for válido
        cla(ax);                     % Apaga o gráfico corrente
        hold(ax,'on');               % Retém o gráfico corrente
        axis(ax,limites_grafico);    % Define os eixos x e y  
        grid(ax,'on');               % Habilita a grade
        plot(ax,[y1(i,1), y1(i,1)+lh*sin(y1(i,2))], [0, lh*cos(y1(i,2))],'blue','LineWidth',2); % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
        rectangle(ax,'Position',[limites_grafico(1,1),-0.2,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.1],'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
        rectangle(ax,'Position',[y1(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
        hold(ax,'off');              % Libera o gráfico corrente
        drawnow();                   % Atualiza a figura com os dados anteriores
        if ~isempty(gcf) & gravar    % Verifica se a figura corrente é válida e se gravar = 1 (true)
          F(i) = getframe(f);        % Obtém a imagem da figura corrente     
          writeVideo(video,F(i));    % Armazena a imagem    
        end
      end
    end  
elseif ~isempty(y1) & ~isempty(y3) & isempty(y2)& isempty(y4)   % Se a animação for para dois sistemas
    f = figure;                          % Cria o conteiner para os gráficos
    f.Position = [200 200 900 400];      % Tamanho e posição configurados para 2 gráficos [x1,y1,x2,y2]
    tl2 = tiledlayout(1,2);              % Layout para dois sistemas com uma linha e duas colunas
    tl2.TileSpacing = 'compact';         % Diminui o espaço entre os gráficos
    tl2.Padding = 'compact';             % Diminui o espaço ao redor dos gráficos
    title(tl2,{titulo;subtitulo});       % Título mestre da figura
    xlabel(tl2,'Posição Horizontal - m') % Texto do eixo x
    ylabel(tl2,'Posição Vertical - m')   % Texto do eixo y
    ax1 = nexttile;                      % Gráfico 1
    title(strcat(titulo1,char(160),char(8658),char(160),'ci =',char(160),ci1)); % Título do gráfico 1
    ax2 = nexttile;                      % Gráfico 2
    title(strcat(titulo3,char(160),char(8658),char(160),'ci =',char(160),ci2)); % Título do gráfico 2
    
    if length(y1) >= length(y3)          % No caso de 2 vetores de tamanhos diferentes especifica o maior e o menor para o de maior tamanho continuar executando enquanto o menor para 
        maior = y1;
        menor = y3;
    else
        maior = y3;
        menor = y1;
    end     
    
        for i = 1:length(maior)                       % Realiza o loop até o tamanho de maior (maior continua e menor para caso tenham tamanhos diferentes)
            if isgraphics(ax1)                        % Verifica se o handle ax1 é valido
                cla(ax1);                             % Limpa a figura anterior
                hold(ax1,'on');                       % Retém o gráfico corrente
                axis(ax1,limites_grafico);            % Define os eixos x e y
                grid(ax1,'on');                       % Habilita a grade
                p1 = plot(ax1,[maior(i,1), maior(i,1)+lh*sin(maior(i,2))], [0, lh*cos(maior(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax1,'Position',[limites_grafico(1,1),-0.2,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.1],'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax1,'Position',[maior(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                
hold(ax1,'off');                      % Retém o gráfico corrente
            end
            if i <= length(menor) & isgraphics(ax2)   % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax2);                             % Limpa a figura anterior
                hold(ax2,'on');                       % Retém o gráfico corrente
                axis(ax2,limites_grafico);            % Define os eixos x e y             
                grid(ax2,'on');                       % Habilita a grade
                p3 = plot(ax2,[menor(i,1), menor(i,1)+lh*sin(menor(i,2))], [0, lh*cos(menor(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax2,'Position',[limites_grafico(1,1),-0.2,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.1],'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax2,'Position',[menor(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax2,'off');                      % Libera o gráfico corrente
            end
            drawnow();                                % Atualiza a figura com os dados anteriores                               
            if ~isempty(gcf) & gravar                 % Verifica se a figura é válida e se gravar = 1 (true)
                F(i) = getframe(gcf);                 % Obtém a imagem da figura corrente
                writeVideo(video,F(i));               % Armazena a imagem
            end
        end  
elseif ~isempty(y2) & ~isempty(y3) & ~isempty(y2) & ~isempty(y4) % Se a animação for para três sistemas
    f = figure;                                      % Cria um figura
    f.Position = [250 50 850 620]; % [x1,y1,x2,y2]   % Tamanho e posição configurados para 3 gráficos [x1,y1,x2,y2]
    tl = tiledlayout(2,2)                            % Cria um layout com 1 linha e três colunas
    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';
    title(tl,{titulo;subtitulo});        % Título mestre da figura
    xlabel(tl,'Posição Horizontal - m'); % Texto do eixo x
    ylabel(tl,'Posição Vertical - m');   % Texto do eixo y

    ax1 = nexttile;                                  % Grafico da medição do sistema não-linear
    title(strcat(titulo1,char(160),char(8658),char(160),'ci =',char(160),ci1));   % Título do grafico 1
    ax2 = nexttile;                                  % Grafico do EKF
    title(strcat(titulo2,char(160),char(8658),char(160),'ci =',char(160),ci2));   % Título do gráfico 2
    ax3 = nexttile;                                  % Grafico da medição do sistema linear
    title(strcat(titulo3,char(160),char(8658),char(160),'ci =',char(160),ci3));   % Título do gráfico 3
    ax4 = nexttile;                                  % Grafico do KF
    title(strcat(titulo4,char(160),char(8658),char(160),'ci =',char(160),ci3));   % Título do gráfico 3
    
    for i = 1:length(y1)                              % Loop até o tamanho de maior
            if isgraphics(ax1)                        % Verifica de ax1 é válido
                cla(ax1);                             % Limpa a figura anterior
                hold(ax1,'on');                       % Retém o gráfico corrente
                axis(ax1,limites_grafico);            % Define os eixos x e y
                grid(ax1,'on');                       % Habilita a grade
                p1 = plot(ax1,[y1(i,1), y1(i,1)+lh*sin(y1(i,2))], [0, lh*cos(y1(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax1,'Position',[limites_grafico(1,1),-0.2,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.1],'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax1,'Position',[y1(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax1,'off');                      % Retém o gráfico corrente
            end
            if  isgraphics(ax2)                       % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax2);                             % Limpa a figura anterior
                hold(ax2,'on');                       % Retém o gráfico corrente
                axis(ax2,limites_grafico);            % Define os eixos x e y             
                grid(ax2,'on');                       % Habilita a grade
                p3 = plot(ax2,[y2(i,1), y2(i,1)+lh*sin(y2(i,2))], [0, lh*cos(y2(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax2,'Position',[limites_grafico(1,1),-0.2,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.1],'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax2,'Position',[y2(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax2,'off');                      % Libera o gráfico corrente
            end
            if isgraphics(ax3)                        % Continua animando o maior e para o menor e verifica se ax3 ainda é válido
                cla(ax3);                             % Limpa a figura anterior
                hold(ax3,'on');                       % Retém o gráfico corrente
                axis(ax3,limites_grafico);            % Define os eixos x e y             
                grid(ax3,'on');                       % Habilita a grade
                p3 = plot(ax3,[y3(i,1), y3(i,1)+lh*sin(y3(i,2))], [0, lh*cos(y3(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax3,'Position',[limites_grafico(1,1),-0.2,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.1],'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax3,'Position',[y3(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax3,'off');                      % Libera o gráfico corrente
            end
            if isgraphics(ax4)                        % Continua animando o maior e para o menor e verifica se ax3 ainda é válido
                cla(ax4);                             % Limpa a figura anterior
                hold(ax4,'on');                       % Retém o gráfico corrente
                axis(ax4,limites_grafico);            % Define os eixos x e y             
                grid(ax4,'on');                       % Habilita a grade
                p3 = plot(ax4,[y4(i,1), y4(i,1)+lh*sin(y4(i,2))], [0, lh*cos(y4(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                rectangle(ax4,'Position',[limites_grafico(1,1),-0.2,abs(limites_grafico(1,1))+abs(limites_grafico(1,2)),0.1],'FaceColor',[0.8500 0.3250 0.0980]) % Carrinho
                rectangle(ax4,'Position',[y4(i,1)-lc/2,-hc,lc,hc],'FaceColor','magenta') % Trilhos
                hold(ax4,'off');                      % Libera o gráfico corrente
            end
            drawnow();                                % Atualiza a figura com os dados anteriores                               
            if ~isempty(gcf) & gravar                 % Verifica se a figura é válida e se gravar = 1 (true)
                F(i) = getframe(gcf);                 % Obtém a imagem da figura corrente
                writeVideo(video,F(i));               % Armazena a imagem
            end
        end                                         
end
%% Fecha o arquivo de vídeo
if gravar % Se gravar = 1 (true)
    close(video);
end
end
