function animar_pendulo(y1,y2,s,l_haste,l_carrinho,h_carrinho,titulo1,titulo2)
    % x -> Estado atual do sistema
    % s -> Escala da animação para aumento ou diminuição das dimensão do desenho
    % l_haste -> Comprimento da haste
    % l_carrinho -> Comprimento do carrinho
    % h_carrinho -> Altura do carrinho

    lh = l_haste*s;            % Comprimento da haste em escala
    lc = l_carrinho*s;         % Comprimento do carrinho em escala
    hc = h_carrinho*s;         % Altura do carrinho em escala
if isempty(y2)
    for i = 1:length(y1)
        clf();                 % Limpa a figura anterior
        hold on;               % Retém o gráfico corrente
        axis([-4,4,-3,3]);     % Define os eixos x e y  
        grid on;               % Habilita a grade
        plot([ y1(i,1), y1(i,1)+lh*sin(y1(i,2))], [0, lh*cos(y1(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
        plot(y1(i,1)+[-lc/2,lc/2,lc/2,-lc/2,-lc/2], [0,0,-hc,-hc,0],'magenta','LineWidth',2); % Desesnha o carrinho a partir da posição corrente 
        title(strcat('Animação -',char(160),titulo1));
        xlabel('Posição do carrinho - m');
        ylabel(strcat('l*Cos ',char(160),char(952)));
        hold off;              % Libera o gráfico corrente
        drawnow();             % Atualiza a figura com os dados anteriores
    end
else
    tiledlayout(1,2)
    ax1 = nexttile;
    title(strcat('Animação -',char(160),titulo1));
    xlabel('Posição do carrinho - m');
    ylabel(strcat('l*Cos',char(160),char(952)));
    ax2 = nexttile;
    title(strcat('Animação -',char(160),titulo2));
    xlabel('Posição do carrinho - m');
    ylabel(strcat('l*Cos ',char(160),char(952)));
    
    if length(y1) >= length(y2) % No caso de 2 vetores de tamanhos diferentes especifica o maior e o menor
        maior = y1
        menor = y2
    else
        maior = y2
        menor = y1
    end     
    
        for i = 1:length(maior)
            if isgraphics(ax1)
                cla(ax1);                             % Limpa a figura anterior
                hold(ax1,'on');                       % Retém o gráfico corrente
                axis(ax1,[-4,4,-3,3]);                % Define os eixos x e y
                grid(ax1,'on');                       % Habilita a grade
                p1 = plot(ax1,[maior(i,1), maior(i,1)+lh*sin(maior(i,2))], [0, lh*cos(maior(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                p2 = plot(ax1,maior(i,1)+[-lc/2,lc/2,lc/2,-lc/2,-lc/2], [0,0,-hc,-hc,0],'magenta','LineWidth',2); % Desesnha o carrinho a partir da posição corrente 
                hold(ax1,'off');                      % Retém o gráfico corrente
            end
            if i <= length(menor) & isgraphics(ax2)   % Continua animando o maior e para o menor e verifica se ax2 ainda é válido
                cla(ax2);                             % Limpa a figura anterior
                hold(ax2,'on');                       % Retém o gráfico corrente
                axis(ax2,[-4,4,-3,3]);                % Define os eixos x e y             
                grid(ax2,'on');                       % Habilita a grade
                p3 = plot(ax2,[menor(i,1), menor(i,1)+lh*sin(menor(i,2))], [0, lh*cos(menor(i,2))],'blue','LineWidth',2);  % Desenha a haste apartir da posição do carrinho e ângulo da haste corrente  
                p4 = plot(ax2,menor(i,1)+[-lc/2,lc/2,lc/2,-lc/2,-lc/2], [0,0,-hc,-hc,0],'magenta','LineWidth',2); % Desesnha o carrinho a partir da posição corrente        
                hold(ax2,'off');                      % Libera o gráfico corrente
            end
            drawnow();                                % Atualiza a figura com os dados anteriores
        end     
end
end
