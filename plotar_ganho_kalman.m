function plotar_ganho_kalman(t,L_k,w)
%%     Função que calcula a norma 2 das duas colunas do Ganho de Kalman e mostra no gráfico.
L_1 = norm(L_k(:,1),2);
L_2 = norm(L_k(:,2),2);
for i = 3:2:length(L_k)
 L_1 = [L_1 norm(L_k(:,i),2)];
 L_2 = [L_2 norm(L_k(:,i+1),2)];
end
plot(t,L_1','LineWidth',w)
hold on
plot(t,L_2','LineWidth',w)
hold off
title('Norma 2 das colunas do Ganho de Kalman vs.Tempo');
ylabel('Ganho de Kalman');
xlabel('Tempo - s');
grid on;
legend('Referente à Saída 1 - Posição do Carrinho','Referente à Saída 2 - Ângular da haste')
end
