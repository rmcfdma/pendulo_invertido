function F = equacoes_nao_lineares(t,x,K)
l = 0.641; % Comprimento da Haste
m = 0.230; % Massa do Pêndulo
M = 0.57; % Massa do Carrinho
g = 9.81; % Aceleração da Gravidade
u = -K*x;
F = zeros(4,1);
F(1) = x(2);
F(2) = (u+m*l*x(4)^2*sin(x(3))-m*g*sin(x(3))*cos(x(3)))/(M+m-m*cos(x(3))^2);
F(3) = x(4);
F(4) = (u*cos(x(3))-(M+m)*g*sin(x(3))+m*l*x(4)^2*sin(x(3))*cos(x(3)))/(m*l*cos(x(3))^2-(M+m)*l);
end

