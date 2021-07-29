

  
% deleta todas as variáveis na memória
clear all

% fecha todas as janelas abertas do sistema
close all

% limpa a tela
clc

% dados da geometria
A=1;
B=1;
E_1=69e9;
E_2=69e9;
A_1=pi/4*(14.3e-3)^2;
A_2=pi/4*(14.3e-3)^2;
% massa da partícula
m=0.5;

%  espaço de trabalho
x=0.65:10e-3:1.5;
y=0.65:10e-3:1.5;
[X,Y]=meshgrid(x,y);

% criação das matrizes dos deslocamentos e frequências
dx=size(X);
dy=size(X);
omega=size(X);

for k=1:size(X,1)
for p=1:size(Y,2)

% Matriz jacobiana
% Js*qs=Je*qe
% Matriz Js

Jx=[2*(X(k,p)-A) 2*Y(k,p);
    2*X(k,p)     2*(Y(k,p)-B)];

% Matriz Je
L_1=sqrt((X(k,p)-A)^2+Y(k,p)^2);
L_2=sqrt((X(k,p))^2+(Y(k,p)-B)^2);
Jq=[2*L_1 0;0 2*L_2];

% Matriz "Qui" rigidez local do atuadores
k_1=E_1*A_1/L_1;
k_2=E_2*A_2/L_2;

Ks=[k_1 0;
    0   k_2];

% J=Jq\Jx;
% K=J'*Ks*J;
% delta=K\[0;1];

% Matriz de compliância para evitar a inversão de uma matriz retangular
Comp=inv(Jx)*Jq*inv(Ks)*transpose(Jq)*inv(transpose(Jx));

% Deslocamento
delta=Comp*[0;1];
% Matriz de rigidez a partir da matriz compliância
K=inv(Comp);

% dados para plotagem da superfície
dx(k,p)=delta(1,1);
dy(k,p)=delta(2,1);

% Cálculo da matriz de rigidez do sistema via conceito de massa concentrada
% do sistema

% Matriz de massa dos elos ativos
M=[m 0;
   0 m];
      
% Cálculo da primeira frequência natural a partir das matrizes de massa e 
% de rigidez
[W] = eig(K,M);
omega(k,p)=(sqrt(min(W)))/(2*pi);
end
end

% subplot(1,2,1)
% surf(Q3,-Q4,dx)
% title('Displacement in x (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dx(mm)')
% subplot(1,2,2)
[~,h] = contour(X*1000,Y*1000,dx*1000,15);
title('Displacement in x Lumped Parameters (mm)')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on

figure
% subplot(1,2,1)
% surf(Q3,-Q4,dy)
% title('Displacement in y (Concentrated Parameters)')
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('dy(mm)')
% subplot(1,2,2)
[~,h] = contour(X*1000,Y*1000,dy*1000,15);
title('Displacement in y Lumped Parameters (mm)')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on

figure
% subplot(1,2,1)
% surf(Q3,-Q4,omega)
% title('First natural frequency (Concentrated Parameters)')
% xlabel('x(m)')
% ylabel('y(m)')
% zlabel('\omega (Hz)')
% subplot(1,2,2)
[C,h] = contour(X*1000,Y*1000,omega,15);
title('First natural frequency Lumped Parameters (all terms considered) (Hz)')
xlabel('x(mm)')
ylabel('y(mm)')
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
grid on
