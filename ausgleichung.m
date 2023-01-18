function [ xdach, v, sx ] = ausgleichung( A, l )

N = A'*A; %Normalgleichungsmatrix
Q = inv(N); %Inversion Normalgleichungsmatrix
n = A'*l; %Absolutglied
xdach = Q*n; %Gekürzter Parametervektor
v = A*xdach-l; %Verbesserungsvektor
W = v'*v; %Verbesserungsquadratsumme
zw = size(A); n = zw(1); u = zw(2); %Freiheitsgrade
s0 = sqrt(W/(n-u)); %empir. Stdabw. der Gewichtseinheit
sx = zeros([u,1]);
for i = 1:u
    sx(i,1) = s0 * sqrt(Q(i,i));
end


end

