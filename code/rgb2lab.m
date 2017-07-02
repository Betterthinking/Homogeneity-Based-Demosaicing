function [ output ] = rgb2lab(fig)
%rgb2lab : rgb2lab(image)
%                @brief: convert a RGB figure to LAB color space   

try
    R = fig(:,:,1);
    G = fig(:,:,2);
    B = fig(:,:,3);
catch
    disp('Error! The input figure must be 3 channle input with RGB');    
end

shape2D = size(R);

r = gamma_lab(R);
g = gamma_lab(G);
b = gamma_lab(B);

M = [0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505];

X = zeros(shape2D);
Y = zeros(shape2D);
Z = zeros(shape2D);
for i = 1:shape2D(1)
    for j = 1:shape2D(2)
        temp = M*[r(i,j);g(i,j);b(i,j)];
        X(i,j) = temp(1);
        Y(i,j) = temp(2);
        Z(i,j) = temp(3);
    end
end

Xn = 95.074;
Yn = 100;
Zn = 108.883;

L = 116.*transfunc(Y./Yn)-16;
A = 500.*(transfunc(X./Xn)-transfunc(Y./Yn));
B = 200.*(transfunc(Y./Yn)-transfunc(Z./Zn));

final(:,:,1) = L;
final(:,:,2) = A;
final(:,:,3) = B;
output = final;
end

function [output] = gamma_lab(X)
    T = double(X)./255;
    T(T>0.04045) = ((T(T>0.04045)+0.055)./1.055).^2.4;
    T(T<=0.04045) = (T(T<=0.04045)/12.92);
    output = T;
end

function [output] = transfunc(X)
    T = X;
    threshold = 6^3/29^3;
     T(T>threshold) = ((T(T>threshold)).^(1/3));
    T(T<=threshold) = (T(T<=threshold).*(29/6)^2./3+4/29);
    output = T;
end