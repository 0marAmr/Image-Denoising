function [u]=speckle_sol(A)
A = im2double(A);
k = 0.5;
[m , n] = size(A);
Med = [];
Med2 = [];
%Modified filter 
for i=2:m-1
    for j=2:n-1
            Med(1) =A(i-1,j) ;
            Med(2) = A(i,j-1);
            Med(3) = A(i,j+1);
            Med(4) = A(i+1,j);        
            Med2 = Med.^2;
            
            M = median(Med);
            h = median(Med2);
            q = h - M^2;
            
             if  abs(A(i,j) - M)^2 > k*q
                A(i,j) = M ;
             end
                
    end    
end
rescaledMatrix = rescale(A, 0, 255);
% Typecasted to uint8 
u = uint8(rescaledMatrix); 
% u = A;