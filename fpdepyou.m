% where u is the noisy input image.
function [frth]=fpdepyou(u,iteration,dt)
% u = Noisy Image
% T - Threshold , Based on the threshold you will get scale sapce images. At a particular 
% value of T you will gwt the converged image. 
% Time step   -  0.25 a good experimental step time

%The laplacian operator is calculated on seperate steps as the MatLab
%cannot calculate it directly.
u=double(u);

for  t=1:iteration
    
    [ux,uy]=gradient(u);                         %get the gradiant of u
    [uxx,~]=gradient(ux);                        %get the gradiant of ux to be uxx
    [~,uyy]=gradient(uy);                        %get the gradiant of uy to be uyy
                                                 %Now we have uxx and uyy
    c=1./(1.+sqrt(uxx.^2+uyy.^2)+0.0000001);     %get g(n) = c(|nabla^2(u)|)nabla^2(u)
    [div1,~]=gradient(c.*uxx);                   %get get the gradiant of g(n) to get gx
    [~,div2]=gradient(c.*uyy);                   %get get the gradiant of g(n) to get gy
    [div11,~]=gradient(div1);                    %get get the gradiant of gx to get gxx
    [~,div22]=gradient(div2);                    %get get the gradiant of gy to get gyy
                                                 %Now we have nabla^2(g(n))
    div=div11+div22;                             %nabla^2[g(n)] 
    
    u=u-(dt.*div);                               %u(n+1) = u(n) - dt*nabla^2[g(n)]
                
end
frth=uint8(u);                                   % Converting to 8 bit image