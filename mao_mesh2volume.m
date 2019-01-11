%given mesh for closed surface, compute volume based on divergence theorem
function [volume area]=mao_mesh2volume(V,F);
area=0;
volume=0;
for j=1:size(F,1)
    V1=V(F(j,1),:);
    V2=V(F(j,2),:);
    V3=V(F(j,3),:);
    a=V2-V1;
    b=V3-V1;
    c=cross(a,b);
    M=(V1+V2+V3)/3;
    n=c/norm(c);
    T=M(1)*n(1)*0.5*norm(c);
    if ~isnan(T)
    volume=volume+T;
    end
    area=area+0.5*norm(c);
end
end