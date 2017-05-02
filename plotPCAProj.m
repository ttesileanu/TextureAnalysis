function [N,vecP] = plotPCAProj(evA,evB,v1,v2)
 
%determine the discriminant that separates ev1 and ev2
%project this discriminant onto the plane spanned by v1 and v2
 
evA = evA(:,2:10);
evB = evB(:,2:10);
 
data   = [evA;evB];
labels = [zeros(size(evA,1),1);ones(size(evB,1),1)];
obj    = fitcdiscr(data,labels);  
 
 
%get discriminant (normal to best separating hyperplane)
N = obj.Coeffs(1,2).Linear'; %normal vector, unit length
N = N/norm(N);
 
 
%get separating point
xL   = -obj.Coeffs(1,2).Const./norm(obj.Coeffs(1,2).Linear);
vec  = [(xL-.5)*N;xL*N;(xL+.5)*N];
proj = [v1,v2];
 
vecP  = vec*proj;
projA = evA*N';
projB = evB*N';
 
%solid line: component A
if mean(projA<xL)
    plot(vecP(1:2,1),vecP(1:2,2),'-o')
    plot(vecP(2:3,1),vecP(2:3,2),'--o')
else
    plot(vecP(1:2,1),vecP(1:2,2),'--o')
    plot(vecP(2:3,1),vecP(2:3,2),'-o')
end

end