function a = simplex_project_dist(x, E)
% SIMPLEX_PROJECT_DIST Project point onto simplex using distance geometry
%
% This function projects the points in x onto the simplex defined by the
% distance matrix E using distance geometry. All points are expressed in
% terms of distances to the vertices. The barycentric coordinates of the
% projected points are returned. The projection is defined as the point 
% within the simplex with minimal L2-distance.
%
% Input:
% x : pxM vector, the squared distances of the M point to project to the vertices
% E : pxp matrix, the squared distance matrix of the p vertices
%
% Output:
% a : pxM vector, barycentric coordinates of the projected points
%
% Author: Rob Heylen
% Visionlab, University of Antwerp, september 2010

delta1=1e-4;
delta2=1e-8;

[p,M]=size(x);
a=zeros(p,M);

vol=factorial(p-1)*calc_simplex_volume(E);
if vol<delta2
    a=0;
    error('Vertices are linearly dependent');
end
if vol<delta1
    warning('Vertices are almost linearly dependent, results may vary');
end

if p==1
    % Special case, we know the answer immediately
    a=ones(1,M);
    return;
end

% Project x onto the simplex plane
x=plane_project(x,E);

if p==2
    d=sqrt(E(1,2));
    x=sqrt(x);
    a([2 1],:)=x/d;
    mask=logical(sum(real(a>1)));
    a(1,mask)=real(a(1,mask)>a(2,mask));
    a(2,mask)=1-a(1,mask);
    return;
end
    
% Find the set of points inside the simplex
s=zeros(p,M);
for i=1:p
    Er=E([1:i-1 i+1:p],[1:i-1 i+1:p]);
    Ei=E([1:i-1 i+1:p],i);
    xr=x([1:i-1 i+1:p],:);
    d=calc_dist2(xr,Ei,Er);
    s(i,:)=(abs(x(i,:)-d(1,:)))<((d(2,:)-d(1,:))/2);
end
s=prod(real(s),1);

% For the points inside the simplex, determine the barycentric coordinates
I=find(s);
xr=x(:,I);
a(:,I)=calc_bary_coor(xr,E);

if length(I)==M
    return;
end

% Calculate the incenter
ac=zeros(p,1);
for i=1:p
    Er=E([1:i-1 i+1:p],[1:i-1 i+1:p]);
    ac(i)=calc_simplex_volume(Er);
end
ac=ac/sum(ac);
c=calc_dist_coor(ac,E);
%xc=calc_dist(x,c*ones(1,M),E);
xc=calc_dist2(x,c,E);

I=find(~s);
Mi=length(I);

for i=1:p
    s=zeros(p-1,Mi);
    si=[1:i-1 i+1:p];
    for j=1:p-1
        sj=si([1:j-1 j+1:p-1]);
        Er=[E(sj,sj) c(sj); c(sj)' 0];
        %Ej=[E(sj,si(j)) ; c(si(j))]*ones(1,Mi);
        Ej=[E(sj,si(j)) ; c(si(j))];
        xr=[x(sj,I); xc(1,I)];
        %d=calc_dist(xr,Ej,Er);
        d=calc_dist2(xr,Ej,Er);
        s(j,:)=(abs(x(si(j),I)-d(1,:)))<((d(2,:)-d(1,:))/2);
    end
    s=prod(real(s),1);
    Ic=I(logical(s));
    Mc=length(Ic);
    if Mc>0
        %Project point onto simplex plane
        x(si,Ic)=plane_project(x(si,Ic),E(si,si));
        a(si,Ic)=simplex_project_dist(x(si,Ic),E(si,si));
    end
    I=setdiff(I,Ic);
    Mi=length(I);
    if Mi==0
        break;
    end
end

% This can happen when there vertices are almost linearly dependent
if Mi~=0
    for i=1:Mi
        a(x(:,I(i))==min(x(:,I(i))),I(i))=1;
    end
end
end





function Dnew=plane_project(X,E)
%PLANE_PROJECT Project points onto hyperplane
%
% Input:
% x: pxM vector of M sets of squared distances to the vertices
% E: pxp squared distance matrix of the p vertices of the simplex.
%
% Output:
% Xnew: vector containing the squared distances of the projected points

[p,M]=size(X);
C=[E ones(p,1); ones(1,p) 0];
Ci=C^(-1);
Xi=[X ; ones(1,M)];
do=Ci*Xi;
do=sum(Xi.*do)/2;
Dnew=X-ones(p,1)*do;
end



function d=calc_dist2(x,y,E)
% CALC_DIST2 Calculate distance between points using distance geometry
%
% This function calculates the unknown distance between all the points in x
% and the single point in y, given the distances to a set of vertices 
% forming a simplex, and the mutual distances between these vertices
%
% Input:
% x: pxM matrix of M sets of squared distances of x to the vertices
% y: px1 vector of the second point y
% E: pxp squared distance matrix of the p vertices of the simplex.
%
% Output:
% d:      vector containing the two different possibilities for the squared
%         distance, sorted


[p,M]=size(x);
D=E^(-1);
q11=(ones(1,p)*(D*ones(p,1)))*ones(1,M);
q21=(D*x);
q31=(ones(1,p)*(D*y))*ones(1,M);
q22=sum(x.*q21);
q32=y'*q21;
q33=(y'*(D*y))*ones(1,M);
q21=ones(1,p)*q21;

A = -q11;
B = 2*q11.*q32 - 2*(q21-1).*(q31-1);
C = q11.*q22.*q33 + 2*q32.*(q21-1).*(q31-1) - q22.*(q31-1).^2 - q33.*(q21-1).^2 - q11.*q32.^2;

T=B.^2-4*A.*C;
d(1,:)=real((-B-sqrt(T))./(2*A));
d(2,:)=real((-B+sqrt(T))./(2*A));
d=sort(d);
end



function a = calc_bary_coor(x, E)
%CALC_BARY_COOR Calculate barycentric coordinates in distance geometry
% This function calculates the barycentric coordinates of the points given
% with distance coordinates in x with respect to the distance matrix E
%
% Input:
% x:  pxM matrix of M sets of distance coordinates of x
% E:  pxp squared distance matrix of the p vertices of the simplex.
%
% Output:
% a:  pxM matrix of the barycentric coordinates of x

[p,M]=size(x);
if p==2
    d=sqrt(E(1,2));
    x=sqrt(x);
    a([2 1],:)=x/d;
    mask=logical(sum(real(a>1)));
    a(1,mask)=a(1,mask).*(a(1,mask)>a(2,mask)) - a(1,mask).*(a(1,mask)<=a(2,mask));
    a(2,mask)=1-a(1,mask);
    return;
end

for i=1:p
    V=calc_simplex_volume(E([1:i-1 i+1:p],[1:i-1 i+1:p]));
    xi=[x([1:i-1 i+1:p],:) ; ones(1,M)];
    C=[E([1:i-1 i+1:p],[1:i-1 i+1:p]) ones(p-1,1); ones(1,p-1) 0];
    Ci=C^(-1);
    do=Ci*xi;
    a(i,:)=V*sqrt(sum(xi.*do));
    
    d=calc_dist2(x([1:i-1 i+1:p],:),E([1:i-1 i+1:p],i),E([1:i-1 i+1:p],[1:i-1 i+1:p]));
    mask=2*((abs(x(i,:)-d(1,:)))<((d(2,:)-d(1,:))/2))-1;
    a(i,:)= a(i,:).*mask;

    
end
a=a./(ones(p,1)*sum(a));
end


function d=calc_bary_dist(ax,ay,E)
% CALC_BARY_DIST Calculate distance between two points given baryc. coord.
%
% This function calculates the unknown distance between x and y, given the
% barycentric coordinates of both with respect to a set of vertices of
% which one only knows the distance matrix.
%
% Input:
% ax, ay: pxM matrices of M sets of barycentric coordinates of x and y
% E:      pxp squared distance matrix of the p vertices of the simplex.
%
% Output:
% d:      vector containing the squared distances


[p,M]=size(ax);
J=eye(p,p)-(1/p)*ones(p,1)*ones(1,p);
A=-0.5*J*E*J;
d=A*(ax-ay);
d=sum((ax-ay).*d);
end


function out=calc_simplex_volume(D)
%SIMPLEXVOLUME Returns simplex volume of distance matrix D

K=size(D,1);
out=abs(det([D ones(K,1); ones(1,K) 0]));
out=sqrt(out/(2^(K-1)*(factorial(K-1)^2)));

end



function d=calc_dist_coor(ax,E)
% CALC_DIST_COOR Calculate distance coordinates given barycentric
% coordinates
%
% This function calculates the distance coordinates of a point with respect
% to the vertices defined by distance matrix E given a set of barycentric
% coordinates.
%
% Input:
% ax: pxM matrix of M sets of barycentric coordinates of x
% E:  pxp squared distance matrix of the p vertices of the simplex.
%
% Output:
% d:  pxM matrix of the distance coordinates of x


[p,M]=size(ax);
J=eye(p,p)-(1/p)*ones(p,1)*ones(1,p);
A=-0.5*J*E*J;
d=zeros(p,M);
for i=1:p
    at=zeros(p,M);
    at(i,:)=ones(1,M);
    ad=(ax-at);
    t=A*ad;
    d(i,:)=sum(ad.*t)';
    
end
end