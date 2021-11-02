function [driftz,drifty] = Dipolhertzquadrandsize(z0,y0,betas,alphas,kas,R0,j)
%interparticle force function for variable particle sizes
%version5 is performance optimized
%the performance limiting steps is now the computation of  d=sqrt(d2) and the calling of the function 

% extract the particle for which we will compute the drift
%z and y are the position of particle zy on which the other particels exert a force
z=z0(j);
y=y0(j);
R=R0(j);
d2 = (z-z0).^2 + (y-y0).^2; % define an array of  the square of the particle distance 
d=sqrt(d2); % define an array of the particle distances

%creat an logic array to later find all intersecting particles
logic= ((d-R0-R)<0);

% defin cos(theta)^2
therm= (z-z0).^2./d2;

 
%compute the drift vector from all particles 

%some performance definitions (reapearing variables are defined here)
  d=d2.*d;
 a=alphas.*(    1*(1-3*therm)./(d));
 d=d.*d2;
%final equation
Fz =(z-z0).*(    (betas)  .*(3-5*therm)./d+          a)                             ;
Fy =(y-y0).*(    (betas)  .*(1-5*therm)./d+          a)                             ; 


% recompute the drift for the particles that intersect
z0=z0(logic);
y0=y0(logic); 
kas=kas(logic);
betas=betas(logic);
alphas=alphas(logic);
R0=R0(logic);
%recomputation is better than finding it 
dd2 = (z-z0).^2 + (y-y0).^2; 
dd=sqrt(dd2);

thermd= (z-z0).^2./dd2;

% compute the repulsion force term
F=abs(dd-R0-R).*sqrt(abs(dd-R0-R));

%replace the forces for the particles that intersected the zy particle
Fz(logic) =(z-z0).*(    (betas)   .*(3-5*thermd)./((R0+R).^4.*dd)+    alphas.*(    1*(1-3*thermd)./((R+R0).^2.*dd) )       +kas.*F./dd)  ;
Fy(logic) =(y-y0).*(    (betas)   .*(1-5*thermd)./((R0+R).^4.*dd)+    alphas.*(    1*(1-3*thermd)./((R+R0).^2.*dd) )      +kas.*F./dd) ;

% eliminate the drift from the particle zy intselfe onto intselfe
Fz(j)=0;
Fy(j)=0;

% add all the drifts together and  and return them
drifty=sum(Fy);


driftz=sum(Fz);


