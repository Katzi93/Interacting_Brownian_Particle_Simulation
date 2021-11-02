% This is a version for variable sized particles 
% Simulation of many interacting brownian particles
% This code was written by Katzi ( Forian Katzmeier)
% Take care this code can produce 8 GB output data
% A quick run can be made be reducing the squrt of particle to 10
% The output data are 1 video file and two txt files containing coordinates
% the consol will print tau which is the scaling of the time dimension
clear all;
close all;
clc;





%set time
time =	40;              %set total time of the simulation in dimensionless time
steps = time*1000;      %define the total number of steps of the  simulation
dt = time/steps;        %define time interval for integration
hop = 25 ;              %define number of integration steps between frames of the video


%set space and particles
sqrtparticles=33;                %initialize the root of particles. The actual particle count will be sqrtparticles^2
particles = sqrtparticles^2;    %compute the actual number of particles
Rscale=0.65*10^-6;                    % scale for the radii can be the average of R0 or anyting you want 
dist= 0.5;                      %distance in muliple particle radii (the center center distance is dist+2*Rmax and Rmax ist the maximum particle radius)

   








%Natural constants/Physical parameter definitions in Si units to compute
%the real time intervall between frames and 


kbT = 1.381*10^-23*298   ;            %Thermal energy as Boltzman konstant times temperature in N*m at 25C
nu = 1.0002*10^-(3) ;               % Dynamic viscosity of Water in (N/(m)^2) *s
epsilon_vac = 8.854*10^-12  ;         % Vacuum permitivity in coulomb^2/((m)^2*N)
epsilon_water = 78.4*epsilon_vac ;
elementarycharge = 1.602*10^-19;    % elementary charge in coulomb
mol = 6.002*10^23;             % Advogradroconstant

tau= 6 *pi*nu*Rscale^3/kbT          %scaling of the time also print it

frame= dt*tau*hop % time between frames in seconds
ttwo= tau*time %real duration of the simulation in seconds
%computation of simulation parameters in µm and seconds


E=(288*1.0/0.017); % electric Field in Volts/m



             % sigma for random number computation


         
                                 %  constant that includes particle repulsion




R0 = wblrnd(0.66,2.1,1,particles)*10^-6; % create an array of random paricle sizes units of wblrnd are in µm
%R0= ones(1,particles)*0.65*10^-6;
% scale for the radii can be the average of R0 or anyting you want                                     
R0=R0/Rscale; % make R0 dimensionless


Rmax=max(R0);                       % find the maximum of the radii
range = (sqrtparticles-1)*(2*Rmax+dist);              %define a range for plotting
actualsize=range*Rscale
dimensionfactor = epsilon_water *6*pi *Rscale^2*E^2/kbT*Rscale; %the  dimensionless scaling factor
sigmas = sqrt(2*(1./R0)*dt);  % compute the standart deviations for stochasitc differential
gamma=0.088; %set gamma
alphas= gamma/2*dimensionfactor*R0.^3; %fill a list with all values for gamma
alphas= [alphas alphas alphas alphas alphas alphas alphas alphas alphas]; % make it periodic for the later defined periodic boundarie condistionss

betas= zeros(particles); % initialize an array for alll the betas
Kdhoch2=0.23; % define Kdhoch 2


tau= 6 *pi*nu*Rscale^3/kbT % set the time scale
frame= dt*tau*hop % time between frames in seconds
ttwo= tau*time
for j=1:particles
betas(j,:)= dimensionfactor*R0.^3.*(R0(j)^2*Kdhoch2- gamma/2*R0.^2-gamma/2*R0(j)^2); % fill the matrix for beta 

end 
betas = [betas betas betas betas betas betas betas betas betas]; % make beta periodic for the boundaries

kold= 75 % define the repulsion constant in the old variables (in the text k s nu)
k=kold*kbT/(Rscale^(3)); % compute the new prefactor which can also be set here
kdimensionfactor= k*Rscale^(3)/kbT; % define the dimension factor of the repulsive force
kas= zeros(particles); % initialize the array for the k's
for j=1:particles
kas(j,:)= kdimensionfactor*(1/R0(j))*(1./R0+1/R0(j)).^(-1./2.); % fill the array for k

end 

kas=[kas kas kas kas kas kas kas kas kas]; % make k periodic 

R0 = [R0 R0 R0 R0 R0 R0 R0 R0 R0]; % make r periodic


    
    
    
    
    
    




% put the particles on a square grid (initial conditions)
 z = -range/2:(2*Rmax+dist):range/2;
 y = -range/2:(2*Rmax+dist):range/2;
 [z,y] = meshgrid(z,y);
 z1 = reshape(z,1,[]).';
 y1 = reshape(y,1,[]).';
 
% initialize the array where the solution of the simulation is saved
 z= zeros(particles, steps);
 z(:,1)=z1;
 y= zeros(particles, steps);
 y(:,1)=y1;
 
% some array which is needed to store particle positions
 zold=zeros(particles);
 yold=zeros(particles);
 zold2=zeros(particles*9);
 yold2=zeros(particles*9);
 
 %define periodic boundary connditions
 % confine partices
 ymin=min(y1)-(2*Rmax+dist)/2; %  the (2*R+dist)/2 is needed that the particles on the edge have still a distance
 ymax= max(y1)+(2*Rmax+dist)/2;
 
 zmin=min(z1)-(2*Rmax+dist)/2;
 zmax= max(z1)+(2*Rmax+dist)/2;


%solution sceme
t=0; % set starting time to 0



%randomize the particle positions 
%run the simulation for 1000 steps with only the repulsion force acting
% i and j are not the same in the skript and here.
% We compute the drift for particle j
for i=1:2000
    %define some convenient variables
    zold=1*z(:,1)+0; % I use a 1 here because I dont want to save the randomization process
    yold=1*y(:,1)+0 ;                    % extend matrix
    
    
    % set the drift to zero for this time step
    dy=zeros(1,particles);    
    dz=dy;
    % make periodic boundaries i.e. copy the array around itselfe
     zold2=[zold  zold          zold                (zold+2*zmax)       (zold+2*zmin)   (zold+2*zmax) (zold+2*zmin) (zold+2*zmax) (zold+2*zmin)];             
     yold2=[yold (yold+2*ymax)  (yold+2*ymin)        yold                yold           (yold+2*ymax) (yold+2*ymin) (yold+2*ymin) (yold+2*ymax)];
   
   %iterate over all particles
   for j=1:particles 
       
          nullen= zeros(1,9*particles);
    % compute the drift for one partikel without alpha and beta 
   
             [dz(j),dy(j)] =Dipolhertzquadrandsize(zold2(:),yold2(:),nullen(:),nullen(:),transpose(kas(j,:)),R0(:),j); % We use R(:) instead of R(:) to get the dimensions right
            
            
        
   end
   % add the random number to the drift (for randomization increase sigma 10 times
   dz=dz*dt+10*sigmas.*randn(1,particles) ;
   dy=dy*dt+10*sigmas.*randn(1,particles) ;
   
    %update particle position 
    z(:,1)=z(:,1) + dz.'; % randomization precess is not saved
    y(:,1)=y(:,1) + dy.';
   
    %reset particle position if they leave the restricted area 
    z((z(:,1)>zmax),1) =  z((z(:,1)>zmax),1)-2*zmax;
    z((z(:,1)<zmin),1) =  z((z(:,1)<zmin),1)+2*zmax;
    y((y(:,1)>ymax),1) =  y((y(:,1)>ymax),1)-2*ymax;
    y((y(:,1)<ymin),1) =  y((y(:,1)<ymin),1)+2*ymax;
    
end










%Here starts the actual simulation 
tic % measure time 

plotlive=301; % some variable to plot live during the simulation
%iterate over all time stepssteps
for i=1:steps 
    %define some convenient variables
    zold=1*z(:,i)+0;
    yold=1*y(:,i)+0 ;                    % extend matrix
    plotlive=plotlive+1;
    
    % live plotting
     if plotlive>300
         figure(1);
         plot(zold,yold,'o')
         
         plotlive=0;
     end
      
    % set the drift to zero for this time step
    dy=zeros(1,particles);   
    dz=dy;
    % make periodic boundaries i.e. copy the array around itselfe
     zold2=[zold  zold          zold                 (zold+2*zmax)   (zold+2*zmin)     (zold+2*zmax) (zold+2*zmin) (zold+2*zmax) (zold+2*zmin)];            % deinfe old positions 
     yold2=[yold (yold+2*ymax)  (yold+2*ymin)        yold                yold           (yold+2*ymax) (yold+2*ymin) (yold+2*ymin) (yold+2*ymax)];
   
   %iterate over all particles 
   for j=1:particles 
             % compute the drift for one partikel  j
             [dz(j),dy(j)] =Dipolhertzquadrandsize(zold2(:),yold2(:),transpose(betas(j,:)),alphas(:),transpose(kas(j,:)),R0(:),j);
                 
   end
   % add a random number vector to the drift vector
   dz=dz*dt+sigmas.*randn(1,particles) ;
   dy=dy*dt+sigmas.*randn(1,particles) ;
   
    %update particle position 
    z(:,i+1)=z(:,i) + dz.';
    y(:,i+1)=y(:,i) + dy.';
   
    %reset particle position if they leave the restricted area 
    z((z(:,i+1)>zmax),i+1) =  z((z(:,i+1)>zmax),i+1)-2*zmax;
    z((z(:,i+1)<zmin),i+1) =  z((z(:,i+1)<zmin),i+1)+2*zmax;
    y((y(:,i+1)>ymax),i+1) =  y((y(:,i+1)>ymax),i+1)-2*ymax;
    y((y(:,i+1)<ymin),i+1) =  y((y(:,i+1)<ymin),i+1)+2*ymax;
    
 end
toc % print time 

% save the data as txt files
zfilnametxt = sprintf("39hoch2_k%f_z_gamma%f_Kdhoch2%f.txt",kdimensionfactor, gamma, Kdhoch2);
writematrix(z,zfilnametxt ,'Delimiter','tab');
yfilnametxt = sprintf("39hoch2_k%f_y_gamma%f_Kdhoch2%f.txt",kdimensionfactor, gamma, Kdhoch2);
writematrix(z,yfilnametxt ,'Delimiter','tab');
%metafilnametxt = sprintf("meta_a%f_b%f.txt", alpha, beta);









% Plotting of the particles. 

b = 0; % some magic constant that I need for reasons


% define the plotting window
set(gca,'units','pixels','position', [200 200 750 750])
set(gcf, 'units', 'centimeters', 'position', [0 0 50 50])
set(gca,'XColor', 'none','YColor','none')

% define last image of the plot
steps=i;

% iterate and plot every hop steps
for i=1:hop:steps   % iterate over  time steps
   cla
   b = b +1;
for j=1:particles  % iterate over all particles

    
% plot a circle for every particle   
xlim([zmin zmax]);
ylim([ymin ymax]);
h(j)=circle(y(j,i),z(j,i),R0(j));
hold on 

end
set(gca,'XColor', 'none','YColor','none') % remove the axis

%save the frame
N(b)= getframe();
delete(h)



end


%create a video file
filname = sprintf("particles%f²frame%fs_actualsize_%f_k%f_gamma%f_Kdhoch2%f.avi",sqrtparticles,frame,actualsize,kdimensionfactor, gamma, Kdhoch2);
% write everting into an avi file
v = VideoWriter(filname);
v.Quality = 100; % this is very important 
open(v)
writeVideo(v,N)
close(v)
