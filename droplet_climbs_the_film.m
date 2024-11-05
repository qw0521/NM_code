%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% This code was written by Zhili Hu and adjusted by Wen Qin, both works in NUAA in 2024.
%% Email: zhili.hu@nuaa.edu.cn
%%
%% This code is based on the phase-field model by Borcia, R. et al. and the original paper is,
%% Static and dynamic contact angles – A phase field modelling. Eur. Phys. J. Spec. Top. 166, 127-131 (2009).
%%
%% We use it to simulate the equilibrium shape profile of a nanodroplet that sits
%% at an atomic step. The droplet can either sit still at the step or climb up the
%% step, depending on its relative size and wetabilities. The exampled parameter
%% set is for a Pt droplet that sits at the edge of a PtSex monolayer, on a MoS2
%% substrate. Users could modify the code for other scensrios.
%%
%% Details of the model can be found in the Methods section in our paper in ?? and
%% a pseudo code is available in the supplementary information.
%%
%% For an easy usage, lines that potential users may need to change is marked with “USER-ATTENTION”
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clear;                          % clear history data

%% USER-ATTENTION: parameters that users may set
Lx=2^6;Ly=2^5;Lz=2^5;           % the X, Y, and Z sizes of the simulation box.
R=10;                     	    % the radius of the Pt droplet. The smaller droplet will stay at step, while larger one will climb the step.
file_name="Climb_films";        % the name for saving simulation results.
lam=180;                        % the coupling constant.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=0.5;                          % the space step size along either X, Y, and Z.
dt=c^2*0.01;                    % the time step
Nx=round(Lx/c);Ny=round(Ly/c);Nz=round(Lz/c); % "Nx", "Ny", and "Nz" are the length of the simulation area regulated by "c".
xi=zeros(Ny,Nx,Nz);             % "the dimensionless supersaturation constant ξ, defind as (vol0-vol)/vol_vapor.
[x y z]=meshgrid([1:Nx]*c-Lx/2,[1:Ny]*c-Ly/2,[1:Nz]*c);
r2=(x+5).^2+y.^2+(z-5).^2;      % USER-ATTENTION: The Pt droplet is initially set at PtSex step.
% r2=(x-15).^2+y.^2+(z-8).^2;   % The Pt droplet is initially set on PtSex top.
rou=sign(R^2-r2)/2+.5;          % the phase-field for the density of the liquid.
                                % rou=1 denotes liquid phase and rou=0 denotes the vapor phase.
K=1;                            % the droplet’s surface tension
mask=1-(z<(atan(10*x)+pi/2)/pi*3+3); % the region allowed for all Pt, either as vapor or droplet
%% USER-ATTENTION: these 2 lines controls the contact angles of the droplet
wall=(1-mask).*(.6323+z*.0);        % We set rou_s=0.6323 such that the contact angle between Pt/PtSex is 67.2
wall=wall-.1259*(wall>0).*(z<3+1);  % We set rou_s=0.5064 (0.6323-0.1259) such that the contact angle between Pt/MoS2 is 88.9
%%
vols=zeros(5,1);                % the volume of the Pt droplet.
rou=rou.*mask;
vol0=sum(rou(:));               % the initial volume of the Pt droplet.
% Initialize video file.
writerObj=VideoWriter(['video_' num2str(file_name) '.avi']);
writerObj.FrameRate = 60;
writerObj.Quality = 95;
open(writerObj);

for t=1:1500/dt % USER-ATTENTION: The total evolution time is set here.
    v=K*lap3(rou,0.5)-rou.*(rou-1).*(2*rou-1)*10+lam*xi.*(rou-rou.^2); % v means d rou/d t
    tf=1;
    v=v./tf.*mask;
    rou=rou+dt*v;
    rou=rou.*mask;
    rou=rou+wall;
    vol=sum(rou(:).*mask(:));           % the condensed volume of all Pt either as vapor or in the droplet.
    xi=(vol0-vol)/sum(mask(:))*mask;    % the dimensionless supersaturation
    if (t-0)/100>0&&mod(t-100,100)==0 %USER-ATTENTION: This controls at which time steps, graphs are plotted.
    % This block is quite lengthy for a better rendering of the images.
    % Users could simply ignore the detailed logic here, if only interested in the physics
        p=rou.*mask;
        ws=sum(1-mask,3);
        for i=1:Ny % Without this block, the plotted droplet will not be very spheric
            for j=1:Nx
                b=false;
                for k=Nz:-1:1
                    if p(i,j,k)>.99
                        for k1=k:-1:1
                            p(i,j,k1)=1;
                        end
                    break
                    end
                end
            end
        end
        pp=(sum(p.*mask,3)+ws);
        p1=pp;
        xx=[1:1:64];yy=[1:1:128];
        xlin = linspace(1,64,640);ylin = linspace(1,128,1280);
        [X,Y] = meshgrid(ylin,xlin);
        Z = interp2(yy,xx,p1,X,Y,'cubic')*10;
        C=sum(p.*mask,3);
        C = interp2(yy,xx,C,X,Y,'cubic');C(C>=2.5)=30;C(C<2.5)=0;
        WS= interp2(yy,xx,ws,X,Y,'cubic');
        C=C+WS;C(C>=12)=30;
        if true
            % video setups
            surf(Z,C,'LineStyle','none','FaceColor','interp','FaceLighting','gouraud')
            axis off;
            axis equal;
            material default
            title(['t=',int2str(round((t-0)/100))])
            light('Position',[30 10 50])
            view(-20,35)
            drawnow;
        end
        frame = getframe(gcf);  % Save the image into the video file.
        writeVideo(writerObj,frame);
        if mod(t-1,round(10/dt))==0
            n=round((t-1)/10*dt)+1;
            vols(round((t-1)/10*dt)+1)=vol;
            toc % This prints the time spent for the simulation
        end
    end
end
close(writerObj);
filenm=['Rt',num2str(file_name),'.mat'];
save (filenm)
saveas(gcf,strcat('Rt',num2str(file_name),'-',num2str(n)),'jpg') % The video is saved.

function y= lap3(x,h) % 3d Laplacian operator
s=size(x);if length(s)<3,s(3)=1;end
a=[1:s(1)];b=[1:s(2)];c=[1:s(3)];
y=(x(mod(a-2,s(1))+1,:,:)+x(mod(a,s(1))+1,:,:)...
    +x(:,mod(b-2,s(2))+1,:)+x(:,mod(b,s(2))+1,:)...
    +x(:,:,mod(c-2,s(3))+1)+x(:,:,mod(c,s(3))+1)...
    -6*x)/h^2;
end

