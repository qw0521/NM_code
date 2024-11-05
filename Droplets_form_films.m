%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% This code was written by Zhili Hu and adjusted by Wen Qin, both works in NUAA in 2024.
%% Email: zhili.hu@nuaa.edu.cn
%%
%% This code is based on our previous code written for paper He, Y. et al. Nat. Commun. 11, 57 (2020).
%% and the phase field model in it is based on Karma, Phys. Rev. Lett. 81, 4444-4447 (1998).
%%
%% We use it to simulate the morphological evolution in our experiments, where
%% Pt nanodroplets drive the formation of earthworm like PtSex nanoribbons, which
%% further stitch into complete films. Users could modify this code for their own purposes.
%%
%% Details of the model can be found in the Methods section in our paper.
%%
%% For an easy usage, lines that potential users may want to change is marked with “USER-ATTENTION”
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic                             %
clear;                          % clear history data

%% USER-ATTENTION: parameters that users may set
Ndroplet=6;                     % Pt droplet #
L=170;                          % simulation box size
c=1;                            % space grid size
dt=c^2/100;                     % time step
N=L/c;                          % the length of the data matrix
thickness=0.011;                % the dimentionless thickness of a PtSex monolayer,
                                % which governs the consumption rate of Pt droplets during motion.
Rball=15;                       % the initial radius of Pt droplets
Rc=35;                          % the critical radius above which Pt droplets solidify hence cease to move
file_name="Form_films";         % the name for saving simulation results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y]=meshgrid([1:N]*c-L/2,[1:N]*c-L/2);
r2=@(x,y)((ones(N,1)*[1:N]-N/2)*c-x).^2+((ones(N,1)*[1:N]-N/2)'*c-y).^2;
phi=-ones(N);                   % the order parameter, whereas  -1 represents the
                                % bare MoS2 substrate, and 1 represents the PtSex single layer.
lambda=10;                      % a dimensionless coupling constant
tau_phi=1;                      %
pbcdist=@(x1,x2)min(abs(x1-x2),L-abs(x1-x2)); % distance along x direction due to periodic doundary condition

%% This block defines functions that use FFT to solve derivatives
q=[0:N/2 (N/2+1:N-1)-N]*2*pi/L;
q2=(ones(N,1)*q.^2+(ones(N,1)*q.^2)');
q=(ones(N,1)*[0:N/2-1 0 (N/2+1:N-1)-N]*2*pi/L)';
nabx=@(f)ifft(fft(f')*1i.*q)'; % x component of a nabla operator
naby=@(f)ifft(fft(f)*1i.*q); % y component
lap=@(f)ifft2(-q2.*fft2(f)); % Laplacian operator
%%%%%%%%%%%%%

r21=@(x1,y1)(pbcdist(x1,x).^2+pbcdist(y1,y).^2);    % distance^2
x2i=@(x)mod(round(x/c+N/2)-1,N)+1;  % convert coordinates of x,y to an index in array droplets

while true                      % initialize the droplets with random positions
  droplets=rand(Ndroplet,3)*L-L/2; % droplets(i,1:2) stores x and y of droplet i
  droplets(:,3)=Rball^3; % droplets(i,3) is proportional to droplet i's volumn
  overlap=0;
  for i=1:Ndroplet % these lines prevent droplets from overlapping
    for j=i+1:Ndroplet
      if (pbcdist(droplets(i,1),droplets(j,1))^2+pbcdist(droplets(i,2),droplets(j,2))^2)^.5...
        <(droplets(i,3)^(1/3)+droplets(j,3)^(1/3))*2
        overlap=1;
      end
    end
  end
  if overlap==0
    break
  end
end
for i=1:Ndroplet                   % PtSex film nucleates at each Pt droplet
  s=rand(1,2)-.5;
  s=s/norm(s)*.1*Rball+droplets(i,1:2);
  phi=phi+sign(Rball^2-r2(s(1),s(2)))+1;
end

%% these 3 lines are for figures
h=imagesc(flipud(phi));
axis equal;axis off;
scale=5;N1=scale*N; % "scale" defines the figure quality

phis=zeros(N1,N1,10);% "phis" stores phi as a function of time
Cs=phis;phis3=phis;

for t=1:600/dt % the morphological evolution of droplets and PtSex
  b=zeros(N);
  for i=1:Ndroplet                 % b represents area covered by the droplets.
    b=b+sign(droplets(i,3)^(2/3)-r21(droplets(i,1),droplets(i,2)))/2+.5;
  end
  xi=0.02*ones(N)+.3*b;            % xi is the concentration field of Se species
  v=2*real(lap(phi))+sin(pi*phi)+lambda*xi.*(1+cos(phi*pi)); % v is proportional to d phi/ d t
  v=real(v).*(.9*b+.1);
  phi=phi+dt/tau_phi*v;
  for i=1:Ndroplet                 % Pt droplets are consumed by the growth of PtSex
    if droplets(i,3)<=0
      continue
    end
    b=sign(droplets(i,3)^(2/3)-r21(droplets(i,1),droplets(i,2)))/2+.5;
    v1=v.*b;
    droplets(i,3)=droplets(i,3)-sum(v1(:))/2*thickness;
  end
  for i=1:Ndroplet                 % Simulate each Pt droplet moving.
    if droplets(i,3)>Rc^3||droplets(i,3)<=0
      continue
    end
    droplets(i,1:2)=droplets(i,1:2)+(rand(1,2)-.5)*80/droplets(i,3)^0.5; % the Browniam motion of droplets.
    phix=real(nabx(phi));
    phiy=real(naby(phi));
    ix0=x2i(droplets(i,1));
    iy0=x2i(droplets(i,2));
    dir=[phix(iy0,ix0),phiy(iy0,ix0)];
    dir=real(dir/norm(dir));
    slope=[];
    ix1=mod(round(ix0+dir(1)/1*[-Rball:Rball])-1,N)+1;
    iy1=mod(round(iy0+dir(2)/1*[-Rball:Rball])-1,N)+1;
    for j=1:length(ix1)
      slope=[slope phi(iy1(j),ix1(j))];
    end
    base=round(min(slope));
    [var,vj]=min(abs(slope-base-1));
    if abs(vj-Rball-1)>20
      vj=Rball+1+2*sign(vj-Rball-1);
    end
    ix=ix1(vj);iy=iy1(vj);
    v_xy=[x(iy,ix) y(iy,ix)]-droplets(i,1:2);
    droplets(i,1:2)=droplets(i,1:2)+v_xy;
    % once a droplet stands on a complete layer, it lays a new nuclus under it.
    b=sign(droplets(i,3)^(2/3)-r21(droplets(i,1),droplets(i,2)))/2+.5;
    phib=(phi+1).*b;area=real(sum(b(:)));phib=sort(phib(:));phib=phib(end-round(area)+1:end);
    if imag(area)~=0
      pause;
    end
    if max(phib)-min(phib)<.1&&droplets(i,3)>300 %tweak norm(dir)<.0001
      s=rand(1,2);
      r=droplets(i,3)^(1/3);
      s=s/norm(s)*.5*r+droplets(i,1:2);
      v1=sign(r^2-r2(s(1),s(2)))+1;
      phi=phi+v1;
      droplets(i,3)=droplets(i,3)-pi*r^2*thickness;
    end
  end
  for i=1:Ndroplet-1       % Simulation of the merging of Pt droplets
    for j=i+1:Ndroplet
      if (pbcdist(droplets(i,1),droplets(j,1))^2+pbcdist(droplets(i,2),droplets(j,2))^2)^.5...
        <(droplets(i,3)^(1/3)+droplets(j,3)^(1/3))*1
        xij=droplets(j,1)-droplets(i,1);
        yij=droplets(j,2)-droplets(i,2);
        xij=mod(xij+L/2,L)-L/2;
        yij=mod(yij+L/2,L)-L/2;
        % the new center is the mass center of the previous two droplets
        droplets(i,1:2)=droplets(i,1:2)+([xij yij]*droplets(j,3))/(droplets(i,3)+droplets(j,3));
        droplets(i,3)=real(droplets(i,3)+droplets(j,3)); % the volumns merge
        droplets(j,3)=0;
      end
    end
  end
  interval=10; % the interval for ploting
  if mod(t-1,interval)==0 % Plotting
    p=real(phi);
    p(1,1:2)=[0 10];
    % we scale up the canvas from var p to p1 by 2D fitting, for a finer image resolution
    % the below line contains function of 2d fit that doesnot work in Octave, so Octave users
    % could comment next a few lines relating to 2D fitting and set p1 =p after them
    surffit = fit([x(:),y(:)],p(:),'cubicinterp','normalize','on');
    [x1 y1]=meshgrid([1:N1]*c/scale-L/2,[1:N1]*c/scale-L/2);
    r22=@(x,y)(pbcdist(x1,x).^2+pbcdist(y1,y).^2);
    if true
        p1=reshape(surffit(x1(:),y1(:)),[N1,N1]);
        C=p1;
        p2=(p1>0)*3;
        for i=1:Ndroplet
            p3=max(droplets(i,3)^(2/3)-r22(droplets(i,1),droplets(i,2)),0).^.5;
            p2=max(p2,p3);
            C=max(C,(p3>0)*9.5);
        end
        surf(x1,y1,real(p2),real(C),'LineStyle','none')
        axis equal;axis off;
        light('Position',[30 10 50])
        phis3(:,:,(t-1)/interval+1)=p2;
        Cs(:,:,(t-1)/interval+1)=C;
    end
    coverage1=numel(find(phi>0.99));coverage2=numel(find(phi>1.99));rate1=coverage1/L^2;rate2=coverage2/L^2;
    title(['rate1=',string(rate1),'rate2=',string(rate2)])
    drawnow
  end
  if mod(t-1,100)==0 % once the whole substrate is covered by PtSex, the code stops
    ['t=' int2str(t*dt) '#' int2str(toc) 's'];
    if sum(phi<.5)==0
      break
    end
    toc
  end
  if sum(droplets(:,3))<=0
      break
  end
end
%% Making video file.
close;
t=0:1:size(phis3,3)-1; % Set video time length and step.
% Initialize video file.
writerObj=VideoWriter(['video_' num2str(file_name) '.avi']);
writerObj.FrameRate = 60;
writerObj.Quality = 95;
open(writerObj);
for i=1:length(t)
    surf(x1,y1,real(phis3(:,:,i)),real(Cs(:,:,i)),'LineStyle','none')
    axis equal;axis off;light('Position',[30 10 50])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end
close(writerObj);