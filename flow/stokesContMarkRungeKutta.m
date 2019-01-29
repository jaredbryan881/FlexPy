% Modified from Taras Gerya's 'Numerical Geodynamics,' 2010.

% Solution of 2D Stokes continuity and advection equations 
% with finite differences and marker-in-cell technique
% on a regular grid using pressure-velocity formulation
% for a medium with variable viscosity
% using of the forth-order accurate in space 
% first order accurate in time 
% Runge-Kutta marker advection scheme

% Clean all variables
clear;
% Clear all figures
clf;

% Numerical model parameters
% Model size, m
xsize   =   1500000; % Horizontal
%ysize   =   1200000; % Vertical
ysize   =   500000;

% Numbers of nodes
xnum    =   51; % Horizontal
%ynum    =   61; % Vertical
ynum    =   31;
% Grid step
xstp    =   xsize/(xnum-1); % Horizontal
ystp    =   ysize/(ynum-1); % Vertical

% Pressure condition in one cell (i==2 && j==3)
p0cell  =   0;

% Gravity acceleration directed downward
gy      =   9.81; % m/s^2

% Making vectors for nodal points positions (basic nodes)
x       =   0:xstp:xsize; % Horizontal
y       =   0:ystp:ysize; % Vertical

% Making vectors for cell centers positions (staggered nodes)
xc      =    xstp/2:xstp:xsize-xstp/2; % Horizontal
yc      =    ystp/2:ystp:ysize-ystp/2; % Vertical

% Maximal timestep, s
timemax=1e+4*365.25*24*3600; %10 kyr
% Maximal marker displacement step, number of gridsteps
markmax=0.5;
% Amount of timesteps
stepmax=100;
% Moving Markers: 
% 0 = not moving at all
% 1 = simple 1-st order advection
% 4 = 4-th order in space  Runge-Kutta
markmove=4;

% Material viscosity, Pa s
etaref=1e+22;
MFLOW(1)=50*etaref;
MFLOW(2)=0.1*etaref;
MFLOW(3)=etaref;
MFLOW(4)=50*etaref;
MFLOW(5)=1e+5*etaref;
MFLOW(6)=1e+13; %Free Surface

% Material density, kg/m^3
MRHO(1)=2900; 
MRHO(2)=3300; 
MRHO(3)=3500;
MRHO(4)=4400;
MRHO(5)=4000; % BLOCK
MRHO(6)=1; % Free Surface

freeSurf = true;
if freeSurf == true
    waterdepth=60000;
else
    waterdepth=0;
end

%waterdepth=0;
blockleftbound=(xsize/2)-40000;
blockrightbound=(xsize/2)+40000;
blockupperbound=160000;
blocklowerbound=220000;
% Making vectors for nodal points positions (basic nodes)
gridx=0:xstp:xsize; % Horizontal
gridy=0:ystp:ysize; % Vertical
% Making vectors for cell centers positions (staggered nodes)
gridcx=xstp/2:xstp:xsize-xstp/2; % Horizontal
gridcy=ystp/2:ystp:ysize-ystp/2; % Vertical

% Defining number of markers and steps between them in the horizontal and vertical direction
mxnum=200; %total number of markers in horizontal direction
mynum=300;  %total number of markers in vertical direction
mxstep=xsize/mxnum; %step between markers in horizontal direction   
mystep=ysize/mynum; %step between markers in vertical direction

% Creating markers arrays
MX=zeros(mynum*mxnum,1);   % X coordinate, m
MY=zeros(mynum*mxnum,1);   % Y coordinate, m
MI=zeros(mynum*mxnum,1);   % Type

% Defining intial position of markers
% Defining lithological structure of the model
% Marker counter
mm1=0;
counter=1;
for xm = 1:1:mxnum
    for ym = 1:1:mynum
        % Update marker counter:
        mm1=mm1+1;
        
        % Coordinates with small random displacement
        MX(mm1)=xm*mxstep-mxstep/2+(rand-0.5)*mxstep;
        MY(mm1)=ym*mystep-mystep/2+(rand-0.5)*mystep;
        % Defining initial rock distribution (marker type)
        if((MY(mm1)<blocklowerbound+waterdepth) && (MY(mm1)>blockupperbound+waterdepth) && (MX(mm1)>blockleftbound) && (MX(mm1)<blockrightbound))
            MI(mm1)=5;
        else
            if(MY(mm1)<=100000+waterdepth)
                MI(mm1)=1;
            elseif((MY(mm1)>100000+waterdepth) && (MY(mm1)<=300000+waterdepth))
                MI(mm1)=2;
            elseif((MY(mm1)>300000+waterdepth) && (MY(mm1)<=660000+waterdepth))
                MI(mm1)=2; %3
            else
                MI(mm1)=2; %4
            end
            
            if(MY(mm1)<=waterdepth)
                MI(mm1)=6;
            end
        end
    end
end

% Save Number of markers
marknum=mm1;

% Density, viscosity arrays
etas = zeros(ynum,xnum);    % Viscosity for shear stress
etan = zeros(ynum,xnum);    % Viscosity for normal stress
rho = zeros(ynum,xnum);     % Density

% keep track of the mm1 coords for the lithosphere markers
lith = find(MI==1);

% find the markers at the 'surface' of the lithosphere. 
% really the minimum depth for makers MI=1 for a given x
surfInd = zeros(mxnum, 1);

for s=1:mxnum
    surfInd(s) = lith(1 + (60 * (s-1)));
end

% free surface coordinates for each time step
surfCoords = zeros(stepmax, mxnum, 2);

% Initial time, s
timesum=0;
% Main Time cycle
for ntimestep=1:1:stepmax
    fprintf('Time step %d \n',ntimestep)
    tic();
    
    % Backup transport properties arrays
    etas0 = etas;
    etan0 = etan;
    rho0 = rho;
    % Clear transport properties arrays
    etas = zeros(ynum,xnum);   % Viscosity for shear stress
    etan = zeros(ynum,xnum);   % Viscosity for normal stress
    rho = zeros(ynum,xnum);    % Density
    % Clear wights for basic nodes
    wtnodes=zeros(ynum,xnum);
    % Clear wights for etas
    wtetas=zeros(ynum,xnum);
    % Clear wights for etan
    wtetan=zeros(ynum,xnum);

    % Interpolating parameters from markers to nodes
    for mm1 = 1:1:marknum

        % Check markers inside the grid
        if (MX(mm1)>=gridx(1) && MX(mm1)<=gridx(xnum) && MY(mm1)>=gridy(1) && MY(mm1)<=gridy(ynum)) 

            %  xn    rho(xn,yn)--------------------rho(xn+1,yn)
            %           ?           ^                  ?
            %           ?           ?                  ?
            %           ?          dy                  ?
            %           ?           ?                  ?
            %           ?           v                  ?
            %           ?<----dx--->o Mrho(xm,ym)       ?
            %           ?                              ?
            %           ?                              ?
            %  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
            %
            %
            % Define indexes for upper left node in the cell where the marker is
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16(MX(mm1)./xstp-0.5))+1;
            yn=double(int16(MY(mm1)./ystp-0.5))+1;
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-1)
                xn=xnum-1;
            end
            if (yn<1)
                yn=1;
            end
            if (yn>ynum-1)
                yn=ynum-1;
            end

            % Define normalized distances from marker to the upper left node;
            dx=(MX(mm1)-gridx(xn))./xstp;
            dy=(MY(mm1)-gridy(yn))./ystp;

            % Define material density and viscosity from marker type
            MRHOCUR=MRHO(MI(mm1)); % Density
            METACUR=MFLOW(MI(mm1)); % Viscosity


            % Add density to 4 surrounding basic nodes
            % Upper-Left node
            rho(yn,xn)=rho(yn,xn)+(1.0-dx).*(1.0-dy).*MRHOCUR;
            wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy);
            % Lower-Left node
            rho(yn+1,xn)=rho(yn+1,xn)+(1.0-dx).*dy.*MRHOCUR;
            wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy;
            % Upper-Right node
            rho(yn,xn+1)=rho(yn,xn+1)+dx.*(1.0-dy).*MRHOCUR;
            wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy);
            % Lower-Right node
            rho(yn+1,xn+1)=rho(yn+1,xn+1)+dx.*dy.*MRHOCUR;
            wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy;

            % Add viscosity etas() to 4 surrounding basic nodes
            % only using markers located at <=0.5 gridstep distances from nodes
            % Upper-Left node
            if(dx<=0.5 && dy<=0.5)
                etas(yn,xn)=etas(yn,xn)+(1.0-dx).*(1.0-dy).*METACUR;
                wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx).*(1.0-dy);
            end
            % Lower-Left node
            if(dx<=0.5 && dy>=0.5)
                etas(yn+1,xn)=etas(yn+1,xn)+(1.0-dx).*dy.*METACUR;
                wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx).*dy;
            end
            % Upper-Right node
            if(dx>=0.5 && dy<=0.5)
                etas(yn,xn+1)=etas(yn,xn+1)+dx.*(1.0-dy).*METACUR;
                wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx.*(1.0-dy);
            end
            % Lower-Right node
            if(dx>=0.5 && dy>=0.5)
                etas(yn+1,xn+1)=etas(yn+1,xn+1)+dx.*dy.*METACUR;
                wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx.*dy;
            end

            % Add viscosity etan() to the center of current cell (pressure node)
            etan(yn+1,xn+1)=etan(yn+1,xn+1)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy)).*METACUR;
            wtetan(yn+1,xn+1)=wtetan(yn+1,xn+1)+(1.0-abs(0.5-dx)).*(1.0-abs(0.5-dy));
        end
    end

    % Computing  Viscosity, density, rock type for nodal points
    for i=1:1:ynum
        for j=1:1:xnum
            % Density
            if (wtnodes(i,j)~=0)
                % Compute new value interpolated from markers
                rho(i,j)=rho(i,j)./wtnodes(i,j);
            else
                % If no new value is interpolated from markers old value is used
                rho(i,j)=rho0(i,j);
            end
            % Viscosity etas() (basic nodes)
            if (wtetas(i,j)~=0)
                % Compute new value interpolated from markers
                etas(i,j)=etas(i,j)./wtetas(i,j);
            else
                % If no new value is interpolated from markers old value is used
                etas(i,j)=etas0(i,j);
            end
            % Viscosity etan() (pressure cells)
             if (wtetan(i,j)~=0)
                % Compute new value interpolated from markers
                etan(i,j)=etan(i,j)./wtetan(i,j);
            else
                % If no new value is interpolated from markers old value is used
                etan(i,j)=etan0(i,j);
            end
        end
    end

    % Matrix of coefficients initialization
    L=sparse(xnum*ynum*3,xnum*ynum*3);
    % Vector of right part initialization
    R=zeros(xnum*ynum*3,1);

    % Computing Kcont and Kbond coefficients 
    etamin=min(min(etas)); % min viscosity in the model
    kcont=2*etamin/(xstp+ystp);
    kbond=4*etamin/(xstp+ystp)^2;

    % Solving x-Stokes, y-Stokes and continuity equations
    % x-Stokes: ETA(d2vx/dx2+d2vx/dy2)-dP/dx=0
    % y-Stokes: ETA(d2vy/dx2+d2vy/dy2)-dP/dy=gy*RHO
    % continuity: dvx/dx+dvy/dy=0
    % Composing matrix of coefficients L()
    % and vector (column) of right parts R()
    % Boundary conditions: free slip
    % Process all Grid points
    for i=1:1:ynum
      for j=1:1:xnum
        % Global index for P, vx, vy in the current node
        inp=((j-1)*ynum+i)*3-2; % P
        invx=inp+1;
        invy=inp+2;


        % Continuity equation
        % Ghost pressure unknowns (i=1, j=1) and boundary nodes (4 corners + one cell)
        if(i==1 || j==1 || (i==2 && j==2) || (i==2 && j==xnum) || (i==ynum && j==2) || (i==ynum && j==xnum) || (i==2 && j==3))
            % Ghost pressure unknowns (i=1, j=1): P(i,j)=0
            if(i==1 || j==1)
                L(inp,inp)          =   1*kbond;            % Coefficient for P(i,j)
                R(inp,1)            =   0;                  % Right part
            end
            % Upper and lower left corners dP/dx=0 => P(i,j)-P(i,j+1)=0
            if((i==2 && j==2) || (i==ynum && j==2))
                L(inp,inp)          =   1*kbond;            % Coefficient for P(i,j) 
                L(inp,inp+ynum*3)   =   -1*kbond;           % Coefficient for P(i,j+1)
                R(inp,1)            =   0;                  % Right part
            end
            % Upper and lower right corners dP/dx=0 => P(i,j)-P(i,j-1)=0
            if((i==2 && j==xnum) || (i==ynum && j==xnum))
                L(inp,inp)          =   1*kbond;            % Coefficient for P(i,j) 
                L(inp,inp-ynum*3)   =   -1*kbond;           % Coefficient for P(i,j-1)
                R(inp,1)            =   0;                  % Right part
            end
            % One cell 
            if (i==2 && j==3)
                L(inp,inp)          =   1*kbond;            % Coefficient for P(i,j)
                R(inp,1)            =   p0cell;             % Right part
            end
        %Internal nodes: dvx/dx+dvy/dy=0
        else
            %dvx/dx=(vx(i-1,j)-vx(i-1,j-1))/dx
            L(inp,invx-3)           =   kcont/xstp;         % Coefficient for vx(i-1,j) 
            L(inp,invx-3-ynum*3)    =   -kcont/xstp;        % Coefficient for vx(i-1,j-1) 
            %dvy/dy=(vy(i,j-1)-vy(i-1,j-1))/dy
            L(inp,invy-ynum*3)      =   kcont/ystp;         % Coefficient for vy(i,j-1) 
            L(inp,invy-3-ynum*3)    =   -kcont/ystp;        % Coefficient for vy(i-1,j-1) 
            % Right part:0
            R(inp,1)=0;
        end

        % x-Stokes equation
        % Ghost vx unknowns (i=ynum) and boundary nodes (i=1, i=ynum-1, j=1, j=xnum)
        if(i==1 || i==ynum-1 || i==ynum || j==1 || j==xnum)
            % Ghost vx unknowns (i=ynum: vx(i,j)=0
            if(i==ynum)
                L(invx,invx)        =   1*kbond; % Coefficient for vx(i,j)
                R(invx,1)           =   0; % Right part
            end
            % Left and Right boundaries (j=1, j=xnum) 
            if((j==1 || j==xnum) && i<ynum)
                % Free slip, No slip: vx(i,j)=0
                L(invx,invx)        =   1*kbond; % Coefficient for vx(i,j)
                R(invx,1)           =   0; % Right part
            end
            % Upper boundary, iner points (i=1, 1<j<xnum)
            if(i==1 && j>1 && j<xnum)
                % Free slip dvx/dy=0: vx(i,j)-vx(i+1,j)=0
                L(invx,invx)        =   1*kbond; % Coefficient for vx(i,j)
                L(invx,invx+3)      =   -1*kbond; % Coefficient for vx(i+1,j)
                R(invx,1)=0; % Right part
    %             % No slip vx=0: vx(i,j)-1/3*vx(i+1,j)=0
    %             L(invx,invx)=1*kbond; % Coefficient for vx(i,j)
    %             L(invx,invx+3)=-1/3*kbond; % Coefficient for vx(i+1,j)
    %             R(invx,1)=0; % Right part
            end
            % Lower boundary, iner points (i=ynum-1, 1<j<xnum)
            if(i==ynum-1 && j>1 && j<xnum)
                % Free slip dvx/dy=0: vx(i,j)-vx(i-1,j)=0
                L(invx,invx)        =   1*kbond; % Coefficient for vx(i,j)
                L(invx,invx-3)      =   -1*kbond; % Coefficient for vx(i-1,j)
                R(invx,1)=0; % Right part
    %             % No slip vx=0: vx(i,j)-1/3*vx(i+1,j)=0
    %             L(invx,invx)=1*kbond; % Coefficient for vx(i,j)
    %             L(invx,invx-3)=-1/3*kbond; % Coefficient for vx(i-1,j)
    %             R(invx,1)=0; % Right part
            end
        %Internal nodes: dSxx/dx+dSxy/dy-dP/dx=0
        else
            %dSxx/dx=2*etan(i+1,j+1)*(vx(i,j+1)-vx(i,j))/dx^2-2*etan(i+1,j)*(vx(i,j)-vx(i,j-1))/dx^2
            L(invx,invx+ynum*3)     =   2*etan(i+1,j+1)/xstp^2;                         % Coefficient for vx(i,j+1)
            L(invx,invx-ynum*3)     =   2*etan(i+1,j)/xstp^2;                           % Coefficient for vx(i,j-1)
            L(invx,invx)            =   -2*etan(i+1,j+1)/xstp^2-2*etan(i+1,j)/xstp^2;   % Coefficient for vx(i,j)

            %dSxy/dy=etas(i+1,j)*((vx(i+1,j)-vx(i,j))/dy^2+(vy(i+1,j)-vy(i+1,j-1))/dx/dy)-
            %         -etas(i,j)*((vx(i,j)-vx(i-1,j))/dy^2+(vy(i,j)-vy(i,j-1))/dx/dy)-
            L(invx,invx+3)          =   etas(i+1,j)/ystp^2;                             % Coefficient for vx(i+1,j)
            L(invx,invx-3)          =   etas(i,j)/ystp^2;                               % Coefficient for vx(i-1,j)
            L(invx,invx)            =   L(invx,invx)-etas(i+1,j)/ystp^2-etas(i,j)/ystp^2; % ADD coefficient for vx(i,j)
            L(invx,invy+3)          =   etas(i+1,j)/xstp/ystp;                          % Coefficient for vy(i+1,j)
            L(invx,invy+3-ynum*3)   =   -etas(i+1,j)/xstp/ystp;                         % Coefficient for vy(i+1,j-1)
            L(invx,invy)            =   -etas(i,j)/xstp/ystp;                           % Coefficient for vy(i,j)
            L(invx,invy-ynum*3)     =   etas(i,j)/xstp/ystp;                            % Coefficient for vy(i,j-1)
            % -dP/dx=(P(i+1,j)-P(i+1,j+1))/dx
            L(invx,inp+3)           =   kcont/xstp;                                     % Coefficient for P(i+1,j)
            L(invx,inp+3+ynum*3)    =   -kcont/xstp;                                    % Coefficient for P(i+1,j+1)
            % Right part:0
            R(invx,1)               =   0;
        end

        % y-Stokes equation
        % Ghost vy unknowns (j=xnum) and boundary nodes (i=1, i=ynum, j=1, j=xnum-1)
        if(i==1 || i==ynum || j==1 || j==xnum-1 || j==xnum)
            % Ghost vy unknowns (j=xnum: vy(i,j)=0
            if(j==xnum)
                L(invy,invy)        =   1*kbond;                                % Coefficient for vy(i,j)
                R(invy,1)           =   0;
            end
            % Upper and lower boundaries (i=1, i=ynum) 
            if((i==1 || i==ynum) && j<xnum)
                % Free slip, No slip: vy(i,j)=0
                L(invy,invy)        =   1*kbond;                                % Coefficient for vy(i,j)
                R(invy,1)           =   0;
            end
            % Left boundary, iner points (j=1, 1<i<ynum)
            if(j==1 && i>1 && i<ynum)
                % Free slip dvy/dx=0: vy(i,j)-vy(i,j+1)=0
                L(invy,invy)        =   1*kbond;                                % Coefficient for vy(i,j)
                L(invy,invy+ynum*3) =   -1*kbond;                               % Coefficient for vy(i,j+1)
                R(invy,1)           =   0;
    %             % No slip vy=0: vy(i,j)-1/3*vy(i,j+1)=0
    %             L(invy,invy)=1*kbond; % Coefficient for vy(i,j)
    %             L(invy,invy+ynum*3)=-1/3*kbond; % Coefficient for vy(i,j+1)
    %             R(invy,1)=0;
            end
            % Right boundary, iner points (j=xnum-1, 1<i<ynum)
            if(j==xnum-1 && i>1 && i<ynum)
                % Free slip dvy/dx=0: vy(i,j)-vy(i,j-1)=0
                L(invy,invy)        =   1*kbond;                                % Coefficient for vy(i,j)
                L(invy,invy-ynum*3) =  -1*kbond;                                % Coefficient for vy(i,j-1)
                R(invy,1)           =   0;
    %             % No slip vy=0: vy(i,j)-1/3*vy(i,j-1)=0
    %             L(invy,invy)=1*kbond; % Coefficient for vy(i,j)
    %             L(invy,invy-ynum*3)=-1/3*kbond; % Coefficient for vy(i,j-1)
    %             R(invy,1)=0;
            end
        %Internal nodes: dSyy/dy+dSxy/dx-dP/dy=-gy*RHO
        else
            %dSyy/dy=2*etan(i+1,j+1)*(vy(i+1,j)-vy(i,j))/dy^2-2*etan(i,j+1)*(vy(i,j)-vy(i-1,j))/dy^2
            L(invy,invy+3)          =   2*etan(i+1,j+1)/ystp^2;                 % Coefficient for vy(i+1,j)
            L(invy,invy-3)          =   2*etan(i,j+1)/ystp^2;                   % Coefficient for vy(i-1,j)
            L(invy,invy)            =   -2*etan(i+1,j+1)/ystp^2-2*etan(i,j+1)/ystp^2; % Coefficient for vy(i,j)

            %dSxy/dx=etas(i,j+1)*((vy(i,j+1)-vy(i,j))/dx^2+(vx(i,j+1)-vx(i-1,j+1))/dx/dy)-
            %         -etas(i,j)*((vy(i,j)-vy(i,j-1))/dx^2+(vx(i,j)-vx(i-1,j))/dx/dy)-
            L(invy,invy+ynum*3)     =   etas(i,j+1)/xstp^2;                     % Coefficient for vy(i,j+1)
            L(invy,invy-ynum*3)     =   etas(i,j)/xstp^2;                       % Coefficient for vy(i,j-1)
            L(invy,invy)            =   L(invy,invy)-etas(i,j+1)/xstp^2-etas(i,j)/xstp^2; % ADD coefficient for vy(i,j)
            L(invy,invx+ynum*3)     =   etas(i,j+1)/xstp/ystp;                  % Coefficient for vx(i,j+1)
            L(invy,invx+ynum*3-3)   =   -etas(i,j+1)/xstp/ystp;                 % Coefficient for vx(i-1,j+1)
            L(invy,invx)            =   -etas(i,j)/xstp/ystp;                   % Coefficient for vx(i,j)
            L(invy,invx-3)          =   etas(i,j)/xstp/ystp;                    % Coefficient for vx(i-1,j)

            % -dP/dy=(P(i,j+1)-P(i+1,j+1))/dx
            L(invy,inp+ynum*3)      =   kcont/ystp;                             % Coefficient for P(i,j+1)
            L(invy,inp+3+ynum*3)    =   -kcont/ystp;                            % Coefficient for P(i+1,j+1)
            % Right part: -RHO*gy
            R(invy,1)               =   -gy*(rho(i,j)+rho(i,j+1))/2;
        end

      end
    end
    
    %Obtaining vector of solutions S()
    S=L\R;
    
    % Reload solutions to 2D p(), vx(), vy() arrays
    % Dimensions of arrays are reduced compared to the basic grid
    p=zeros(ynum,xnum);
    vy=zeros(ynum,xnum);
    vx=zeros(ynum,xnum);
    % Process all Grid points
    for i=1:1:ynum
      for j=1:1:xnum
        % Global index for P, vx, vy in S()
        inp=((j-1)*ynum+i)*3-2; % P
        invx=inp+1;
        invy=inp+2;
        % P
        p(i,j)=S(inp)*kcont;
        % vx
        vx(i,j)=S(invx);
        % vy
        vy(i,j)=S(invy);
      end
    end

    % Compute vx,vy for internal nodes
    vx1=zeros(ynum,xnum);
    vy1=zeros(ynum,xnum);
    
    exx=zeros(ynum,xnum);
    eyy=zeros(ynum,xnum);
    exy=zeros(ynum,xnum);
    
    sigmaxx=zeros(ynum,xnum);
    sigmayy=zeros(ynum,xnum);
    sigmaxy=zeros(ynum,xnum);
    
    % Process internal Grid points
    for i=2:1:ynum-1
      for j=2:1:xnum-1
          % Velocity components
          vx1(i,j)=(vx(i-1,j)+vx(i,j))/2;
          vy1(i,j)=(vy(i,j-1)+vy(i,j))/2;
        
          exy(i,j)=(1/2)*((vx(i,j)-vx(i-1,j))/ystp+(vy(i,j)-vy(i,j-1))/xstp);
          sigmaxy(i,j)=2*etas(i,j)*exy(i,j);
      end
    end

    for i=2:1:ynum
        for j=2:1:xnum
            exx(i,j)=(vx(i-1,j)-vx(i-1,j-1))/xstp;
            eyy(i,j)=(-exx(i,j));
            sigmaxx(i,j)=2*etan(i,j)*exx(i,j);
            sigmayy(i,j)=(-sigmaxx(i,j));
        end
    end
        
    lith_vert_norm_stress = zeros(stepmax, xnum); % vertical normal stress at the base of the lithosphere
    lith_shear_stresses = zeros(stepmax, round((100000+waterdepth)/(ysize/ynum)+1), xnum);
    for count=1:1:xnum
        lith_vert_norm_stress(ntimestep, count) = sigmayy(round((100000+waterdepth)/(ysize/ynum)+1), count);
        lith_shear_stresses(ntimestep, :, :) = sigmaxy(1:round((100000+waterdepth)/(ysize/ynum)+1), :);
    end
    
    % write vertical normal stress at base of lithosphere to disk
    %dlmwrite('vertNormStress_1Myr_free.txt', lith_vert_norm_stress(ntimestep, :), '-append');
    % write shear stresses in lithosphere to disk
    for depth=1:1:round((waterdepth+100000)/(ysize/ynum)+1)
        %dlmwrite(sprintf('lithShearStress_1Myr_free_%d.txt', depth), lith_shear_stresses(ntimestep, depth, :), '-append');
    end
    % write normal stresses at the base of the lithosphere to disk
    if ntimestep == 50
        colRange = linspace(1, 50, 50);
        xlswrite(sprintf('lithVertNormStress.xls'), transpose(colRange(1:49)), 'A1:A50');
        xlswrite(sprintf('lithVertNormStress.xls'), transpose(-lith_vert_norm_stress(50, 2:50)), 'B1:B50');
        dlmwrite('vertNormStress.txt', -lith_vert_norm_stress(50, :), 'precision', 10);
    end

    log_exx=zeros(ynum,xnum);
    log_neg_exx=zeros(ynum,xnum);
    log_exy=zeros(ynum,xnum);
    log_neg_exy=zeros(ynum,xnum);
    
    log_sigmaxx=zeros(ynum,xnum);
    log_neg_sigmaxx=zeros(ynum,xnum);
    log_sigmaxy=zeros(ynum,xnum);
    log_neg_sigmaxy=zeros(ynum,xnum);
    
    for i=2:1:ynum-1
        for j=2:1:xnum-1
            % Log strain rates
            if(exx(i,j)==0)
                log_exx(i,j)=NaN;
                log_neg_exx(i,j)=NaN;
            elseif(exx(i,j)<0)
                log_neg_exx(i,j)=(sign(exx(i,j))*(log10(abs(exx(i,j)))));
                log_exx(i,j)=NaN;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %log_exx(i,j)=(sign(exx(i,j))*(log10(abs(exx(i,j)))));
            else
                log_exx(i,j)=log10(exx(i,j));
                log_neg_exx(i,j)=NaN;
            end
            if(exy(i,j)==0)
                log_exy(i,j)=NaN;
                log_neg_exy(i,j)=NaN;
            elseif(exy(i,j)<0)
                log_neg_exy(i,j)=(sign(exy(i,j))*(log10(abs(exy(i,j)))));
                log_exy(i,j)=NaN;
            else
                log_exy(i,j)=log10(exy(i,j));
                log_neg_exy(i,j)=NaN;
            end
            % Log stresses
            if(sigmaxx(i,j)==0)
                log_sigmaxx(i,j)=NaN;
                log_neg_sigmaxy(i,j)=NaN;
            elseif(sigmaxx(i,j)<0)
                log_neg_sigmaxx(i,j)=(sign(sigmaxx(i,j))*(log10(abs(sigmaxx(i,j)))));
                log_sigmaxx(i,j)=NaN;
            else
                log_sigmaxx(i,j)=log10(sigmaxx(i,j));
                log_neg_sigmaxx(i,j)=NaN;
            end
            if(sigmaxy(i,j)==0)
                log_sigmaxy(i,j)=NaN;
                log_neg_sigmaxy(i,j)=NaN;
            elseif(sigmaxy(i,j)<0)
                log_neg_sigmaxy(i,j)=(sign(sigmaxy(i,j))*(log10(abs(sigmaxy(i,j)))));
                log_sigmaxy(i,j)=NaN;
            else
                log_sigmaxy(i,j)=log10(sigmaxy(i,j));
                log_neg_sigmaxy(i,j)=NaN;
            end
        end
    end
    
    low_bound_epsilon = min([(min(min(exx))) (min(min(exy))) (min(min(eyy)))]);
    up_bound_epsilon = max([(max(max(exx))) (max(max(exy))) (max(max(eyy)))]);
    low_bound_sigma = min([(min(min(sigmaxx))) (min(min(sigmayy))) (min(min(sigmaxy)))]);
    up_bound_sigma = max([(max(max(sigmaxx))) (max(max(sigmayy))) (max(max(sigmaxy)))]);
    
    low_bound_lith_sigmaxy = -1e+6;
    up_bound_lith_sigmaxy = 1e+6;
    
    % set location of free surface before advection of markers
    surfCoords(ntimestep, :, 1) = MX(surfInd(:));
    surfCoords(ntimestep, :, 2) = MY(surfInd(:));
    if mod(ntimestep, 50) == 0
        figure(counter), clf;
        scatter(surfCoords(ntimestep, :, 1), surfCoords(ntimestep, :, 2))
        ylim([waterdepth-5000, waterdepth+5000])
        axis ij;
        counter = counter + 1;
    elseif ntimestep == 1
        figure(counter), clf;
        scatter(surfCoords(ntimestep, :, 1), surfCoords(ntimestep, :, 2))
        ylim([waterdepth-5000, waterdepth+5000])
        axis ij;
        counter = counter + 1;
    end
%{
%% FIGURE 1 - DENSITY AND VISCOSITY w/ VELOCITY FIELD
    figure(counter), clf;
    % Plotting log viscosity as colormap
    subplot(1,2,1)
    pcolor(x/1000,y/1000,log10(etas));
    shading interp;     % making smooth transitions between colors
    colorbar;           % showing a colorbar for the map
    hold on;            % continuing plotting on the colormap
    % Plotting velocity vector as arrows using internal nodes only
    quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
    hold off;           % stop plotting on the colormap
    box on;             % making a box around the plot
    title('log_{10} viscosity (color, Pa s), velocity (arrows)'); % title for the plot
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits

    % Plotting density as colormap
    subplot(1,2,2);
    pcolor(x/1000,y/1000,rho);      % making a colormap
    shading interp;     % making smooth transitions between colors
    colorbar;           % showing a colorbar for the map
    hold on;            % continuing plotting on the colormap
    % Plotting velocity vector as arrows using internal nodes only
    quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
    hold off;           % stop plotting on the colormap
    box on;             % making a box around the plot
    title(['Density (kg/m^3) Step=',num2str(ntimestep),' Myr=',num2str(timesum*1e-6/(365.25*24*3600))]); 
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    drawnow % draw above figure without delay
    % saveas(figure(counter), sprintf('visc_den_vel_%d', ntimestep), 'png');
    counter = counter + 1;

% FIGURE 2 - VELOCITY COMPONENTS   
    % Plotting velocity componenets
    figure(counter), clf;
    subplot(1,2,1);
    pcolor(x/1000,y/1000,vx1);
    shading interp;
    colorbar;
    hold on;
    quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
    hold off;           % stop plotting on the colormap
    title('v_x_1');
    box on;
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    
    subplot(1,2,2);
    pcolor(x/1000,y/1000,vy1);
    shading interp;
    colorbar;
    hold on;
    quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
    hold off;           % stop plotting on the colormap
    title('v_x_1');
    box on;
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    % saveas(figure(counter), sprintf('vel_comps_%d', ntimestep), 'png');
    counter = counter + 1;

% FIGURE 3 - STRAIN RATES    
    % Plotting strain rates and stresses
    figure(counter), clf;
    
    subplot(1,2,1);
    pcolor(x/1000,y/1000,exx);
    shading interp;
    colorbar;
    box on;
    title('\epsilon_x_x strain rate');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([low_bound_epsilon up_bound_epsilon]);

    subplot(1,2,2);
    pcolor(x/1000,y/1000,exy);
    shading interp;
    colorbar;
    box on;
    title('\epsilon_x_y strain rate');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([low_bound_epsilon up_bound_epsilon]);
    saveas(figure(counter), sprintf('strainrate_%d', ntimestep), 'png');
    % saveas(figure(counter), fullfile('/Users/jaredbryan/Desktop/MatlabFigs/strainrate/' ,sprintf('strainrate_%d', ntimestep)), 'png');
    counter = counter + 1;


% FIGURE 4 - STRESSES
    figure(counter), clf;
    % Stresses
    subplot(3,1,1);
    pcolor(x/1000,y/1000,sigmaxx);
    shading interp;
    colorbar;
    box on;
    title('\sigma''_x_x');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([low_bound_sigma up_bound_sigma]);

    subplot(3,1,2);
    pcolor(x/1000,y/1000,sigmaxy);
    shading interp;
    colorbar;
    box on;
    title('\sigma''_x_y');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([low_bound_sigma up_bound_sigma]);
    % saveas(figure(counter), sprintf('stress_%d', ntimestep), 'png');
    counter = counter + 1;
    % Scale to cm/year with *100*60*60*24*365.25

    subplot(3,1,3);
    pcolor(x/1000,y/1000,sigmayy);
    shading interp;
    colorbar;
    box on;
    title('\sigma''_y_y');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([low_bound_sigma up_bound_sigma]);
%}
%{   
    figure(counter), clf;
    subplot(3,1,1)
    pcolor(x/1000,y(1:round((100000+waterdepth)/(ysize/ynum)+1))/1000,sigmaxx(1:round((100000+waterdepth)/(ysize/ynum)+1), :));
    shading interp;
    colorbar;
    box on;
    title('\sigma''_x_x');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 (100000+waterdepth)/1000]); % Making axes limits
    caxis([low_bound_sigma up_bound_sigma]);
    
    subplot(3,1,2)
    pcolor(x/1000,y(1:round((100000+waterdepth)/(ysize/ynum)+1))/1000,sigmaxy(1:round((100000+waterdepth)/(ysize/ynum)+1), :));
    shading interp;
    colorbar;
    box on;
    title('\sigma''_x_y');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 (100000+waterdepth)/1000]); % Making axes limits
    caxis([low_bound_sigma up_bound_sigma]);
    
        subplot(3,1,2)
    pcolor(x/1000,y(1:round((100000+waterdepth)/(ysize/ynum)+1))/1000,sigmayy(1:round((100000+waterdepth)/(ysize/ynum)+1), :));
    shading interp;
    colorbar;
    box on;
    title('\sigma''_y_y');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 (100000+waterdepth)/1000]); % Making axes limits
    caxis([low_bound_sigma up_bound_sigma]);

    log_low_bound_sigma = sign(low_bound_sigma) * log10(abs(low_bound_sigma));
    log_low_bound_epsilon = sign(low_bound_epsilon) * log10(abs(low_bound_epsilon));
    log_up_bound_sigma = log10(up_bound_sigma);
    log_up_bound_epsilon = log10(up_bound_epsilon);
%}
%{  
% FIGURE 5 - LOG STRAIN RATES    
    % Plotting log strain rates and stresses
    figure(counter), clf;
    subplot(2,2,1);
    pcolor(x/1000,y/1000,log_exx);
    shading interp;
    colorbar;
    box on;
    title('log \epsilon_x_x strain rate');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    %caxis([log_low_bound_epsilon log_up_bound_epsilon]);
    
    subplot(2,2,2);
    pcolor(x/1000,y/1000,log_exy);
    shading interp;
    colorbar;
    box on;
    title('log \epsilon_x_y strain rate');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    %caxis([log_low_bound_epsilon log_up_bound_epsilon]);
    
    subplot(2,2,3);
    pcolor(x/1000,y/1000,log_neg_exx);
    shading interp;
    colorbar;
    box on;
    title('log \epsilon_x_x strain rate');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    %caxis([log_low_bound_epsilon log_up_bound_epsilon]);
    
    subplot(2,2,4);
    pcolor(x/1000,y/1000,log_neg_exy);
    shading interp;
    colorbar;
    box on;
    title('log \epsilon_x_y strain rate');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    %caxis([log_low_bound_epsilon log_up_bound_epsilon]);
    
    % saveas(figure(counter), sprintf('logstrainrate%d', ntimestep), 'png');
    counter = counter + 1;
   
% FIGURE 6 - LOG STRESSES
    % Log Stresses
    figure(counter), clf;
    subplot(2,2,1);
    pcolor(x/1000,y/1000,log_sigmaxx);
    shading interp;
    colorbar;
    box on;
    title('log \sigma''_x_x');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([0 log_up_bound_sigma]);

    subplot(2,2,2);
    pcolor(x/1000,y/1000,log_sigmaxy);
    shading interp;
    colorbar;
    box on;
    title('log \sigma''_x_y');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([0 log_up_bound_sigma]);
    
    subplot(2,2,3);
    pcolor(x/1000,y/1000,log_neg_sigmaxx);
    shading interp;
    colorbar;
    box on;
    title('log neg \sigma''_x_x');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([log_low_bound_sigma 0]);

    subplot(2,2,4);
    pcolor(x/1000,y/1000,log_neg_sigmaxy);
    shading interp;
    colorbar;
    box on;
    title('log neg \sigma''_x_y');
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
    caxis([log_low_bound_sigma 0]);

    % saveas(figure(counter), sprintf('logstress%d', ntimestep), 'png');
    counter = counter + 1;
 
    figure(counter); clf;
    subplot(3,1,1);
    pcolor(x/1000,y(1:12)/1000,sigmaxx(1:12, :));
    shading interp;
    colorbar;
    box on;
    title(['\sigma''_x_x step ', num2str(ntimestep)]);
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 200]); % Making axes limits
    %caxis([low_bound_lith_sigmaxy up_bound_lith_sigmaxy])
    subplot(3,1,2);
    pcolor(x/1000,y(1:12)/1000,sigmaxy(1:12, :));
    shading interp;
    colorbar;
    box on;
    title(['\sigma''_x_y step ', num2str(ntimestep)]);
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 200]); % Making axes limits
    subplot(3,1,3);
    pcolor(x/1000,y(1:12)/1000,sigmayy(1:12, :));
    shading interp;
    colorbar;
    box on;
    title(['\sigma''_y_y step ', num2str(ntimestep)]);
    xlabel('x, km');        % title for the horizontal axis
    ylabel('y, km');        % title for the vertical axis
    axis ij image ;     % directing vertical axis downward, making proper dimensions
    axis([0 xsize/1000 0 200]); % Making axes limits
    %caxis([low_bound_lith_sigmaxy up_bound_lith_sigmaxy])
    %caxis([low_bound_lith_sigmaxy up_bound_lith_sigmaxy])
    %saveas(figure(counter), sprintf('shearStressesLith%d', ntimestep), 'png');
    counter = counter + 1;
%}    
    
    % Check maximal velocity
    vxmax=max(max(abs(vx)));
    vymax=max(max(abs(vy)));
    % Set initial value for the timestep
    timestep=timemax;
    % Check marker displacement step
    if (vxmax>0)
        if (timestep>markmax*xstp/vxmax)
            timestep=markmax*xstp/vxmax;
        end
    end
    if (vymax>0)
        if (timestep>markmax*ystp/vymax)
            timestep=markmax*ystp/vymax;
        end
    end

    % Moving Markers by velocity field
    if(markmove>0)
        % Create arrays for velocity of markers
        vxm=zeros(4,1);
        vym=zeros(4,1);
        % Marker cycle
        for mm1 = 1:1:marknum

            % Check markers inside the grid
            if (MX(mm1)>=gridx(1) && MX(mm1)<=gridx(xnum) && MY(mm1)>=gridy(1) && MY(mm1)<=gridy(ynum)) 

                % Save marker coordinates
                xcur=MX(mm1);
                ycur=MY(mm1);
                % Defining number of Runge-Kutta cycles
                for rk=1:1:markmove

                    %  xn    V(xn,yn)--------------------V(xn+1,yn)
                    %           ?           ^                  ?
                    %           ?           ?                  ?
                    %           ?          dy                  ?
                    %           ?           ?                  ?
                    %           ?           v                  ?
                    %           ?<----dx--->o Mrho(xm,ym)       ?
                    %           ?                              ?
                    %           ?                              ?
                    %  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)

                    % Define indexes for upper left node in the Vx-cell where the marker is
                    % Vx-cells are displaced downward for 1/2 of vertical gridstep
                    % !!! SUBTRACT 0.5 since int16(0.5)=1
                    xn=double(int16(xcur./xstp-0.5))+1;
                    yn=double(int16((ycur-ystp/2.0)./ystp-0.5))+1;
                    % Check indexes:
                    % vertical index for upper left Vx-node must be between 1 and ynum-2
                    if (xn<1)
                        xn=1;
                    end
                    if (xn>xnum-1)
                        xn=xnum-1;
                    end
                    if (yn<1)
                        yn=1;
                    end
                    if (yn>ynum-2)
                        yn=ynum-2;
                    end

                    % Define and check normalized distances from marker to the upper left Vx-node;
                    dx=(xcur-gridx(xn))./xstp;
                    dy=(ycur-gridcy(yn))./ystp;

                    % Calculate Marker velocity from four surrounding Vx nodes
                    vxm(rk)=0;
                    vxm(rk)=vxm(rk)+(1.0-dx).*(1.0-dy).*vx(yn,xn);
                    vxm(rk)=vxm(rk)+(1.0-dx).*dy.*vx(yn+1,xn);
                    vxm(rk)=vxm(rk)+dx.*(1.0-dy).*vx(yn,xn+1);
                    vxm(rk)=vxm(rk)+dx.*dy.*vx(yn+1,xn+1);

                    % Define indexes for upper left node in the Vy-cell where the marker is
                    % Vy-cells are displaced rightward for 1/2 of horizontal gridstep
                    % !!! SUBTRACT 0.5 since int16(0.5)=1
                    xn=double(int16((xcur-xstp/2.0)./xstp-0.5))+1;
                    yn=double(int16(ycur./ystp-0.5))+1;
                    % Check indexes:
                    % horizontal index for upper left Vy-node must be between 1 and xnum-2
                    if (xn<1)
                        xn=1;
                    end
                    if (xn>xnum-2)
                        xn=xnum-2;
                    end
                    if (yn<1)
                        yn=1;
                    end
                    if (yn>ynum-1)
                        yn=ynum-1;
                    end

                    % Define and check normalized distances from marker to the upper left Vy-node;
                    dx=(xcur-gridcx(xn))./xstp;
                    dy=(ycur-gridy(yn))./ystp;

                    % Calculate Marker velocity from four surrounding nodes
                    vym(rk)=0;
                    vym(rk)=vym(rk)+(1.0-dx).*(1.0-dy).*vy(yn,xn);
                    vym(rk)=vym(rk)+(1.0-dx).*dy.*vy(yn+1,xn);
                    vym(rk)=vym(rk)+dx.*(1.0-dy).*vy(yn,xn+1);
                    vym(rk)=vym(rk)+dx.*dy.*vy(yn+1,xn+1);

                    % Update coordinates for the next Runge-Kutta cycle
                    if(rk<4)
                        if (rk<3)
                            xcur=MX(mm1)+timestep/2*vxm(rk);
                            ycur=MY(mm1)+timestep/2*vym(rk);
                        else
                            xcur=MX(mm1)+timestep*vxm(rk);
                            ycur=MY(mm1)+timestep*vym(rk);
                        end
                    end
                end
                % Recompute velocity using 4-th order Runge_Kutta
                if (markmove==4)
                    vxm(1)=(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4))/6;
                    vym(1)=(vym(1)+2*vym(2)+2*vym(3)+vym(4))/6;
                end

                % Displacing Marker according to its velocity
                MX(mm1)=MX(mm1)+timestep*vxm(1);
                MY(mm1)=MY(mm1)+timestep*vym(1);
                
            end
        end
    end
    
    % Advance in time
    timesum=timesum+timestep;
    toc()
    
end
