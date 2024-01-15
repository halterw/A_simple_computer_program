% Iterative finite difference program for 2D viscous deformation, Matlab version
% This code is free software under the Creative Commons CC-BY-NC-ND license
% Published in W.R.Halter, E.Macherel, and S.M.Schmalholz (2022) JSG, https://doi.org/10.1016/j.jsg.2022.104617
clear all, close all, clc, tic, load('colormaps\vik.mat'), colormap(vik) % Colormap from Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
ps          = 1;        % ps = 1 models pure shear; ps = 0 models simple shear
etaRatio    = 1e-3;     % Viscosity ratio (matrix/inclusion)
nout        = 1e3;
% Numerical parameters
nx          = 301;              ny          = 186;                  % Numerical resolution
Lx          = 1;                Ly          = Lx/nx*ny;             % Model dimension
dx          = Lx/(nx-1);        dy          = Ly/(ny-1);            % Grid spacing
x           = [-Lx/2:dx:Lx/2];  y           = [-Ly/2:dy:Ly/2];      % Coordinate vectors
x_vx        = [x(1)-dx/2,( x(1:end-1) + x(2:end) )/2,x(end)+dx/2];  % Horizontal vector for Vx which is one more than basic grid
y_vy        = [y(1)-dy/2,( y(1:end-1) + y(2:end) )/2,y(end)+dy/2];  % Vertical   vector for Vy which is one more than basic grid
% 3 numerical grids due to staggered grid
[X,Y]       = ndgrid(x,y);      [X_vx,Y_vx] = ndgrid(x_vx,y); [X_vy,Y_vy] = ndgrid(x,y_vy);
% Physical parameters
eta_B       = 1;                n_exp       = 1;
s_ref       = 1.5;              D_B         =-1;
% Initialization
P         	= zeros(nx,  ny);   ETA       	= eta_B*ones(nx, ny);
% Define garnet inclusion
G = [
    -0.0924021300036714, -0.0665524544007828;
     0.0198515118129514, -0.1222508262945290;
     0.1346758477169800, -0.0125678785653086;
     0.0172808177255486,  0.1176806218631410;
    -0.0906883339454028,  0.0508425755906477;
];
inside = inpolygon(X,Y, G(:,1), G(:,2));
ETA(inside) = eta_B/etaRatio;   etamin      = 0.1*min(ETA(:));
for smo=1:2; Iix  = [2:nx-1]; Iiy  = [2:ny-1];                      % Smoothing of the initial viscosity field
    ETA(Iix,:)    = ETA(Iix,:) + 0.4*(ETA(Iix+1,:)-2*ETA(Iix,:)+ETA(Iix-1,:));
    ETA(:,Iiy)    = ETA(:,Iiy) + 0.4*(ETA(:,Iiy+1)-2*ETA(:,Iiy)+ETA(:,Iiy-1));
end
ETA_L           = ETA;              ETA_XY          = n2c(ETA);     ETA_PL       = ETA;
RES_VX_relaxed  = zeros(nx-1,ny-2); RES_VY_relaxed  = zeros(nx-2,ny-1);
% Boundary condition
if     ps==1,   VX      = -D_B*X_vx;        VY      = D_B*Y_vy;
elseif ps==0,   VX      =  D_B*Y_vx;        VY      =   0*Y_vy; end
% Parameters for pseudo-transient iterations
tol             = 1e-6;             err_absolute    = 1;            err_relative = 1;
%%% Numerical parameters for iterative solver %%%
CFLV            = 1/5;              CFLP            = CFLV/nx/20;   % Numerical parameter for iterative solver
dpt_P           = zeros(nx,ny);
dpt_Vx          = CFLV./(max(ETA(1:end-1,2:end-1),ETA(2:end,2:end-1))/dx^2 + max(ETA_XY(:,1:end-1),ETA_XY(:,2:end))/dx/dy);
dpt_Vy          = CFLV./(max(ETA(2:end-1,1:end-1),ETA(2:end-1,2:end))/dy^2 + max(ETA_XY(1:end-1,:),ETA_XY(2:end,:))/dx/dy);
dpt_V_dx2       = max(max(dpt_Vx(1:end-1,:),dpt_Vx(2:end,:))/dx^2, max(dpt_Vy(:,1:end-1),dpt_Vy(:,2:end))/dy^2);
dpt_P(2:end-1,2:end-1) = CFLP./dpt_V_dx2;
dpt_P([1 end],:)       = dpt_P([2 end-1],:);  dpt_P(:,[1 end])  = dpt_P(:,[2 end-1]);
%%% Numerical parameters for iterative solver %%%
iter        = 0;
while err_absolute(end)>tol||err_relative(end)>tol; iter = iter+1;  % START of iteration loop
    DXX                 = diff(VX,1,1)/dx;                          % Eq (1)
    DYY                 = diff(VY,1,2)/dy;                          % Eq (2)
    DXY                 = 1/2*( diff(VX(2:end-1,:),1,2)/dy ...
                              + diff(VY(:,2:end-1),1,1)/dx );       % Eq (3)
    TXX                 = 2.*ETA.*DXX;                              % Eq (4)
    TYY                 = 2.*ETA.*DYY;                              % Eq (5)
    TXY                 = 2.*n2c(ETA).*DXY;                         % Eq (6)
    SXX                 = -P(:,2:end-1)+TXX(:,2:end-1);             % Eq (7)
    SYY                 = -P(2:end-1,:)+TYY(2:end-1,:);             % Eq (8)
    RES_P               = -( diff(VX,1,1)/dx + diff(VY,1,2)/dy );   % Eq (16)
    RES_VX              = diff(SXX,1,1)/dx + diff(TXY,1,2)/dy;      % Eq (17)
    RES_VY              = diff(SYY,1,2)/dy + diff(TXY,1,1)/dx;      % Eq (18)
    RES_VX_relaxed      = RES_VX_relaxed*(1-5/nx) + RES_VX.*dpt_Vx; % Relaxation on residual      
    RES_VY_relaxed      = RES_VY_relaxed*(1-5/ny) + RES_VY.*dpt_Vy; % Relaxation on residual      
    P                   = P                   + dpt_P.*RES_P;       % Eq (22)
    VX(2:end-1,2:end-1) = VX(2:end-1,2:end-1) + RES_VX_relaxed;     % Eq (23)
    VY(2:end-1,2:end-1) = VY(2:end-1,2:end-1) + RES_VY_relaxed;     % Eq (24)
    TII                 = sqrt(0.5*(TXX.^2 + TYY.^2 + 2*c2n(TXY).^2)); % Eq (13)
    if n_exp>1 % Power-law viscosity
        ETA_PL_it       = ETA_PL;           % Viscosity of previous iteration step
        ETA_PL          = ETA_L.*(TII/s_ref).^(1-n_exp);            % Eq (12)
        ETA_PL          = exp(log(ETA_PL)*0.5+log(ETA_PL_it)*0.5);
        ETA             = 1./( 1./ETA_L + 1./ETA_PL );              % Eq (14)
        ETA(inside)     = ETA_L(inside);    % Power-law viscosity only applied to matrix
        ETA(ETA<etamin) = etamin;           % Minimum viscosity
    end
    err_absolute(iter) = max([max(abs(RES_VX(:))), max(abs(RES_VY(:))), max(abs(RES_P(:)))]);
    err_relative(iter) = max([max(abs(RES_VX_relaxed)), max(abs(RES_VY_relaxed)), max(abs(dpt_P.*RES_P))]);
    if mod(iter,nout)==0 % Visualization during calculation
        subplot(221),pcolor(X,Y,(P-mean(P(:)))/(2*eta_B*abs(D_B))),shading interp,axis equal,axis tight,
        caxis([min((P(:)-mean(P(:)))/(2*eta_B*abs(D_B))) max((P(:)-mean(P(:)))/(2*eta_B*abs(D_B)))]/1.5), colorbar
        title({'A) Pressure, P / (2\eta_BD_B)'}),ylabel('Height [ ]')
        subplot(222),pcolor(X,Y,log10(ETA/eta_B)), shading interp, hold on
        title(['B) Log_{10} of Viscosity, \eta / \eta_B, and velocity arrows'])
        st = 15; VXC = (VX(1:end-1,:)+VX(2:end,:))/2; VYC = (VY(:,1:end-1)+VY(:,2:end))/2;
        quiver(X(1:st:end,1:st:end),Y(1:st:end,1:st:end),VXC(1:st:end,1:st:end),VYC(1:st:end,1:st:end),'w');
        axis equal,axis tight,colorbar,axis([-Lx/2,Lx/2,-Ly/2,Ly/2]),hold off
        subplot(223),pcolor(X,Y,TII/(2*eta_B*abs(D_B))),shading interp,axis equal,axis tight,colorbar
        title(['C) 2nd stress invariant, T_{II} / (2\eta_BD_B)']), xlabel('Width [ ]'), ylabel('Height [ ]')
        subplot(224),loglog(1:length(err_absolute),err_absolute,1:length(err_relative),err_relative),title('D) Error (residual) convergence'),xlabel('Iterations'),ylabel('Residual')
        axis([1 iter+1e4 tol/10 1e6]), hold on,plot([1 iter+1e4],[tol tol],'--k'),text(1e3,tol*2,'Tolerance'), legend('absolute error','relative error',''), hold off
        set(gcf,'position',[281 47.4000 975.2000 724]), drawnow
    end
end, toc % END of iteration loop
% Additional functions perfoming interpolations on the numerical grid
function A1 = n2c(A0) % Interpolation of nodal points to center points
A1      = (A0(2:end,:) + A0(1:end-1,:))/2;
A1      = (A1(:,2:end) + A1(:,1:end-1))/2;
end
function A2 = c2n(A0) % Interpolation of center points to nodal points
A1    	= zeros(size(A0,1)+1,size(A0,2));
A1(:,:)	= [1.5*A0(1,:)-0.5*A0(2,:); (A0(2:end,:)+A0(1:end-1,:))/2; 1.5*A0(end,:)-0.5*A0(end-1,:)];
A2   	= zeros(size(A1,1),size(A1,2)+1);
A2(:,:)	= [1.5*A1(:,1)-0.5*A1(:,2), (A1(:,2:end)+A1(:,1:end-1))/2, 1.5*A1(:,end)-0.5*A1(:,end-1)];
end