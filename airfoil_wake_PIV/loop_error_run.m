SourceFolder = '../../Experiment_wing3/OUT_TRPIV/';
SaveFolder = 'amic_loop';
ROI = [41 110 41 150];         % subdomain of PIV field
                               % [first-last row, first-last column]
Mm_per_px_ratio = 0.1197;      % spatial resolution of photos
Sample_rate = 30;              % sample rate of PIV
Vector_Spacing = 10;           % vector per pixel in PIV
Nu = 1e-6;                     % kinetic viscosity coefficient of water
Rho = 1e3;                     % density of water

AFrame_F = 11:430;             % frames of filtering
AFrame_P = 111:310;            % frames of pressure calculation
InitLoop = 0;                  % initial loop of AMIC
FinalLoop = 50;                % final loop of AMIC
Lambda = [0.1 0.05];           % iteration coefficients

load([SourceFolder, 'Grid_Wing.mat'], 'X','Y');
% entire domain coordinates
X_etr = X;                     Y_etr = Y;
% map of subdomain
map = false(size(X_etr));
map(ROI(1):ROI(2), ROI(3):ROI(4)) = true;
% subdomain coordinates normalized by the chord length
X = X_etr(ROI(1):ROI(2),ROI(3):ROI(4));
Y = Y_etr(ROI(1):ROI(2),ROI(3):ROI(4));
X = (X - X(1)).*Mm_per_px_ratio./80 + 0.70;
Y = (Y - Y(1)).*Mm_per_px_ratio./80 - 0.53;
dx = Vector_Spacing*Mm_per_px_ratio*1e-3;
dy = dx;
dt = 1/Sample_rate;
mkdir(SaveFolder);
tic

disp('loading data...');
U_array = zeros(numel(X_etr), length(AFrame_F));   V_array = U_array;
iCount = 0;
for iFrame = AFrame_F
    iCount = iCount + 1;
    tmp_file = [SourceFolder,'Wing_',num2str(iFrame,'%06u'),'.mat'];
    load(tmp_file, 'U','V');
    U_array(:, iCount) = V(:);    V_array(:, iCount) = U(:);
    % swapping U and V
end
U_array = U_array.*Sample_rate.*Mm_per_px_ratio.*1e-3;
V_array = V_array.*Sample_rate.*Mm_per_px_ratio.*1e-3;
U_array0 = U_array;            V_array0 = V_array;
clear iFrame iCount U V
if InitLoop == 0
else
    load([SaveFolder,'/',num2str(InitLoop, '%03u'),'.mat']);
end
Enabled = ones(size(U_array));
toc


% processing
Vrms = rms(U_array,'all') + rms(V_array,'all');
P_AMIC = zeros(numel(X), length(AFrame_P));
for iLoop = InitLoop:FinalLoop
    % turbulence frozen filter
    % sequence frame P+1 P-1 P+2 P-2 P+3 P-3...
    U_cor = zeros([size(U_array) length(Lambda)*2]);
    V_cor = U_cor;
    for iFrame = 1:length(AFrame_F)
        u0    = reshape(U_array(:,iFrame), size(X_etr));
        v0    = reshape(V_array(:,iFrame), size(X_etr));
        mask0 = reshape(Enabled(:,iFrame), size(X_etr));
        u = u0;   v = v0;   mask = mask0;
        for iStep = 1:length(Lambda)
            [u, v] = MRK4Prop2D(u, v, X_etr, Y_etr, dt);
            mask = interp2(X_etr, Y_etr, mask, X_etr-u*dt, Y_etr-v*dt);
            mask(isnan(mask)) = 0; mask = imgaussfilt(mask, 5);
            iSub = iFrame + iStep;
            if iSub < length(AFrame_F)
                U_cor(:,iSub,iStep*2-1) = (u(:)-U_array(:,iSub)).*mask(:);
                V_cor(:,iSub,iStep*2-1) = (v(:)-V_array(:,iSub)).*mask(:);
            end
        end
        u = u0;   v = v0;   mask = mask0;
        for iStep = 1:length(Lambda)
            [u, v] = MRK4Prop2D(u, v, X_etr, Y_etr,-dt);
            mask = interp2(X_etr, Y_etr, mask, X_etr+u*dt, Y_etr+v*dt);
            mask(isnan(mask)) = 0; mask = imgaussfilt(mask, 5);
            iSub = iFrame - iStep;
            if iSub > 0
                U_cor(:,iSub,iStep*2) = (u(:)-U_array(:,iSub)).*mask(:);
                V_cor(:,iSub,iStep*2) = (v(:)-V_array(:,iSub)).*mask(:);
            end
        end
    end
    if iLoop == 0
        fprintf(['iLoop = ',num2str(iLoop),'\n']);
        rms_PIV = 0;
    else
        % ratio of update
        idx_update = (rms(U_cor(U_cor~=0))+rms(V_cor(V_cor~=0)))/Vrms;
        % correcting
        for iStep = 1:length(Lambda)
            U_array = U_array...
                + Lambda(iStep)*(U_cor(:,:,iStep*2-1)+U_cor(:,:,iStep*2));
            V_array = V_array...
                + Lambda(iStep)*(V_cor(:,:,iStep*2-1)+V_cor(:,:,iStep*2));
        end
        % RMS to the PIV field
        rms_PIV =...
            sqrt(mean(([U_array;V_array]-[U_array0;V_array0]).^2,'all'));
        fprintf(['iLoop = ',num2str(iLoop),...
            '\tvelocity correction = ',num2str(idx_update,'%.5g'),'\n']);
        fprintf(['RMS to PIV field = ',num2str(rms_PIV,'%.5g'),'\n']);
    end
    % pressure calculation
    parfor iFrame = 1:length(AFrame_P)
        p_iFrame = iFrame - AFrame_F(1) + AFrame_P(1);
        u1 = reshape(U_array(map(:),p_iFrame-1), size(X));
        v1 = reshape(V_array(map(:),p_iFrame-1), size(X));
        u2 = reshape(U_array(map(:),p_iFrame  ), size(X));
        v2 = reshape(V_array(map(:),p_iFrame  ), size(X));
        u3 = reshape(U_array(map(:),p_iFrame+1), size(X));
        v3 = reshape(V_array(map(:),p_iFrame+1), size(X));
        ut = (u3 - u1)./(2*dt);
        vt = (v3 - v1)./(2*dt);
        P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,dx,dy,0*X);
        P_AMIC(:,iFrame) = P(:);
    end
    [xcor0, ycor0] = meshgrid(1:size(X,2),1:size(X,1));
    for iFrame = 2:length(AFrame_P)
        % the position of grid point at frame 0.5 moving to frame 0
        p_iFrame = iFrame - AFrame_F(1) + AFrame_P(1) - 1;
        u = reshape(U_array(map(:), p_iFrame), size(X));
        v = reshape(V_array(map(:), p_iFrame), size(X));
        u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
        v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
        xcor1 = xcor0 - u./2./Vector_Spacing;
        ycor1 = ycor0 - v./2./Vector_Spacing;
        % the position of grid point at frame 0.5 moving to frame 1
        p_iFrame = iFrame - AFrame_F(1) + AFrame_P(1);
        u = reshape(U_array(map(:), p_iFrame), size(X));
        v = reshape(V_array(map(:), p_iFrame), size(X));
        u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
        v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
        xcor2 = xcor0 + u./2./Vector_Spacing;
        ycor2 = ycor0 + v./2./Vector_Spacing;
        % minimize dp/dt over the available grid points
        p1 = interp2(reshape(P_AMIC(:,iFrame-1),size(xcor1)), xcor1,ycor1);
        p2 = interp2(reshape(P_AMIC(:,iFrame),  size(xcor2)), xcor2,ycor2);
        P_AMIC(:,iFrame) =...
            P_AMIC(:,iFrame)-mean(p2(:)-p1(:),'all','omitnan');
    end
    res_P = ResDP2D(U_array(map(:),AFrame_P), V_array(map(:),AFrame_P),...
        P_AMIC, Rho, Nu, X, dx, dt);
    fprintf(['pressure dynamic error = ',...
        num2str(mean(res_P),'%.5g'),'\n']);
    save([SaveFolder,'/',num2str(iLoop, '%03u'),'.mat'],...
        'AFrame_F', 'AFrame_P', 'U_array', 'V_array', 'P_AMIC',...
        'res_P', 'rms_PIV');
    toc
end

%% functions
function [ut, vt] = TaylorDt2D(u, v, Increx, Increy)
% ut = -(u_mean\cdot\nabla)u_fluc
u_mean = imgaussfilt(u, 7);
v_mean = imgaussfilt(v, 7);
u_fluc = u - u_mean;
v_fluc = v - v_mean;
[ux, uy] = gradient(u_fluc); ux = ux./Increx; uy = uy./Increy;
[vx, vy] = gradient(v_fluc); vx = vx./Increx; vy = vy./Increy;
ut = -u_mean.*ux - v_mean.*uy;
vt = -u_mean.*vx - v_mean.*vy;
ut = imgaussfilt3(ut, 0.25);
vt = imgaussfilt3(vt, 0.25);
end

function [u, v] = MRK4Prop2D(u0, v0, x, y, dt)
% Runge-Kutta propogation of 2D flow fields
% u0, v0, current velocity
% x, y, corrodinates
% dt time increments
% u, v, velocity after dt
dx = x(1,2) - x(1);            dy = y(2,1) - y(1);
[ut1, vt1] = TaylorDt2D(u0,            v0,            dx, dy);
[ut2, vt2] = TaylorDt2D(u0+0.5*ut1*dt, v0+0.5*vt1*dt, dx, dy);
[ut3, vt3] = TaylorDt2D(u0+0.5*ut2*dt, v0+0.5*vt2*dt, dx, dy);
[ut4, vt4] = TaylorDt2D(u0+    ut3*dt, v0+    vt3*dt, dx, dy);
u = u0 + (ut1+2*ut2+2*ut3+ut4)*dt/6;
v = v0 + (vt1+2*vt2+2*vt3+vt4)*dt/6;
end

function D = LaplaceOperation2(f, dx, dy)
% 2D Laplacian operator on f
% calculates the values on the edges by linearly extrapolating the second
% differences from the interior
fx = 0.*f;
fx(:,2:end-1) = (f(:,1:end-2)+f(:,3:end)-2*f(:,2:end-1))./dx^2;
fx(:,1) = 2*fx(:,2) - fx(:,3);
fx(:,end) = 2*fx(:,end-1) - fx(:,end-2);
fy = 0.*f;
fy(2:end-1,:) = (f(1:end-2,:)+f(3:end,:)-2*f(2:end-1,:))./dy^2;
fy(1,:) = 2*fy(2,:) - fy(3,:);
fy(end,:) = 2*fy(end-1,:) - fy(end-2,:);
D = fx + fy;
end

function P = PIterSolver2(u,v,ut,vt,Nu,Rho,Increx,Increy,Pinit)
% Iterative Method for Pressure Estimation
% Pressure solver for time-resolved 2D2C velocity field
% need function LaplaceOperation2
% input
% u, v - velocity in current frame, arranged in meshgrid(X,Y)
% ut, vt - time derivative of u, v
% Nu, Rho - properties of fluid
% Increx, Increy - increment of meshgrid, uniform in the field
% Pinit - initial value of pressure, can be set to 0 in the first frame
[ux, uy] = gradient(u); ux = ux./Increx; uy = uy./Increy;
[vx, vy] = gradient(v); vx = vx./Increx; vy = vy./Increy;
delta_u = LaplaceOperation2(u, Increx, Increy);
delta_v = LaplaceOperation2(v, Increx, Increy);
px = -Rho.*(ut + u.*ux + v.*uy - Nu.*delta_u);
py = -Rho.*(vt + u.*vx + v.*vy - Nu.*delta_v);

P = zeros(size(u) + 2);
P(2:end-1, 2:end-1) = Pinit;
P(1,:) = P(2,:);   P(end,:) = P(end-1,:);
P(:,1) = P(:,2);   P(:,end) = P(:,end-1);
[xq, yq] = meshgrid([1.5:size(px,2)-0.5], 1:size(px,1));
px = interp2(px, xq, yq, 'spline');
px = [px(1,:);px;px(end,:)];
px1 = [0*px(:,1:2),px,0*px(:,1)];
px2 = [0*px(:,1),px,0*px(:,1:2)];
[xq, yq] = meshgrid(1:size(py,2), [1.5:size(py,1)-0.5]);
py = interp2(py, xq, yq, 'spline');
py = [py(:,1),py,py(:,end)];
py1 = [0*py(1:2,:);py;0*py(1,:)];
py2 = [0*py(1,:);py;0*py(1:2,:)];
% Lambda = 0.1;
for iCount = 1:1e5
    Lambda = 0.2/sqrt(iCount) + 0.05; % Adaptive relaxation coefficient
    PD = Increx*px1 - (P - circshift(P, 1, 2))...
        -Increx*px2 - (P - circshift(P,-1, 2))...
        +Increy*py1 - (P - circshift(P, 1, 1))...
        -Increy*py2 - (P - circshift(P,-1, 1));
    P = P + Lambda.*PD;
    P(1,:) = P(2,:); P(end,:) = P(end-1,:);
    P(:,1) = P(:,2); P(:,end) = P(:,end-1);
    if mean(abs(PD(2:end-1,2:end-1)),'all') < 1e-8
    % if mean(abs(PD),'all') < 1e-8
        break;
    end
end
P = P(2:end-1,2:end-1);
end

function res = ResDP2D(U, V, P, Rho, Nu, X, dx, dt)
nframe = size(U,2);
res = zeros(nframe-2,1);
for iframe = 2:nframe-1
    u1 = reshape(U(:,iframe-1), size(X));
    v1 = reshape(V(:,iframe-1), size(X));
    p1 = reshape(P(:,iframe-1), size(X));
    u2 = reshape(U(:,iframe  ), size(X));
    v2 = reshape(V(:,iframe  ), size(X));
    p2 = reshape(P(:,iframe  ), size(X));
    u3 = reshape(U(:,iframe+1), size(X));
    v3 = reshape(V(:,iframe+1), size(X));
    p3 = reshape(P(:,iframe+1), size(X));
    ut = (u3-u1)./(2*dt);      utt = (u3-2*u2+u1)./dt./dt;
    vt = (v3-v1)./(2*dt);      vtt = (v3-2*v2+v1)./dt./dt; 
    pt = (p3-p1)./(2*dt);
    [ux,  uy]  = gradient(u2, dx, dx);
    [vx,  vy]  = gradient(v2, dx, dx);
    [utx, uty] = gradient(ut, dx, dx);
    [vtx, vty] = gradient(vt, dx, dx);
    [ptx, pty] = gradient(pt, dx, dx);
    delta_ut   = LaplaceOperation2(ut, dx, dx);
    delta_vt   = LaplaceOperation2(vt, dx, dx);
    resx = ptx+Rho*(utt+ut.*ux+vt.*uy+u2.*utx+v2.*uty-Nu.*delta_ut);
    resy = pty+Rho*(vtt+ut.*vx+vt.*vy+u2.*vtx+v2.*vty-Nu.*delta_vt);
    res(iframe-1) = rms([resx;resy], 'all');
end
end