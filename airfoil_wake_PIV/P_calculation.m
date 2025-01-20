% woii EAP UC3M junwei.chen@uc3m.es
% 20241017

SourceFolder = '../../Experiment_wing3/OUT_TRPIV/';
OutputFile = 'Pressure.mat';

ROI = [41 110 41 150];         % subdomain of PIV field
                               % [first-last row, first-last column]
Mm_per_px_ratio = 0.1197;      % spatial resolution of photos
Sample_rate = 30;              % sample rate of PIV
Vector_Spacing = 10;           % vector per pixel in PIV
Nu = 1e-6;                     % kinetic viscosity coefficient of water
Rho = 1e3;                     % density of water

AFrame = 11:530;               % frames to be processed
AFrame_P = 21:260;             % frames to calculate pressure

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
Increx = Vector_Spacing*Mm_per_px_ratio*1e-3;
Increy = Increx;
Incret = 1/Sample_rate;

if ~exist(OutputFile, 'file')
    save(OutputFile, 'AFrame_P');
else
    save(OutputFile, 'AFrame_P', '-append');
end

%% pressure from original PIV field
tic
V_PIV.U = zeros(numel(X_etr), length(AFrame));   V_PIV.V = V_PIV.U;
iCount = 0;
for iFrame = AFrame
    iCount = iCount + 1;
    load([SourceFolder,'Wing_',num2str(iFrame,'%06u'),'.mat'], 'U','V');
    V_PIV.U(:, iCount) = V(:);    V_PIV.V(:, iCount) = U(:);
    % swapping U and V
end
V_PIV.U = V_PIV.U.*Sample_rate.*Mm_per_px_ratio.*1e-3;
V_PIV.V = V_PIV.V.*Sample_rate.*Mm_per_px_ratio.*1e-3;
clear iFrame iCount U V
P_PIV = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u1 = reshape(V_PIV.U(map(:),p_iFrame-1), size(X));
    v1 = reshape(V_PIV.V(map(:),p_iFrame-1), size(X));
    u2 = reshape(V_PIV.U(map(:),p_iFrame  ), size(X));
    v2 = reshape(V_PIV.V(map(:),p_iFrame  ), size(X));
    u3 = reshape(V_PIV.U(map(:),p_iFrame+1), size(X));
    v3 = reshape(V_PIV.V(map(:),p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*Incret);
    vt = (v3 - v1)./(2*Incret);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,Increx,Increy,0*X);
    P_PIV(:,iFrame) = P(:);
end
[xcor0, ycor0] = meshgrid(1:size(X,2),1:size(X,1));
for iFrame = 2:length(AFrame_P)
    % % minimize dp/dt over the whole domain
    % p1 = reshape(P_PIV(:,iFrame-1), size(X));
    % p2 = reshape(P_PIV(:,iFrame  ), size(X));
    % p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    % u1 = reshape(V_PIV.U(map(:),p_iFrame-1), size(X));
    % v1 = reshape(V_PIV.V(map(:),p_iFrame-1), size(X));
    % u2 = reshape(V_PIV.U(map(:),p_iFrame  ), size(X));
    % v2 = reshape(V_PIV.V(map(:),p_iFrame  ), size(X));
    % [p1x, p1y] = gradient(p1); p1x = p1x./Increx; p1y = p1y./Increy;
    % [p2x, p2y] = gradient(p2); p2x = p2x./Increx; p2y = p2y./Increy;
    % P_PIV(:,iFrame) = P_PIV(:,iFrame) -...
    %     mean(p2-p1+(u1.*p1x+v1.*p1y+u2.*p2x+v2.*p2y).*Incret./2, 'all');

    % the position of grid point at frame 0.5 moving to frame 0
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1) - 1;
    u = reshape(V_PIV.U(map(:), p_iFrame), size(X));
    v = reshape(V_PIV.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor1 = xcor0 - u./2./Vector_Spacing;
    ycor1 = ycor0 - v./2./Vector_Spacing;
    % the position of grid point at frame 0.5 moving to frame 1
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u = reshape(V_PIV.U(map(:), p_iFrame), size(X));
    v = reshape(V_PIV.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor2 = xcor0 + u./2./Vector_Spacing;
    ycor2 = ycor0 + v./2./Vector_Spacing;
    % minimize dp/dt over the available grid points
    p1 = interp2(reshape(P_PIV(:,iFrame-1),size(xcor1)), xcor1, ycor1);
    p2 = interp2(reshape(P_PIV(:,iFrame),  size(xcor2)), xcor2, ycor2);
    P_PIV(:,iFrame) = P_PIV(:,iFrame)-mean(p2(:)-p1(:),'all','omitnan');
end
save(OutputFile, 'P_PIV', '-append');
toc

%% pressure from Savitzky-Golay filtered field
tic
V_SG = load('V_SG.mat');
P_SG = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - V_SG.AFrame(1) + AFrame_P(1);
    u1 = reshape(V_SG.U(map(:),p_iFrame-1), size(X));
    v1 = reshape(V_SG.V(map(:),p_iFrame-1), size(X));
    u2 = reshape(V_SG.U(map(:),p_iFrame  ), size(X));
    v2 = reshape(V_SG.V(map(:),p_iFrame  ), size(X));
    u3 = reshape(V_SG.U(map(:),p_iFrame+1), size(X));
    v3 = reshape(V_SG.V(map(:),p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*Incret);
    vt = (v3 - v1)./(2*Incret);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,Increx,Increy,0*X);
    P_SG(:,iFrame) = P(:);
end
[xcor0, ycor0] = meshgrid(1:size(X,2),1:size(X,1));
for iFrame = 2:length(AFrame_P)
    % the position of grid point at frame 0.5 moving to frame 0
    p_iFrame = iFrame - V_SG.AFrame(1) + AFrame_P(1) - 1;
    u = reshape(V_SG.U(map(:), p_iFrame), size(X));
    v = reshape(V_SG.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor1 = xcor0 - u./2./Vector_Spacing;
    ycor1 = ycor0 - v./2./Vector_Spacing;
    % the position of grid point at frame 0.5 moving to frame 1
    p_iFrame = iFrame - V_SG.AFrame(1) + AFrame_P(1);
    u = reshape(V_SG.U(map(:), p_iFrame), size(X));
    v = reshape(V_SG.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor2 = xcor0 + u./2./Vector_Spacing;
    ycor2 = ycor0 + v./2./Vector_Spacing;
    % minimize dp/dt over the available grid points
    p1 = interp2(reshape(P_SG(:,iFrame-1),size(xcor1)), xcor1, ycor1);
    p2 = interp2(reshape(P_SG(:,iFrame),  size(xcor2)), xcor2, ycor2);
    P_SG(:,iFrame) = P_SG(:,iFrame)-mean(p2(:)-p1(:),'all','omitnan');
end
save(OutputFile, 'P_SG', '-append');
toc

%% pressure from POD processed field
tic
V_POD = load('V_POD.mat');
P_POD = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - V_POD.AFrame(1) + AFrame_P(1);
    u1 = reshape(V_POD.U(map(:),p_iFrame-1), size(X));
    v1 = reshape(V_POD.V(map(:),p_iFrame-1), size(X));
    u2 = reshape(V_POD.U(map(:),p_iFrame  ), size(X));
    v2 = reshape(V_POD.V(map(:),p_iFrame  ), size(X));
    u3 = reshape(V_POD.U(map(:),p_iFrame+1), size(X));
    v3 = reshape(V_POD.V(map(:),p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*Incret);
    vt = (v3 - v1)./(2*Incret);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,Increx,Increy,0*X);
    P_POD(:,iFrame) = P(:);
end
[xcor0, ycor0] = meshgrid(1:size(X,2),1:size(X,1));
for iFrame = 2:length(AFrame_P)
    % the position of grid point at frame 0.5 moving to frame 0
    p_iFrame = iFrame - V_POD.AFrame(1) + AFrame_P(1) - 1;
    u = reshape(V_POD.U(map(:), p_iFrame), size(X));
    v = reshape(V_POD.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor1 = xcor0 - u./2./Vector_Spacing;
    ycor1 = ycor0 - v./2./Vector_Spacing;
    % the position of grid point at frame 0.5 moving to frame 1
    p_iFrame = iFrame - V_POD.AFrame(1) + AFrame_P(1);
    u = reshape(V_POD.U(map(:), p_iFrame), size(X));
    v = reshape(V_POD.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor2 = xcor0 + u./2./Vector_Spacing;
    ycor2 = ycor0 + v./2./Vector_Spacing;
    % minimize dp/dt over the available grid points
    p1 = interp2(reshape(P_POD(:,iFrame-1),size(xcor1)), xcor1, ycor1);
    p2 = interp2(reshape(P_POD(:,iFrame),  size(xcor2)), xcor2, ycor2);
    P_POD(:,iFrame) = P_POD(:,iFrame)-mean(p2(:)-p1(:),'all','omitnan');
end
save(OutputFile, 'P_POD', '-append');
toc

%% pressure from Auto Encoder processed field
tic
V_AE = load('V_AE.mat');
P_AE = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - V_AE.AFrame(1) + AFrame_P(1);
    u1 = reshape(V_AE.U(:,p_iFrame-1), size(X));
    v1 = reshape(V_AE.V(:,p_iFrame-1), size(X));
    u2 = reshape(V_AE.U(:,p_iFrame  ), size(X));
    v2 = reshape(V_AE.V(:,p_iFrame  ), size(X));
    u3 = reshape(V_AE.U(:,p_iFrame+1), size(X));
    v3 = reshape(V_AE.V(:,p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*Incret);
    vt = (v3 - v1)./(2*Incret);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,Increx,Increy,0*X);
    P_AE(:,iFrame) = P(:);
end
[xcor0, ycor0] = meshgrid(1:size(X,2),1:size(X,1));
for iFrame = 2:length(AFrame_P)
    % the position of grid point at frame 0.5 moving to frame 0
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1) - 1;
    u = reshape(V_AE.U(:, p_iFrame), size(X));
    v = reshape(V_AE.V(:, p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor1 = xcor0 - u./2./Vector_Spacing;
    ycor1 = ycor0 - v./2./Vector_Spacing;
    % the position of grid point at frame 0.5 moving to frame 1
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u = reshape(V_AE.U(:, p_iFrame), size(X));
    v = reshape(V_AE.V(:, p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor2 = xcor0 + u./2./Vector_Spacing;
    ycor2 = ycor0 + v./2./Vector_Spacing;
    % minimize dp/dt over the available grid points
    p1 = interp2(reshape(P_AE(:,iFrame-1),size(xcor1)), xcor1, ycor1);
    p2 = interp2(reshape(P_AE(:,iFrame),  size(xcor2)), xcor2, ycor2);
    P_AE(:,iFrame) = P_AE(:,iFrame)-mean(p2(:)-p1(:),'all','omitnan');
end
save(OutputFile, 'P_AE', '-append');
toc

%% pressure from Advection Multiframe Iterative Correction processed field
tic
V_AMIC = load('V_AMIC.mat');
P_AMIC = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - V_AMIC.AFrame(1) + AFrame_P(1);
    u1 = reshape(V_AMIC.U(map(:),p_iFrame-1), size(X));
    v1 = reshape(V_AMIC.V(map(:),p_iFrame-1), size(X));
    u2 = reshape(V_AMIC.U(map(:),p_iFrame  ), size(X));
    v2 = reshape(V_AMIC.V(map(:),p_iFrame  ), size(X));
    u3 = reshape(V_AMIC.U(map(:),p_iFrame+1), size(X));
    v3 = reshape(V_AMIC.V(map(:),p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*Incret);
    vt = (v3 - v1)./(2*Incret);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,Increx,Increy,0*X);
    P_AMIC(:,iFrame) = P(:);
end
[xcor0, ycor0] = meshgrid(1:size(X,2),1:size(X,1));
for iFrame = 2:length(AFrame_P)
    % the position of grid point at frame 0.5 moving to frame 0
    p_iFrame = iFrame - V_AMIC.AFrame(1) + AFrame_P(1) - 1;
    u = reshape(V_AMIC.U(map(:), p_iFrame), size(X));
    v = reshape(V_AMIC.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor1 = xcor0 - u./2./Vector_Spacing;
    ycor1 = ycor0 - v./2./Vector_Spacing;
    % the position of grid point at frame 0.5 moving to frame 1
    p_iFrame = iFrame - V_AMIC.AFrame(1) + AFrame_P(1);
    u = reshape(V_AMIC.U(map(:), p_iFrame), size(X));
    v = reshape(V_AMIC.V(map(:), p_iFrame), size(X));
    u = u./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    v = v./(Sample_rate.*Mm_per_px_ratio.*1e-3);
    xcor2 = xcor0 + u./2./Vector_Spacing;
    ycor2 = ycor0 + v./2./Vector_Spacing;
    % minimize dp/dt over the available grid points
    p1 = interp2(reshape(P_AMIC(:,iFrame-1),size(xcor1)), xcor1, ycor1);
    p2 = interp2(reshape(P_AMIC(:,iFrame),  size(xcor2)), xcor2, ycor2);
    P_AMIC(:,iFrame) = P_AMIC(:,iFrame)-mean(p2(:)-p1(:),'all','omitnan');
end
save(OutputFile, 'P_AMIC', '-append');
toc




%% functions
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
