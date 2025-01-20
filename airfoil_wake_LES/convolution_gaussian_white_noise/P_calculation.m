% 240617 computing pressure
OutputFile = 'Pressure.mat';
AFrame_P = 11:310;             % frames to calculate pressure
tic

% loading flow parameters
file = '../../../airfoilLES/airfoilLES_parameters.h5';
Re = h5read(file, '/Re');
dt = h5read(file, '/dt');
% Ma = 0.3;
% Re = 23000;
Rho = 1.2;                     % assume to be so
Nu = 1/Re;                     % non-dimensional far field velocity is 1 
                               % and chord length is 1

% load velocity field reference and noise
load('V_NOISE.mat');
dx = xb(2) - xb(1);
[X, Y] = meshgrid(xb, yb);
AFrame = 1:size(U_clean, 2);

if ~exist(OutputFile, 'file')
    save(OutputFile, 'AFrame_P');
else
    save(OutputFile, 'AFrame_P', '-append');
end

%% pressure from original LES field
disp('computing reference pressure...');
P_REF = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u1 = reshape(U_clean(:,p_iFrame-1), size(X));
    v1 = reshape(V_clean(:,p_iFrame-1), size(X));
    u2 = reshape(U_clean(:,p_iFrame  ), size(X));
    v2 = reshape(V_clean(:,p_iFrame  ), size(X));
    u3 = reshape(U_clean(:,p_iFrame+1), size(X));
    v3 = reshape(V_clean(:,p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*dt);
    vt = (v3 - v1)./(2*dt);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,dx,dt,0*X);
    P_REF(:,iFrame) = P(:);
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
    u = reshape(U_clean(:, p_iFrame), size(X));
    v = reshape(V_clean(:, p_iFrame), size(X));
    xcor1 = xcor0 - u.*dt./2;
    ycor1 = ycor0 - v.*dt./2;
    % the position of grid point at frame 0.5 moving to frame 1
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u = reshape(U_clean(:, p_iFrame), size(X));
    v = reshape(V_clean(:, p_iFrame), size(X));
    xcor2 = xcor0 + u.*dt./2;
    ycor2 = ycor0 + v.*dt./2;
    % minimize dp/dt over the available grid points
    p1 = interp2(reshape(P_REF(:,iFrame-1),size(xcor1)), xcor1, ycor1);
    p2 = interp2(reshape(P_REF(:,iFrame),  size(xcor2)), xcor2, ycor2);
    P_REF(:,iFrame) = P_REF(:,iFrame)-mean(p2(:)-p1(:),'all','omitnan');
end
save(OutputFile, 'P_REF', '-append');
toc

%% pressure from noisy field
disp('computing pressure from noisy field...');
P_NOISY = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u1 = reshape(U_noisy(:,p_iFrame-1), size(X));
    v1 = reshape(V_noisy(:,p_iFrame-1), size(X));
    u2 = reshape(U_noisy(:,p_iFrame  ), size(X));
    v2 = reshape(V_noisy(:,p_iFrame  ), size(X));
    u3 = reshape(U_noisy(:,p_iFrame+1), size(X));
    v3 = reshape(V_noisy(:,p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*dt);
    vt = (v3 - v1)./(2*dt);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,dx,dt,0*X);
    P_NOISY(:,iFrame)...
        = P(:)- mean(P, 'all') + mean(P_REF(:,iFrame), 'all');
end
save(OutputFile, 'P_NOISY', '-append');
toc

%% pressure from Savitzky-Golay filtered field
disp('computing pressure from Savitzky-Golay filtered field...');
load('V_SG.mat');
P_SG = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u1 = reshape(U(:,p_iFrame-1), size(X));
    v1 = reshape(V(:,p_iFrame-1), size(X));
    u2 = reshape(U(:,p_iFrame  ), size(X));
    v2 = reshape(V(:,p_iFrame  ), size(X));
    u3 = reshape(U(:,p_iFrame+1), size(X));
    v3 = reshape(V(:,p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*dt);
    vt = (v3 - v1)./(2*dt);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,dx,dt,0*X);
    P_SG(:,iFrame)...
        = P(:)- mean(P, 'all') + mean(P_REF(:,iFrame), 'all');
end
save(OutputFile, 'P_SG', '-append');
toc

%% pressure from POD truncated field
disp('computing pressure from POD truncated field...');
load('V_POD.mat');
P_POD = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u1 = reshape(U(:,p_iFrame-1), size(X));
    v1 = reshape(V(:,p_iFrame-1), size(X));
    u2 = reshape(U(:,p_iFrame  ), size(X));
    v2 = reshape(V(:,p_iFrame  ), size(X));
    u3 = reshape(U(:,p_iFrame+1), size(X));
    v3 = reshape(V(:,p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*dt);
    vt = (v3 - v1)./(2*dt);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,dx,dt,0*X);
    P_POD(:,iFrame)...
        = P(:)- mean(P, 'all') + mean(P_REF(:,iFrame), 'all');
end
save(OutputFile, 'P_POD', '-append');
toc

%% pressure from advection-based multiframe iterativly corrected field
disp('computing pressure from AMIC processed field...');
load('V_AMIC.mat');
P_AMIC = zeros(numel(X), length(AFrame_P));
parfor iFrame = 1:length(AFrame_P)
    p_iFrame = iFrame - AFrame(1) + AFrame_P(1);
    u1 = reshape(U(:,p_iFrame-1), size(X));
    v1 = reshape(V(:,p_iFrame-1), size(X));
    u2 = reshape(U(:,p_iFrame  ), size(X));
    v2 = reshape(V(:,p_iFrame  ), size(X));
    u3 = reshape(U(:,p_iFrame+1), size(X));
    v3 = reshape(V(:,p_iFrame+1), size(X));
    ut = (u3 - u1)./(2*dt);
    vt = (v3 - v1)./(2*dt);
    P  = PIterSolver2(u2,v2,ut,vt,Nu,Rho,dx,dt,0*X);
    P_AMIC(:,iFrame)...
        = P(:)- mean(P, 'all') + mean(P_REF(:,iFrame), 'all');
end
save(OutputFile, 'P_AMIC', '-append');
toc

%% test play
% figure;
% for iFrame = 1:300
% % pcolor(X, Y, reshape(P_REF(:,iFrame),size(X)));
% % pcolor(X, Y, reshape(P_NOISY(:,iFrame),size(X)));
% % pcolor(X, Y, reshape(P_SG(:,iFrame),size(X)));
% pcolor(X, Y, reshape(P_POD(:,iFrame),size(X)));
% % pcolor(X, Y, reshape(P_AMIC(:,iFrame),size(X)));
% shading flat; colormap jet; caxis([-0.1 0.1]); pause(0.1);
% end

%% error
err = std(P_NOISY-P_REF, 1, 'all');
disp(['std error of noisy field: ', num2str(err)]);
err = std(P_SG   -P_REF, 1, 'all');
disp(['std error of SG filtered field: ', num2str(err)]);
err = std(P_POD  -P_REF, 1, 'all');
disp(['std error of POD truncated field: ', num2str(err)]);
err = std(P_AMIC -P_REF, 1, 'all');
disp(['std error of AMIC corrected field: ', num2str(err)]);

%% cosine similarity
simi = mean(sum(P_REF.*P_NOISY)./sqrt(sum(P_REF.^2).*sum(P_NOISY.^2)));
disp(['cosine similarity of noisy field: ', num2str(simi)]);
simi = mean(sum(P_REF.*P_SG   )./sqrt(sum(P_REF.^2).*sum(P_SG   .^2)));
disp(['cosine similarity of SG filtered field: ', num2str(simi)]);
simi = mean(sum(P_REF.*P_POD  )./sqrt(sum(P_REF.^2).*sum(P_POD .^2)));
disp(['cosine similarity of POD truncated field: ', num2str(simi)]);
simi = mean(sum(P_REF.*P_AMIC )./sqrt(sum(P_REF.^2).*sum(P_AMIC .^2)));
disp(['cosine similarity of AMIC corrected field: ', num2str(simi)]);

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

function P = PIterSolver2Compressible(u,v,ut,vt,Nu,Rho,dx,dy,Pinit)
% Iterative Method for Pressure Estimation
% Pressure solver for time-resolved 2D2C velocity field
% compressible
% need function LaplaceOperation2
% input
% u, v - velocity in current frame, arranged in meshgrid(X,Y)
% ut, vt - time derivative of u, v
% Nu, Rho - properties of fluid
% Increx, Increy - increment of meshgrid, uniform in the field
% Pinit - initial value of pressure, can be set to 0 in the first frame
[ux, uy] = gradient(u); ux = ux./dx; uy = uy./dy;
[vx, vy] = gradient(v); vx = vx./dx; vy = vy./dy;
delta_u = LaplaceOperation2(u, dx, dy);
delta_v = LaplaceOperation2(v, dx, dy);
[divx, divy] = gradient(divergence(u, v)./dx); % supposing dx = dy
divx = divx./dx; divy = divy./dy;
px = -Rho.*(ut + u.*ux + v.*uy - Nu.*(delta_u+1/3.*divx));
py = -Rho.*(vt + u.*vx + v.*vy - Nu.*(delta_v+1/3.*divy));

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
for iCount = 1:5e5
    Lambda = 0.2/sqrt(iCount) + 0.05; % Adaptive relaxation coefficient
    PD = dx*px1 - (P - circshift(P, 1, 2))...
        -dx*px2 - (P - circshift(P,-1, 2))...
        +dy*py1 - (P - circshift(P, 1, 1))...
        -dy*py2 - (P - circshift(P,-1, 1));
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