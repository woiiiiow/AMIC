%% parameters
AFrame = 11:110;
Nu = 5e-5;                     Rho = 1;
dt = 0.0065;   dx = 0.0114;
xb = 0:dx:1;   yb = 0:dx:1;   zb = 0:dx:0.5;   zb = zb(1:12);
[X,Y,Z] = meshgrid(xb,yb,zb);
OutputFile = 'Pressure.mat';
dev = 'GPU'; % 'GPU'/'gpu': use GPU to calculate pressure
if isempty(gcp('nocreate')), parpool(8); end
tic

%% copy reference pressure
disp('copying reference pressure...');
SourceFolder = '../../Channel/volumn_edge/Fields_testing/';
ROI = [1 88 1 88 1 12]; % [begin end begin end begin end] of the 1/2/3 dim
P_ref = zeros(numel(X), length(AFrame));
for iFrame = 1:length(AFrame)
    sub_name = ['Field_', num2str(AFrame(iFrame),'%06u'), '.mat'];
    % loading field
    load([SourceFolder, sub_name], 'p');
    % croping
    p = p(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    P_ref(:,iFrame) = p(:);
end
if exist(OutputFile, 'file')
    save(OutputFile, 'P_ref', '-append');
else
    save(OutputFile, 'P_ref');
end
toc

%% pressure from DNS field as varification
% disp('calculating pressure from field without noise...');
% P_vari = zeros(numel(X), length(AFrame));
% for iFrame = 1:length(AFrame)
%     sub_name = ['Field_', num2str(AFrame(iFrame)-1,'%06u'), '.mat'];
%     load([SourceFolder, sub_name], 'u', 'v', 'w');
%     u1 = permute(u(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     v1 = permute(v(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     w1 = permute(w(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     sub_name = ['Field_', num2str(AFrame(iFrame),  '%06u'), '.mat'];
%     load([SourceFolder, sub_name], 'u', 'v', 'w');
%     u2 = permute(u(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     v2 = permute(v(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     w2 = permute(w(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     sub_name = ['Field_', num2str(AFrame(iFrame)+1,'%06u'), '.mat'];
%     load([SourceFolder, sub_name], 'u', 'v', 'w');
%     u3 = permute(u(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     v3 = permute(v(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     w3 = permute(w(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6)), [2 1 3]);
%     ut = (u3-u1)./(2*dt); vt = (v3-v1)./(2*dt); wt = (w3-w1)./(2*dt);
%     P = PIterSolver3(u2,v2,w2,ut,vt,wt,Nu,Rho,dx,dx,dx,0*u2,dev);
%     P = permute(P, [2 1 3]);
%     if 1 % exist('P_ref', 'var')
%         P = P - mean(P,'all') + mean(P_ref(:,iFrame),'all');
%     end
%     P_vari(:,iFrame) = P(:);
% end
% save(OutputFile, 'P_vari', '-append');
% toc

%% pressure from field with noise
disp('calculating pressure from field with noise...');
load('V_NOISE.mat', 'U_array','V_array','W_array');
P_noise = zeros(numel(X), length(AFrame));
parfor iFrame = 1:length(AFrame)
    u1 = permute(reshape(U_array(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    v1 = permute(reshape(V_array(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    w1 = permute(reshape(W_array(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    u2 = permute(reshape(U_array(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    v2 = permute(reshape(V_array(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    w2 = permute(reshape(W_array(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    u3 = permute(reshape(U_array(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    v3 = permute(reshape(V_array(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    w3 = permute(reshape(W_array(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    ut = (u3-u1)./(2*dt); vt = (v3-v1)./(2*dt); wt = (w3-w1)./(2*dt);
    P = PIterSolver3(u2,v2,w2,ut,vt,wt,Nu,Rho,dx,dx,dx,0*u2,dev);
    P = permute(P, [2 1 3]);
    if 1 % exist('P_ref', 'var')
        P = P - mean(P,'all') + mean(P_ref(:,iFrame),'all');
    end
    P_noise(:,iFrame) = P(:);
end
save(OutputFile, 'P_noise', '-append');
toc

%% pressure from Savitzky-Golay filtered field
disp('calculating pressure from Savitzky-Golay filtered field...');
load('V_SG.mat', 'U','V','W');
P_sg = zeros(numel(X), length(AFrame));
parfor iFrame = 1:length(AFrame)
    u1 = permute(reshape(U(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    v1 = permute(reshape(V(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    w1 = permute(reshape(W(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    u2 = permute(reshape(U(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    v2 = permute(reshape(V(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    w2 = permute(reshape(W(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    u3 = permute(reshape(U(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    v3 = permute(reshape(V(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    w3 = permute(reshape(W(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    ut = (u3-u1)./(2*dt); vt = (v3-v1)./(2*dt); wt = (w3-w1)./(2*dt);
    P = PIterSolver3(u2,v2,w2,ut,vt,wt,Nu,Rho,dx,dx,dx,0*u2,dev);
    P = permute(P, [2 1 3]);
    if 1 % exist('P_ref', 'var')
        P = P - mean(P,'all') + mean(P_ref(:,iFrame),'all');
    end
    P_sg(:,iFrame) = P(:);
end
save(OutputFile, 'P_sg', '-append');
toc

%% pressure from POD filtered field
disp('calculating pressure from POD filtered field...');
load('V_POD.mat', 'U','V','W');
P_pod = zeros(numel(X), length(AFrame));
parfor iFrame = 1:length(AFrame)
    u1 = permute(reshape(U(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    v1 = permute(reshape(V(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    w1 = permute(reshape(W(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    u2 = permute(reshape(U(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    v2 = permute(reshape(V(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    w2 = permute(reshape(W(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    u3 = permute(reshape(U(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    v3 = permute(reshape(V(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    w3 = permute(reshape(W(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    ut = (u3-u1)./(2*dt); vt = (v3-v1)./(2*dt); wt = (w3-w1)./(2*dt);
    P = PIterSolver3(u2,v2,w2,ut,vt,wt,Nu,Rho,dx,dx,dx,0*u2,dev);
    P = permute(P, [2 1 3]);
    if 1 % exist('P_ref', 'var')
        P = P - mean(P,'all') + mean(P_ref(:,iFrame),'all');
    end
    P_pod(:,iFrame) = P(:);
end
save(OutputFile, 'P_pod', '-append');
toc

%% pressure from AMIC filtered field
disp('calculating pressure from AMIC filtered field...');
load('V_AMIC.mat', 'U','V','W');
P_amic = zeros(numel(X), length(AFrame));
parfor iFrame = 1:length(AFrame)
    u1 = permute(reshape(U(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    v1 = permute(reshape(V(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    w1 = permute(reshape(W(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
    u2 = permute(reshape(U(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    v2 = permute(reshape(V(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    w2 = permute(reshape(W(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
    u3 = permute(reshape(U(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    v3 = permute(reshape(V(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    w3 = permute(reshape(W(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
    ut = (u3-u1)./(2*dt); vt = (v3-v1)./(2*dt); wt = (w3-w1)./(2*dt);
    P = PIterSolver3(u2,v2,w2,ut,vt,wt,Nu,Rho,dx,dx,dx,0*u2,dev);
    P = permute(P, [2 1 3]);
    if 1 % exist('P_ref', 'var')
        P = P - mean(P,'all') + mean(P_ref(:,iFrame),'all');
    end
    P_amic(:,iFrame) = P(:);
end
save(OutputFile, 'P_amic', '-append');
toc

%% pressure from autocoder filtered field
% disp('calculating pressure from AutoEncoder filtered field...');
% load('V_AE.mat', 'U','V','W');
% P_ae = zeros(numel(X), length(AFrame));
% parfor iFrame = 1:length(AFrame)
%     u1 = permute(reshape(U(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
%     v1 = permute(reshape(V(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
%     w1 = permute(reshape(W(:,AFrame(iFrame)-1), size(X)), [2 1 3]);
%     u2 = permute(reshape(U(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
%     v2 = permute(reshape(V(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
%     w2 = permute(reshape(W(:,AFrame(iFrame)  ), size(X)), [2 1 3]);
%     u3 = permute(reshape(U(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
%     v3 = permute(reshape(V(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
%     w3 = permute(reshape(W(:,AFrame(iFrame)+1), size(X)), [2 1 3]);
%     ut = (u3-u1)./(2*dt); vt = (v3-v1)./(2*dt); wt = (w3-w1)./(2*dt);
%     P = PIterSolver3(u2,v2,w2,ut,vt,wt,Nu,Rho,dx,dx,dx,0*u2,dev);
%     P = permute(P, [2 1 3]);
%     if 1 % exist('P_ref', 'var')
%         P = P - mean(P,'all') + mean(P_ref(:,iFrame),'all');
%     end
%     P_ae(:,iFrame) = P(:);
% end
% save(OutputFile, 'P_ae', '-append');
% toc

%%
function D = LaplaceOperation3(f, dx, dy, dz)
% 3D Laplacian operator on f
% calculates the values on the edges by linearly extrapolating the second
% differences from the interior
fx = zeros(size(f)); fy = fx; fz = fx;
fx(:,2:end-1,:) = (f(:,1:end-2,:)-2*f(:,2:end-1,:)+f(:,3:end,:))./dx^2;
fx(:,1,:)   = 2*fx(:,2,:)     - fx(:,3,:);
fx(:,end,:) = 2*fx(:,end-1,:) - fx(:,end-2,:);
fy(2:end-1,:,:) = (f(1:end-2,:,:)-2*f(2:end-1,:,:)+f(3:end,:,:))./dy^2;
fy(1,:,:)   = 2*fy(2,:,:)     - fy(3,:,:);
fy(end,:,:) = 2*fy(end-1,:,:) - fy(end-2,:,:);
fz(:,:,2:end-1) = (f(:,:,1:end-2)-2*f(:,:,2:end-1)+f(:,:,3:end))./dz^2;
fz(:,:,1)   = 2*fz(:,:,2)     - fz(:,:,3);
fz(:,:,end) = 2*fz(:,:,end-1) - fz(:,:,end-2);
D = fx + fy + fz;
end

function P = PIterSolver3(u,v,w,ut,vt,wt,Nu,Rho,dx,dy,dz,Pinit,dev)
% Iterative Method for Pressure Estimation
% Pressure solver for time-resolved 3D3C velocity field
% need function LaplaceOperation3
% data arranged in MATLAB default meshgrid(Y, X, Z)
% u, v, w - velocity in current frame
% ut, vt, wt - time derivative of u, v, w
% Nu, Rho - properties of fluid
% Increx, Increy, Increz - increment of meshgrid, uniform in the field
% Pinit - initial value of pressure, can be set to 0 in the first frame
% dev - computing platform, can be 'CPU' or 'GPU'
[ux,uy,uz] = gradient(u, dx,dy,dz);
[vx,vy,vz] = gradient(v, dx,dy,dz);
[wx,wy,wz] = gradient(w, dx,dy,dz);
delta_u = LaplaceOperation3(u, dx,dy,dz);
delta_v = LaplaceOperation3(v, dx,dy,dz);
delta_w = LaplaceOperation3(w, dx,dy,dz);
px = -Rho.*(ut + u.*ux + v.*uy + w.*uz - Nu.*delta_u);
py = -Rho.*(vt + u.*vx + v.*vy + w.*vz - Nu.*delta_v);
pz = -Rho.*(wt + u.*wx + v.*wy + w.*wz - Nu.*delta_w);
[ly,lx,lz] = size(u);
px = interp3(px,(1.5:lx-0.5)',1:ly,1:lz);
py = interp3(py,1:lx,(1.5:ly-0.5)',1:lz);
pz = interp3(pz,1:lx,1:ly,(1.5:lz-0.5)');
P = zeros(size(u) + 2);
if strcmp(dev,'GPU')|strcmp(dev,'gpu')
    px = gpuArray(px); py = gpuArray(py); pz = gpuArray(pz);
    P = gpuArray(P);
end
px = cat(1,px(1,:,:),px,px(end,:,:));
px = cat(3,px(:,:,1),px,px(:,:,end));
px1= cat(2,0*px(:,1:2,:),px,0*px(:,1,:));
px2= cat(2,0*px(:,1,:),px,0*px(:,1:2,:));
py = cat(2,py(:,1,:),py,py(:,end,:));
py = cat(3,py(:,:,1),py,py(:,:,end));
py1= cat(1,0*py(1:2,:,:),py,0*py(1,:,:));
py2= cat(1,0*py(1,:,:),py,0*py(1:2,:,:));
pz = cat(1,pz(1,:,:),pz,pz(end,:,:));
pz = cat(2,pz(:,1,:),pz,pz(:,end,:));
pz1= cat(3,0*pz(:,:,1:2),pz,0*pz(:,:,1));
pz2= cat(3,0*pz(:,:,1),pz,0*pz(:,:,1:2));
P(2:end-1, 2:end-1, 2:end-1) = Pinit;
P(1,:,:) = P(2,:,:);   P(end,:,:) = P(end-1,:,:);
P(:,1,:) = P(:,2,:);   P(:,end,:) = P(:,end-1,:);
P(:,:,1) = P(:,:,2);   P(:,:,end) = P(:,:,end-1);
for iCount = 1:1e5
    % Adaptive relaxation coefficient
    Lambda = 0.17/sqrt((iCount-1)/10+1) + 0.08;
    PD = dx*px1 - (P - circshift(P, 1, 2))...
        -dx*px2 - (P - circshift(P,-1, 2))...
        +dy*py1 - (P - circshift(P, 1, 1))...
        -dy*py2 - (P - circshift(P,-1, 1))...
        +dz*pz1 - (P - circshift(P, 1, 3))...
        -dz*pz2 - (P - circshift(P,-1, 3));
    P = P + Lambda.*PD;
    P(1,:,:) = P(2,:,:); P(end,:,:) = P(end-1,:,:);
    P(:,1,:) = P(:,2,:); P(:,end,:) = P(:,end-1,:);
    P(:,:,1) = P(:,:,2); P(:,:,end) = P(:,:,end-1);
    if mean(abs(PD(2:end-1,2:end-1,2:end-1)),'all') < 1e-8
        break;
    end
end
P = P(2:end-1,2:end-1,2:end-1);
end