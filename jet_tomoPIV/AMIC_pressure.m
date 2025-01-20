% woii EAP UC3M junwei.chen@uc3m.es
% 20241017
% setting of loading and saving data
RootFolder = './';
SubSet     = 'JetAR1_02';
PreFix     = '/GPU_vec';
SufFix     = '.dat';
SaveFolder = './';
AFrame     = 7:4:5995;          % the frames to process (both AMIC and pressure computation)
% geometries
dt         = 1/500;             % original set 2000 Hz
% setting of AMIC
Lambda     = [0.15 0.15];       % coefficient of AMIC smoothing
NLoop      = 64;                % maximum number of loop
ThVCI      = 0.9;               % iteration stops when the velocity
                                % correction index exceeds the threshold
                                % default: 0.9
% setting of pressure computation
FIncre     = 1;                 % multiples of time increment for pressure
Nu         = 1e-6;              % kinetic viscosity coefficient of water
Rho        = 1e3;               % density of water
dev        = 'GPU';             % 'GPU'/'gpu': use GPU to calculate pressure
ParNo      = 4;                 % number of parallel pool

disp(['save to file "',SaveFolder,SubSet,'_conAMIC.mat"']);
tic

%% loading data
disp('Loading data...');
A = importdata([RootFolder,SubSet,...
    PreFix,num2str(AFrame(1),'%05u'),SufFix]);
tmp = regexp(A.textdata{2}, '\d+', 'match');
Size= [str2num(tmp{1}) str2num(tmp{2}) str2num(tmp{3})];
X = reshape(A.data(:,2), Size);
Y = reshape(A.data(:,1), Size);
Z = reshape(A.data(:,3), Size);
dx = (X(1,2,1) - X(1)).*1e-3;
dy = (Y(2,1,1) - Y(1)).*1e-3;
dz = (Z(1,1,2) - Z(1)).*1e-3;
U_array = zeros(numel(X), length(AFrame));
V_array = zeros(numel(X), length(AFrame));
W_array = zeros(numel(X), length(AFrame));
Enabled = false(numel(X), length(AFrame)); % enabled vectors
for iFrame = 1:length(AFrame)
    A = importdata([RootFolder,SubSet,PreFix,...
        num2str(AFrame(iFrame),'%05u'),SufFix]);
    U = reshape(A.data(:,5), Size);
    V = reshape(A.data(:,4), Size);
    W =-reshape(A.data(:,6), Size);
    E = ~(V > 0.0076 & V < 0.00765);
    U_array(:,iFrame) = U(:);
    V_array(:,iFrame) = V(:);
    W_array(:,iFrame) = W(:);
    Enabled(:,iFrame) = E(:);
end
U_noisy = U_array;      V_noisy = V_array;      W_noisy = W_array;
toc

%% AMIC (advection-based multiframe iterative correction) for velocity
disp('Smoothing the velocity field using AMIC...');
disp(['AMIC (advection-based multiframe iterative correction)',...
    ' is amicable!']);
A_Vdis = zeros(1,NLoop); % L-2 distance of velocity to the noisy field
U_noise = U_array;  V_noise = V_array;  W_noise = W_array;
v_cor_idx = 0;                 % velocity correction index
for iLoop = 1:NLoop
    % turbulence frozen filter
    % sequence frame P+1 P-1 P+2 P-2 P+3 P-3...
    U_cor = zeros([size(U_array) length(Lambda)*2]);
    V_cor = U_cor;             W_cor = U_cor;
    for iFrame = 1:length(AFrame)
        u0 = reshape(U_array(:,iFrame), size(X));
        v0 = reshape(V_array(:,iFrame), size(X));
        w0 = reshape(W_array(:,iFrame), size(X));
        mask0 = double(reshape(Enabled(:,iFrame), size(X)));
        u = u0;   v = v0;   w = w0;   mask = mask0;
        for iStep = 1:length(Lambda)
            [u, v, w] = MRK4Prop3D(u, v, w, dx, dy, dz, dt);
            mask = interp3(X,Y,Z, mask, X-u*dt,Y-v*dt,Z-w*dt);
            mask(isnan(mask)) = 0; mask = imgaussfilt3(mask, 5);
            if iFrame + iStep < length(AFrame)
                U_cor(:,iFrame+iStep,iStep*2-1)...
                    = (u(:) - U_array(:,iFrame+iStep)).*mask(:);
                V_cor(:,iFrame+iStep,iStep*2-1)...
                    = (v(:) - V_array(:,iFrame+iStep)).*mask(:);
                W_cor(:,iFrame+iStep,iStep*2-1)...
                    = (w(:) - W_array(:,iFrame+iStep)).*mask(:);
            end
        end
        u = u0;   v = v0;   w = w0;   mask = mask0;
        for iStep = 1:length(Lambda)
            [u, v, w] = MRK4Prop3D(u, v, w, dx, dy, dz,-dt);
            mask = interp3(X,Y,Z, mask, X+u*dt,Y+v*dt,Z+w*dt);
            mask(isnan(mask)) = 0; mask = imgaussfilt3(mask, 5);
            if iFrame - iStep >0
                U_cor(:,iFrame-iStep,iStep*2)...
                    = (u(:) - U_array(:,iFrame-iStep)).*mask(:);
                V_cor(:,iFrame-iStep,iStep*2)...
                    = (v(:) - V_array(:,iFrame-iStep)).*mask(:);
                W_cor(:,iFrame-iStep,iStep*2)...
                    = (w(:) - W_array(:,iFrame-iStep)).*mask(:);
            end
        end
    end
    % ratio of update
    index_cth =...
        rms(U_cor(U_cor~=0))+rms(V_cor(V_cor~=0))+rms(W_cor(W_cor~=0));
    % correcting
    for iStep = 1:length(Lambda)
        U_array = U_array...
            + Lambda(iStep)*(U_cor(:,:,iStep*2-1)+U_cor(:,:,iStep*2));
        V_array = V_array...
            + Lambda(iStep)*(V_cor(:,:,iStep*2-1)+V_cor(:,:,iStep*2));
        W_array = W_array...
            + Lambda(iStep)*(W_cor(:,:,iStep*2-1)+W_cor(:,:,iStep*2));
    end
    % % divergence free filter
    % if iLoop >= SLoopDiv
    %     div = 0;
    %     for iFrame = 1:length(AFrame)
    %         u0 = reshape(U_array(:,iFrame), size(X));
    %         v0 = reshape(V_array(:,iFrame), size(X));
    %         w0 = reshape(W_array(:,iFrame), size(X));
    %         [u, v, w] = DivFilter3D(u0, v0, w0, RateDiv, NItDiv);
    %         U_array(:,iFrame) = u(:);
    %         V_array(:,iFrame) = v(:);
    %         W_array(:,iFrame) = w(:);
    %         div = div+rms(u-u0,'all')+rms(v-v0,'all')+rms(w-w0,'all');
    %     end
    %     index_cdiv = div/3/length(AFrame)/Vrms;
    % end
    A_Vdis(iLoop) =...
        rms([U_array-U_noise;V_array-V_noise;W_array-W_noise],'all');
    fprintf(['iLoop = ',num2str(iLoop),...
        '\tcorrection advection = ',num2str(index_cth,'%.5g'),'\n']);
    if iLoop > 2
        v_cor_idx = (A_Vdis(iLoop)-A_Vdis(iLoop-1))/...
            (A_Vdis(iLoop-1)-A_Vdis(iLoop-2));
        fprintf(['\b\tvelocity correction index = ',...
            num2str(v_cor_idx,'%.5g'),'\n']);
    end
    if v_cor_idx > ThVCI
        break
    end
    toc
end

%% computing pressure
disp('Computing pressure field...');
P_array = zeros(numel(X), length(AFrame));
if isempty(gcp('nocreate')), parpool(ParNo); end
parfor iFrame = 1+FIncre:length(AFrame)-FIncre
    disp(iFrame);
    u1 = reshape(U_array(:,iFrame-FIncre), size(X));
    v1 = reshape(V_array(:,iFrame-FIncre), size(X));
    w1 = reshape(W_array(:,iFrame-FIncre), size(X));
    e1 = reshape(Enabled(:,iFrame-FIncre), size(X));
    u2 = reshape(U_array(:,iFrame       ), size(X));
    v2 = reshape(V_array(:,iFrame       ), size(X));
    w2 = reshape(W_array(:,iFrame       ), size(X));
    e2 = reshape(Enabled(:,iFrame       ), size(X));
    u3 = reshape(U_array(:,iFrame+FIncre), size(X));
    v3 = reshape(V_array(:,iFrame+FIncre), size(X));
    w3 = reshape(W_array(:,iFrame+FIncre), size(X));
    e3 = reshape(Enabled(:,iFrame+FIncre), size(X));
    % e = e1&e2&e3;
    e = imclose(e1,ones(3)) & imclose(e2,ones(3)) & imclose(e3,ones(3));
    ut = (u3-u1)./(2*dt*FIncre);
    vt = (v3-v1)./(2*dt*FIncre);
    wt = (w3-w1)./(2*dt*FIncre);
    P = PIterSolver3(u2,v2,w2,ut,vt,wt,e,Nu,Rho,dx,dy,dz,0*u2,dev);
    P_array(:,iFrame) = P(:);
end
% correcting interframe pressure
[xcor0, ycor0, zcor0] = meshgrid(1:size(X,2),1:size(X,1),1:size(X,3));
for iFrame = 1+FIncre:length(AFrame)-FIncre+1
    % the position of grid point at frame 0.5 moving to frame 0
    p_iFrame = iFrame - 1;
    % velocity unit switch from m/s to edge length/frame
    u = reshape(U_array(:, p_iFrame), size(X)).*(dt/dx);
    v = reshape(V_array(:, p_iFrame), size(X)).*(dt/dx);
    w = reshape(W_array(:, p_iFrame), size(X)).*(dt/dx);
    xcor1 = xcor0 - u.*dt./2;
    ycor1 = ycor0 - v.*dt./2;
    zcor1 = zcor0 - w.*dt./2;
    % the position of grid point at frame 0.5 moving to frame 1
    p_iFrame = iFrame;
    u = reshape(U_array(:, p_iFrame), size(X)).*(dt/dx);
    v = reshape(V_array(:, p_iFrame), size(X)).*(dt/dx);
    w = reshape(W_array(:, p_iFrame), size(X)).*(dt/dx);
    xcor2 = xcor0 + u.*dt./2;
    ycor2 = ycor0 + v.*dt./2;
    zcor2 = zcor0 + w.*dt./2;
    % minimize dp/dt over the available grid points
    p1 = interp3(reshape(P_array(:,iFrame-1),size(X)), xcor1,ycor1,zcor1);
    p2 = interp3(reshape(P_array(:,iFrame),  size(X)), xcor2,ycor2,zcor2);
    P_array(:,iFrame)...
        = P_array(:,iFrame) - mean(p2(:)-p1(:),'all','omitnan');
end
toc

%% saving data
save([SaveFolder,SubSet,'_conAMIC.mat'], 'X', 'Y', 'Z',...
    'U_noisy', 'V_noisy', 'W_noisy', 'Enabled',...
    'U_array', 'V_array', 'W_array', 'P_array', '-v7.3');
toc

%% functions
function [ut, vt, wt] = TaylorDt3D(u, v, w, dx, dy, dz)
% ut = -(u_mean\cdot\nabla)u_fluc
u_mean = imgaussfilt3(u, 7);
v_mean = imgaussfilt3(v, 7);
w_mean = imgaussfilt3(w, 7);
u_fluc = u - u_mean;
v_fluc = v - v_mean;
w_fluc = w - w_mean;
[ux, uy, uz] = gradient(u_fluc);
[vx, vy, vz] = gradient(v_fluc);
[wx, wy, wz] = gradient(w_fluc);
ux = ux./dx;    uy = uy./dy;    uz = uz./dz;
vx = vx./dx;    vy = vy./dy;    vz = vz./dz;
wx = wx./dx;    wy = wy./dy;    wz = wz./dz;
ut = - u_mean.*ux - v_mean.*uy - w_mean.*uz;
vt = - u_mean.*vx - v_mean.*vy - w_mean.*vz;
wt = - u_mean.*wx - v_mean.*wy - w_mean.*wz;
ut = imgaussfilt3(ut, 0.75);
vt = imgaussfilt3(vt, 0.75);
wt = imgaussfilt3(wt, 0.75);
end

function [u, v, w] = MRK4Prop3D(u0, v0, w0, dx, dy, dz, dt)
% Runge-Kutta propogation of 3D flow fields
% u0, v0, w0, current velocity
% x, y, z corrodinates
% dt time increments
% u, v, w velocity after dtw
[ut1, vt1, wt1] =...
    TaylorDt3D(u0,            v0,            w0,            dx,dy,dz);
[ut2, vt2, wt2] =...
    TaylorDt3D(u0+0.5*ut1*dt, v0+0.5*vt1*dt, w0+0.5*wt1*dt, dx,dy,dz);
[ut3, vt3, wt3] =...
    TaylorDt3D(u0+0.5*ut2*dt, v0+0.5*vt2*dt, w0+0.5*wt2*dt, dx,dy,dz);
[ut4, vt4, wt4] =...
    TaylorDt3D(u0+    ut3*dt, v0+    vt3*dt, w0+    wt3*dt, dx,dy,dz);
u = u0 + (ut1+2*ut2+2*ut3+ut4)*dt/6;
v = v0 + (vt1+2*vt2+2*vt3+vt4)*dt/6;
w = w0 + (wt1+2*wt2+2*wt3+wt4)*dt/6;
end

function [u, v, w] = DivFilter3D(u0, v0, w0, Rate0, NLoop)
% 3D divergence-free filter using gradient descend optimization
% minimizing ||u0 - u||^2 under the constrain div u = 0
% Rate0: initial rate (much smaller than std(u0[v0,w0]))
% Nloop: number of loop
u = u0; v = v0; w = w0; d = divergence(u, v, w);
lambda = Rate0*d;
for iLoop = 1:NLoop
    rate = Rate0*(0.9/sqrt(iLoop/5+0.8)+0.1);
    tmp = 2*(u - u0) + circshift(lambda, 1,2) - circshift(lambda,-1,2);
    u(2:end-1,2:end-1,2:end-1)...
        = u(2:end-1,2:end-1,2:end-1) - rate*tmp(2:end-1,2:end-1,2:end-1);
    tmp = 2*(v - v0) + circshift(lambda, 1,1) - circshift(lambda,-1,1);
    v(2:end-1,2:end-1,2:end-1)...
        = v(2:end-1,2:end-1,2:end-1) - rate*tmp(2:end-1,2:end-1,2:end-1);
    tmp = 2*(w - w0) + circshift(lambda, 1,3) - circshift(lambda,-1,3);
    w(2:end-1,2:end-1,2:end-1)...
        = w(2:end-1,2:end-1,2:end-1) - rate*tmp(2:end-1,2:end-1,2:end-1);
    d = divergence(u, v, w);
    lambda = lambda + rate*d;
end
end

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

function P = PIterSolver3(u,v,w,ut,vt,wt,Mask,Nu,Rho,dx,dy,dz,Pinit,dev)
% Iterative Method for Pressure Estimation
% Pressure solver for time-resolved 3D3C velocity field
% need function LaplaceOperation3
% data arranged in MATLAB default meshgrid(Y, X, Z)
% u, v, w - velocity in current frame
% ut, vt, wt - time derivative of u, v, w
% Mask - 1 for activated nodes and 0 for others
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
mask = zeros(size(Mask) + 2);
mask(2:end-1,2:end-1,2:end-1) = Mask;
mask1 = circshift(mask, 1, 2).*mask;
mask2 = circshift(mask,-1, 2).*mask;
mask3 = circshift(mask, 1, 1).*mask;
mask4 = circshift(mask,-1, 1).*mask;
mask5 = circshift(mask, 1, 3).*mask;
mask6 = circshift(mask,-1, 3).*mask;
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
    PD = (+dx*px1 - (P - circshift(P, 1, 2))).*mask1 +...
         (-dx*px2 - (P - circshift(P,-1, 2))).*mask2 +...
         (+dy*py1 - (P - circshift(P, 1, 1))).*mask3 +...
         (-dy*py2 - (P - circshift(P,-1, 1))).*mask4 +...
         (+dz*pz1 - (P - circshift(P, 1, 3))).*mask5 +...
         (-dz*pz2 - (P - circshift(P,-1, 3))).*mask6;
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
