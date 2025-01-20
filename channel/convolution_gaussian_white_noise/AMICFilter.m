SourceFile = 'V_NOISE.mat';
OutputFile = 'V_AMIC.mat';
AFrame = 1:120;
% grid
dt = 0.0065;   dx = 0.0114;
xb = 0:dx:1;   yb = 0:dx:1;   zb = 0:dx:0.5;   zb = zb(1:12);
[X,Y,Z] = meshgrid(xb,yb,zb);
% smoothing parameters
NLoop = 8;
Lambda = [0.15 0.15];
SLoopDiv = 1;
RateDiv = 0.4;
NItDiv = 0;
tic

% disp('loading fields...')
load(SourceFile, 'U_array','V_array','W_array','Enabled');
U_noise = U_array;  V_noise = V_array;  W_noise = W_array;

% smoothing
disp('smoothing fields using physics-based multifold iterative method...');
Vrms = rms(U_array,'all')+rms(V_array,'all')+rms(W_array,'all');
A_Vdis = zeros(1,NLoop); % L-2 distance to noisy field, velocity
for iLoop = 1:NLoop
    % turbulence frozen filter
    % sequence frame P+1 P-1 P+2 P-2 P+3 P-3...
    U_cor = zeros([size(U_array) length(Lambda)*2]);
    V_cor = U_cor;             W_cor = U_cor;
    for iFrame = 1:length(AFrame)
        u0 = reshape(U_array(:,iFrame), size(X));
        v0 = reshape(V_array(:,iFrame), size(X));
        w0 = reshape(W_array(:,iFrame), size(X));
        mask0 = reshape(Enabled(:,iFrame), size(X));
        u = u0;   v = v0;   w = w0;   mask = mask0;
        for iStep = 1:length(Lambda)
            [u, v, w] = MRK4Prop3D(u, v, w, X, Y, Z, dt);
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
            [u, v, w] = MRK4Prop3D(u, v, w, X, Y, Z,-dt);
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
    index_cth = (rms(U_cor(U_cor~=0))+rms(V_cor(V_cor~=0))...
        +rms(W_cor(W_cor~=0)))/Vrms;
    % correcting
    for iStep = 1:length(Lambda)
        U_array = U_array...
            + Lambda(iStep)*(U_cor(:,:,iStep*2-1)+U_cor(:,:,iStep*2));
        V_array = V_array...
            + Lambda(iStep)*(V_cor(:,:,iStep*2-1)+V_cor(:,:,iStep*2));
        W_array = W_array...
            + Lambda(iStep)*(W_cor(:,:,iStep*2-1)+W_cor(:,:,iStep*2));
    end

    % divergence free filter
    if iLoop >= SLoopDiv
        div = 0;
        for iFrame = 1:length(AFrame)
            u0 = reshape(U_array(:,iFrame), size(X));
            v0 = reshape(V_array(:,iFrame), size(X));
            w0 = reshape(W_array(:,iFrame), size(X));
            [u, v, w] = DivFilter3D(u0, v0, w0, RateDiv, NItDiv);
            U_array(:,iFrame) = u(:);
            V_array(:,iFrame) = v(:);
            W_array(:,iFrame) = w(:);
            div = div+rms(u-u0,'all')+rms(v-v0,'all')+rms(w-w0,'all');
        end
        index_cdiv = div/3/length(AFrame)/Vrms;
    end

    A_Vdis(iLoop) =...
        std([U_array-U_noise;V_array-V_noise;W_array-W_noise],1, 'all');

    fprintf(['iLoop = ',num2str(iLoop),...
        '\tcorrection advection = ',num2str(index_cth,'%.5g'),...
        '\tcorrection div = ',num2str(index_cdiv,'%.5g'),'\n']);
    if iLoop > 2
        fprintf(['\b\tcorrection index = ',...
            num2str((A_Vdis(iLoop)-A_Vdis(iLoop-1))/...
            (A_Vdis(iLoop-1)-A_Vdis(iLoop-2)),'%.5g'),'\n']);
    end
    toc
end

% save
U = U_array(:,AFrame);   V = V_array(:,AFrame);   W = W_array(:,AFrame);
save(OutputFile, 'xb','yb','zb','dt', 'U','V','W');

%% functions
function [ut, vt, wt] = TaylorDt3D(u, v, w, Increx, Increy, Increz)
% ut = -(u_mean\cdot\nabla)u_fluc
u_mean = imgaussfilt3(u, 7);
v_mean = imgaussfilt3(v, 7);
w_mean = imgaussfilt3(w, 7);

u_profile = mean(u, [1,3]);
for i = 73:88
    u_mean(:,i,:) = u_mean(:,i,:)...
        - mean(u_mean(:,i,:),'all') + u_profile(i);
end

u_fluc = u - u_mean;
v_fluc = v - v_mean;
w_fluc = w - w_mean;
[ux, uy, uz] = gradient(u_fluc);
[vx, vy, vz] = gradient(v_fluc);
[wx, wy, wz] = gradient(w_fluc);
ux = ux./Increx; uy = uy./Increy; uz = uz./Increz;
vx = vx./Increx; vy = vy./Increy; vz = vz./Increz;
wx = wx./Increx; wy = wy./Increy; wz = wz./Increz;
ut = - u_mean.*ux - v_mean.*uy - w_mean.*uz;
vt = - u_mean.*vx - v_mean.*vy - w_mean.*vz;
wt = - u_mean.*wx - v_mean.*wy - w_mean.*wz;

ut = imgaussfilt3(ut, 0.75);
vt = imgaussfilt3(vt, 0.75);
wt = imgaussfilt3(wt, 0.75);
end

function [u, v, w] = MRK4Prop3D(u0, v0, w0, x, y, z, dt)
% Runge-Kutta propogation of 3D flow fields
% u0, v0, w0, current velocity
% x, y, z corrodinates
% dt time increments
% u, v, w velocity after dtw
dx = x(1,2,1) - x(1); dy = y(2,1,1) - y(1); dz = z(1,1,2) - z(1);
% mask = interp3(x, y, z, ones(size(x)), x-u0*dt, y-v0*dt, z-w0*dt);
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

