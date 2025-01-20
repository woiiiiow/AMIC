% Woii 240619
% for additional wing foil LES

PathField = 'V_NOISE.mat';
PathParameter = '../../../airfoilLES/airfoilLES_parameters.h5';
OutputFile = 'V_AMIC.mat';
NLoop = 11;                    % Total loops of AMIC
Lambda = [0.15 0.1];           % iteration coefficients
tic

disp('loading data...');
load(PathField);
[X, Y] = meshgrid(xb, yb);
dt = h5read(PathParameter, '/dt');
toc

% smoothing
disp('smoothing fields using AMIC...');
U_array = U_noisy;             V_array = V_noisy;
Enabled = ones(size(U_array));
Vrms = rms(U_array,'all') + rms(V_array,'all');
A_Vdis = zeros(1,NLoop); % L-2 distance to noisy field, velocity
for iLoop = 1:NLoop
    % turbulence frozen filter
    % sequence frame P+1 P-1 P+2 P-2 P+3 P-3...
    U_cor = zeros([size(U_array) length(Lambda)*2]);
    V_cor = U_cor;
    for iFrame = 1:length(AFrame)
        u0 = reshape(U_array(:,iFrame), size(X));
        v0 = reshape(V_array(:,iFrame), size(X));
        mask0 = reshape(Enabled(:,iFrame), size(X));
        u = u0;   v = v0;   mask = mask0;
        for iStep = 1:length(Lambda)
            [u, v] = MRK4Prop2D(u, v, X, Y, dt);
            mask = interp2(X,Y, mask, X-u*dt,Y-v*dt);
            mask(isnan(mask)) = 0;
            mask = imgaussfilt(imerode(mask,ones(3)), 5);
            iSub = iFrame + iStep;
            if iSub < length(AFrame)
                U_cor(:,iSub,iStep*2-1) = (u(:)-U_array(:,iSub)).*mask(:);
                V_cor(:,iSub,iStep*2-1) = (v(:)-V_array(:,iSub)).*mask(:);
            end
        end
        u = u0;   v = v0;   mask = mask0;
        for iStep = 1:length(Lambda)
            [u, v] = MRK4Prop2D(u, v, X, Y,-dt);
            mask = interp2(X,Y, mask, X+u*dt,Y+v*dt);
            mask(isnan(mask)) = 0;
            mask = imgaussfilt(imerode(mask,ones(3)), 5);
            iSub = iFrame - iStep;
            if iSub > 0
                U_cor(:,iSub,iStep*2) = (u(:)-U_array(:,iSub)).*mask(:);
                V_cor(:,iSub,iStep*2) = (v(:)-V_array(:,iSub)).*mask(:);
            end
        end
    end
    % ratio of update
    index_cth = (rms(U_cor(U_cor~=0))+rms(V_cor(V_cor~=0)))/Vrms;
    % correcting
    for iStep = 1:length(Lambda)
        U_array = U_array...
            + Lambda(iStep)*(U_cor(:,:,iStep*2-1)+U_cor(:,:,iStep*2));
        V_array = V_array...
            + Lambda(iStep)*(V_cor(:,:,iStep*2-1)+V_cor(:,:,iStep*2));
    end

    A_Vdis(iLoop) = std([U_array-U_noisy;V_array-V_noisy],1, 'all');
    fprintf(['iLoop = ',num2str(iLoop),...
        '\tcorrection Taylor = ',num2str(index_cth,'%.5g'),'\n']);
    if iLoop > 2
        fprintf(['\b\tcorrection index = ',...
            num2str((A_Vdis(iLoop)-A_Vdis(iLoop-1))/...
            (A_Vdis(iLoop-1)-A_Vdis(iLoop-2)),'%.5g'),'\n']);
    end
    toc
end

% save
U = U_array;                   V = V_array;
save(OutputFile, 'U','V','AFrame');

%% display
figure; iFrame = 24;
subplot(3,2,1); pcolor(X, Y, reshape(U_clean(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([0.6 1.1]);
subplot(3,2,3); pcolor(X, Y, reshape(U_noisy(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([0.6 1.1]);
subplot(3,2,5); pcolor(X, Y, reshape(U      (:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([0.6 1.1]);
subplot(3,2,2); pcolor(X, Y, reshape(V_clean(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([-0.25 0.25]);
subplot(3,2,4); pcolor(X, Y, reshape(V_noisy(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([-0.25 0.25]);
subplot(3,2,6); pcolor(X, Y, reshape(V      (:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([-0.25 0.25]);

%% functions
function [ut, vt] = TaylorDt2D(u, v, dx, dy)
% ut = -(u_mean\cdot\nabla)u_fluc
u_mean = imgaussfilt(u, 7);
v_mean = imgaussfilt(v, 7);
u_fluc = u - u_mean;
v_fluc = v - v_mean;
[ux, uy] = gradient(u_fluc); ux = ux./dx; uy = uy./dy;
[vx, vy] = gradient(v_fluc); vx = vx./dx; vy = vy./dy;
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
% mask = interp3(x, y, z, ones(size(x)), x-u0*dt, y-v0*dt, z-w0*dt);
[ut1, vt1] = TaylorDt2D(u0,            v0,            dx, dy);
[ut2, vt2] = TaylorDt2D(u0+0.5*ut1*dt, v0+0.5*vt1*dt, dx, dy);
[ut3, vt3] = TaylorDt2D(u0+0.5*ut2*dt, v0+0.5*vt2*dt, dx, dy);
[ut4, vt4] = TaylorDt2D(u0+    ut3*dt, v0+    vt3*dt, dx, dy);
u = u0 + (ut1+2*ut2+2*ut3+ut4)*dt/6;
v = v0 + (vt1+2*vt2+2*vt3+vt4)*dt/6;
end