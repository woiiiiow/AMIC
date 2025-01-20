% Woii, junwei.chen@uc3m.es, EAP, UC3M
% AMIC to the 2D PIV dataset of airfoil's wake

SourceFolder    = '../../Experiment_wing3/OUT_TRPIV/';
OutputFile      = 'V_AMIC.mat';
ROI = [41 110 41 150];          % subdomain of PIV field
                                % [first-last row, first-last column]
Mm_per_px_ratio = 0.1197;       % spatial resolution of photos
Sample_rate     = 30;           % sample rate of PIV
Vector_Spacing  = 10;           % vector per pixel in PIV

AFrame          = 11:530;       % frames to be processed
NLoop           = 100;          % Total loops of AMIC
Lambda          = [0.1 0.05];   % iteration coefficients
ThVCI           = 0.9;          % iteration stops when the velocity
                                % correction index exceeds the threshold
                                % default: 0.9
InitLoopDiv     = 100;          % first loop to start divergence-free filter
RateDiv         = 0.4;          % initial rate of divergence-free filter
NItDiv          = 1;            % loops of divergence-free filter everytime

load([SourceFolder, 'Grid_Wing.mat'], 'X','Y');
% entire domain coordinates
X_etr = X;                      Y_etr = Y;
dt = 1/Sample_rate;
tic

disp('loading data...');
U_array = zeros(numel(X_etr), length(AFrame));   V_array = U_array;
iCount = 0;
for iFrame = AFrame
    iCount = iCount + 1;
    load([SourceFolder,'Wing_',num2str(iFrame,'%06u'),'.mat'], 'U','V');
    U_array(:, iCount) = V(:);    V_array(:, iCount) = U(:);
    % swapping U and V
end
U_array = U_array.*Sample_rate.*Mm_per_px_ratio.*1e-3;
V_array = V_array.*Sample_rate.*Mm_per_px_ratio.*1e-3;
clear iFrame iCount U V
Enabled = ones(size(U_array));
toc

% smoothing
disp('smoothing fields using AMIC...');
Vrms = rms(U_array,'all') + rms(V_array,'all');
A_Vdis = zeros(1,NLoop); % L-2 distance of velocity to the noisy field
U_noise = U_array;  V_noise = V_array;
v_cor_idx = 0;                 % velocity correction index
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
            mask(isnan(mask)) = 0; mask = imgaussfilt(mask, 5);
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
            mask(isnan(mask)) = 0; mask = imgaussfilt(mask, 5);
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

    % divergence-free filter
    index_cdiv = 0;
    if iLoop >= InitLoopDiv
        div = 0;
        for iFrame = 1:length(AFrame)
            u0 = reshape(U_array(:,iFrame), size(X));
            v0 = reshape(V_array(:,iFrame), size(X));
            [u, v] = DivFilter2D(u0, v0, RateDiv, NItDiv);
            U_array(:,iFrame) = u(:);
            V_array(:,iFrame) = v(:);
            div = div+rms(u-u0,'all')+rms(v-v0,'all');
        end
        index_cdiv = div/2/length(AFrame)/Vrms;
    end

    fprintf(['iLoop = ',num2str(iLoop),...
        '\tcorrection Taylor = ',num2str(index_cth,'%.5g'),...
        '\tcorrection Div = ',num2str(index_cdiv,'%.5g'),'\n']);
    A_Vdis(iLoop) = rms([U_array-U_noise;V_array-V_noise],'all');
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

% save
U = U_array;                   V = V_array;
save(OutputFile, 'U','V','AFrame');

%% functions
function enabled = BMedian2D(u, FiltSize, E, Th)
u_med = medfilt2(u, [FiltSize, FiltSize]);
r     = abs(u - u_med);
r_med = medfilt2(r, [FiltSize, FiltSize]);
enabled = r./(r_med+E) < Th;
end

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
% mask = interp3(x, y, z, ones(size(x)), x-u0*dt, y-v0*dt, z-w0*dt);
[ut1, vt1] = TaylorDt2D(u0,            v0,            dx, dy);
[ut2, vt2] = TaylorDt2D(u0+0.5*ut1*dt, v0+0.5*vt1*dt, dx, dy);
[ut3, vt3] = TaylorDt2D(u0+0.5*ut2*dt, v0+0.5*vt2*dt, dx, dy);
[ut4, vt4] = TaylorDt2D(u0+    ut3*dt, v0+    vt3*dt, dx, dy);
u = u0 + (ut1+2*ut2+2*ut3+ut4)*dt/6;
v = v0 + (vt1+2*vt2+2*vt3+vt4)*dt/6;
end

function [u, v] = DivFilter2D(u0, v0, Rate0, NLoop)
% 3D divergence-free filter using gradient descend optimization
% minimizing ||u0 - u||^2 under the constrain div u = 0
% Rate0: initial rate (much smaller than std(u0) or std(v0))
% Nloop: number of loop
u = u0; v = v0; d = divergence(u, v);
lambda = Rate0*d;
for iLoop = 1:NLoop
    rate = Rate0*(0.9/sqrt(iLoop/5+0.8)+0.1);
    tmp = 2*(u - u0) + circshift(lambda, 1,2) - circshift(lambda,-1,2);
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) - rate*tmp(2:end-1,2:end-1);
    tmp = 2*(v - v0) + circshift(lambda, 1,1) - circshift(lambda,-1,1);
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) - rate*tmp(2:end-1,2:end-1);
    d = divergence(u, v);
    lambda = lambda + rate*d;
end
end
