PathField = 'V_NOISE.mat';
OutputFile = 'V_SG.mat';
tic

% loading
load(PathField);
[~,M_SG,~] = MatSG(5,1,1);     % Local Polynomial (Savitzky-Golay) filter
[X, Y] = meshgrid(xb, yb);
U_3D_array = reshape(U_noisy, length(yb), length(xb), length(AFrame));
V_3D_array = reshape(V_noisy, length(yb), length(xb), length(AFrame));

% smoothing
disp('smoothing fields using Savitzky-Golay filter...');
weight = convn(ones(size(U_3D_array)), M_SG, 'same');
U_3D_array = convn(U_3D_array, M_SG, 'same')./weight;
V_3D_array = convn(V_3D_array, M_SG, 'same')./weight;
U = reshape(U_3D_array, [numel(X), length(AFrame)]);
V = reshape(V_3D_array, [numel(X), length(AFrame)]);
toc

% saving
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
function [MM,M,Mdt] = MatSG(w,dx,dt)
% w: window size; dx: spatial increment in grid; dt: temporal ...

polynomial =@(i,j,k) [1,...
                      i*dx, j*dx, k*dt,...
                      i*j*dx^2, i*k*dx*dt, j*k*dx*dt,...
                      i*i*dx^2, j*j*dx^2, k*k*dt^2];

w2=floor(w/2);
range = -w2:w2;
M = zeros(numel(range)^3,numel(polynomial(1,1,1)));

cont=0;
for k=-range
    for j=range
        for i=range
            cont=cont+1;
            M(cont,:) = polynomial(i,j,k);
        end
    end
end

MM  = (M'*M)\(M');%inv(M'*M)*M';

convcoef =@(col,w) reshape(col,[w w w]);

M    = convcoef(MM(1,:),w);
Mdt  = convcoef(MM(4,:),w);
Mddx = convcoef(MM(9,:),w);
Mddy = convcoef(MM(8,:),w);
end