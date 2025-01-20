SourceFile = 'V_NOISE.mat';
OutputFile = 'V_SG.mat';
AFrame = 1:120;
% grid
dt = 0.0065;   dx = 0.0114;
xb = 0:dx:1;   yb = 0:dx:1;   zb = 0:dx:0.5;   zb = zb(1:12);
[X,Y,Z] = meshgrid(xb,yb,zb);
% Local Polynomial (Savitzky-Golay) filter settings
[~, M_SG, ~] = MatSG4D(5, 1, 1);
tic

load(SourceFile, 'U_array','V_array','W_array');
disp('smoothing fields using Savitzky-Golay filter...');
U_4D_array = reshape(U_array, [size(X), size(U_array,2)]);
V_4D_array = reshape(V_array, [size(X), size(V_array,2)]);
W_4D_array = reshape(W_array, [size(X), size(W_array,2)]);
weight = convn(ones(size(U_4D_array)), M_SG, 'same');
U_4D_array = convn(U_4D_array, M_SG, 'same')./weight;
V_4D_array = convn(V_4D_array, M_SG, 'same')./weight;
W_4D_array = convn(W_4D_array, M_SG, 'same')./weight;
U = reshape(U_4D_array(:,:,:,AFrame), [numel(X), length(AFrame)]);
V = reshape(V_4D_array(:,:,:,AFrame), [numel(X), length(AFrame)]);
W = reshape(W_4D_array(:,:,:,AFrame), [numel(X), length(AFrame)]);
toc

% save
save(OutputFile, 'xb','yb','zb','dt', 'U','V','W');

%% functions
function [MM, M, Mdt] = MatSG4D(w, dx, dt)
% w: window size; dx: spatial increment in grid; dt: temporal ...
polynomial =@(i,j,k,l) [1, i*dx, j*dx, k*dx, l*dt,...
    i*i*dx*dx, i*j*dx*dx, i*k*dx*dx, i*l*dx*dt, j*j*dx*dx, j*k*dx*dx,...
    j*l*dx*dt, k*k*dx*dx, k*l*dx*dt, l*l*dt*dt];
w2=floor(w/2);
range = -w2:w2;
M = zeros(numel(range)^4,numel(polynomial(1,1,1,1)));
cont=0;
for l=range
    for k=range
        for j=range
            for i=range
                cont=cont+1;
                M(cont,:) = polynomial(i,j,k,l);
            end
        end
    end
end
MM  = (M'*M)\(M');%inv(M'*M)*M';
convcoef =@(col,w) reshape(col,[w w w w]);
M    = convcoef(MM(1,:),w);
Mdt  = convcoef(MM(5,:),w);
end