function SDE_plot_inf()

data = importdata('samples/ftrack_data_output.mat');
track_max = 100;
n_max = 7; %data.n_max;
savepath = 'samples/inf_free';
l_size = 16;

% Find best model (choose one)
for i = 1:4
    %evi(i) = data(i).Z_norm;
end
%[~,best] = max(evi);
best = 1;

MM = [0 0; 1 0; 0 1; 1 1];
SDE_params(data.results(best).maxLpar,MM(best,:))
fprintf('\nInformation for model %i: %2.3f \n',best,data.results(best).H(1))

H_star = zeros(track_max,n_max);
H_ref = zeros(track_max,n_max);
frac = zeros(1,n_max);
prob = zeros(1,n_max);
post = [data.results(best).samples.post];
post_cum = cumsum(post);

for n=1:n_max

   scaled_obs = SDE_scaling(data.data,n);
    for i=1:track_max
        % Use rand to draw random sample and generate a track 
        point = rand;
        draw = find(post_cum > point,1);
        theta = SDE_params(data.results(best).samples(draw).theta,MM(best,:));
        new_obs = SDE_replicate(scaled_obs,theta,n); % 
        
        % Calculate and store information
        logl = SDE_logl_m(new_obs,theta,n); % 
        H_star(i,n) = logl - data.results(best).logZ(n); %
        logl_ref = SDE_logl_m(scaled_obs,theta,n);
        H_ref(i,n) = logl_ref - data.results(best).logZ(n);
        
    end
    
    frac(n) = mean(H_star(:,n) > data.results(best).H(n));
    prob(n) = mean(H_star(:,n) > H_ref(:,n));
end

for i = 1:n_max;
    fprintf('\nInformation for model %i (scaled by %i)     : %2.3f \n',best,i,data.results(best).H(i))
    fprintf('Fraction of tracks with H* < H (scaled by %i) :  %1.2f  \n',i,frac(i))
    fprintf('Probability of H* < H_ref (scaled by %i)      :   %.2f \n',i,prob(i))
end
% Plot information along with original

figure(3);
inp = scatter(H_star(:,1),ones(1,track_max),'d');
hold on
scatter(H_ref(:,1),ones(1,track_max) + 0.1,'Black','x')
txt = sprintf('$p(H^* > H) =  %.2f$',prob(1));
text(80,1,txt,'Interpreter','Latex','FontSize',l_size)
title('Information per step for scaled tracks','Interpreter','Latex','FontSize',l_size)
ylabel('$n$','Interpreter','Latex','FontSize',l_size,'Rotation',0)
DH = -130;
dh = 1.4;
n_max = n_max - 2;
scatter(DH-1,n_max + 2,700,'Black','x','Linewidth',3)
text(DH + 4*dh,n_max + 2,'Scaled data ','Interpreter','Latex','FontSize',l_size)
scatter(DH-1.0 - dh,n_max + 1,'d')
scatter(DH-1.0 - 2*dh,n_max + 1.2,'d')
scatter(DH-1 + dh,n_max + 1.4,'d')
scatter(DH-1 + dh,n_max + 1.,'d')
scatter(DH-1 + 2*dh,n_max + 1.2,'d')
scatter(DH-1 - dh,n_max + 1.4,'d')
scatter(DH-1,n_max + 1.2,'d')
text(DH + 4*dh,n_max + 1.2,'Simulated tracks ','Interpreter','Latex','FontSize',l_size)
text(DH + 13*dh,n_max + .7,'($f = 0$)','Interpreter','Latex','FontSize',l_size)
scatter(data.results(best).H(1) - 0.1,1,700,'Black','x','Linewidth',3)
for n = 2:n_max + 2 
    scatter(H_star(:,n),n * ones(1,track_max) - 0.1 ,'d')
    scatter(H_ref(:,n),n * ones(1,track_max) + 0.1,'Black','x')
    scatter(data.results(best).H(n),n,700,'Black','x','Linewidth',3)
    txt = sprintf('$p(H^* > H) =  %.2f$',prob(n));
text(80,n,txt,'Interpreter','Latex','FontSize',l_size)
end
ylim([0 7.9])
xlim([-150 155])
xt = get(gca, 'XTick');
set(gca, 'FontSize', l_size)
hold off
print(savepath,'-depsc');

