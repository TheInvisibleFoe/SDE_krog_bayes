function SDE_plot_inf()

data = importdata('samples/ftrack_output.mat');
track_max = 100;
nlist = data.options.nlist;
savepath = '../draft/inf_force';
l_size = 16;

% Find best model (choose one)
for i = 1:4
    %evi(i) = data(i).Z_norm;
end
%[~,best] = max(evi);
best = 4;

MM = [0 0; 1 0; 0 1; 1 1];
fprintf('\nInformation for model %i: %2.3f \n',best,data.results(best).H(1))

N = length(nlist);
H_star = zeros(track_max,N);
H_ref = zeros(track_max,N);
prob = zeros(1,N);
post = [data.results(best).samples.post];
post_cum = cumsum(post);
H_check = zeros(1,length(nlist));


for n=1:N
    n2=nlist(n);
    scaled_obs = SDE_scaling(data.data,n2);
    H_check(n) = -data.results(best).logZ(n);
    for j = 1:length(post)
        tht = SDE_params(data.results(best).samples(j).theta,MM(best,:));
        H_check(n) = H_check(n) + SDE_logl_m(scaled_obs,tht,n2)*data.results(best).samples(j).post;
    end
        
   for i=1:track_max
        % Use rand to draw random sample and generate a track 
        point = rand;
        draw = find(post_cum > point,1);
        theta = SDE_params(data.results(best).samples(draw).theta,MM(best,:));
        new_obs = SDE_replicate(scaled_obs,theta,n2); % 
        
        % Calculate and store information
        logl = SDE_logl_m(new_obs,theta,n2); % 
        H_star(i,n) = logl - data.results(best).logZ(n); %
        logl_ref = SDE_logl_m(scaled_obs,theta,n2);
        H_ref(i,n) = logl_ref - data.results(best).logZ(n);
        
   end
   prob(n) = mean(H_star(:,n) > H_ref(:,n));
end

%for i = 1:N;
  %  n2 ) nlist(n);
  %  fprintf('\nInformation for model %i (scaled by %i)     : %2.3f \n',best,i,data.results(best).H(i))
  %  fprintf('Fraction of tracks with H* < H (scaled by %i) :  %1.2f  \n',i,frac(i))
  %  fprintf('Probability of H* < H_ref (scaled by %i)      :   %.2f \n',i,prob(i))
%end
% Plot information along with original

figure(3);
scatter(H_star(:,1),ones(1,track_max),'d');
hold on
scatter(H_ref(:,1),ones(1,track_max) + 0.3,'Black','x')
txt = sprintf('$p(h^* > h) =  %.2f$',prob(1));
text(80,1,txt,'Interpreter','Latex','FontSize',l_size)
title('Information per step for scaled tracks','Interpreter','Latex','FontSize',l_size)
ylabel('$n$','Interpreter','Latex','FontSize',l_size,'Rotation',0)
DH = -130;
dh = 1.4;
scatter(DH-1,N + 2,700,'Black','x','Linewidth',2)
text(DH + 4*dh,N + 2,'Mean for scaled data ','Interpreter','Latex','FontSize',l_size)
scatter(DH-1.0 - dh,N + 1,'d')
scatter(DH-1.0 - 2*dh,N + 1.2,'d')
scatter(DH-1 + dh,N + 1.4,'d')
scatter(DH-1 + dh,N + 1.,'d')
scatter(DH-1 + 2*dh,N + 1.2,'d')
scatter(DH-1 - dh,N + 1.4,'d')
scatter(DH-1,N + 1.2,'d')
text(DH + 4*dh,N + 1.2,'Simulated tracks ','Interpreter','Latex','FontSize',l_size)
text(DH + 13*dh,N + .7,'($\sigma_{mn} \neq 0$)','Interpreter','Latex','FontSize',l_size)
scatter(mean(H_ref(:,1)),1,700,'Black','x','LineWidth',2)
for n = 2:N;
    n2 = nlist(n); 
    scatter(H_star(:,n),n2 * ones(1,track_max) - 0. ,'d')
    scatter(H_ref(:,n),n2 * ones(1,track_max) + 0.3,'Black','x')
    scatter(mean(H_ref(:,n)),n2,700,'Black','x','LineWidth',2)
    txt = sprintf('$p(h^* > h) =  %.2f$',prob(n));
text(80,n2,txt,'Interpreter','Latex','FontSize',l_size)
end
ylim([0 N+3])
xlim([-150 155])
xt = get(gca, 'XTick');
set(gca, 'FontSize', l_size)
hold off
%print(savepath,'-depsc');

    

