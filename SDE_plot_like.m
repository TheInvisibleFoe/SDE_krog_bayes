function SDE_plot_like()

output = importdata('~/Dropbox/bayesspt/SDE/code/samples/noisetrack_output.mat');
%output = importdata('~/Dropbox/bayesspt/SDE/code/samples/ftrack_output.mat');
for i = 1:4
    theta(i) = output(3).param_mean(i);
end
theta
l_size = 16;
%track = importdata('~/Dropbox/bayesspt/SDE/code/samples/ftrack.mat');
track = importdata('~/Dropbox/bayesspt/SDE/code/samples/noisetrack.mat');

obs = track.obs;
clean = track.clean;

%obs = track;
%clean =obs;

figure(1)

tp = subplot('Position',[0.1 0.2 0.35 0.55]);
plot(clean,'-')
hold on
plot(obs,'--o','MarkerSize',3)
hold off
xlabel('$t$','Interpreter','latex','FontSize',l_size)
ylabel('$x$','Interpreter','latex','Rotation',0,'Position',[-17 4],'FontSize',l_size)
tp.FontSize = l_size - 4;
grid on

D(1:1000) = 0.001/2 * (1:1000);
alpha(1:2001) = 0.002 * (-1000:1000);
f(1:2001) = -0.05 * 2/100 * (-1000:1000);

locerr(1:100) = max(theta(4),0.3) * (2/100) * (1:100);
for i = 1:length(f) 
   ll(i) = SDE_logl(obs,[theta(1) theta(2) f(i) theta(4)]);
end
like = exp(ll-max(ll)); 
tl = exp(SDE_logl(obs,[theta(1) theta(2) 0 0])-max(ll));
ax2 = subplot('Position',[0.55 0.17 0.15 0.25]);
plot(f,like)
line([-0. -0.],[0 1],'Color',[0 0.6 0],'Linewidth',1,'LineStyle','--'); 
ax2.XTick = [-1 -0.5 0 0.5 1];
grid on
xlabel('$f$','Interpreter','latex','FontSize',l_size)
ylabel('$P(data|\theta)$','Interpreter','latex','Rotation',90,'FontSize',l_size)
text(-0.95,2.6,'Normalized likelihood','Interpreter','latex','FontSize',l_size+2)
text(-6.2,2.6,'Simulated data','Interpreter','latex','FontSize',l_size+2)
ax2.YTick = [];
for i = 1:length(alpha) 
   ll(i) = SDE_logl(obs,[theta(1) alpha(i) theta(3) theta(4)]);
end
like = exp(ll-max(ll)); 
ax3 = subplot('Position',[0.75 0.52 0.15 0.25]);
plot(alpha,like)
line([1 1],[0 1],'Color',[0 0.6 0],'Linewidth',1,'LineStyle','--');    
ax3.XTick = [-2 -1 0 1 2];
grid on
xlabel('$\alpha$','Interpreter','latex','FontSize',l_size)
ylabel('$P(data|\theta)$','Interpreter','latex','Rotation',90,'FontSize',l_size)
ax3.YTick = [];
for i = 1:length(D) 
   llD(i) = SDE_logl(obs,[D(i) theta(2) theta(3) theta(4)]);
end
likeD = exp(llD-max(llD)); 
ax4 = subplot('Position',[0.55 0.52 0.15 0.25]);
plot(D,likeD)
line([.2 .2],[0 1],'Color',[0 0.6 0],'Linewidth',1,'LineStyle','--'); 
ax4.XTick = [0 0.2 0.4];
grid on
xlabel('$D_0$','Interpreter','latex','FontSize',l_size)
ylabel('$P(data|\theta)$','Interpreter','latex','Rotation',90,'FontSize',l_size)
ax4.YTick = [];
for i = 1:length(locerr) 
   lle(i) = SDE_logl(obs,[theta(1) theta(2) theta(3) locerr(i)]);
end
likee = exp(lle-max(lle)); 
ax5 = subplot('Position',[0.75 0.17 0.15 0.25]);
plot(locerr,likee)
%line([40 40],[0 1],'Color',[0 0.6 0],'Linewidth',1,'LineStyle','--'); 
line([2 2],[0 1],'Color',[0 0.6 0],'Linewidth',1,'LineStyle','--'); 
grid on
%ax5.XTick = [0 30 60 90];
ax5.XTick = [0 2 4 6];
xlabel('$\sigma_{mn}^2$','Interpreter','latex','FontSize',l_size)
ylabel('$P(data|\theta)$','Interpreter','latex','Rotation',90,'FontSize',l_size)
axis([0 inf 0 1])
ax5.YTick = [];
ax1.FontSize=l_size-4;
ax2.FontSize=l_size-4;
ax3.FontSize=l_size-4;
ax4.FontSize=l_size-4;
ax5.FontSize=l_size-4;

fig=gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 14/2 10/2];
print('~samples/samplefig_noise','-depsc')
