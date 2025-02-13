clc
clear all
close all
data = importdata( 'frag900.txt' );
xraw      = data(:,1);
yraw  = data(:,2);
praw     = data(:,5);

ft = fittype( 'Emax*( CA/IC50A + CB/IC50B + alpha*( CA/IC50A ) * ( CB/IC50B ) )^n /(( CA/IC50A + CB/IC50B + alpha*( CA/IC50A ) * ( CB/IC50B ) )^n  + 1 )', ...
    'independent', {'CA', 'CB'}, 'dependent', 'z', 'problem', 'Emax' )

Emax = 1.01;

opts = fitoptions( ft );

%opts.Lower =[0, 0,0, -0]
opts.Lower = [0, 0, 0, -0];
%opts.Lower = [-20, 5, 3, -10];
%opts.Lower =[-20, -20, -4.1, 1]
opts.Robust = 'LAR';
%opts.StartPoint = [0.0089, 0.706, 1.0, 0.746];
opts.StartPoint = [0.0089, 0.706, 1.0, 0.746];


[f, gof] = fit( [xraw, yraw], praw, ft,...
    opts, 'problem', Emax )

plot( f, [xraw, yraw], praw );
zlim([0 1.03])

%shading interp
 xlabel('Sa(T1=0.94 sec) [g]')
% xh = get(gca,'XLabel'); % Handle of the x label
% set(xh, 'Units', 'Normalized')
% pos = get(xh, 'Position');
% set(xh, 'Position',pos.*[1,1,1],'Rotation',12)
% 
 ylabel('Momentum Flux(m^3/s)')
 zlabel('Failure Probability')
IC50A =   f.IC50A;
IC50B =   f.IC50B;
alpha =   f.alpha;
%alpha=5
n =       f.n ;
%n=6
%..........................................

SA=0:0.07:4; % set surge height range for fitted fragily
Moment=0:50:1600;  % set wave height range for fitted fragily


k=1;
for i=1:length(SA)
    for j=1:length(Moment)        
            CA=SA(i);
            CB=Moment(j);
            X(k)=SA(i);
            Y(k)=Moment(j);
            Z2(k) = Emax*( CA/IC50A + CB/IC50B + alpha*( CA/IC50A ) * ( CB/IC50B ...
                                ) )^n /(( CA/IC50A + CB/IC50B + alpha*( CA/IC50A ) ...
                                * ( CB/IC50B ) )^n  + 1 );
        k=k+1;    
end
end
figure(2)

fragility=[X' Y' Z2'];
Z=Z2;
% triangulate and plot
tri = delaunay(X, Y);
trisurf(tri, X, Y, Z);
% optional, could help make the plot look nicer
%shading interp
 xlabel('Sa(T1=0.94 sec) [g]')
% xh = get(gca,'XLabel'); % Handle of the x label
% set(xh, 'Units', 'Normalized')
% pos = get(xh, 'Position');
% set(xh, 'Position',pos.*[1,1,1],'Rotation',12)
% 
 ylabel('Moment Flux(m^3/s)')
 zlabel('Cumulative Probability')
% 
% yh = get(gca,'YLabel'); % Handle of the y label
% set(yh, 'Units', 'Normalized')
% pos = get(yh, 'Position');
% set(yh, 'Position',pos.*[1,1,1],'Rotation',-30)
zlim([0 1.03])
ylim([0 1500])
zlim([0 1.03])


%new one
figure(3)

plot3( xraw , yraw , praw,'o','Color','w','MarkerSize',6.5,...
    'MarkerFaceColor','#0000FF')   %'#0000FF'
ylabel('Moment Flux(m^3/s)')
zlabel('Cumulative Probability')
xlabel('Sa(T1=0.94 sec) [g]')