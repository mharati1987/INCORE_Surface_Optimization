clc
clear all
close all
data = importdata( 'frag400.txt' );
xraw      = data(:,1);
yraw  = data(:,2);
praw     = data(:,5);

%this is the first damage state

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

IC50A =   f.IC50A;
IC50B =   f.IC50B;
alpha =   f.alpha;
%alpha=5
n =       f.n ;
%n=6
%..........................................

SA=0:0.1:4; % set surge height range for fitted fragily
Moment=0:50:2500;  % set wave height range for fitted fragily


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
figure(1)

fragility=[X' Y' Z2'];
Z=Z2;
% triangulate and plot
tri = delaunay(X, Y);
trisurf(tri, X, Y, Z,'FaceColor','r');
% optional, could help make the plot look nicer
%shading interp
 xlabel('Sa (T1)[g]')
% xh = get(gca,'XLabel'); % Handle of the x label
% set(xh, 'Units', 'Normalized')
% pos = get(xh, 'Position');
% set(xh, 'Position',pos.*[1,1,1],'Rotation',12)
% 
 ylabel('Momentum Flux(m^3/s)')
 zlabel('Cumulative Probability')

% 
% yh = get(gca,'YLabel'); % Handle of the y label
% set(yh, 'Units', 'Normalized')
% pos = get(yh, 'Position');
% set(yh, 'Position',pos.*[1,1,1],'Rotation',-30)

hold on




%this is the next damage state (ds2)    %%%%%%%%%%%%%%%%%%%%%%%%%%%


xraw      = data(:,1);
yraw  = data(:,2);
praw     = data(:,4);

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

% plot( f, [xraw, yraw], praw );

IC50A =   f.IC50A;
IC50B =   f.IC50B;
alpha =   f.alpha;
%alpha=5
n =       f.n ;
%n=6
%..........................................
% 
% SA=0:0.07:6; % set surge height range for fitted fragily
% Moment=0:10:600;  % set wave height range for fitted fragily


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

ax1=figure(1)
fragility=[X' Y' Z2'];
Z=Z2;
% triangulate and plot
tri = delaunay(X, Y);
trisurf(tri, X, Y, Z,'FaceColor','g');
%colormap(ax1,pink)
hold on 



%this is the first damage state (ds3)   %%%%%%%%%%%%%%%%%%%%%%%

xraw      = data(:,1);
yraw  = data(:,2);
praw     = data(:,3);

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

% plot( f, [xraw, yraw], praw );

IC50A =   f.IC50A;
IC50B =   f.IC50B;
alpha =   f.alpha;
%alpha=5
n =       f.n ;
%n=6
%..........................................

% SA=0:0.07:6; % set surge height range for fitted fragily
% Moment=0:10:600;  % set wave height range for fitted fragily


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

fragility=[X' Y' Z2'];
Z=Z2;
% triangulate and plot

%ax2=figure(2)

tri = delaunay(X, Y);
trisurf(tri, X, Y, Z,'FaceColor','b');
%colormap(Z2,summer)

%colormap winter
hold on 
zlim([0 1.03])
ylim([0 1500])
xlim([0 4])
zlabel('Cumulative Probability')
