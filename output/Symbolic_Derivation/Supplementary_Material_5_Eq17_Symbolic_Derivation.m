%% Supplementary Material 5: MATLAB symbolic derivation of Eq. (17)
% Rows of K:    sigma_11, sigma_22, sigma_12
% Columns of K: n = 1, 2, 3, 4
%
% The notation follows Appendix B of the revised manuscript:
%   Sigma_1, Sigma_2 : remote normal stresses
%   S, D             : mean and deviatoric remote stresses
%   u = a^2/r^2      : dimensionless radial coordinate
%   s_d              : Cartesian disturbance stress
%   J(i,n)           : angular integral
%   K(i,n)           : coefficient in Eqs. (16) and (17)
%
% Required product: Symbolic Math Toolbox

clear;
clc;

%% 1. Symbols and Kirsch solution
syms a r theta u Sigma_1 Sigma_2 real
assume(a > 0);
assume(r >= a);
assumeAlso(u > 0);
assumeAlso(u <= 1);

S = (Sigma_1 + Sigma_2)/2;
D = (Sigma_1 - Sigma_2)/2;
q = (a/r)^2;

% Polar stress components in Eq. (1).
s_rr = S*(1-q) + D*(1-4*q+3*q^2)*cos(2*theta);
s_tt = S*(1+q) - D*(1+3*q^2)*cos(2*theta);
s_rt = -D*(1+2*q-3*q^2)*sin(2*theta);

%% 2. Cartesian disturbance stress
% Cartesian stress components in Eq. (2).
s_H = [ ...
    (s_rr+s_tt)/2 + (s_rr-s_tt)/2*cos(2*theta) ...
        - s_rt*sin(2*theta); ...
    (s_rr+s_tt)/2 - (s_rr-s_tt)/2*cos(2*theta) ...
        + s_rt*sin(2*theta); ...
    (s_rr-s_tt)/2*sin(2*theta) + s_rt*cos(2*theta)];

% Disturbance stress in Eq. (3).
s_d = simplify(s_H - [Sigma_1; Sigma_2; 0]);

% Compact form used in Appendix B.
s_d_u = [ ...
    -(S+2*D)*u*cos(2*theta) + D*u*(3*u-2)*cos(4*theta); ...
     (S-2*D)*u*cos(2*theta) - D*u*(3*u-2)*cos(4*theta); ...
    -S*u*sin(2*theta)       + D*u*(3*u-2)*sin(4*theta)];

% Check that the compact form equals the transformed Kirsch field.
check_field = simplify(expand(s_d - subs(s_d_u,u,q)));
if ~isequal(check_field,sym(zeros(3,1)))
    error('The compact disturbance field does not match the Kirsch solution.');
end

%% 3. Angular and radial integrations
% From u=a^2/r^2 and A0=pi*a^2,
% K(i,n) = 1/(2*pi) * int_0^1 u^(-2) J(i,n) du,
% J(i,n) = int_0^(2*pi) [s_d_u(i)]^n dtheta.
n = 1:4;
J = sym(zeros(3,4));
K = sym(zeros(3,4));

for i = 1:3
    for j = 1:4
        J(i,j) = simplify(int(expand(s_d_u(i)^n(j)),theta,0,2*pi));
        K(i,j) = simplify(int(J(i,j)/u^2,u,0,1)/(2*pi));
    end
end
K = simplify(expand(K));

%% 4. Eq. (17) and exact symbolic check
K_eq17 = sym(zeros(3,4));

K_eq17(1,2) = 5*Sigma_1^2/4-Sigma_1*Sigma_2+Sigma_2^2/4;
K_eq17(2,2) = Sigma_1^2/4-Sigma_1*Sigma_2+5*Sigma_2^2/4;
K_eq17(3,2) = Sigma_1^2/4+Sigma_2^2/4;

K_eq17(1,4) = 335*Sigma_1^4/448-323*Sigma_1^3*Sigma_2/280 ...
    +797*Sigma_1^2*Sigma_2^2/1120-57*Sigma_1*Sigma_2^3/280 ...
    +51*Sigma_2^4/2240;
K_eq17(2,4) = 51*Sigma_1^4/2240-57*Sigma_1^3*Sigma_2/280 ...
    +797*Sigma_1^2*Sigma_2^2/1120-323*Sigma_1*Sigma_2^3/280 ...
    +335*Sigma_2^4/448;
K_eq17(3,4) = 51*Sigma_1^4/2240+3*Sigma_1^3*Sigma_2/140 ...
    +41*Sigma_1^2*Sigma_2^2/1120+3*Sigma_1*Sigma_2^3/140 ...
    +51*Sigma_2^4/2240;

check_eq17 = simplify(expand(K-K_eq17));
if ~isequal(check_eq17,sym(zeros(3,4)))
    error('The derived matrix K does not match Eq. (17).');
end

%% 5. Concise output and formula visualization
fprintf('\nEq. (17) verified: K-K_eq17 is the zero matrix.\n');
fprintf('Rows: sigma_11, sigma_22, sigma_12; columns: n=1,2,3,4.\n\n');
disp('K =');
disp(K)

folder = fileparts(mfilename('fullpath'));
mat_file = fullfile(folder,'Supplementary_Material_5_Eq17_Symbolic_Results.mat');
txt_file = fullfile(folder,'Supplementary_Material_5_Eq17_Symbolic_Results.txt');
png_file = fullfile(folder,'Supplementary_Material_5_Eq17_Symbolic_Matrix.png');

% Store only the quantities needed to reproduce Eq. (17).
save(mat_file,'K','K_eq17','s_d_u','Sigma_1','Sigma_2','S','D','u');

fid = fopen(txt_file,'w','n','UTF-8');
if fid == -1
    error('Unable to create %s.',txt_file);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'Eq. (17): rows = sigma_11, sigma_22, sigma_12; columns = n=1,2,3,4\n\n');
fprintf(fid,'K =\n%s\n\n',char(K));
fprintf(fid,'LaTeX:\n%s\n',latex(K));
clear cleaner  % Close the text file before generating the visualization.

showEq17Matrix(K,png_file);
fprintf('Saved: %s\n',mat_file);
fprintf('Saved: %s\n',txt_file);
fprintf('Saved: %s\n',png_file);

%% Local function: display and export Eq. (17), following the original script
function showEq17Matrix(K,png_file)
    fig = figure('Color','w','Name','Equation (17)', ...
        'NumberTitle','off','MenuBar','none','Position',[100 100 1900 720], ...
        'Visible','off');
    ax = axes(fig,'Position',[0.02 0.03 0.96 0.94]);
    axis(ax,[0 1 0 1]);
    axis(ax,'off');
    hold(ax,'on');

    text(ax,0.5,0.95,'$[K_{in}(\Sigma)]$', ...
        'Interpreter','latex','FontSize',18,'FontWeight','bold', ...
        'HorizontalAlignment','center');

    x = [0.17 0.31 0.43 0.72];
    y = [0.68 0.43 0.18];
    row_name = {'$\sigma_{11}$','$\sigma_{22}$','$\sigma_{12}$'};
    for j = 1:4
        text(ax,x(j),0.84,sprintf('$n=%d$',j),'Interpreter','latex', ...
            'FontSize',14,'HorizontalAlignment','center');
    end
    for i = 1:3
        text(ax,0.055,y(i),row_name{i},'Interpreter','latex', ...
            'FontSize',15,'HorizontalAlignment','center');
        for j = 1:4
            font_size = 13;
            if j == 4
                font_size = 10;
            end
            text(ax,x(j),y(i),['$',latex(K(i,j)),'$'], ...
                'Interpreter','latex','FontSize',font_size, ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end

    % Light table guides make the row/column correspondence unambiguous.
    plot(ax,[0.10 0.98],[0.78 0.78],'k-','LineWidth',0.8);
    plot(ax,[0.10 0.98],[0.555 0.555],'Color',[0.75 0.75 0.75]);
    plot(ax,[0.10 0.98],[0.305 0.305],'Color',[0.75 0.75 0.75]);
    plot(ax,[0.10 0.98],[0.055 0.055],'k-','LineWidth',0.8);
    plot(ax,[0.10 0.10],[0.055 0.78],'k-','LineWidth',0.8);

    print(fig,png_file,'-dpng','-r200');
    if usejava('desktop')
        set(fig,'Visible','on');
    else
        close(fig);
    end
end
