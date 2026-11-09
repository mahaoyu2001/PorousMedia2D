function results = verify_pair_prefactor(outputDirectory)
% VERIFY_PAIR_PREFACTOR Independent spatial FFT check of C_i,0 = -1.
% No Abaqus solves, FEM fields, or statistical fits are used.
% Run: verify_pair_prefactor
% The exact result holds for the infinite zero-extended Kirsch kernel with
% isotropic step pair correlation and exclusion distance D >= 2a.
% Finite-box, periodic FFT results are checks, not finite-density corrections.
if nargin < 1
    outputDirectory = fileparts(mfilename('fullpath'));
end
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
settings = [16,1/32; 32,1/32; 32,1/64; 64,1/32];
diameters = [2,2.04];
rows = [];
for n = 1:size(settings,1)
    R = settings(n,1); dx = settings(n,2);
    N = round(2*R/dx);
    x = ((0:N-1)+0.5)*dx-R;
    [X,Y] = meshgrid(x,x);
    r2 = X.^2+Y.^2;
    u = (r2 >= 1)./r2;
    c2 = (X.^2-Y.^2)./r2; s2 = 2*X.*Y./r2;
    c4 = 2*c2.^2-1; s4 = 2*s2.*c2;
    B = -0.5*u.*(3*u-2);
    clear X Y r2
    q = ceil(max(diameters)/dx)+1;
    shifts = -q:q; indices = mod(shifts,N)+1;
    [sx,sy] = meshgrid(shifts*dx,shifts*dx);
    weights = cell(size(diameters));
    offsets = ((0:7)+0.5)/8-0.5;
    for d = 1:numel(diameters)
        w = zeros(size(sx));
        for ox = offsets
            for oy = offsets
                w = w+double((sx+ox*dx).^2+(sy+oy*dx).^2 < diameters(d)^2)/64;
            end
        end
        weights{d} = w;
    end
    for i = 1:3
        switch i
            case 1, h = 0.5*u.*c2+B.*c4;
            case 2, h = 1.5*u.*c2-B.*c4;
            case 3, h = -0.5*u.*s2+B.*s4;
        end
        I = sum(h.^2,'all')*dx^2;
        exactI = pi*[0.25,1.25,0.25];
        F = fft2(h);
        correlation = real(ifft2(abs(F).^2))*dx^2;
        local = correlation(indices,indices);
        for d = 1:numel(diameters)
            coefficient = -sum(local.*weights{d},'all')*dx^2/(pi*I);
            rows(end+1,:) = [R,dx,N,i,diameters(d),I,exactI(i),100*(I/exactI(i)-1),coefficient,100*abs(coefficient+1)]; %#ok<AGROW>
        end
        clear h F correlation
    end
    fprintf('R/a=%g, dx/a=%g completed.\n',R,dx);
    clear u c2 s2 c4 s4 B
end
results = array2table(rows,'VariableNames', ...
    {'R_over_a','dx_over_a','N','component_index','D_over_a', ...
    'I_numeric','I_exact','I_error_percent','C0_numeric','C0_error_percent'});
assert(all(results.C0_error_percent < 3), ...
    'FFT verification differs from the exact result by more than 3 percent.');
writetable(results,fullfile(outputDirectory,'Pair_Prefactor_FFT_Convergence.csv'));
disp(results);
end
