%% contraction_invariance_check.m
% Numerical verification that the forward and backward MPCI sweep operators
% (as built in odeMPCI.m / odeMPCIrev.m) are similar up to sign, per
% contraction_invariance_analysis.md Section 6. Operator-construction lines
% below are copied verbatim from odeMPCI.m:123-139 / odeMPCIrev.m:135-152 so
% the test matches production code, not an idealized reconstruction.
%
% Outputs:
%   contraction_invariance_check.csv  -- Summary + Eigenvalues + SingularValues sections
%   contraction_invariance_check.tex  -- same three sections as LaTeX tables
% plus a console summary.

order = 16;      % options.N as used throughout the pipeline (N_cheb)
N = order - 1;   % internal Chebyshev degree, matching the N=N-1 decrement in both files

tau = fliplr(cos((0:N)*(pi/N)));  % CGL nodes

%% Shared operator pieces (identical in odeMPCI.m and odeMPCIrev.m)
W = eye(length(tau));
W(1,1) = 0.5;
W(end,end) = 0.5;
T = cos((0:N-1)'*acos(tau))';    % truncated, (N+1) x N -- feeds A
s_ = 1./(4:2:2*N);
S_3 = zeros(N+1);
S_3(2:end,2:end) = diag([-.5,-s_(1:end-1)],1);
S_2 = diag([1,s_],-1);
S_1 = S_2 + S_3;
S = S_1(:,1:N);
S(1,:) = [1/4, zeros(1,N-1)];      % matches Woollands; 1/4 at 0-indexed (0,0)
A = (T'*W*T)\T'*W;
T_full = cos((0:N)'*acos(tau))';  % full, (N+1) x (N+1) -- for evaluation

%% Forward vs backward: only the boundary row differs
Tm1_fwd = cos((0:N).*pi());      % T_k(-1) = (-1)^k
Tm1_bwd = ones(1,N+1);           % T_k(+1) = 1 for all k
L_fwd = [Tm1_fwd; zeros(N,N+1)];
L_bwd = [Tm1_bwd; zeros(N,N+1)];
P_fwd = (eye(N+1) - L_fwd) * S;
P_bwd = (eye(N+1) - L_bwd) * S;

Kf = T_full * P_fwd * A;
Kb = T_full * P_bwd * A;

%% Checks (contraction_invariance_analysis.md Section 6)
R = flipud(eye(size(Kf,1)));
similarity_residual = norm(Kb + R*Kf*R, 'fro') / norm(Kf, 'fro');   % (a)

eig_f = eig(Kf);
eig_b = eig(Kb);
rho_f  = max(abs(eig_f));                                          % (c)
rho_b  = max(abs(eig_b));
% Note: cond(Kf)/cond(Kb) are NOT reported. Kf and Kb are singular by
% construction -- the boundary-row projection (I-L) zeroes out one degree
% of freedom in the (N+1)-dimensional space, giving each an exact zero
% eigenvalue (confirmed numerically: eig index 1 is ~1e-16). cond() then
% divides by that near-zero value, producing floating-point noise (not a
% real forward/backward discrepancy), so it's omitted rather than reported.
%
% Instead, compare singular values directly (no division involved). Since
% Kb = -R*Kf*R with R orthogonal, Kb'*Kb = R*(Kf'*Kf)*R is an orthogonal
% similarity transform of Kf'*Kf, so eig(Kb'*Kb) = eig(Kf'*Kf) exactly --
% svd(Kf) and svd(Kb) should match to machine precision, individually.
sv_f = svd(Kf);
sv_b = svd(Kb);
[sv_f_sorted, ~] = sort(sv_f, 'descend');
[sv_b_sorted, ~] = sort(sv_b, 'descend');
sv_abs_diff = abs(sv_f_sorted - sv_b_sorted);

%% Pair eigenvalues for reporting: eig_f should match -eig_b (possibly reordered)
neg_eig_b = -eig_b;
[eig_f_sorted, ~] = sortrows([real(eig_f), imag(eig_f)]);
[eig_b_sorted, ~] = sortrows([real(neg_eig_b), imag(neg_eig_b)]);
eig_abs_diff = sqrt((eig_f_sorted(:,1) - eig_b_sorted(:,1)).^2 + ...
                     (eig_f_sorted(:,2) - eig_b_sorted(:,2)).^2);

%% Console report
fprintf('=== Contraction Invariance Check (N=%d) ===\n\n', order);
fprintf('SimilarityResidual : %.3e  (expect ~1e-16)\n', similarity_residual);
fprintf('SpectralRadius     : forward=%.10f  backward=%.10f  |diff|=%.3e\n', ...
    rho_f, rho_b, abs(rho_f - rho_b));
fprintf('MaxSingularValue   : forward=%.10f  backward=%.10f  |diff|=%.3e\n', ...
    sv_f_sorted(1), sv_b_sorted(1), sv_abs_diff(1));
fprintf('\nFull eigenvalue/singular-value comparisons written to .csv / .tex.\n\n');

%% Write CSV (Summary + Eigenvalues)
fid_csv = fopen('contraction_invariance_check.csv', 'w');
fprintf(fid_csv, 'Summary\n');
fprintf(fid_csv, 'Metric,Forward,Backward,AbsDiff\n');
fprintf(fid_csv, 'SpectralRadius,%.10f,%.10f,%.3e\n', rho_f, rho_b, abs(rho_f - rho_b));
fprintf(fid_csv, 'MaxSingularValue,%.10f,%.10f,%.3e\n', sv_f_sorted(1), sv_b_sorted(1), sv_abs_diff(1));
fprintf(fid_csv, 'SimilarityResidual,%.3e,,\n', similarity_residual);
fprintf(fid_csv, '\n');
fprintf(fid_csv, 'Eigenvalues\n');
fprintf(fid_csv, 'Index,Re(eig_f),Im(eig_f),Re(-eig_b),Im(-eig_b),AbsDiff\n');
for k = 1:length(eig_f_sorted)
    fprintf(fid_csv, '%d,%.10f,%.10f,%.10f,%.10f,%.3e\n', k, ...
        eig_f_sorted(k,1), eig_f_sorted(k,2), eig_b_sorted(k,1), eig_b_sorted(k,2), eig_abs_diff(k));
end
fprintf(fid_csv, '\n');
fprintf(fid_csv, 'SingularValues\n');
fprintf(fid_csv, 'Index,sigma_f,sigma_b,AbsDiff\n');
for k = 1:length(sv_f_sorted)
    fprintf(fid_csv, '%d,%.10f,%.10f,%.3e\n', k, sv_f_sorted(k), sv_b_sorted(k), sv_abs_diff(k));
end
fclose(fid_csv);
fprintf('CSV saved to contraction_invariance_check.csv\n');

%% Write LaTeX tables (Summary + Eigenvalues)
fid_tex = fopen('contraction_invariance_check.tex', 'w');
fprintf(fid_tex, '%% Summary (N=%d)\n', order);
fprintf(fid_tex, '\\begin{tabular}{l r r r}\n\\hline\n');
fprintf(fid_tex, 'Metric & Forward & Backward & $|$Diff$|$ \\\\\n\\hline\n');
fprintf(fid_tex, 'Spectral radius & %.10f & %.10f & %.3e \\\\\n', rho_f, rho_b, abs(rho_f - rho_b));
fprintf(fid_tex, 'Max singular value & %.10f & %.10f & %.3e \\\\\n', sv_f_sorted(1), sv_b_sorted(1), sv_abs_diff(1));
fprintf(fid_tex, 'Similarity residual & %.3e & --- & --- \\\\\n', similarity_residual);
fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

fprintf(fid_tex, '%% Eigenvalue comparison (N=%d)\n', order);
fprintf(fid_tex, '\\begin{tabular}{l r r r r r}\n\\hline\n');
fprintf(fid_tex, 'Index & Re($\\lambda_f$) & Im($\\lambda_f$) & Re($-\\lambda_b$) & Im($-\\lambda_b$) & $|\\Delta|$ \\\\\n\\hline\n');
for k = 1:length(eig_f_sorted)
    fprintf(fid_tex, '%d & %.10f & %.10f & %.10f & %.10f & %.3e \\\\\n', k, ...
        eig_f_sorted(k,1), eig_f_sorted(k,2), eig_b_sorted(k,1), eig_b_sorted(k,2), eig_abs_diff(k));
end
fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

fprintf(fid_tex, '%% Singular value comparison (N=%d)\n', order);
fprintf(fid_tex, '\\begin{tabular}{l r r r}\n\\hline\n');
fprintf(fid_tex, 'Index & $\\sigma_f$ & $\\sigma_b$ & $|\\Delta|$ \\\\\n\\hline\n');
for k = 1:length(sv_f_sorted)
    fprintf(fid_tex, '%d & %.10f & %.10f & %.3e \\\\\n', k, sv_f_sorted(k), sv_b_sorted(k), sv_abs_diff(k));
end
fprintf(fid_tex, '\\hline\n\\end{tabular}\n');
fclose(fid_tex);
fprintf('LaTeX tables saved to contraction_invariance_check.tex\n');
