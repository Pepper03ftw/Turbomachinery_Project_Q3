%% RUN_MEANLINE_DESIGN
clear; clc; close;

in = default_fan_inputs();

% Duty coefficients and all loss-model switches are defined in
% default_fan_inputs.m. Override them here only if you intentionally want a
% different study case.

out = meanline_fan_design(in);

disp(' ')
disp('=== PRELIMINARY FAN MEANLINE DESIGN ===')
fprintf('Solve mode        : %s\n', out.geometry.solve_mode);
fprintf('Mode note         : %s\n', out.geometry.solve_mode_note);
fprintf('phi (initial/used): %.4f / %.4f\n', out.coeffs.phi_initial, out.coeffs.phi_used);
fprintf('phi shift         : %+ .4f\n', out.coeffs.phi_used - out.coeffs.phi_initial);
fprintf('psi               : %.4f\n', out.coeffs.psi);
fprintf('R                 : %.4f\n', out.coeffs.R);
fprintf('Estimated eta_tt  : %.4f\n', out.performance.eta_tt_est);
fprintf('Estimated beta_tt : %.4f\n', out.coeffs.beta_tt_est);
fprintf('Power [MW]        : %.4f\n', out.performance.power_W/1e6);
fprintf('Mean radius [m]   : %.4f\n', out.geometry.rm);
fprintf('Hub/tip radii [m] : %.4f / %.4f\n', out.geometry.rh, out.geometry.rt);
fprintf('Annulus area [m^2]: %.4f\n', out.geometry.A1);
fprintf('Rotor/stator count: %d / %d\n', out.geometry.rotor_blades, out.geometry.stator_vanes);

fprintf('\nAngles [deg]\n')
disp(struct2table(out.angles_deg,'AsArray',true))

fprintf('Howell turning checks\n')
fprintf('Rotor : required = %.2f deg, limit = %.2f deg, residual = %+ .2f deg --> %s\n', ...
    out.howell.rotor_required_deg, out.howell.rotor_limit_deg, out.howell.rotor_residual_deg, passfail(out.howell.rotor_ok));
fprintf('Stator: required = %.2f deg, limit = %.2f deg, residual = %+ .2f deg --> %s\n', ...
    out.howell.stator_required_deg, out.howell.stator_limit_deg, out.howell.stator_residual_deg, passfail(out.howell.stator_ok));

fprintf('\nDiffusion-factor checks\n')
fprintf('Rotor : DF = %.3f, limit = %.3f, residual = %+ .3f --> %s\n', ...
    out.diffusion.rotor_lieblein, out.diffusion.limit, out.diffusion.rotor_residual, passfail(out.diffusion.rotor_ok));
fprintf('Stator: DF = %.3f, limit = %.3f, residual = %+ .3f --> %s\n', ...
    out.diffusion.stator_lieblein, out.diffusion.limit, out.diffusion.stator_residual, passfail(out.diffusion.stator_ok));

fprintf('\nLieblein optimum incidence/deviation estimates\n')
if isfield(out,'lieblein') && isfield(out.lieblein,'enabled') && out.lieblein.enabled
    print_lieblein_row('Rotor ', out.lieblein.rotor);
    print_lieblein_row('Stator', out.lieblein.stator);
else
    if isfield(out,'lieblein') && isfield(out.lieblein,'error')
        fprintf('Lieblein evaluation unavailable: %s\n', out.lieblein.error);
    else
        fprintf('Lieblein evaluation disabled.\n');
    end
end

fprintf('\nConstraint residual vector c <= 0\n')
fprintf('  Howell rotor  [deg] : %+ .4f\n', out.constraints.residuals.howell_rotor_deg);
fprintf('  Howell stator [deg] : %+ .4f\n', out.constraints.residuals.howell_stator_deg);
fprintf('  DF rotor      [-]   : %+ .4f\n', out.constraints.residuals.DF_rotor);
fprintf('  DF stator     [-]   : %+ .4f\n', out.constraints.residuals.DF_stator);
fprintf('Overall feasibility   : %s\n', passfail(out.constraints.is_feasible));

fprintf('\nLoss breakdown\n')
fprintf('loss_model.type label : %s\n', out.loss.type_label);
fprintf('Rotor  zeta total  : %.4f\n', out.loss.zetaR);
fprintf('  profile          : %.4f\n', out.loss.breakdown.rotor.profile);
fprintf('  trailing edge    : %.4f\n', out.loss.breakdown.rotor.trailing_edge);
fprintf('  shock            : %.4f\n', out.loss.breakdown.rotor.shock);
fprintf('  tip leakage      : %.4f\n', out.loss.breakdown.rotor.tip);
fprintf('  endwall          : %.4f\n', out.loss.breakdown.rotor.endwall);
fprintf('Stator zeta total  : %.4f\n', out.loss.zetaS);
fprintf('  profile          : %.4f\n', out.loss.breakdown.stator.profile);
fprintf('  trailing edge    : %.4f\n', out.loss.breakdown.stator.trailing_edge);
fprintf('  shock            : %.4f\n', out.loss.breakdown.stator.shock);
fprintf('  tip leakage      : %.4f\n', out.loss.breakdown.stator.tip);
fprintf('  endwall          : %.4f\n', out.loss.breakdown.stator.endwall);

% Optional parametric sweep over duty coefficients.
phi_vec = linspace(0.35,1.15,81);
psi_vec = linspace(0.20,0.80,81);

ETA  = nan(numel(psi_vec), numel(phi_vec));   % efficiency shown only for feasible points
GRES = nan(numel(psi_vec), numel(phi_vec));   % max inequality residual
FEAS = false(numel(psi_vec), numel(phi_vec)); % feasibility mask

% Store a reference design point separately
in_ref  = default_fan_inputs();
out_ref = meanline_fan_design(in_ref);

for i = 1:numel(psi_vec)
    for j = 1:numel(phi_vec)
        in = default_fan_inputs();
        in.phi = phi_vec(j);   % initial guess / selected phi
        in.psi = psi_vec(i);

        try
            out = meanline_fan_design(in);

            % Only keep points that returned finite values
            if isfinite(out.performance.eta_tt_est) && ...
               all(isfinite(out.constraints.residuals.all_ineq))

                GRES(i,j) = max(out.constraints.residuals.all_ineq);

                % Feasible means all inequality residuals <= 0
                if GRES(i,j) <= 0
                    FEAS(i,j) = true;
                    ETA(i,j)  = out.performance.eta_tt_est;
                end
            end

        catch
            % Leave ETA/GRES as NaN for failed or nonphysical points
            ETA(i,j)  = NaN;
            GRES(i,j) = NaN;
        end
    end
end

figure;
contourf(phi_vec, psi_vec, ETA, 40, 'LineColor', 'none');
cb = colorbar;
cb.Label.String = 'Estimated total-to-total efficiency, \eta_{tt} [-]';
hold on;

% Plot feasibility boundary only where residual field exists
if any(isfinite(GRES(:)))
    contour(phi_vec, psi_vec, GRES, [0 0], 'k', 'LineWidth', 1.5);
end

% Plot the reference design point
plot(out_ref.coeffs.phi_used, in_ref.psi, 'wo', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 7);

xlabel('\phi');
ylabel('\psi');
title('Feasible-region efficiency map');

% Optional: tighten color scaling to actual feasible range
eta_min = min(ETA(:), [], 'omitnan');
eta_max = max(ETA(:), [], 'omitnan');
if isfinite(eta_min) && isfinite(eta_max) && eta_max > eta_min
    clim([eta_min eta_max]);
end

grid on;
box on;

% Helpful console output
fprintf('Feasible points: %d / %d\n', nnz(FEAS), numel(FEAS));
fprintf('ETA min/max over feasible points: %.4f / %.4f\n', ...
    min(ETA(:),[],'omitnan'), max(ETA(:),[],'omitnan'));
fprintf('GRES min/max: %.4f / %.4f\n', ...
    min(GRES(:),[],'omitnan'), max(GRES(:),[],'omitnan'));

% figure;
% contourf(phi_vec, psi_vec, ETA, 20, 'LineColor', 'none');
% cb = colorbar;
% cb.Label.String = 'Estimated total-to-total efficiency, \eta_{tt} [-]';
% hold on;
% contour(phi_vec, psi_vec, double(feasible), [1 1], 'k', 'LineWidth', 1.2);
% plot(out.coeffs.phi_used, in.psi, 'wo', 'MarkerFaceColor', 'k');
% xlabel('\phi'); ylabel('\psi');
% title('Estimated meanline efficiency with feasibility boundary');
% grid on;

% Velocity triangles for the final, converged design.
plot_velocity_triangles(out);

function s = passfail(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
end
end

function print_lieblein_row(label, r)
if r.valid
    fprintf('%s: sigma = %.3f, t/c = %.4f, i = %+ .3f deg, delta = %+ .3f deg, camber = %.3f deg, stagger = %.3f deg\n', ...
        label, r.sigma, r.tmax_over_c, r.incidence_deg, r.deviation_deg, r.camber_deg, r.stagger_deg);
else
    fprintf('%s: invalid or outside digitized chart range. beta_in = %.3f deg, sigma = %.3f, t/c = %.4f\n', ...
        label, r.chi_in_deg, r.sigma, r.tmax_over_c);
end
end

function plot_velocity_triangles(out)
Vx = out.velocities.Vx;
U  = 2*pi*out.input.N_rpm/60 * out.geometry.rm;
Vt1 = out.velocities.Vt1;
Vt2 = out.velocities.Vt2;
Wt1 = out.velocities.Wt1;
Wt2 = out.velocities.Wt2;
Vt3 = out.velocities.Vt3;

figure;
hold on; axis equal; grid on;
quiver(0,0,Vx,Vt1,0,'LineWidth',1.5,'MaxHeadSize',0.5);
quiver(0,0,0,U,0,'LineWidth',1.5,'MaxHeadSize',0.5);
quiver(0,U,Vx,Wt1,0,'LineWidth',1.5,'MaxHeadSize',0.5);
quiver(0,0,Vx,Vt2,0,'LineWidth',1.5,'MaxHeadSize',0.5);
quiver(0,U,Vx,Wt2,0,'LineWidth',1.5,'MaxHeadSize',0.5);
quiver(0,0,Vx,Vt3,0,'LineWidth',1.5,'MaxHeadSize',0.5);
legend('V_1','U','W_1','V_2','W_2','V_3','Location','best');
xlabel('Axial component [m/s]');
ylabel('Tangential component [m/s]');
title('Final-design velocity triangles');
end
