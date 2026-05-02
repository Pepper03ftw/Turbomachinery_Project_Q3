%% RUN_MEANLINE_DESIGN
clear; clc; close all;

in = default_fan_inputs();
out = meanline_fan_design(in);

disp(' ')
disp('=== PRELIMINARY FAN MEANLINE DESIGN ===')
fprintf('Solve mode        : %s\n', out.geometry.solve_mode);
fprintf('Mode note         : %s\n', out.geometry.solve_mode_note);
fprintf('phi (initial/used): %.4f / %.4f\n', out.coeffs.phi_initial, out.coeffs.phi_used);
fprintf('phi shift         : %+ .4f\n', out.coeffs.phi_used - out.coeffs.phi_initial);
fprintf('psi               : %.4f\n', out.coeffs.psi);
fprintf('R                 : %.4f\n', out.coeffs.R);
fprintf('R_R recovery ratio : %.4f\n', out.coeffs.RR);
if isfield(out,'R_search')
    fprintf('R search feasible : %s (%d/%d feasible candidates)\n', ...
        yesno(out.R_search.found_feasible), nnz(out.R_search.feasible), numel(out.R_search.feasible));
end
if isfield(out,'beta_tt_search')
    fprintf('beta_tt search mode : %s (%d/%d feasible candidates)\n', ...
        out.beta_tt_search.mode, nnz(out.beta_tt_search.feasible), numel(out.beta_tt_search.feasible));
    fprintf('selected beta_tt    : %.4f\n', out.beta_tt_search.selected_beta_tt);
    fprintf('beta objective J    : %.6f\n', out.beta_tt_search.objective);
end
if isfield(out,'sweep')
    fprintf('Sweep selected    : %s (%d/%d feasible candidates)\n', ...
        out.sweep.selected_source, out.sweep.feasible_count, out.sweep.total_count);
    fprintf('Sweep objective J : %.6f\n', out.sweep.objective);
end
fprintf('Estimated eta_tt  : %.4f\n', out.performance.eta_tt_est);
fprintf('Estimated beta_tt : %.4f\n', out.coeffs.beta_tt_est);
fprintf('Power [MW]        : %.4f\n', out.performance.power_W/1e6);
fprintf('Mean radius [m]   : %.4f\n', out.geometry.rm);
fprintf('Hub/tip radii [m] : %.4f / %.4f\n', out.geometry.rh, out.geometry.rt);
if isfield(out,'design_basis') && isfield(out.design_basis,'hub_to_tip_actual')
    fprintf('Hub-to-tip ratio   : %.4f (input/fallback %.4f)\n', ...
        out.design_basis.hub_to_tip_actual, out.design_basis.hub_to_tip_input);
end
fprintf('Annulus area [m^2]: %.4f\n', out.geometry.A1);
fprintf('Rotor/stator count: %d / %d\n', out.geometry.rotor_blades, out.geometry.stator_vanes);

fprintf('\nMach / work / mass-flow design-basis checks\n')
fprintf('Mx target/actual        : %.4f / %.4f\n', out.design_basis.Mx_target, out.design_basis.Mx_actual);
fprintf('Mrel,tip target/actual  : %.4f / %.4f\n', out.design_basis.Mrel_tip_target, out.design_basis.Mrel_tip_actual);
fprintf('eta work guess          : %.4f\n', out.design_basis.eta_work_guess);
fprintf('computed phi range      : %.4f <= %.4f <= %.4f\n', out.design_basis.phi_min, out.coeffs.phi_used, out.design_basis.phi_max);
fprintf('computed psi range      : %.4f <= %.4f <= %.4f\n', out.design_basis.psi_min, out.coeffs.psi, out.design_basis.psi_max);
fprintf('meanline mdot target/actual [kg/s]: %.3f / %.3f  (error %+ .3f kg/s, %+ .3f%%)\n', ...
    out.design_basis.mdot_target, out.design_basis.mdot_meanline, ...
    out.design_basis.mdot_error, 100*out.design_basis.mdot_rel_error);
fprintf('fixed tolerance rel      : %.3f%%\n', 100*out.design_basis.fixed_tolerance_rel);
fprintf('hub-to-tip target/actual: %.4f / %.4f\n', out.design_basis.hub_to_tip_input, out.design_basis.hub_to_tip_actual);
fprintf('R_R soft minimum/actual  : %.4f / %.4f --> %s\n', out.stability.RR_soft_min, out.stability.RR, passfail(out.stability.RR_soft_ok));

fprintf('\nAngles [deg]\n')
disp(struct2table(out.angles_deg,'AsArray',true))

fprintf('Howell turning checks\n')
fprintf('Rotor : required = %.2f deg, limit = %.2f deg, residual = %+ .2f deg --> %s\n', ...
    out.howell.rotor_required_deg, out.howell.rotor_limit_deg, out.howell.rotor_residual_deg, passfail(out.howell.rotor_ok));
fprintf('Stator: required = %.2f deg, limit = %.2f deg, residual = %+ .2f deg --> %s\n', ...
    out.howell.stator_required_deg, out.howell.stator_limit_deg, out.howell.stator_residual_deg, passfail(out.howell.stator_ok));

fprintf('\nDiffusion-factor checks\n')
fprintf('Rotor : DF = %.3f, range = [%.3f, %.3f], high residual = %+ .3f, low residual = %+ .3f --> %s\n', ...
    out.diffusion.rotor_lieblein, out.diffusion.minimum, out.diffusion.limit, ...
    out.diffusion.rotor_residual, out.diffusion.rotor_low_residual, passfail(out.diffusion.rotor_ok));
fprintf('Stator: DF = %.3f, range = [%.3f, %.3f], high residual = %+ .3f, low residual = %+ .3f --> %s\n', ...
    out.diffusion.stator_lieblein, out.diffusion.minimum, out.diffusion.limit, ...
    out.diffusion.stator_residual, out.diffusion.stator_low_residual, passfail(out.diffusion.stator_ok));

fprintf('\nLieblein optimum incidence/deviation and signed stagger estimates\n')
if isfield(out,'lieblein') && isfield(out.lieblein,'enabled') && out.lieblein.enabled
    print_lieblein_row('Rotor ', out.lieblein.rotor);
    print_lieblein_row('Stator', out.lieblein.stator);
    if isfield(out.loss,'loss_model_used') && isfield(out.loss.loss_model_used,'lieblein_stagger_auto_applied')
        fprintf('Auto-fed to loss-model stagger modes: %s\n', yesno(out.loss.loss_model_used.lieblein_stagger_auto_applied));
    end
else
    if isfield(out,'lieblein') && isfield(out.lieblein,'error')
        fprintf('Lieblein evaluation unavailable: %s\n', out.lieblein.error);
    else
        fprintf('Lieblein evaluation disabled.\n');
    end
end

fprintf('\nAxial chord projections for MEANGEN input\n')
if isfield(out.geometry,'rotor_axial_chord')
    fprintf('Rotor : c = %.5f m, Cx = %.5f m, stagger used = %+ .3f deg\n', ...
        out.geometry.rotor_chord, ...
        out.geometry.rotor_axial_chord, ...
        out.geometry.rotor_stagger_for_axial_chord_deg);
end
if isfield(out.geometry,'stator_axial_chord')
    fprintf('Stator: c = %.5f m, Cx = %.5f m, stagger used = %+ .3f deg\n', ...
        out.geometry.stator_chord, ...
        out.geometry.stator_axial_chord, ...
        out.geometry.stator_stagger_for_axial_chord_deg);
end
if isfield(out.geometry,'rotor_axial_chord') && isfield(out.geometry,'stator_axial_chord')
    fprintf('MEANGEN axial chord line: %.6f %.6f\n', ...
        out.geometry.rotor_axial_chord, out.geometry.stator_axial_chord);
end

fprintf('\nConstraint and diagnostic residuals\n')
fprintf('  Howell rotor  [deg] : %+ .4f\n', out.constraints.residuals.howell_rotor_deg);
fprintf('  Howell stator [deg] : %+ .4f\n', out.constraints.residuals.howell_stator_deg);
fprintf('  DF rotor high [-]   : %+ .4f\n', out.constraints.residuals.DF_rotor);
fprintf('  DF stator high[-]   : %+ .4f\n', out.constraints.residuals.DF_stator);
fprintf('  DF rotor low  [-]   : %+ .4f\n', out.constraints.residuals.DF_rotor_low);
fprintf('  DF stator low [-]   : %+ .4f\n', out.constraints.residuals.DF_stator_low);
fprintf('  phi min/max   [-]   : %+ .4f / %+ .4f\n', out.constraints.residuals.phi_min, out.constraints.residuals.phi_max);
fprintf('  psi min/max   [-]   : %+ .4f / %+ .4f\n', out.constraints.residuals.psi_min, out.constraints.residuals.psi_max);
fprintf('  mdot abs      [kg/s]: %+ .4f\n', out.constraints.residuals.mdot_abs);
fprintf('  mdot rel      [-]   : %+ .4f\n', out.constraints.residuals.mdot_rel);
if isfield(out.constraints.residuals,'beta_tt_min')
    fprintf('  beta_tt min/max [-]: %+ .4f / %+ .4f\n', ...
        out.constraints.residuals.beta_tt_min, out.constraints.residuals.beta_tt_max);
end
fprintf('  Lieblein rotor [-]  : %+ .4f\n', out.constraints.residuals.lieblein_rotor);
fprintf('  Lieblein stator[-]  : %+ .4f\n', out.constraints.residuals.lieblein_stator);
fprintf('  Mx residual   [-]   : %+ .4e\n', out.constraints.residuals.Mx_design);
fprintf('  Mrel tip res. [-]   : %+ .4e\n', out.constraints.residuals.Mrel_tip_design);
fprintf('Overall hard-aero feasibility : %s\n', passfail(out.constraints.is_feasible));
if isfield(out.constraints.residuals,'project_ineq')
    fprintf('Project matching diagnostics are soft by default. Max project residual: %+ .4e\n', max(out.constraints.residuals.project_ineq));
end
if isfield(out.constraints.residuals,'beta_ineq')
    fprintf('Beta_tt range diagnostics are soft by default. Max beta residual   : %+ .4e\n', max(out.constraints.residuals.beta_ineq));
end

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

if isfield(out,'sweep')
    fprintf('\nDuty-coefficient sweep was executed by meanline_fan_design.\n');
elseif isfield(in,'plot') && isfield(in.plot,'enable_duty_sweep') && in.plot.enable_duty_sweep
    fprintf('\nPlot-only duty sweep flag is enabled, but the design sweep is controlled by in.design.solve_mode.\n');
else
    fprintf('\nNo additional plot-only duty sweep requested.\n');
end

plot_velocity_triangles(out);

function print_lieblein_row(label, r)
if ~isfield(r,'valid') || ~r.valid
    fprintf('%s: invalid / outside digitized chart range.\n', label);
    return
end
fprintf('%s: sigma = %.3f, t/c = %.4f, i = %+ .3f deg, delta = %+ .3f deg, camber = %.3f deg\n', ...
    label, r.sigma, r.tmax_over_c, r.incidence_deg, r.deviation_deg, r.camber_deg);
if isfield(r,'stagger_signed_deg')
    fprintf('       signed metal in/out = %+ .3f / %+ .3f deg, signed stagger = %+ .3f deg\n', ...
        r.metal_in_signed_deg, r.metal_out_signed_deg, r.stagger_signed_deg);
    fprintf('       magnitude metal in/out = %.3f / %.3f deg, magnitude stagger = %.3f deg\n', ...
        r.metal_in_magnitude_deg, r.metal_out_magnitude_deg, r.stagger_magnitude_deg);
else
    fprintf('       metal in/out = %.3f / %.3f deg, stagger = %.3f deg\n', ...
        r.metal_in_deg, r.metal_out_deg, r.stagger_deg);
end
end

function s = yesno(tf)
if tf
    s = 'YES';
else
    s = 'NO';
end
end

function s = passfail(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
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
quiver(0,U,Vx,Wt2,0,'LineWidth',1.5,'MaxHeadSize',0.5);
quiver(0,0,Vx,Vt2,0,'LineWidth',1.5,'MaxHeadSize',0.5);
quiver(0,0,Vx,Vt3,0,'LineWidth',1.5,'MaxHeadSize',0.5);
legend('V_1','U','W_1','W_2','V_2','V_3','Location','bestoutside');
xlabel('Axial component [m/s]');
ylabel('Tangential component [m/s]');
title('Velocity Triangles at Mean Radius');
end
