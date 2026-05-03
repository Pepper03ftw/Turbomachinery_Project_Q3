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
if isfield(out.coeffs,'RR_signed')
    fprintf('R_R signed diagnostic: %.4f  (not used for stability ranking unless RR_angle_convention=''signed'')\n', out.coeffs.RR_signed);
    fprintf('R_R angle convention : %s\n', out.coeffs.RR_angle_convention);
    fprintf('R_R angles used [deg]: alpha1=%.3f, beta2_RR=%.3f\n', out.coeffs.alpha1_RR_deg, out.coeffs.beta2_RR_deg);
end
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
if isfield(out,'geometry') && isfield(out.geometry,'rm_continuity') && ...
        isfield(out.geometry,'fixed_mean_radius_duty_active') && out.geometry.fixed_mean_radius_duty_active
    fprintf('Fixed-radius duty mode: continuity-reference r_m = %.4f m, selected r_m = %.4f m (shift %+ .2f%%)\n', ...
        out.geometry.rm_continuity, out.geometry.rm, 100*(out.geometry.rm - out.geometry.rm_continuity)/max(out.geometry.rm_continuity,1e-12));
end
if isfield(out,'design_basis') && isfield(out.design_basis,'tipmach_radius_match') && out.design_basis.tipmach_radius_match.enabled
    tm = out.design_basis.tipmach_radius_match;
    fprintf('Continuity r_mean [m]: %.4f  --> tip-Mach matched r_mean [m]: %.4f  (shift %+ .2f%%, status: %s)\n', ...
        tm.rm_continuity, tm.rm_selected, 100*tm.relative_shift, tm.status);
    fprintf('Tip-Mach match residual: target %.4f, matched %.4f, rel. error %+ .3f%%\n', ...
        tm.Mrel_tip_target, tm.Mrel_tip_selected, 100*tm.Mrel_tip_rel_error);
end
fprintf('Hub/tip radii [m] : %.4f / %.4f\n', out.geometry.rh, out.geometry.rt);
if isfield(out,'design_basis') && isfield(out.design_basis,'rho_radius_eval')
    fprintf('Radius eval rho [kg/m^3] : %.5f\n', out.design_basis.rho_radius_eval);
    fprintf('Radius eval T/p [K/Pa]   : %.3f / %.1f\n', out.design_basis.T_radius_eval, out.design_basis.p_radius_eval);
    fprintf('Radius eval Vx/A [m/s,m^2]: %.3f / %.5f\n', out.design_basis.Vx_radius_eval, out.design_basis.A_radius_eval);
end
if isfield(out,'design_basis') && isfield(out.design_basis,'continuity_ambient_direct')
    ca = out.design_basis.continuity_active;
    cd = out.design_basis.continuity_ambient_direct;
    fprintf('Continuity radius check: active rm=%.4f m using rho=%.5f, Vx=%.2f; ambient-direct rm=%.4f m using rho=%.5f, Vx=%.2f\n', ...
        ca.rm, ca.rho, ca.Vx, cd.rm, cd.rho, cd.Vx);
end
if isfield(out,'design_basis') && isfield(out.design_basis,'hub_to_tip_actual')
    fprintf('Hub-to-tip ratio   : %.4f (input/fallback %.4f)\n', ...
        out.design_basis.hub_to_tip_actual, out.design_basis.hub_to_tip_input);
end
fprintf('Annulus area [m^2]: %.4f\n', out.geometry.A1);
if isfield(out,'blade_count_table')
    fprintf('Rotor/stator reference count: %d / %d (not objective-selected)\n', out.geometry.rotor_blades, out.geometry.stator_vanes);
else
    fprintf('Rotor/stator count: %d / %d\n', out.geometry.rotor_blades, out.geometry.stator_vanes);
end
if isfield(out.geometry,'chord_closure')
    fprintf('Chord closure     : %s\n', out.geometry.chord_closure.mode);
end
if isfield(out.geometry,'rotor_t_over_c')
    fprintf('Rotor/stator t/c : %.4f / %.4f\n', out.geometry.rotor_t_over_c, out.geometry.stator_t_over_c);
end
fprintf('Rotor/stator pitch [m]: %.5f / %.5f\n', out.geometry.rotor_pitch, out.geometry.stator_pitch);
fprintf('Rotor/stator solidity : %.4f / %.4f\n', out.geometry.rotor_solidity, out.geometry.stator_solidity);
if isfield(out,'diffusion') && isfield(out.diffusion,'rotor_solidity_design')
    if isfield(out.diffusion.rotor_solidity_design,'selected_reason')
        fprintf('Rotor solidity sizing : %s, sigma_DF=%.4f, DF=%.4f\n', ...
            out.diffusion.rotor_solidity_design.selected_reason, ...
            out.diffusion.rotor_solidity_design.sigma_initial_from_DF_target, ...
            out.diffusion.rotor_solidity_design.DF);
    end
    if isfield(out.diffusion.stator_solidity_design,'selected_reason')
        fprintf('Stator solidity sizing: %s, sigma_DF=%.4f, DF=%.4f\n', ...
            out.diffusion.stator_solidity_design.selected_reason, ...
            out.diffusion.stator_solidity_design.sigma_initial_from_DF_target, ...
            out.diffusion.stator_solidity_design.DF);
    end
end

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
if isfield(out.stability,'RR_soft_residual')
    fprintf('R_R soft residual        : %+ .4f  (positive means below soft minimum)\n', out.stability.RR_soft_residual);
end
if isfield(out.stability,'RR_signed')
    fprintf('R_R signed diagnostic    : %.4f using signed alpha1=%.3f deg, beta2=%.3f deg\n', ...
        out.stability.RR_signed, out.stability.alpha1_signed_deg, out.stability.beta2_signed_deg);
    fprintf('R_R metric angle basis   : %s; alpha1=%.3f deg, beta2_RR=%.3f deg\n', ...
        out.stability.RR_angle_convention, out.stability.alpha1_RR_deg, out.stability.beta2_RR_deg);
    if isfield(out.stability,'beta2_RR_positive_ok')
        fprintf('R_R beta2 convention ok  : %s  (residual %+ .4e rad)\n', ...
            passfail(out.stability.beta2_RR_positive_ok), out.stability.beta2_RR_positive_residual);
    end
end

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
    fprintf('Rotor : c = %.5f m, t = %.5f m, Cx = %.5f m, stagger used = %+ .3f deg\n', ...
        out.geometry.rotor_chord, ...
        out.geometry.rotor_thickness, ...
        out.geometry.rotor_axial_chord, ...
        out.geometry.rotor_stagger_for_axial_chord_deg);
end
if isfield(out.geometry,'stator_axial_chord')
    fprintf('Stator: c = %.5f m, t = %.5f m, Cx = %.5f m, stagger used = %+ .3f deg\n', ...
        out.geometry.stator_chord, ...
        out.geometry.stator_thickness, ...
        out.geometry.stator_axial_chord, ...
        out.geometry.stator_stagger_for_axial_chord_deg);
end
if isfield(out.geometry,'rotor_axial_chord') && isfield(out.geometry,'stator_axial_chord')
    fprintf('MEANGEN axial chord line: %.6f %.6f\n', ...
        out.geometry.rotor_axial_chord, out.geometry.stator_axial_chord);
end

if isfield(out,'blade_count_table')
    fprintf('\nBlade-count geometry table\n')
    fprintf('Blade count was not used to select the meanline objective; table is post-processing only.\n');
    Tbc = out.blade_count_table;
    fprintf('Generated %d rotor/stator count combinations.\n', height(Tbc));
    fprintf('Rotor chord range [m] : %.5f to %.5f\n', min(Tbc.chord_R), max(Tbc.chord_R));
    fprintf('Stator chord range [m]: %.5f to %.5f\n', min(Tbc.chord_S), max(Tbc.chord_S));
    fprintf('Rotor Cx range [m]    : %.5f to %.5f\n', min(Tbc.Cx_R), max(Tbc.Cx_R));
    fprintf('Stator Cx range [m]   : %.5f to %.5f\n', min(Tbc.Cx_S), max(Tbc.Cx_S));
    nshow = min(10,height(Tbc));
    disp(Tbc(1:nshow,:));
    try
        writetable(Tbc,'blade_count_geometry_table.csv');
        fprintf('Saved full blade-count table to blade_count_geometry_table.csv\n');
    catch ME
        warning('Could not save blade-count table: %s', ME.message);
    end
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
if isfield(out.constraints.residuals,'mdot_excess_rel')
    fprintf('  mdot excess beyond tol [-]: %+ .4f\n', out.constraints.residuals.mdot_excess_rel);
end
if isfield(out.constraints.residuals,'beta_tt_min')
    fprintf('  beta_tt min/max [-]: %+ .4f / %+ .4f\n', ...
        out.constraints.residuals.beta_tt_min, out.constraints.residuals.beta_tt_max);
end
fprintf('  Lieblein rotor [-]  : %+ .4f\n', out.constraints.residuals.lieblein_rotor);
fprintf('  Lieblein stator[-]  : %+ .4f\n', out.constraints.residuals.lieblein_stator);
fprintf('  Mx residual   [-]   : %+ .4e\n', out.constraints.residuals.Mx_design);
fprintf('  Mrel tip res. [-]   : %+ .4e\n', out.constraints.residuals.Mrel_tip_design);
if isfield(out.constraints.residuals,'Mrel_tip_excess_rel')
    fprintf('  Mrel tip excess beyond tol [-]: %+ .4e\n', out.constraints.residuals.Mrel_tip_excess_rel);
end
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
if isfield(out,'loss') && isfield(out.loss,'eta_penalty')
    fprintf('Estimated eta penalty from endwall terms : %.5f (%.3f percentage points)\n', ...
        out.loss.eta_penalty_endwall, 100*out.loss.eta_penalty_endwall);
    fprintf('Eta_tt if only endwall terms are removed : %.4f\n', out.performance.eta_tt_no_endwall_est);
    fprintf('Beta_tt if only endwall terms removed   : %.4f\n', out.performance.beta_tt_no_endwall_est);
end

if isfield(out,'sweep')
    fprintf('\nDuty-coefficient sweep was executed by meanline_fan_design.\n');
elseif isfield(in,'plot') && isfield(in.plot,'enable_duty_sweep') && in.plot.enable_duty_sweep
    fprintf('\nPlot-only duty sweep flag is enabled, but the design sweep is controlled by in.design.solve_mode.\n');
else
    fprintf('\nNo additional plot-only duty sweep requested.\n');
end

if isfield(in,'plot') && isfield(in.plot,'enable_beta_mdot_map') && in.plot.enable_beta_mdot_map
    plot_beta_mdot_from_sweep(out, in);
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

function plot_beta_mdot_from_sweep(out, in)
if isfield(out,'speedline') && isfield(out.speedline,'computed') && ~out.speedline.computed
    fprintf('beta_tt-mdot speedline not computed: %s\n', out.speedline.reason);
    return
end
if ~isfield(out,'sweep') || ~isfield(out.sweep,'mdot') || ~isfield(out.sweep,'beta_tt')
    fprintf('beta_tt-mdot plot skipped: no sweep history available.\n');
    return
end
mdot = out.sweep.mdot(:);
beta = out.sweep.beta_tt(:);
mask = isfinite(mdot) & isfinite(beta);
if isfield(in.plot,'beta_mdot_only_feasible') && in.plot.beta_mdot_only_feasible && isfield(out.sweep,'feasible')
    mask = mask & logical(out.sweep.feasible(:));
end
if nnz(mask) < 2
    fprintf('beta_tt-mdot plot skipped: fewer than two valid sweep points.\n');
    return
end
mdot = mdot(mask);
beta = beta(mask);
[mdot, idx] = sort(mdot);
beta = beta(idx);

figure;
hold on; grid on; box on;
plot(mdot, beta, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
plot(out.design_basis.mdot_meanline, out.coeffs.beta_tt_est, 's', 'MarkerSize', 8, 'LineWidth', 1.5);
if isfield(out.input,'mdot')
    xline(out.input.mdot, '--', 'target \dot{m}');
end
if isfield(out.input,'design') && isfield(out.input.design,'beta_tt_min')
    yline(out.input.design.beta_tt_min, '--', '\beta_{tt,min}');
end
if isfield(out.input,'design') && isfield(out.input.design,'beta_tt_max')
    yline(out.input.design.beta_tt_max, '--', '\beta_{tt,max}');
end
xlabel('mass flow \dot{m} [kg/s]');
ylabel('\beta_{tt} [-]');
title('\beta_{tt} vs \dot{m} from preliminary sweep history');
legend('sweep candidates','selected candidate','Location','best');

if max(mdot)-min(mdot) < 1e-6
    warning(['beta_tt-mdot sweep history has almost no mdot spread. ', ...
             'With continuity-fixed radius and fixed Mx, mdot is fixed by construction. ', ...
             'A true speedline/choke check requires fixed-geometry off-design points, e.g. MULTALL runs.']);
end

if isfield(in.plot,'beta_mdot_save_png') && in.plot.beta_mdot_save_png
    fname = 'beta_tt_vs_mdot_preliminary_sweep.png';
    if isfield(in.plot,'beta_mdot_png_name') && ~isempty(in.plot.beta_mdot_png_name)
        fname = in.plot.beta_mdot_png_name;
    end
    try
        saveas(gcf, fname);
        fprintf('Saved beta_tt-mdot diagnostic plot to %s\n', fname);
    catch ME
        warning('Could not save beta_tt-mdot plot: %s', ME.message);
    end
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
