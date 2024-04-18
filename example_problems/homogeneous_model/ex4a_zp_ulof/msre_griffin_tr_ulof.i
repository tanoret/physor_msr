# ================================================================================================================
# Model description
# Molten Salt Reactor Experiment (MSRE) Model
# Core Neutronics Model
# Integrates:
# - Doppler-Temperature feedback with interpolation from tabulated cross sections
# - Density-Temperature feedback with field functions for density
# - DNPC-Drift feedback
# MSRE: reference plant design based on 10 MW of MSRE Experiment.
# ================================================================================================================
# Author(s): Dr. Mustafa K. Jaradat
# ================================================================================================================
# ================================================================================================================
# MODEL PARAMETERS
# ================================================================================================================
total_power          = 10.0E+2 # Total reactor Power (W)
T_Salt_initial       = 923.0   # K
Salt_Density_initial = 2263.0  # kg/m3
# ================================================================================================================
solid_blocks         = 'core core_barrel'
salt_blocks          = 'core lower_plenum upper_plenum down_comer riser pump elbow'
non_solid_blocks     = 'lower_plenum upper_plenum down_comer riser pump elbow'
all_blocks           = 'core lower_plenum upper_plenum down_comer core_barrel riser pump elbow'
# ================================================================================================================
# GLOBAL PARAMETERS
# ================================================================================================================
[GlobalParams]
  library_file       = '../xs_msre_micro.xml'
  library_name       = 'MSRE-Simplified'
  is_meter           = true
  grid_names         = 'Tfuel'
  grid_variables     = 'T_salt'
  plus               = true
[]
# ================================================================================================================
# TRANSPORT SYSTEM
# ================================================================================================================
[TransportSystems]
  particle                       = neutron
  equation_type                  = transient
  G                              = 16
  ReflectingBoundary             = 'left'
  VacuumBoundary                 = 'bottom right top loop_boundary top_core_barrel'
  [transport]
    scheme                       = CFEM-Diffusion
    family                       = LAGRANGE
    order                        = FIRST
    n_delay_groups               = 6
    assemble_scattering_jacobian = true
    assemble_fission_jacobian    = true
    external_dnp_variable        = 'dnp'
    fission_source_aux           = true
  []
[]
# ================================================================================================================
# GEOMETRY AND MESH
# ================================================================================================================
[Mesh]
  [fmg]
    type     = FileMeshGenerator
    file     = '../mesh_msre_in.e'
  []
  coord_type = 'RZ'
[]
# ================================================================================================================
# AUXVARIABLES AND AUXKERNELS
# ================================================================================================================
[Functions]
  [dt_max_fn]
    type = PiecewiseLinear
    x = '  0   25   100  1000   2000'
    y = '1.0  1.0   5.0    50    500'
  []
[]
[AuxVariables]
  [vel_x]
    family            = MONOMIAL
    order             = CONSTANT
    block             = ${salt_blocks}
  []
  [vel_y]
    family            = MONOMIAL
    order             = CONSTANT
    block             = ${salt_blocks}
  []
  [T_salt]
    family            = MONOMIAL
    order             = CONSTANT
    block             = ${all_blocks}
  []
  [T_solid]
    family            = MONOMIAL
    order             = CONSTANT
    block             = ${solid_blocks} 
  []
  [c1]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [c2]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [c3]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [c4]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [c5]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [c6]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [dnp]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
    components        = 6
  []

  [ad_U235]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_U238]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_Be9]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_Li7]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_F9]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_Zr90]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_Zr91]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_Zr92]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_Zr94]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_Zr96]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${salt_blocks}
  []
  [ad_C12]
    order             = CONSTANT
    family            = MONOMIAL
    block             = ${solid_blocks}
  []
[]

[AuxKernels]
  # ---------------------------------------------------------------------------------------------
  # Dealyed Neutron Source
  # ---------------------------------------------------------------------------------------------
  [build_dnp]
    type                = BuildArrayVariableAux
    variable            = dnp
    component_variables = 'c1 c2 c3 c4 c5 c6'
    execute_on          = 'initial timestep_begin final'
  []
  # ---------------------------------------------------------------------------------------------
  # Updating Non-Solid Regions Atom Densities
  # ---------------------------------------------------------------------------------------------
  [update_ad_f_U235]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_U235
    coupled_variables   = 'T_salt'
    expression          = '9.692763E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_U238]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_U238
    coupled_variables   = 'T_salt'
    expression          = '1.967924E-04 *(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_Be9]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_Be9
    coupled_variables   = 'T_salt'
    expression          = '9.496943E-03*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_Li7]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_Li7
    coupled_variables   = 'T_salt'
    expression          = '2.121311E-02 *(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_F9]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_F9
    coupled_variables   = 'T_salt'
    expression          = '4.790899E-02*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_Zr90]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_Zr90
    coupled_variables   = 'T_salt'
    expression          = '8.395505E-04*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_Zr91]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_Zr91
    coupled_variables   = 'T_salt'
    expression          = '1.830855E-04*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_Zr92]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_Zr92
    coupled_variables   = 'T_salt'
    expression          = '2.798495E-04*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_Zr94]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_Zr94
    coupled_variables   = 'T_salt'
    expression          = '2.836027E-04*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_f_Zr96]
    type                = ParsedAux
    block               = ${non_solid_blocks}
    variable            = ad_Zr96
    coupled_variables   = 'T_salt'
    expression          = '4.568967E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  # ---------------------------------------------------------------------------------------------
  # Updating Core (Homogenized Salt & Moderator) Atom Densities
  # ---------------------------------------------------------------------------------------------
  [update_ad_c_U235]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_U235
    coupled_variables   = 'T_salt'
    expression          = '2.159856E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_U238]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_U238
    coupled_variables   = 'T_salt'
    expression          = '4.385162E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_Be9]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_Be9
    coupled_variables   = 'T_salt'
    expression          = '2.116221E-03*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_Li7]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_Li7
    coupled_variables   = 'T_salt'
    expression          = '4.726958E-03*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_F9]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_F9
    coupled_variables   = 'T_salt'
    expression          = '1.067565E-02*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_Zr90]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_Zr90
    coupled_variables   = 'T_salt'
    expression          = '1.870786E-04*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_Zr91]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_Zr91
    coupled_variables   = 'T_salt'
    expression          = '4.079728E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_Zr92]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_Zr92
    coupled_variables   = 'T_salt'
    expression          = '6.235937E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_c_Zr94]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_Zr94
    coupled_variables   = 'T_salt'
    expression          = '6.319571E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
  [update_ad_Zr96]
    type                = ParsedAux
    block               = 'core'
    variable            = ad_Zr96
    coupled_variables   = 'T_salt'
    expression          = '1.018111E-05*(1.0-0.4798*(T_salt-${T_Salt_initial})/${Salt_Density_initial})'
    execute_on          = 'INITIAL timestep_end'
  []
[]
# ================================================================================================================
# MATERIALS
# ================================================================================================================
[PowerDensity]
  power                   = '${fparse total_power}'
  power_density_variable  = power_density
  family                  = MONOMIAL
  order                   = CONSTANT
[]
[Materials]
  [activeCore]
    type           = CoupledFeedbackNeutronicsMaterial
    isotopes       = '    C12     U235     U238      BE9      LI7     F19
                         ZR90     ZR91     ZR92     ZR94     ZR96'
    densities      = ' ad_C12  ad_U235  ad_U238   ad_Be9   ad_Li7   ad_F9
                      ad_Zr90  ad_Zr91  ad_Zr92  ad_Zr94  ad_Zr96'
    material_id    = 1
    block          = 'core'
  []
  [flow_loop]
    type           = CoupledFeedbackNeutronicsMaterial
    isotopes       = '   U235     U238      BE9      LI7     F19
                         ZR90     ZR91     ZR92     ZR94     ZR96'
    densities      = 'ad_U235  ad_U238   ad_Be9   ad_Li7   ad_F9
                      ad_Zr90  ad_Zr91  ad_Zr92  ad_Zr94  ad_Zr96'
    material_id    = 2
    block          = ${non_solid_blocks}
  []
  [core_barrel]
    type           = CoupledFeedbackNeutronicsMaterial
    isotopes       = '   C12'
    densities      = 'ad_C12'
    material_id    = 1
    block          = 'core_barrel'
  []
[]
# ================================================================================================================
# POSTPROCESSORS
# ================================================================================================================
[Postprocessors]
  [dt]
    type             = TimestepSize
  []
  [dt_max_pp]
    type             = FunctionValuePostprocessor
    function         = dt_max_fn
    execute_on       = 'INITIAL TIMESTEP_END'
  []
  [power_total]
    type             = ElementIntegralVariablePostprocessor
    block            = ${salt_blocks}
    variable         = power_density
    execute_on       = 'initial timestep_end'
  []
  [power_avg]
    type             = ElementAverageValue
    variable         = power_density
    block            = ${salt_blocks}
    execute_on       = 'initial timestep_end'
  []
  [core_vol]
    type             = VolumePostprocessor
    block            = 'core'
    execute_on       = 'initial timestep_end'
  []
  [Tmax_fuel]
    type             = ElementExtremeValue
    value_type       = max
    variable         = T_salt
    block            = ${salt_blocks}
    execute_on       = 'initial timestep_end'
  []
  [Tavg_fuel]
    type             = ElementAverageValue
    variable         = T_salt
    block            = ${salt_blocks}
    execute_on       = 'initial timestep_end'
  []
  [Tmax_core_fuel]
    type             = ElementExtremeValue
    value_type       = max
    variable         = T_salt
    block            = 'core'
    execute_on       = 'initial timestep_end'
  []
  [Tavg_core_fuel]
    type             = ElementAverageValue
    variable         = T_salt
    block            = 'core'
    execute_on       = 'initial timestep_end'
  []
  [Tmax_mod]
    type             = ElementExtremeValue
    value_type       = max
    variable         = T_solid
    block            = 'core'
    execute_on       = 'initial timestep_end'
  []
  [Tavg_mod]
    type             = ElementAverageValue
    variable         = T_solid
    block            = 'core'
    execute_on       = 'initial timestep_end'
  []
[]
# ================================================================================================================
# EXECUTION PARAMETERS
# ================================================================================================================
[Preconditioning]
  [SMP]
    type               = SMP
    full               = true
  []
[]
[Executioner]
  type                 = Transient
  solve_type           = PJFNK
  petsc_options_iname  = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value  = 'hypre boomeramg 50'
  
  start_time = 0.0
  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 1.0
    timestep_limiting_postprocessor = dt_max_pp
    optimal_iterations = 20
    iteration_window   = 2
    growth_factor      = 2
    cutback_factor     = 0.5
  []
  end_time             = 2000.0
  auto_advance         = true

  l_max_its            = 200
  nl_abs_tol           = 1e-5
  fixed_point_min_its  = 1
  fixed_point_max_its  = 50
  fixed_point_rel_tol  = 1e-5
  fixed_point_abs_tol  = 1e-5
[]
# ================================================================================================================
# MULTIAPPS AND TRANSFERS
# ================================================================================================================
[MultiApps]
  [flow_dnp]
    type                         = TransientMultiApp
    input_files                  = 'msre_ph_tr_ulof.i'
    execute_on                   = 'timestep_end'
    max_procs_per_app            = 48
    keep_solution_during_restore = true
    catch_up                     = true
    app_type                     = PronghornApp
    library_name                 = ''
  []
[]
[Transfers]
  [power_density]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app           = flow_dnp
    source_variable        = 'power_density'
    variable               = 'power_density'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [fission_source]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app           = flow_dnp
    source_variable        = 'fission_source'
    variable               = 'fission_source'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [c1]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'c1'
    variable               = 'c1'
    execute_on             = 'timestep_end'
    search_value_conflicts = false

  []
  [c2]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'c2'
    variable               = 'c2'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [c3]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'c3'
    variable               = 'c3'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [c4]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'c4'
    variable               = 'c4'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [c5]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'c5'
    variable               = 'c5'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [c6]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'c6'
    variable               = 'c6'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [T_salt]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'T_fluid'
    variable               = 'T_salt'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [T_graph]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'T_solid'
    variable               = 'T_solid'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [vel_x]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'superficial_vel_x'
    variable               = 'vel_x'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
  [vel_y]
    type                   = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app         = flow_dnp
    source_variable        = 'superficial_vel_y'
    variable               = 'vel_y'
    execute_on             = 'timestep_end'
    search_value_conflicts = false
  []
[]
# ================================================================================================================
# OUTPUTS & DEBUG
# ================================================================================================================
[Debug]
  show_var_residual_norms       = false
[]
[Outputs]
  file_base                     = msre_tr_ulof_out
  csv                           = true
  exodus                        = true
  perf_graph                    = true
  print_linear_converged_reason = false
  print_linear_residuals        = false
  execute_on                    = 'INITIAL FINAL TIMESTEP_END'
[]
# ================================================================================================================
# USER OBJECTS FOR RESTART
# ================================================================================================================
[UserObjects]
  [transport_solution]
    type             = TransportSolutionVectorFile
    transport_system = transport
    writing          = false
    execute_on       = 'INITIAL'
    scale_with_keff  = false
    folder           = '../ex3a_zp_ss/gold/'
  []
  [auxvar_solution]
    type             = SolutionVectorFile
    var              = '  vel_x   vel_y  T_salt  T_solid     c1     c2     c3     c4     c5     c6     ad_C12
                        ad_U235 ad_U238  ad_Be9   ad_Li7  ad_F9  ad_Zr90  ad_Zr91  ad_Zr92  ad_Zr94   ad_Zr96'
    loading_var      = '  vel_x   vel_y  T_salt  T_solid     c1     c2     c3     c4     c5     c6     ad_C12
                        ad_U235 ad_U238  ad_Be9   ad_Li7  ad_F9  ad_Zr90  ad_Zr91  ad_Zr92  ad_Zr94   ad_Zr96'
    writing          = false
    execute_on       = 'INITIAL'
    folder           = '../ex3a_zp_ss/gold/'
  []
[]