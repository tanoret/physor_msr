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
total_power          = 10.0E+6 # Total reactor Power (W)
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
  library_file       = '../../xs_msre_micro.xml'
  library_name       = 'MSRE-Simplified'
  is_meter           = true
  grid_names         = 'Tfuel'
  grid_variables     = 'T_salt'
  plus               = true
  scalar_flux        = 'sflux_g0   sflux_g1   sflux_g2   sflux_g3   sflux_g4   sflux_g5   sflux_g6   sflux_g7
                        sflux_g8   sflux_g9  sflux_g10  sflux_g11  sflux_g12  sflux_g13  sflux_g14  sflux_g15'
[]
# ================================================================================================================
# TRANSPORT SYSTEM
# ================================================================================================================
[TransportSystems]
  particle                       = neutron
  equation_type                  = eigenvalue
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
  []
[]
# ================================================================================================================
# GEOMETRY AND MESH
# ================================================================================================================
[Mesh]
  [fmg]
    type     = FileMeshGenerator
    file     = '../../mesh_msre_in.e'
  []
  coord_type = 'RZ'
[]
# ================================================================================================================
# AUXVARIABLES AND AUXKERNELS
# ================================================================================================================
[AuxVariables]
  [Peaking_factor]
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
  [nuFission_RR]
    family            = MONOMIAL
    order             = CONSTANT
    block             = ${salt_blocks}
  []
  [Fission_RR]
    family            = MONOMIAL
    order             = CONSTANT
    block             = ${salt_blocks}
  []
  [Absorption_RR]
    family            = MONOMIAL
    order             = CONSTANT
    block             = ${all_blocks}
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
  # Reaction Rates
  # ---------------------------------------------------------------------------------------------
  [nuFission_RR]
    type                = VectorReactionRate
    block               = ${salt_blocks}
    variable            = nuFission_RR
    cross_section       = nu_sigma_fission
    scale_factor        = power_scaling
  []
  [Fission_RR]
    type                = VectorReactionRate
    block               = ${salt_blocks}
    variable            = Fission_RR
    cross_section       = sigma_fission
    scale_factor        = power_scaling
  []
  [Absorption_RR]
    type                = VectorReactionRate
    block               = ${all_blocks} 
    variable            = Absorption_RR
    cross_section       = sigma_absorption
    scale_factor        = power_scaling
  []
  [Peaking_factor]
    type                = NormalizationAux
    variable            = Peaking_factor
    source_variable     = power_density
    normalization       = power_avg
    execute_on          = 'initial timestep_end'
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
  [power_peak]
    type             = ElementExtremeValue
    variable         = Peaking_factor
    value_type       = max
    block            = ${salt_blocks}
    execute_on       = 'initial timestep_end'
  []
  [core_vol]
    type             = VolumePostprocessor
    block            = 'core'
    execute_on       = 'initial timestep_end'
  []
  [RR_Production]
    type             = ElementIntegralVariablePostprocessor
    variable         = nuFission_RR
    block            = ${salt_blocks}
    execute_on       = 'initial timestep_end'
  []
  [RR_Fission]
    type             = ElementIntegralVariablePostprocessor
    variable         = Fission_RR
    block            = ${salt_blocks}
    execute_on       = 'initial timestep_end'
  []
  [RR_Absorption]
    type             = ElementIntegralVariablePostprocessor
    variable         = Absorption_RR
    block            = ${all_blocks} 
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
  [Leakage]
    type             = PartialSurfaceCurrent
    boundary         = 'right bottom top'
    transport_system = transport
    execute_on       = 'initial timestep_end'
  []
[]
# ================================================================================================================
# EXECUTION PARAMETERS
# ================================================================================================================
[Preconditioning]
  [SMP]
    type              = SMP
    full              = true
  []
[]
[Executioner]
  type                = Eigenvalue
  solve_type          = PJFNKMO
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 50'
  l_max_its           = 200
  nl_abs_tol          = 1e-6
[]
# ================================================================================================================
# OUTPUTS & DEBUG
# ================================================================================================================
[Debug]
  show_var_residual_norms       = false
[]
[Outputs]                       
  file_base                     = msre_ss_nodnp_out
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
  []
  [auxvar_solution]
    type             = SolutionVectorFile
    var              = 'T_salt  T_solid   ad_C12 ad_U235 ad_U238  ad_Be9   ad_Li7  ad_F9  ad_Zr90  ad_Zr91  ad_Zr92  ad_Zr94   ad_Zr96'
    loading_var      = 'T_salt  T_solid   ad_C12 ad_U235 ad_U238  ad_Be9   ad_Li7  ad_F9  ad_Zr90  ad_Zr91  ad_Zr92  ad_Zr94   ad_Zr96'
    writing          = false
    execute_on       = 'INITIAL'
  []
[]
