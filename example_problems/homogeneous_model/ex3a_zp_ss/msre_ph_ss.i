# ================================================================================================================
# Model description
# Molten Salt Reactor Experiment (MSRE) Model - Steady-State Model
# Primary Loop Thermal Hydraulics Model
# Integrates:
# - Porous media model for reactor primary loop
# - Weakly compressible, turbulent flow formulation
# MSRE: reference plant design based on 10.0 MW of MSRE Experiment.
# ================================================================================================================
# Author(s): Dr. Mauricio Tano & Dr. Mustafa K. Jaradat
# ================================================================================================================
# MODEL PARAMETERS
# ================================================================================================================
# Geometry
# ----------------------------------------------------------------------------------------------------------------
core_radius           = 0.69793684
graph_heat_frac       = 0.05467
# ----------------------------------------------------------------------------------------------------------------
# Material Thermal Properties
# ----------------------------------------------------------------------------------------------------------------
cp_graph              = 1757.3              # (J/(kg.K)) specific heat of graphite
rho_graph             = 1860.0              # (kg/(m3)) density of graphite
k_graph               = 40.1                # (W/(m.k)) density of graphite
# ----------------------------------------------------------------------------------------------------------------
cp_steel              = 500.0               # (J/(kg.K)) specific heat of steel
rho_steel             = 8000.0              # (kg/(m3)) density of steel
k_steel               = 15.0                # (W/(m.k)) density of steel
# ----------------------------------------------------------------------------------------------------------------
# Porosity 
# ----------------------------------------------------------------------------------------------------------------
core_porosity         = 0.222831853         # core porosity salt VF=0.222831853, Graphite VF=0.777168147
down_comer_porosity   = 1.0                 # downcomer porosity
lower_plenum_porosity = 0.5                 # lower pelnum porosity
upper_plenum_porosity = 1.0                 # upper pelnum porosity
riser_porosity        = 1.0                 # riser porosity
pump_porosity         = 1.0                 # pump porosity
elbow_porosity        = 1.0                 # elbow porosity
# ----------------------------------------------------------------------------------------------------------------
# Thermal-Hydraulic diameters
# ----------------------------------------------------------------------------------------------------------------
D_H_fuel_channel      = 0.0191334114              # Hydraulic diameter of bypass
D_H_downcomer         = 0.045589414               # Hydraulic diameter of riser
D_H_pipe              = '${fparse 5*0.0254}'      # Riser Hydraulic Diameter
D_H_plena             = '${fparse 2*core_radius}' # Hydraulic diameter of riser
# ----------------------------------------------------------------------------------------------------------------
# Operational Parameters
# ----------------------------------------------------------------------------------------------------------------
T_inlet_hx            = 908.15              # Salt inlet temperature (K)
bulk_htc              = 20000.0             # (W/(m3.K)) core bulk volumetric heat exchange coefficient (already callibrated)
p_outlet              = 1.50653E+05         # 1.01325e+05 # Reactor outlet pressure (Pa)
T_Salt_initial        = 908.15              # inital salt temperature (will change in steady-state)
pump_force            = 5.5E+06            # pump force functor (set to get a loop circulation time of ~25 seconds)
vol_hx                = 1.0E+10             # (W/(m3.K)) volumetric heat exchange coefficient for heat exchanger
                                            # Note: vol_hx need to be tuned to match intermediate HX performance for transients
# ----------------------------------------------------------------------------------------------------------------
# Delayed Neutron Data
# ----------------------------------------------------------------------------------------------------------------
Sc_t                  = 100.0               # turbulent Schmidt number
lambda_1              = 0.013336            # 0.0133104
lambda_2              = 0.032738999         # 0.0305427
lambda_3              = 0.120779999         # 0.115179
lambda_4              = 0.302780002         # 0.301152
lambda_5              = 0.849489987         # 0.879376
lambda_6              = 2.85299993          # 2.91303
beta_1                = 0.000024446         # 8.42817e-05
beta_2                = 0.001257016         # 0.000684616
beta_3                = 0.001200046         # 0.000479796
beta_4                = 0.002674925         # 0.00103883
beta_5                = 0.001102124         # 0.000549185
beta_6                = 0.000460964         # 0.000184087
# ----------------------------------------------------------------------------------------------------------------
fluid_blocks          = 'core lower_plenum upper_plenum down_comer riser pump elbow'
solid_blocks          = 'core core_barrel'
non_solid_blocks      = 'lower_plenum upper_plenum down_comer riser pump elbow'
# ================================================================================================================
# GLOBAL PARAMETERS
# ================================================================================================================
[GlobalParams]
  fp                     = fluid_properties_obj
  porosity               = 'porosity'
  rhie_chow_user_object  = 'pins_rhie_chow_interpolator'
  u                      = superficial_vel_x
  v                      = superficial_vel_y
  advected_interp_method = 'upwind'
  velocity_interp_method = 'rc'
  mixing_length          = 'mixing_length'
[]
# ================================================================================================================
# GEOMETRY AND MESH
# ================================================================================================================
[Mesh]
  coord_type             = 'RZ'
  [fmg]
    type = FileMeshGenerator
    file                 = '../mesh_msre_in.e'
  []
[]
[Problem]
  kernel_coverage_check = false
[]
# ================================================================================================================
# FV VARIABLES
# ================================================================================================================
[Variables]
  [superficial_vel_x]
    type              = PINSFVSuperficialVelocityVariable
    initial_condition = 1e-8
    block             = ${fluid_blocks}
  []
  [superficial_vel_y]
    type              = PINSFVSuperficialVelocityVariable
    initial_condition = 1e-8
    block             = ${fluid_blocks}
  []
  [pressure]
    type              = INSFVPressureVariable
    initial_condition = ${p_outlet}
    block             = ${fluid_blocks}
  []
  [T_fluid]
    type              = INSFVEnergyVariable
    initial_condition = ${T_Salt_initial}
    block             = ${fluid_blocks}
  []
  [T_solid]
    type              = INSFVEnergyVariable
    initial_condition = ${T_Salt_initial}
    block             = ${solid_blocks}
  []
  [c1]
    type              = MooseVariableFVReal
    block             = ${fluid_blocks}
  []
  [c2]
    type              = MooseVariableFVReal
    block             = ${fluid_blocks}
  []
  [c3]
    type              = MooseVariableFVReal
    block             = ${fluid_blocks}
  []
  [c4]
    type              = MooseVariableFVReal
    block             = ${fluid_blocks}
  []
  [c5]
    type              = MooseVariableFVReal
    block             = ${fluid_blocks}
  []
  [c6]
    type              = MooseVariableFVReal
    block             = ${fluid_blocks}
  []
[]
# ================================================================================================================
# THERMAL-HYDRAULICS PROBLEM SETUP
# ================================================================================================================
[FluidProperties]
  [fluid_properties_obj]
    type                             = SimpleFluidProperties
    density0                         = 2705.8554   # kg/m^3
    thermal_expansion                = 0.000177319 # K^{-1}
    cp                               = 1868.0      # J/kg·K
    viscosity                        = 0.008268    # Pa-s11
    thermal_conductivity             = 1.4         # W/m·K
  []
[]
[Modules]
  [NavierStokesFV]
    # Basic settings - weakly-compressible, turbulent flow with buoyancy
    block                            = ${fluid_blocks}
    compressibility                  = 'weakly-compressible'
    porous_medium_treatment          = true
    add_energy_equation              = true
    gravity                          = '0.0 -9.81 0.0'
    
    # Variable naming                
    velocity_variable                = 'superficial_vel_x superficial_vel_y'
    pressure_variable                = 'pressure'
    fluid_temperature_variable       = 'T_fluid'
    
    # Numerical schemes
    pressure_face_interpolation      = average
    momentum_advection_interpolation = upwind
    mass_advection_interpolation     = upwind
    energy_advection_interpolation   = upwind
    velocity_interpolation           = rc
    
    # Porous & Friction treatement
    use_friction_correction          = true
    friction_types                   = 'darcy forchheimer'
    friction_coeffs                  = 'Darcy_coefficient Forchheimer_coefficient'
    consistent_scaling               = 100.0
    porosity_smoothing_layers        = 2
    turbulence_handling              = 'mixing-length'

    # fluid properties
    density                          = 'rho'
    dynamic_viscosity                = 'mu'
    thermal_conductivity             = 'kappa'
    specific_heat                    = 'cp'

    # Energy source-sink
    external_heat_source             = 'power_density_fuel'

    # Boundary Conditions
    wall_boundaries                  = 'left      top      bottom   right    loop_boundary '
    momentum_wall_types              = 'symmetry  slip     noslip   noslip   noslip'
    energy_wall_types                = 'heatflux  heatflux heatflux heatflux heatflux'
    energy_wall_function             = '0         0        0        0        0'

    # Constrain Pressure
    pin_pressure                     = true
    pinned_pressure_value            = ${p_outlet}
    pinned_pressure_point            = '0.0 2.13859 0.0'
    pinned_pressure_type             = point-value-uo

    # Passive Scalar -- solved separetely to integrate porosity jumps
    add_scalar_equation              = false

    #Scaling -- used mainly for nonlinear solves
    momentum_scaling                 = 1e-3
    mass_scaling                     = 10
  []
[]
[FVKernels]
  [energy_storage]
    type                  = PINSFVEnergyTimeDerivative
    variable              = T_solid
    rho                   = rho_s
    cp                    = cp_s
    is_solid              = true
  []
  [solid_energy_diffusion_core]
    type                  = PINSFVEnergyAnisotropicDiffusion
    variable              = T_solid
    kappa                 = 'effective_thermal_conductivity'
    effective_diffusivity = true
    porosity              = 1
  []
  [heat_source]
    type                  = FVCoupledForce
    variable              = T_solid
    v                     = power_density_graph
    block                 = 'core'
  []
  # ----------------------------------------------------------------------------------------------------------------
  [pump_x]
    type                  = INSFVBodyForce
    variable              = superficial_vel_x
    functor               = ${pump_force}
    block                 = 'pump'
    momentum_component    = 'x'
    rhie_chow_user_object = 'pins_rhie_chow_interpolator'
  []
  [pump_y]
    type                  = INSFVBodyForce
    variable              = superficial_vel_y
    functor               = ${pump_force}
    block                 = 'pump'
    momentum_component    = 'y'
    rhie_chow_user_object = 'pins_rhie_chow_interpolator'
  []
  [convection_fluid_hx]
    type                  = NSFVEnergyAmbientConvection
    variable              = T_fluid
    T_ambient             = ${T_inlet_hx}
    alpha                 = ${vol_hx}
    block                 = 'pump'
  []
  # ----------------------------------------------------------------------------------------------------------------
  [convection_core]
    type                  = PINSFVEnergyAmbientConvection
    variable              = T_solid
    T_fluid               = T_fluid
    T_solid               = T_solid
    is_solid              = true
    h_solid_fluid         = ${bulk_htc}
    block                 = 'core'
  []
  [convection_core_completmeent]
    type                  = PINSFVEnergyAmbientConvection
    variable              = T_fluid
    T_fluid               = T_fluid
    T_solid               = T_solid
    is_solid              = false
    h_solid_fluid         = ${bulk_htc}
    block                 = 'core'
  []
  # ----------------------------------------------------------------------------------------------------------------
  # Kernels for solve of delayed neutron precursor transport
  # ----------------------------------------------------------------------------------------------------------------
  [c1_time]
    type                 = FVFunctorTimeKernel
    variable             = 'c1'
  []
  [c2_time]
    type                 = FVFunctorTimeKernel
    variable             = 'c2'
  []
  [c3_time]
    type                 = FVFunctorTimeKernel
    variable             = 'c3'
  []
  [c4_time]
    type                 = FVFunctorTimeKernel
    variable             = 'c4'
  []
  [c5_time]
    type                 = FVFunctorTimeKernel
    variable             = 'c5'
  []
  [c6_time]
    type                 = FVFunctorTimeKernel
    variable             = 'c6'
  []
  [c1_advection]
    type                 = PINSFVMassAdvection
    variable             =  c1
    rho                  = 'c1_porous'
    block                = ${fluid_blocks}
  []
  [c2_advection]
    type                 = PINSFVMassAdvection
    variable             =  c2
    rho                  = 'c2_porous'
    block                = ${fluid_blocks}
  []
  [c3_advection]
    type                 = PINSFVMassAdvection
    variable             =  c3
    rho                  = 'c3_porous'
    block                = ${fluid_blocks}
  []
  [c4_advection]
    type                 = PINSFVMassAdvection
    variable             =  c4
    rho                  = 'c4_porous'
    block                = ${fluid_blocks}
  []
  [c5_advection]
    type                 = PINSFVMassAdvection
    variable             =  c5
    rho                  = 'c5_porous'
    block                = ${fluid_blocks}
  []
  [c6_advection]
    type                 = PINSFVMassAdvection
    variable             =  c6
    rho                  = 'c6_porous'
    block                = ${fluid_blocks}
  []
  [c1_turb_diffusion]
    type                 = INSFVMixingLengthScalarDiffusion
    schmidt_number       = ${Sc_t}
    variable             = c1
    block                = ${fluid_blocks}
  []
  [c2_turb_diffusion]
    type                 = INSFVMixingLengthScalarDiffusion
    schmidt_number       = ${Sc_t}
    variable             = c2
    block                = ${fluid_blocks}
  []
  [c3_turb_diffusion]
    type                 = INSFVMixingLengthScalarDiffusion
    schmidt_number       = ${Sc_t}
    variable             = c3
    block                = ${fluid_blocks}
  []
  [c4_turb_diffusion]
    type                 = INSFVMixingLengthScalarDiffusion
    schmidt_number       = ${Sc_t}
    variable             = c4
    block                = ${fluid_blocks}
  []
  [c5_turb_diffusion]
    type                 = INSFVMixingLengthScalarDiffusion
    schmidt_number       = ${Sc_t}
    variable             = c5
    block                = ${fluid_blocks}
  []
  [c6_turb_diffusion]
    type                 = INSFVMixingLengthScalarDiffusion
    schmidt_number       = ${Sc_t}
    variable             = c6
    block                = ${fluid_blocks}
  []
  [c1_src]
    type                 = FVCoupledForce
    variable             = c1
    v                    = fission_source
    coef                 = ${beta_1}
    block                = ${fluid_blocks}
  []
  [c2_src]
    type                 = FVCoupledForce
    variable             = c2
    v                    = fission_source
    coef                 = ${beta_2}
    block                = ${fluid_blocks}
  []
  [c3_src]
    type                 = FVCoupledForce
    variable             = c3
    v                    = fission_source
    coef                 = ${beta_3}
    block                = ${fluid_blocks}
  []
  [c4_src]
    type                 = FVCoupledForce
    variable             = c4
    v                    = fission_source
    coef                 = ${beta_4}
    block                = ${fluid_blocks}
  []
  [c5_src]
    type                 = FVCoupledForce
    variable             = c5
    v                    = fission_source
    coef                 = ${beta_5}
    block                = ${fluid_blocks}
  []
  [c6_src]
    type                 = FVCoupledForce
    variable             = c6
    v                    = fission_source
    coef                 = ${beta_6}
    block                = ${fluid_blocks}
  []
  [c1_decay]
    type                 = FVReaction
    variable             = c1
    rate                 = ${lambda_1}
    block                = ${fluid_blocks}
  []
  [c2_decay]
    type                 = FVReaction
    variable             = c2
    rate                 = ${lambda_2}
    block                = ${fluid_blocks}
  []
  [c3_decay]
    type                 = FVReaction
    variable             = c3
    rate                 = ${lambda_3}
    block                = ${fluid_blocks}
  []
  [c4_decay]
    type                 = FVReaction
    variable             = c4
    rate                 = ${lambda_4}
    block                = ${fluid_blocks}
  []
  [c5_decay]
    type                 = FVReaction
    variable             = c5
    rate                 = ${lambda_5}
    block                = ${fluid_blocks}
  []
  [c6_decay]
    type                 = FVReaction
    variable             = c6
    rate                 = ${lambda_6}
    block                = ${fluid_blocks}
  []
[]
[FVInterfaceKernels]
  # Conjugated heat transfer with core barrel
  [convection]
    type                = FVConvectionCorrelationInterface
    variable1           = T_fluid
    variable2           = T_solid
    boundary            = 'core_barrel'
    h                   = ${bulk_htc}
    T_solid             = T_solid
    T_fluid             = T_fluid
    subdomain1          = 'core down_comer lower_plenum upper_plenum'
    subdomain2          = 'core_barrel'
    wall_cell_is_bulk   = true
  []
[]
# ================================================================================================================
# AUXVARIABLES & AUXKERNELS
# ================================================================================================================
[AuxVariables]
  [power_density]
    type              = MooseVariableFVReal
    initial_condition = 3.049766e+06
  []
  [power_density_fuel]
    type              = MooseVariableFVReal
    initial_condition = 0.0
  []
  [power_density_graph]
    type              = MooseVariableFVReal
    initial_condition = 0.0
  []
  [fission_source]
    type              = MooseVariableFVReal
    initial_condition = 1.0
  []
  [porosity_var]
    type              = MooseVariableFVReal
    block             = ${fluid_blocks}
  []
  [rho_var]
    type              = MooseVariableFVReal
    initial_condition = 1.0
    block             = ${fluid_blocks}
  []
[]
[AuxKernels]
  [porosity_var_aux]
    type                = FunctorAux
    variable            = porosity_var
    functor             = 'porosity'
    block               = ${fluid_blocks}
  []
  [rho_var_aux]
    type                = FunctorAux
    variable            = 'rho_var'
    functor             = 'rho'
    block               = ${fluid_blocks}
  []
  [fuel_power_density_core]
    type                = ParsedAux
    variable            = power_density_fuel
    coupled_variables   = 'power_density'
    expression          = 'power_density * (1.0-${graph_heat_frac})'
    execute_on          = 'INITIAL timestep_end'
    block               = 'core'
  []
  [fuel_power_density_others]
    type                = ParsedAux
    variable            = power_density_fuel
    coupled_variables   = 'power_density'
    expression          = 'power_density * 1.0'
    execute_on          = 'INITIAL timestep_end'
    block               = ${non_solid_blocks}
  []
  [graph_power_density]
    type                = ParsedAux
    variable            = power_density_graph
    coupled_variables   = 'power_density'
    expression          = 'power_density * ${graph_heat_frac}'
    execute_on          = 'INITIAL timestep_end'
    block               = 'core'
  []
[]
# ================================================================================================================
# MATERIALS
# ================================================================================================================
[FunctorMaterials]
  # ----------------------------------------------------------------------------------------------------------------
  # Setting up material porosities at fluid blocks
  # ----------------------------------------------------------------------------------------------------------------
  [porosity]
    type                    = ADPiecewiseByBlockFunctorMaterial
    prop_name               = 'porosity'
    subdomain_to_prop_value = 'core             ${core_porosity}
                               lower_plenum     ${lower_plenum_porosity}
                               upper_plenum     ${upper_plenum_porosity}
                               down_comer       ${down_comer_porosity}
                               riser            ${riser_porosity}
                               pump             ${pump_porosity}
                               elbow            ${elbow_porosity}
                               core_barrel      0'
  []
  # ----------------------------------------------------------------------------------------------------------------
  # Setting up hydraulic diameters at fluid blocks
  # ----------------------------------------------------------------------------------------------------------------
  [hydraulic_diameter]
    type                    = PiecewiseByBlockFunctorMaterial
    prop_name               = 'characteristic_length'
    subdomain_to_prop_value = 'core             ${D_H_fuel_channel}
                               lower_plenum     ${D_H_plena}
                               upper_plenum     ${D_H_plena}
                               down_comer       ${D_H_downcomer}
                               riser            ${D_H_pipe}
                               pump             ${D_H_pipe}
                               elbow            ${D_H_pipe}'
    block                   = ${fluid_blocks}
  []
  # ----------------------------------------------------------------------------------------------------------------
  # Setting up Fluid & Solid properties 
  # ----------------------------------------------------------------------------------------------------------------
  [fluid_props_to_mat_props]
    type                    = GeneralFunctorFluidProps
    pressure                = 'pressure'
    T_fluid                 = 'T_fluid'
    speed                   = 'speed'
    characteristic_length   = characteristic_length
    block                   = ${fluid_blocks}
  []
  [core_moderator]
    type                    = ADGenericFunctorMaterial
    prop_names              = 'rho_s   cp_s   k_s'
    prop_values             = '${rho_graph} ${cp_graph} ${k_graph}'
    block                   = 'core'
  []
  [core_barrel_steel]
    type                    = ADGenericFunctorMaterial
    prop_names              = 'rho_s   cp_s   k_s'
    prop_values             = '${rho_steel} ${cp_steel} ${k_steel}'
    block                   = 'core_barrel'
  []
  [effective_fluid_thermal_conductivity]
    type                    = ADGenericVectorFunctorMaterial
    prop_names              = 'kappa'
    prop_values             = 'k k k'
    block                   = ${fluid_blocks}
  []
  [effective_solid_thermal_conductivity]
    type                    = ADGenericVectorFunctorMaterial
    prop_names              = 'effective_thermal_conductivity'
    prop_values             = 'k_s k_s k_s'
    block                   =  ${solid_blocks}
  []
  # Drag correlations per block
  [isotropic_drag_core]
    type                    = FunctorChurchillDragCoefficients
    multipliers             = '100000 100 100000'
    block                   = 'core'
  []
  [drag_lower_plenum]
    type                    = FunctorChurchillDragCoefficients
    multipliers             = '10 1 10'
    block                   = 'upper_plenum'
  []
  [drag_upper_plenum]
    type                    = FunctorChurchillDragCoefficients
    multipliers             = '1 1 1'
    block                   = 'lower_plenum'
  []
  [drag_downcomer]
    type                    = FunctorChurchillDragCoefficients
    multipliers             = '1 1 1'
    block                   = 'down_comer'
  []
  [drag_piping]
    type                    = FunctorChurchillDragCoefficients
    multipliers             = '0 0 0'
    block                   = 'riser pump elbow'
  []
  # ----------------------------------------------------------------------------------------------------------------
  # Materials for computing corrected DNP advection
  # ----------------------------------------------------------------------------------------------------------------
  [c1_mat]
    type                    = ADParsedFunctorMaterial
    expression              = 'c1 / porosity'
    functor_names           = 'c1 porosity'
    functor_symbols         = 'c1 porosity'
    property_name           = 'c1_porous'
  []
  [c2_mat]
    type                    = ADParsedFunctorMaterial
    expression              = 'c2 / porosity'
    functor_names           = 'c2 porosity'
    functor_symbols         = 'c2 porosity'
    property_name           = 'c2_porous'
  []
  [c3_mat]
    type                    = ADParsedFunctorMaterial
    expression              = 'c3 / porosity'
    functor_names           = 'c3 porosity'
    functor_symbols         = 'c3 porosity'
    property_name           = 'c3_porous'
  []
  [c4_mat]
    type                    = ADParsedFunctorMaterial
    expression              = 'c4 / porosity'
    functor_names           = 'c4 porosity'
    functor_symbols         = 'c4 porosity'
    property_name           = 'c4_porous'
  []
  [c5_mat]
    type                    = ADParsedFunctorMaterial
    expression              = 'c5 / porosity'
    functor_names           = 'c5 porosity'
    functor_symbols         = 'c5 porosity'
    property_name           = 'c5_porous'
  []
  [c6_mat]
    type                    = ADParsedFunctorMaterial
    expression              = 'c6 / porosity'
    functor_names           = 'c6 porosity'
    functor_symbols         = 'c6 porosity'
    property_name           = 'c6_porous'
  []
[]
# ================================================================================================================
# POSTPROCESSORS
# ================================================================================================================
[Postprocessors]
  [pressure_outlet]
    type                    = SideAverageValue
    variable                = pressure
    boundary                = 'pump_inlet'
  []
  [pressure_inlet]
    type                    = SideAverageValue
    variable                = 'pressure'
    boundary                = 'downcomer_outlet'
  []
  [pressure_core_delta]
    type                   = ParsedPostprocessor
    expression             = 'pressure_inlet - pressure_outlet'
    pp_names               = 'pressure_inlet pressure_outlet'
    execute_on             = 'initial timestep_end'
  []
  [T_inlet]
    type                    = SideAverageValue
    variable                = 'T_fluid'
    boundary                = 'downcomer_outlet'
  []
  [T_outlet]
    type                    = SideAverageValue
    variable                = 'T_fluid'
    boundary                = 'riser_inlet'
  []
  [T_core_inlet]
    type                    = SideAverageValue
    variable                = 'T_fluid'
    boundary                = 'core_in'
  []
  [T_core_outlet]
    type                    = SideAverageValue
    variable                = 'T_fluid'
    boundary                = 'core_out'
  []
  [v_core_inlet]
    type                    = SideAverageValue
    variable                = 'superficial_vel_y'
    boundary                = 'core_in'
  []
  [v_core_outlet]
    type                    = SideAverageValue
    variable                = 'superficial_vel_y'
    boundary                = 'core_out'
  []
  [T_core_delta]
    type                   = ParsedPostprocessor
    expression             = 'T_core_outlet - T_core_inlet'
    pp_names               = 'T_core_outlet T_core_inlet'
    execute_on             = 'initial timestep_end'
  []
  [area_pp_downcomer_inlet]
    type                    = AreaPostprocessor
    boundary                = 'downcomer_inlet'
    execute_on              = 'INITIAL'
  [] 
  [vfr_downcomer]
    type                    = VolumetricFlowRate
    vel_x                   = superficial_vel_x
    vel_y                   = superficial_vel_y
    advected_quantity       = 1.0
    boundary                = 'downcomer_inlet'
  []
  [vfr_pump]
    type                    = VolumetricFlowRate
    vel_x                   = superficial_vel_x
    vel_y                   = superficial_vel_y
    advected_quantity       = 1.0
    boundary                = 'pump_outlet'
  []
  [mfr_core_inlet]
    type                    = VolumetricFlowRate
    vel_x                   = superficial_vel_x
    vel_y                   = superficial_vel_y
    advected_quantity       = rho
    boundary                = 'downcomer_outlet'
  []
  [mfr_core_outlet]
    type                    = VolumetricFlowRate
    vel_x                   = superficial_vel_x
    vel_y                   = superficial_vel_y
    advected_quantity       = rho
    boundary                = 'pump_inlet'
  []
  [core_vol]
    type                    = VolumePostprocessor
    block                   = 'core'
    execute_on              = 'initial timestep_end'
  []                        
  [loop_vol]                
    type                    = VolumePostprocessor
    block                   = ${non_solid_blocks}
    execute_on              = 'initial timestep_end'
  []                        
  [Tmax_fuel]               
    type                    = ElementExtremeValue
    value_type              = max
    variable                = T_fluid
    block                   = ${fluid_blocks}
    execute_on              = 'initial timestep_end'
  []                        
  [Tavg_fuel]               
    type                    = ElementAverageValue
    variable                = T_fluid
    block                   = ${fluid_blocks}
    execute_on              = 'initial timestep_end'
  []                        
  [Tmax_core_fuel]          
    type                    = ElementExtremeValue
    value_type              = max
    variable                = T_fluid
    block                   = 'core'
    execute_on              = 'initial timestep_end'
  []                        
  [Tavg_core_fuel]          
    type                    = ElementAverageValue
    variable                = T_fluid
    block                   = 'core'
    execute_on              = 'initial timestep_end'
  []                        
  [Tmax_mod]                
    type                    = ElementExtremeValue
    value_type              = max
    variable                = T_solid
    block                   = 'core'
    execute_on              = 'initial timestep_end'
  []                        
  [Tavg_mod]                
    type                    = ElementAverageValue
    variable                = T_solid
    block                   = 'core'
    execute_on              = 'initial timestep_end'
  []                        
  [power_total]             
    type                    = ElementIntegralVariablePostprocessor
    variable                = power_density
    execute_on              = 'initial timestep_end'
  []                        
  [power_avg]               
    type                    = ElementAverageValue
    variable                = power_density
    execute_on              = 'initial timestep_end'
  []
  [power_fuel_total]             
    type                    = ElementIntegralVariablePostprocessor
    variable                = power_density_fuel
    execute_on              = 'initial timestep_end'
  [] 
  [power_ghrap_total]             
    type                    = ElementIntegralVariablePostprocessor
    variable                = power_density_graph
    execute_on              = 'initial timestep_end'
  []
  [power_total_2]
    type                   = ParsedPostprocessor
    expression             = 'power_ghrap_total + power_fuel_total'
    pp_names               = 'power_ghrap_total power_fuel_total'
    execute_on             = 'initial timestep_end'
  []
[]
# ================================================================================================================
# EXECUTION PARAMETERS
# ================================================================================================================
[Executioner]
  type                             = Transient
  solve_type                       = NEWTON
  petsc_options_iname              = '-pc_type -sub_pc_factor_shift_type'
  petsc_options_value              = ' lu       NONZERO'
  automatic_scaling                = true
  nl_abs_tol                       = 1e-6
  nl_max_its                       = 10
  [TimeStepper]
    type                           = IterationAdaptiveDT
    dt                             = 0.1
    optimal_iterations             = 20
    iteration_window               = 2
    growth_factor                  = 2
    cutback_factor                 = 0.5
  []
  end_time                         = 1e10
  steady_state_detection           = true
  steady_state_tolerance           = 1e-16
[]
# ================================================================================================================
# OUTPUTS & DEBUG
# ================================================================================================================
[Debug]
  show_var_residual_norms          = false
[]
[Outputs]
  csv                              = true
  exodus                           = true
  print_linear_converged_reason    = false
  print_linear_residuals           = false
  print_nonlinear_converged_reason = false
[]
