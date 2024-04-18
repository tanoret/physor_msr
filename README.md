# Repository for MSR Physor Workshop 2024
Repository for MSR training workshop at Physor 2024

# Running instructions
Copy the project directory into your local directory: `cp -r /projects/physor_molten_salt_training_2024/ .`
Load the environment variable for the binaries:
1. Clean your modules: `module purge`
2. Load ParaView for post-processing: `moduleload paraview`
3. Loaf griffin-opt and pronghorn-opt: `source /projects/physor_molten_salt_training_2024/env.sh`

# Example run command
`griffin-opt –i mesh_msre.i --mesh-only`
