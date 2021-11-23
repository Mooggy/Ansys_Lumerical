
## Package overview:
 - Ansys-Lumerical integration is dependent-package of SiEPIC-Tools, and it offers a layout-first design methodology. This allows the user to directly simulate from the layout, without needing to create the schematic first. 


## Package features:

- Circuit simulations:
  - Add INTERCONNECT IOs (text labels)
  - Netlist generation: Creating a Spice netlist suitable for for circuit simulations. This includes extracting the waveguide length (wg_length) for all waveguides.
  - Menu item "Import circuit to INTC" will automatically: generate the netlist, launch Lumerical INTERCONNECT to perform the circuit simulations. It also allows users to reuse their testbenches by importing their testbench files.
- Environmental setup: 
  - Set up project directory to store temporary gds and netlist files
  - Set up Lumerical installation path to ensure a stable integration with Lumerical simulation softwares


