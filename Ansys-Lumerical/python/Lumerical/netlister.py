# Copyright (c) 2021, Ansys, Inc. All rights reserved.


#This package is adopted from circuit simulation part of SiEPIC Tools
#This contains most functions needed for generating netlist. Some core function (e,g. find_pins(),find_nets()) are from SiEPIC Tools.

import pya


def has_duplicates(value):
    if len(value) != len(set(value)):
        return True
    else:
        return False

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def trim_netlist_ansys(nets, components, selected_component, verbose=None):
    selected = selected_component
    #>17        <2
    # nets[0].pins[0].component.idx
    trimmed_net = []
    net_idx = [[each.pins[0].component.idx, each.pins[1].component.idx] for each in nets]
    len_net_idx = len(net_idx)
    count = 0
    trimmed_nets, trimmed_components = [], []
    while count < (len_net_idx - 1):
        for i in range(count + 1, len_net_idx):  # i keep track of nets from next net to last net
            # first set is formed of elements from current to backwards
            first_set = set(net_idx[count])
            # second set is formed of elements from current + 1 to forward
            second_set = set(net_idx[i])
            if len(first_set.intersection(second_set)) > 0:  # if there are common elements between two sets
                net_idx.pop(i)  # remove the nets from the list
                net_idx.pop(count)  # remove the count net as well
                # merged them and add to the list
                net_idx.append(list(first_set.union(second_set)))
                len_net_idx -= 1  # 2 removed 1 added so reduce 1
                count -= 1  # readjust count as the elements have shifted to left
                break
        count += 1
    for net in net_idx:
        if(selected.idx in net):
            trimmed_components = [each for each in components if each.idx in net]
            trimmed_nets = [each for each in nets if (
                each.pins[0].component.idx in net or each.pins[1].component.idx in net)]
            if verbose:
                print("success - netlist trimmed")

    return trimmed_nets, trimmed_components

class port_info:

    def __init__(self, detector_net, port_name):
        self.detector_net = detector_net
        self.port_name = port_name
        
def get_LumericalINTERCONNECT_analyzers_ansys(self, components, verbose=None):
    """
    Find - LumericalINTERCONNECT_Port
    get their parameters
    determine which OpticalIO they are connected to, and find their nets
    Assume that the detectors and laser are on the topcell (not subcells); don't perform transformations.

    returns: parameters, nets in order

    usage:
    laser_net, detector_nets, wavelength_start, wavelength_stop, wavelength_points, ignoreOpticalIOs = get_LumericalINTERCONNECT_analyzers(topcell, components)
    """

    topcell = self

    #from . import _globals
    # Define enumeration for pins
    PIN_TYPES = enum('OPTICALIO', 'OPTICAL', 'ELECTRICAL')
    PIN_LENGTH = 100  # 0.1 micron
    from SiEPIC.utils import select_paths, get_technology
    from SiEPIC.core import Net
    TECHNOLOGY = get_technology()

    layout = topcell.layout()
    LayerLumericalN = self.layout().layer(TECHNOLOGY['Lumerical'])

    # default is the 1st polarization
    orthogonal_identifier = 1

    # Find the laser and detectors in the layout.
    iter1 = topcell.begin_shapes_rec(LayerLumericalN)
    n_IO = 0
    ports_info = []
    port_names = ''
    laser_net = None
    #wavelength_start, wavelength_stop, wavelength_points, orthogonal_identifier, ignoreOpticalIOs = 0, 0, 0, 0, 0
    while not(iter1.at_end()):
        subcell = iter1.cell()             # cell (component) to which this shape belongs
        if iter1.shape().is_box():
            box = iter1.shape().box.transformed(iter1.itrans())
            #changed to search for Lumerical_INTERCONNECT_Port for netlisting
            if iter1.cell().basic_name() == ("Lumerical_INTERCONNECT_Port"):
                n_IO += 1
                # *** todo read parameters from Text labels rather than PCell:
                #change to search for Port name
                port_name = subcell.pcell_parameters_by_name()["name"]
                if verbose:
                    print("%s: Port {%s} %s, box -- %s; %s" %
                          (n_IO, subcell.basic_name(), port_name, box.p1, box.p2))
                #opticalIO is extracted from here
                port_names += ' ' + port_name
                # find components which have an IO pin inside the Lumerical box:
                components_IO = [c for c in components if any(
                    [box.contains(p.center) for p in c.pins if p.type == PIN_TYPES.OPTICALIO])]
                if len(components_IO) > 1:
                    raise Exception("Error - more than 1 optical IO connected to the port.")
                if len(components_IO) == 0:
                    print("Warning - No optical IO connected to the port.")
#          raise Exception("Error - 0 optical IO connected to the port.")
                else:
                    p = [p for p in components_IO[0].pins if p.type == PIN_TYPES.OPTICALIO]
                    # change hard coded 'detector' string to pcell property (port.name)
                    p[0].pin_name = port_name
                    p[0].net = Net(idx=p[0].pin_name, pins=p)
                    ports_info.append(port_info(p[0].net, port_name))
                    if verbose:
                        print(" - pin_name: %s" % (p[0].pin_name))

            if iter1.cell().basic_name() == ("Lumerical_INTERCONNECT_Laser"):
                n_IO += 1
                # *** todo read parameters from Text labels rather than PCell:
                wavelength_start = subcell.pcell_parameters_by_name()["wavelength_start"]
                wavelength_stop = subcell.pcell_parameters_by_name()["wavelength_stop"]
                wavelength_points = subcell.pcell_parameters_by_name()["npoints"]
                orthogonal_identifier = subcell.pcell_parameters_by_name()["orthogonal_identifier"]
                ignoreOpticalIOs = subcell.pcell_parameters_by_name()["ignoreOpticalIOs"]
                if verbose:
                    print("%s: Laser {%s}, box -- %s; %s" %
                          (n_IO, subcell.basic_name(), box.p1, box.p2))
                # find components which have an IO pin inside the Lumerical box:
                components_IO = [c for c in components if any(
                    [box.contains(p.center) for p in c.pins if p.type == PIN_TYPES.OPTICALIO])]
                if len(components_IO) > 1:
                    raise Exception("Error - more than 1 optical IO connected to the laser.")
                if len(components_IO) == 0:
                    print("Warning - No optical IO connected to the laser.")
#          raise Exception("Error - 0 optical IO connected to the laser.")
                else:
                    p = [p for p in components_IO[0].pins if p.type == PIN_TYPES.OPTICALIO]
                    p[0].pin_name += '_laser' + str(n_IO)
                    laser_net = p[0].net = Net(idx=p[0].pin_name, pins=p)
                    if verbose:
                        print(" - pin_name: %s" % (p[0].pin_name))

        iter1.next()

    # Sort the detectors:
    ports_info2 = sorted(ports_info, key=lambda d: d.port_name)

    # output:
    detector_nets = []
    for d in ports_info2:
        detector_nets.append(d.detector_net)

    return laser_net, detector_nets, port_names

def get_LumericalINTERCONNECT_analyzers_from_text_label_ansys(self, components, verbose=True):
    from SiEPIC.utils import get_technology, get_technology_by_name, select_instances, find_automated_measurement_labels, get_layout_variables
    from Lumerical.netlister import identify_nets_ansys
    from SiEPIC.extend import find_components, find_pins
    from SiEPIC.core import Pin
    #TECHNOLOGY, lv, ly, topcell = get_layout_variables()
    PIN_TYPES = enum('OPTICALIO', 'OPTICAL', 'ELECTRICAL')
    PIN_LENGTH = 100  # 0.1 micron

    topcell = self
    layout = topcell.layout()
    TECHNOLOGY = get_technology()

    Text_1 = topcell.layout().layer(TECHNOLOGY['Text'])  

    # find all text labels starting with INTC_IO and opt_in
    # find labels and generate port_names
    port_names = ''
    nets, components = topcell.identify_nets_ansys(verbose = True)

    pins = topcell.find_pins()

    components_sorted = []
    #iterate all labels touching this pin.
    iter2 = topcell.begin_shapes_rec(Text_1)
    while not(iter2.at_end()):
        text = iter2.shape().text
        #if text.string.find("INTC_IO") > -1:
        if text.string.find("INTC_IO") > -1 or text.string.find("opt_in") > -1:  
            print( "\n- -label text: %s" %text.string)
            port_names = port_names + " " + text.string
            pin_path = text.trans
            t = text
            components_sorted = sorted([c for c in components],key=lambda x: x.trans.disp.to_p().distance(pya.Point(t.x, t.y).to_dtype(1)))
            # KLayout issue: waveguide location is always 0,0. So sorting distance won't work for waveguides. Ff the closest component is a waveguide, do the following. 
            # Get a list of all pins in the layout
            all_pin = []
            for c in components:
              all_pin += c.pins
            print("All pins in the layout: %s" %all_pin)
            #for p in all_pin, get the closest pin to the text label. Modify the pin properties 
            pins_wg = sorted([p for p in all_pin],key=lambda x: x.center.distance(pya.Point(t.x, t.y).to_dtype(1))) 
            print("Closest pin of the waveguide ") 
            pins_wg[0].display()
            pins_wg[0].pin_name = text.string
            pins_wg[0].type = 0
            print("connected component:")
            components_sorted[0].display()
        iter2.next()
    print(port_names)


    return port_names, components_sorted


def get_LumericalINTERCONNECT_analyzers_from_opt_in_ansys(self, components, verbose=None, opt_in_selection_text=[]):
    """
    From the opt_in label, find the trimmed circuit, and assign a laser and detectors

    returns: parameters, nets in order

    usage:
    laser_net, detector_nets, wavelength_start, wavelength_stop, wavelength_points, ignoreOpticalIOs, detector_list = get_LumericalINTERCONNECT_analyzers_from_opt_in(topcell, components)
    """
    #from . import _globals
    from SiEPIC.core import Net

    PIN_TYPES = enum('OPTICALIO', 'OPTICAL', 'ELECTRICAL')
    PIN_LENGTH = 100  # 0.1 micron

    from SiEPIC.utils import load_DFT
    DFT = load_DFT()
    if not DFT:
        if verbose:
            print(' no DFT rules available.')
        return False, False, False, False, False, False, False, False

    from SiEPIC.scripts import user_select_opt_in
    opt_in_selection_text, opt_in_dict = user_select_opt_in(
        verbose=verbose, option_all=False, opt_in_selection_text=opt_in_selection_text)
    if not opt_in_dict:
        if verbose:
            print(' no opt_in selected.')
        return False, False, False, False, False, False, False, False

    # find closest GC to opt_in (pick the 1st one... ignore the others)
    t = opt_in_dict[0]['Text']
    components_sorted = sorted([c for c in components if [p for p in c.pins if p.type == PIN_TYPES.OPTICALIO]],
                               key=lambda x: x.trans.disp.to_p().distance(pya.Point(t.x, t.y).to_dtype(1)))
    if not(components_sorted):
        warning = pya.QMessageBox()
        warning.setStandardButtons(pya.QMessageBox.Ok)
        warning.setText("To run a simulation, you need to have optical IO in the layout." )
        pya.QMessageBox_StandardButton(warning.exec_())
        return False, False, False, False, False, False, False, False
        
    dist_optin_c = components_sorted[0].trans.disp.to_p().distance(pya.Point(t.x, t.y).to_dtype(1))
    if verbose:
        print(" - Found opt_in: %s, nearest GC: %s.  Locations: %s, %s. distance: %s" % (opt_in_dict[0][
              'Text'], components_sorted[0].instance,  components_sorted[0].center, pya.Point(t.x, t.y), dist_optin_c))
    if dist_optin_c > float(DFT['design-for-test']['opt_in']['max-distance-to-grating-coupler']) * 1000:
        warning = pya.QMessageBox()
        warning.setStandardButtons(pya.QMessageBox.Ok)
        warning.setText("To run a simulation, you need to have an opt_in label with %s microns from the nearest grating coupler" % int(
            DFT['design-for-test']['opt_in']['max-distance-to-grating-coupler']))
        pya.QMessageBox_StandardButton(warning.exec_())
        return False, False, False, False, False, False, False, False
    # starting with the opt_in label, identify the sub-circuit, then GCs
    detector_GCs = [c for c in components if [p for p in c.pins if p.type == PIN_TYPES.OPTICALIO] if (
        c.trans.disp - components_sorted[0].trans.disp).to_p() != pya.DPoint(0, 0)]
    if verbose:
        print("   N=%s, detector GCs: %s" %
              (len(detector_GCs), [c.display() for c in detector_GCs]))
    vect_optin_GCs = [(c.trans.disp - components_sorted[0].trans.disp).to_p()
                      for c in detector_GCs]

    # Laser at the opt_in GC:
    p = [p for p in components_sorted[0].pins if p.type == PIN_TYPES.OPTICALIO]
    p[0].pin_name += '_laser'
    laser_net = p[0].net = Net(idx=p[0].pin_name, pins=p)
    if verbose:
        print(" - pin_name: %s" % (p[0].pin_name))

    if DFT['design-for-test']['tunable-laser']['wavelength'] == opt_in_dict[0]['wavelength']:
        wavelength_start, wavelength_stop, wavelength_points = float(DFT['design-for-test']['tunable-laser']['wavelength-start']), float(
            DFT['design-for-test']['tunable-laser']['wavelength-stop']), int(DFT['design-for-test']['tunable-laser']['wavelength-points'])
    else:
        warning = pya.QMessageBox()
        warning.setStandardButtons(pya.QMessageBox.Ok)
        warning.setText("No laser at %s nm is available. Tunable laser definition is in the technology's DFT.xml file." %
                        opt_in_dict[0]['wavelength'])
        pya.QMessageBox_StandardButton(warning.exec_())
        return False, False, False, False, False, False, False, False

    if opt_in_dict[0]['pol'] == 'TE':
        orthogonal_identifier = 1
    elif opt_in_dict[0]['pol'] == 'TM':
        orthogonal_identifier = 2
    else:
        warning = pya.QMessageBox()
        warning.setStandardButtons(pya.QMessageBox.Ok)
        warning.setText("Unknown polarization: %s." % opt_in_dict[0]['pol'])
        pya.QMessageBox_StandardButton(warning.exec_())
        return False, False, False, False, False, False, False, False
    ignoreOpticalIOs = False

    # find the GCs in the circuit and connect detectors based on DFT rules
    ports_info = []
    port_name = 0
    detector_lookuptable = {1: 1, -1: 2, -2: 3}
    detector_list = []
    for d in list(range(int(DFT['design-for-test']['grating-couplers']['detectors-above-laser']) + 0, 0, -1)) + list(range(-1, -int(DFT['design-for-test']['grating-couplers']['detectors-below-laser']) - 1, -1)):
        if pya.DPoint(0, d * float(DFT['design-for-test']['grating-couplers']['gc-pitch']) * 1000) in vect_optin_GCs:
            port_name += 1
            detector_list += [detector_lookuptable[d]]
            index = vect_optin_GCs.index(pya.DPoint(
                0, d * float(DFT['design-for-test']['grating-couplers']['gc-pitch']) * 1000))
            # detector_GCs[index] # component

            p = [p for p in detector_GCs[index].pins if p.type == PIN_TYPES.OPTICALIO]
            p[0].pin_name += '_detector' + str(port_name)
            p[0].net = Net(idx=p[0].pin_name, pins=p)
            ports_info.append(port_info(p[0].net, port_name))
            if verbose:
                print(" - pin_name: %s" % (p[0].pin_name))

    # Sort the detectors:
    ports_info2 = sorted(ports_info, key=lambda d: d.port_name)

    # output:
    detector_nets = []
    for d in ports_info2:
        detector_nets.append(d.detector_net)

    return laser_net, detector_nets, wavelength_start, wavelength_stop, wavelength_points, orthogonal_identifier, ignoreOpticalIOs, detector_list


def identify_nets_ansys(self, verbose=False):
    # function to identify all the nets in the cell layout
    # use the data in Optical_pin, Optical_waveguide to find overlaps
    # and save results in components
    PIN_TYPES = enum('OPTICALIO', 'OPTICAL', 'ELECTRICAL')
    PIN_LENGTH = 100  # 0.1 micron
    if verbose:
        print("SiEPIC.extend.identify_nets():")

    #from SiEPIC import _globals
    from SiEPIC.core import Net

    # output: array of Net[]
    nets = []
    from SiEPIC.extend import find_components, find_pins
    # find components and pins in the cell layout
    components = self.find_components()
    pins = self.find_pins()

    # Optical Pins:
    optical_pins = [p for p in pins if p.type == PIN_TYPES.OPTICAL]

    # Loop through all pairs components (c1, c2); only look at touching components
    for c1 in components:
        for c2 in components[c1.idx + 1: len(components)]:
            if verbose:
                print(" - Components: [%s-%s], [%s-%s].  Pins: %s, %s" % (c1.component, c1.idx, c2.component, c2.idx, c1.pins, c2.pins))

            if c1.polygon.bbox().overlaps(c2.polygon.bbox()) or c1.polygon.bbox().touches(c2.polygon.bbox()):
                # Loop through all the pins (p1) in c1
                # - Compare to all other pins, find other overlapping pins (p2) in c2
                for p1 in [p for p in c1.pins if p.type == PIN_TYPES.OPTICAL]:
                    for p2 in [p for p in c2.pins if p.type == PIN_TYPES.OPTICAL]:
                        if verbose:
                            print(" - Components, pins: [%s-%s, %s, %s, %s], [%s-%s, %s, %s, %s]; difference: %s"
                                  % (c1.component, c1.idx, p1.pin_name, p1.center, p1.rotation, c2.component, c2.idx, p2.pin_name, p2.center, p2.rotation, p1.center - p2.center))
                        # check that pins are facing each other, 180 degree
                        check1 = ((p1.rotation - p2.rotation) % 360) == 180

                        # check that the pin centres are perfectly overlapping
                        # (to avoid slight disconnections, and phase errors in simulations)
                        check2 = (p1.center == p2.center)

                        if check1 and check2:  # found connected pins:
                            # make a new optical net index
                            net_idx = len(nets)
                            # optical net connects two pins; keep track of the pins, Pin[] :
                            nets.append(
                                Net(idx=net_idx, pins=[p1, p2], _type=PIN_TYPES.OPTICAL))
                            # assign this net number to the pins
                            p1.net = nets[-1]
                            p2.net = nets[-1]

                            if verbose:
                                print(" - pin-pin, net: %s, component, pin: [%s-%s, %s, %s, %s], [%s-%s, %s, %s, %s]"
                                      % (net_idx, c1.component, c1.idx, p1.pin_name, p1.center, p1.rotation, c2.component, c2.idx, p2.pin_name, p2.center, p2.rotation))

    return nets, components
    
    
    



def check_components_models_ansys():

    # Check if all the components in the cell have compact models loaded in INTERCONNECT

    # test for Component.has_compactmodel()
    from .utils import get_layout_variables
    TECHNOLOGY, lv, ly, cell = get_layout_variables()

    print("* find_components()")
    components = cell.find_components()
    print("* Display list of components")

    if not all([c.has_model() for c in components]):
        # missing models, find which one
        components_havemodels = [[c.has_model(), c.component, c.instance] for c in components]
        missing_models = []
        for c in components_havemodels:
            if c[0] == False:
                missing_models.append([c[1], c[2]])
        missing = ("We have %s component(s) missing models, as follows: %s" %
                   (len(missing_models), missing_models))
        v = pya.MessageBox.warning("Errors", missing, pya.MessageBox.Ok)
    else:
        print('check_components_models(): all models are present.')
        v = pya.MessageBox.warning(
            "All ok", "All components have models. Ok to simulate the circuit.", pya.MessageBox.Ok)


def spice_netlist_export_ansys(self, verbose=True, opt_in_selection_text=[],i=1,loc_x=[], loc_y =[]):
    import SiEPIC

    #from . import _globals

    # get coordinates of the circuits in the layout
    #circuit_number = number
    cell_name = self.basic_name()
    print("Generating netlist for %s" %cell_name)
    # Define enumeration for pins
    PIN_TYPES = enum('OPTICALIO', 'OPTICAL', 'ELECTRICAL')
    PIN_LENGTH = 100  # 0.1 micron
    from time import strftime
    from SiEPIC.utils import eng_str

    from SiEPIC.utils import get_technology
    TECHNOLOGY = get_technology()
    if not TECHNOLOGY['technology_name']:
        v = pya.MessageBox.warning("Errors", "Lumerical-Tools requires a technology to be chosen.  \n\nThe active technology is displayed on the bottom-left of the KLayout window, next to the T. \n\nChange the technology using KLayout File | Layout Properties, then choose Technology and find the correct one (e.g., EBeam, GSiP).", pya.MessageBox.Ok)
        return '', '', 0, []

    # get the netlist from the entire layout
    nets, components = self.identify_nets_ansys(verbose=True)
    if not components:
        v = pya.MessageBox.warning("Errors", "No components found.", pya.MessageBox.Ok)
        return '', '', 0, []

    # Get information about the laser and detectors:
    # this updates the Optical IO Net
    laser_net, detector_nets, port_names= \
        get_LumericalINTERCONNECT_analyzers_ansys(self, components, verbose=verbose)
    detector_list = []

    
    # if Laser and Detectors are not defined
    #if not laser_net or not detector_nets: Ansys: detector_nets changed to port_names
    if not port_names:
        # Use opt_in labels
        port_names, components = get_LumericalINTERCONNECT_analyzers_from_text_label_ansys(self,components, verbose=verbose)

    if not components:
        pya.MessageBox.warning("Error: netlist extraction",
                               "Error: netlist extraction. No components found. Please make sure to select the correct Technology.", pya.MessageBox.Ok)
        return '', '', 0, []
    verbose = True
    if verbose:
        print("* Display list of components:")
        [c.display() for c in components]
        print("* Display list of nets:")
        [n.display() for n in nets]

    text_main = '* Spice output from KLayout Lumerical Package, %s.\n\n' % (
        strftime("%Y-%m-%d %H:%M:%S"))
    text_subckt = text_main

    # convert KLayout GDS rotation/flip to Lumerical INTERCONNECT
    # KLayout defines mirror as an x-axis flip, whereas INTERCONNECT does y-axis flip
    # KLayout defines rotation as counter-clockwise, whereas INTERCONNECT does clockwise
    # input is KLayout Rotation,Flip; output is INTERCONNECT:
    KLayoutInterconnectRotFlip = \
        {(0, False): [0, False],
         (90, False): [270, False],
         (180, False): [180, False],
         (270, False): [90, False],
         (0, True): [180, True],
         (90, True): [90, True],
         (180, True): [0, True],
         (270, True): [270, False]}

    # Determine the Layout-to-Schematic (x,y) coordinate scaling
    # Find the distances between all the components, in order to determine scaling
    sch_positions = [o.Dcenter for o in components]
    sch_distances = []
    for j in range(len(sch_positions)):
        for k in range(j + 1, len(sch_positions)):
            dist = (sch_positions[j] - sch_positions[k]).abs()
            sch_distances.append(dist)
    sch_distances.sort()
    if verbose:
        print("Distances between components: %s" % sch_distances)
    # remove any 0 distances:
    while 0.0 in sch_distances:
        sch_distances.remove(0.0)
    # scaling based on nearest neighbour:
    Lumerical_schematic_scaling = 0.6 / min(sch_distances)
    print("Scaling for Lumerical INTERCONNECT schematic: %s" % Lumerical_schematic_scaling)
    # but if the layout is too big, limit the size
    MAX_size = 0.05 * 1e3
    if max(sch_distances) * Lumerical_schematic_scaling > MAX_size:
        Lumerical_schematic_scaling = MAX_size / max(sch_distances)
    print("Scaling for Lumerical INTERCONNECT schematic: %s" % Lumerical_schematic_scaling)

    # find electrical IO pins
    electricalIO_pins = ""
    DCsources = ""  # string to create DC sources for each pin
    Vn = 1
    SINGLE_DC_SOURCE = 2
    # (1) attach all electrical pins to the same DC source
    # (2) or to individual DC sources
    # (3) or choose based on number of DC sources, if > 5, use single DC source

    # create individual sources:
    for c in components:
        for idx,p in enumerate(c.pins):
            if p.type == PIN_TYPES.ELECTRICAL:
                if idx > 0:  
                    if p.pin_name == c.pins[idx - 1].pin_name: continue  # Skip pins that have exactly the same name (assume they are internally connected in the component)
                # just want port name other than component_port name in the netlist
                NetName = " " + c.component + '_' + str(c.idx) + '_' + p.pin_name
                electricalIO_pins += NetName
                DCsources += "N" + \
                    str(Vn) + NetName + \
                    " dcsource amplitude=0 sch_x=%s sch_y=%s\n" % (-2 - Vn / 3., -2 + Vn / 8.)
                Vn += 1
    electricalIO_pins_subckt = electricalIO_pins

    # create 1 source
    if (SINGLE_DC_SOURCE == 1) or ((SINGLE_DC_SOURCE == 3) and (Vn > 5)):
        electricalIO_pins_subckt = ""
        for c in components:
            for p in c.pins:
                if p.type == PIN_TYPES.ELECTRICAL:
                    # to avoid same pin names in multiple circuits generation, all the "N$" are replaced by cell_name
                    #NetName = " N$"
                    NetName = " " + cell_name
                    electricalIO_pins_subckt += NetName
                    DCsources = "N1" + NetName + " dcsource amplitude=0 sch_x=-2 sch_y=0\n"

    # find optical IO pins
    opticalIO_pins = ''
    for c in components:
        for p in c.pins:
            if p.type == PIN_TYPES.OPTICALIO:
                # only enbale lumerical_ports to be added to the netlist as opticalIOs
                NetName = ' ' + p.pin_name
                print(p.pin_name)
                opticalIO_pins += NetName


    circuit_name = self.name.replace('.', '')  # remove "."
    if '_' in circuit_name[0]:
        circuit_name = ''.join(circuit_name.split('_', 1))  # remove leading _
    
    # only use lumerical_ports as optical IO
    opticalIO_pins_1 = port_names
    # create the top subckt:
    text_subckt += '.subckt %s%s%s\n' % (circuit_name, electricalIO_pins, opticalIO_pins_1)

    all_nets = ''
    for c in components:
        # Check pins to see if explicitly ordered numerically - requires first character in pin name to be a number (Stefan Preble, RIT)
        explicit_ordering = False
        for p in c.pins:
            pinname1 = p.pin_name[0]
            if pinname1.isdigit():
                explicit_ordering = True
            else:
                explicit_ordering = False  # all pins must be numbered
                break
               
        nets_str = ''        
        if explicit_ordering:   # Order the pins numerically (Stefan Preble, RIT)
            for idx, p in enumerate(c.pins):
                if idx > 0:  
                    if p.pin_name == c.pins[idx - 1].pin_name: continue  # Skip pins that have exactly the same name (assume they are internally connected in the component)
                if p.type == PIN_TYPES.ELECTRICAL:
                    nets_str += " " + c.component + '_' + str(c.idx) + '_' + p.pin_name
                if p.type == PIN_TYPES.OPTICALIO:
                    nets_str += " " + str(p.net.idx)    
                if p.type == PIN_TYPES.OPTICAL:
                    nets_str += " " + " N$" + str(p.net.idx)
                else:
                    nets_str += " " + " N$" + str(p.net.idx)
        else:
            
            # do one loop instead of 3 loops for ordering nets
            for p in c.pins:
                if p.type == PIN_TYPES.ELECTRICAL:
                    nets_str += " " + c.component + '_' + str(c.idx) + '_' + p.pin_name   
                if p.type == PIN_TYPES.OPTICALIO:
                    # To change nets_str for text label netlisting
                    nets_str += " " + p.pin_name                    
                if p.type == PIN_TYPES.OPTICAL:
                    import random
                    
                    if p.net.idx or p.net.idx == 0:
                        nets_str += " " + cell_name + str(p.net.idx)
                    else: 
                        import random
                        num = random.randint(1,99)
                        str_1 = " " + " N$$" + str(num)
                        # check if pin net has repeataions 
                        while all_nets.__contains__(str_1):
                            num = random.randint(1,99)
                            str_1 = " " + " N$$" + str(num)
                            print("Unique pin: %s" %str_1)
                            print("nets_str: %s" %all_nets)
                        nets_str += " " + " N$$" + str(num)
                        all_nets += " " + " N$$" + str(num) 
                        #nets_str += " " + " N$" + str(random.randint(1,99))
              
                    
            '''
            # optical nets: must be ordered electrical, optical IO, then optical
            for p in c.pins:
                if p.type == PIN_TYPES.ELECTRICAL:
                    nets_str += " " + c.component + '_' + str(c.idx) + '_' + p.pin_name
            
            for p in c.pins:
                if p.type == PIN_TYPES.OPTICALIO:
                    #nets_str += " " + str(p.net.idx)
                    #Ansys: to change nets_str for text label netlisting
                    nets_str += " " + p.pin_name
            for p in c.pins:
                if p.type == PIN_TYPES.OPTICAL:
                    #nets_str += " N$" + str(p.net.idx)
                    nets_str += " " + cell_name + str(p.net.idx)
            else:
                nets_str += " " + cell_name + str(p.net.idx)
            '''
        trans = KLayoutInterconnectRotFlip[(c.trans.angle, c.trans.is_mirror())]

        flip = ' sch_f=true' if trans[1] else ''
        if trans[0] > 0:
            rotate = ' sch_r=%s' % str(trans[0])
        else:
            rotate = ''

        # Check to see if this component is an Optical IO type.
        pinIOtype = any([p for p in c.pins if p.type == PIN_TYPES.OPTICALIO])

        #Ansys: remove igonreOpticalIOs
        ignoreOpticalIOs = False
        if ignoreOpticalIOs and pinIOtype:
            # Replace the Grating Coupler or Edge Coupler with a 0-length waveguide.
            component1 = "ebeam_wg_strip_1550"
            params1 = "wg_length=0u wg_width=0.500u"
        else:
            component1 = c.component
            params1 = c.params
        
        # Remove "$N" from component's name for cell instance arrays of the same name     
        if "$" in component1:
            component1 = component1[:component1.find("$")]

        text_subckt += ' %s %s %s ' % (component1.replace(' ', '_') +
                                       "_" + str(c.idx), nets_str, component1.replace(' ', '_'))
        if c.library != None:
            text_subckt += 'library="%s" ' % c.library
        x, y = c.Dcenter.x, c.Dcenter.y
        text_subckt += '%s lay_x=%s lay_y=%s sch_x=%s sch_y=%s %s%s\n' % \
            (params1,
             eng_str(x * 1e-6), eng_str(y * 1e-6),
             eng_str(x * Lumerical_schematic_scaling), eng_str(y * Lumerical_schematic_scaling),
             rotate, flip)

    text_subckt += '.ends %s\n\n' % (circuit_name)

    if laser_net:
        text_main += '* Optical Network Analyzer:\n'
        text_main += '.ona input_unit=wavelength input_parameter=start_and_stop\n  + minimum_loss=80\n  + analysis_type=scattering_data\n  + multithreading=user_defined number_of_threads=1\n'
        text_main += '  + orthogonal_identifier=%s\n' % orthogonal_identifier
        text_main += '  + start=%4.3fe-9\n' % wavelength_start
        text_main += '  + stop=%4.3fe-9\n' % wavelength_stop
        text_main += '  + number_of_points=%s\n' % wavelength_points
        for i in range(0, len(detector_nets)):
            text_main += '  + input(%s)=%s,%s\n' % (i + 1, circuit_name, detector_nets[i].idx)
        text_main += '  + output=%s,%s\n' % (circuit_name, laser_net.idx)
    

    # main circuit
    if loc_x:
        text_subckt += '%s %s %s %s sch_x=%f sch_y=%f ' %\
            (circuit_name, electricalIO_pins_subckt, opticalIO_pins_1, circuit_name, float(loc_x[i-1])/200 , float(loc_y[i-1])/200)
    else:
        text_subckt += '%s %s %s %s sch_x=-1 sch_y=-1 ' %\
            (circuit_name, electricalIO_pins_subckt, opticalIO_pins_1, circuit_name)
    if len(DCsources) > 0:
        text_subckt += 'sch_r=270\n\n'
    else:
        text_subckt += '\n\n'

    text_main += DCsources

    return text_subckt, text_main, len(detector_nets), detector_list, opticalIO_pins_1



def selected_circuits_netlists_ansys(verbose = False, simulate = False):


    from SiEPIC.utils import get_technology, get_technology_by_name, select_instances
    from SiEPIC.extend import find_components, find_pins
    from SiEPIC.utils import get_layout_variables
    import os
    import random

    #TECHNOLOGY, lv, ly, cell = get_layout_variables()
  
    print("=======================================================")
    print("Netlisting for multiple circuits")
  
    verbose = True
    cells = []
  
    # Get the index of the selected cells 
    selection = select_instances()
    print("num of selected inst: %s" %len(selection))
    #selection = lv.object_selection
    for i in range(len(selection)):
        index = selection[i].cell_index()
        cell_name = selection[i].layout().cell(index).basic_name()
        print("cell number %s: [cell_index number: %s; cell name: %s]" %((i+1),index, cell_name))
        cells.append(selection[i].layout().cell(index))
    

    # get names of the selected cells
    selected_names = []
    for i in range(len(selection)):
        index = selection[i].cell_index()
        cell = selection[i].layout().cell(index)
        cell_name = cell.basic_name()
        print("index: %s, cell name: %s"%(index, cell_name))
        selected_names.append(cell_name)

    print("selected names: %s" %selected_names)



    # define a list containing all optical pins to check duplicates
    all_optical_pins = []
    TECHNOLOGY, lv, ly, cell = get_layout_variables()
    top = ly.top_cell()
    loc_x = []
    loc_y = []
    iter = pya.RecursiveInstanceIterator.new(ly,top)
    while not(iter.at_end()):
        if iter.inst_cell().basic_name() in selected_names:
            r,t,loc = iter.inst_dtrans().to_s().split()
            x,y = loc.split(",")
            loc_x.append(x)
            loc_y.append(y)
        iter.next()

    # export netlist only
    if not simulate:
        # get saving dir and check name conflic
        file_dir = pya.FileDialog.ask_save_file_name('Please specify the save diretory for your netlist','default','.spi')
        if file_dir != None:
            dirname = os.path.dirname(file_dir)
            base_name = os.path.basename(file_dir)
            filename_subckt_1 = os.path.join(dirname,  '%s.spi' % base_name)
            if os.path.isfile(filename_subckt_1):
                raise Exception('There already exists a file with this name. Please choose a different file name!')
            file = open(filename_subckt_1, 'w')


            # Get netlists for all selected circuits and write into user defined spi file
            for i in range(len(selection)):
                text_Spice, text_Spice_main, num_detectors, detector_list, optical_pins = cells[i].spice_netlist_export_ansys(verbose=False, opt_in_selection_text=[], i=i, loc_x = loc_x, loc_y = loc_y)
                if not text_Spice:
                    raise Exception("No netlist available. Cannot run simulation.")
                all_optical_pins.extend(optical_pins.split())
                if verbose:   
                    print(text_Spice)
                    print('-------------------------')
                file.write(text_Spice)
            file.close()
            print(all_optical_pins)
            if has_duplicates(all_optical_pins):
                raise Exception("Optical ports with the same names are found. Please rename your optical ports so each optical port has an unique name!")

    # import netlist to INTC
    if simulate:
        file_path_1 = os.path.join(pya.Application.instance().inst_path(),'working_dir.txt')
        file_path_2 = os.path.join(pya.Application.instance().klayout_path()[0], 'working_dir.txt')
        if(os.path.isfile(file_path_1)==True):
            f = open(file_path_1, "r")
            path = f.read()
        elif(os.path.isfile(file_path_2)==True):
            f = open(file_path_2, "r")
            path = f.read()
        else:
            raise Exception("Warning: The program could not find your project directory. Please specify your project directory manually (Ansys Lumerical > Setup Project Directory).")

        filename_subckt = os.path.join(path, '%s_circuits_selected_%i.spi' %(len(selection), random.randint(1,101)))
        file = open(filename_subckt, 'w')

        
        for i in range(len(selection)):
            text_Spice, text_Spice_main, num_detectors, detector_list, optical_pins = cells[i].spice_netlist_export_ansys(verbose=False, opt_in_selection_text=[], i=i, loc_x = loc_x, loc_y = loc_y)
            if not text_Spice:
                raise Exception("No netlist available. Cannot run simulation.")
            all_optical_pins.extend(optical_pins.split())
            if verbose:   
                print(text_Spice)
                print('-------------------------')
            file.write (text_Spice)
        file.close()

        if has_duplicates(all_optical_pins):
            raise Exception("Optical ports with the same names are found. Please rename your optical ports so each optical port has an unique name!")
            return

        import os
        from Lumerical.interconnect_ansys import INTC_commandline
        netlist_path = r'%s' % filename_subckt
        INTC_commandline(netlist_path)

    '''
    cell_1.spice_netlist_export_ansys(verbose=True, opt_in_selection_text=[])
    text_Spice, text_Spice_main, num_detectors, detector_list = cell_1.spice_netlist_export_ansys(verbose=True, opt_in_selection_text=[])
    
    if not text_Spice:
        raise Exception("No netlist available. Cannot run simulation.")
    if verbose:   
        print(text_Spice)
    
    circuit_name = cell_1.name.replace('.','') # remove "."
    if '_' in circuit_name[0]:
        circuit_name = ''.join(circuit_name.split('_', 1))  # remove leading _
    
    if file_dir != None:
        dirname = os.path.dirname(file_dir)
        base_name = os.path.basename(file_dir)
        filename_subckt_1 = os.path.join(dirname,  '%s.spi' % base_name)
        if os.path.isfile(filename_subckt_1):
        raise Exception('There already exists a file with this name. Please choose a different file name!')
        else:
        file = open(filename_subckt_1, 'w')
        file.write (text_Spice)
        file.close()
    '''



pya.Cell.identify_nets_ansys = identify_nets_ansys
pya.Cell.get_LumericalINTERCONNECT_analyzers_ansys = get_LumericalINTERCONNECT_analyzers_ansys
pya.Cell.get_LumericalINTERCONNECT_analyzers_from_opt_in_ansys = get_LumericalINTERCONNECT_analyzers_from_opt_in_ansys
pya.Cell.spice_netlist_export_ansys = spice_netlist_export_ansys
pya.Cell.get_LumericalINTERCONNECT_analyzers_from_text_label_ansys = get_LumericalINTERCONNECT_analyzers_from_text_label_ansys