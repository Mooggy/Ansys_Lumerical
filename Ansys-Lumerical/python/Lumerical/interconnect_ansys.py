# Copyright (c) 2021, Ansys, Inc. All rights reserved.

'''
Circuit simulations using Lumerical INTERCONNECT and a Compact Model Library

- run_INTC: run INTERCONNECT using Python integration
- INTC_commandline: invoke INTC via the command line, with an lsf file as input.
- Setup_Lumerical_KLayoutPython_integration
    Configure PATH env, import lumapi, run interconnect, 
    Install technology CML, read CML elements
- circuit_simulation: netlist extract and run simulation

usage:
 import Lumerical.interconnect_ansys


################################################################################
'''



import pya
from Lumerical.netlister import *


class INTC_GUI(pya.QDialog):
  def getInputs(self):
    return (str(self.target.text), str(self.le.text))
    print("Selected file: " + test_name)
    print("Selected target: "+ target_name)
    
  def cancel(self):
    self.btnlistener = 'cancel'
    self.close()


  def ok(self):
    self.btnlistener = 'ok'
    self.getInputs()
    self.close()
    test_name = str(self.le.text)
    target_name = str(self.target.text)
    return (str(self.target.text), str(self.le.text))
       
  def getfile(self):
    fname = pya.QFileDialog.getOpenFileName(self,'Load Test File for INTERCONNECT','c:\\',"Test files(*.icp *.lsf)")
    self.le.setText(str(fname))
    pya.Application.instance().set_config('test_file_path', str(fname))
    return(str(fname))
    
  def btnListener(self):
    return self.btnlistener
    
  def __init__(self, parent = None):
    """ Dialog constructor """
    
    super(INTC_GUI, self).__init__()

    self.setWindowTitle("INTERCONNECT Testbench Setup (Optional)")

    mainLayout = pya.QVBoxLayout(self)
    self.setLayout(mainLayout)
    
    self.resize(600, 150)
    self.btnlistener = ''
    #define v
    layout_h1 = pya.QHBoxLayout(self)
    layout_h2 = pya.QHBoxLayout(self)
    #layout_h3 = pya.QHBoxLayout(self)

    
    # add netlist target editor  
    self.target_ins = pya.QLabel("Compound Element Name:", self)
    self.target_ins.setFont(pya.QFont('Times', 12))
    self.target_ins.setFixedWidth(200)
    layout_h1.addWidget(self.target_ins)
    self.target = pya.QLineEdit("Default", self)
    self.target.setFont(pya.QFont('Times', 12))
    #self.target.setFixedWidth(450)
    layout_h1.addWidget(self.target)
  
    
    # add image
    self.image = pya.QLabel("Test File (.lsf or .icp):", self)
    self.image.setFont(pya.QFont('Times', 12))
    self.image.setFixedWidth(200)
    layout_h2.addWidget(self.image)
    
    # show selected test file
    try:
      self.fname = pya.Application.instance().get_config('test_file_path')
      self.le = pya.QLineEdit(str(self.fname),self)
    except:
      self.le = pya.QLineEdit("",self)
    self.le.setFont(pya.QFont('Times', 12))
    layout_h2.addWidget(self.le)  
    
    #button = pya.QPushButton('Load .lsf or .icp files', self)
    button = pya.QPushButton('Load', self)
    button.setFont(pya.QFont('Times', 12, pya.QFont.Bold))
    #button.setFont(pya.QFont('Times', 12))
    #self.button.setFixedWidth(100)
    layout_h2.addWidget(button)
    

    
    # add display test file button
    '''
    self.la = pya.QLabel("Selected test file:")
    self.la.setFont(pya.QFont('Times', 12))
    layout_h3.addWidget(self.la) 
    '''     



    # OK Cancel
    layout_h4 = pya.QHBoxLayout(self);
    ok = pya.QPushButton("OK",self)
    ok.setFont(pya.QFont('Times', 12))
    fname = ok.clicked(self.ok)
    cancel = pya.QPushButton("Cancel",self)
    cancel.setFont(pya.QFont('Times', 12))
    cancel.clicked(self.cancel)
    layout_h4.addWidget(ok)
    layout_h4.addWidget(cancel)
    
    
    #set layout

    mainLayout.addLayout(layout_h2)
    #mainLayout.addLayout(layout_h3)
    mainLayout.addLayout(layout_h1)
    mainLayout.addLayout(layout_h4)

    #self.setLayout(layout_v)
    
    # get file action once clicked
    self.fname = button.clicked(self.getfile)





def run_INTC(verbose=False):
  from SiEPIC import _globals
  from SiEPIC.lumerical import load_lumapi
  lumapi = _globals.LUMAPI
  if not lumapi:
    print("SiEPIC.lumerical.interconnect.run_INTC: lumapi not loaded; reloading load_lumapi.")
    import sys
    if sys.version_info[0] == 3:
        if sys.version_info[1] < 4:
            from imp import reload
        else:
            from importlib import reload
    elif sys.version_info[0] == 2:
        from imp import reload    
    reload(load_lumapi)

  if not lumapi:
    print("SiEPIC.lumerical.interconnect.run_INTC: lumapi not loaded")
    pya.MessageBox.warning("Cannot load Lumerical Python integration.", "Cannot load Lumerical Python integration. \nSome SiEPIC-Tools Lumerical functionality will not be available.", pya.MessageBox.Cancel)
    return
  
  if verbose:
    print(_globals.INTC)  # Python Lumerical INTERCONNECT integration handle
  
  if not _globals.INTC: # Not running, start a new session
    _globals.INTC = lumapi.open('interconnect')
    if verbose:
      print(_globals.INTC)  # Python Lumerical INTERCONNECT integration handle
  else: # found open INTC session
    try:
      lumapi.evalScript(_globals.INTC, "?'KLayout integration test.\n';\n")
    except: # but can't communicate with INTC; perhaps it was closed by the user
      _globals.INTC = lumapi.open('interconnect') # run again.
      if verbose:
        print(_globals.INTC)  # Python Lumerical INTERCONNECT integration handle
  try: # check again
    lumapi.evalScript(_globals.INTC, "?'KLayout integration test.\n';\n")
  except:
    raise Exception ("Can't run Lumerical INTERCONNECT via Python integration.")


def Setup_Lumerical_KLayoutPython_integration(verbose=False):
  import sys, os, string, pya

  from ..utils import get_technology, get_technology_by_name
  # get current technology
  TECHNOLOGY = get_technology(query_activecellview_technology=True)  
  # load more technology details (CML file location)
  TECHNOLOGY = get_technology_by_name(TECHNOLOGY['technology_name'])

  # location for the where the CMLs will locally be installed:
  dir_path = os.path.join(pya.Application.instance().application_data_path(), 'Lumerical_CMLs')

  try:
    libraries = [n for n in pya.Library.library_names() if (pya.Library.library_by_name(n).technology == TECHNOLOGY['technology_name'])]
    for n in [pya.Library.library_by_name(l) for l in libraries]:
      print(n.layout().meta_info_value("path"))
  except:
    pass

  question = pya.QMessageBox()
  question.setStandardButtons(pya.QMessageBox.Yes | pya.QMessageBox.No)
  question.setDefaultButton(pya.QMessageBox.Yes)
  question.setText("SiEPIC-Tools will install the Compact Model Library (CML) in Lumerical INTERCONNECT for the currently active technology. \nThis includes the libraries %s.  \nProceed?" % libraries)
  informative_text = "\nTechnology: %s\n" % TECHNOLOGY['technology_name']
  for i in range(0,len(TECHNOLOGY['INTC_CMLs_path'])):
    informative_text += "Source CML file {}: {}\n".format(i+1, TECHNOLOGY['INTC_CMLs_path'][i])
  informative_text += "Install location: %s" % dir_path
  question.setInformativeText(informative_text)
  if(pya.QMessageBox_StandardButton(question.exec_()) == pya.QMessageBox.No):
    return

  
  ##################################################################
  # Load Lumerical API: 
  from .. import _globals
  run_INTC()
  lumapi = _globals.LUMAPI
  if not lumapi:
    print('SiEPIC.lumerical.interconnect.Setup_Lumerical_KLayoutPython_integration: lumapi not loaded')
    return

  import os 
  # Read INTC element library
  lumapi.evalScript(_globals.INTC, "out=library;")
  _globals.INTC_ELEMENTS=lumapi.getVar(_globals.INTC, "out")

  # Install technology CML if missing in INTC
  # check if the latest version of the CML is in KLayout's tech
  if not ("design kits::"+TECHNOLOGY['technology_name'].lower()+"::"+TECHNOLOGY['INTC_CML_version'].lower().replace('.cml','').lower()) in _globals.INTC_ELEMENTS:
    # install CML
    print("Lumerical INTC, installdesignkit ('%s', '%s', true);" % (TECHNOLOGY['INTC_CML_path'], dir_path ) )
    lumapi.evalScript(_globals.INTC, "installdesignkit ('%s', '%s', true);" % (TECHNOLOGY['INTC_CML_path'], dir_path ) )
    
  # Install other CMLs within technology
  for i in range(0,len(TECHNOLOGY['INTC_CMLs_name'])):
    if not ("design kits::"+TECHNOLOGY['INTC_CMLs_name'][i].lower()+"::"+TECHNOLOGY['INTC_CMLs_version'][i].lower().replace('.cml','').lower()) in _globals.INTC_ELEMENTS:
        # install CML
        print("Lumerical INTC, installdesignkit ('%s', '%s', true);" % (TECHNOLOGY['INTC_CMLs_path'][i], dir_path ) )
        lumapi.evalScript(_globals.INTC, "installdesignkit ('%s', '%s', true);" % (TECHNOLOGY['INTC_CMLs_path'][i], dir_path ) )
        
  # Re-Read INTC element library
  lumapi.evalScript(_globals.INTC, "out=library;")
  _globals.INTC_ELEMENTS=lumapi.getVar(_globals.INTC, "out")
  # Close INTERCONNECT so that the library information is saved, then re-open
  lumapi.close(_globals.INTC)
  run_INTC()

  # Save INTC element library to KLayout application data path
  if not os.path.exists(dir_path):
    os.makedirs(dir_path)
  fh = open(os.path.join(dir_path,"Lumerical_INTC_CMLs.txt"), "w")
  fh.writelines(_globals.INTC_ELEMENTS)
  fh.close()
  
  integration_success_message = "message('KLayout-Lumerical INTERCONNECT integration successful, CML library/libraries:\n"
  for cml_name in TECHNOLOGY['INTC_CMLs_name']:
    integration_success_message += "design kits::"+cml_name.lower()+"\n"
  integration_success_message += "');switchtodesign;\n"
  lumapi.evalScript(_globals.INTC, integration_success_message)

  # instantiate all library elements onto the canvas
  question = pya.QMessageBox()
  question.setStandardButtons(pya.QMessageBox.Yes | pya.QMessageBox.No)
  question.setDefaultButton(pya.QMessageBox.Yes)
  question.setText("Do you wish to see all the components in %s library?" % TECHNOLOGY['technology_name'])
#  question.setInformativeText("Do you wish to see all the components in the library?")
  if(pya.QMessageBox_StandardButton(question.exec_()) == pya.QMessageBox.No):
    # lumapi.evalScript(_globals.INTC, "b=0:0.01:10; plot(b,sin(b),'Congratulations, Lumerical is now available from KLayout','','Congratulations, Lumerical is now available from KLayout');")
    return
  intc_elements = _globals.INTC_ELEMENTS.split('\n')
#  tech_elements = [ e.split('::')[-1] for e in intc_elements if "design kits::"+TECHNOLOGY['technology_name'].lower()+"::" in e ]
  tech_elements = [ e for e in intc_elements if "design kits::"+TECHNOLOGY['technology_name'].lower()+"::" in e ]
  i, x, y, num = 0, 0, 0, len(tech_elements)
  for i in range(0, num):
    lumapi.evalScript(_globals.INTC, "a=addelement('%s'); setposition(a,%s,%s); " % (tech_elements[i],x,y) )
    y += 250
    if (i+1) % int(num**0.5) == 0:
      x += 250
      y = 0

def INTC_commandline(filename2, exct = False):
  print ("Running Lumerical INTERCONNECT using the command interface.")
  import sys, os, string
  
  if sys.platform.startswith('linux'):
    import subprocess
    # Linux-specific code here...
    print("Running INTERCONNECT")
    # Location of INTERCONNECT program (this found from RPM installation)
    file_path = '/opt/lumerical/interconnect/bin/interconnect'
    subprocess.Popen([file_path, '-run', filename2])
      
  
  elif sys.platform.startswith('darwin'):
    # OSX specific
    import sys
    if int(sys.version[0]) > 2:
      import subprocess
      subprocess.Popen(['/usr/bin/open -n /Applications/Lumerical/INTERCONNECT/INTERCONNECT.app', '-run', '--args -run %s' % filename2])          
    else:
      import commands
      print("Running INTERCONNECT")
      runcmd = ('source ~/.bash_profile; /usr/bin/open -n /Applications/Lumerical/INTERCONNECT/INTERCONNECT.app --args -run %s' % filename2)
      print("Running in shell: %s" % runcmd)
      a=commands.getstatusoutput(runcmd)
      print(a)

  
  elif sys.platform.startswith('win'):
    # Windows specific code here
    import subprocess
    print("Running INTERCONNECT")
    #check Interconnect installation directory
    #Ansys: file_path_c the version is hardcoded and needs improvement 
    file_path_1 = os.path.join(pya.Application.instance().inst_path(),'lum_path.txt')
    file_path_2 = os.path.join(pya.Application.instance().klayout_path()[0], 'lum_path.txt')
    if(os.path.isfile(file_path_1)==True):
      f = open(file_path_1, "r")
      path = f.read()
      file_path_c = r'%s/bin/interconnect.exe'%path
      if(os.path.isfile(file_path_c)==True):
        subprocess.Popen(args=[file_path_c, '-run', filename2], shell=exct)
      else:
        raise Exception("Warning: The program could not find INTERCONNECT. Please specify the the location of INTERCONNECT manually.")
    elif(os.path.isfile(file_path_2)==True):
      f = open(file_path_2, "r")
      path = f.read()
      file_path_c = r'%s/bin/interconnect.exe'%path
      if(os.path.isfile(file_path_c)==True):
        subprocess.Popen(args=[file_path_c, '-run', filename2], shell=exct)
      else:
        raise Exception("Warning: The program could not find INTERCONNECT. Please specify the the location of INTERCONNECT manually.")
    
    else:
      raise Exception("Warning: The program could not find INTERCONNECT. Please specify the the location of INTERCONNECT manually.")
    '''  
    if(os.path.isfile(file_path_a)==True):
      subprocess.Popen(args=[file_path_a, '-run', filename2], shell=exct)
    elif(os.path.isfile(file_path_b)==True):
      subprocess.Popen(args=[file_path_b, '-run', filename2], shell=exct)
    elif(os.path.isfile(file_path_c)==True):
      subprocess.Popen(args=[file_path_c, '-run', filename2], shell=exct)
    else:
      warning_window = pya.QMessageBox()
      warning_window.setText("Warning: The program could not find INTERCONNECT.")
      warning_window.setInformativeText("Do you want to specify it manually?")
      warning_window.setStandardButtons(pya.QMessageBox.Yes | pya.QMessageBox.Cancel);
      warning_window.setDefaultButton(pya.QMessageBox.Yes)
      response = warning_window.exec_()        
      if(response == pya.QMessageBox.Yes):
        dialog = pya.QFileDialog()
        path = str(dialog.getOpenFileName())
        path = path.replace('/', '\\')
        subprocess.Popen(args=[path, '-run', filename2], shell=exct)
    '''

def component_simulation(verbose=False, simulate=True):
  import sys, os, string
  from .. import _globals

  # get selected instances
  from ..utils import select_instances
  selected_instances = select_instances()

  from ..utils import get_layout_variables
  TECHNOLOGY, lv, ly, cell = get_layout_variables()
    
  
  # check that it is one or more:
  error = pya.QMessageBox()
  error.setStandardButtons(pya.QMessageBox.Ok )
  if len(selected_instances) == 0:
    error.setText("Error: Need to have a component selected.")
    return
  warning = pya.QMessageBox()
  warning.setStandardButtons(pya.QMessageBox.Yes | pya.QMessageBox.Cancel)
  warning.setDefaultButton(pya.QMessageBox.Yes)
  if len(selected_instances) > 1 :
    warning.setText("Warning: More than one component selected.")
    warning.setInformativeText("Do you want to Proceed?")
    if(pya.QMessageBox_StandardButton(warning.exec_()) == pya.QMessageBox.Cancel):
      return
  
  # Check if the component has a compact model loaded in INTERCONNECT
  # Loop if more than one component selected
  for obj in selected_instances:
# *** not working. .returns Flattened.
#    c = obj.inst().cell.find_components()[0]
    if verbose:
      print("  selected component: %s" % obj.inst().cell )
    c = cell.find_components(cell_selected=[obj.inst().cell])
    if c:
      c=c[0]
    else:
      continue
    
    if not c.has_model():
      if len(selected_instances) == 0:
        error.setText("Error: Component '%s' does not have a compact model. Cannot perform simulation." % c)
        continue

    # GUI to ask which pin to inject light into
    pin_names = [p.pin_name for p in c.pins if p.type == _globals.PIN_TYPES.OPTICAL or p.type == _globals.PIN_TYPES.OPTICALIO]
    if not pin_names:
      continue
    pin_injection = pya.InputDialog.ask_item("Pin selection", "Choose one of the pins in component '%s' to inject light into." % c.component, pin_names, 0)
    if not pin_injection:
      return
    if verbose:
      print("Pin selected from InputDialog = %s, for component '%s'." % (pin_injection, c.component) )
    
    # Write spice netlist and simulation script
    from ..utils import get_technology
    TECHNOLOGY = get_technology()  # get current technology
    import SiEPIC
    from time import strftime 
    text_main = '* Spice output from KLayout Lumerical Package v%s, %s technology (SiEPIC.lumerical.interconnect.component_simulation), %s.\n\n' % (SiEPIC.__version__, TECHNOLOGY['technology_name'], strftime("%Y-%m-%d %H:%M:%S") )



    # find electrical IO pins
    electricalIO_pins = ""
    DCsources = "" # string to create DC sources for each pin
    Vn = 1
    # (2) or to individual DC sources
    # create individual sources:
    for p in c.pins:
      if p.type == _globals.PIN_TYPES.ELECTRICAL:
        NetName = " " + c.component +'_' + str(c.idx) + '_' + p.pin_name
        electricalIO_pins += NetName
        DCsources += "N" + str(Vn) + NetName + " dcsource amplitude=0 sch_x=%s sch_y=%s\n" % (-2-Vn/3., -2+Vn/8.)
        Vn += 1
    electricalIO_pins_subckt = electricalIO_pins


    # component nets: must be ordered electrical, optical IO, then optical
    nets_str = ''
    DCsources = "" # string to create DC sources for each pin
    Vn = 1
    for p in c.pins: 
      if p.type == _globals.PIN_TYPES.ELECTRICAL:
        if not p.pin_name:
          continue
        NetName = " " + c.component +'_' + str(c.idx) + '_' + p.pin_name
        nets_str += NetName
        DCsources += "N" + str(Vn) + NetName + " dcsource amplitude=0 sch_x=%s sch_y=%s\n" % (-2-Vn/3., -2+Vn/8.)
        Vn += 1
      if p.type == _globals.PIN_TYPES.OPTICAL or p.type == _globals.PIN_TYPES.OPTICALIO:
        nets_str += " " + str(p.pin_name)


        
    # *** todo: some other way of getting this information; not hard coded.
    # GUI? Defaults from PCell?
    orthogonal_identifier=1
    wavelength_start=1500
    wavelength_stop=1600
    wavelength_points=2000
    text_main += '* Optical Network Analyzer:\n'
    text_main += '.ona input_unit=wavelength input_parameter=start_and_stop\n  + minimum_loss=80\n  + analysis_type=scattering_data\n  + multithreading=user_defined number_of_threads=1\n' 
    text_main += '  + orthogonal_identifier=%s\n' % orthogonal_identifier
    text_main += '  + start=%4.3fe-9\n' % wavelength_start
    text_main += '  + stop=%4.3fe-9\n' % wavelength_stop
    text_main += '  + number_of_points=%s\n' % wavelength_points
    for i in range(0,len(pin_names)):
      text_main += '  + input(%s)=SUBCIRCUIT,%s\n' % (i+1, pin_names[i])
    text_main += '  + output=SUBCIRCUIT,%s\n\n' % (pin_injection)

    text_main += DCsources

    text_main += 'SUBCIRCUIT %s SUBCIRCUIT sch_x=-1 sch_y=-1 \n\n' % (nets_str)
    text_main += '.subckt SUBCIRCUIT %s\n' % (nets_str)
    text_main += ' %s %s %s ' % ( c.component.replace(' ', '_') +"_1", nets_str, c.component.replace(' ', '_') ) 
    if c.library != None:
      text_main += 'library="%s" %s ' % (c.library, c.params)
    text_main += '\n.ends SUBCIRCUIT\n'

    from .. import _globals
    tmp_folder = _globals.TEMP_FOLDER
    import os    
    filename = os.path.join(tmp_folder, '%s_main.spi' % c.component)
    filename2 = os.path.join(tmp_folder, '%s.lsf' % c.component)
    filename_icp = os.path.join(tmp_folder, '%s.icp' % c.component)

    # Write the Spice netlist to file
    file = open(filename, 'w')
    file.write (text_main)
    file.close()
    if verbose:
      print(text_main)

    '''
    # Ask user whether to start a new visualizer, or use an existing one.
    opt_in_labels = [o['opt_in'] for o in opt_in]
    opt_in_labels.insert(0,'All opt-in labels')
    opt_in_selection_text = pya.InputDialog.ask_item("opt_in selection", "Choose one of the opt_in labels, to fetch experimental data.",  opt_in_labels, 0)
    if not opt_in_selection_text: # user pressed cancel
      pass
    '''    

    # Write the Lumerical INTERCONNECT start-up script.
    text_lsf =  'switchtolayout;\n'
    text_lsf +=  "cd('%s');\n" % tmp_folder
    text_lsf += 'deleteall;\n'
    text_lsf += "importnetlist('%s');\n" % filename
    text_lsf += "save('%s');\n" % filename_icp
    text_lsf += 'run;\n'
    if 0:
      for i in range(0, len(pin_names)):
        text_lsf += 'h%s = haveresult("ONA_1", "input %s/mode 1/gain");\n' % (i+1, i+1)
        text_lsf += 'if (h%s>0) { visualize(getresult("ONA_1", "input %s/mode 1/gain")); } \n' % (i+1, i+1)
    if 1:
      text_lsf += 't = "";\n'
      for i in range(0, len(pin_names)):
        text_lsf += 'h%s = haveresult("ONA_1", "input %s/mode 1/gain");\n' % (i+1, i+1)
        text_lsf += 'if (h%s>0) { t%s = getresult("ONA_1", "input %s/mode 1/gain"); t=t+"t%s,"; } \n' % (i+1, i+1, i+1, i+1)
      text_lsf += 't = substring(t, 1, length(t) - 1);\n'
      text_lsf += 'eval("visualize(" + t + ");");\n'

    file = open(filename2, 'w')
    file.write (text_lsf)
    file.close()
    if verbose:
      print(text_lsf)
    
    if simulate:
      # Run using Python integration:
      try: 
        from .. import _globals
        run_INTC()
        # Run using Python integration:
        lumapi = _globals.LUMAPI
        lumapi.evalScript(_globals.INTC, "cd ('" + tmp_folder + "');")
        lumapi.evalScript(_globals.INTC, "feval('"+ c.component + "');\n")
      except:
        from .. import scripts
        scripts.open_folder(tmp_folder)
        INTC_commandline(filename)
    else:
      from .. import scripts
      scripts.open_folder(tmp_folder)

def circuit_simulation_toolbar():
  circuit_simulation(verbose=False,opt_in_selection_text=[], matlab_data_files=[], simulate=True)

def circuit_simulation_ansys(verbose=True,opt_in_selection_text=[], matlab_data_files=[], simulate=True):
  from Lumerical.netlister import spice_netlist_export_ansys
  from Lumerical.netlister import selected_circuits_netlists_ansys
  from SiEPIC.utils import get_technology, get_technology_by_name, select_instances
  from Lumerical.netlister import has_duplicates
  from SiEPIC import _globals
  from SiEPIC.lumerical import load_lumapi
  import lumapi  

  print ('*** circuit_simulation(), opt_in: %s' % opt_in_selection_text)
  if verbose:
    print('*** circuit_simulation()')
  
  # check for supported operating system, tested on:
  # Windows 7, 10
  # OSX Sierra, High Sierra
  # Linux
  import sys
  if not any([sys.platform.startswith(p) for p in {"win","linux","darwin"}]):
    raise Exception("Unsupported operating system: %s" % sys.platform)
  
  #from SiEPIC import _globals
  from SiEPIC.utils import get_layout_variables
  TECHNOLOGY, lv, layout, topcell = get_layout_variables()
  
  # Save the layout prior to running simulations, if there are changes.
  mw = pya.Application.instance().main_window()
  if mw.manager().has_undo():
    mw.cm_save()
  layout_filename = mw.current_view().active_cellview().filename()
  if len(layout_filename) == 0:
    raise Exception("Please save your layout before running the simulation")


  # get target name and test file as type(str)
  gui = INTC_GUI()
  target_name, test_name = gui.ok()
  print("Selected file: " + test_name)
  print("Selected target: "+ target_name)

  # *** todo    
  #   Add the "disconnected" component to all disconnected pins
  #  optical_waveguides, optical_components = terminate_all_disconnected_pins()
  simulate = simulate
  # check if there is any cell selection
  selection = select_instances()


  #check if the selection is a circuit or not. 
  l = 0
  for i in range(len(selection)):
    if selection[i].layout().cell(selection[i].cell_index()).is_library_cell():
      l=l+1

  if selection and l == 0:    #selection is circuit from top cell
    selected_circuits_netlists_ansys(verbose=False, simulate = simulate)

  else:
    # Output the Spice netlist:
    text_Spice, text_Spice_main, num_detectors, detector_list, optical_pins = \
      topcell.spice_netlist_export_ansys(verbose=verbose, opt_in_selection_text=opt_in_selection_text)
    if not text_Spice:
      raise Exception("No netlist available. Cannot run simulation.")
      return
    if verbose:   
      print(text_Spice)
      print("optical pins: %s" %(optical_pins))
    if has_duplicates(optical_pins.split()):
      raise Exception("Optical ports with the same names are found. Please rename your optical ports so each optical port has an unique name!")


    circuit_name = topcell.name.replace('.','') # remove "."
    if '_' in circuit_name[0]:
      circuit_name = ''.join(circuit_name.split('_', 1))  # remove leading _
    
    import tempfile
    tmp_folder = tempfile.mkdtemp()
    import os
    filename = os.path.join(tmp_folder, '%s_main.spi' % circuit_name)
    filename_subckt = os.path.join(tmp_folder,  '%s.spi' % circuit_name)
    filename2 = os.path.join(tmp_folder, '%s.lsf' % circuit_name)
    filename_icp = os.path.join(tmp_folder, '%s.icp' % circuit_name)
    
    text_Spice_main += '.INCLUDE "%s"\n\n' % (filename_subckt)
    
    # Write the Spice netlist to file
    #file = open(filename, 'w')
    #file.write (text_Spice_main)
    #file.close()
    #Ansys: let user choose the directory where the netlist is saved if export netlist is selected
    if simulate == False:
      file_dir = pya.FileDialog.ask_save_file_name('Please specify the save diretory for your netlist','','.spi')
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
        

    #Ansys: if simulate means if import netlist to INTC, this will save netlist to working dir and directly opens INTC and import netlist

    if simulate: 
      # find project dir
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
      
      # get the temp netlist path
      filename_subckt = os.path.join(path, '%s.spi' % circuit_name)
      file = open(filename_subckt, 'w')
      file.write (text_Spice)
      file.close()
      import os
      netlist_path = r'%s' % filename_subckt  
      a= "net='%s';"%netlist_path

      '''
      # ask if import lsf for testbench setup
      p = pya.FileDialog.get_open_file_name("Optinal: Import an lsf script file to INTERCONNECT", "", "All Files (*)")
      lsf = p.to_s()
      print(lsf)
      '''
      gui = INTC_GUI()
      returnCode=gui.exec_()
      btn = gui.btnListener()
      if btn == "ok":
        target_name, test_name = gui.ok()
        print("Selected file: " + test_name)
        print("Selected target: "+ target_name)
        print("button:" + btn)

        # if a file selected, check if there is a target_name
        if test_name != '':

          #trim the last line of the netlist in order to import it into a compound
          with open(filename_subckt, 'r+') as fp:
              # read an store all lines into list
              lines = fp.readlines()
              # move file pointer to the beginning of a file
              fp.seek(0)
              # truncate the file
              fp.truncate()
              # start writing lines except the last line
              # lines[:-1] from line 0 to the third last line
              fp.writelines(lines[:-2])     


          # if test file is .lsf
          # Case NO.1.1: if lsf + target
          if test_name.endswith(".lsf"):
            netlist_command = "importnetlist('COMPOUND_1','%s');" %filename_subckt
            lsf_command = "lsf='%s';"%test_name 
            _globals.INTC_1 = lumapi.open('interconnect')
            lumapi.evalScript(_globals.INTC_1, "createcompound;")
            lumapi.evalScript(_globals.INTC_1, "%s" %netlist_command)
            lumapi.evalScript(_globals.INTC_1, "%s" %lsf_command)
            lumapi.evalScript(_globals.INTC_1, "feval(lsf);")
          # if test file is .icp
          if test_name.endswith(".icp"):
            # Case NO.1.2: if icp + target
            if target_name != 'Default' and target_name != '':
              icp_command = "icp='%s';"%test_name  
              netlist_command = "importnetlist('%s','%s');" %(target_name,filename_subckt) 
              _globals.INTC_1 = lumapi.open('interconnect')
              lumapi.evalScript(_globals.INTC_1, "%s" %icp_command)
              lumapi.evalScript(_globals.INTC_1, "load(icp);")
              lumapi.evalScript(_globals.INTC_1, "switchtodesign;")
              lumapi.evalScript(_globals.INTC_1, "%s" %netlist_command)
            # Case NO.2:  if there is no target_name but only file name
            else:
              pya.MessageBox.warning(
                  "Warning", "Please provide a correct compound name!", pya.MessageBox.Ok)
              return          
        # Case NO.3: if no file or target selected, proceed with ''import to INTC'
        else:
          INTC_commandline(netlist_path)

      else:
        print("button:" +btn)
        return

    if verbose:
      print('Done Lumerical INTERCONNECT circuit simulation.')
