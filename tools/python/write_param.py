"""
This script provides a menu-driven interface to write out slug
parameter files.
"""

import sys
import curses
import curses.textpad
if sys.version_info[0] < 3:
    # Python 2 version
    from curses.wrapper import wrapper
else:
    # Python 3 version
    from curses import wrapper
import os.path as osp
from collections import OrderedDict
import textwrap

# An enum to represent types of menu items
class ItemType:
    Str, Bool, Real, PosInt, PosReal, Vals, PosRealList = range(7)

# This class defines a menu option
class Menu_opt(object):
    def __init__(self, menu_str, start_val, menu_loc, 
                 param_type, allowed_vals = None,
                 optional = False, required_if = None, 
                 only_if_num = None, only_if_not = None,                  
                 only_if = None, only_if_empty = None,
                 trigger_redraw = False, help_txt = None):
        self.menu_str = menu_str
        self.val = start_val
        self.menu_loc = menu_loc
        self.param_type = param_type
        if self.param_type == ItemType.Vals:
            self.allowed_vals = allowed_vals
        self.optional = optional
        self.required_if = required_if
        self.only_if_num = only_if_num
        self.only_if_not = only_if_not
        self.only_if = only_if
        self.only_if_empty = only_if_empty
        self.trigger_redraw = True
        self.help_txt = help_txt

class Menu(object):

    # Initialization
    def __init__(self, stdscr):

        # Store stdscr
        self.stdscr = stdscr

        # Get the window size; make sure it is large enough
        self.winy, self.winx = self.stdscr.getmaxyx()
        if self.winy < 27 or self.winx < 80:
            raise RuntimeError("Terminal must be at least 80 x 27")

        # Set up the text wrapper utility class
        self.wrapper = textwrap.TextWrapper(
            width = self.winx-4, initial_indent='* ',
            subsequent_indent='  ')

        # Set the location of the options column
        self.offset = 0.5

        # Store default parameter values and menu keys
        self.params = OrderedDict([
            ('file_name', Menu_opt(
                'Parameter file', 'default.param',
                'Basic', ItemType.Str, help_txt = 
                'name of parameter file to be written')),
            ('model_name', Menu_opt(
                'Model name', 'SLUG_DEF', 'Home', ItemType.Str,
                help_txt = 
                'name of model (base name of SLUG output files)')),
            ('out_dir', Menu_opt(
                'Output directory', '', 'Basic', 
                ItemType.Str, optional = True, help_txt = 
                'directory in which to place SLUG output files;'
                ' default = current working directory')),
            ('verbosity', Menu_opt(
                'Verbosity level', '1', 'Basic', ItemType.Vals,
                allowed_vals = ['0', '1', '2'], help_txt = 
                '0 = SLUG runs silently, 1 = SLUG produces '
                + 'moderate output, 2 = SLUG produces verbose '
                + 'output')),
            ('sim_type', Menu_opt(
                'Simulation Type', 'galaxy', 'Control', ItemType.Vals,
                allowed_vals = ['galaxy', 'cluster'],
                trigger_redraw = True, help_txt =
                'type of stellar popualtion to simulate; cluster = ' +
                'simple stellar population, galaxy = '+
                'composite stellar population')),
            ('n_trials', Menu_opt(
                'Number of Trials', 1, 'Control', ItemType.PosInt,
                help_txt = 'number of Monte Carlo simulations to run')),
            ('log_time', Menu_opt(
                'Logarithmic Timestepping', False, 'Control', 
                ItemType.Bool, only_if_num = ['time_step'],
                only_if_empty = ['output_times'],
                trigger_redraw = True,
                help_txt =
                'true = logarithmically-spaced outputs, false = ' +
                'linearly-spaced outputs')),
            ('time_step', Menu_opt(
                'Time Step [yr, dex, or file name]', '1.0e6',
                'Control', ItemType.Str, trigger_redraw = True,
                only_if_empty = ['output_times'],
                help_txt = 'time step in yr (for linear time step) ' +
                'or dex (for logarithmic time step), or name of PDF file')),
            ('start_time', Menu_opt(
                'Start Time [yr]', '', 'Control', ItemType.PosReal,
                only_if_num = ['time_step'],
                only_if_empty = ['output_times'],
                optional = True, required_if = [('log_time', True)],
                help_txt = 'time of first output; if left ' +
                'unspecified, and time stepping is linear, defaults '+
                'to the time step')),
            ('end_time', Menu_opt(
                'End Time [yr]', '1.0e8', 'Control', ItemType.PosReal,
                only_if_num = ['time_step'],
                only_if_empty = ['output_times'],
                help_txt = 'time at which to end simulation')),
            ('output_times', Menu_opt(
                'Output Times [yr]', '', 'Control', ItemType.PosRealList,
                optional=True, trigger_redraw = True,
                help_txt = 'comma-separated list '
                'of exact output times, in yr; if set, overrides '
                'all other timestepping options')),
            ('sfr', Menu_opt(
                'Star Formation Rate [Msun/yr or "sfh"]', '1.0',
                'Control', ItemType.Str,
                only_if = [('sim_type', 'galaxy')],
                help_txt = 'numerical values give ' +
                'constant SFR; "sfh" causes the SFR to be read ' +
                'from a file')),
            ('sfh', Menu_opt(
                'Star Formation History', '', 'Control', 
                ItemType.Str, optional = True,
                only_if = [('sim_type', 'galaxy'),
                           ('sfr', 'sfh')],
                help_txt = 'name of a PDF file specifying the star ' +
                'formation history')),
            ('cluster_mass', Menu_opt(
                'Cluster Mass [Msun or file name]', '1.0e3', 'Control', 
                ItemType.Str, only_if=[('sim_type', 'cluster')], help_txt =
                'numerical values give the ' +
                'cluster mass; non-numerical values specify a '
                'PDF file')), 
            ('redshift', Menu_opt(
                'Redshift', '', 'Control', ItemType.PosReal,
                optional=True, help_txt = 
                'redshift of the system; all output spectra and ' +
                'photometry will be computed in the observed frame;' +
                ' defaults to 0')),
            ('out_integrated', Menu_opt(
                'Write integrated properties', True, 'Output', 
                ItemType.Bool, only_if = [('sim_type', 'galaxy')], 
                help_txt = 'if true, write out physical properties ' +
                'of whole galaxy')),
            ('out_integrated_phot', Menu_opt(
                'Write integrated photometry', True, 'Output', 
                ItemType.Bool, only_if = [('sim_type', 'galaxy')],
                help_txt = 'if true, write out photometric ' +
                'properties of galaxy as a whole')),
            ('out_integrated_sn', Menu_opt(
                'Write integrated SN numbers', True,
                'Output', ItemType.Bool,
                only_if = [('sim_type', 'galaxy')],
                help_txt = 'if true, write out number of SNe in ' +
                'galaxy as a whole')),
            ('out_integrated_spec', Menu_opt(
                'Write integrated spectra', True,
                'Output', ItemType.Bool,
                only_if = [('sim_type', 'galaxy')],
                help_txt = 'if true, write out spectra of ' +
                'galaxy as a whole')),
            ('out_integrated_yield', Menu_opt(
                'Write integrated yield', True,
                'Output', ItemType.Bool,
                only_if = [('sim_type', 'galaxy')],
                help_txt = 'if true, write out chemical yield of ' +
                'galaxy as a whole')),
            ('out_cluster', Menu_opt(
                'Write cluster properties', True,
                'Output', ItemType.Bool,
                help_txt = 'if true, write out physical ' +
                'properties of each star cluster')),
            ('out_cluster_phot', Menu_opt(
                'Write cluster photometry', True,
                'Output', ItemType.Bool,
                help_txt = 'if true, write out photometric ' +
                'properties of each star cluster')),
            ('out_cluster_sn', Menu_opt(
                'Write cluster SN numbers', True,
                'Output', ItemType.Bool,
                help_txt = 'if true, write out number of SN in '+
                'eacsh cluster')),
            ('out_cluster_spec', Menu_opt(
                'Write cluster spectra', True,
                'Output', ItemType.Bool,
                help_txt = 'if true, write out spectra ' +
                'of each star cluster')),
            ('out_cluster_yield', Menu_opt(
                'Write cluster yields', True,
                'Output', ItemType.Bool,
                help_txt = 'if true, write out chemical yields ' +
                'of each star cluster')),
            ('output_mode', Menu_opt(
                'Output format', 'ascii', 'Output', ItemType.Vals,
                allowed_vals = ['ascii', 'binary', 'fits'],
                help_txt = 'format of output files: ASCII, binary, '+
                'or FITS')),
            ('imf', Menu_opt(
                'Initial mass function', 
                osp.join('lib', 'imf', 'chabrier.imf'),
                'Stellar', ItemType.Str,
                help_txt = 'name of file specifying the stellar ' +
                'IMF')),
            ('cmf', Menu_opt(
                'Cluster mass function', 
                osp.join('lib', 'cmf', 'slug_default.cmf'),
                'Stellar', ItemType.Str, help_txt = 
                'name of file specifying cluster ' +
                'mass function')),
            ('clf', Menu_opt(
                'Cluster lifetime function', 
                osp.join('lib', 'clf', 'slug_default.clf'),
                'Stellar', ItemType.Str, help_txt = 
                'file giving the cluster ' +
                'lifetime distribution')),
            ('tracks', Menu_opt(
                'Stellar evolution tracks', 
                'geneva_2013_vvcrit_00',
                'Stellar', ItemType.Str, help_txt =
                'name of the stellar evolution ' +
                'track set or file')),
            ('atmospheres', Menu_opt(
                'Stellar atmosphere directory',
                osp.join('lib', 'atmospheres')+osp.sep,
                'Stellar', ItemType.Str, help_txt = 
                'directory containing the stellar ' +
                'atmospheres')),
            ('specsyn_mode', Menu_opt(
                'Spectral synthesis mode', 'SB99',
                'Stellar', ItemType.Vals,
                allowed_vals = ['planck', 'Kurucz', 'Kurucz+Hillier',
                                'Kurucz+Pauldrach', 'SB99'],
                help_txt = 'method to use for stellar spectrum '+
                'calculation')),
            ('clust_frac', Menu_opt(
                'Cluster fraction', '1.0', 'Stellar', ItemType.Real,
                help_txt = 'fraction of stars formed in clusters')),
            ('min_stoch_mass', Menu_opt(
                'Minimum stochastic mass [Msun]', '0.0',
                'Stellar', ItemType.Real, help_txt = 
                'minimum star mass to treat stochastically')),
            ('metallicity', Menu_opt(
                'Metallicity [Zsun]', '',
                'Stellar', ItemType.PosReal, optional=True,
                help_txt = 'metallicity normalized to Solar')),
            ('A_V', Menu_opt(
                'A_V [mag or file name]',
                '',
                'Extinction', ItemType.Str, optional=True,
                trigger_redraw=True, help_txt =
                'visual extinction; numeric values give a constant '+
                'extinction, non-numeric values give a file name '+
                'containing the PDF of extinction; if omitted, no '+
                'extinction is included')),
            ('extinction_curve', Menu_opt(
                'Extinction curve', osp.join('lib', 'extinct',
                                             'SB_ATT_SLUG.dat'),
                'Extinction', ItemType.Str, optional=True,
                only_if_not = [('A_V', '')], help_txt =
                'file containing the shape of the extinction curve')),
            ('nebular_extinction_factor', Menu_opt(
                'A_V,neb/A_V [num or name]',
                '',
                'Extinction', ItemType.Str, optional=True,
                only_if_not = [('A_V', '')],
                help_txt =
                'ratio of nebular to stellar extinction; numeric values ' +
                'give a constant '+
                'extinction, non-numeric values give a file name '+
                'containing the PDF of extinction; if omitted, no '+
                'extinction is included')),
            ('compute_nebular', Menu_opt(
                'Compute nebular emission', True,
                'Nebular', ItemType.Bool, help_txt =
                'true = compute nebular contribution for clusters < '+
                '10 Myr old; false = do not comput nebular ' +
                'contribution')),
            ('atomic_data', Menu_opt(
                'Atomic data directory', 
                osp.join('lib', 'atomic')+osp.sep,
                'Nebular', ItemType.Str,
                only_if = [('compute_nebular', True)], help_txt =
                'directory containing atomic data used in nebular '+
                'computations')),
            ('nebular_no_metals', Menu_opt(
                'Omit nebular metal lines', False,
                'Nebular', ItemType.Bool, 
                only_if = [('compute_nebular', True)], help_txt =
                'true = treat nebula as pure H; false = nebular of '+
                'scaled Solar composition')),
            ('nebular_den', Menu_opt(
                'Nebular density [cm^-3]', '1.0e2',
                'Nebular', ItemType.PosReal,
                only_if = [('compute_nebular', True)], help_txt =
                'density of nebula to be used in computing '+
                'continuum emission')),
            ('nebular_temp', Menu_opt(
                'Nebular temp. [K, or -1 for auto]',
                '-1', 'Nebular', ItemType.Real,
                only_if = [('compute_nebular', True)], help_txt =
                'temperature of nebular to be used in computing '+
                'continuum emission; if set to a value < 0, '+
                'temperature is determined automatically from the ' +
                'value of log U')),
            ('nebular_logU', Menu_opt(
                'log ionization parameter (U)',
                '-3', 'Nebular', ItemType.Vals,
                allowed_vals = ['-3', '-2.5', '-2'],
                only_if = [('compute_nebular', True)], help_txt =
                'value of log U to use in computing the nebular '+
                'emission')),
            ('nebular_phi', Menu_opt(
                'phi (ioniz. photon absorption frac.)',
                '0.73', 'Nebular', ItemType.Real,
                only_if = [('compute_nebular', True)], help_txt =
                'fraction of ionizing photons absorbed by H ' +
                'within the observational aperture and thus ' +
                'able to contribute to nebular emission')),
            ('phot_mode', Menu_opt(
                'Photometry output mode', 'L_nu',
                'Photometry', ItemType.Vals,
                allowed_vals = ['L_nu', 'L_lambda', 'AB', 'STMAG',
                                'VEGA'], help_txt = 
                'photometric system to use for output; L_nu = '+
                'frequency-averaged luminosity (erg/s/Hz), '+
                'L_lambda = wavelength-averaed luminosity '+
                '(erg/s/Ang), AB = AB magnitudes, STMAG = '+
                'ST magnitudes, VEGA = Vega magnitudes')),
            ('filters', Menu_opt(
                'Filter library directory', 
                osp.join('lib', 'filters') + osp.sep,
                'Photometry', ItemType.Str, help_txt = 
                'directory containing photometric filter data')),
            ('phot_bands', Menu_opt(
                'Photometric bands', 'Lbol, QH0',
                'Photometry', ItemType.Str, help_txt =
                'comma-separated list of filters / bands for ' +
                'which photometry should be computed; see ' +
                osp.join('lib', 'filters', 'FILTER_LIST') +
                ' for a full list of available filters')),
            ('yield_dir', Menu_opt(
                'Yields directory',
                osp.join('lib', 'yields')+osp.sep,
                'Yields', ItemType.Str,
                help_txt =
                'directory containing yield tables')),
            ('yield_mode', Menu_opt(
                'Yield computation mode',
                'sukhbold16+karakas16+doherty14',
                'Yields',
                ItemType.Vals,
                allowed_vals = [
                    'sukhbold16+karakas16+doherty14',
                    'karakas16+doherty14',
                    'sukhbold16'],
                help_txt = 'yield tables to use')),
            ('no_decay_isotopes', Menu_opt(
                'Disable isotope decay',
                False,
                'Yields', ItemType.Bool,
                help_txt = 'if true, no radioactive decay is'
                ' applied to isotope masses')),
            ('isotopes_included', Menu_opt(
                'Isotopes included',
                'intersection',
                'Yields',
                ItemType.Vals,
                allowed_vals = [ 'intersection', 'union' ],
                help_txt = 'when combining yield tables, include '
                'only isotopes present in all tables (intersection)'
                ' or present in any table (union)'))
        ])

        # List of menus and names
        self.menus = OrderedDict([
            ('Home', 'Main menu'),
            ('Basic', 'Basic keywords'),
            ('Control', 'Simulation control keywords'),
            ('Output', 'Output control keywords'),
            ('Stellar', 'Stellar model keywords'),
            ('Extinction', 'Extinction keywords'),
            ('Nebular', 'Nebular emission keywords'),
            ('Photometry', 'Photometry keywords'),
            ('Yields', 'Yield keywords')
        ])

        # Hide the cursor
        curses.curs_set(0)

        # Draw the home menu
        self.draw_menu('Home')
        self.show_help()

        # Main loop to read user input
        self.done = False
        while not self.done:
            c = self.stdscr.getch()
            if c == curses.KEY_DOWN:
                self.navigate_updown(1)
                self.show_help()
            elif c == curses.KEY_UP:
                self.navigate_updown(-1)
                self.show_help()
            elif c == curses.KEY_RIGHT:
                self.navigate_leftright(1)
                self.show_help()
            elif c == curses.KEY_LEFT:
                self.navigate_leftright(-1)
                self.show_help()
            elif c == curses.KEY_HOME:
                self.reset_cursor()
                self.show_help()
            elif c == ord('\n'):
                self.navigate_select()


    # Function to write out a parameter file
    def write_file(self):

        # Open file
        fp = open(self.params['file_name'].val, 'w')
        fp.write('# SLUG2 parameter file auto-generated by write_param\n')
        for k in self.params.keys():

            # Skip parameters that should be omitted
            if k == 'file_name':
                continue
            if self.params[k].only_if is not None:
                include = True
                for con in self.params[k].only_if:
                    if self.params[con[0]].val != con[1]:
                        include = False
                if include == False:
                    continue
            if self.params[k].only_if_not is not None:
                include = True
                for con in self.params[k].only_if_not:
                    if self.params[con[0]].val == con[1]:
                        include = False
                if include == False:
                    continue
            if self.params[k].only_if_num is not None:
                include = True
                for con in self.params[k].only_if_num:
                    try:
                        numval = float(self.params[con].val)
                    except ValueError:
                        include = False
                if include == False:
                    continue
            if self.params[k].only_if_empty is not None:
                include = True
                for con in self.params[k].only_if_empty:
                    if len(self.params[con].val) > 0:
                        include = False
                if include == False:
                    continue

            # Skip optional parameters that have been left blank
            if self.params[k].optional and self.params[k].val == '':
                continue

            # Write
            fp.write(k + '     ')
            if (self.params[k].param_type == ItemType.Str) or \
               (self.params[k].param_type == ItemType.Vals) or \
               (self.params[k].param_type == ItemType.PosReal) or \
               (self.params[k].param_type == ItemType.Real) or \
               (self.params[k].param_type == ItemType.PosRealList):
                fp.write(self.params[k].val + '\n')
            elif self.params[k].param_type == ItemType.PosInt:
                fp.write(str(self.params[k].val) + '\n')
            elif self.params[k].param_type == ItemType.Bool:
                if self.params[k].val == True:
                    fp.write('1\n')
                else:
                    fp.write('0\n')

        # Close and report status
        fp.close()
        self.message('Wrote '+self.params['file_name'].val)


    # Function to handle when the user hits return
    def navigate_select(self):

        # Get current position
        posy, posx = self.stdscr.getyx()

        # Are we in the left column?
        if posx == 1:

            # Is this a command line or an option line? If the former,
            # execute; if the latter, go to the right
            if posy in self.opt_lines:
                self.navigate_leftright(1)
            else:
                com = self.stdscr.instr(
                    int(self.offset*self.winx)-1).strip()
                if sys.version_info[0] > 2:
                    # Deal with python 3 string bullshit
                    com = com.decode('ascii')
                if com == 'Exit without writing':
                    self.done = True
                elif com == 'Write and exit':
                    self.write_file()
                    self.done = True
                elif com == 'Write':
                    self.write_file()
                else:
                    for k, v in zip(self.menus.keys(),
                                    self.menus.values()):
                        if v == com:
                            self.draw_menu(k)
                            break

        else:

            # We're in the right column, on an option. See what type
            # of option this is.
            p = self.cur_menu_params[self.opt_lines.index(posy)]

            # Act based on the type of option this is
            old_val = p.val
            if (p.param_type == ItemType.Str) or \
               (p.param_type == ItemType.PosInt) or \
               (p.param_type == ItemType.PosReal) or \
               (p.param_type == ItemType.Real) or \
               (p.param_type == ItemType.PosRealList):

                # String or number option: go into text editing mode

                # Set up new window to hold textpad
                ncols = self.winx-int(self.offset*self.winx)-1
                win = curses.newwin(1, ncols,
                                    posy, int(self.offset*self.winx))
                win.addstr(str(old_val))
                tb = curses.textpad.Textbox(win)

                # Turn cursor back on and edit; make sure that the
                # user has entered a valid value before allowing them
                # to exit editing
                curses.curs_set(1)
                edit_done = False
                while not edit_done:
                    newtxt = tb.edit().strip()
                    edit_done = True
                    if len(newtxt) == 0 and not p.optional:
                        edit_done = False
                    elif p.param_type == ItemType.PosInt:
                        try:
                            numval = int(newtxt)
                        except ValueError:
                            numval = 0
                        if numval <= 0:
                            self.message('Error: ' + p.menu_str + 
                                         ' must be a positive' +
                                         ' integer')
                            edit_done = False
                    elif p.param_type == ItemType.PosReal:
                        try:
                            numval = float(newtxt)
                        except ValueError:
                            numval = 0
                            newtxt = ''
                        if numval <= 0 and not p.optional:
                            self.message('Error: ' + p.menu_str + 
                                         ' must be a positive' +
                                         ' real number')
                            edit_done = False
                    elif p.param_type == ItemType.Real:
                        try:
                            numval = float(newtxt)
                        except ValueError:
                            newtxt = ''
                            if not p.optional:
                                self.message('Error: ' + p.menu_str + 
                                             ' must be a real number')
                                edit_done = False
                    elif p.param_type == ItemType.PosRealList:
                        try:
                            listval = newtxt.split(',')
                            if len(listval) == 0:
                                raise ValueError
                            for i in range(len(listval)):
                                listval[i] = float(listval[i])
                        except ValueError:
                            newtxt = ''
                            if not p.optional:
                                self.message(
                                    'Error: ' + p.menu_str +
                                    ' must be a comma-separated'
                                    ' list of positive real numbers')
                                edit_done = False
                if p.param_type == ItemType.PosInt:
                    p.val = numval
                elif (p.param_type == ItemType.Str) or \
                     (p.param_type == ItemType.PosReal) or \
                     (p.param_type == ItemType.Real) or \
                     (p.param_type == ItemType.PosRealList):
                    if len(newtxt) == 0:
                        newtxt = ''
                    p.val = newtxt

                # Turn cursor off and redraw
                curses.curs_set(0)
                self.clear_message()
                self.stdscr.addnstr(posy, posx, newtxt, len(newtxt),
                                    curses.A_REVERSE)
                if len(newtxt) < ncols:
                    self.stdscr.addnstr(posy, posx+len(newtxt),
                                        ' '*(ncols-len(newtxt)),
                                        ncols-len(newtxt))

            elif p.param_type == ItemType.Vals:

                # List of values -- just cycle through the list
                curidx = p.allowed_vals.index(p.val)
                curidx = (curidx + 1) % len(p.allowed_vals)
                p.val = p.allowed_vals[curidx]
                self.stdscr.addnstr(posy, posx, p.val, len(p.val),
                                    curses.A_REVERSE)
                nspace = self.winx - int(self.offset*self.winx) \
                         - len(p.val) - 1
                if nspace > 0:
                    self.stdscr.addnstr(posy, posx+len(p.val), 
                                        ' '*nspace, nspace)

            # Move cursor back and redraw screen
            self.stdscr.move(posy, posx)
            self.stdscr.refresh()

            # Redraw if required
            if old_val != p.val and p.trigger_redraw:
                self.draw_menu(self.cur_menu, start_opt=p.menu_str,
                               start_right=True)

            # Put help message back
            self.show_help()

    # Function to navigate the menu left and right
    def navigate_leftright(self, dx):

        # Get current position
        posy, posx = self.stdscr.getyx()

        # If requested to move right, and we can do so, move the
        # cursor and highlighting
        if (dx > 0) and (posx == 1) and (posy in self.opt_lines):
            self.highlight(False)
            self.stdscr.move(posy, int(self.offset*self.winx))
            self.highlight(True)

        elif (dx < 0):

            # Move right
            self.highlight(False)
            self.stdscr.move(posy, 1)
            self.highlight(True)

    # Function to navigate the menu up and down
    def navigate_updown(self, dy):

        # See if we're on the left or right
        linenum, posx = self.stdscr.getyx()
        if posx == 1:

            # We're on the left

            # Unhighlight current line
            self.highlight(False)

            # Move to next line if possible
            idx = (self.selectable_lines.index(linenum) + dy) \
                  % len(self.selectable_lines)
            linenum = self.selectable_lines[idx]
            self.stdscr.move(linenum, 1)

            # Highlight new line
            self.highlight(True)

        else:

            # We're on the right. If this is a boolean, specified
            # value, or interger parameter, cycle the value
            p = self.cur_menu_params[self.opt_lines.index(linenum)]
            old_val = p.val
            if p.param_type == ItemType.Bool:
                # Cycle boolean values
                p.val = not p.val
                if p.val:
                    self.stdscr.addnstr(linenum, posx, 'True', 4,
                                        curses.A_REVERSE)
                    self.stdscr.addnstr(linenum, posx+4, ' ', 1)
                else:
                    self.stdscr.addnstr(linenum, posx, 'False', 5,
                                        curses.A_REVERSE)
            elif p.param_type == ItemType.Vals:
                # Cycle through pre-defined options
                curidx = p.allowed_vals.index(p.val)
                curidx = (curidx + dy) % len(p.allowed_vals)
                p.val = p.allowed_vals[curidx]
                self.stdscr.addnstr(linenum, posx, p.val, len(p.val),
                                    curses.A_REVERSE)
                nspace = self.winx - int(self.offset*self.winx) \
                         - len(p.val) - 1
                if nspace > 0:
                    self.stdscr.addnstr(linenum, posx+len(p.val), 
                                        ' '*nspace, nspace)
            elif p.param_type == ItemType.PosInt:
                # Add or subtract
                p.val = p.val + dy
                if p.val == 0:
                    p.val = 1
                self.stdscr.addnstr(linenum, posx, str(p.val), 
                                    len(str(p.val)), curses.A_REVERSE)
                nspace = self.winx - int(self.offset*self.winx) \
                         - len(str(p.val)) - 1
                if nspace > 0:
                    self.stdscr.addnstr(linenum, posx+len(str(p.val)), 
                                        ' '*nspace, nspace)

            # Reset cursor
            self.stdscr.move(linenum, posx)                

            # Redraw menu if needed
            if old_val != p.val and p.trigger_redraw:
                self.draw_menu(self.cur_menu, start_opt=p.menu_str,
                               start_right=True)

        # Redraw
        self.stdscr.refresh()

    # Show the help text
    def show_help(self):

        # Read the line we're on
        posy, posx = self.stdscr.getyx()

        # Show the help text if available
        if posy in self.opt_lines:
            name = self.cur_menu_params[
                self.opt_lines.index(posy)].menu_str
            txt = self.cur_menu_params[
                self.opt_lines.index(posy)].help_txt
            if txt is not None:
                self.message(name + ': ' + txt)
        elif posy in self.menu_lines:
            self.clear_message()
        elif posy in self.com_lines:
            idx = self.com_lines.index(posy)
            if idx == 0:
                self.message('Write out ' +
                             self.params['file_name'].val)
            elif idx == 1:
                self.message('Write out ' +
                             self.params['file_name'].val +
                             ' and exit')
            elif idx == 2:
                self.message('Exit without writing output')


    # Function to highlight the currently-selected option
    def highlight(self, on=True):

        # Save cursor position
        posy, posx = self.stdscr.getyx()

        # Read the current text
        cur_txt \
            = self.stdscr.instr(int(self.offset*self.winx)-1).strip()
        if len(cur_txt) == 0:
            cur_txt = ' '

        # Write it back with the requested highlighting
        if on:
            self.stdscr.addnstr(posy, posx, cur_txt, len(cur_txt),
                                curses.A_REVERSE)
        else:
            self.stdscr.addnstr(posy, posx, cur_txt, len(cur_txt))

        # Move cursor back to starting location
        self.stdscr.move(posy, posx)

    # Function to reset the cursor
    def reset_cursor(self):

        # Move cursor and highlight
        self.highlight(False)
        self.stdscr.move(4, 1)
        self.highlight(True)

    # Function to write an option
    def write_opt(self, name, default, linenum=None, attr=None):

        # Get current y position if needed
        if linenum is None:
            linenum, posx = self.stdscr.getyx()

        # Write the option name
        if attr is None:
            self.stdscr.addnstr(linenum, 1, name, len(name))
        else:
            self.stdscr.addnstr(linenum, 1, name, len(name), attr)

        # Write the default
        if default is not None:
            if attr is None:
                self.stdscr.addnstr(linenum, int(self.offset*self.winx), 
                                    str(default), 
                                    len(str(default)))
            else:
                self.stdscr.addnstr(linenum, int(self.offset*self.winx),
                                    str(default), 
                                    len(str(default)), attr)

    # Function to draw the a menu
    def draw_menu(self, name, start_opt=None, start_right=False):

        # Clear
        self.stdscr.erase()

        # Draw a border and a title
        self.stdscr.border()
        titstr = 'SLUG2 Parameter File Writer'
        self.stdscr.addnstr(1, int(self.winx/2-len(titstr)/2), 
                            titstr, len(titstr))

        # Record which lines contain options and which contain commands
        self.opt_lines = []
        self.menu_lines = []
        self.com_lines = []

        # Figure out which parameters go on this menu
        self.cur_menu = name
        self.cur_menu_params = []
        for p in self.params.values():
            if p.menu_loc == name:

                # Skip if conditions are not met
                if p.only_if is not None:
                    show = True
                    for con in p.only_if:
                        if self.params[con[0]].val != con[1]:
                            show = False
                    if show == False:
                        continue
                if p.only_if_not is not None:
                    show = True
                    for con in p.only_if_not:
                        if self.params[con[0]].val == con[1]:
                            show = False
                    if show == False:
                        continue
                if p.only_if_num is not None:
                    show = True
                    for con in p.only_if_num:
                        try:
                            numval = float(self.params[con].val)
                        except ValueError:
                            show = False
                    if show == False:
                        continue
                if p.only_if_empty is not None:
                    show = True
                    for con in p.only_if_empty:
                        if len(self.params[con].val) > 0:
                            show = False
                    if show == False:
                        continue

                # Add to menu
                self.cur_menu_params.append(p)

        # Draw the part of the menu containing options
        self.write_opt('Option', 'Value', 3, 
                  attr=curses.A_BOLD | curses.A_UNDERLINE)
        linenum = 4
        start_line = 0
        for p in self.cur_menu_params:
            menu_str = p.menu_str
            if p.optional:
                opt = True
                if p.required_if is not None:
                    for con in p.required_if:
                        if self.params[con[0]].val == con[1]:
                            opt = False
                if opt:
                    menu_str = menu_str + ' [optional]'
            self.write_opt(menu_str, p.val, linenum)
            self.opt_lines.append(linenum)
            if start_opt is not None:
                if start_opt == p.menu_str:
                    start_line = linenum
            linenum = linenum+1

        # Draw the part of the menu containing links to other menus
        linenum = linenum+1
        self.write_opt('More parameters', None, linenum,
                       attr=curses.A_BOLD | curses.A_UNDERLINE)
        linenum = linenum+1
        for m in self.menus.keys():
            if m != name:
                self.write_opt(self.menus[m], None, linenum)
                self.menu_lines.append(linenum)
                linenum = linenum+1

        # Draw the part of the menu for writing out files, if we're on
        # the home menu
        if name == 'Home':
            linenum = linenum+1
            self.write_opt('Commands', None, linenum,
                           attr=curses.A_BOLD | curses.A_UNDERLINE)
            linenum = linenum+1
            self.write_opt('Write', None, linenum)
            self.com_lines.append(linenum)
            linenum = linenum+1
            self.write_opt('Write and exit', None, linenum)
            self.com_lines.append(linenum)
            linenum = linenum+1
            self.write_opt('Exit without writing', None, linenum)
            self.com_lines.append(linenum)
            linenum = linenum+1

        # Get all lines
        self.selectable_lines = self.com_lines + self.menu_lines \
                                + self.opt_lines
        self.selectable_lines.sort()

        # Set the initial cursor line
        if start_line == 0:
            self.reset_cursor()
        else:
            if start_right:
                self.stdscr.move(start_line,
                                 int(self.offset*self.winx))
            else:
                self.stdscr.move(start_line, 1)
            self.highlight(True)            

        # Store the number of the line after the end of the text; this
        # is used to issue messages
        linenum = linenum+1
        self.msg_top_line = linenum

        # Redraw
        self.stdscr.refresh()

    # Function to write a message at the bottom of the screen
    def message(self, msg):
        self.clear_message()
        posy, posx = self.stdscr.getyx()
        msg_break = self.wrapper.wrap(msg)
        for i, m in enumerate(msg_break):
            self.stdscr.addnstr(self.winy-1-len(msg_break)+i, 
                                2, m, len(m))
        self.stdscr.move(posy, posx)
        self.stdscr.refresh()

    # Function to clear the message at the bottom of the screen
    def clear_message(self):
        posy, posx = self.stdscr.getyx()
        for i in range(self.msg_top_line, self.winy-1):
            self.stdscr.addnstr(i, 1, ' '*(self.winx-2), 
                                self.winx-2)
        self.stdscr.move(posy, posx)
        self.stdscr.refresh()


# Main function to drive the menus
def runwindow(stdscr):

    menu = Menu(stdscr)

# Run the function
wrapper(runwindow)
