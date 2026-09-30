import os
import numpy as np
from CSXCAD.CSProperties import CSPropProbeBox, CSPropDumpBox
from CSXCAD.Utilities import CheckNyDir
from openEMS import utilities
def _load_ui_file(filepath):
    """Read an openEMS probe output file in a single pass.

    Returns (data, col_names) where data is a float64 array of shape (N, ncols)
    and col_names is the list of column-name strings from the last header line
    (e.g. ['t/s', 'voltage', 'mode_purity']), or None if no header was found.
    """
    comments = []
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('%'):
                comments.append(line[1:].strip())
            else:
                s = line.strip()
                if s:
                    rows.append(s.split())
    col_names = comments[-1].split() if comments else None
    if col_names is not None and col_names[0] != 't/s':
        raise ValueError('{}: unexpected first column "{}", expected "t/s"'.format(filepath, col_names[0]))
    return np.array(rows, dtype=np.double), col_names
class UI_data:
    def __init__(self, fns, path, freq, signal_type='pulse', **kw):
        self.path = path
        if type(fns)==str:
            fns = [fns]
        self.fns  = fns

        if np.isscalar(freq):
            freq = [freq]
        self.freq = freq

        self.ui_time        = []
        self.ui_val         = []
        self.ui_f_val       = []
        self.col_names      = []  # column names per file, from the file header
        self.ui_mode_purity = []  # mode-purity time series or None

        for fn in fns:
            data, col_names = _load_ui_file(os.path.join(path, fn))
            self.col_names.append(col_names)
            self.ui_time.append(data[:, 0])
            self.ui_val.append(data[:, 1])
            self.ui_f_val.append(utilities.DFT_time2freq(data[:, 0], data[:, 1], freq, signal_type=signal_type))
            has_purity = (col_names is not None and len(col_names) > 1
                          and col_names[-1] == 'mode_purity')
            self.ui_mode_purity.append(data[:, -1] if has_purity else None)
class Port(object):
    """
    The port base class.

    :param CSX: Continuous Structure
    :param port_nr: int -- port number, must be unique among all ports with the
        same PortNamePrefix, as openEMS writes each port probe to a file named
        after the port number
    :param R: float -- port reference impedance, e.g. 50 (Ohms)
    :param start, stop: (3,) array -- Start/Stop box coordinates
    :param p_dir: int -- port direction
    :param excite: float -- port excitation amplitude
    :param priority: int -- priority of all contained primtives
    :param PortNamePrefix: str -- a prefix for all ports-names
    :param delay: float -- a positive delay value to e.g. emulate a phase shift
    """
    def __init__(self, CSX, port_nr, start, stop, excite, **kw):
        self.CSX      = CSX
        self.number   = port_nr
        self.excite   = excite
        self.start    = np.array(start, np.double)
        self.stop     = np.array(stop, np.double)
        self.Z_ref    = None
        self.U_filenames = kw.get('U_filenames', [])
        self.I_filenames = kw.get('I_filenames', [])
        self.port_props = []

        self.priority = 0
        if 'priority' in kw:
            self.priority = kw['priority']

        self.prefix = ''
        if 'PortNamePrefix' in kw:
            self.prefix = kw['PortNamePrefix']
        self.delay = 0

        if 'delay' in kw:
            self.delay = kw['delay']

        self.lbl_temp = self.prefix + 'port_{}' +  '_{}'.format(self.number)

    def _AddProbe(self, CSX, name, **kw):
        # openEMS writes each probe to a file of its name, two probes with the
        # same name would corrupt each other's file
        for prop in CSX.GetPropertiesByName(name):
            if isinstance(prop, CSPropProbeBox) and not isinstance(prop, CSPropDumpBox):
                raise ValueError('port {}: a probe named "{}" already exists, port numbers '
                                 'must be unique (or use a different PortNamePrefix)'.format(self.number, name))
        return CSX.AddProbe(name, **kw)

    def SetEnabled(self, val):
        from CSXCAD.CSProperties import CSPropExcitation
        found_any = False
        for prop in self.port_props:
            if type(prop) == CSPropExcitation:
                prop.SetEnabled(val)
                found_any = True
        if not found_any and val:
            # if we attempt to activate this port and it does not have any excitation set, raise an exception!
            raise Exception('Unable to enable port! No excitation found!')

    def ReadUIData(self, sim_path, freq, signal_type ='pulse'):
        self.u_data = UI_data(self.U_filenames, sim_path, freq, signal_type )
        self.uf_tot = 0
        self.ut_tot = 0
        for n in range(len(self.u_data.fns)):
            self.uf_tot += self.u_data.ui_f_val[n]
            self.ut_tot += self.u_data.ui_val[n]
        self.u_time = self.u_data.ui_time[0]

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type )
        self.if_tot = 0
        self.it_tot = 0
        for n in range(len(self.i_data.fns)):
            self.if_tot += self.i_data.ui_f_val[n]
            self.it_tot += self.i_data.ui_val[n]
        self.i_time = self.i_data.ui_time[0]

        # mode purity: extra column written by ProcessModeMatch (index 0 = purity)
        self.u_mode_purity = self.u_data.ui_mode_purity
        self.i_mode_purity = self.i_data.ui_mode_purity


    def CalcPort(self, sim_path, freq, ref_impedance=None, ref_plane_shift=None, signal_type='pulse'):
        self.ReadUIData(sim_path, freq, signal_type)

        if ref_impedance is not None:
            self.Z_ref = ref_impedance
        if self.Z_ref is None:
            raise Exception('Port Z_ref should not be None!')

        if ref_plane_shift is not None:
            if not hasattr(self, 'beta'):
                raise Exception('Port has no beta attribute!')
            shift = ref_plane_shift
            if self.measplane_shift:
                shift -= self.measplane_shift
            shift *= self.CSX.GetGrid().GetDeltaUnit()
            phase = np.real(self.beta)*shift
            uf_tot = self.uf_tot * np.cos(-phase) + 1j * self.if_tot * self.Z_ref * np.sin(-phase)
            if_tot = self.if_tot * np.cos(-phase) + 1j * self.uf_tot / self.Z_ref * np.sin(-phase)
            self.uf_tot = uf_tot
            self.if_tot = if_tot

        self.uf_inc = 0.5 * ( self.uf_tot + self.if_tot * self.Z_ref )
        self.if_inc = 0.5 * ( self.if_tot + self.uf_tot / self.Z_ref )
        self.uf_ref = self.uf_tot - self.uf_inc
        self.if_ref = self.if_inc - self.if_tot

        if type(self.Z_ref) in [int, float]:
            self.ut_inc = 0.5 * ( self.ut_tot + self.it_tot * self.Z_ref )
            self.it_inc = 0.5 * ( self.it_tot + self.ut_tot / self.Z_ref )
            self.ut_ref = self.ut_tot - self.ut_inc
            self.it_ref = self.it_inc - self.it_tot

        # calc some more port parameter
        # incoming power
        self.P_inc = 0.5*np.real(self.uf_inc*np.conj(self.if_inc))
        # reflected power
        self.P_ref = 0.5*np.real(self.uf_ref*np.conj(self.if_ref))
        # accepted power (incoming - reflected)
        self.P_acc = 0.5*np.real(self.uf_tot*np.conj(self.if_tot))
class CPWPort(Port):
    """Coplanar waveguide (CPW) port.

    Creates the CPW metal, gap-spanning voltage probes (left and right gaps),
    current probes, and an optional symmetric excitation across the two gaps.

    Parameters
    ----------
    metal_prop : CSProperties
        Metal property for the CPW conductor.
    prop_dir : int or str
        Direction of wave propagation (0/1/2 or 'x'/'y'/'z').
    exc_dir : int or str
        E-field direction across the gaps (0/1/2 or 'x'/'y'/'z'), i.e. the
        width direction of the CPW. The CPW plane is normal to the cross
        product of ``prop_dir`` and ``exc_dir``.
    gap_width : float
        Width of each CPW gap in drawing units.
    excite : bool or float, optional
        Enable port excitation.
    FeedShift : float, optional
        Excitation shift from ``start``.
    Feed_R : float, optional
        Lumped resistance in Ohms (applied to each gap as 2*R).
    MeasPlaneShift : float, optional
        Measurement plane distance from ``start``.

    See Also
    --------
    Port, MSLPort, StripLinePort
    """

    def __init__(self, CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir,
                 gap_width, excite=0, **kw):
        super(CPWPort, self).__init__(CSX, port_nr=port_nr, start=start, stop=stop, excite=excite, **kw)

        self.prop_ny = CheckNyDir(prop_dir)

        # Width direction = E-field direction across the gaps; height direction
        # (normal of the CPW plane) = cross product, as in AddCPWPort.m
        self.width_ny = CheckNyDir(exc_dir)
        if self.width_ny == self.prop_ny:
            raise Exception('CPWPort: exc_dir must differ from prop_dir')
        self.height_ny = 3 - self.prop_ny - self.width_ny

        if start[self.height_ny] != stop[self.height_ny]:
            raise Exception('CPWPort: start/stop in height direction must be equal')

        self.direction = np.sign(stop[self.prop_ny] - start[self.prop_ny])
        if self.direction == 0:
            raise Exception('CPWPort: start/stop in prop direction must differ')

        feed_shift = kw.get('FeedShift', 0)
        feed_R     = kw.get('Feed_R', np.inf)

        nstart = np.minimum(start, stop)
        nstop  = np.maximum(start, stop)

        # CPW metal layer
        metal_prop.AddBox(np.array(start), np.array(stop), priority=self.priority)
        self.port_props.append(metal_prop)

        # Measurement plane
        measplane_pos = 0.5 * (nstart[self.prop_ny] + nstop[self.prop_ny])
        if 'MeasPlaneShift' in kw:
            measplane_pos = start[self.prop_ny] + self.direction * kw['MeasPlaneShift']

        mesh       = CSX.GetGrid()
        prop_lines = mesh.GetLines(self.prop_ny)
        idx = np.argmin(np.abs(prop_lines - measplane_pos))
        idx = max(1, min(idx, len(prop_lines) - 2))
        meshlines = prop_lines[idx-1:idx+2]
        if self.direction < 0:
            meshlines = meshlines[::-1]

        self.measplane_shift = abs(meshlines[1] - start[self.prop_ny])
        self.U_delta = np.diff(meshlines)
        i_pos = meshlines[:2] + np.diff(meshlines) / 2.0
        self.I_delta = np.diff(i_pos)

        # Half-width and gap offsets in width direction
        w_center = 0.5 * (nstart[self.width_ny] + nstop[self.width_ny])
        half_w   = 0.5 * (nstop[self.width_ny] - nstart[self.width_ny])

        w_add_start = np.zeros(3)
        w_add_stop  = np.zeros(3)
        w_add_start[self.width_ny] = half_w
        w_add_stop[self.width_ny]  = half_w + gap_width

        # Voltage probes: pairs (left/right gap) at three prop positions
        suffix_pairs = [('A1', 'A2'), ('B1', 'B2'), ('C1', 'C2')]
        self.U_filenames = []
        for n, (s1, s2) in enumerate(suffix_pairs):
            v_pt = np.zeros(3)
            v_pt[self.prop_ny]   = meshlines[n]
            v_pt[self.width_ny]  = w_center
            v_pt[self.height_ny] = start[self.height_ny]

            for s, sign in [(s1, -1), (s2, +1)]:
                u_name = self.lbl_temp.format('ut') + s
                self.U_filenames.append(u_name)
                u_probe = self._AddProbe(CSX, u_name, p_type=0, weight=0.5)
                u_probe.AddBox(v_pt + sign*w_add_start, v_pt + sign*w_add_stop,
                               priority=self.priority)
                self.port_props.append(u_probe)

        # Current probes: span CPW + gaps width, one cell in height
        height_lines = mesh.GetLines(self.height_ny)
        width_lines  = mesh.GetLines(self.width_ny)
        h_idx = np.argmin(np.abs(height_lines - start[self.height_ny]))
        h_idx = max(2, min(h_idx, len(height_lines) - 3))

        w_idx_start = np.argmin(np.abs(width_lines - nstart[self.width_ny]))
        w_idx_start = max(1, min(w_idx_start, len(width_lines) - 1))
        w_idx_stop  = np.argmin(np.abs(width_lines - nstop[self.width_ny]))
        w_idx_stop  = max(0, min(w_idx_stop,  len(width_lines) - 2))

        i_base_start = np.zeros(3)
        i_base_start[self.width_ny]  = 0.5 * (width_lines[w_idx_start-1] + width_lines[w_idx_start])
        i_base_start[self.height_ny] = 0.5 * (height_lines[h_idx-2] + height_lines[h_idx-1])
        i_base_stop = np.zeros(3)
        i_base_stop[self.width_ny]   = 0.5 * (width_lines[w_idx_stop]   + width_lines[w_idx_stop+1])
        i_base_stop[self.height_ny]  = 0.5 * (height_lines[h_idx+1] + height_lines[h_idx+2])

        self.I_filenames = []
        for n, s in enumerate(['A', 'B']):
            i_start = i_base_start.copy()
            i_stop  = i_base_stop.copy()
            i_start[self.prop_ny] = 0.5 * (meshlines[n]   + meshlines[n+1])
            i_stop[self.prop_ny]  = i_start[self.prop_ny]

            i_name = self.lbl_temp.format('it') + s
            self.I_filenames.append(i_name)
            i_probe = self._AddProbe(CSX, i_name, p_type=1, weight=self.direction, norm_dir=self.prop_ny)
            i_probe.AddBox(i_start, i_stop)
            self.port_props.append(i_probe)

        # Excitation: two gap boxes, E-field in the width (gap) direction
        if excite != 0:
            feed_idx = np.argmin(np.abs(prop_lines - (start[self.prop_ny] + feed_shift * self.direction)))
            ex_pt = np.zeros(3)
            ex_pt[self.prop_ny]   = prop_lines[feed_idx]
            ex_pt[self.width_ny]  = w_center
            ex_pt[self.height_ny] = start[self.height_ny]

            for lbl_s, sign in [('excite_1', -1), ('excite_2', +1)]:
                exc_val = np.zeros(3)
                exc_val[self.width_ny] = sign*excite
                exc = CSX.AddExcitation(self.lbl_temp.format(lbl_s), exc_type=0,
                                        exc_val=exc_val, delay=self.delay)
                exc.AddBox(ex_pt + sign*w_add_start, ex_pt + sign*w_add_stop,
                           priority=self.priority)
                self.port_props.append(exc)

        # Termination resistance at start — centred on w_center, one box per gap
        r_pt = np.zeros(3)
        r_pt[self.prop_ny]   = start[self.prop_ny]
        r_pt[self.width_ny]  = w_center
        r_pt[self.height_ny] = start[self.height_ny]
        if feed_R > 0 and not np.isinf(feed_R):
            lumped = CSX.AddLumpedElement(self.lbl_temp.format('resist'),
                                          ny=self.width_ny, R=2*feed_R)
            for sign in [-1, +1]:
                lumped.AddBox(r_pt + sign*w_add_start, r_pt + sign*w_add_stop, priority=self.priority)
            self.port_props.append(lumped)
        elif np.isinf(feed_R):
            pass
        elif feed_R == 0:
            for sign in [-1, +1]:
                metal_prop.AddBox(r_pt + sign*w_add_start, r_pt + sign*w_add_stop, priority=self.priority)
        else:
            raise Exception('CPWPort: Feed_R must be >= 0')

    def ReadUIData(self, sim_path, freq, signal_type='pulse'):
        all_u = UI_data(self.U_filenames, sim_path, freq, signal_type)

        # Sum paired (left+right gap) probes at each of the three positions
        uf_A = all_u.ui_f_val[0] + all_u.ui_f_val[1]
        uf_B = all_u.ui_f_val[2] + all_u.ui_f_val[3]
        uf_C = all_u.ui_f_val[4] + all_u.ui_f_val[5]
        ut_B = all_u.ui_val[2]   + all_u.ui_val[3]

        self.uf_tot = uf_B
        self.ut_tot = ut_B

        self.i_data = UI_data(self.I_filenames, sim_path, freq, signal_type)
        self.if_tot = 0.5 * (self.i_data.ui_f_val[0] + self.i_data.ui_f_val[1])
        self.it_tot = 0.5 * (self.i_data.ui_val[0]   + self.i_data.ui_val[1])

        unit = self.CSX.GetGrid().GetDeltaUnit()
        Et   = uf_B
        dEt  = (uf_C - uf_A) / (np.sum(np.abs(self.U_delta)) * unit)
        Ht   = self.if_tot
        dHt  = (self.i_data.ui_f_val[1] - self.i_data.ui_f_val[0]) / (np.abs(self.I_delta[0]) * unit)

        beta = np.sqrt(-dEt * dHt / (Ht * Et))
        beta[np.real(beta) < 0] *= -1
        self.beta  = beta
        self.Z_ref = np.sqrt(Et * dEt / (Ht * dHt))
def AddCPWPort(self, port_nr, metal_prop, start, stop, prop_dir, exc_dir, gap_width, excite=0, **kw):
        """ AddCPWPort(port_nr, metal_prop, start, stop, prop_dir, exc_dir, gap_width, excite=0, **kw)

        Add a coplanar waveguide port.

        See Also
        --------
        openEMS.ports.CPWPort
        """
        if self.__CSX is None:
            raise Exception('AddCPWPort: CSX is not set!')
        return CPWPort(self.__CSX, port_nr, metal_prop, start, stop, prop_dir, exc_dir, gap_width, excite, **kw)
