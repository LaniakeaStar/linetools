""" Module for XSpecGui
"""
from __future__ import print_function, absolute_import, division, unicode_literals
import sys
import numpy as np
import pdb
from qtpy import QtGui
from qtpy.QtWidgets import QWidget, QLabel, QPushButton, QMainWindow
from qtpy.QtWidgets import QVBoxLayout, QHBoxLayout
from qtpy import QtCore
from qtpy.QtWidgets import QListWidget
from matplotlib import rcParams
from astropy.units import Quantity
from astropy import units as u
from linetools.guis import utils as ltgu
from linetools.guis import line_widgets as ltgl
from linetools.guis import spec_widgets as ltgsp
from linetools import utils as ltu
from linetools.analysis import voigt as lav
from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm import utils as ltiu
from linetools.lists.linelist import LineList
from linetools.spectra.xspectrum1d import XSpectrum1D
class XSpecGui(QMainWindow):
    """ GUI to replace XIDL x_specplot (which simulated a GUI by T. Barlow)
    """
    def __init__(self, ispec, guessfile=None, parent=None, zsys=None, norm=None, exten=None,
                 rsp_kwargs={}, unit_test=False, screen_scale=1.,
                 **kwargs):
        QMainWindow.__init__(self, parent)
        """
        ispec = str, XSpectrum1D or tuple of arrays
          Input spectrum or spectrum filename.  If tuple then (wave,
          fx), (wave, fx, sig) or (wave, fx, sig, co)
        guessfile : str, optional
          name of the .json file generated with igmguesses GUI in Pyigm (see https://github.com/pyigm/pyigm/blob/master/docs/igmguesses.rst)
          if not None - overplot fitted line profiles from igmguesses
        parent : Widget parent, optional
        zsys : float, optional
          intial redshift
        exten : int, optional
          extension for the spectrum in multi-extension FITS file
        norm : bool, optional
          True if the spectrum is normalized
        screen_scale : float, optional
          Scale the default sizes for the gui size
        """
        #reload(ltgl)
        #reload(ltgsp)
        # INIT
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        if isinstance(ispec, list):
            self.multispec_mode = True
            self.ispec_list = ispec
            self.exten_list = exten if isinstance(exten, list) else [0] * len(ispec)
        else:
            self.multispec_mode = False
            self.ispec_list = [ispec]
            self.exten_list = [exten if exten is not None else 0]
        self.scale = screen_scale
        # Needed to avoid crash in large spectral files
        rcParams['agg.path.chunksize'] = 20000
        rcParams['axes.formatter.useoffset'] = False  # avoid scientific notation in axes tick labels
        # Build a widget combining several others
        self.main_widget = QWidget()
        # Status bar
        self.create_status_bar()
        # Grab the pieces and tie together
        self.pltline_widg = ltgl.PlotLinesWidget(status=self.statusBar,
            init_z=zsys, screen_scale=self.scale)
        self.pltline_widg.setMaximumWidth(int(300*self.scale))
        ## Abs sys
        abs_sys = None
        voigtsfit = None
        if guessfile is not None:
            # Load
            ism = LineList('ISM')
            igm_guess = ltu.loadjson(guessfile)
            comps = []
            for key in igm_guess['cmps'].keys():
                comp = AbsComponent.from_dict(igm_guess['cmps'][key], chk_vel=False, linelist=ism)
                comps.append(comp)
            abs_sys = ltiu.build_systems_from_components(comps,
                                                         vsys=500. * u.km / u.s)  # ,chk_z=False)  ### 100000.*u.km/u.s   ok
            ### voigt fit - added
            # Spectrum
            spec, spec_fil = ltgu.read_spec(ispec, exten=exten, norm=norm,
                                            rsp_kwargs=rsp_kwargs)
            voigtsfit = np.asarray([0] * len(spec.wavelength))
            alllines = []
            for iabs_sys in abs_sys:
                lines = iabs_sys.list_of_abslines()
                alllines = alllines + lines
            if len(alllines) > 0:
                voigtsfit = lav.voigt_from_abslines(spec.wavelength, alllines, fwhm=3.).flux.value
            if not norm:
                voigtsfit = voigtsfit * spec.co
        self.spec_list = []
        for i, f in enumerate(self.ispec_list):
            if f.endswith('.csv'):
                data = np.genfromtxt(f, delimiter=',', names=True)
                names = data.dtype.names
                wave = data[names[0]]
                flux = data[names[1]]
                # Aquí definimos una unidad válida de flujo
                flux_unit = u.Unit("erg / (s cm2 Angstrom)")
                if len(names) >= 3:
                    sig = data[names[2]]
                    spec = XSpectrum1D(wave * u.AA, flux * flux_unit, sig * flux_unit)
                else:
                    spec = XSpectrum1D(wave * u.AA, flux * flux_unit)
            else:
                spec, _ = ltgu.read_spec(f, exten=self.exten_list[i], norm=norm, rsp_kwargs=rsp_kwargs)
            self.spec_list.append(spec)
        # Hook the spec widget to Plot Line
        self.spec_widg = ltgsp.ExamineSpecWidget(self.spec_list[0], guessfile=guessfile, voigtsfit=voigtsfit,
                                         status=self.statusBar, parent=self, llist=self.pltline_widg.llist,
                                         zsys=zsys, norm=norm, exten=self.exten_list[0],
                                         abs_sys=abs_sys, screen_scale=self.scale,
                                         rsp_kwargs=rsp_kwargs, **kwargs)
        # Reset redshift from spec
        if zsys is None:
            if hasattr(self.spec_widg.spec, 'z'):
                self.pltline_widg.setz(str(self.spec_widg.spec.z[self.spec_widg.select]))
        # Auto set line list if spec has proper object type
        if hasattr(self.spec_widg.spec, 'stypes'):
            if self.spec_widg.spec.stypes[self.spec_widg.select].lower() == 'galaxy':
                self.pltline_widg.llist = ltgu.set_llist('Galaxy',in_dict=self.pltline_widg.llist)
            elif self.spec_widg.spec.stypes[self.spec_widg.select].lower() == 'absorber':
                self.pltline_widg.llist = ltgu.set_llist('Strong',in_dict=self.pltline_widg.llist)
            self.pltline_widg.llist['Plot'] = True
            idx = self.pltline_widg.lists.index(self.pltline_widg.llist['List'])
            self.pltline_widg.llist_widget.setCurrentRow(idx)
        #
        self.pltline_widg.spec_widg = self.spec_widg

        self.pltline_widg.llist_widget.currentRowChanged.connect(self._on_llist_changed)

        # Multi spec
        self.mspec_widg = ltgsp.MultiSpecWidget(self.spec_widg)
        #self.mspec_widg.clear()
        #for f in self.ispec_list:
        #    self.mspec_widg.addItem(str(f))
        #self.mspec_widg.currentRowChanged.connect(self.change_spectrum)
        self.current_spec_index = 0  # keep track of current spectrum index
        self._per_spec_z = {}
        self.spec_selector = QListWidget()
        self.spec_selector.addItems([str(f) for f in self.ispec_list])
        self.spec_selector.currentRowChanged.connect(self.change_spectrum)
        self.spec_widg.canvas.mpl_connect('button_press_event', self.on_click)
        # Layout
        # Extras
        extras = QWidget()
        extras.setMinimumWidth(int(180*self.scale))
        extras.setMaximumWidth(int(280*self.scale))
        vbox = QVBoxLayout()
        qbtn = QPushButton(self)
        qbtn.setText('Quit')
        qbtn.clicked.connect(self.quit)
        #vbox.addWidget(self.pltline_widg)
        #vbox.addWidget(self.mspec_widg)
        #vbox.addWidget(qbtn)
        vbox.addWidget(self.pltline_widg)
        vbox.addWidget(self.spec_selector)  # new widget for spectrum selection
        vbox.addWidget(qbtn)
        extras.setLayout(vbox)
        # Main window
        hbox = QHBoxLayout()
        hbox.addWidget(self.spec_widg)
        hbox.addWidget(extras)
        self.main_widget.setLayout(hbox)
        # Point MainWindow
        self.setCentralWidget(self.main_widget)
        if unit_test:
            self.quit()
    def create_status_bar(self):
        """ Status bar for the GUI
        """
        self.status_text = QLabel("XSpec")
        self.statusBar().addWidget(self.status_text, 1)
    def on_click(self, event):
        """ Over-loads click events
        """
        if event.button == 3: # Set redshift
            if event.xdata is None:  # Mac bug [I think]
                return
            if self.pltline_widg.llist['List'] is None:
                return
            self.select_line_widg = ltgl.SelectLineWidget(
                self.pltline_widg.llist[self.pltline_widg.llist['List']]._data,
                scale=self.scale)
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            if line.strip() == 'None':
                return
            #
            quant = line.split('::')[1].lstrip()
            spltw = quant.split(' ')
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            wrest = Quantity(float(spltw[0]), unit=self.pltline_widg.llist[
                self.pltline_widg.llist['List']]._data['wrest'].unit) # spltw[1])  [A bit risky!]
            z = event.xdata/wrest.value - 1.
            self.pltline_widg.llist['z'] = z
            self._per_spec_z[self._current_spec_index] = z     
            self.pltline_widg.zbox.setText('{:.5f}'.format(z))  

            # redraw
            self.spec_widg.on_draw()

# New
    def _finite_vals(self, q):
        """Returns finite values of a Quantity or array."""
        vals = np.asarray(q.value if hasattr(q, 'value') else q, dtype=float)
        return vals[np.isfinite(vals)]

    def _autoscale_from_spec(self, spec, pad_x=0.02, pad_y=0.05, pct_low=1.0, pct_high=99.0, is_norm=None):
        """Adjust the plot limits based on the spectrum"""

        ax = getattr(self.spec_widg.canvas, 'ax', None)
        if ax is None:
            ax = self.spec_widg.canvas.figure.gca()

        wv = self._finite_vals(spec.wavelength)
        fx = self._finite_vals(spec.flux)

        if wv.size == 0 or fx.size == 0:
            return  # nothing to scale

        # 3) X-lims per percentiles + padding
        x1, x2 = np.nanpercentile(wv, [pct_low, pct_high])
        xr = (x2 - x1)
        if xr <= 0:
            xr = max(1.0, abs(x1))  # fallback
        ax.set_xlim(x1 - pad_x * xr, x2 + pad_x * xr)

        # 4) Y-lims: if normalized centered ~1; if not, percentiles
        if is_norm is None:
            is_norm = getattr(self, 'norm', None)
            if is_norm is None:
                med = np.nanmedian(fx)
                is_norm = (0.5 < med < 2.0)

        if is_norm:
            med = np.nanmedian(fx)
            std = np.nanstd(fx)
            lo = med - 5*std
            hi = med + 5*std
            if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
                lo, hi = np.nanpercentile(fx, [pct_low, pct_high])
        else:
            lo, hi = np.nanpercentile(fx, [pct_low, pct_high])

        yr = hi - lo
        if not np.isfinite(yr) or yr <= 0:
            lo = np.nanmin(fx)
            hi = np.nanmax(fx)
            yr = hi - lo if np.isfinite(hi - lo) and (hi - lo) > 0 else 1.0

        ax.set_ylim(lo - pad_y * yr, hi + pad_y * yr)

        # Redraw
        self.spec_widg.canvas.draw_idle()  



    # Quit
    def quit(self):
        self.close()

    def _sync_llist(self, keep_plot=True):
        """
        Ensure both widgets share the same llist dict
        """
        # nombre de lista actual y z
        sel = self.pltline_widg.llist.get('List', None)
        zval = self.pltline_widg.llist.get('z', 0.0)

        if sel is not None:
            # IMPORTANTÍSIMO: no reemplazar el dict; refrescar su contenido in-place
            # usando set_llist pero sobre el MISMO in_dict
            ltgu.set_llist(sel, in_dict=self.pltline_widg.llist)  # actualiza claves del MISMO objeto
            self.pltline_widg.llist['z'] = zval
            if keep_plot:
                self.pltline_widg.llist['Plot'] = True

            # aseguramos que ambos widgets compartan EXACTAMENTE el mismo objeto
            self.spec_widg.llist = self.pltline_widg.llist
            self.pltline_widg.spec_widg = self.spec_widg

            # reflejar selección visual en el listwidget (por si acaso)
            try:
                idx = self.pltline_widg.lists.index(sel)
                self.pltline_widg.llist_widget.setCurrentRow(idx)
            except Exception:
                pass

    def _ensure_valid_llist(self):
        """
        Asegura que llist tenga un nombre de lista válido y sincroniza referencias
        entre PlotLinesWidget y ExamineSpecWidget para evitar KeyError.
        """
        ll = self.pltline_widg.llist  # MISMA referencia
        # nombre seleccionado en el panel
        sel_name = None
        try:
            row = self.pltline_widg.llist_widget.currentRow()
            if row is not None and 0 <= row < len(self.pltline_widg.lists):
                sel_name = self.pltline_widg.lists[row]
        except Exception:
            pass

        curr = ll.get('List', None)

        # Si el usuario eligió "None", no pintamos líneas
        if sel_name == 'None':
            ll['Plot'] = False
            # Evita KeyError dentro de on_draw fijando un nombre “dummy”
            if curr is None:
                ltgu.set_llist('Galaxy', in_dict=ll)
            # Sincroniza referencias
            self.spec_widg.llist = ll
            self.pltline_widg.spec_widg = self.spec_widg
            return

        # Si vamos a pintar, aseguremos nombre válido
        if ll.get('Plot', True):  # por defecto, pinta
            name = sel_name or curr or 'Galaxy'
            ltgu.set_llist(name, in_dict=ll)  # **in-place**, sin crear dict nuevo
        else:
            if curr is None:
                ltgu.set_llist('Galaxy', in_dict=ll)

        # Sincroniza referencias entre widgets
        self.spec_widg.llist = ll
        self.pltline_widg.spec_widg = self.spec_widg

        

    def _fix_limits_current_spec(self):
        """Fija xlim/ylim al espectro actualmente cargado y sincroniza ventana interna."""
        ax = getattr(self.spec_widg.canvas, 'ax', None)
        if ax is None:
            ax = self.spec_widg.canvas.figure.gca()

        spec = self.spec_widg.spec
        try:
            w = np.asarray(spec.wavelength.value, dtype=float)
            f = np.asarray(spec.flux.value, dtype=float)
        except Exception:
            return

        w = w[np.isfinite(w)]
        f = f[np.isfinite(f)]
        if w.size == 0 or f.size == 0:
            return

        x1, x2 = np.nanpercentile(w, [1, 99])
        y1, y2 = np.nanpercentile(f, [1, 99])
        xr = max(x2 - x1, 1.0)
        yr = max(y2 - y1, 1.0)
        ax.set_xlim(x1 - 0.02 * xr, x2 + 0.02 * xr)
        ax.set_ylim(y1 - 0.05 * yr, y2 + 0.05 * yr)

        new_xlim = ax.get_xlim()
        # sincroniza posibles nombres de “ventana” usados por tu versión
        for attr in ('xlim', 'wvlim', 'wvmnx', 'wlim', 'x_mnx'):
            try:
                setattr(self.spec_widg, attr, new_xlim)
            except Exception:
                pass

        self.spec_widg.canvas.draw_idle()

    def _on_llist_changed(self, *_):
        """
        Cambiar la line list desde el panel, evitando 'List'==None cuando redibujamos.
        """
        ll = self.pltline_widg.llist
        try:
            row = self.pltline_widg.llist_widget.currentRow()
            name = self.pltline_widg.lists[row] if (row is not None and row >= 0) else 'Galaxy'
        except Exception:
            name = 'Galaxy'

        if name == 'None':
            ll['Plot'] = False
            if ll.get('List') is None:
                ltgu.set_llist('Galaxy', in_dict=ll)
        else:
            ltgu.set_llist(name, in_dict=ll)
            ll['Plot'] = True

        self.spec_widg.llist = ll
        self.pltline_widg.spec_widg = self.spec_widg

        # Redibuja y reimpone límites del espectro ACTUAL
        self.spec_widg.on_draw()
        self._fix_limits_current_spec()

    def change_spectrum(self, index):
        """Switch spectrum shown in the main plot"""
        if index < 0 or index >= len(self.spec_list):
            return

        self._current_spec_index = index
        new_spec = self.spec_list[index]

        # Limpiar eje para evitar artistas “fantasma”
        ax = getattr(self.spec_widg.canvas, 'ax', None)
        if ax is None:
            ax = self.spec_widg.canvas.figure.gca()
        ax.cla()

        # Cambiar espectro
        self.spec_widg.set_spectrum(new_spec)

        # Fijar límites/ventana al NUEVO espectro (antes de dibujar)
        self._fix_limits_current_spec()

        # Asegurar que la line list activa sea válida y compartida
        # (no crea dict nuevo; actualiza in-place y comparte ref)
        ll = self.pltline_widg.llist
        try:
            row = self.pltline_widg.llist_widget.currentRow()
            name = self.pltline_widg.lists[row] if (row is not None and row >= 0) else 'Galaxy'
        except Exception:
            name = 'Galaxy'

        if name == 'None':
            ll['Plot'] = False
            if ll.get('List') is None:
                ltgu.set_llist('Galaxy', in_dict=ll)  # evita KeyError aunque no se pinte
        else:
            ltgu.set_llist(name, in_dict=ll)
            ll['Plot'] = True

        # z del NUEVO espectro: usa el recordado para este índice, si existe;
        # si el XSpectrum trae .z lo respetamos como fallback; si no, 0.0
        zval = self._per_spec_z.get(index, None)
        if zval is None:
            zval = 0.0
            if hasattr(new_spec, 'z'):
                try:
                    zval = float(new_spec.z[getattr(self.spec_widg, 'select', 0)])
                except Exception:
                    try:
                        zval = float(new_spec.z)
                    except Exception:
                        zval = 0.0
        ll['z'] = float(zval)
        self.pltline_widg.zbox.setText('{:.6f}'.format(float(zval)))  # NO usamos setz()

        # Compartir la MISMA referencia de llist entre widgets
        self.spec_widg.llist = ll
        self.pltline_widg.spec_widg = self.spec_widg

        # Redibujo final (espectro + líneas si Plot=True)
        self.spec_widg.on_draw()
        self.spec_widg.canvas.draw_idle()


def main(args, **kwargs):
    from qtpy.QtWidgets import QApplication
    from linetools.spectra.xspectrum1d import XSpectrum1D
    if not isinstance(args,(XSpectrum1D,tuple,str)):
        raise IOError("Bad input")
    # Run
    app = QApplication(sys.argv)
    gui = XSpecGui(args, **kwargs)
    gui.show()
    app.exec_()
    return