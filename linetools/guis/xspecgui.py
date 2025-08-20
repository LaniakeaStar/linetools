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

        self.norm = norm

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

        self.pltline_widg.llist = ltgu.set_llist('ISM', in_dict=self.pltline_widg.llist)
        self.pltline_widg.llist['Plot'] = True

        # Refleja en el widget visual
        try:
            idx = self.pltline_widg.lists.index(self.pltline_widg.llist['List'])
            self.pltline_widg.llist_widget.setCurrentRow(idx)
        except Exception:
            pass

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

            self.spec_list.append(self._ensure_vacuum_angstrom(spec))



        # Hook the spec widget to Plot Line
        self.spec_widg = ltgsp.ExamineSpecWidget(self.spec_list[0], guessfile=guessfile, voigtsfit=voigtsfit,
                                         status=self.statusBar, parent=self, llist=self.pltline_widg.llist,
                                         zsys=zsys, norm=norm, exten=self.exten_list[0],
                                         abs_sys=abs_sys, screen_scale=self.scale,
                                         rsp_kwargs=rsp_kwargs, **kwargs)
        
        self.spec_widg.llist = self.pltline_widg.llist
        
        # Reset redshift from spec
        if zsys is None:
            if hasattr(self.spec_widg.spec, 'z'):
                self.pltline_widg.setz(str(self.spec_widg.spec.z[self.spec_widg.select]))
        # Auto set line list if spec has proper object type
        if hasattr(self.spec_widg.spec, 'stypes'):
            st = str(self.spec_widg.spec.stypes[self.spec_widg.select]).lower()
            if st == 'galaxy':
                self.pltline_widg.llist = ltgu.set_llist('Galaxy', in_dict=self.pltline_widg.llist)
            elif st == 'absorber':
                self.pltline_widg.llist = ltgu.set_llist('Strong', in_dict=self.pltline_widg.llist)

            self.pltline_widg.llist['Plot'] = True
            idx = self.pltline_widg.lists.index(self.pltline_widg.llist['List'])
            self.pltline_widg.llist_widget.setCurrentRow(idx)
            # vuelve a compartir el dict actualizado
            self.spec_widg.llist = self.pltline_widg.llist
        #
        
        self.pltline_widg.spec_widg = self.spec_widg
        # Multi spec
        self.mspec_widg = ltgsp.MultiSpecWidget(self.spec_widg)

        #self.mspec_widg.clear()
        #for f in self.ispec_list:
        #    self.mspec_widg.addItem(str(f))

        #self.mspec_widg.currentRowChanged.connect(self.change_spectrum)

        self.spec_selector = QListWidget()
        self.spec_selector.addItems([str(f) for f in self.ispec_list])
        self.spec_selector.currentRowChanged.connect(self.change_spectrum)


        self.spec_widg.canvas.mpl_connect('button_press_event', self.on_click)
        self.spec_widg.canvas.mpl_connect('draw_event', self._ensure_visible_after_draw)

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
                self._ensure_llist_ready()
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
            self._ensure_llist_ready()

            # Accede de forma segura al objeto de la lista y su unidad
            ll = self.pltline_widg.llist
            llist_name = ll.get('List', 'ISM') or 'ISM'
            # Si por alguna razón no está poblada, repoblarla
            self.pltline_widg.llist = ltgu.set_llist(llist_name, in_dict=ll)
            ll = self.pltline_widg.llist  # refresca ref local

            llist_obj = ll[ll['List']]
            wrest_col = llist_obj._data['wrest']
            wrest_unit = getattr(wrest_col, 'unit', u.AA)  # fallback a Å si no hay unidad

            # wrest del item seleccionado (en la UI viene como string, tomamos el número)
            try:
                wrest_val = float(spltw[0])
            except Exception:
                # último recurso: quitar posibles letras/unidades del número
                wrest_val = float(''.join(ch for ch in spltw[0] if (ch.isdigit() or ch in '.-+eE')))

            # Normaliza a Å
            wrest = Quantity(wrest_val, wrest_unit).to(u.AA)

            # xdata debe ser finito
            if (event.xdata is None) or (not np.isfinite(event.xdata)):
                return

            # Calcula z
            z = event.xdata / wrest.value - 1.0
            print("z={:.5f}".format(z))
            self.statusBar().showMessage('z = {:f}'.format(z))

            # Usa la API oficial (refresca overlays y widgets internos)
            self.pltline_widg.setz('{:.5f}'.format(z))

            # Asegura que ambos sigan compartiendo el MISMO dict
            self.spec_widg.llist = self.pltline_widg.llist

            # Redibuja de forma segura (reintenta si la lista quedó momentáneamente inconsistente)
            try:
                self.spec_widg.on_draw()
            except KeyError:
                # p.ej. llist['List'] fue None por un instante → recupera y reintenta
                self._ensure_llist_ready()
                self.spec_widg.llist = self.pltline_widg.llist
                self.spec_widg.on_draw()


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

        # --- AIR/VAC utils (Edlén/Peck, suficiente para GUI) ---
    def _air_to_vacuum_angstrom(self, wav_air):
        w = np.asarray(wav_air, dtype=float)  # Å en aire
        s2 = (1e4 / w)**2
        n = 1.0 + 0.0000834254 + 0.02406147/(130.0 - s2) + 0.00015998/(38.9 - s2)
        return w * n  # Å en vacío

    def _ensure_vacuum_angstrom(self, spec):
        """Ensure wavelength is in Å. Convert AIR->VAC if header says AIR or instrument is KCWI."""
        try:
            unit = getattr(spec.wavelength, 'unit', None)
            if unit is None:
                return spec

            wav = spec.wavelength.to(u.AA).value

            meta = getattr(spec, 'meta', {}) or {}
            hdr = (str(meta.get('WAVEFRAME', '')).lower() + ' ' +
                   str(meta.get('AIRORVAC', '')).lower() + ' ' +
                   str(meta.get('CTYPE1', '')).lower() + ' ' +
                   str(meta.get('WCSNAME', '')).lower())

            instr = str(meta.get('INSTRUME', '')).lower()

            is_kcwi = ('kcwi' in instr)   
            says_air = ('air' in hdr) or ('air' in instr)

            if is_kcwi or says_air:
                wav = self._air_to_vacuum_angstrom(wav)

            spec.wavelength = (wav * u.AA)
        except Exception:
            pass
        return spec
    
    def _ensure_llist_ready(self):
        """ Ensure the line list is set and ready to plot"""
        try:
            # Nombre actual o ISM por defecto
            cur = self.pltline_widg.llist.get('List', 'ISM') or 'ISM'
            # Repoblar el dict (esto garantiza la clave self.pltline_widg.llist[cur] con el objeto LineList)
            self.pltline_widg.llist = ltgu.set_llist(cur, in_dict=self.pltline_widg.llist)

            # Sincroniza selección visual
            try:
                row = self.pltline_widg.lists.index(self.pltline_widg.llist['List'])
                self.pltline_widg.llist_widget.setCurrentRow(row)
            except Exception:
                pass

            # Plot ON y comparte el MISMO dict con spec_widg
            self.pltline_widg.llist['Plot'] = True
            self.spec_widg.llist = self.pltline_widg.llist
        except Exception:
            pass



    # Quit
    def quit(self):
        self.close()
    
    def change_spectrum(self, index):
        """Switch spectrum shown in the main plot y autocalibra ejes."""
        if index < 0 or index >= len(self.spec_list):
            return

        # --- Asegura que la calibración de longitudes de onda esté en vacío ---
        new_spec = self._ensure_vacuum_angstrom(self.spec_list[index])

        # --- (2) Refresca la line list, asegurando que no quede 'None' ---
        self._ensure_llist_ready()

        # --- (3) Comparte el MISMO dict con spec_widg (importantísimo) ---
        self.spec_widg.llist = self.pltline_widg.llist

        # Cambia el espectro
        self.spec_widg.set_spectrum(new_spec)

        # Sincroniza z del espectro con la line list
        self._sync_line_list_to_spec(new_spec)

        # Limpia el eje
        ax = getattr(self.spec_widg.canvas, 'ax', None)
        if ax is None:
            ax = self.spec_widg.canvas.figure.gca()
        ax.clear()

        # Redibuja
        self.spec_widg.on_draw()
        self._autoscale_from_spec(new_spec)



    def _sync_line_list_to_spec(self, spec):
        """Sincroniza z del LineList con el del espectro y refresca líneas correctamente."""
        # 1) Obtiene z del XSpectrum1D (si existe)
        z = None
        try:
            z_attr = getattr(spec, 'z', None)
            if z_attr is not None:
                if np.ndim(z_attr) == 0:
                    z = float(z_attr)
                else:
                    sel = getattr(self.spec_widg, 'select', 0)
                    z = float(z_attr[sel])
        except Exception:
            z = None

        if z is None or not np.isfinite(z):
            z = 0.0

        # 2) Usa la API del widget (dispara refresco de overlays)
        try:
            self.pltline_widg.setz(f"{z:.5f}")
        except Exception:
            # fallback
            self.pltline_widg.llist['z'] = z
            try:
                self.pltline_widg.zbox.setText(f"{z:.5f}")
            except Exception:
                pass

        # 3) Re‑compartir el MISMO dict por si algún método lo reemplazó
        self.spec_widg.llist = self.pltline_widg.llist

    def _ensure_visible_after_draw(self, event=None):
        try:
            ax = getattr(self.spec_widg.canvas, 'ax', None)
            if ax is None:
                ax = self.spec_widg.canvas.figure.gca()
            x1, x2 = ax.get_xlim()
            wv = np.asarray(self.spec_widg.spec.wavelength.to(u.AA).value, dtype=float)
            fx = np.asarray(self.spec_widg.spec.flux.value if hasattr(self.spec_widg.spec.flux,'value') else self.spec_widg.spec.flux, dtype=float)
            m = np.isfinite(wv) & np.isfinite(fx) & (wv >= min(x1,x2)) & (wv <= max(x1,x2))
            if m.sum() < 10:  # if less than 10 points visible, autoscale
                self._autoscale_from_spec(self.spec_widg.spec)
        except Exception:
            pass


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

