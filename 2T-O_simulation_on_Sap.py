import xrayutilities as xu
from xrayutilities.materials import elements as el
from matplotlib.pylab import *
from xrayutilities.materials.material import (Crystal, HexagonalElasticTensor)
from xrayutilities.materials.spacegrouplattice import SGLattice
import os
import tkinter
from tkinter import filedialog

GaN = Crystal("GaN",
              SGLattice(186, 3.189, 5.186, atoms=[el.Ga, el.N],
                        pos=[('2b', 0), ('2b', 3 / 8.)]),
              HexagonalElasticTensor(390.e9, 145.e9, 106.e9, 398.e9, 105.e9),
              thetaDebye=600)
AlN = Crystal('AlN',
              SGLattice(186, 3.1130, 4.9816, atoms=[el.Al, el.N],
                        pos=[('2b', 0), ('2b', 3 / 8.)]),
              xu.materials.HexagonalElasticTensor(410.e9, 149.e9, 99.e9, 389.e9, 125.e9),
              thetaDebye=1150)
InN = Crystal('AlN',
              SGLattice(186, 3.5380, 5.7020, atoms=[el.In, el.N],
                        pos=[('2b', 0), ('2b', 3 / 8.)]),
              xu.materials.HexagonalElasticTensor(190.e9, 104.e9, 121.e9, 182.e9, 10.e9),
              thetaDebye=660)
AlGaN = xu.materials.material.Alloy(GaN, AlN, x=0.25)

wavelength = xu.wavelength('CuKa1')
offset = 0

sub = xu.simpack.Layer(xu.materials.Al2O3, 6300)
lay1 = xu.simpack.Layer(GaN, 2300)
lay2 = xu.simpack.Layer(AlN, 1, relaxation=0)
lay3 = xu.simpack.Layer(AlGaN, 25, relaxation=0)
epi = xu.simpack.LayerStack('list',  lay1 + lay2 + lay3)

thetaMono = arcsin(wavelength / (2 * xu.materials.Ge.planeDistance(2, 2, 0)))
Cmono = cos(2 * thetaMono)
dyn = xu.simpack.DynamicalModel(epi, I0=1.5e9, background=0,
                                resolution_width=2e-3, polarization='both',
                                Cmono=Cmono)
fitmdyn = xu.simpack.FitModel(dyn)
fitmdyn.set_param_hint('GaN_thickness', vary=True)
fitmdyn.set_param_hint('GaN_a', vary=True)
fitmdyn.set_param_hint('GaN_c', vary=True)
fitmdyn.set_param_hint('AlN_thickness', vary=True)
fitmdyn.set_param_hint('AlN_a', vary=True,)
fitmdyn.set_param_hint('AlN_c', vary=True)
fitmdyn.set_param_hint('GaN_0_75_AlN_0_25__thickness', vary=True)
fitmdyn.set_param_hint('GaN_0_75_AlN_0_25__a', vary=True)
fitmdyn.set_param_hint('GaN_0_75_AlN_0_25__c', vary=True)
fitmdyn.set_param_hint('GaN_0_75_AlN_0_25__at0_Ga_2b_occupation', vary=True)
fitmdyn.set_param_hint('GaN_0_75_AlN_0_25__at2_Al_2b_occupation', vary=True)
fitmdyn.set_param_hint('resolution_width', vary=True)
params = fitmdyn.make_params()

# plot experimental data
f = figure(figsize=(7, 5))

root = tkinter.Tk()
root.withdraw()
fTyp = [("xrdファイル", "*.xrdml")]
iDir = os.path.abspath(os.path.dirname(r"\\Sirius\carbon-nas\SR4000\003_XRD"))
filepath = filedialog.askopenfilename(filetypes=fTyp, initialdir=iDir)
filename = os.path.basename(filepath)
dirname = os.path.dirname(filepath)
d = xu.io.panalytical_xml.XRDMLFile(filename,path=dirname)
scan = d.scans[-1]
tt = scan.scanmot - offset
semilogy(tt, scan.int, 'o-', ms=3, label='data')

# perform fit and plot the result
fitmdyn.lmodel.set_hkl((0, 0, 2))
ai = (d.scans[-1].scanmot - offset) / 2
fitr = fitmdyn.fit(d.scans[-1].int, params, ai)
print(fitr.fit_report())  # for older lmfit use: lmfit.report_fit(fitr)


plt.show()