#!/usr/bin/python
#
## Xiaofeng Qian @ mit, 2010
#
import numpy, sys

#-----------------------------------------------------------#
## User defined parameters

e_strength      = 0.01           # 1/Angstrom
dt              = 2              # atto-second
E_begin         = 0              # eV
E_end           = 20             # eV
E_begin_plot    = 0              # eV
E_end_plot      = 20             # eV
NEgrid          = 3000           # number of energy grids
damping_factor  = 0.02           # damping factor
molecule        = 'sih4'         # molecule name

s = 'xyz'
if len(sys.argv) > 1:
    pol = sys.argv[1][0].lower()
    e_direction = s.find(pol) + 1
    dipole_filename = 'dip_'+pol+'.dat'
    oas_figure_filename = molecule+'_'+pol+'.png'
    oas_output_filename = molecule+'_'+pol+'.dat'
else:
    print 'usage: ./plot_optical_absorption.py  x | y | z'
    sys.exit(1)


## Do not change the lines below
#-----------------------------------------------------------#
## constants in SI unit
EV_IN_J = 1.60217733e-19
HBAR    = 1.054572669125102e-34
FEMTOSECOND = 1.0e-18
ELECTRON_REST_MASS_IN_KG = 9.1093897e-31
ELEMENTARY_CHARGE_IN_C   = 1.60217733e-19
A_IN_M = 1.0e-10

#-----------------------------------------------------------#
## read dipole file
dp_file = file(dipole_filename, "r")
data=[]
for line in dp_file.readlines():
    if not line.startswith('DIP'): continue
    data.append(map(float, line.split()[3:]))

dp_file.close()
dp      = numpy.array(data)
dp_raw  = dp[:, e_direction-1]
dp_diff = dp_raw - dp_raw[0]

## setup time and energy grid
nstep   = dp.shape[0]
dt      = dt * FEMTOSECOND * EV_IN_J / HBAR
nt      = numpy.linspace(0, (nstep-1), nstep) 
t       = nt * dt
dE      = (E_end - E_begin) / numpy.float(NEgrid-1)
dp_damping = numpy.exp( - t**2 * damping_factor )


## initialize energy axis and dipole strength function
S       = numpy.zeros(NEgrid)
E_axis  = numpy.zeros(NEgrid)

## calculate dipole strength function in energy
for it in range(NEgrid):
    Eit     = E_begin + dE * it
    E_axis[it] = Eit
    dp_fft = numpy.sum( dp_diff[range(nstep)] * dp_damping * numpy.sin( Eit * t))
    S[it]  = dp_fft* ( dE*it ) * 2 * dt /e_strength /numpy.pi

S  = S * ELECTRON_REST_MASS_IN_KG * ELEMENTARY_CHARGE_IN_C * HBAR**(-2) * A_IN_M**2

## check sum_rule = total number of valence electrons
SumRule  = numpy.sum( S ) * dE
print "Sum rule = ", SumRule

## print spectrum in output
f = open(oas_output_filename, "wt")
for i in xrange(len(E_axis)):
    f.write("%f %f\n" % (E_axis[i],S[i]))
f.close()


## plot out optical absorption spectrum
try:
    import pylab
except:
    print "pylab not available: not plotting now"
    sys.exit(1)

params = {'axes.labelsize':  16,
          'text.fontsize':   16,
          'legend.fontsize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex':     True}
pylab.rcParams.update(params)
pylab.rc("axes", linewidth=2.0)
pylab.rc("lines", linewidth=3.0)
pylab.plot(E_axis, S)
pylab.xlabel('energy (eV)')
pylab.ylabel('dipole strength function ($eV^{-1}$)')
pylab.title('optical absorption spectrum')
pylab.grid(True)
pylab.axis([E_begin_plot, E_end_plot, 0.0, numpy.max(S)*1.2])
ax = pylab.gca()
ax.title.set_fontsize(24)
pylab.savefig(oas_figure_filename)

pylab.show()


