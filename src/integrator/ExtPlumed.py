#  Copyright (C) 2017,2018
#      Max Planck Institute for Polymer Research
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""
******************************
espressopp.integrator.ExtPlumed
******************************

To call PLUMED (PLUgin for MolEcular Dynamics), set environment variable PLUMED_KERNEL to the full path of
libplumedKernel.so, and add the path of the directory of libplumed.so to environment variable LD_LIBRARY_PATH, when
PLUMED is linked through runtime linking. If PLUMED is linked by dynamic linking, i.e., it is linked as a shared
library, PLUMED_KERNEL needs not to be set. If PLUMED is linked statically, neither environment variable needs to be
changed. By default, PLUMED is linked by runtime linking. Error messages related to "undefined symbol" can be ignored,
as long as the simulation runs.

usage:

plumed = espressopp.integrator.ExtPlumed(system, "plumed.dat", "log.plumed", 0.005)
plumed.setNaturalUnits()
plumed.Init()
integrator.addExtension(plumed)

or:

plumed = espressopp.integrator.ExtPlumed(system, "plumed.dat", "log.plumed", 0.005)
plumed.setNaturalUnits()
plumed.setRestart(1)
plumed.Init()
integrator.addExtension(plumed)

.. function:: espressopp.integrator.ExtPlumed(system, cmd, log, dt)

		:param system: The Espresso++ system object.
                :type system: espressopp.System
                :param cmd: input file for PLUMED
                :type cmd: ``str``
                :param log: log file for PLUMED
                :type log:  ``str``
                :param dt:  time step
                :type dt: ``float`` (default: 0.005)

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ExtPlumed
import mpi4py.MPI as MPI

class ExtPlumedLocal(ExtensionLocal, integrator_ExtPlumed):

    def __init__(self, system, cmd, log, dt):
        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, integrator_ExtPlumed, system, pmi._PMIComm.getMPIsubcomm(), cmd, log, dt)
            else:
                pass
        else:
            cxxinit(self, integrator_ExtPlumed, system, pmi._MPIcomm, cmd, log, dt)

if pmi.isController :
    class ExtPlumed(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.ExtPlumedLocal',
            pmicall = ['getBias', 'setNaturalUnits', 'setTimeUnit', 'setEnergyUnit', 'setLengthUnit', 'setKbT', 'setRealPrecison', 'setMDChargeUnit', 'setMDMassUnit', 'setRestart', 'readInputLine', 'Init']
            )
