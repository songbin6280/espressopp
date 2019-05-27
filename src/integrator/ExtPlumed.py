#  Copyright (C) 2018
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
*************************************************
espressopp.integrator.ExtPlumed
*************************************************

This extension serves as interface between ESPResSo++ and `PLUMED <http://www.plumed.org/home>`_. PLUMED can
be used to run realtime analysis and bias the simulation along chosen collective variables. Details of using
PLUMED and its input file can be found on its website.

To call PLUMED, set environment variable :envvar:`PLUMED_KERNEL` to the full path of libplumedKernel.so, and add
the path of the directory of libplumed.so to environment variable :envvar:`LD_LIBRARY_PATH`, when PLUMED is linked through runtime linking.
If PLUMED is linked by dynamic linking, i.e., it is linked as a shared library to ESPResSo++, :envvar:`PLUMED_KERNEL` needs not to be set.
If PLUMED is linked statically, neither environment variable needs to be provided.
Error messages regarding to "undefined symbol" can be ignored, as long as the simulation runs.

By default, PLUMED is linked by runtime linking. The user can change the way of PLUMED linking to ESPResSo++ by changing CMake variable :makevar:`PLUMED_LINK_TYPE`
to either "shared" or "static" from the default "runtime". If your distribution of PLUMED is not installed in a default system directory, you can define the root
path of PLUMED installation as CMake variable :makevar:`PLUMED_HOME`. If you keep the default "runtime" linking, your installation of PLUMED is not searched. Hence, you
do not need to worry about whether your PLUMED is not installed in a default system directory or setting :makevar:`PLUMED_HOME`.

usage:

.. code:: python

    plumed = espressopp.integrator.ExtPlumed(system, "plumed.dat", "log.plumed", 0.005)
    plumed.setNaturalUnits()
    plumed.Init()
    integrator.addExtension(plumed)

or:

.. code:: python

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

    def __init__(self, system, cmd, log, dt=0.005):
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
            pmicall = ['getBias', 'setNaturalUnits', 'setTimeUnit', 'setEnergyUnit', 'setLengthUnit', 'setKbT', 'setRealPrecison', 'setMDChargeUnit', 'setMDMassUnit', 'setRestart', 'Init']
            )
