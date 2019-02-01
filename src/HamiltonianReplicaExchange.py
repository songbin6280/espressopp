#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
*************************************
espressopp.HamiltonianReplicaExchange
*************************************

"""

from espressopp import MultiSystem
from espressopp import pmi
import random
from math import exp

class HamiltonianReplicaExchange(object):

    def __init__(self, NumOfSys = 4, RNG= None):
        pass

    def startDefiningSystem(self, n):
        if not (n in range(0,self._nsystems)):
            print "ERROR: System number must be between 0 and ",self._nsystems
        else:
            pmi.activate(self._comm[n])
            self._multisystem.beginSystemDefinition()

    def endDefiningSystem(self, n):
        if not (n in range(0,self._nsystems)):
            print "ERROR: System number must be between 0 and ",self._nsystems
        else:
            pmi.deactivate(self._comm[n])

    def getNumberOfSystems(self):
        return self._nsystems

    def getNumberOfCPUsPerSystem(self):
        return self._ncpuspersystem

    def setIntegrator(self, integrator, thermostat):
        self._multisystem.setIntegrator(integrator)
        self._thermostat.append(thermostat)

    def setAnalysisE(self, analysisE):
        self._multisystem.setAnalysisPotential(analysisE)

    def setAnalysisT(self, analysisT):
        self._multisystem.setAnalysisTemperature(analysisT)

    def setAnalysisNPart(self, analysisNPart):
        self._multisystem.setAnalysisNPart(analysisNPart)

    def setDumpConfXYZ(self, dumpconf):
        set._multisystem.setDumpConfXYZ(dumpconf)

    def runDumpConfXYZ(self):
        set._multisystem.runDumpConfXYZ()

    def run(self, nsteps):
        self._multisystem.runIntegrator(nsteps)

    def exchange(self):
        pass
