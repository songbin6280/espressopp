/*
  Copyright (C) 2018
  Max Planck Institute for Polymer Research

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// ESPP_CLASS
#ifndef _INTEGRATOR_ExtPlumed_HPP
#define _INTEGRATOR_ExtPlumed_HPP

#include "python.hpp"
#include "types.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "mpi.hpp"
#include "Particle.hpp"
#include "Plumed.h"

namespace espressopp {
  namespace integrator {

    /** ExtPlumed */
    class ExtPlumed : public Extension {

    public:
      ExtPlumed(shared_ptr < System >, python::object, std::string, std::string, real);
      real getBias();
      void setNaturalUnits();
      void setRealPrecision(int);
      void setMDChargeUnits(real);
      void setMDMassUnits(real);
      void setRestart(int);
      void setTimeUnit(real);
      void setEnergyUnit(real);
      void setLengthUnit(real);
      void setKbT(real);
      void Init();

      virtual ~ExtPlumed();
      /** Register this class so it can be used from Python. */
      static void registerPython();

    private:
      PLMD::Plumed * p;
      std::string dat;
      std::string units;
      std::string log;
      real dt;
      int step;

      longint nreal; // total number of atoms (real & ghost) on the processor
      longint natoms; // total number of atoms
      int *  gatindex;
      real * masses;
      real * charges;
      real * pos;
      real * f;
      real bias;
      int plumedNeedsEnergy;
      bool dat_is_file;

      boost::signals2::connection _runInit, _aftCalcF, _aftIntP;
      void connect();
      void disconnect();
      void setStep();
      void updateForces();
      void updateStep();
    };
  }
}

#endif
