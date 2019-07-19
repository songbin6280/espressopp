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

#include <string>
#include "python.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "interaction/Interaction.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"
#include "ExtPlumed.hpp"
#include "mpi.hpp"

namespace espressopp {

  using namespace iterator;
  using std::string;

  namespace integrator {

    ExtPlumed::ExtPlumed(shared_ptr<System> _system, python::object _pyobj, string _dat, string _log, real _dt, bool _restart):
      Extension(_system),
      dt(_dt),
      step(0),
      nreal(0),
      gatindex(NULL),
      masses(NULL),
      f(NULL),
      pos(NULL),
      charges(NULL),
      particlesChanged(false)
    {
      p=new PLMD::Plumed;
      int real_precision = sizeof(real);
      p->cmd("setRealPrecision",&real_precision);
      p->cmd("setMDEngine","ESPResSo++");
      MPI_Comm comm = MPI_Comm(*_system->comm);
      p->cmd("setMPIComm", &comm);
      bool dat_is_file = (_dat.find_first_of("\n") > _dat.size()); // test if input is a multline string.
      if (dat_is_file) p->cmd("setPlumedDat", _dat.c_str());
      p->cmd("setLogFile", _log.c_str());
      p->cmd("setTimestep",&dt);
      longint tmp = _system->storage->getNRealParticles();
      boost::mpi::all_reduce(*_system->comm, tmp, natoms, std::plus<longint>());
      p->cmd("setNatoms",&natoms);
      if (_restart) {
        int res = 1;
        p->cmd("setRestart", &res);
      }
      p->cmd("init");
      if (!dat_is_file) {
        std::istringstream iss(_dat);
        std::string token;
        iss >> std::ws; // remove leading white spaces in a line
        while(std::getline(iss, token)) {
          iss >> std::ws; // remove leading white spaces in a line
          if (token[0] == '#') continue;
          p->cmd("readInputLine", token.c_str());
          token.clear();
        }
      }
      _onParticlesChanged = _system->storage->onParticlesChanged.connect(boost::bind(&ExtPlumed::onParticlesChanged, this));
    }

    ExtPlumed::~ExtPlumed() {
      _onParticlesChanged.disconnect();
      delete [] f;
      delete [] pos;
      delete [] charges;
      delete [] gatindex;
      delete [] masses;
      delete p;
    }

    real ExtPlumed::getBias() {
      return bias;
    }

    void ExtPlumed::onParticlesChanged() {
      particlesChanged = true;
    }

    void ExtPlumed::disconnect() {
      _runInit.disconnect();
      _aftCalcF.disconnect();
      _aftIntP.disconnect();
    }

    void ExtPlumed::connect() {
      _runInit = integrator->runInit.connect( boost::bind(&ExtPlumed::setStep, this));
      _aftCalcF = integrator->aftCalcF.connect( boost::bind(&ExtPlumed::updateForces, this));
      _aftIntP = integrator->aftIntP.connect( boost::bind(&ExtPlumed::updateStep, this));
    }

    void ExtPlumed::setStep() {
      step = integrator->getStep();
    }

    void ExtPlumed::updateForces() {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      if (nreal!=system.storage->getNRealParticles()) {
        nreal = system.storage->getNRealParticles();
        if(charges) delete [] charges;
        if(masses) delete [] masses;
        if(gatindex) delete [] gatindex;
        if(pos) delete [] pos;
        if(f) delete [] f;
        gatindex = new int [nreal];
        masses = new real [nreal];
        charges = new real [nreal];
        pos = new real[nreal*3];
        f = new real[nreal*3];

        for(auto tp=std::make_pair(0, CellListIterator(realCells));
            tp.first<nreal && !tp.second.isDone();
            ++tp.first, ++tp.second)
          {
            gatindex[tp.first] = static_cast<int>(tp.second->id())-1;
          }
        particlesChanged = false;

      } else if (particlesChanged) {
        for(auto tp=std::make_pair(0, CellListIterator(realCells));
            tp.first<nreal && !tp.second.isDone();
            ++tp.first, ++tp.second)
          {
            gatindex[tp.first] = static_cast<int>(tp.second->id())-1;
          }
        particlesChanged = false;
      }

      for(auto tp=std::make_pair(0, CellListIterator(realCells));
          tp.first<nreal && !tp.second.isDone();
          ++tp.first, ++tp.second)
        {
          masses[tp.first] = tp.second->mass();
          charges[tp.first] = tp.second->q();
          std::copy(tp.second->force().begin(), tp.second->force().end(), &f[tp.first*3]);
          std::copy(tp.second->position().begin(), tp.second->position().end(), &pos[tp.first*3]);
        }

      p->cmd("setStep",&step);
      p->cmd("setAtomsNlocal",&nreal);
      p->cmd("setAtomsGatindex",&gatindex[0]);

      real box[3][3];
      for(int i=0;i<3;i++) for(int j=0;j<3;j++) box[i][j]=0.0;
      Real3D L = system.bc->getBoxL();
      box[0][0]=L[0];
      box[1][1]=L[1];
      box[2][2]=L[2];

      real virial[3][3];
      for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) virial[i][j]=0.0;

      p->cmd("setPositions", pos);
      p->cmd("setForces", f);
      p->cmd("setBox",&box[0][0]);
      p->cmd("setMasses",masses);
      p->cmd("setCharges",charges);
      p->cmd("setVirial", &virial[0][0]);
      p->cmd("getBias",&bias);
      p->cmd("prepareCalc");

      int plumedNeedsEnergy = 0;
      p->cmd("isEnergyNeeded", &plumedNeedsEnergy);
      if (plumedNeedsEnergy) {
        real pot_energy = 0.;
        const interaction::InteractionList& srIL = system.shortRangeInteractions;
        for (size_t j =0; j < srIL.size(); ++j) {
          pot_energy += srIL[j]->computeEnergy();
        }
        pot_energy /= system.comm->size(); // PLUMED defines PE this way.
        p->cmd("setEnergy", &pot_energy);
      }
      p->cmd("performCalc");

      for(auto tp=std::make_pair(0, CellListIterator(realCells));
          tp.first<nreal && !tp.second.isDone();
          ++tp.first, ++tp.second)
        {
          std::copy(&f[tp.first*3], &f[tp.first*3+3], tp.second->force().begin());
        }
    }

    void ExtPlumed::updateStep() {
      step = integrator->getStep()+1;
    }

    /****************************************************
     ** REGISTRATION WITH PYTHON
     ****************************************************/

    void ExtPlumed::registerPython() {

      using namespace espressopp::python;

      class_<ExtPlumed, shared_ptr<ExtPlumed>, bases<Extension> >

        ("integrator_ExtPlumed", init< shared_ptr< System >, python::object, string, string, real, bool>())
        .def("getBias", &ExtPlumed::getBias)
        .def("connect", &ExtPlumed::connect)
        .def("disconnect", &ExtPlumed::disconnect)
        ;
    }
  }
}
