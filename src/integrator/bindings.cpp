#include "bindings.hpp"

#include "VelocityVerlet.hpp"

namespace espresso {
  namespace integrator {
    void registerPython() {
      VelocityVerlet::registerPython();
    }
  }
}