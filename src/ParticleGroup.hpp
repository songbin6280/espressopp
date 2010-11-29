// ESPP_CLASS
#ifndef __ESPRESSOPP_PARTICLEGROUP_H
#define	__ESPRESSOPP_PARTICLEGROUP_H

#include "Particle.hpp"
#include "log4espp.hpp"
#include "types.hpp"
#include <list>
#include <boost/signals2.hpp>

namespace espresso {

  /**
   * \brief group of particles
   *
   * This part contains a list of particles to e.g. organize the system into
   * molecules. The particles are stored as ids, in addition a list of active
   * particles on the processor (or will it be cell) is generated by connecting
   * to the communicator signals.
   *
   * Tis is a first try. Further extensions might be to put groups into groups.
   * 
   */
  class ParticleGroup
  {
  public:
    ParticleGroup(shared_ptr< storage::Storage > _storage);
    ~ParticleGroup();

    void add(longint pid);

    // for debugging purpose
    void print();

    static void registerPython();
    
  protected:
    std::map<longint, Particle*> active;
    // replace this by a better thing
    std::map<longint, longint> particles;
    shared_ptr< storage::Storage > storage;

    // some signalling stuff to keep track of the particles in cell
    boost::signals2::connection con_send, con_recv, con_changed;
    void beforeSendParticles(ParticleList& pl,
			     class OutBuffer& buf);
    void afterRecvParticles(ParticleList& pl,
			    class InBuffer& buf);
    void onParticlesChanged();
    
    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif	/* PARTICLEGROUP_H */

