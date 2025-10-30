//-----------------------------------------------------------------------------
//  Jul 26 2001 P.Murat  : STNTUPLE utility routines
//-----------------------------------------------------------------------------
#include "Stntuple/mod/StntupleUtilities.hh"

namespace stntuple {

//-----------------------------------------------------------------------------
// returns proper time of the mothers, in ns
//-----------------------------------------------------------------------------
  double get_proper_time(const mu2e::SimParticle* Simp) {
    double dt(0.);
    const bool verbose(false);

    // The mu2ePrimary code means that G4 track was created by our
    // PrimaryGeneratorAction.  If the particle is "mu2ePrimary" but
    // still has a parent, it is a continuation of a particle from the
    // previous simulation stage, and their proper times should be
    // combined.

    const mu2e::SimParticle* part  (Simp);
    const mu2e::SimParticle* parent(nullptr);
    const auto process = mu2e::ProcessCode::mu2ePrimary; // code to check for

    while(part->parent().isNonnull()) {
      if(verbose) std::cout << __func__ << ": part PDG = " << part->pdgId() << " code = " << part->creationCode() << " endProperTime = " << part->endProperTime()
                            << " dt = " << dt << std::endl;
      if (part->creationCode() == process) {
//-----------------------------------------------------------------------------
// The current particle is a continuation from the previous stage, 
// not a physically different particle.  We need to find its record in 
// the StepPointMC collection to get the correct proper time.
//-----------------------------------------------------------------------------
        parent = part->parent().get();

        if (parent->pdgId() == part->pdgId()) {
          dt  += parent->endProperTime();
        }
        else {
          break;
        }
      }
      else {
        break;
      }
      part = parent;
    }
    if(verbose) std::cout << __func__ << ": --> final dt = " << dt << std::endl;

    return dt;
  }
}
