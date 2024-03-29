///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/gui/TMu2eBField.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/BFieldManagerMaker.hh"
#include "Offline/GeometryService/inc/BFieldConfigMaker.hh"
#include "Offline/GeometryService/inc/BeamlineMaker.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "CLHEP/Vector/ThreeVector.h"

#include "TBranch.h"
#include "TNtuple.h"

ClassImp(stntuple::TMu2eBField)

namespace stntuple {
//-----------------------------------------------------------------------------
TMu2eBField::TMu2eBField(): TEveMagField() {

  bool _allowReplacement    (true);
  bool _messageOnReplacement(false);
  bool _messageOnDefault    (false);

  std::string input_file = "Mu2eG4/test/geom_01.txt";

  mu2e::SimpleConfig* _config = new mu2e::SimpleConfig(input_file,
						       _allowReplacement,
						       _messageOnReplacement,
						       _messageOnDefault );

  std::unique_ptr<mu2e::Beamline> tmp = mu2e::BeamlineMaker::make(*_config);
  fBeamline = std::move(tmp).get();
  tmp.release();

  std::unique_ptr<mu2e::BFieldConfig> bfc(mu2e::BFieldConfigMaker(*_config, *fBeamline).getBFieldConfig());
  fBfc      = bfc.get();
  bfc.release();

  fBfmm     = new mu2e::BFieldManagerMaker(*fBfc);
  std::unique_ptr<mu2e::BFieldManager> bfm = fBfmm->getBFieldManager();
  fBmgr     = bfm.get();
  bfm.release();
}

//-----------------------------------------------------------------------------
TMu2eBField::~TMu2eBField() {
  delete fBeamline;
  delete fBfc;
  delete fBfmm;
  delete fBmgr;
}

//-----------------------------------------------------------------------------
// for some reason, TEveTrackPropagator RungeKutta method explicitly inverts
// the magnetic field - add one more inversion to get it right
//-----------------------------------------------------------------------------
TEveVector TMu2eBField::GetField(float X, float Y, float Z) const {

  double bx(0), by(0), bz(0);
  CLHEP::Hep3Vector field = fBmgr->getBField(CLHEP::Hep3Vector(X*10,Y*10,Z*10));

  bx = field.x();
  by = field.y();
  bz = field.z();
  return TEveVector(-bx,-by,-bz);
}
}
