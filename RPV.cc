#include "RPV.h"
RPV rpv;
namespace hak {
	const float &met() { return rpv.met(); }
	const float &metPhi() { return rpv.metPhi(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_p4() { return rpv.jets_p4(); }
	const vector<float> &jets_p4Correction() { return rpv.jets_p4Correction(); }
	const int &type() { return rpv.type(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ll_p4() { return rpv.ll_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lt_p4() { return rpv.lt_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &total_p4() { return rpv.total_p4(); }
	const int &ll_id() { return rpv.ll_id(); }
	const int &lt_id() { return rpv.lt_id(); }
	const int &ll_charge() { return rpv.ll_charge(); }
	const int &lt_charge() { return rpv.lt_charge(); }
	const int &ll_index() { return rpv.ll_index(); }
	const int &lt_index() { return rpv.lt_index(); }
	const int &eventNumber() { return rpv.eventNumber(); }
	const int &runNumber() { return rpv.runNumber(); }
	const int &lumiBlock() { return rpv.lumiBlock(); }
	const float &scale_1fb() { return rpv.scale_1fb(); }
	const vector<float> &btagDiscriminant() { return rpv.btagDiscriminant(); }
	const int &numEvents() { return rpv.numEvents(); }
}
