using namespace std;

const double PI = 3.14159265;

const double kEtheta2Cut = 0.004; //GeV rad^2
const double kRazzleCutVal = 0.6;
const vector<double> kPrismCentroid {-74.,0.}; //x,y [cm]
const double kDistanceFromBNB = 110e2; // [cm]
const double nt = 1.53e30; //Number of argon targets in active volume
const double GeV2perm2 = 2.56819e31;

const double kElectronMass   =  5.109989461e-04;        // GeV
const double kMuonMass       =  1.056583745e-01;         // GeV
const double kTauMass        =  1.77686e+00;             // GeV
const double kPionMass       =  1.3957018e-01;          // GeV
const double kPi0Mass        =  1.349766e-01;           // GeV
const double kProtonMass     =  9.38272081e-01;           // GeV
const double kNeutronMass    =  9.39565413e-01;           // GeV
const double kPhotonMass    =  0;           // GeV

map<int,double> pdgmass{
  {11,kElectronMass},
  {13,kMuonMass},
  {15,kTauMass},
  {211,kPionMass},
  {111,kPi0Mass},
  {2212,kProtonMass},
  {2112,kNeutronMass},
  {22,kPhotonMass},
  {0,0}, //No id, set mass to 0
  {5,0}, //No id, set mass to 0
  {-5,0} //No id, set mass to 0
};

const std::vector<int> colors = {kRed,kGreen,kBlue,kYellow,kMagenta,kCyan,kOrange,kSpring,kTeal,kAzure,kViolet,kPink};



