# Reaction Aware (R)DataFrames

For processing specific electron scattering reactions with podio style data.

Goals :

1. Simplifying analysis code for general final states.
2. Same code for different types of data, e.g. HepMC, ePIC.
3. No additional dependencies (ROOT only), no build (header only, parsed by root at runtime).
4. Simply define the final state particles then use standardised functions to add columns/branches to output, 1 line of code per branch.
5. Hide boilerplate and C++isms from user.
6. Automate MC matching and calculation of equivalent truth variables.
7. Automate combinitorial analysis (!!! To be done)

To install just download the code from git and add the path to ROOT_INCLUDE_PATH
If you do not have the base rad code already installed you can add it as a submodule

      git clone --recurse-submodules https://github.com/dglazier/epic-rad.git
      setenv EPICRAD /to/where/is/epic-rad
      setenv RAD ${EPICRAD}/rad

If RAD is installed already you do not need to download the submodule

      git clone https://github.com/dglazier/epic-rad.git
      setenv EPICRAD /to/where/is/epic-rad
      setenv RAD /to/where/is/rad

In either case you then need to add the include path to ROOT_INCLUDE_PATH so the files are visibel in an interactive root session.

      setenv ROOT_INCLUDE_PATH ${RAD}/include:${EPICRAD}/include
      or
      setenv ROOT_INCLUDE_PATH ${ROOT_INCLUDE_PATH}:${RAD}/include:${EPICRAD}/include
      
  Example code :

        // create an epic reaction
        rad::config::ePICReaction epic{"events", "data_file.root");
        //choose processing scheme i.e. match reconstructed and generated events
        epic.AliasColumnsAndMatchWithMC(false);
        //Assign particles names and indices
        //indicing comes from ordering in hepmc file as we matched Sim and Rec.
        epic.setBeamIonIndex(rad::beams::BeamEleFix());
        epic.setBeamElectronIndex(rad::beams::BeamIonFix());
        epic.setScatElectronIndex(1);
        //give final state hadrons names,
        //if we give a PDG code it will generate el_OK branches etc
        //el_OK = 1 if electron reconstructed with right PDG
        epic.setParticleIndex("el",6,11);
        epic.setParticleIndex("po",7,-11);
        epic.setParticleIndex("p",5,2212);
        
        //Group particles into top and bottom vertices
        //aka Meson and Baryon components
        //this is required for calcualting reaction kinematics
        //e.g. t distributions
        epic.setBaryonParticles({"p"});
        epic.setMesonParticles({"el","po"});

        //must call this after all particles are configured
        epic.makeParticleMap();

        //create column for invariant mass of e+ and e-
        rad::rdf::Mass(epic,"IMass","{el,po}");

        // we can now add further columns, make a snapshot or draw a histogram
        // draw a histogram. Not I must prepend rec_ or tru_ to get the reconstructed or truth variable
        auto df0 = epic.CurrFrame(); //get the current dataframe node. Now operate like regular RDataFrame
        auto hInvMassRec = df0.Histo1D({"InvMassRec","Recon M(e-,e+) [GeV]",100,.3,5.},"rec_IMass");
        auto hInvMassTru = df0.Histo1D({"InvMassTru","Truth M(e-,e+) [GeV]",100,.3,5.},"tru_IMass");
        hInvMassRec->DrawCopy();
        hInvMassTrue->DrawCopy("same");

The matching generated with reconstructed is the simplest analysis for simulated data. However to be more 
realistic you need to add algorithms for choosing which particle is associated with your defined final state particles.
Some examples of this are given in the examples! e.g. choose the first electron in ReconstructedParticles for the scattered electron, 
or choose the electron with the highest momentum. Ultimately this will require full combinitiral analysis to be implemented.

To find the mcmatch index, you need to know the order of the particles in the hepmc3 file. This can be found by checking the MCParticles branch in the reconstructed tree. Open the file in root and get events tree,

      events->Scan("MCParticles.PDG:MCParticles.generatorStatus")

      ***********************************************
      *    Row   * Instance * MCParticl * MCParticl *
      ***********************************************
      *        0 *        0 *        11 *         4 *
      *        0 *        1 *      2212 *         4 *
      *        0 *        2 *        11 *         1 *
      *        0 *        3 *       -11 *         1 *
      *        0 *        4 *       211 *         1 *
      *        0 *        5 *      2112 *         1 *
      *        0 *        6 *        11 *         1 *
      *        0 *        7 *        11 *         0 *
      *        0 *        8 *        22 *         0 *
      *        0 *        9 *        11 *         0 *
      *        0 *       10 *        11 *         0 *
      *        0 *       11 *        22 *         0 *
      *        0 *       12 *        22 *         0 *
      *        0 *       13 *        22 *         0 *


The particles with generatorStatus = 4 are the beams. generatorStatus=1 are the final state particles which have been thrown in genat4. generatorStatus=0 are secondaries which should be ignored. You can match the Status=1 PDG values with the particles in your reaction.



## Developing your own column calculations

To create the user-friendly function rad::rdf::Mass etc, requires 2 steps. Currently this is organised in 2 seperate files.
One for the raw C++ calculation, the other to interface this to RDataFrame via a Define call. When developing your own 
calculations you should try and group them in physics processes, for example a file for compton scattering kinematic calculations.
Lets look at an example, MissMass : given some final state particles, these are subtracted from the sum of the beams and the 
resulting mass is returned. First I must define the c++ function (see include/ReactionKinematics.h),

    template<typename Tp, typename Tm>
    Tp MissMass(const config::RVecIndexMap& react,const RVecI &ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    { 
      auto psum = beams::BeamIonFourVector(react[names::BeamIonIdx()][0],px,py,pz,m);
      psum+=beams::BeamEleFourVector(react[names::BeamEleIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.M();
    }

Here we see some C++ stuff that we want to hide from users, like templating the vectors, this protects against their types changing 
from float to double for example. The type Tp and Tm are deduced at run time and the types of momentum and mass arrays can be different.
To sum the beams we start with the ion beams:: means this function is defined in Beams.h, BeamIonFourVector returns either a fixed 
4-vector which you must have defined at the start of your script, or if you define your beam with an indice, the value given for that 
event, this can be useful for processing simulated or generated data.
The function SubtractFourVector is part of BasicKinematics.h and it just subtracted the four-momentum components indiced in ineg,
which the user will define in their script.
Note the object react which is RVecIndex map allows you to find the indice for specific parts of your final state. The approriate functions to use as a key in the map are given in [DefineNames](https://github.com/dglazier/rad/blob/master/include/DefineNames.h) . So you can replace BeamIonIdx with any of the other Idx functions to get the indices for those particles.

Second in ReactionKinemticsRDF.h we interface to RDataFrame :

    void MissMass(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::MissMass(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }

Note components_p4 => px,py,pz,pm. Using components_p4 allows DefineForAllTypes to switch in the approriate componenets for rec or truth. 

Here use of DefineForAllTypes adds columns for both rec_ and tru_ variables. "name" will be the name of the new column, and neg
is the list of particles to be subtracted e.g. "{el,po}" . rad is able to use this to find the actual index of the electron and 
positron for the event and subtract those particles.

When creating these 2 files you should adhere to the namespacing convention. c++ functions are in namespace rad, Rdataframe 
interfaces are in namespace rad::rdf.
 

# Snapshot Tree

If you create a snapshot via something like

      epic.Snapshot("output.root");

Then the tree will contain all aliased or nnewly defined columns. In particular momentum components and any calculation which was requested. The tree will just be flat single entry per calculation, and multi-entries for momentum components, one for each particle. If you are using MCMatching there will be both a truth branch and reconstructed, allowing you to determine resolutions of all quantities etc. To access a particular particles components you just need to index by the name you gave it e.g.

      //plot the reconstructed momentum of the electron
      rad_tree->Draw("rec_pmag[el]>>p(100,0,20)");
      //plot the truth W
      rad_tree->Draw("tru_W>>w(100,0,50)");



## Matching detector information

For ePIC analysis we can match the detector information to each particle.
To do this you must use the ePICDetectorReaction rather than just ePICReaction.
To specify what information you want to use (for example to save linked info synchronised to tree) :

       //Add some detector associations
       //From TaggerTrackerTracks use "momentum.x"
       epic.AssociateTracks({"TaggerTrackerTracks"},
   		       {"momentum.x"});
       //From various calorimeters use "energy"	       
       epic.AssociateClusters({"EcalBarrelClusters","EcalBarrelImagingClusters",
       	"EcalBarrelScFiClusters",
      	"EcalEndcapNClusters","EcalEndcapPClusters","EcalEndcapPInsertClusters",
      	"HcalBarrelClusters","HcalEndcapNClusters","LFHCALClusters",
      	"EcalFarForwardZDCClusters","HcalFarForwardZDCClusters"},
    	{"energy"});
 

