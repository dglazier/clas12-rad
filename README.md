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

To run on ifarm it is simplest to use my build

      module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
      module load clas12
      module unload hipo/4.2.0

      setenv CLAS12RAD /work/clas12/dglazier/clas12-rad
      setenv RAD ${CLAS12RAD}/rad
      setenv HIPO ${CLAS12RAD}/hipo/install
      setenv ROOT_INCLUDE_PATH ${RAD}/include:${CLAS12RAD}/include:${HIPO}/include
      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HIPO}/lib

To install on laptop etc, just download the code from git and add the path to ROOT_INCLUDE_PATH
If you do not have the base rad code already installed you can add it via a submodule

      git clone --recurse-submodules https://github.com/dglazier/clas12-rad.git
      setenv CLAS12RAD /to/where/is/clas12-rad
      setenv RAD ${CLAS12RAD}/rad

If RAD is installed already you do not need to download the submodule

      git clone https://github.com/dglazier/epic-rad.git
      setenv CLAS12RAD /to/where/is/clas12-rad
      setenv RAD /to/where/is/rad

We also need access to the hipo4 library. A submodule has also been added for this in case you do not already have access. In principle the hipo aspects should already be availbale on ifarm.

To build hipo see

https://github.com/gavalian/hipo?tab=readme-ov-file#installing-the-package-c

      setenv HIPO ${CLAS12RAD}/hipo/install


In either case you then need to add the include path to ROOT_INCLUDE_PATH so the files are visible in an interactive root session.

      setenv ROOT_INCLUDE_PATH ${RAD}/include:${CLAS12RAD}/include:${HIPO}/include
      or
      setenv ROOT_INCLUDE_PATH ${ROOT_INCLUDE_PATH}:${RAD}/include:${CLAS12RAD}/include:${HIPO}/include

And the hipo libdir

      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HIPO}/lib

      
  Example code :

        // create an clas12 reaction
 	rad::clas12::CLAS12Reaction c12{files};	
	//choose processing scheme i.e. match reconstructed and generated events
        c12.AliasColumnsAndMatchWithMC();
	//Set beam and target energy
	c12.FixBeamElectronMomentum(0,0,10.4); //default e- mass
  	c12.FixBeamIonMomentum(0,0,0); //default p mass
        //Assign particles names and indices
        //indicing comes from ordering in hepmc file as we matched Sim and Rec.

	//give final state hadrons names,
	//note ScatElectron name is "scat_ele"
        //if we give a PDG code it will generate el_OK branches etc
        //scat_ele_OK = 1 if electron reconstructed with right PDG
	c12.setScatElectronIndex(rad::indice::useNthOccurance(1,11),{pidtype});//n=1, pid == 11
  	c12.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{pidtype},211);
 	c12.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{pidtype},-211);
 	c12.setParticleIndex("proton",rad::indice::useNthOccurance(1,2212),{pidtype},2212);
         
        //Group particles into top and bottom vertices
        //aka Meson and Baryon components
        //this is required for calcualting reaction kinematics
        //e.g. t distributions
    	c12.setBaryonParticles({"proton"}); //recoil proton
 	c12.setMesonParticles({"pip","pim"}); //intermediate meson

        //must call this after all particles are configured
 	c12.makeParticleMap();
	
        //create column for invariant mass of pi+ and pi-
        rad::rdf::Mass(c12,"IMass","{pip,pim}");

        // we can now add further columns, make a snapshot or draw a histogram
        // draw a histogram. Not I must prepend rec_ or tru_ to get the reconstructed or truth variable
        auto df0 = c12.CurrFrame(); //get the current dataframe node. Now operate like regular RDataFrame
        auto hInvMassRec = df0.Histo1D({"InvMassRec","Recon M(#pi-,#pi+) [GeV]",100,.3,2.},"rec_IMass");
        auto hInvMassTru = df0.Histo1D({"InvMassTru","Truth M(#pi-,#pi+) [GeV]",100,.3,2.},"tru_IMass");
        hInvMassRec->DrawCopy();
        hInvMassTrue->DrawCopy("same");

	//or save all columns to a root tree file
	c12.Snapshot("tree.root");


The matching generated with reconstructed is the simplest analysis for simulated data. However to be more 
realistic you need to add algorithms for choosing which particle is associated with your defined final state particles.
e.g. choose the first electron in ReconstructedParticles for the scattered electron, 
or choose the electron with the highest momentum. Ultimately this will require full combinitorial analysis to be implemented.

Note some helpful branches are added : rec_pmag , rec_theta and rec_phi and if truth matching is on the same with tru_ and res_, where the latter give the difference between rec and tru.

If you snapshot a tree you can access particular particle elements using their name.
e.g

      rad_tree->Draw("rec_pmag[scat_ele]");
      rad_tree->Draw("rec_theta[pip]");

## Developing your own column calculations

To create the user-friendly function rad::rdf::Mass etc, requires 2 steps. Currently this is organised in 2 seperate files.
One for the raw C++ calculation, the other to interface this to RDataFrame via a Define call. When developing your own 
calculations you should try and group them in physics processes, for example a file for compton scattering kinematic calculations.
Lets look at an example, MissMass : given some final state particles, these are subtracted from the sum of the beams and the 
resulting mass is returned. First I must define the c++ function (see rad/include/ReactionKinematics.h),

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

