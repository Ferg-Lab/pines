/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023 of Nicholas Herringer and Siva Dasetty.

The PINES module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The PINES module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithVirtualAtom.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"
#include "tools/Stopwatch.h"

// -- SD header file for both ANN function and PINES
#include "PINES.h"

#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace std;

namespace PLMD
{
namespace PINES
{

//+PLUMEDOC PINESMOD_COLVAR PINES
/*

******************WARNING: ONLYDIRECT FUNCTIONALITY HAS ONLY BEEN TESTED FOR AR13**********************************

\par Examples

\plumedfile
PINES ...
LABEL=PINES
PRECISION=1000
NLIST
ONLYCROSS
REF_FILE=md.pdb
PINESATOMS=5
ATOMTYPES=NA,CL,OW,HW1,HW2
SFACTOR=1.0,1.0,1.0,1.0,1.0
SORT=1.0,1.0,1.0,1.0,1.0
SWITCH1={RATIONAL R_0=0.38 MM=12 NN=6}
SWITCH2={RATIONAL R_0=0.35 MM=12 NN=6}
SWITCH3={RATIONAL R_0=0.41 MM=12 NN=6}
SWITCH4={RATIONAL R_0=0.40 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.28 MM=12 NN=6}
NL_CUTOFF=5.0,1.5,1.5,1.5,1.5
NL_STRIDE=10.0,10.0,10.0,10.0,10.0
NL_SKIN=0.1,0.1,0.1,0.1,0.1
PINESREP
NL_CONSTANT_SIZE=10
... PINES
PRINT STRIDE=1 ARG=PINES FILE=PINES_trj
\plumedfile

This module is designed for compatability with the ANN module. The following is an example of a plumed.dat file using both modules:
\plumedfile
PINES ...
LABEL=PINES
PRECISION=1000
NLIST
ONLYCROSS
REF_FILE=md.pdb
PINESATOMS=5
ATOMTYPES=NA,CL,OW,HW1,HW2
SFACTOR=1.0,1.0,1.0,1.0,1.0
SORT=1.0,1.0,1.0,1.0,1.0
SWITCH1={RATIONAL R_0=0.38 MM=12 NN=6}
SWITCH2={RATIONAL R_0=0.35 MM=12 NN=6}
SWITCH3={RATIONAL R_0=0.41 MM=12 NN=6}
SWITCH4={RATIONAL R_0=0.40 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.28 MM=12 NN=6}
NL_CUTOFF=5.0,1.5,1.5,1.5,1.5
NL_STRIDE=10.0,10.0,10.0,10.0,10.0
NL_SKIN=0.1,0.1,0.1,0.1,0.1
PINESREP
NL_CONSTANT_SIZE=10
... PINES

ANN ...
LABEL=ann
ARG=PINES.*
NUM_LAYERS=6
NUM_NODES=61,64,32,16,3,3
ACTIVATIONS=BNTanh,BNTanh,BNTanh,BNTanh,Linear
WEIGHTS0=                                                                                                        1,1           Top
WEIGHTS1=
WEIGHTS2=
BIASES0=
BIASES1=
BIASES2=
... ANN

PBMETAD ...
ARG=ann.node-0,ann.node-1,ann.node-2 SIGMA=0.1,0.1,0.1 HEIGHT=1.2
TEMP=300 PACE=500 BIASFACTOR=2 LABEL=pb
FILE=hills_pc0.out,hills_pc1.out,hills_pc2.out
GRID_MIN=-2.5,-2.5,-2.5 GRID_MAX=2.5,2.5,2.5 GRID_BIN=5000,5000,5000
... PBMETAD
PRINT STRIDE=1 ARG=ann.node-0,ann.node-1,ann.node-2,pb.bias FILE=colvar.out
\plumedfile


*/
//+ENDPLUMEDOC


PLUMED_REGISTER_ACTION(PINES,"PINES")

void PINES::registerKeywords( Keywords& keys )
{
  Colvar::registerKeywords( keys );
  keys.add("numbered","SWITCH","The switching functions parameter."
           "You should specify a Switching function for all PINES blocks."
           "Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add("compulsory","PRECISION","the precision for approximating reals with integers in sorting.");
  keys.add("compulsory","REF_FILE","PDB file name that contains the information about system connectivity and labels.");
  keys.add("compulsory","PINESATOMS","Number of atoms to use for PINES.");
  keys.add("compulsory","SORT","Whether to sort or not the PINES block.");
  keys.add("compulsory","ATOMTYPES","The atom types to use for PINES.");
  keys.add("optional","SFACTOR","Scale the PINES-distance by such block-specific factor");
  keys.add("optional","VOLUME","Scale atom-atom distances by the cubic root of the cell volume. The input volume is used to scale the R_0 value of the switching function. ");
  keys.add("optional","UPDATEPINES","Frequency (in steps) at which the PINES is updated.");
  keys.addFlag("TEST",false,"Print the actual and reference PINES and exit");
  keys.addFlag("ONLYCROSS",false,"Use only cross-terms (A-B, A-C, B-C, ...) in PINES");
  keys.addFlag("ONLYDIRECT",false,"Use only direct-terms (A-A, B-B, C-C, ...) in PINES");
  keys.addFlag("DERIVATIVES",false,"Activate the calculation of the PINES for every class (needed for numerical derivatives).");
  keys.addFlag("NLIST",false,"Use a neighbor list for distance calculations.");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("TIMER",false,"Perform timing analysis on heavy loops.");
  keys.add("optional","NL_CONSTANT_SIZE","Fix the number of elements in all blocks to be constant. Blocks which have a total number of possible elements less than the chosen constant size will not be affected.");
  keys.add("optional","NL_CUTOFF","Neighbor lists cutoff.");
  keys.add("optional","NL_STRIDE","Update neighbor lists every NL_STRIDE steps. When post-processing a trajectory a stride of 1 should always be used.");
  keys.add("optional","NL_SKIN","The maximum atom displacement tolerated for the neighbor lists update.");
  // -- SD Flag for writing PINES values in a single file when using plumed driver.
  keys.addFlag("WRITEPINESTRAJ",false,"Flag to enable or disable writing PINES_representation when using plumed driver.");
  // -- SD Variables to control frequency of writing PINES values and ANN PINES derivatives during simulation.
  keys.add("optional","WRITEPINESSTRIDE","STRIDE to write PINES_representation.");
  componentsAreNotOptional(keys);
  // Changing "COMPONENTS" to "default" and slightly modifying the name. Added components for ANN_SUM_DERIV
  keys.addOutputComponent("ELEMENT", "default", "Elements of the PINES block. The position in the N choose 2 interactions (i) and the neighbor in the neighbor list (j) is given as PINES-i-j.");
  //keys.addOutputComponent("ANNSUMDERIV", "default", "2D array of PINES element partial derivatives (used with ANN module).");
  keys.reset_style("SWITCH","compulsory");
}

PINES::PINES(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  timer(false),
  NL_const_size(0),
  updatePINES(1),
  Nprec(1000),
  Natm(1),
  Nlist(1),
  NLsize(1),
  Fvol(1.),
  Vol0(0.),
  m_PINESdistance(0.),
  solv_blocks(1),
  scaling(std:: vector<double>(Nlist)),
  r00(std:: vector<double>(Nlist)),
  nl_skin(std:: vector<double>(Nlist)),
  fmass(std:: vector<double>(Nlist)),
  dosort(std:: vector<bool>(Nlist)),
  compos(std:: vector<Vector>(NLsize)),
  sw(std:: vector<string>(Nlist)),
  nl(std:: vector<NeighborList *>(Nlist)),
  nl_small(std:: vector<NeighborList *>(Nlist)),
  nlcom(std:: vector<NeighborList *>(NLsize)),
  m_deriv(std:: vector<Vector>(1)),
  ds_array(std:: vector<double>(1)),
  ann_deriv(std:: vector<std:: vector<Vector> >(1)),
  PINES_Pair0(std:: vector<int>(1)),
  PINES_Pair1(std:: vector<int>(1)),
  Plist(std:: vector<std:: vector<AtomNumber> >(1)),
  listall(std:: vector<AtomNumber>(1)),
  AtomToResID_Dict(std:: vector<unsigned>(1)),
  NList_OW_blocks(std:: vector<unsigned>(1)),
  NList_HW_blocks(std:: vector<unsigned>(1)),
  listreduced(std:: vector<AtomNumber>(1)),
  listnonwater(std:: vector<AtomNumber>(1)),
  atype(std:: vector<string>(1)),
  nl_cut(std:: vector<double>(Nlist)),
  nl_st(std:: vector<int>(Nlist)),
  Svol(false),
  cross(true),
  direct(true),
  doneigh(false),
  test(false),
  CompDer(false),
  // SD -- local variables corresponding to user defined flags.
  writePINEStraj(false),
  writestride(false),
  writePINESstride(-1),
  cart2PINES(true),
  // SD -- used in prepare function.
  invalidateList(true),
  firsttime(true)
{
  log << "Starting PINES Constructor\n";

  // Precision on the real-to-integer transformation for the sorting
  parse("PRECISION",Nprec);
  if(Nprec<2) error("Precision must be => 2");

  // PBC
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc) {
    log << "Using Periodic Boundary Conditions\n";
  } else  {
    log << "Isolated System (NO PBC)\n";
  }

  // SERIAL/PARALLEL
  parseFlag("SERIAL",serial);
  if(serial) {
    log << "Serial PINES construction\n";
  } else     {
    log << "Parallel PINES construction\n";
  }

  // Derivatives
  parseFlag("DERIVATIVES",CompDer);
  if(CompDer) log << "Computing Derivatives\n";

  // Timing
  parseFlag("TIMER",timer);
  if(timer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }

  // Test
  parseFlag("TEST",test);

  // Constant Neighbor List Size
  if(keywords.exists("NL_CONSTANT_SIZE")) {
    parse("NL_CONSTANT_SIZE",NL_const_size);
  }

  // UPDATEPINES
  if(keywords.exists("UPDATEPINES")) {
    parse("UPDATEPINES",updatePINES);
  }

  // Volume Scaling
  parse("VOLUME",Vol0);
  if (Vol0>0) {
    Svol=true;
  }

  // PINES direct and cross blocks
  bool oc=false,od=false;
  parseFlag("ONLYCROSS",oc);
  parseFlag("ONLYDIRECT",od);
  if (oc&&od) {
    error("ONLYCROSS and ONLYDIRECT are incompatible options!");
  }
  if(oc) {
    direct=false;
    log << "Using only CROSS-PINES blocks\n";
  }
  if(od) {
    cross=false;
    log << "Using only DIRECT-PINES blocks\n";
  }

  // Atoms for PINES
  parse("PINESATOMS",Natm);
  atype.resize(Natm);
  parseVector("ATOMTYPES",atype);
  solv_blocks=0;
  if (std::find(atype.begin(), atype.end(), "OW") != atype.end()) {
    solv_blocks=1;
  }
  if ( (std::find(atype.begin(), atype.end(), "HW1") != atype.end()) && (std::find(atype.begin(), atype.end(), "HW2") != atype.end()) ) {
    solv_blocks=2;
  }
  if ( (std::find(atype.begin(), atype.end(), "OW") != atype.end()) && (std::find(atype.begin(), atype.end(), "HW1") != atype.end()) && (std::find(atype.begin(), atype.end(), "HW2") != atype.end()) ) {
    solv_blocks=3;
  }
  // If included, HW1 and HW2 atomtypes are combined into a single HW
  // atomtype so Natm must be decreased by 1 to account for this
  if (solv_blocks > 1) {
    Natm = Natm - 1; 
  }
  // Specific ordering of solvent atomtypes in input file is necessary
  // to avoid complicating the code. This checks to make sure that the
  // expected order is being followed in the input file.
  if (solv_blocks == 1) {
    if (atype[atype.size()-1] != "OW") {
      error("Water atoms must be listed last for ATOMTYPES and must be listed in order as (OW), (HW1,HW2), or (OW,HW1,HW2)");
    }
  } else if (solv_blocks == 2) {
    if ((atype[atype.size()-2] != "HW1") || (atype[atype.size()-1] != "HW2")) {
      error("Water atoms must be listed last for ATOMTYPES and must be listed in order as (OW), (HW1,HW2), or (OW,HW1,HW2)");
    }
  } else if (solv_blocks == 3) {
    if ((atype[atype.size()-3] != "OW") || (atype[atype.size()-2] != "HW1") || (atype[atype.size()-1] != "HW2")) {
      error("Water atoms must be listed last for ATOMTYPES and must be listed in order as (OW), (HW1,HW2), or (OW,HW1,HW2)");
    }
  }

  // Reference PDB file
  parse("REF_FILE",ref_file);
  PDB mypdb;
  FILE* fp=fopen(ref_file.c_str(),"r");
  if (fp!=NULL) {
    log<<"Opening PDB file with reference frame: "<<ref_file.c_str()<<"\n";
    mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    fclose (fp);
  } else {
    error("Error in reference PDB file");
  }

  // Build COM/Atom lists of AtomNumbers (this might be done in PBC.cpp)
  // Atomlist or Plist used to build pair lists
  Plist.resize(Natm);
  // Atomlist used to build list of atoms for each COM
  std:: vector<std:: vector<AtomNumber> > comatm(1);
  // NLsize is the number of atoms in the pdb cell
  NLsize=mypdb.getAtomNumbers().size();
  // In the following P stands for Point (either an Atom or a COM)
  unsigned resnum=0;
  // Presind (array size: number of residues) contains the contains the residue number
  //   this is because the residue numbers may not always be ordered from 1 to resnum
  std:: vector<unsigned> Presind;
  // Build Presind
  for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
    unsigned rind=mypdb.getResidueNumber(mypdb.getAtomNumbers()[i]);
    bool oldres=false;
    for (unsigned j=0; j<Presind.size(); j++) {
      if(rind==Presind[j]) {
        oldres=true;
      }
    }
    if(!oldres) {
      Presind.push_back(rind);
    }
  }
  resnum=Presind.size();

  // Pind0 is the atom/COM used in Nlists (for COM Pind0 is the first atom in the pdb belonging to that COM)
  unsigned Pind0size;

  Pind0size=NLsize;
  
  std:: vector<unsigned> Pind0(Pind0size);

  comatm.resize(NLsize);

  log << "Total COM/Atoms: " << Natm*resnum << " \n";
  // Build lists of Atoms/COMs for NLists
  //   comatm filled also for non_COM calculation for analysis purposes
  unsigned countIndex = 0;
  for (unsigned j=0; j<Natm; j++) {
    unsigned oind;
    for (unsigned i=0; i<Pind0.size(); i++) {
      Pind0[i]=0;
    }
    for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
      // Residue/Atom AtomNumber: used to build NL for COMS/Atoms pairs.
      AtomNumber anum=mypdb.getAtomNumbers()[i];
      // ResidueName/Atomname associated to atom
      string rname=mypdb.getResidueName(anum);
      string aname=mypdb.getAtomName(anum);
      // Index associated to residue/atom: used to separate COM-lists
      unsigned rind=mypdb.getResidueNumber(anum);
      unsigned aind=anum.index();
      // This builds lists for NL
      string Pname;
      unsigned Pind;

      Pname=aname;
      Pind=aind;
      
      if(Pname==atype[j]) {
        if(Pind0[Pind]==0) {
          // adding the atomnumber to the atom/COM list for pairs
          // SD local variable of type AtomNumber. Its value is set using countIndex.
          AtomNumber ati;
          ati.setIndex(countIndex);
          Plist[j].push_back(ati); //anum) -- SD (previously, it is same as atom number in PDB file).;
          Pind0[Pind]=aind+1;
          oind=Pind;
          countIndex += 1;
        }
        // adding the atomnumber to list of atoms for every COM/Atoms
        comatm[Pind0[Pind]-1].push_back(anum);
      }
      // This block accounts for the combination of HW1 and HW2 so that
      // HW2 atom IDs are still included
      else if( (solv_blocks > 1) && (j == Natm-1) ) {
        if(Pname==atype[atype.size()-1]) {
          if(Pind0[Pind]==0) {
            // adding the atomnumber to the atom/COM list for pairs
            // SD local variable of type AtomNumber. Its value is set using countIndex.
            AtomNumber ati;
            ati.setIndex(countIndex);
            Plist[j].push_back(ati); //anum) -- SD (previously, it is same as atom number in PDB file).;
            Pind0[Pind]=aind+1;
            oind=Pind;
            countIndex += 1;
          }
          // adding the atomnumber to list of atoms for every COM/Atoms
          comatm[Pind0[Pind]-1].push_back(anum);
        }
      }
    }
    // Output Lists
    log << "  Groups of type  " << j << ": " << Plist[j].size() << " \n";
    string gname;
    unsigned gsize;

    gname=mypdb.getAtomName(comatm[Pind0[oind]-1][0]);
    gsize=1;

    log.printf("    %6s %3s %13s %10i %6s\n", "type  ", gname.c_str(),"   containing ",gsize," atoms");
  }

  // SD This is to build the list with the atoms required for PINES.
  // std:: vector<AtomNumber> listall;
  listall.clear();
  // AtomToResID is used to ensure that all relevant atoms from a water
  // molecule are considered by matching their residue ID. AtomToResID
  // uses the same indexing as listall for easy comparison between atom
  // IDs and their corresponding residue ID.
  AtomToResID_Dict.clear();
  if (solv_blocks > 1) {
    for (unsigned j=0; j<Natm+1; j++) {
      for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
        // SD --- including only user defined atom types;
        AtomNumber at_num = mypdb.getAtomNumbers()[i];
        unsigned rind=mypdb.getResidueNumber(mypdb.getAtomNumbers()[i]);                                                                        
        // ResidueName/Atomname associated to atom                                                                        
        string at_name = mypdb.getAtomName(at_num);                                                                             
        if(at_name == atype[j]) {
          // -- SD listall should contain the actual atom numbers in the PDB file.
          listall.push_back(at_num);
          AtomToResID_Dict.push_back(rind);
        }                                                                                                             
      }                                                                                                                 
    }
  } else {
    for (unsigned j=0; j<Natm; j++) {
      for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
        // SD --- including only user defined atom types;
        AtomNumber at_num = mypdb.getAtomNumbers()[i];                                                                        
        // ResidueName/Atomname associated to atom                                                                        
        string at_name = mypdb.getAtomName(at_num);                                                                             
        if(at_name == atype[j]) {
          // -- SD listall should contain the actual atom numbers in the PDB file.
          listall.push_back(at_num);
        } 
      }                                                                                                                 
    }
  }

  // PINES blocks and Neighbour Lists
  Nlist=0;
  // Direct adds the A-A ad B-B blocks (N)
  if(direct) {
    Nlist=Nlist+unsigned(Natm);
  }
  // Cross adds the A-B blocks (N*(N-1)/2)
  if(cross) {
    Nlist=Nlist+unsigned(double(Natm*(Natm-1))/2.);
    // If water oxygens and water hydrogens are included
    // there will be an undesirable interaction block corresponding to
    // OW-HW. This removes that block from the total block count (Nlist)
    if (solv_blocks == 3) {
      Nlist = Nlist - 1;
    }
  }

  unsigned ncnt=0;
  // NL_const_size limits the size of water oxygen blocks (also can be
  // thought of as the number of water molecules), but there should be
  // twice as many hydrogen considered as oxygen (2:1 hydrogen/oxygen
  // ratio in water). This helps with this later by specifying which
  // block numbers out of Nlist correspond to which block types 
  NList_OW_blocks.clear();
  NList_HW_blocks.clear();
  if(cross){
    for (unsigned j=0; j<Natm; j++) {
      for (unsigned i=j+1; i<Natm; i++) {
        if (ncnt < Nlist) {
          if (solv_blocks == 1) {
            if (i == Natm-1) {
              NList_OW_blocks.push_back(ncnt);
            }
          } else if (solv_blocks == 2) {
            if (i == Natm-1) {
              NList_HW_blocks.push_back(ncnt);
            }
          } else if (solv_blocks == 3) {
            if (i == Natm-1) {
              NList_HW_blocks.push_back(ncnt);
            } else if (i == Natm-2) {
              NList_OW_blocks.push_back(ncnt);
            }
          }
        }
        ncnt+=1;
      }
    }
  }
  // PINES scaled option
  scaling.resize(Nlist);
  for(unsigned j=0; j<Nlist; j++) {
    scaling[j]=1.;
  }

  if(keywords.exists("SFACTOR")) {
    parseVector("SFACTOR",scaling);
  }

  // Added STRIDE to write PINES representation and ANN sum derivatives -- SD
  if(keywords.exists("WRITEPINESTRAJ")){
    parseFlag("WRITEPINESTRAJ",writePINEStraj);
  }
  if(keywords.exists("WRITEPINESSTRIDE")) { 
    parse("WRITEPINESSTRIDE",writePINESstride);
  }
  if (writePINESstride != -1) {
    writestride=true;
  }

  // Neighbour Lists option
  parseFlag("NLIST",doneigh);
  nl.resize(Nlist);
  nl_small.resize(Nlist);
  nl_skin.resize(Nlist);
  nl_cut.resize(Nlist,0.);
  nl_st.resize(Nlist,0);
  if(doneigh) {
    parseVector("NL_CUTOFF",nl_cut);
    parseVector("NL_STRIDE",nl_st);
    parseVector("NL_SKIN",nl_skin);
    for (unsigned j=0; j<Nlist; j++) {
      if(nl_cut[j]<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
      if(nl_st[j]<=0) error("NL_STRIDE should be explicitly specified and positive");
      if(nl_skin[j]<=0.) error("NL_SKIN should be explicitly specified and positive");
      nl_cut[j]=nl_cut[j]+nl_skin[j];
    }
    log << "Creating Neighbor Lists \n";
    // SD -- nlall is a neighbor list created using list all. nl_cut[0] and nl_st[0] are probably not needed.
    // WARNING: is nl_cut meaningful here?
    nlall= new NeighborList(listall,true,pbc,getPbc(),comm,nl_cut[0],nl_st[0]);
    unsigned ncnt=0;
    // Direct blocks AA, BB, CC, ...
    if(direct) {
      for (unsigned j=0; j<Natm; j++) {
        nl[ncnt]= new NeighborList(Plist[j],true,pbc,getPbc(),comm,nl_cut[j],nl_st[j]);
        ncnt+=1;
      }
    }

    // Cross blocks AB, AC, BC, ...
    if(cross) {
      // -- SD No changes here. nl depends on Plist for each j. Plist for each j = [0, Num of atoms of type j]
      // -- SD example: if j=0 corresponds to C1, Plist[0] = [0] because there is only one C1; if is OW, it is [0, 1, 2, ... Nwaters]
      for (unsigned j=0; j<Natm; j++) {
        for (unsigned i=j+1; i<Natm; i++) {
          // This accounts for the case where Nlist != Natm*(Natm-1)/2 (when solv_blocks = 3)
          if (ncnt < Nlist) {
            nl[ncnt]= new NeighborList(Plist[i],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[ncnt],nl_st[ncnt]);
          }
          ncnt+=1;
        }
      }

    }

    if(getStep()==0) {
      std::vector<vector<AtomNumber>> snn_list;
      std::vector<vector<AtomNumber>> Hydrogen_nn_list;
      std::vector<double> snn_mags;
      std::vector<AtomNumber> all_ids;
      std::vector<double> all_mags;
      std::vector<AtomNumber> total_waterOx_list;
      std::vector<AtomNumber> total_waterHydrogen_list;
      int Nlist_count;

      // Sized for Natm - solvent atom 
      if (solv_blocks == 3) {
        Hydrogen_nn_list.resize(Natm-2);
        snn_list.resize(Natm-2);
      } else {
        Hydrogen_nn_list.resize(Natm-1);
        snn_list.resize(Natm-1);
      }

      Vector ddist;
      double smallest_val;
      Nlist_count=0;

      for(unsigned j=0; j<Natm; j++) {
        for(unsigned i=j+1; i<Natm; i++) {
          if( (atype[i] == "OW") || ((solv_blocks == 2) && (atype[i] == "HW1")) ) {
            all_mags.clear();
            all_ids.clear();
            
            for(unsigned k=0; k<nl[Nlist_count]->size(); k++) {
              unsigned id0=(nl[Nlist_count]->getClosePairAtomNumber(k).first).index();
              unsigned id1=(nl[Nlist_count]->getClosePairAtomNumber(k).second).index();
              Vector Pos0,Pos1;
              double mag;
              
              Pos0=mypdb.getPosition(listall[id0]);
              Pos1=mypdb.getPosition(listall[id1]);

              ddist=pbcDistance(Pos0,Pos1);
              mag=ddist.modulo();

              all_mags.push_back(mag);
              all_ids.push_back(listall[id0]);
            }
            for(unsigned x=0; x<NL_const_size+10; x++) {
              smallest_val = all_mags[0];
              int nl_pos = 0;
              for(unsigned y=0; y<all_mags.size(); y++) {
                if(smallest_val > all_mags[y]) {
                  smallest_val = all_mags[y];
                  nl_pos = y;
                }
              }
              snn_mags.push_back(smallest_val);
              snn_list[j].push_back(all_ids[nl_pos]);

              // If water hydrogen are included, 
              // find them by residue ID
              unsigned atom_indx = listall.size();
              if (solv_blocks > 1) {
                for(unsigned k=0; k<listall.size(); k++) {
                  if (listall[k] == all_ids[nl_pos]) {
                    atom_indx = k;
                  }
                }
                for(unsigned k=0; k<listall.size(); k++) {
                  if (AtomToResID_Dict[atom_indx] == AtomToResID_Dict[k]) {
                    if( (solv_blocks == 3) && (atom_indx != k) ){
                      Hydrogen_nn_list[j].push_back(listall[k]);
                    } else if ( (solv_blocks == 2) ) {
                      Hydrogen_nn_list[j].push_back(listall[k]);
                    }
                  }
                }
              }

              all_mags[nl_pos] = 99999;
              if(total_waterOx_list.empty()) {
                total_waterOx_list.push_back(all_ids[nl_pos]);
              } else {
                bool id_present=false;
                for(unsigned y=0; y<total_waterOx_list.size(); y++) {
                  if(all_ids[nl_pos] == total_waterOx_list[y]) {
                    id_present=true;
                  }
                }
                if(!id_present) {
                  total_waterOx_list.push_back(all_ids[nl_pos]);
                }
              }
              // If hydrogen IDs aren't in the total ID list, add them
              if (solv_blocks > 1) {
                if (std::find(total_waterHydrogen_list.begin(), total_waterHydrogen_list.end(), Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-2]) == total_waterHydrogen_list.end()) {
                  total_waterHydrogen_list.push_back(Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-2]);
                }
                if (std::find(total_waterHydrogen_list.begin(), total_waterHydrogen_list.end(), Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-1]) == total_waterHydrogen_list.end()) {
                  total_waterHydrogen_list.push_back(Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-1]);
                }
              }
            }
          }
          Nlist_count += 1;
        }
      }
      // // fclose(debugging_file);
      if (solv_blocks > 0) {
        std::sort(total_waterOx_list.begin(), total_waterOx_list.end());
      }
      listreduced.clear();
      listnonwater.clear();
      for (unsigned j=0; j<Natm; j++) {
        for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
          // SD --- including only user defined atom types;
          AtomNumber at_num = mypdb.getAtomNumbers()[i];                                                                        
          // ResidueName/Atomname associated to atom                                                                        
          string at_name = mypdb.getAtomName(at_num);                                                                             
          if( (at_name == atype[j]) && ( (at_name != "OW") && (at_name != "HW1") && (at_name != "HW2") ) ) {
            // -- SD listall should contain the actual atom numbers in the PDB file.
            listnonwater.push_back(at_num);
          }
          if( (direct) && ( (at_name == "C12") || (at_name == "C43") || (at_name == "C56") ) ) {
            listnonwater.push_back(at_num);
          }                                                                                                               
        }                                                                                                                 
      }
      listreduced = listnonwater;
      // Add solvent IDs to the total ID list
      if (solv_blocks == 1) {
        for(unsigned i=0; i<total_waterOx_list.size(); i++) {
          listreduced.push_back(total_waterOx_list[i]);
        }
      } else if (solv_blocks == 2) {
        for(unsigned i=0; i<total_waterHydrogen_list.size(); i++) {
          listreduced.push_back(total_waterHydrogen_list[i]);
        }
      } else if (solv_blocks == 3) {
        for(unsigned i=0; i<total_waterOx_list.size(); i++) {
          listreduced.push_back(total_waterOx_list[i]);
          listreduced.push_back(total_waterHydrogen_list[i*2]);
          listreduced.push_back(total_waterHydrogen_list[i*2+1]);
        }
      }
      // Match atom IDs to positional indexing in listreduced    
      if (solv_blocks == 3) { 
        for(unsigned j=0; j<Natm-2; j++) {
          for(unsigned k=0; k<NL_const_size+10; k++) {
            for(unsigned i=0; i<listreduced.size(); i++) {
              if(snn_list[j][k]==listreduced[i]){
                AtomNumber atom_id;
                snn_list[j][k]=atom_id.setIndex(i);
              }
              if(Hydrogen_nn_list[j][2*k]==listreduced[i]){
                AtomNumber atom_id;
                Hydrogen_nn_list[j][2*k]=atom_id.setIndex(i);
              }
              if(Hydrogen_nn_list[j][2*k+1]==listreduced[i]){
                AtomNumber atom_id;
                Hydrogen_nn_list[j][2*k+1]=atom_id.setIndex(i);
              }
            }
          }
        }
      } else {
        for(unsigned j=0; j<Natm-1; j++) {
          for(unsigned k=0; k<NL_const_size+10; k++) {
            for(unsigned i=0; i<listreduced.size(); i++) {
              if (solv_blocks == 1) {
                if(snn_list[j][k]==listreduced[i]){
                  AtomNumber atom_id;
                  snn_list[j][k]=atom_id.setIndex(i);
                }
              }
              if (solv_blocks == 2) {
                if(Hydrogen_nn_list[j][2*k]==listreduced[i]){
                  AtomNumber atom_id;
                  Hydrogen_nn_list[j][2*k]=atom_id.setIndex(i);
                }
                if(Hydrogen_nn_list[j][2*k+1]==listreduced[i]){
                  AtomNumber atom_id;
                  Hydrogen_nn_list[j][2*k+1]=atom_id.setIndex(i);
                }
              }
            }
          }
        }
      }
      int count_nl_loop=0;
      if(cross) {
        for(unsigned j=0; j<Natm-1; j++) {
          for(unsigned i=j+1; i<Natm; i++) {
            if (count_nl_loop < Nlist) {
              if(atype[i] == "OW") {
                nl_small[count_nl_loop] = new NeighborList(snn_list[j],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[count_nl_loop],nl_st[count_nl_loop]);
              } else if (atype[i] == "HW1") {
                nl_small[count_nl_loop] = new NeighborList(Hydrogen_nn_list[j],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[count_nl_loop],nl_st[count_nl_loop]);
              } else {
                nl_small[count_nl_loop] = new NeighborList(Plist[i],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[ncnt],nl_st[ncnt]);
              }
            }
            count_nl_loop += 1;
          }
        }
      }
      else if(direct) {
        for(unsigned j=0; j<Plist.size(); j++) {
          nl_small[count_nl_loop] = new NeighborList(Plist[j],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[ncnt],nl_st[ncnt]);
        }
      }
      nlreduced= new NeighborList(listreduced,true,pbc,getPbc(),comm,nl_cut[0],nl_st[0]);
    }
  } else {
    log << "WARNING: Neighbor List not activated this has not been tested!!  \n";
    nlall= new NeighborList(listall,true,pbc,getPbc(),comm);
    for (unsigned j=0; j<Nlist; j++) {
      nl[j]= new NeighborList(Plist[j],Plist[j],true,true,pbc,getPbc(),comm);
    }
  }
  // Output Nlist
  log << "Total Nlists: " << Nlist << " \n";
  for (unsigned j=0; j<Nlist; j++) {
    log << "  list " << j+1 << "   size " << nl[j]->size() << " \n";
  }


  // Sorting
  dosort.resize(Nlist);
  std:: vector<int> ynsort(Nlist);
  parseVector("SORT",ynsort);
  if(cart2PINES) {
    for (unsigned i=0; i<Nlist; i++) {
      if(ynsort[i]==0) {
        dosort[i]=false;
      } else {
        dosort[i]=true;
      }
    }
  } else {
    for (unsigned i=0; i<Nlist; i++) {
      if(ynsort[i]==0||CompDer) {
        dosort[i]=false;
      } else {
        dosort[i]=true;
      }
    }
  }
  //build box vectors and correct for pbc
  log << "Building the box from PDB data ... \n";
  Tensor Box=mypdb.getBoxVec();
  log << "  Done! A,B,C vectors in Cartesian space:  \n";
  log.printf("  A:  %12.6f%12.6f%12.6f\n", Box[0][0],Box[0][1],Box[0][2]);
  log.printf("  B:  %12.6f%12.6f%12.6f\n", Box[1][0],Box[1][1],Box[1][2]);
  log.printf("  C:  %12.6f%12.6f%12.6f\n", Box[2][0],Box[2][1],Box[2][2]);
  log << "Changing the PBC according to the new box \n";
  Pbc mypbc;
  mypbc.setBox(Box);
  log << "The box volume is " << mypbc.getBox().determinant() << " \n";

  //Compute scaling factor
  if(Svol) {
    Fvol=cbrt(Vol0/mypbc.getBox().determinant());
    log << "Scaling atom distances by  " << Fvol << " \n";
  } else {
    log << "Using unscaled atom distances \n";
  }

  r00.resize(Nlist);
  sw.resize(Nlist);
  for (unsigned j=0; j<Nlist; j++) {
    if( !parseNumbered( "SWITCH", j+1, sw[j] ) ) break;
  }
  if(CompDer) {
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    sfs.resize(Nlist);
    std::string errors;
    for (unsigned j=0; j<Nlist; j++) {
      if(Svol) {
        double r0;
        vector<string> data=Tools::getWords(sw[j]);
        data.erase(data.begin());
        Tools::parse(data,"R_0",r0);
        std::string old_r0; Tools::convert(r0,old_r0);
        r0*=Fvol;
        std::string new_r0; Tools::convert(r0,new_r0);
        std::size_t pos = sw[j].find("R_0");
        sw[j].replace(pos+4,old_r0.size(),new_r0);
      }
      sfs[j].set(sw[j],errors);
      std::string num;
      Tools::convert(j+1, num);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
      r00[j]=sfs[j].get_r0();
      log << "  Swf: " << j << "  r0=" << (sfs[j].description()).c_str() << " \n";
    }
  }


  // build the rPINES distances (transformation and sorting is done afterwards)
  if(CompDer) {
    log << "  PINES  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
  }
  checkRead();
  // Create components of PINES

  // cPINES hasn't been created yet so it can't be used in the loop.
  // The loop is set up generally as N(N-1)/2. 
  // The if/else statement accounts for there being >1 elements in the solute-solvent blocks, 
  // and expects that the interaction with solvent is the last block of the each solute atom's interactions
  // Total count keeps a running tally of the elements in the entire PINES so that there will be an equal number of components.
  unsigned total_count=0;
  int count_nl_loop=0;
  if(cross) {
  for(int j = 0; j < Natm; j++) {
    for(int i= j+1; i < Natm; i++) {
      if (count_nl_loop < Nlist) {
        // Add elements for the various solvent blocks if needed
        if(atype[i] == "OW") {
          for(int n = 0; n < NL_const_size; n++) {
            string comp = "ELEMENT-" + to_string(total_count);
            addComponentWithDerivatives(comp); 
            componentIsNotPeriodic(comp);
            total_count += 1;
          }
        } else if(atype[i] == "HW1") {
          for(int n = 0; n < 2*NL_const_size; n++) {
            string comp = "ELEMENT-" + to_string(total_count);
            addComponentWithDerivatives(comp); 
            componentIsNotPeriodic(comp);
            total_count += 1;
          }
        } else {
          string comp = "ELEMENT-" + to_string(total_count);
          addComponentWithDerivatives(comp); 
          componentIsNotPeriodic(comp);
          total_count +=1;
        }
      }
      count_nl_loop +=1;
    }
  }
  }
  else if(direct) {
    for(int j = 0; j < 78; j++) {
      string comp = "ELEMENT-" + to_string(total_count);
      addComponentWithDerivatives(comp);
      componentIsNotPeriodic(comp);
      total_count +=1;
    }
  }


  requestAtoms(nlreduced->getFullAtomList());


  if(cross){
    ann_deriv.resize(listreduced.size());
  }
  if(direct){
    ann_deriv.resize(listreduced.size()*(listreduced.size()-1)*0.5);
  }

 // printf("TOTAL COUNT :%d\n\n\n", total_count);
  for (int i = 0; i < ann_deriv.size(); i++) {
    ann_deriv[i].resize(total_count);
  }

  ds_array.resize(total_count);
}

// The following deallocates pointers
PINES::~PINES()
{
  for (unsigned j=0; j<Nlist; j++) {
    delete nl[j];
  }
  for(unsigned j=0; j<Nlist; j++) {
    delete nl_small[j];
  }
  delete nlall;
  delete nlreduced;
}


// SD request atoms in every frame.
void PINES::prepare() {
  if(nlall->getStride()>0) {
    if((getStep()+1)%nlall->getStride()==0) {
      requestAtoms(nlall->getFullAtomList());

      if(cross){
        ann_deriv.resize(listall.size());
      }
      if (direct){
	      ann_deriv.resize(listall.size()*(listall.size()-1)*0.5);
      }

      int total_count=0;
      // Adjust total_count based on the solvent atoms considered.
      // Total_count should be the total number of elements in the
      // PINES block.

      if (solv_blocks == 3) {
        for(unsigned j=0; j<Natm-2; j++) {
          total_count += j;
        }
        total_count += NL_const_size*(Natm-2)*3;
      } else if (solv_blocks == 2) {
        for(unsigned j=0; j<Natm-1; j++) {
          total_count += j;
        }
        total_count += NL_const_size*(Natm-1)*2;
      } else if (solv_blocks == 1) {
        for(unsigned j=0; j<Natm-1; j++) {
          total_count += j;
        }
        total_count += NL_const_size*(Natm-1);
      } else {
        for(unsigned j=0; j<Natm; j++) {
          total_count += j;
        }
      }

      if(direct){
	      total_count = 78;
      }

      for(unsigned i=0; i < ann_deriv.size(); i++) {
        ann_deriv[i].resize(total_count);
      }
      ds_array.resize(total_count);
    } else if((getStep()%nlall->getStride()==0) && (getStep()!=0)) {
      
      std::vector<vector<AtomNumber>> snn_list;
      std::vector<vector<AtomNumber>> Hydrogen_nn_list;
      std::vector<double> snn_mags;
      std::vector<AtomNumber> all_ids;
      std::vector<double> all_mags;
      std::vector<AtomNumber> total_waterOx_list;
      std::vector<AtomNumber> total_waterHydrogen_list;

      int buffer=10;
      int Nlist_count;
      snn_list.clear();
      Hydrogen_nn_list.clear();
      total_waterOx_list.clear();
      total_waterHydrogen_list.clear();

      if (solv_blocks == 3) {
        Hydrogen_nn_list.resize(Natm-2);
        snn_list.resize(Natm-2);
      } else {
        Hydrogen_nn_list.resize(Natm-1);
        snn_list.resize(Natm-1);
      }

      Vector ddist;
      double smallest_val;
      Nlist_count=0;

      Nlist_count=0;
      for(unsigned j=0; j<Natm; j++) {
        for(unsigned i=j+1; i<Natm; i++) {
          if( (atype[i] == "OW") || ((solv_blocks == 2) && (atype[i] == "HW1")) ) {
            all_mags.clear();
            all_ids.clear();

            for(unsigned k=0; k<nl[Nlist_count]->size(); k++) {
              unsigned id0=(nl[Nlist_count]->getClosePairAtomNumber(k).first).index();
              unsigned id1=(nl[Nlist_count]->getClosePairAtomNumber(k).second).index();
              Vector Pos0,Pos1;
              double mag;
              
              Pos0=getPosition(id0);
              Pos1=getPosition(id1);

              ddist=pbcDistance(Pos0,Pos1);
              mag=ddist.modulo();

              all_mags.push_back(mag);
              all_ids.push_back(listall[id0]);
            }
            for(unsigned x=0; x<NL_const_size+10; x++) {
              smallest_val = all_mags[0];
              int nl_pos = 0;
              for(unsigned y=0; y<all_mags.size(); y++) {
                if(smallest_val > all_mags[y]) {
                  smallest_val = all_mags[y];
                  nl_pos = y;
                }
              }

              snn_mags.push_back(smallest_val);
              snn_list[j].push_back(all_ids[nl_pos]);

              unsigned atom_indx = listall.size();
              if (solv_blocks > 1) {
                for(unsigned k=0; k<listall.size(); k++) {
                  if (listall[k] == all_ids[nl_pos]) {
                    atom_indx = k;
                  }
                }
                for(unsigned k=0; k<listall.size(); k++) {
                  if (AtomToResID_Dict[atom_indx] == AtomToResID_Dict[k]) {
                    if( (solv_blocks == 3) && (atom_indx != k) ){
                      Hydrogen_nn_list[j].push_back(listall[k]);
                    } else if ( (solv_blocks == 2) ) {
                      Hydrogen_nn_list[j].push_back(listall[k]);
                    }
                  }
                }
              }


              all_mags[nl_pos] = 99999.;
              if(total_waterOx_list.empty()) {
                total_waterOx_list.push_back(all_ids[nl_pos]);
              } else {
                bool id_present=false;
                for(unsigned y=0; y<total_waterOx_list.size(); y++) {
                  if(all_ids[nl_pos] == total_waterOx_list[y]) {
                    id_present=true;
                  }
                }
                if(!id_present) {
                  total_waterOx_list.push_back(all_ids[nl_pos]);
                }
              }
              if (solv_blocks > 1) {
                if (std::find(total_waterHydrogen_list.begin(), total_waterHydrogen_list.end(), Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-2]) == total_waterHydrogen_list.end()) {
                  total_waterHydrogen_list.push_back(Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-2]);
                }
                if (std::find(total_waterHydrogen_list.begin(), total_waterHydrogen_list.end(), Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-1]) == total_waterHydrogen_list.end()) {
                  total_waterHydrogen_list.push_back(Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-1]);
                }
              }
            }
          }
          Nlist_count += 1;
        }
      }
      if (solv_blocks > 0) {
        std::sort(total_waterOx_list.begin(), total_waterOx_list.end());
      }
      listreduced.clear();
      listreduced = listnonwater;

      if (solv_blocks == 1) {
        for(unsigned i=0; i<total_waterOx_list.size(); i++) {
          listreduced.push_back(total_waterOx_list[i]);
        }
      } else if (solv_blocks == 2) {
        for(unsigned i=0; i<total_waterHydrogen_list.size(); i++) {
          listreduced.push_back(total_waterHydrogen_list[i]);
        }
      } else if (solv_blocks == 3) {
        for(unsigned i=0; i<total_waterOx_list.size(); i++) {
          listreduced.push_back(total_waterOx_list[i]);
          listreduced.push_back(total_waterHydrogen_list[i*2]);
          listreduced.push_back(total_waterHydrogen_list[i*2+1]);
        }
      }
  
      if (solv_blocks == 3) { 
        for(unsigned j=0; j<Natm-2; j++) {
          for(unsigned k=0; k<NL_const_size+10; k++) {
            for(unsigned i=0; i<listreduced.size(); i++) {
              if(snn_list[j][k]==listreduced[i]){
                AtomNumber atom_id;
                snn_list[j][k]=atom_id.setIndex(i);
              }
              if(Hydrogen_nn_list[j][2*k]==listreduced[i]){
                AtomNumber atom_id;
                Hydrogen_nn_list[j][2*k]=atom_id.setIndex(i);
              }
              if(Hydrogen_nn_list[j][2*k+1]==listreduced[i]){
                AtomNumber atom_id;
                Hydrogen_nn_list[j][2*k+1]=atom_id.setIndex(i);
              }
            }
          }
        }
      } else {
        for(unsigned j=0; j<Natm-1; j++) {
          for(unsigned k=0; k<NL_const_size+10; k++) {
            for(unsigned i=0; i<listreduced.size(); i++) {
              if (solv_blocks == 1) {
                if(snn_list[j][k]==listreduced[i]){
                  AtomNumber atom_id;
                  snn_list[j][k]=atom_id.setIndex(i);
                }
              }
              if (solv_blocks == 2) {
                if(Hydrogen_nn_list[j][2*k]==listreduced[i]){
                  AtomNumber atom_id;
                  Hydrogen_nn_list[j][2*k]=atom_id.setIndex(i);
                }
                if(Hydrogen_nn_list[j][2*k+1]==listreduced[i]){
                  AtomNumber atom_id;
                  Hydrogen_nn_list[j][2*k+1]=atom_id.setIndex(i);
                }
              }
            }
          }
        }
      }
      int count_nl_loop=0;
      for(unsigned j=0; j<Natm-1; j++) {
        for(unsigned i=j+1; i<Natm; i++) {
          if (count_nl_loop < Nlist) {
            // Only update solvent blocks (included IDs are dynamic)
            if(atype[i] == "OW") {
              delete nl_small[count_nl_loop];
              nl_small[count_nl_loop] = new NeighborList(snn_list[j],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[count_nl_loop],nl_st[count_nl_loop]);
            } else if (atype[i] == "HW1") {
              delete nl_small[count_nl_loop];
              nl_small[count_nl_loop] = new NeighborList(Hydrogen_nn_list[j],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[count_nl_loop],nl_st[count_nl_loop]);
            }
          }
          count_nl_loop += 1;
        }
      }

      
      delete nlreduced;
      nlreduced= new NeighborList(listreduced,true,pbc,getPbc(),comm,nl_cut[0],nl_st[0]);

      requestAtoms(nlreduced->getFullAtomList());
      if(cross){
        ann_deriv.resize(listreduced.size());
      }
      if(direct){
	      ann_deriv.resize(listreduced.size()*(listreduced.size()-1)*0.5);
      }
     // printf("ANN size %d",ann_deriv.size());
      int total_count=0;

      if(cross){
        if (solv_blocks == 3) {
          for(unsigned j=0; j<Natm-2; j++) {
            total_count += j;
          }
          total_count += NL_const_size*(Natm-2)*3;
        } else if (solv_blocks == 2) {
          for(unsigned j=0; j<Natm-1; j++) {
            total_count += j;
          }
          total_count += NL_const_size*(Natm-1)*2;
        } else if (solv_blocks == 1) {
          for(unsigned j=0; j<Natm-1; j++) {
            total_count += j;
          }
          total_count += NL_const_size*(Natm-1);
        } else {
          for(unsigned j=0; j<Natm; j++) {
            total_count += j;
          }
        }
      }
      if(direct){
        total_count = 78;
      }

      for(unsigned i=0; i < ann_deriv.size(); i++) {
        ann_deriv[i].resize(total_count);
      }
      ds_array.resize(total_count);
    }
  }
}

void PINES::calculate()
{

  // Local variables
  // The following are probably needed as static arrays
  static int prev_stp=-1;
  static int init_stp=1;
  static std:: vector<std:: vector<Vector> > prev_pos(Nlist);
  static std:: vector<std:: vector<double> > cPINES(Nlist);
  static std:: vector<std:: vector<int> > Atom0(Nlist);
  static std:: vector<std:: vector<int> > Atom1(Nlist);
  std:: vector<std:: vector<int> > A0(Nprec);
  std:: vector<std:: vector<int> > A1(Nprec);
  size_t stride=1;
  unsigned rank=0;

  if(!serial) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  } else {
    stride=1;
    rank=0;
  }

  // Transform (and sort) the rPINES before starting the dynamics
  if (((prev_stp==-1) || (init_stp==1)) &&!CompDer) {
    if(prev_stp!=-1) {init_stp=0;}
    // Calculate the volume scaling factor
    if(Svol) {
      Fvol=cbrt(Vol0/getBox().determinant());
    }
    //Set switching function parameters
    log << "\n";
    log << "REFERENCE PDB # " << prev_stp+2 << " \n";
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    sfs.resize(Nlist);
    std::string errors;
    for (unsigned j=0; j<Nlist; j++) {
      if(Svol) {
        double r0;
        vector<string> data=Tools::getWords(sw[j]);
        data.erase(data.begin());
        Tools::parse(data,"R_0",r0);
        std::string old_r0; Tools::convert(r0,old_r0);
        r0*=Fvol;
        std::string new_r0; Tools::convert(r0,new_r0);
        std::size_t pos = sw[j].find("R_0");
        sw[j].replace(pos+4,old_r0.size(),new_r0);
      }
      sfs[j].set(sw[j],errors);
      std::string num;
      Tools::convert(j+1, num);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
      r00[j]=sfs[j].get_r0();
      log << "  Swf: " << j << "  r0=" << (sfs[j].description()).c_str() << " \n";
    }
  }
  // Do the sorting only once per timestep to avoid building the PINES N times for N rPINES PDB structures!
  if ((getStep()>prev_stp&&getStep()%updatePINES==0)||CompDer) {
    if (CompDer) log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PINES \n";
    //
    // build COMs from positions if requested
    // update neighbor lists when an atom moves out of the Neighbor list skin
    if (doneigh && ((getStep()+1)%nlall->getStride()==0)) {
      bool doupdate=false;
      // For the first step build previous positions = actual positions
      if (prev_stp==-1) {
        bool docom=com;
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if(docom) {
              Pos=compos[i];
            } else {
              Pos=getPosition(nl[j]->getFullAtomList()[i].index());
            }
            prev_pos[j].push_back(Pos);
          }
        }
        doupdate=true;
      }
      // Decide whether to update lists based on atom displacement, every stride
      std:: vector<std:: vector<Vector> > tmp_pos(Nlist);
      if (getStep() % nlall->getStride() ==0) {
        bool docom=com;
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if(docom) {
              Pos=compos[i];
            } else {
              Pos=getPosition(nl[j]->getFullAtomList()[i].index());
            }
            tmp_pos[j].push_back(Pos);
            if (pbcDistance(tmp_pos[j][i],prev_pos[j][i]).modulo()>=nl_skin[j]) {
              doupdate=true;
            }
          }
        }
      }
      // Update Nlists if needed
      if (doupdate==true) {
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            prev_pos[j][i]=tmp_pos[j][i];
          }
          nl[j]->update(prev_pos[j]);
          log << " Step " << getStep() << "  Neighbour lists updated " << nl[j]->size() << " \n";
        }
      }
    }
    // Calculate the volume scaling factor
    if(Svol) {
      Fvol=cbrt(Vol0/getBox().determinant());
    }
    Vector ddist;
    // Global to local variables
    bool doserial=serial;

    // Build "Nlist" PINES blocks
    for(unsigned j=0; j<Nlist; j++) {
      if(dosort[j]) {
        // from global to local variables to speedup the for loop with if statements
        bool docom=com;
        bool dopbc=pbc;
        // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
        std:: vector<int> OrdVec(Nprec,0);
        cPINES[j].resize(0);
        Atom0[j].resize(0);
        Atom1[j].resize(0);
        // Building distances for the PINES vector at time t
        if(timer) stopwatch.start("1 Build cPINES");
        if((getStep()+1)%nlall->getStride()==0) {
          for(unsigned i=rank; i<nl[j]->size(); i+=stride) {
            unsigned i0=(nl[j]->getClosePairAtomNumber(i).first).index();
            unsigned i1=(nl[j]->getClosePairAtomNumber(i).second).index();\
            Vector Pos0,Pos1;
            if(docom) {
              Pos0=compos[i0];
              Pos1=compos[i1];
            } else {
              Pos0=getPosition(i0);
              Pos1=getPosition(i1);
            }
            if(dopbc) {
              ddist=pbcDistance(Pos0,Pos1);
            } else {
              ddist=delta(Pos0,Pos1);
            }
            double df=0.;
            //Integer sorting ... faster!
            //Transforming distances with the Switching function + real to integer transformation
            int Vint=int(sfs[j].calculate(ddist.modulo()*Fvol, df)*double(Nprec-1)+0.5);
            // Enables low precision with standard PINES sizes.
            if(cart2PINES) {
              if(Vint == 0) {
                Vint = 1;
              }
            }

            //Integer transformed distance values as index of the Ordering Vector OrdVec
            OrdVec[Vint]+=1;
            //Keeps track of atom indices for force and virial calculations
            A0[Vint].push_back(i0);
            A1[Vint].push_back(i1);
          }
          int sb_count=0;
          int discards=0;
          // Account for twice as many hydrogen as oxygen
          int max_solv_atoms = NL_const_size;
          if (std::find(NList_HW_blocks.begin(), NList_HW_blocks.end(), j) != NList_HW_blocks.end()) {
            max_solv_atoms = 2*NL_const_size;
          }

          if(nl[j]->size() >= max_solv_atoms) {
            // Limit solvent blocks to NL_constant_size
            for(unsigned i=Nprec-1; i<0; i--) {
              if(OrdVec[i] != 0) {
                sb_count += OrdVec[i];
              }
              if(sb_count > max_solv_atoms) {
                OrdVec[i] = 0;
                if(sb_count - A0[i].size() < max_solv_atoms) {
                  discards = sb_count - max_solv_atoms;
                  A0[i].resize(A0[i].size()-discards);
                  A1[i].resize(A1[i].size()-discards);
                } else {
                  A0[i].resize(1,0);
                  A1[i].resize(1,0);
                }
              }
            }
          }
        } else {
          for(unsigned i=rank; i<nl_small[j]->size(); i+=stride) {
            unsigned i0=(nl_small[j]->getClosePairAtomNumber(i).first).index();
            unsigned i1=(nl_small[j]->getClosePairAtomNumber(i).second).index();
            Vector Pos0,Pos1;
            if(docom) {
              Pos0=compos[i0];
              Pos1=compos[i1];
            } else {
              Pos0=getPosition(i0);
              Pos1=getPosition(i1);
            }
            if(dopbc) {
              ddist=pbcDistance(Pos0,Pos1);
            } else {
              ddist=delta(Pos0,Pos1);
            }
            double df=0.;
            //Integer sorting ... faster!
            //Transforming distances with the Switching function + real to integer transformation
            int Vint=int(sfs[j].calculate(ddist.modulo()*Fvol, df)*double(Nprec-1)+0.5);
            // Enables low precision with standard PINES sizes.
            if(cart2PINES) {
              if(Vint == 0) {
                Vint = 1;
              }
            }

            //Integer transformed distance values as index of the Ordering Vector OrdVec
            OrdVec[Vint]+=1;
            //Keeps track of atom indices for force and virial calculations
            A0[Vint].push_back(i0);
            A1[Vint].push_back(i1);
          }
          int sb_count=0;
          int discards=0;

          int max_solv_atoms = NL_const_size;
          if (std::find(NList_HW_blocks.begin(), NList_HW_blocks.end(), j) != NList_HW_blocks.end()) {
            max_solv_atoms = 2*NL_const_size;
          }

          if(nl_small[j]->size() >= max_solv_atoms) {
            // Limit solvent blocks to NL_constant_size
            for(unsigned i=Nprec-1; i<0; i--) {
              if(OrdVec[i] != 0) {
                sb_count += OrdVec[i];
              }
              if(sb_count > max_solv_atoms) {
                OrdVec[i] = 0;
                if(sb_count - A0[i].size() < max_solv_atoms) {
                  discards = sb_count - max_solv_atoms;
                  A0[i].resize(A0[i].size()-discards);
                  A1[i].resize(A1[i].size()-discards);
                } else {
                  A0[i].resize(1,0);
                  A1[i].resize(1,0);
                }
              }
            }
          }
        }
        if(timer) stopwatch.stop("1 Build cPINES");
        if(timer) stopwatch.start("2 Sort cPINES");
        if(!doserial && comm.initialized()) {
          // Vectors keeping track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          std:: vector<int> Vdim(stride,0);
          std:: vector<int> Vpos(stride,0);
          // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
          std:: vector<int> OrdVecAll(stride*Nprec);
          // Big vectors containing all Atom indexes for every occupancy (Atom0O(Nprec,n) and Atom1O(Nprec,n) matrices in one vector)
          std:: vector<int> Atom0F;
          std:: vector<int> Atom1F;
          // Vector used to reconstruct arrays
          std:: vector<unsigned> k(stride,0);
          // Zeros might be many, this slows down a lot due to MPI communication
          // Avoid passing the zeros (i=1) for atom indices
          for(unsigned i=1; i<Nprec; i++) {
            // Building long vectors with all atom indexes for occupancies ordered from i=1 to i=Nprec-1
            // Can this be avoided ???
            Atom0F.insert(Atom0F.end(),A0[i].begin(),A0[i].end());
            Atom1F.insert(Atom1F.end(),A1[i].begin(),A1[i].end());
            A0[i].resize(0);
            A1[i].resize(0);
          }
          // Resize partial arrays to fill up for the next PINES block
          A0[0].resize(0);
          A1[0].resize(0);
          A0[Nprec-1].resize(0);
          A1[Nprec-1].resize(0);
          // Avoid passing the zeros (i=1) for atom indices
          OrdVec[0]=0;
          OrdVec[Nprec-1]=0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();

          // pass the array sizes before passing the arrays
          int dim=Atom0F.size();
          // Vdim and Vpos keep track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          comm.Allgather(&dim,1,&Vdim[0],1);

          // TO BE IMPROVED: the following may be done by the rank 0 (now every rank does it)
          int Fdim=0;
          for(unsigned i=1; i<stride; i++) {
            Vpos[i]=Vpos[i-1]+Vdim[i-1];
            Fdim+=Vdim[i];
          }
          Fdim+=Vdim[0];
          // build big vectors for atom pairs on all ranks for all ranks
          std:: vector<int> Atom0FAll(Fdim);
          std:: vector<int> Atom1FAll(Fdim);
          // TO BE IMPROVED: Allgathers may be substituted by gathers by proc 0
          //   Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PINES
          comm.Allgather(&OrdVec[0],Nprec,&OrdVecAll[0],Nprec);
          // Gather the vectors of atom pairs to keep track of the idexes for the forces
          comm.Allgatherv(&Atom0F[0],Atom0F.size(),&Atom0FAll[0],&Vdim[0],&Vpos[0]);
          comm.Allgatherv(&Atom1F[0],Atom1F.size(),&Atom1FAll[0],&Vdim[0],&Vpos[0]);

          // Reconstruct the full vectors from collections of Allgathered parts (this is a serial step)
          // This is the tricky serial step, to assemble together PINES and atom-pair info from head-tail big vectors
          // Loop before on l and then on i would be better but the allgather should be modified
          // Loop on blocks
          // Loop on Ordering Vector size excluding zeros (i=1)
          if(timer) stopwatch.stop("2 Sort cPINES");
          if(timer) stopwatch.start("3 Reconstruct cPINES");
          for(unsigned i=1; i<Nprec; i++) {
            // Loop on the ranks
            for(unsigned l=0; l<stride; l++) {
              // Loop on the number of head-to-tail pieces
              for(unsigned m=0; m<OrdVecAll[i+l*Nprec]; m++) {
                // cPINES is the current PINES at time t
                cPINES[j].push_back(double(i)/double(Nprec-1));
                Atom0[j].push_back(Atom0FAll[k[l]+Vpos[l]]);
                Atom1[j].push_back(Atom1FAll[k[l]+Vpos[l]]);
                k[l]+=1;
              }
            }
          }
          if(timer) stopwatch.stop("3 Reconstruct cPINES");
        } else {
          for(unsigned i=1; i<Nprec; i++) {
            for(unsigned m=0; m<OrdVec[i]; m++) {
              cPINES[j].push_back(double(i)/double(Nprec-1));
              Atom0[j].push_back(A0[i][m]);
              Atom1[j].push_back(A1[i][m]);
            }
          }
        }
      }
    }
  }
  Vector distance;
  double dfunc=0.;
  // Calculate volume scaling factor
  if(Svol) {
    Fvol=cbrt(Vol0/getBox().determinant());
  }


                                                    
  FILE *PINES_rep_file = NULL;
  if (writestride) {
    if (getStep() % writePINESstride == 0) {                                                                                      
      string PINES_rep_fileName = "PINES_representation_" + to_string(getStep()) + ".dat";                                    
      PINES_rep_file = fopen(PINES_rep_fileName.c_str(), "w+"); 
    }
  }
  
  FILE *PINES_rep_file_traj = NULL;
  if (writePINEStraj) {
    string PINES_rep_fileName_traj = "PINES_representation_traj.dat";
    PINES_rep_file_traj = fopen(PINES_rep_fileName_traj.c_str(), "a");
  }

  // SD countLoopLimit for debugging.
  int countLoopLimit = 0;
  for(unsigned j=0; j<Nlist; j++) {
    bool dosorting=dosort[j];
    unsigned limit=0;
    // Set limit to size of PINES block. Solute-solute blocks will be 1 (for tetracosane) 
    // and much larger for solute-solvent blocks (likely hundreds)
    limit = cPINES[j].size();
    // Allow for non-constant block sizes if desired
    if(NL_const_size > 0) {
      // Solute-solvent blocks have more neighbors than necessary so that padding is not necessary.
      // The solute-solvent blocks are already sorted so the last elements of the block (size NL_const_size) are the desired interactions to include
      // i.e. the closest solute-solvent interactions. This is irrelevant for the solute-solute interactions. 
      int start_val=0;
      // This sets the start value to be NL_const_size away from the end of the sorted block to choose the desired interactions.
      
      int max_solv_atoms = NL_const_size;
      if (std::find(NList_HW_blocks.begin(), NList_HW_blocks.end(), j) != NList_HW_blocks.end()) {
        max_solv_atoms = 2*NL_const_size;
      }

      if (cross){
        if(limit > max_solv_atoms) {
          start_val = limit - max_solv_atoms;
        }
      }
      if (writePINEStraj) {
        for(unsigned i=start_val; i<limit; i++) {
          fprintf(PINES_rep_file_traj, "%8.6f\t", cPINES[j][i]);
        }
      }
      if (writestride) {
        if ( getStep() % writePINESstride == 0) {
          for(unsigned i=start_val; i<limit; i++) {
            fprintf(PINES_rep_file, "%8.6f\t", cPINES[j][i]);
            countLoopLimit += 1;
          }
        }
      }
    } else {
      // Prints out in the same PINES block element format as TEST
      if (writePINEStraj) {
        for(unsigned i=0; i<limit; i++) {
          fprintf(PINES_rep_file_traj, "%8.6f\t", cPINES[j][i]);
        } 
      }
      if (writestride) {
        if ( getStep() % writePINESstride == 0) {
          for(unsigned i=0; i<limit; i++) {
            fprintf(PINES_rep_file, "%8.6f\t", cPINES[j][i]);
            countLoopLimit += 1;
          }
        }
      }
    }
  }

  if (writePINEStraj) {
    fprintf(PINES_rep_file_traj, "\n#END OF FRAME\n");
    fclose(PINES_rep_file_traj);
  }
  if (writestride) {
    if ( getStep() % writePINESstride == 0) {
      fprintf(PINES_rep_file, "\n#END OF FRAME: %d \n", getStep());
      fclose(PINES_rep_file);
    }
  }

  if(timer) stopwatch.start("4 Build For Derivatives");
  // non-global variables Nder and Scalevol defined to speedup if structures in cycles
  bool Nder=CompDer;
  bool Scalevol=Svol;
  if(cart2PINES) {
    for(unsigned j=0; j<ann_deriv.size(); j++) {
      for(unsigned i=0; i<ann_deriv[j].size(); i++) {
        for(unsigned k=0; k<3; k++) {ann_deriv[j][i][k]=0.;}
      }
    }
    for(unsigned j=0; j<3; j++) {
      for(unsigned k=0; k<3; k++) {
        m_virial[j][k]=0.;
      }
    }
    // resize vectors to the appropriate sizes and set starting values to zero --NH

    PINES_Pair0.resize(ds_array.size());
    PINES_Pair1.resize(ds_array.size());

    unsigned PINES_element=0;
    for(unsigned j=0; j<Nlist; j++) {
      unsigned limit=0;
      limit = cPINES[j].size();
      int start_val=0;

      int max_solv_atoms = NL_const_size;
      if (std::find(NList_HW_blocks.begin(), NList_HW_blocks.end(), j) != NList_HW_blocks.end()) {
        max_solv_atoms = 2*NL_const_size;
      }
      
      if (cross) {
      	if(limit > max_solv_atoms) {
          start_val = limit - max_solv_atoms;
        }
      }
      for(unsigned i=start_val; i<limit; i++) {
        unsigned i0=0;
        unsigned i1=0;
        // Atom0 and Atom1 are lists that index atoms for PINES elements
        i0=Atom0[j][i];
        // Record the atom IDs for the PINES elements of interest --NH
        PINES_Pair0[PINES_element] = i0;
        i1=Atom1[j][i];
        PINES_Pair1[PINES_element] = i1;
        // Pos0 and Pos1 are 1x3 vectors that hold the xyz coordinates of the indexed atoms
        Vector Pos0,Pos1;
        Pos0=getPosition(i0);
        Pos1=getPosition(i1);
        // distance is also a 1x3 vector of the xyz distances between the two atoms after consideration of the pbc
        distance=pbcDistance(Pos0,Pos1);
        dfunc=0.;
        // dm is a scalar value that is the magnitude of the distance between the atoms
        double dm=distance.modulo();
        // sfs[j] is the parameters for the switching function, which can be chosen to be different for different blocks
        // In this case, all blocks use the same switching function so all sfs[j] are the same function.
        // Used with .calculate(dm*Fvol, dfunc), the PINES element value is returned and the derivative stored in dfunc
        double tPINES = sfs[j].calculate(dm*Fvol, dfunc);

        double ds_element=0.;
        // Create the ds_array one element at a time --NH
        ds_element = scaling[j]*Fvol*Fvol*dfunc*dm;
        ds_array[PINES_element] = ds_element;
        // Create 1x3 vector of (dr/dx,dr/dy,dr/dz) --NH
        Vector dr_dcoord = distance/dm;
        
        // the xyz components of the distance between atoms is scaled by tmp and added or subtracted to reflect
        // that distance is calculated as Pos1 - Pos0
        // Calculate ann_deriv values for the current PINES element in the loop --NH
        ann_deriv[i0][PINES_element] = -ds_element*dr_dcoord;
        ann_deriv[i1][PINES_element] =  ds_element*dr_dcoord;

        // This m_virial is likely not correct but has been included in case it is necessary to test the code --NH
        m_virial    -= ds_element*Tensor(distance,distance); // Question
        PINES_element += 1;

      }
    }
    
    if (!serial && comm.initialized() ) {
      int count = 0;
      for(unsigned j=0; j<Nlist; j++) {
          for(unsigned i=0; i<cPINES[j].size(); i++) {
              count += 1;
          }
      }
      
      comm.Barrier();
      // SD -- This probably works because cPINES[j] size is variable for each j.
      for (unsigned j=0; j< Nlist; j++) {
        for (unsigned k=0; k<cPINES[j].size(); k++) {
          comm.Sum(cPINES[j][k]);
          cPINES[j][k] /= comm.Get_size();
        }
      }

      // SD -- This probably works because comm.Sum cannot handle 3D vectors.
      if(!ann_deriv.empty()) {
        for (unsigned i=0;  i<ann_deriv.size(); i++) {
          for (unsigned j=0; j<ann_deriv[j].size(); j++) {
            for (unsigned k=0; k<3; k++) {
              comm.Sum(ann_deriv[i][j][k]);
              ann_deriv[i][j][k] /= comm.Get_size();
            }
          }
        }
      }
      // SD -- this is probably not needed.
      comm.Sum(&m_virial[0][0],9);
    }
  }
  prev_stp=getStep();

  //Timing
  if(timer) stopwatch.stop("4 Build For Derivatives");
  if(timer) {
    log.printf("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log<<stopwatch;
  }
  unsigned total_count=0;
  for (int j = 0; j < Nlist; j++) {
    unsigned limit=0;
    limit = cPINES[j].size();
    int start_val=0;

    int max_solv_atoms = NL_const_size;
    if (std::find(NList_HW_blocks.begin(), NList_HW_blocks.end(), j) != NList_HW_blocks.end()) {
      max_solv_atoms = 2*NL_const_size;
    }

    if (cross) {
      if( (limit > max_solv_atoms) && (solv_blocks != 0) ){
        start_val = limit - max_solv_atoms;
      }
    }
    for (int i = start_val; i < limit; i++) {
      string comp = "ELEMENT-" + to_string(total_count);
      Value* valueNew=getPntrToComponent(comp);
      valueNew -> set(cPINES[j][i]);
      // Pass the 3D array to the plumed core --NH
      // A 2D array is passed for each PINES element (component) --NH
      for(unsigned k=0; k<ann_deriv.size(); k++) {
        setAtomsDerivatives(valueNew, k, ann_deriv[k][total_count]);
      }
      total_count += 1;
    }
  }

}
//Close Namespaces at the very beginning
}
}
