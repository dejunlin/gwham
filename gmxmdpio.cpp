#include "typedefs.hpp"
#include "gmxmdpio.hpp"
#include "fileio.hpp"
#include <string>
#include <map>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iterator>

using namespace std;

mdp2pullpot::mdp2pullpot(const bitset<MAXNRST>& _rcmask, const fileio& _fio) :
  rcmask(_rcmask), 
  fio(_fio)
  {}

linecounter mdp2pullpot::operator() (const string& fname, map<uint,vector<umbrella*> >& funct, vector<Hamiltonian<RSTXLAMBDAsgl>* >& V) {
  fio.fopen(fname);
  //There are 2 possible modes: pull-group or contact-group
  //for pull-group mode, we take k and init from pull-group parameters
  //for contact-group, we take them from contact-group parameters
  //pull-group mode is activated when pull_ncontactgroups is not defined
  //contact-group mode is activated when pull_ncontactgroups is defined
  //contact-group mode will overwrite pull-group mode
  valtype T; //temperature
  vector<valtype> rstlambdas, kA, kB, initA, initB;
  vector<vector<valtype> > contactlambdas;
  bitset<MAXNRST> ckA, ckB, cinitA, cinitB;
  uint current_lambda = 0;
  uint npullgrps = 0;
  uint ncontactgrps = 0;
  linecounter nlines = 0;
  while(fio.readaline()) {
    vector<string> tmpstr = fio.line2str();
    //ignore empty lines
    if(!tmpstr.size()) { continue; }
    ++nlines;
    if(tmpstr[0].compare("init-lambda-state") == 0) {
      current_lambda = atoi(tmpstr[2].c_str());
    } else if(tmpstr[0].compare("restraint-lambdas") == 0) {
      for(uint i = 2; i < tmpstr.size(); ++i) {
	rstlambdas.push_back(atof(tmpstr[i].c_str()));
      }
    } else if(tmpstr[0].compare("umbrella-contact1-lambdas") == 0) {
      if(contactlambdas.size()) {
	cerr << "umbrella-contact1-lambdas in mdp file must be specified before umbrella-contact2-lambdas\n";
	exit(-1);
      }
      vector<valtype> lambdas;
      for(uint i = 2; i < tmpstr.size(); ++i) {
	lambdas.push_back(atof(tmpstr[i].c_str()));
      }
      contactlambdas.push_back(lambdas);
    } else if(tmpstr[0].compare("umbrella-contact2-lambdas") == 0) {
      vector<valtype> lambdas;
      for(uint i = 2; i < tmpstr.size(); ++i) {
	lambdas.push_back(atof(tmpstr[i].c_str()));
      }
      contactlambdas.push_back(lambdas);
    } else if(tmpstr[0].compare("ref_t") == 0) {
      T = atof(tmpstr[2].c_str());
    } else if(tmpstr[0].compare("pull") == 0) {
      if(tmpstr[2].compare("umbrella") != 0 && tmpstr[2].compare("contact") != 0) {
	cerr << "The pulling potential in a mdp file is not umbrella/contact but you're asking it to read umbrella/contact data" << endl;
	exit(-1);
      }
    } else if(tmpstr[0].compare("pull_ngroups") == 0) {
      npullgrps = atoi(tmpstr[2].c_str());
      cout << "#Reading pull potential for " << npullgrps << " pull-groups" << endl;
      //NOTE: Here we resize the parameter arrays
      if(!ncontactgrps) {
        kA.resize(npullgrps); 
        kB.resize(npullgrps); 
        initA.resize(npullgrps); 
        initB.resize(npullgrps); 
      }
    } else if(tmpstr[0].compare("pull_ncontactgroups") == 0) {
      ncontactgrps = atoi(tmpstr[2].c_str());
      cout << "#Reading pull potential for " << ncontactgrps << " contact-groups" << endl;
      //NOTE: Here we resize the parameter arrays 
      kA.resize(ncontactgrps); 
      kB.resize(ncontactgrps); 
      initA.resize(ncontactgrps); 
      initB.resize(ncontactgrps); 
    } else if(tmpstr[0].compare(0,string("pull_group").length(),"pull_group") == 0 && !ncontactgrps) {
      const uint pullgrpid = getpullgrpid("pull_group",tmpstr[0]);
      cout << "#Processing pull-group " << pullgrpid << endl;
    } else if(tmpstr[0].compare(0,string("pull_init").length(),"pull_init") == 0 && !ncontactgrps) {
      if(tmpstr[0].compare(0,string("pull_initB").length(),"pull_initB") == 0) {
        const uint pullgrpid = getpullgrpid("pull_initB",tmpstr[0]);
        initB[pullgrpid-1] = atof(tmpstr[2].c_str());
        cinitB |= (bitunit << pullgrpid-1);
      } else {
        const uint pullgrpid = getpullgrpid("pull_init",tmpstr[0]);
        initA[pullgrpid-1] = atof(tmpstr[2].c_str());
        cinitA |= (bitunit << pullgrpid-1);
      }
    } else if(tmpstr[0].compare(0,string("pull_k").length(),"pull_k") == 0 && !ncontactgrps) {
      if(tmpstr[0].compare(0,string("pull_kB").length(),"pull_kB") == 0) {
        const uint pullgrpid = getpullgrpid("pull_kB",tmpstr[0]);
        kB[pullgrpid-1] = atof(tmpstr[2].c_str());
        ckB |= (bitunit << pullgrpid-1);
      } else {
        const uint pullgrpid = getpullgrpid("pull_k",tmpstr[0]);
        kA[pullgrpid-1] = atof(tmpstr[2].c_str());
        ckA |= (bitunit << pullgrpid-1);
      }
    } else if(tmpstr[0].compare(0,string("pull_contactgroup").length(),"pull_contactgroup") == 0 && ncontactgrps) {
      const uint contactgrpid = getcontactgrpid(tmpstr[0]);
      char contactgrpidstr[MAXNDIGWIN];
      sprintf(contactgrpidstr,"%d",contactgrpid);
      string kbstr = string("pull_contactgroup") + contactgrpidstr + "_kB";
      string kstr = string("pull_contactgroup") + contactgrpidstr + "_k";
      string nc0bstr = string("pull_contactgroup") + contactgrpidstr + "_nc0B";
      string nc0str = string("pull_contactgroup") + contactgrpidstr + "_nc0";
      if(tmpstr[0].compare(0,kbstr.length(), kbstr) == 0) {
        cout << "#Processing contact-group " << contactgrpid << endl;
        kB[contactgrpid] = atof(tmpstr[2].c_str());
        ckB |= (bitunit << contactgrpid);
      } else if(tmpstr[0].compare(0,kstr.length(), kstr) == 0){
        kA[contactgrpid] = atof(tmpstr[2].c_str());
        ckA |= (bitunit << contactgrpid);
      } else if(tmpstr[0].compare(0,nc0bstr.length(), nc0bstr) == 0) {
        initB[contactgrpid] = atof(tmpstr[2].c_str());
        cinitB |= (bitunit << contactgrpid);
      } else if(tmpstr[0].compare(0,nc0str.length(), nc0str) == 0){
        initA[contactgrpid] = atof(tmpstr[2].c_str());
        cinitA |= (bitunit << contactgrpid);
      }
    } 
  }
  //Check if we read the correct number of parameters
  if(!ncontactgrps) {
    if(ckA.count() != npullgrps || cinitA.count() != npullgrps) {
      cerr << "Number of state-A parameters is less than the number pull-groups " << npullgrps << endl;
      exit(-1);
    }
  } else {
    if(ckA.count() != ncontactgrps || cinitA.count() != ncontactgrps) {
      cerr << "Number of state-A parameters is less than the number contact-groups " << ncontactgrps << endl;
      exit(-1);
    }
  }
  //cout << "# pullgrpid kA kB initA initB lambda" << endl;
  const uint maxi = ncontactgrps ? ncontactgrps : npullgrps;
  for(uint i = 0; i < maxi;  ++i) {
    if(kA[i] == kB[i] && kA[i] == 0 && (rcmask & (bitunit << i)).none()) { continue; }
    //Occasionaly, the B-state parameters are not available and we'll assume they're the same as the A-state ones
    if((ckA & (bitunit << i)).any() && (ckB & (bitunit << i)).none()) { kB[i] = kA[i]; } 
    if((cinitA & (bitunit << i)).any() && (cinitB & (bitunit << i)).none()) { initB[i] = initA[i]; }
    valtype l;
    if(!rstlambdas.size() && !contactlambdas.size()) { l = 0;}
    else if (ncontactgrps) { 
      l = contactlambdas[i][current_lambda]; 
    }

    //cout << "# " << i << " " << kA[i] << " " << kB[i] << " " << initA[i] << " " << initB[i] << " " << l << endl;
    //funct.insert(pair<uint,umbrella>(i, umbrella(kA[i],kB[i],initA[i],initB[i],l)));
    funct[i].push_back(new umbrella(kA[i],kB[i],initA[i],initB[i],l));
    //cout << "# " << i << " " << funct[i].getk() << " " << funct[i].getinit() << endl; 
  }
  //Construct Hamiltonian
  vector<valtype> k, init;
  vector<uint> RC_grpids;
  for(uint i = 0; i < maxi; ++i) {
    //For RCs that are interesting, we build RST hamiltonian on them
    //NOTE that even for RCs that are virtually NOT restrained, we still build RST 
    //hamiltonian on them but the parameters are implicitly zero
    if( (rcmask & (bitunit << i)).any() ) {
      if(funct.find(i) == funct.end()) {
	if(ncontactgrps) {
	  cerr << "Contact group " << i << " is interesting but there's no umbrella restraint on it.\n";
	} else {
	  cerr << "Pull group " << i << " is interesting but there's no umbrella restraint on it.\n";
	}
	exit(-1);
      }
      k.push_back(funct[i][funct[i].size()-1]->getk());
      init.push_back(funct[i][funct[i].size()-1]->getinit());
      RC_grpids.push_back(i);
    }
  }
  //RC_grpids.size() is how many RC we're interested in, which will be the first RC_grpids.size() elements of the histogram data point
  //funct.begin()->size() tells how man states we've processed, including the current one
  //so the potential due to the restrained RCs we're not interested will be combined as the whichlambda'th element in one histogram data point
  if(!funct.size()) { cerr << "There's no any restraint applied on any of the RCs you're interested in. Check the input files again\n"; exit(-1); }
  //dumRCsize is how many restrained RC we're NOT interested in
  const uint dumRCsize = funct.size() - RC_grpids.size();
  const uint whichlambda = dumRCsize ? (RC_grpids.size() + (funct.begin()->second).size() - 1) : 0;
  vector<valtype> params;
  params.push_back(BoltzmannkJ);
  params.push_back(T);
  params.insert(params.end(),k.begin(),k.end());
  params.insert(params.end(),init.begin(),init.end());
  params.push_back(dumRCsize ? 1.0 : 0.0); //should just weight it to 0 if dumRCsize == 0
  params.push_back(whichlambda); 
  V.push_back(new Hamiltonian<RSTXLAMBDAsgl>(params));
  cout << "#mdp2pullpot processed " << RC_grpids.size() << " RCs you're interested in, whose id (indexed begins from zero) are: ";
  copy(RC_grpids.begin(),RC_grpids.end(),ostream_iterator<uint>(cout," "));
  cout << " and they'll be binned into histogram.\n";
  cout << "#The rest " << dumRCsize << " restrained RCs will be combined into their total bias energy, which will be binned into histogram\n";
  cout << "#We're now at the " << (funct.begin()->second).size() << "'th hamiltonian\n";
  return nlines;
}

uint mdp2pullpot::getpullgrpid(const string& directive, const string& str) const {
  const uint size_dir = directive.length();
  const uint size_tot = str.length();
  const uint size_num = size_tot - size_dir;
  const string numstr = str.substr(size_dir,size_num);
  /*cout << "dir = " << directive << endl;
  cout << "str = " << str << endl;
  cout << "size_num = " << size_num << endl;*/
  return atoi(numstr.c_str());
}

uint mdp2pullpot::getcontactgrpid(const string& str) const {
  vector<string> out;
  split<string>(out, str, "_");
  return getpullgrpid("contactgroup",out[1]);
}

linecounter mdp2pullpot::operator() (const string& fname, map<uint,vector<umbrella_fb*> >& funct, vector<Hamiltonian<RST_fbXLAMBDAsgl>* >& V) {
  fio.fopen(fname);
  linecounter nlines = 0;
  vector<valtype> lambdas, k0A, k0B, k1A, k1B, initA, initB, r0A, r0B, r1A, r1B;
  // counters to tell if the number of parameters read in is correct
  bitset<MAXNRST> ck0A, ck0B, ck1A, ck1B, cinitA, cinitB, cr0A, cr0B, cr1A, cr1B;
  valtype T; //temperature
  uint current_lambda = 0;
  uint npullgrps = 0;
  while(fio.readaline()) {
    vector<string> tmpstr = fio.line2str();
    //ignore empty lines
    if(!tmpstr.size()) { continue; }
    ++nlines;
    if(tmpstr[0].compare("init-lambda-state") == 0) {
      current_lambda = atoi(tmpstr[2].c_str());
    } else if(tmpstr[0].compare("restraint-lambdas") == 0) {
      for(uint i = 2; i < tmpstr.size(); ++i) {
	lambdas.push_back(atof(tmpstr[i].c_str()));
      }
    } else if(tmpstr[0].compare("ref_t") == 0) {
      T = atof(tmpstr[2].c_str());
    } else if(tmpstr[0].compare("pull") == 0) {
      if(tmpstr[2].compare("umbrella-flat-bottom") != 0) {
	cerr << "The pulling potential in a mdp file is not umbrella-flat-bottom but you're asking it to read umbrella-flat-bottom data" << endl;
	exit(-1);
      }
    } else if(tmpstr[0].compare("pull_ngroups") == 0) {
      npullgrps = atoi(tmpstr[2].c_str());
      cout << "#Reading pull potential for " << npullgrps << " pull-groups" << endl;
      //NOTE: Here we resize the parameter arrays 
      k0A.resize(npullgrps,0); 
      k0B.resize(npullgrps,0); 
      k1A.resize(npullgrps,0); 
      k1B.resize(npullgrps,0); 
      initA.resize(npullgrps,0); 
      initB.resize(npullgrps,0); 
      r0A.resize(npullgrps,0); 
      r0B.resize(npullgrps,0); 
      r1A.resize(npullgrps,0); 
      r1B.resize(npullgrps,0); 
    } else if(tmpstr[0].compare(0,string("pull_group").length(),"pull_group") == 0) {
      const uint pullgrpid = getpullgrpid("pull_group",tmpstr[0]);
      if(!pullgrpid) { cout << endl << "#Processing pull-group "; }
      cout <<  pullgrpid << " ";
    } else if(tmpstr[0].compare(0,string("pull_umb_k0").length(),"pull_umb_k0") == 0) {
      if(tmpstr[0].compare(0,string("pull_umb_k0B").length(),"pull_umb_k0B") == 0) {
        const uint pullgrpid = getpullgrpid("pull_umb_k0B",tmpstr[0]);
        k0B[pullgrpid-1] = atof(tmpstr[2].c_str());
        ck0B |= (bitunit << pullgrpid-1);
	//cout << "# " << pullgrpid-1 << " k0B = " << k0B[pullgrpid-1] << endl;
      } else {
        const uint pullgrpid = getpullgrpid("pull_umb_k0",tmpstr[0]);
        k0A[pullgrpid-1] = atof(tmpstr[2].c_str());
        ck0A |= (bitunit << pullgrpid-1);
	//cout << "# " << pullgrpid-1 << " k0A = " << k0A[pullgrpid-1] << endl;
      }
    } else if(tmpstr[0].compare(0,string("pull_umb_k1").length(),"pull_umb_k1") == 0) {
      if(tmpstr[0].compare(0,string("pull_umb_k1B").length(),"pull_umb_k1B") == 0) {
        const uint pullgrpid = getpullgrpid("pull_umb_k1B",tmpstr[0]);
        k1B[pullgrpid-1] = atof(tmpstr[2].c_str());
        ck1B |= (bitunit << pullgrpid-1);
	//cout << "# " << pullgrpid-1 << " k1B = " << k1B[pullgrpid-1] << endl;
      } else {
        const uint pullgrpid = getpullgrpid("pull_umb_k1",tmpstr[0]);
        k1A[pullgrpid-1] = atof(tmpstr[2].c_str());
        ck1A |= (bitunit << pullgrpid-1);
	//cout << "# " << pullgrpid-1 << " k1A = " << k1A[pullgrpid-1] << endl;
      }
    } else if(tmpstr[0].compare(0,string("pull_umb_r0").length(),"pull_umb_r0") == 0) {
      if(tmpstr[0].compare(0,string("pull_umb_r0B").length(),"pull_umb_r0B") == 0) {
        const uint pullgrpid = getpullgrpid("pull_umb_r0B",tmpstr[0]);
        r0B[pullgrpid-1] = atof(tmpstr[2].c_str());
        cr0B |= (bitunit << pullgrpid-1);
      } else {
        const uint pullgrpid = getpullgrpid("pull_umb_r0",tmpstr[0]);
        r0A[pullgrpid-1] = atof(tmpstr[2].c_str());
        cr0A |= (bitunit << pullgrpid-1);
      }
    } else if(tmpstr[0].compare(0,string("pull_umb_r1").length(),"pull_umb_r1") == 0) {
      if(tmpstr[0].compare(0,string("pull_umb_r1B").length(),"pull_umb_r1B") == 0) {
        const uint pullgrpid = getpullgrpid("pull_umb_r1B",tmpstr[0]);
        r1B[pullgrpid-1] = atof(tmpstr[2].c_str());
        cr1B |= (bitunit << pullgrpid-1);
      } else {
        const uint pullgrpid = getpullgrpid("pull_umb_r1",tmpstr[0]);
        r1A[pullgrpid-1] = atof(tmpstr[2].c_str());
        cr1A |= (bitunit << pullgrpid-1);
      }
    } else if(tmpstr[0].compare(0,string("pull_init").length(),"pull_init") == 0) {
      if(tmpstr[0].compare(0,string("pull_initB").length(),"pull_initB") == 0) {
        const uint pullgrpid = getpullgrpid("pull_initB",tmpstr[0]);
        initB[pullgrpid-1] = atof(tmpstr[2].c_str());
        cinitB |= (bitunit << pullgrpid-1);
      } else {
        const uint pullgrpid = getpullgrpid("pull_init",tmpstr[0]);
        initA[pullgrpid-1] = atof(tmpstr[2].c_str());
        cinitA |= (bitunit << pullgrpid-1);
      }
    }
  }
  cout << endl << endl;
  //Check if we read the correct number of parameters
  if(ck0A.count() != npullgrps || ck1A.count() != npullgrps || cinitA.count() != npullgrps || cr0A.count() != npullgrps || cr1A.count() != npullgrps) {
    cerr << "Number of state-A parameters is less than the number pull-groups " << npullgrps << endl;
    cout << "number of k0: " << ck0A.count() << endl;
    cout << "number of k1: " << ck1A.count() << endl;
    cout << "number of init: " << cinitA.count() << endl;
    cout << "number of r0: " << cr0A.count() << endl;
    cout << "number of r1: " << cr1A.count() << endl;
    exit(-1);
  }   
  //Construct the functor for each pull-group
  for(uint i = 0; i < npullgrps; ++i) {
    if(k0A[i] == k0B[i] && k0A[i] == 0 && k1A[i] == k1B[i] && k1A[i] ==0 && (rcmask & (bitunit << i)).none()) { 
      //printf("%d %f %f %f %f %f %f %f %f %f %f %f skipped\n",i,k0A[i],k0B[i],k1A[i],k1B[i],initA[i],initB[i],r0A[i],r0B[i],r1A[i],r1B[i],lambdas[current_lambda]);
      continue; 
    }
    //Occasionaly, the B-state parameters are not available and we'll assume they're the same as the A-state ones
    if((ck0A & (bitunit << i)).any() && (ck0B & (bitunit << i)).none()) { k0B[i] = k0A[i]; } 
    if((ck1A & (bitunit << i)).any() && (ck1B & (bitunit << i)).none()) { k1B[i] = k1A[i]; } 
    if((cinitA & (bitunit << i)).any() && (cinitB & (bitunit << i)).none()) { initB[i] = initA[i]; } 
    if((cr0A & (bitunit << i)).any() && (cr0B & (bitunit << i)).none()) { r0B[i] = r0A[i]; } 
    if((cr1A & (bitunit << i)).any() && (cr1B & (bitunit << i)).none()) { r1B[i] = r1A[i]; } 
    //printf("%d %f %f %f %f %f %f %f %f %f %f %f\n",i,k0A[i],k0B[i],k1A[i],k1B[i],initA[i],initB[i],r0A[i],r0B[i],r1A[i],r1B[i],lambdas[current_lambda]);
    valtype l;
    if(!lambdas.size()) { l = 0;}
    else { l = lambdas[current_lambda]; }
    //funct.insert(pair<uint,umbrella_fb>(i,umbrella_fb(k0A[i],k0B[i],k1A[i],k1B[i],initA[i],initB[i],r0A[i],r0B[i],r1A[i],r1B[i],l)));
    funct[i].push_back(new umbrella_fb(k0A[i],k0B[i],k1A[i],k1B[i],initA[i],initB[i],r0A[i],r0B[i],r1A[i],r1B[i],l));
  }
  //Construct Hamiltonian
  vector<valtype> k0,k1,init,r0,r1;
  vector<uint> RC_pullgrpids;
  for(uint i = 0; i < npullgrps; ++i) {
    //For RCs that are interesting, we build RST hamiltonian on them
    //NOTE that even for RCs that are virtually NOT restrained, we still build RST 
    //hamiltonian on them but the parameters are implicitly zero
    if( (rcmask & (bitunit << i)).any() ) {
      if(funct.find(i) == funct.end()) {
	cerr << "Pull group " << i << " is interesting but there's no umbrella restraint on it.\n";
	exit(-1);
      }
      k0.push_back(funct[i][funct[i].size()-1]->getk0());
      k1.push_back(funct[i][funct[i].size()-1]->getk1());
      init.push_back(funct[i][funct[i].size()-1]->getinit());
      r0.push_back(funct[i][funct[i].size()-1]->getr0());
      r1.push_back(funct[i][funct[i].size()-1]->getr1());
      RC_pullgrpids.push_back(i);
    }
  }
  //RC_pullgrpids.size() is how many RC we're interested in, which will be the first RC_pullgrpids.size() elements of the histogram data point
  //funct.begin()->size() tells how man states we've processed, including the current one
  //so the potential due to the restrained RCs we're not interested will be combined as the whichlambda'th element in one histogram data point
  if(!funct.size()) { cerr << "There's no any restraint applied on any of the RCs you're interested in. Check the input files again\n"; exit(-1); }
  //dumRCsize is how many restrained RC we're NOT interested in
  const uint dumRCsize = funct.size() - RC_pullgrpids.size();
  const uint whichlambda = dumRCsize ? (RC_pullgrpids.size() + (funct.begin()->second).size() - 1) : 0;
  vector<valtype> params;
  params.push_back(BoltzmannkJ);
  params.push_back(T);
  params.insert(params.end(),k0.begin(),k0.end());
  params.insert(params.end(),k1.begin(),k1.end());
  params.insert(params.end(),init.begin(),init.end());
  params.insert(params.end(),r0.begin(),r0.end());
  params.insert(params.end(),r1.begin(),r1.end());
  params.push_back(dumRCsize ? 1.0 : 0.0);  //should just weight it to 0 if dumRCsize == 0
  params.push_back(whichlambda); 
  V.push_back(new Hamiltonian<RST_fbXLAMBDAsgl>(params));
  cout << "#mdp2pullpot processed " << RC_pullgrpids.size() << " RCs you're interested in, whose id (indexed begins from zero) are: ";
  copy(RC_pullgrpids.begin(),RC_pullgrpids.end(),ostream_iterator<uint>(cout," "));
  cout << " and they'll be binned into histogram.\n";
  cout << "#The rest " << dumRCsize << " restrained RCs will be combined into their total bias energy, which will be binned into histogram\n";
  cout << "#We're now at the " << (funct.begin()->second).size() << "'th hamiltonian\n";
  return nlines;
}
