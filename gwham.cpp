#include "typedefs.hpp"
#include "fileio_utils.hpp"
#include "ensemble.hpp"
#include "fileio.hpp"
#include "gnarray.hpp"
#include "trjsubtrj.hpp"
#include <vector>

using namespace std;

typedef gnarray<coordtype, valtype, valtype> narray;
typedef gnarray<coordtype, valtype, valtype> histogram;

int main(int argc, char* argv[]) {
  if (argc != 9) {
    cerr << "Usage: gwham nwin rcsample rcprint Temperature nbins hv lv stride\n";
    exit(-1);
  }
  uint k = 1;
  const uint nwin = atoi(argv[k++]); //# of replica
  const string rcsmpstr = string(argv[k++]); //which pull-group/RC you want to be histogramed on, where the pull-group id starts from 0 (1st pull-group not the reference group)
  const string rcprtstr = string(argv[k++]); //which pull-group/RC you want the PMF on. 
  					   //Since the mpd2rst and xxvg2hist classes will re-order the RCs so that 
					   //they are the 1st N dimension of the histogram (N is the number of RCs you histogramed on), 
					   //rcrpt must be in the range of [0,N-1]. We could potentially allow this range to extend over N-1 but
					   //that will require we know if all the restrained RCs are the RCs we're interested in
  const valtype T = atof(argv[k++]); //what temperature to use to calculate PMF
  const string nbinstr = string(argv[k++]); //# of bins in each dimension
  const string hvstr = string(argv[k++]); //upper bounds in each dimension
  const string lvstr = string(argv[k++]); //lower bounds in each dimension
  const uint stride = atoi(argv[k++]); //Only read every stride lines (excluding comment-lines) in the x.xvg files

  vector<uint> rcsmp, rcprt;
  parser<uint>(rcsmp,rcsmpstr);
  parser<uint>(rcprt,rcprtstr);
  vector<uint> nbins;
  vector<valtype> hv, lv;
  parser<uint>(nbins,nbinstr);
  parser<valtype>(hv,hvstr);
  parser<valtype>(lv,lvstr);

  //Check if the histogram parameters (nbins, hv, lv) make sense
  if(nbins.size() != hv.size() || hv.size() != lv.size()) {
    cerr << "The size of nbins, hv and lv don't match one another" << endl;
    exit(-1);
  }
  for(uint i = 0; i < hv.size(); ++i) {
    if(hv[i] <= lv[i]) {
      cerr << "The upper bound of dimension " << i << " is no larger than the lower bound\n";
      cerr << "hv = " << hv[i] << endl;
      cerr << "lv = " << lv[i] << endl;
      exit(-1);
    }
  }

  //Turn rc into a bitmask
  bitset<MAXNRST> rcsmpmask;
  for(uint i = 0; i < rcsmp.size(); ++i) {
    rcsmpmask |= (bitunit << rcsmp[i]);
  }

  //Read the files
  const NVT nvt(Boltzmannkcal,T);
  trjsubtrj<histogram> trj(stride,rcsmpmask,nvt);
  const fileio<trjsubtrj<histogram>, histogram, map<histogram::hist_coord, vector<valtype> > > reader(trj,std::fstream::in);
  const uint ndim = nbins.size();
  vector<histogram> hists(nwin,histogram(ndim,nbins,hv,lv));
  vector<map<histogram::hist_coord, vector<valtype> > > subtrjs(nwin);
  for(uint i = 0; i < nwin; ++i) {
    char winid[MAXNDIGWIN];
    sprintf(winid,"%d",i);
    string fname = winid;
    reader(fname,hists[i],subtrjs[i]);
  }
}
