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
  vector<uint> N(nwin,0);
  for(uint i = 0; i < nwin; ++i) {
    char winid[MAXNDIGWIN];
    sprintf(winid,"%d",i);
    string fname = winid;
    N[i] += reader(fname,hists[i],subtrjs[i]);
    hists[i].print(valtype(N[i]));
  }

  //Calculate variance of the mean at each bin from each traj
  vector<map<histogram::hist_coord, valtype> > variances(nwin);
  map<histogram::hist_coord, vector<uint> > record;
  //cout << "#coord hist mean sqrmean variance trjsize\n";
  for(uint i = 0; i < nwin; ++i) {
    map<histogram::hist_coord, valtype>& variance = variances[i];
    map<histogram::hist_coord, vector<valtype> >& subtrj = subtrjs[i];
    map<histogram::hist_coord, vector<valtype> >::const_iterator it;
    for(it = subtrj.begin(); it != subtrj.end(); ++it) {
      const histogram::hist_coord coord = it->first;
      const vector<valtype>& trj = it->second;
      //first calculate the mean
      valtype mean = 0.0;
      valtype sqrmean = 0.0;
      vector<valtype>::const_iterator trjit;
      for(trjit = trj.begin(); trjit != trj.end(); ++trjit) {
	const valtype x = *trjit;
	mean += x;
	sqrmean += x*x;
      }
      mean /= N[i]; //not trj.size() here since we're estimating the histogram function
      sqrmean /= N[i];
      //then calculate the variance of the mean
      //NOTE: this is actually just variance of the sample not of the mean
      //since we'll have to multiple by N[i]**2 later anyway, we don't 
      //divide by N[i] but instead multiply by N[i] (instead of N[i]**2) later 
      variance[coord] = (sqrmean - mean*mean);

      /*copy(coord.begin(),coord.end(),ostream_iterator<uint>(cout, " "));
      cout << i << " " << mean << " " << sqrmean << " " << variance[coord] << " " << trj.size() << endl;*/

    }

    //build records of what bins are contributed from which histograms btw
    histogram::const_iterator histit;
    for(histit = hists[i].begin(); histit != hists[i].end(); ++histit) {
      const histogram::hist_coord coord = histit->first;
      if(record.find(coord) != record.end()) { record[coord].push_back(i); }
      else { record[coord] = vector<uint>(1,i); }
    }
  }

  //Weight each histogram using the inverse of the variance at each bin
  narray rho(ndim,nbins,hv,lv);
  map<histogram::hist_coord,vector<uint> >::const_iterator it;
  for(it = record.begin(); it != record.end(); ++it) {
    const histogram::hist_coord coord = it->first;
    const vector<uint>& histids = it->second;
    vector<uint>::const_iterator hit;
    valtype num = 0.0;
    valtype denom = 0.0;
    for(hit = histids.begin(); hit != histids.end(); ++hit) {
      const uint histid = *hit;
      const valtype weight = 1/(variances[histid][coord]*N[histid]);
      denom += weight;
      num += hists[histid][coord]*weight;
    }
    rho[coord] =num/denom;
  }
  //Normalize the probability (just for comparision with other programs)
  printf("#%10s%30s%30s%30s\n","Bin","Vals","PMF","RhoNormalized");
  valtype sum = 0.0;
  narray::iterator itmax = rho.begin();
  for(narray::iterator it = rho.begin(); it != rho.end(); ++it) {
    if(it->second > itmax->second) { itmax = it;}
    sum += it->second;
  }
  for(narray::iterator it = rho.begin(); it != rho.end(); ++it) {
    const coordtype bin = it->first;
    const vector<valtype> val = rho.coord2val(bin);
    const valtype pmf = -Boltzmannkcal*T*log(it->second/itmax->second);
    const valtype rhonorm = it->second/sum;
    for(uint i = 0; i < bin.size(); ++i) { printf("%10d",bin[i]);}
    for(uint i = 0; i < val.size(); ++i) { printf("%30.15lf",val[i]);}
    printf("%30.15lf%30.15lf\n",pmf,rhonorm);
  }
}
