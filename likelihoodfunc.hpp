#ifndef LIKELIHOODFUNC_HPP
#define LIKELIHOODFUNC_HPP
/*
 * =====================================================================================
 *
 *       Filename:  likelyhoodfunc.hpp
 *
 *    Description:  The WHAM likelyhood functor
 *
 *        Version:  1.0
 *        Created:  22/09/14 11:16:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <cmath>
#include <vector>
#include <iterator>
#include <array>
#include <list>
#include "exception.hpp"
using namespace std;

//! Build the branch of a ftree, return a vector of nodes in this branch
//including the branch node
template < class FTree >
vector<uint> buildbranch(FTree& tree, map<uint,bool>& seen, const vector<vector<uint>>& nbnodes, const uint& i, vector<typename FTree::Data>& f) {
  if(seen.find(i) != seen.end()) { return vector<uint>{}; }
  //we always start from f[0]
  const typename FTree::Node headnode = f.begin();
  tree.addnode(headnode+i);
  vector<uint> branch{tree.nodesize()-1}; 
  seen[i] = true;
  for(const auto& j : nbnodes[i]) {
    //return if reach a dead end
    if(j == i) return branch;
    //skip if has been seen before
    if(seen.find(j) != seen.end()) { 
      continue;
    }
    //coming all the way here means a valid edge
    tree.addedge({headnode+i, headnode+j});
    const auto hubid = tree.edgesize()-1;
    //recurse over all the branches
    auto subbranch = buildbranch(tree, seen, nbnodes, j, f);
    branch.insert(branch.end(), subbranch.begin(), subbranch.end());
    tree.addport(hubid, subbranch);
  }
  return branch;
}

//! build a ftree based on the neighbor list for each node
template < class FTree >
vector<FTree> buildftree(const vector<vector<uint>>& nbnodes, vector<typename FTree::Data>& f) {
  if(nbnodes.size() != f.size()) {
    throw(General_Exception("input size of list of neighboring node is not the same as the size of f array"));
  }
  vector<FTree> ans;
  FTree tree;
  map<uint, bool> seen;
  for(uint i = 0; i < nbnodes.size(); ++i) {
    if(seen.find(i) != seen.end()) { continue; }
    buildbranch(tree, seen, nbnodes, i, f);
    if(tree.nodesize()) { 
      ans.emplace_back(tree); 
      tree = FTree{};
    } 
  }
  return ans;
}

//! This class represents the tree structure of f (free energy)
//and deltaf (difference between two f) with each node being 
//one unique f and each edge being one unique deltaf
template < class V >
class ftree {
  public:
    typedef V Data;
    typedef typename vector<V>::iterator Node;
    typedef vector<Node> Nodes;
  private:
    typedef ftree<V> ThisType;
    typedef array<Node, 2> Edge;
    typedef vector<Edge> Edges;
  friend vector<ThisType> buildftree<ThisType>(const vector<vector<uint>>&, vector<V>&);
  friend vector<uint> buildbranch<ThisType>(ThisType&, map<uint,bool>&, const vector<vector<uint>>&, const uint&, vector<V>&);
  public:
    //! Construct the tree 
    ftree() {};

    //! All the nodes in this tree
    Nodes nodes;

    //! All the edges
    /** NOTE the 1st element of each Edge
     * is the minuend and the 2nd one is the 
     * subtrahend
     */
    Edges edges;

    //! For each edge, the indices of the
    //nodes that depend on it to connect to 
    //the tip of the tree
    vector<vector<uint> > ports;

    //! add a node
    inline void addnode(const Node& node) { nodes.emplace_back(node); }
    //! add an edge
    inline void addedge(const Edge& edge) { edges.emplace_back(edge); }
    //! add an port
    inline void addport(const uint& hubid, const vector<uint>& port) { 
      ports.resize(edges.size());
      ports[hubid] = port;
    }
    //! number of nodes
    inline uint nodesize() const { return nodes.size(); }
    //! number of edges
    inline uint edgesize() const { return edges.size(); }
    //! return tne node list
    inline const Nodes& getnodes() const { return nodes; } 
    //! return tne edge list
    inline const Edges& getedges() const { return edges; } 
    //! return tne port list
    inline const vector<vector<uint>>& getports() const { return ports; } 
};

template < class DOS, class NARRAY >
class LikeliHoodFunc {
  typedef typename NARRAY::iterator narrit;
  typedef typename NARRAY::const_iterator narrcit;
  typedef typename DOS::iterator dosit;
  typedef typename DOS::const_iterator doscit;
  typedef decltype( declval<narrit>()->second ) valtype; 
  typedef unsigned int uint;
  typedef unsigned long ulong;
  public:
    LikeliHoodFunc(DOS& _dos, const vector<ulong>& _N, const NARRAY& _C, const vector<NARRAY>& _NgexpmH, vector<valtype>& _f, vector<valtype>& _expf) 
      : dos(_dos), N(_N), C(_C), NgexpmH(_NgexpmH), x(_N.size()-1, numeric_limits<valtype>::max()), f(_f), expf(_expf) 
    {
      if(N.size() != NgexpmH.size() || N.size() != f.size() || N.size() != expf.size()) {
	throw(General_Exception("Size of array N, NgexpmH, f and expf are not the same"));
      }
    };

    valtype operator() (const vector<valtype>& deltaf) {
      //if this operator or this->df() has been called using exactly the 
      //same argument, we just retrieve results from last 
      //call without any update; otherwise, we need to update 
      //dos and x
      const bool update = (deltaf != x);
      //here we constraint f[0] to 0.0
      f[0] = 0;
      expf[0] = 1;
      
      valtype fret = 0;
      narrcit itC = C.begin();
      if(update) {
        x = deltaf;
        for(uint i = 0; i < deltaf.size(); ++i) {
  	  f[i+1] = f[i] + deltaf[i];
  	  const auto e = exp(f[i+1]);
  	  expf[i+1] = e;
  	  fret -= N[i+1]*f[i+1];
        }
        vector<narrcit> itsNgexpmH;
        for_each(NgexpmH.begin(), NgexpmH.end(), [&itsNgexpmH](const NARRAY& narr) { itsNgexpmH.emplace_back(narr.begin()); });
        for(dosit itdos = dos.begin(); itdos != dos.end(); ++itdos) {
          const auto& c = itC->second;
          const auto d = dos(itC, itsNgexpmH, expf);
          itdos->second = d;
          fret -= c*log(d);
        }
      } else {
        for(uint i = 0; i < deltaf.size(); ++i) {
  	  fret -= N[i+1]*f[i+1];
        }
        for(dosit itdos = dos.begin(); itdos != dos.end(); ++itdos ) {
  	  const auto& c = itC->second;
	  const auto& d = itdos->second;
	  ++itC;
	  fret -= c*log(d);
	}
      }
      return fret;
    }
    
    void df (const vector<valtype>& deltaf, vector<valtype>& graddf) {
      //if this operator or this->df() has been called using exactly the 
      //same argument, we just retrieve results from last 
      //call without any update; otherwise, we need to update 
      //dos and x
      const bool update = (deltaf != x);
      vector<valtype> gradf(deltaf.size(), 0.0);
      narrcit itC = C.begin();
      if(update) {
	x = deltaf;
        //here we constrain f[0] to 0.0
        f[0] = 0;
        expf[0] = 1.0;

        for(uint i = 0; i < deltaf.size(); ++i) {
          f[i+1] = f[i] + deltaf[i];
          const auto e = exp(f[i+1]);
          expf[i+1] = e;
          //since f[0] is constrained to 0, gradf[0] is actually 
          //gradient with respect to f[1]
          gradf[i] -= N[i+1];
        }
        vector<narrcit> itsNgexpmH;
        for_each(NgexpmH.begin(), NgexpmH.end(), [&itsNgexpmH](const NARRAY& narr) { itsNgexpmH.emplace_back(narr.begin()); });
        for(dosit itdos = dos.begin(); itdos != dos.end(); ++itdos ) {
          vector<narrcit> itsNgexpmH_buff{itsNgexpmH};
          const auto d = dos(itC, itsNgexpmH, expf);
	  itdos->second = d;
          for(uint i = 0; i < gradf.size(); ++i) {
            gradf[i] += d*itsNgexpmH_buff[i+1]->second*expf[i+1];
            ++itsNgexpmH_buff[i+1];
          }
        }
      } else {
        for(uint i = 0; i < deltaf.size(); ++i) {
          gradf[i] -= N[i+1];
	}
        vector<narrcit> itsNgexpmH_buff;
        for_each(NgexpmH.begin(), NgexpmH.end(), [&itsNgexpmH_buff](const NARRAY& narr) { itsNgexpmH_buff.emplace_back(narr.begin()); });
        for(dosit itdos = dos.begin(); itdos != dos.end(); ++itdos ) {
          const auto d = itdos->second;
          for(uint i = 0; i < deltaf.size(); ++i) {
            gradf[i] += d*itsNgexpmH_buff[i+1]->second*expf[i+1];
            ++itsNgexpmH_buff[i+1];
          }
        }
      }
      graddf.assign(gradf.size(), 0.0);
      for(uint i = 0; i < gradf.size(); ++i) {
	const auto& gf = gradf[i];
	for(uint j = 0; j <= i; ++j) {
	  graddf[j] += gf;
	}
      }
    }
     
  private:
    DOS& dos;
    const vector<ulong>& N;
    const NARRAY& C;
    const vector<NARRAY>& NgexpmH;
    vector<valtype> x;
    vector<valtype>& f;
    vector<valtype>& expf;
};


#endif
