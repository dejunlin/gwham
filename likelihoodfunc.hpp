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
#include "exception.hpp"
using namespace std;

//! merge the trees by recursively base on the nearest neighbor
//NOTE that the earlier trees in the array will
//absorb those come latter and destroy their topology, i.e.,
//the merged tree will persist consistency with the earlier ones
template < class FTree >
void mergebranchNNB(vector<FTree>& trees) {
  if(trees.size() <= 1) { return; }
  const uint Ntreesb = trees.size();
  for(uint i = 0; i < trees.size(); ++i) {
    auto& root = trees[i];
    for(uint j = i+1; j < trees.size(); ) {
      auto& branch = trees[j];
      auto& rroutes = root.getroutes(); 
      auto& broutes = branch.getroutes();
      bool ifmerged = false;
      for(auto itr = rroutes.begin(); itr != rroutes.end(); ++itr) {
        const auto rnode = itr->first;
        auto itb = broutes.find(rnode);
        if(itb != broutes.end()) {
          //loop through all edges in branch and add them to root
	  //TODO: test if we need to first connect the joint node
          const auto& bedges = branch.getedges();
          for(const auto& bedge : bedges) {
            root.addedge(bedge);
          }
          trees.erase(trees.begin()+j);
          ifmerged = true;
          //we only establish one connection
          break;
        }
      }
      if(!ifmerged) { ++j; }
    }
  }
  //recurse until all possible branches are merged
  if(trees.size() > 1 && trees.size() != Ntreesb) {
    mergebranchNNB(trees);
  }

  //reaching here means all branches are merged
}

//! merge the trees based on 2nd connection
template < class FTree >
void mergebranch2nd(const vector<vector<uint>>& nbnodes, vector<FTree>& trees, vector<typename FTree::Data>& f) {
  if(trees.size() <= 1) { return; }
  const auto headnode = f.begin();
  const uint Ntreeb = trees.size();
  for(uint i = 0; i < trees.size(); ++i) {
    auto& root = trees[i];
    for(uint j = i+1; j < trees.size(); ) {
      auto& branch = trees[j];
      //for each node in the root, lookup the 2nd neighbors and see 
      //if any of them can be found in the branch; if so, merge them
      auto& rroutes = root.getroutes(); 
      auto& broutes = branch.getroutes();
      bool ifmerged = false;
      for(auto itr = rroutes.begin(); itr != rroutes.end(); ++itr) {
        const auto& rnblist = nbnodes[itr->first-headnode];
	for(const auto& id : rnblist) {
	  const auto node = headnode + id;
	  auto itb = broutes.find(node);
	  if(itb != broutes.end()) {
	    const typename FTree::Edge newedge = {itr->first ,node};
	    //recursively add the branch to the root
	    root.append(branch, newedge);
	    ifmerged = true;
	    break;
	  }
	}
	if(ifmerged) { break; }
      }
      if(ifmerged) { 
	trees.erase(trees.begin()+j);
	continue;
      }
      ++j;
    }
  }
  if(trees.size() > 1 && trees.size() != Ntreeb) {
    mergebranch2nd(nbnodes, trees, f);
  }
}

//! merge the trees
//NOTE that the earlier trees in the array will
//absorb those come latter and destroy their topology, i.e.,
//the merged tree will persist consistency with the earlier ones
template < class FTree >
void mergebranch(const vector<vector<uint>>& nbnodes, vector<FTree>& trees, vector<typename FTree::Data>& f) {
  //first merge trees using nearest neighbors
  mergebranchNNB(trees);
  if(trees.size() <= 1) { return; }
  //then connect the primary trees
  mergebranch2nd(nbnodes, trees, f);
}

//! build a ftree based on the neighbor list for each node
//the first neighbor of each node is the nearest neighbor
//and we'll first build N trees using the nearest neighbors 
//of each node and merge the tree recursively
template < class FTree >
vector<FTree> buildftree(const vector<vector<uint>>& nbnodes, vector<typename FTree::Data>& f) {
  if(nbnodes.size() != f.size()) {
    throw(General_Exception("input size of list of neighboring node is not the same as the size of f array"));
  }
  typedef typename FTree::Node Node;
  const uint Ntree = nbnodes.size();
  vector<FTree> ans(Ntree);
  const Node headnode = f.begin();
  for(uint i = 0; i < Ntree; ++i) {
    const Node node1 = headnode + i;
    //first neighbor of each state is the nearest neighbor
    const Node node2 = headnode + nbnodes[i][0];
    //assume node1 is closer to the headnode
    ans[i].addedge({node1, node2});
  }
  //! eliminate duplicate pair -- this need to be done only once
  for(uint i = 0; i < ans.size(); ++i) {
    const auto treei = ans[i];
    for(uint j = i+1; j < ans.size();) {
      if(treei.nodeset_equal(ans[j])) {
	ans.erase(ans.begin() + j);
	continue;
      }
      ++j;
    }
  }
  mergebranch(nbnodes, ans, f);

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
    typedef array<Node, 2> Edge;
  private:
    typedef ftree<V> ThisType;
    typedef vector<Edge> Edges;
  friend vector<ThisType> buildftree<ThisType>(const vector<vector<uint>>&, vector<V>&);
  public:
    //! Construct the tree 
    ftree() {};

    //! The indices of all the edges touching a node
    map<Node, vector<uint> > joints;

    //! The indices of edges that connect each node to the headnode 
    map<Node, vector<uint> > routes;

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
    vector<Nodes> ports;

    //! if the 2 trees have the same set of nodes
    bool nodeset_equal(const ThisType& rhs) const {
      bool ans = true;
      const auto& rhsroutes = rhs.getroutes();
      for(auto it = routes.begin(); it!= routes.end(); ++it) {
	const auto& node = it->first;
	if(rhsroutes.find(node) == rhsroutes.end()) {
	  ans = false;
	  break;
	}
      }
      return ans;
    }

    //! add a joint
    inline void addjoint(const Node& node, const uint& edgeid) { 
      auto it = joints.find(node);
      if(it != joints.end()) { it->second.push_back(edgeid); }
      else { joints[node] = {edgeid}; }
    }

    //! add a node
    inline void addnode(const Node& node) { nodes.emplace_back(node); }
    //! add an port
    inline void addport(const uint& hubid, const Nodes& newport) { 
      ports.resize(edges.size());
      auto& port = ports[hubid];
      port.insert(port.end(), newport.begin(), newport.end());
    }
    //! add a route
    inline void addroute(const Node& node, const uint& edgeid) {
      auto it = routes.find(node);
      if(it != routes.end()) {
	it->second.push_back(edgeid);
      } else {
        routes[node] = {edgeid};
      }
      addport(edgeid, {node});
      //also update this node's route to former edges
      const auto& edge = edges[edgeid];
      const auto& itfnode = routes.find(edge[0]);
      //if such former edge exist and the joint node is not the headnode
      if(itfnode != routes.end() && itfnode->second.size()) {
	//we assume that in the array of edges each node map to, the smaller
	//the corresponding index is, the closer such edge is to the node
	addroute(node, (itfnode->second)[0]);
      }
    }
    //! add an edge
    /** NOTE this will re-establish the polarity of the edge
     */
    inline void addedge(const Edge& edge) {
      const auto head = edge[0];
      const auto tail = edge[1];
      const auto ithead = routes.find(head);
      const auto ittail = routes.find(tail);
      //if both nodes in this edge exist in the routes, just skip
      if(ithead != routes.end() && ittail != routes.end()) {
	return;
	//otherwise, whichever end exists is assumed closer to the headnode
      } else if(ithead != routes.end()) {
	//head exists but tail doesn't
        addnode(tail);
        edges.emplace_back(edge);
        addroute(tail, edges.size()-1);
      } else if(ittail != routes.end()) {
	//tail exists but head doesn't
        addnode(head);
	//re-establish polarity of this edge
	Edge newedge = {tail, head};
        edges.emplace_back(newedge);
        addroute(head, edges.size()-1);
      } else {
	//if neither exists, this has to been a new tree otherwise we bail
	if(!this->empty()) {
	  throw(General_Exception("Can't add an isolated edge to non-empty ftree"));
	}
	//assume the first node in this head become head node
        addnode(head);
        addnode(tail);
        edges.emplace_back(edge);
	routes[head] = vector<uint>{};
	addroute(tail, edges.size()-1);
      }
      //map both nodes to this new edge
      addjoint(head, edges.size()-1);
      addjoint(tail, edges.size()-1);
    }

    //! append another tree to this one via the input joint node
    inline ThisType& append(const ThisType& src, const Edge& junction) {
      const auto& node1 = junction[0];
      const auto& node2 = junction[1];
      auto it1 = routes.find(node1);
      auto it2 = routes.find(node2);
      //if both nodes belong to this tree already, just return;
      if(it1 != routes.end() && it2 != routes.end()) {
	return *this;
      } else if(!this->empty() &&  it1 == routes.end() && it2 == routes.end()) {
	throw(General_Exception("Can't append to a non-empty ftree without junction"));
      }
      //NOTE: if this is an tempty tree, doesn't matter which one becomes head
      const Node newnode = (it1 == routes.end()) ?  node1 : node2;

      addedge(junction);

      //recursively add all the source edges 
      const auto& srcjoints = src.getjoints();
      const auto itsrcjoint = srcjoints.find(newnode); 
      if(itsrcjoint == srcjoints.end()) {
	throw(General_Exception("Junction node doesn't belong to the ftree to be appended"));
      }
      const auto& srcedges = src.getedges();
      for(const auto& srcedgeid : itsrcjoint->second) {
	this->append(src, srcedges[srcedgeid]); 
      }
      return *this;
    }
    
    //! if this is an empty tree
    inline bool empty() const { return !(nodes.size() || edges.size() || routes.size() || ports.size()); } 
    //! number of nodes
    inline uint nodesize() const { return nodes.size(); }
    //! number of edges
    inline uint edgesize() const { return edges.size(); }
    //! return tne node list
    inline const Nodes& getnodes() const { return nodes; } 
    //! return tne edge list
    inline const Edges& getedges() const { return edges; } 
    //! return tne port list
    inline const vector<Nodes>& getports() const { return ports; } 
    //! return tne joint list
    inline const map<Node, vector<uint> >& getjoints() const { return joints; } 
    //! return tne route list
    inline map<Node, vector<uint> >& getroutes() { return routes; } 
    inline const map<Node, vector<uint> >& getroutes() const { return routes; } 
};

template < class DOS, class NARRAY, class FTree >
class LikeliHoodFunc {
  typedef typename NARRAY::iterator narrit;
  typedef typename NARRAY::const_iterator narrcit;
  typedef typename DOS::iterator dosit;
  typedef typename DOS::const_iterator doscit;
  typedef decltype( declval<narrit>()->second ) valtype; 
  typedef unsigned int uint;
  typedef unsigned long ulong;
  public:
    LikeliHoodFunc(DOS& _dos, 
                   const vector<ulong>& _N, 
		   const NARRAY& _C, 
		   const vector<NARRAY>& _NgexpmH, 
		   vector<valtype>& _f, 
		   vector<valtype>& _expf, 
		   FTree& _tree
		  ) 
      : dos(_dos), N(_N), C(_C), NgexpmH(_NgexpmH), 
        x(_N.size()-1, numeric_limits<valtype>::max()), fret(0), gradf(x.size(), 0), graddf(x.size(), 0), 
	f(_f), expf(_expf), tree(_tree) 
    {
      if(N.size() != NgexpmH.size() || N.size() != f.size() || N.size() != expf.size()) {
	throw(General_Exception("Size of array N, NgexpmH, f and expf are not the same"));
      }
      if(tree.nodesize() != f.size()) {
        throw(General_Exception("Number of nodes in the tree isn't the same \
        as the number of states. Check the implementation of buildtree function"));
      }
      if(tree.nodesize() != tree.edgesize()+1) {
        throw(General_Exception("Number of nodes in the tree isn't the same \
        as the number of edges + 1. Check the implementation of buildtree function"));
      }
      const auto& edges = tree.getedges();
      const auto& headnode = f.begin();
      cout << "#Edges of tree: \n";
      for(const auto& edge : edges) {
	const auto& nodei = edge[0];
	const auto& nodej = edge[1];
	cout << "# " << nodei-headnode << "----" << nodej-headnode << endl;
      }

    };
    
    //calculate both the value and gradient of the likelihood
    //function
    valtype operator() (const vector<valtype>& deltaf) {
      //if this operator or this->df() has been called using exactly the 
      //same argument, we just retrieve results from last 
      //call without any update; otherwise, we need to update 
      //dos and x
      const bool update = (deltaf != x);
      
      fret = 0;
      narrcit itC = C.begin();
      if(update) {
        x = deltaf;
        //here we constraint f[0] to 0.0
        f[0] = 0;
        expf[0] = 1;

	const auto& edges = tree.getedges();
	const auto headnode = f.begin();
        for(uint i = 0; i < deltaf.size(); ++i) {
	  const auto& edge = edges[i];
	  //iterators to the f array
	  const auto& itfi = edge[1];
	  const auto& itfj = edge[0];
	  *itfi = *itfj + deltaf[i];
	  //id is the index in the f array
	  const uint id = itfi-headnode;
	  const auto& fi = *itfi;
  	  const auto e = exp(fi);
  	  expf[id] = e;
  	  fret -= N[id]*fi;
          //since f[0] is constrained to 0, gradf[0] is actually 
          //gradient with respect to f[1]; Here we need to make 
	  //sure gradf[i] corresponds to f[i-1]
          gradf[id-1] = -valtype(N[id]);
        }
        vector<narrcit> itsNgexpmH;
        for_each(NgexpmH.begin(), NgexpmH.end(), [&itsNgexpmH](const NARRAY& narr) { itsNgexpmH.emplace_back(narr.begin()); });
        for(dosit itdos = dos.begin(); itdos != dos.end(); ++itdos) {
          vector<narrcit> itsNgexpmH_buff{itsNgexpmH};
          const auto& c = itC->second;
          const auto d = dos(itC, itsNgexpmH, expf);
          itdos->second = d;
          fret -= c*log(d);
          for(uint i = 0; i < gradf.size(); ++i) {
            gradf[i] += d*itsNgexpmH_buff[i+1]->second*expf[i+1];
            ++itsNgexpmH_buff[i+1];
          }
        }
        graddf.assign(gradf.size(), 0.0);
        const auto& ports = tree.getports();
        const auto& nodes = tree.getnodes();
        for(uint i = 0; i < graddf.size(); ++i) {
          auto& gdf = graddf[i];
          const auto& port = ports[i];
          for(auto& node : port) { 
            const auto id = node - headnode; 
            gdf += gradf[id-1];
          }
        }
      } 
      return fret;
    }
    
    void df (const vector<valtype>& deltaf, vector<valtype>& _graddf) {
      //if this operator or this->df() has been called using exactly the 
      //same argument, we just retrieve results from last 
      //call without any update; otherwise, we need to update 
      //dos and x
      const bool update = (deltaf != x);
      if(update) {
	this->operator()(deltaf);
      }
      _graddf = graddf;
    }
     
  private:
    DOS& dos;
    const vector<ulong>& N;
    const NARRAY& C;
    const vector<NARRAY>& NgexpmH;
    vector<valtype> x;
    valtype fret;
    vector<valtype> gradf;
    vector<valtype> graddf;
    vector<valtype>& f;
    vector<valtype>& expf;
    FTree& tree;
};


#endif
