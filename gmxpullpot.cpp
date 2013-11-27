#include "gmxpullpot.hpp"
#include <iostream>

umbrella::umbrella(const valtype& _kA, const valtype& _kB, const valtype& _initA, const valtype& _initB, const valtype& _L) :
  kA(_kA),
  kB(_kB),
  initA(_initA),
  initB(_initB),
  L(_L),
  k((1-L)*kA + L*kB),
  init((1-L)*initA + L*initB)
  {
  }

umbrella::umbrella() :
  kA(0),
  kB(0),
  initA(0),
  initB(0),
  L(0),
  k(0),
  init(0)
  {}

valtype umbrella::operator() (const valtype& r) const {
  return 0.5*k*(r-init)*(r-init);
}

umbrella_fb::umbrella_fb (const valtype& _k0A, const valtype& _k0B, 
                          const valtype& _k1A, const valtype& _k1B,
                          const valtype& _initA, const valtype& _initB,
                          const valtype& _r0A, const valtype& _r0B,
                          const valtype& _r1A, const valtype& _r1B,
		          const valtype& _L
                         ) :
  k0A(_k0A),
  k0B(_k0B),
  k1A(_k1A),
  k1B(_k1B),
  initA(_initA),
  initB(_initB),
  r0A(_r0A),
  r0B(_r0B),
  r1A(_r1A),
  r1B(_r1B),
  L(_L),
  k0((1-L)*k0A + L*k0B),
  k1((1-L)*k1A + L*k1B),
  init((1-L)*initA + L*initB),
  r0((1-L)*r0A + L*r0B),
  r1((1-L)*r1A + L*r1B)
  {
  }

umbrella_fb::umbrella_fb () :
  k0A(0),
  k0B(0),
  k1A(0),
  k1B(0),
  initA(0),
  initB(0),
  r0A(0),
  r0B(0),
  r1A(0),
  r1B(0),
  L(0),
  k0(0),
  k1(0),
  init(0),
  r0(0),
  r1(0)
  {
  }

valtype umbrella_fb::operator() (const valtype& r) const {
  const valtype dev = r - init;
  if( dev >= r0 && dev <= r1 ) { return 0.0; }
  else if ( dev < r0 ) { return 0.5*k0*(dev-r0)*(dev-r0); }
  else { return 0.5*k1*(dev-r1)*(dev-r1); }
}
