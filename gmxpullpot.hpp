#if !defined(GMXPULLPOTIO_HPP)
#define GMXPULLPOTIO_HPP

#include "typedefs.hpp"

//! Basic harmonic potential functor used in umbrella sampling
/** This only handles 1d restraint; TODO: we should extend this 
 * to handle all 3 dimension
 */
class umbrella {
  public:
    //! Cosntructor
    umbrella();
    //! Cosntructor
    umbrella(const valtype& _kA, const valtype& _kB, const valtype& _initA, const valtype& _initB, const valtype& _L);
    //! calculate potential given a distance
    valtype operator() (const valtype& r) const;
    valtype getk() const { return k; }
    valtype getinit() const { return init; }
  private:
    //! force constant in state A
    const valtype kA;
    //! force constant in state B
    const valtype kB;
    //! equilibrium position in state A
    const valtype initA;
    //! equilibrium position in state B
    const valtype initB;
    //! Thermodynamic state (Lambda parameters that switches from 0 to 1)
    const valtype L;
    //! force constant in current state
    const valtype k;
    //! equilibrium position in current state
    const valtype init;
};

//! Flat-bottom harmonic potential functor
/** This only handles 1d restraint; TODO: we should extend this 
 * to handle all 3 dimension
 */
class umbrella_fb {
  public:
    //! Constructor
    umbrella_fb(); 
    //! Constructor
    umbrella_fb(const valtype& _k0A, const valtype& _k0B, 
                const valtype& _k1A, const valtype& _k1B,
                const valtype& _initA, const valtype& _initB,
		const valtype& _r0A, const valtype& _r0B,
		const valtype& _r1A, const valtype& _r1B,
		const valtype& _L
               );
    //! calculate potential given a distance
    valtype operator() (const valtype& r) const;
  private:
    //! force constant of the left wing in state A
    const valtype k0A;
    //! force constant of the left wing in state B
    const valtype k0B;
    //! force constant of the right wing in state A
    const valtype k1A;
    //! force constant of the right wing in state B
    const valtype k1B;
    //! equilibrium position in state A
    const valtype initA;
    //! equilibrium position in state B
    const valtype initB;
    //! width of the left wing in state A
    const valtype r0A;
    //! width of the left wing in state B
    const valtype r0B;
    //! width of the right wing in state A
    const valtype r1A;
    //! width of the right wing in state B
    const valtype r1B;
    //! Thermodynamic state (Lambda parameters that switches from 0 to 1)
    const valtype L;

    //! force constant of the left wing in current state
    const valtype k0;
    //! force constant of the right wing in current state
    const valtype k1;
    //! equilibrium position in current state
    const valtype init;
    //! width of the left wing in current state
    const valtype r0;
    //! width of the right wing in current state
    const valtype r1;
};

#endif
