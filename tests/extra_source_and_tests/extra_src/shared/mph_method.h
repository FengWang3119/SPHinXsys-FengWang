#ifndef MPH_METHOD_H
#define MPH_METHOD_H

#include "sphinxsys.h"

namespace SPH
{
int test_ccc = 1;

class KernelOfMPH : public Kernel
{
public:
    KernelOfMPH(Real h, Real particle_spacing);

    /**
     * Calculates the kernel value for
     * the given distance of two particles
     */
    virtual Real W_1D(const Real q) const override;
    virtual Real W_2D(const Real q) const override;
    virtual Real W_3D(const Real q) const override;

    virtual Real dW_1D(const Real q) const override;
    virtual Real dW_2D(const Real q) const override;
    virtual Real dW_3D(const Real q) const override;

    virtual Real d2W_1D(const Real q) const override;
    virtual Real d2W_2D(const Real q) const override;
    virtual Real d2W_3D(const Real q) const override;
};

}
#endif // MPH_METHOD_H