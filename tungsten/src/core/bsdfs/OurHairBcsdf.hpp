#include "PrecomputedAzimuthalLobe.hpp"

#include "bsdfs/Bsdf.hpp"

#include "math/Angle.hpp"

#include <memory>

namespace Tungsten {

class OurHairBscdf
{
    
    const float alpha = 1.0f;

    // float _scaleAngleDeg;
    // float _melaninRatio;
    // float _melaninConcentration;
    // bool _overridesSigmaA;
    // Vec3f _sigmaA;
    // float _roughness;

    // float _scaleAngleRad;
    // std::unique_ptr<PrecomputedAzimuthalLobe> _nR, _nTT, _nTRT;
    // float _betaR, _betaTT, _betaTRT;
    // float _vR, _vTT, _vTRT;

    // static float I0(float x);
    // static float logI0(float x);

    // static float g(float beta, float theta);
    // static float D(float beta, float phi);

    // static float Phi(float gammaI, float gammaT, int p);

    // float ourM(float v, float sinThetaI, float sinThetaO, float cosThetaI, float cosThetaO) const;

    // float NrIntegrand(float beta, float wiDotWo, float phi, float h) const;
    // Vec3f NpIntegrand(float beta, float cosThetaD, float phi, int p, float h) const;

    // float sampleM(float v, float sinThetaI, float cosThetaI, float xi1, float xi2) const;

    // void precomputeAzimuthalDistributions();

public:
    static float csch (float theta);
    static float thetaConeR(float alpha, float thetaI);
    static float thetaConeTT(float alpha, float thetaI);
    static float thetaConeTRT(float alpha, float thetaI);
    static float OurM(float v, float sinThetaC, float sinThetaO, float cosThetaC, float cosThetaO, float IO);
};

}
