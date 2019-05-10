#include "HairBcsdf.hpp"

#include "sampling/PathSampleGenerator.hpp"

#include "bsdfs/Fresnel.hpp"

#include "math/GaussLegendre.hpp"

#include "io/JsonObject.hpp"

namespace Tungsten {

HairBcsdf::HairBcsdf()
: _scaleAngleDeg(2.0f),
  _melaninRatio(0.5f),
  _melaninConcentration(0.25f),
  _overridesSigmaA(false),
  _sigmaA(0.0f),
  _roughness(0.1f)
{
    _lobes = BsdfLobes(BsdfLobes::GlossyLobe | BsdfLobes::AnisotropicLobe);
}

float HairBcsdf::T(float mu, float h, float etaPrime) const {
    float gammaT = asin(h/etaPrime);
	return exp(-2 * mu * (1 + cos(2 * etaPrime)));
}

float HairBcsdf::A(float p, float h, float cosThetaT, float cosThetaD) const {
    float mu_a_prime = muA / cosThetaT;
    float f = Fresnel::dielectricReflectance(1.0f/Eta, cosThetaD*trigInverse(h));
    float etaPrime = std::sqrt(Eta*Eta - (1.0f - cosThetaD*cosThetaD))/cosThetaD;
    // DOUBLE CHECK REFLECTION CASE
    if (p == 1.0f)
        return (1.0f - f)*(1.0f - f)*T(mu_a_prime, h, etaPrime);
	return pow(1 - f, 2) * pow(f, p - 1) * pow(T(mu_a_prime, h, etaPrime), p);
}

float HairBcsdf::G(Vec2f U) const {
    return sqrt(-2.0f * log(U.x())) * cos(2 * PI * U.y());
}

float HairBcsdf::u(float x, float v) const {
    return v * std::log(exp(1.0f/v) - 2 * x * std::sinh(1.0f / v));
}

float HairBcsdf::OurSampleM(float v, float thetaCone, float x1, float x2) const {
    float theta_prime = PI / 2.0f - thetaCone;
    float u_x1 = u(x1, v);
    return u_x1 * cos(theta_prime) + sqrt(1 - pow(u_x1, 2)) * cos(2.0f * PI * x2) * sin(theta_prime);
}

float HairBcsdf::csch (float theta) const {
    return 1.0f / sinh(theta);
}

float HairBcsdf::OurM(float v, float sinThetaC, float sinThetaO, float cosThetaC, float cosThetaO) const {
    float a = csch(1.0f / v) / (2.0f * v);
    float b = exp(sinThetaC * sinThetaO / v);
    float c = cosThetaC * cosThetaO / v;
    float d = 1.0f / cosThetaO;
    float io = I0(c);
    return a * b * io * c * d;
}

float HairBcsdf::Phi(float gammaI, float gammaT, int p)
{
    return 2.0f*p*gammaT - 2.0f*gammaI + p*PI;
}

bool HairBcsdf::sample(SurfaceScatterEvent &event) const
{
    if (!event.requestedLobe.test(BsdfLobes::GlossyLobe))
        return false;
    // Random samples
    Vec2f x1 = event.sampler->next2D();
    Vec2f x2 = event.sampler->next2D();

    Vec3f wi = event.wi;

    float sinThetaI = wi.y();
    float thetaI = std::asin(sinThetaI);
    float cosThetaI = std::cos(thetaI);
    float thetaT = asin(sinThetaI / Eta);

    // Theta cone values for different lobes
    float thetaCR   = -thetaI + 2.0f*_scaleAngleRad;
    float thetaCTT  = -thetaI -      _scaleAngleRad;
    float thetaCTRT = -thetaI - 4.0f*_scaleAngleRad;

    // Lobe selection
    float v;
    float thetaC;
    float h = 2.0f * event.sampler->next1D() - 1.0f;
    Vec3f wo_spec = Vec3f(-event.wi.x(), -event.wi.y(), event.wi.z());
    float sinThetaO_spec = wo_spec.y();
    float thetaO_spec = asin(sinThetaO_spec);
    float thetaD_spec = (thetaO_spec - thetaI) / 2.0f;
    float cosThetaT = cos(thetaT);
    float cosThetaD_spec = cos(thetaD_spec);
    float wR   = A(pR, h, cosThetaT, cosThetaD_spec);
    float wTT  = A(pTT, h, cosThetaT, cosThetaD_spec);
    float wTRT = A(pTRT, h, cosThetaT, cosThetaD_spec);
    
    float p;
    float target = x2.x()*(wR + wTT + wTRT);
    if (target < wR) {
        v = _vR;
        thetaC = thetaCR;
        p = pR;
    } else if (target < wR + wTT) {
        v = _vTT;
        thetaC = thetaCTT;
        p = pTT;
    } else {
        v = _vTRT;
        thetaC = thetaCTRT;
        p = pTRT;
    }

    //Longitudinal sampling
    float sinThetaO = OurSampleM(v, thetaC, x1.x(), x1.y());
    float cosThetaO = trigInverse(sinThetaO);

    float thetaO = std::asin(sinThetaO);
    float thetaD = (thetaO - thetaI)*0.5f;
    float cosThetaD = std::cos(thetaD);

    // Azimuthal Sampling
    //lobe->sample(cosThetaD, x2.y(), phi, phiPdf);
    
    float etaPrime = std::sqrt(Eta*Eta - (1.0f - cosThetaD*cosThetaD))/cosThetaD;
    float gammaI = std::asin(h);
    float gammaT = std::asin(h/etaPrime);
    float phi = Phi(gammaI, gammaT, p) + G(event.sampler->next2D()); 
    //float phiPdf = PhiPdf();

    float sinPhi = std::sin(phi);
    float cosPhi = std::cos(phi);

    event.wo = Vec3f(sinPhi*cosThetaO, sinThetaO, cosPhi*cosThetaO);
    event.pdf = pdf(event);
    event.weight = eval(event)/event.pdf;
    event.sampledLobe = BsdfLobes::GlossyLobe;

    return true;
}
}
