#include "OurHairBcsdf.hpp"

#include "sampling/PathSampleGenerator.hpp"

#include "bsdfs/Fresnel.hpp"

#include "math/GaussLegendre.hpp"

#include "io/JsonObject.hpp"

namespace Tungsten {

    namespace OurHairBcsdf {

        static float csch (float theta) {
            //2/(e^x - e^(-x))
            return 2.0f / (exp(theta) - exp(-theta));
        }

        static float thetaConeR (float alpha, float thetaI) {
            return -thetaI + 2.0f*alpha;
        }

        static float thetaConeTT (float alpha, float thetaI) {
            return -thetaI - alpha;
        }

        static float thetaConeTRT (float alpha, float thetaI) {
            return -thetaI - 4.0f*alpha;
        }

        static float OurM(float v, float sinThetaC, float sinThetaO, float cosThetaC, float cosThetaO, float IO) {
            float a = csch(1.0f / v) / (2.0f * v);
            float b = exp(sinThetaC * sinThetaO / v);
            float c = cosThetaC * sinThetaO / v;
            float d = 1.0f / cosThetaO;
            return a * b * IO * c * d;
        }
    }

}
