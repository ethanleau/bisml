
// accelerators/lighttreesampler.cpp*
#include "lighttreesampler.h"
#include "sampler.h"
#include "reflection.h"
#include "interaction.h"
#include "scene.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

static PBRT_CONSTEXPR BxDFType bsdfFlags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
static PBRT_CONSTEXPR Float DISTANCE_EPSILON = 1e-9;

// 3 sampling techniques.
#define NAIVE_LIGHT_SELECTION                   // traditional light sampling technique
#define BRDF_ORIENTED_LIGHT_SELECTION           // BRDF-oriented light sampling technique
#define NAIVE_BRDF_SAMPLING                     // traditional BRDF sampling

// Adaptive sample allocation
#define ADAPTIVE_SAMPLE_ALLOCATION

// Light cluster sampling (sample one cluster according to pdf) or light cluster summing (consider every cluster and sum contributions of them)
//#define NAIVE_LIGHT_SELECTION_CLUSTER_SAMPLING        // for tradtional light sampling
#define BRDF_ORIENTED_LIGHT_SELECTION_CLUSTER_SAMPLING  // for BRDF-oriented light sampling

LightTreeSampler::LightTreeSampler(const Scene &scene, Sampler *sampler, const SurfaceInteraction &isect, const Float splitThreshold,
                                   const Float brdfThresold)
    : scene(scene), lightNodes(scene.lightTree->nodes), sampler(sampler), isect(isect), irradianceSplitThreshold(splitThreshold),
      brdfSplitThreshold(brdfThresold)
{
    success = LTC::GetBRDFRecord(isect, record);

#ifdef BRDF_ORIENTED_LIGHT_SELECTION
    willBRDFSelection = record.Ks > 0;
#else
    willBRDFSelection = false;
#endif
}

Spectrum LightTreeSampler::Evaluate()
{
    if (!success)
        return Spectrum(0.f);

    Spectrum Ld(0.f);

    /*--------------------------------------Light Splitting To Get Light Clusters----------------------------------*/
    /*-----Traditional Light Sampling-----*/
    {
#ifdef NAIVE_LIGHT_SELECTION
        ProfilePhase _(Prof::IrradiaceOrientedLightSampling);
        RecursiveSplit(0, BRDF_NONE);   // search a light cut
#ifdef NAIVE_LIGHT_SELECTION_CLUSTER_SAMPLING
        std::vector<Float> probs(cuts[BRDF_NONE].size());
        for (int i = 0; i < cuts[BRDF_NONE].size(); ++i)
            probs[i] = ClusterImportance(lightNodes[cuts[BRDF_NONE][i]], BRDF_NONE);
        distribs[BRDF_NONE] = std::unique_ptr<Distribution1D>(new Distribution1D(&probs[0], probs.size()));
#else
        numSamples[BRDF_NONE] = 1;  // for cluster summing, sample count is always 1
#endif
#endif
    }
    /*-----BRDF-Oriented Light Sampling-----*/
    if (willBRDFSelection)
    {
        ProfilePhase _(Prof::BRDFOrientedLightSampling);
        RecursiveSplit(0, BRDF_ONLY);  // search a light cut
#ifdef BRDF_ORIENTED_LIGHT_SELECTION_CLUSTER_SAMPLING
        std::vector<Float> probs(cuts[BRDF_ONLY].size());
        for (int i = 0; i < cuts[BRDF_ONLY].size(); ++i)
            probs[i] = ClusterImportance(lightNodes[cuts[BRDF_ONLY][i]], BRDF_ONLY);
        distribs[BRDF_ONLY] = std::unique_ptr<Distribution1D>(new Distribution1D(&probs[0], probs.size()));
#else
        numSamples[BRDF_ONLY] = 1;  // for cluster summing, sample count is always 1
#endif
    }

    /*----------------------------------------Direct Illumination Sampling------------------------------------------*/
#if defined(ADAPTIVE_SAMPLE_ALLOCATION) && defined(NAIVE_BRDF_SAMPLING) // adaptive sample allocation
    if (willBRDFSelection)
    {
        Float integral = 0.f;
        {
            // Compute the upper bound of BRDF integral (only consider clusters' to improve efficiency)
            ProfilePhase _(Prof::SampleAllocationForBRDFAware);
            for (int i = 0; i < cuts[BRDF_ONLY].size(); ++i)
            {
                Float e = brdfIntegrals[i];
                if (e < 0)
                    e = LTC::EvaluatePivotIntegral(isect, record, lightNodes[cuts[BRDF_ONLY][i]].boundingSphere);
                if (e > integral)
                    integral = e;
            }
            integral = Clamp(integral, 0.f, 1.f);
        }
        int totalNum = integral * 50;   // get total sample count of BRDF-aware methods

        // Sample allocation between BRDF-oriented light sampling and BRDF sampling
        int Nfail = 1, Nsucess = 1;
        // at least once for each sampling technqie
        BRDFOrientedLightSelection();
        if (NaiveBRDFSampling())
            Nsucess += 10;
        else
            ++Nfail;
        numSamples[BRDF_ONLY] = numBRDFSamples = 1;
        // allocate recursively
        Float u = sampler->Get1D();
        for (int step = 0; step < totalNum; ++step)
        {
            Float p = static_cast<Float>(Nfail) / static_cast<Float>(Nfail + Nsucess);
            if (u < p)
            {
                BRDFOrientedLightSelection();
                ++numSamples[BRDF_ONLY];
                u /= p;
            }
            else
            {
                if (NaiveBRDFSampling())
                    Nsucess += 10;
                else
                    ++Nfail;
                ++numBRDFSamples;
                u = (u - p) / (1.f - p);
            }
        }

        //// Visualization of sample allocation result
        // int totalSamples = numBRDFSamples + numSamples[BRDF_ONLY];
        // Float ratio = static_cast<Float>(numBRDFSamples) / totalSamples;
        // return ScalarVisualization(ratio);
    }
    else
    {
        numSamples[BRDF_ONLY] = 0;
        for (int step = 0; step < numBRDFSamples; ++step)
            NaiveBRDFSampling();
    }
#else
    /*-----BRDF-Oriented Light Sampling-----*/
    if (willBRDFSelection)
    {
        for (int step = 0; step < numSamples[BRDF_ONLY]; ++step)
            BRDFOrientedLightSelection();
    }
    else
        numSamples[BRDF_ONLY] = 0;
    /*-----Traditional BRDF Sampling-----*/
    {
#ifdef NAIVE_BRDF_SAMPLING
        for (int step = 0; step < numBRDFSamples; ++step)
            NaiveBRDFSampling();
#else
        numBRDFSamples = 0;
#endif
    }
#endif
    /*-----Traditional Light Sampling-----*/
    {
#ifdef NAIVE_LIGHT_SELECTION
        for (int step = 0; step < numSamples[BRDF_NONE]; ++step)
            NaiveLightSelection();
#else
        numSamples[BRDF_NONE] = 0;
#endif
    }

    /*----------------------------------------Direct Illumination Computing using MIS------------------------------------------*/
    if (!evaluateInfo[BRDF_NONE].empty())
    {
        ProfilePhase _(Prof::IrradiaceOrientedLightSampling);
        Spectrum L0(0.f);
        for (auto &entry : evaluateInfo[BRDF_NONE])
        {
            Float pdf2 = 0.f;
            Float scatteringPdf = 0.f;
            if (numSamples[BRDF_ONLY] > 0)
            {
                pdf2 = Pdf(entry.nodeID, BRDF_ONLY) * entry.lightPdf;
            }
            if (!entry.delta)
            {
                scatteringPdf = entry.scatteringPdf;
            }
            Float weight = PowerHeuristic(numSamples[BRDF_NONE], entry.pdf, numSamples[BRDF_ONLY], pdf2, numBRDFSamples, scatteringPdf);
            L0 += entry.L * weight / entry.pdf;
        }
        Ld += L0 / numSamples[BRDF_NONE];
    }
    if (!evaluateInfo[BRDF_ONLY].empty())
    {
        ProfilePhase _(Prof::BRDFOrientedLightSampling);
        Spectrum L1(0.f);
        for (auto &entry : evaluateInfo[BRDF_ONLY])
        {
            Float pdf2 = 0.f;
            Float scatteringPdf = 0.f;
            if (numSamples[BRDF_NONE] > 0)
            {
                pdf2 = Pdf(entry.nodeID, BRDF_NONE) * entry.lightPdf;
            }
            if (!entry.delta)
            {
                scatteringPdf = entry.scatteringPdf;
            }
            Float weight = PowerHeuristic(numSamples[BRDF_ONLY], entry.pdf, numSamples[BRDF_NONE], pdf2, numBRDFSamples, scatteringPdf);
            L1 += entry.L * weight / entry.pdf;
        }
        Ld += L1 / numSamples[BRDF_ONLY];
    }
    if (!evaluateInfo[2].empty())
    {
        ProfilePhase _(Prof::NaiveBRDFSampling);
        Spectrum L2(0.f);
        for (auto &entry : evaluateInfo[2])
        {
            Float pdf1 = 0.f;
            Float pdf2 = 0.f;
            if (numSamples[BRDF_NONE] > 0)
            {
                pdf1 = Pdf(entry.nodeID, BRDF_NONE) * entry.lightPdf;
            }
            if (numSamples[BRDF_ONLY] > 0)
            {
                pdf2 = Pdf(entry.nodeID, BRDF_ONLY) * entry.lightPdf;
            }
            Float weight = PowerHeuristic(numBRDFSamples, entry.pdf, numSamples[BRDF_NONE], pdf1, numSamples[BRDF_ONLY], pdf2);
            L2 += entry.L * weight / entry.pdf;
        }
        Ld += L2 / numBRDFSamples;
    }
    return Ld;
}

void LightTreeSampler::NaiveLightSelection()
{
    ProfilePhase _(Prof::IrradiaceOrientedLightSampling);
#ifdef NAIVE_LIGHT_SELECTION_CLUSTER_SAMPLING
    Float clusterPdf;
    int num = distribs[BRDF_NONE]->SampleDiscrete(sampler->Get1D(), &clusterPdf);   // sample a cluster
    RecursiveSample(cuts[BRDF_NONE][num], sampler->Get1D(), clusterPdf, BRDF_NONE); // sample a light in cluster
#else
    for (auto nodeID : cuts[BRDF_NONE]) // for each cluster
        RecursiveSample(nodeID, sampler->Get1D(), 1.f, BRDF_NONE);  // sample a light in cluster
#endif
}

void LightTreeSampler::BRDFOrientedLightSelection()
{
    ProfilePhase _(Prof::BRDFOrientedLightSampling);
#ifdef BRDF_ORIENTED_LIGHT_SELECTION_CLUSTER_SAMPLING
    Float clusterPdf;
    int num = distribs[BRDF_ONLY]->SampleDiscrete(sampler->Get1D(), &clusterPdf);   // sample a cluster
    RecursiveSample(cuts[BRDF_ONLY][num], sampler->Get1D(), clusterPdf, BRDF_ONLY); // sample a light in cluster
#else
    for (auto nodeID : cuts[BRDF_ONLY]) // for each cluster
        RecursiveSample(nodeID, sampler->Get1D(), 1.f, BRDF_ONLY);  // sample a light in cluster
#endif
}

bool LightTreeSampler::NaiveBRDFSampling()
{
    ProfilePhase _(Prof::NaiveBRDFSampling);
    Float scatteringPdf = 0;
    Vector3f wi;
    BxDFType sampledType;
    bool sampledSpecular = false;
    Spectrum f = isect.bsdf->Sample_f(isect.wo, &wi, sampler->Get2D(), &scatteringPdf, bsdfFlags, &sampledType);
    f *= AbsDot(wi, isect.shading.n);
    sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
    if (!f.IsBlack() && scatteringPdf > 0)
    {
        // Find intersection
        SurfaceInteraction lightIsect;
        Ray ray = isect.SpawnRay(wi);
        bool foundSurfaceInteraction = scene.Intersect(ray, &lightIsect);
        // Add light contribution
        if (foundSurfaceInteraction)
        {
            const Light *light = lightIsect.primitive->GetAreaLight();
            Spectrum Li = lightIsect.Le(-wi);
            if (!Li.IsBlack())
            {
                if (!sampledSpecular)
                {
                    evaluateInfo[2].push_back(EvaluateEntry(false, light->nodeID, f * Li, scatteringPdf, light->Pdf_Li(isect, wi)));
                }
                else
                    evaluateInfo[2].push_back(EvaluateEntry(true, light->nodeID, f * Li, scatteringPdf));

                evaluateSum[1] += Spectrum(Li * f / scatteringPdf).y();
                return true;
            }
        }
    }
    return false;
}

void LightTreeSampler::RecursiveSplit(int nodeID, BRDF_TYPE brdfType)
{
    ProfilePhase _(Prof::SearchLightCut);
    if (lightNodes[nodeID].isLeaf())
    {
        cuts[brdfType].push_back(nodeID);
        if (brdfType == BRDF_ONLY)
            brdfIntegrals.push_back(-1.f);
    }
    else
    {
        const BoundingSphere &boundingSphere = lightNodes[nodeID].boundingSphere;
        Float d2 = (boundingSphere.c - isect.p).LengthSquared();
        if (brdfType == BRDF_ONLY)
        {
            Float integral = -1.f;
            // Splitting metric of BRDF-oriented light selection
            if (d2 <= boundingSphere.r * boundingSphere.r ||
                (integral = LTC::EvaluatePivotIntegral(isect, record, boundingSphere)) > brdfSplitThreshold)
            {
                RecursiveSplit(nodeID + 1, brdfType);
                RecursiveSplit(lightNodes[nodeID].secondChildOffset, brdfType);
            }
            else
            {
                cuts[brdfType].push_back(nodeID);
                brdfIntegrals.push_back(integral);
            }
        }
        else
        {
            // Splitting metric of traditional light sampling
            if (d2 <= boundingSphere.r * boundingSphere.r ||
                ClusterVariance(lightNodes[nodeID], boundingSphere.c, std::sqrt(d2), boundingSphere.r) < irradianceSplitThreshold)
            {
                RecursiveSplit(nodeID + 1, brdfType);
                RecursiveSplit(lightNodes[nodeID].secondChildOffset, brdfType);
            }
            else  // Traverse
                cuts[brdfType].push_back(nodeID);
        }
    }
}

void LightTreeSampler::RecursiveSample(int nodeID, Float u, Float pdf, BRDF_TYPE brdfType)
{
    if (lightNodes[nodeID].isLeaf())
    {
        // Evaluate
        EvaluateLeafNode(lightNodes[nodeID], pdf, brdfType);
        return;
    }

    int leftID = nodeID + 1;
    int rightID = lightNodes[nodeID].secondChildOffset;
    Float leftMeasure = ClusterImportance(lightNodes[leftID], brdfType);
    Float rightMeasure = ClusterImportance(lightNodes[rightID], brdfType);
    Float totalMeasure = leftMeasure + rightMeasure;
    Float p;
    if (totalMeasure > 0.f)
        p = leftMeasure / totalMeasure;
    else
        return;

    if (u < p)
    {
        RecursiveSample(leftID, u / p, pdf * p, brdfType);
    }
    else
        RecursiveSample(rightID, (u - p) / (1.f - p), pdf * (1.f - p), brdfType);
}

void LightTreeSampler::EvaluateLeafNode(const LinearLightNode &lightNode, Float pdf, BRDF_TYPE brdfType)
{
    if (pdf == 0)
        return;

    // Assumption: each leaf node only has one light
    CHECK_EQ(lightNode.nLights, 1);
    // pdf *= 1.0f / lightNode.nLights;
    // int lightNumber = lightNode.lightsOffset + std::min(int(sampler->Get1D() * lightNode.nLights), lightNode.nLights - 1);
    const std::shared_ptr<Light> &light = scene.lights[lightNode.lightsOffset];

    // Sample a point on the light
    {
        Float lightPdf = 0;
        Vector3f wi;
        VisibilityTester visibility;
        Spectrum Li = light->Sample_Li(isect, sampler->Get2D(), &wi, &lightPdf, &visibility);
        if (lightPdf > 0 && !Li.IsBlack())
        {
            // Evaluate BSDF
            Spectrum f = isect.bsdf->f(isect.wo, wi, bsdfFlags) * AbsDot(wi, isect.shading.n);
            if (!f.IsBlack() && visibility.Unoccluded(scene))
            {
                if (IsDeltaLight(light->flags))
                {
                    evaluateInfo[brdfType].push_back(EvaluateEntry(true, light->nodeID, Li * f, pdf * lightPdf, lightPdf));
                }
                else
                {
                    evaluateInfo[brdfType].push_back(
                        EvaluateEntry(false, light->nodeID, Li * f, pdf * lightPdf, lightPdf, isect.bsdf->Pdf(isect.wo, wi, bsdfFlags)));
                }
                if (brdfType == BRDF_ONLY)
                    evaluateSum[0] += Spectrum(Li * f / (pdf * lightPdf)).y();
            }
        }
    }
}

Float LightTreeSampler::ClusterImportance(const LinearLightNode &lightNode, BRDF_TYPE brdfType) const
{
    ProfilePhase _(Prof::ComputeImportaceOfLightNode);

    Float cosX, cosY;

    const BoundingSphere &boundingSphere = lightNode.boundingSphere;
    Vector3f c = boundingSphere.c - isect.p;
    Float d2 = c.LengthSquared();
    Float r2 = boundingSphere.r * boundingSphere.r;
    if (d2 < r2)
    {
        d2 = std::max(d2, r2 * 0.25f);
        cosX = cosY = 1.f;
    }
    else
    {
        Float d = std::sqrt(d2);
        c = c / d;
        // Compute cosX
        Float sinTheta_u = boundingSphere.r / d;
        Float cosTheta_u = std::sqrt(1 - sinTheta_u * sinTheta_u);
        Float cosTheta_i = Clamp(Dot(record.ng, c), -1, 1);
        if (cosTheta_i >= cosTheta_u)
            cosX = 1.f;
        else
        {
            Float sinTheta_i = std::sqrt(1 - cosTheta_i * cosTheta_i);
            cosX = cosTheta_i * cosTheta_u + sinTheta_i * sinTheta_u;
            if (cosX <= 0)
                return 0;
        }
        // Compute cosY
        Float cosTheta = Clamp(Dot(lightNode.cone.axis, -c), -1, 1);
        Float cosTheta_o = std::cos(lightNode.cone.theta_o);
        if (cosTheta >= cosTheta_o)
            cosY = 1.f;
        else
        {
            Float sinTheta = std::sqrt(1 - cosTheta * cosTheta);
            Float sinTheta_o = std::sqrt(1 - cosTheta_o * cosTheta_o);
            Float cosThetaMinusTheta_o = cosTheta * cosTheta_o + sinTheta * sinTheta_o;
            if (cosThetaMinusTheta_o >= cosTheta_u)
                cosY = 1.f;
            else
            {
                Float sinThetaMinusTheta_o = sinTheta * cosTheta_o - cosTheta * sinTheta_o;
                cosY = cosThetaMinusTheta_o * cosTheta_u + sinThetaMinusTheta_o * sinTheta_u;
                if (cosY <= std::cos(lightNode.cone.theta_e))
                    return 0;
            }
        }
    }

    Float e = lightNode.power * /*cosX **/ cosY / d2;
    if (brdfType == BRDF_ONLY)
    {
        //Float diffuse = record.Kd * InvPi * cosX;
        Float glossy;
        glossy = LTC::EvaluatePivot(isect, record, boundingSphere);

        e *= glossy;
    }
    else
        e *= cosX;

    return e;
}

Float LightTreeSampler::ClusterVariance(const LinearLightNode &lightNode, const Point3f &clusterCenter, Float d, Float clusterRadius) const
{
    Float meanPower = lightNode.power / lightNode.nLights;
    Float varPower = lightNode.nLights > 1 ?  1 / (lightNode.nLights - 1) * (lightNode.power2 - meanPower * meanPower * lightNode.nLights) : 0;
    // Compute mean & var of geometric term
    Float a = std::max(DISTANCE_EPSILON, d - clusterRadius);
    Float b = d + clusterRadius;
    Float meanGeo = 1.0f / (a * b);
    Float a2 = a * a;
    Float b2 = b * b;
    Float a3 = a2 * a;
    Float b3 = b2 * b;
    Float varGeo = (b3 - a3) / (3 * (b - a) * a3 * b3) - 1.0f / (a2 * b2);
    // compute sigma2
    Float sigma2 =
        (varPower * varGeo + varPower * meanGeo * meanGeo + meanPower * meanPower * varGeo) * lightNode.nLights * lightNode.nLights;
    return std::pow(1.0f / (1.0f + std::sqrt(sigma2)), 1.0f / 4);
}

Float LightTreeSampler::Pdf(int leafNodeID, BRDF_TYPE brdfType) const
{
    ProfilePhase _(Prof::ComputeLightTreeSamplingPdf);

    int num = cuts[brdfType].size() - 1;
    // Find the cluster this leaf node belongs to
    for (int i = 0; i < cuts[brdfType].size(); ++i)
    {
        if (leafNodeID < cuts[brdfType][i])
        {
            num = i - 1;
            break;
        }
    }
    int clusterNodeID = cuts[brdfType][num];
    uint64_t bit = lightNodes[leafNodeID].bit >> lightNodes[clusterNodeID].depth;

    Float pdf = 1.f;
    if (brdfType == BRDF_NONE)
    {
#ifdef NAIVE_LIGHT_SELECTION_CLUSTER_SAMPLING
        pdf = distribs[brdfType]->DiscretePDF(num);
#endif
    }
    else
    {
#ifdef BRDF_ORIENTED_LIGHT_SELECTION_CLUSTER_SAMPLING
        pdf = distribs[brdfType]->DiscretePDF(num);
#endif
    }

    return PdfSample(clusterNodeID, bit, pdf, brdfType);
}

Float LightTreeSampler::PdfSample(int nodeID, uint64_t bit, Float pdf, BRDF_TYPE brdfType) const
{
    if (lightNodes[nodeID].isLeaf())
        return pdf;
    else
    {
        int leftID = nodeID + 1;
        int rightID = lightNodes[nodeID].secondChildOffset;
        Float leftMeasure = ClusterImportance(lightNodes[leftID], brdfType);
        Float rightMeasure = ClusterImportance(lightNodes[rightID], brdfType);
        Float totalMeasure = leftMeasure + rightMeasure;
        if (totalMeasure == 0)
            return 0;
        Float p = leftMeasure / totalMeasure;
        if ((bit & 0x1) == 0)
        {
            nodeID = leftID;
            pdf *= p;
        }
        else
        {
            nodeID = rightID;
            pdf *= 1 - p;
        }
        bit = bit >> 1;

        return PdfSample(nodeID, bit, pdf, brdfType);
    }
}

Spectrum LightTreeSampler::ScalarVisualization(Float e) const
{
    Float rgb[3];
    Float a = (1 - e) / 0.25f;
    int X = std::floor(a);
    float Y = a - X;
    switch (X)
    {
    case 0:
        rgb[0] = 1;
        rgb[1] = Y;
        rgb[2] = 0;
        break;
    case 1:
        rgb[0] = 1 - Y;
        rgb[1] = 1;
        rgb[2] = 0;
        break;
    case 2:
        rgb[0] = 0;
        rgb[1] = 1;
        rgb[2] = Y;
        break;
    case 3:
        rgb[0] = 0;
        rgb[1] = 1 - Y;
        rgb[2] = 1;
        break;
    case 4:
        rgb[0] = 0;
        rgb[1] = 0;
        rgb[2] = 1;
        break;
    }
    return Spectrum::FromRGB(rgb);
}
}  // namespace pbrt
