
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_LIGHTTREESAMPLER_H
#define PBRT_ACCELERATORS_LIGHTTREESAMPLER_H

// accelerators/lighttreesampler.h*
#include "pbrt.h"
#include "lighttree.h"
#include "materials/ltc.h"

namespace pbrt {

enum BRDF_TYPE 
{ 
    BRDF_NONE = 0, // traditional light sampling
    BRDF_ONLY = 1  // BRDF-oriented light sampling
};

class LightTreeSampler
{
  public:
    LightTreeSampler(const Scene &scene, Sampler *sampler, const SurfaceInteraction &isect, const Float splitThreshold,
                     const Float brdfThresold);
    Spectrum Evaluate();
    inline int GetLightCutSize() const { return cuts[0].size(); }
    inline int GetBRDFCutSize() const { return cuts[1].size(); }
    inline int GetBRDFWiseSampleCount() const { return numSamples[BRDF_ONLY] + numBRDFSamples; }
    inline bool DidBRDFClustering() const { return willBRDFSelection; }

  private:
    void NaiveLightSelection();
    void BRDFOrientedLightSelection();
    bool NaiveBRDFSampling();
    void RecursiveSplit(int nodeID, BRDF_TYPE brdfType);    // split light node
    void RecursiveSample(int nodeID, Float u, Float pdf, BRDF_TYPE brdfType);   // sample light node
    void EvaluateLeafNode(const LinearLightNode &lightNode, Float pdf, BRDF_TYPE brdfType);
    Float ClusterImportance(const LinearLightNode &lightNode, BRDF_TYPE brdfType) const;
    Float ClusterVariance(const LinearLightNode &lightNode, const Point3f &clusterCenter, Float d, Float clusterRadius) const;
    Float Pdf(int leafNodeID, BRDF_TYPE brdfType) const;
    Float PdfSample(int nodeID, uint64_t bit, Float pdf, BRDF_TYPE brdfType) const;
    Spectrum ScalarVisualization(Float e) const;

    const Scene &scene;
    LinearLightNode *lightNodes;        // flatten light tree
    Sampler *sampler;
    const SurfaceInteraction &isect;
    BRDFRecord record;                  // record BRDFs information of current shading point
    const Float irradianceSplitThreshold;
    const Float brdfSplitThreshold;

    bool success;
    std::vector<int> cuts[2];           // light cut. 0: traditional light sampling    1: BRDF-oriented light sampling
    std::vector<Float> brdfIntegrals;
    std::unique_ptr<Distribution1D> distribs[2];

    bool willBRDFSelection;             // whether use BRDF-oriented light sampling technique
    int numSamples[2] = {1, 1};         // sample count of 0: traditional light sampling    1: BRDF-oriented light sampling
    int numBRDFSamples = 1;             // sample count of traditional BRDF sampling

    struct EvaluateEntry
    {
        bool delta;
        int nodeID;
        Spectrum L;
        Float pdf;
        Float lightPdf;
        Float scatteringPdf;

        EvaluateEntry(bool delta, int nodeID, Spectrum L, Float pdf = 0, Float lightPdf = 0, Float scatteringPdf = 0)
            : delta(delta), nodeID(nodeID), L(L), pdf(pdf), lightPdf(lightPdf), scatteringPdf(scatteringPdf)
        {
        }
    };
    std::vector<EvaluateEntry> evaluateInfo[3];
    Float evaluateSum[2];
};

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_LIGHTTREESAMPLER_H