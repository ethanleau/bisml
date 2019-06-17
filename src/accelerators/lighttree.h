
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_LIGHTTREE_H
#define PBRT_ACCELERATORS_LIGHTTREE_H

// accelerators/lighttree.h*
#include "pbrt.h"
#include "light.h"
#include "primitive.h"

namespace pbrt {
struct LightNode;

struct LightInfo;

struct LightCone
{
    LightCone() : axis(Vector3f(0.f, 0.f, 0.f)), theta_o(0.f), theta_e(0.f) {}
    LightCone(const Vector3f &a, Float o, Float e) : axis(a), theta_o(o), theta_e(e) {}
    LightCone(const LightCone &c2)
    {
        axis = c2.axis;
        theta_o = c2.theta_o;
        theta_e = c2.theta_e;
    }
    Float Measure() const
    {
        if (axis == Vector3f(0.f, 0.f, 0.f))
            return 0.f;
        Float cos_theta_o = std::cos(theta_o);
        Float sin_theta_o = std::sin(theta_o);
        Float theta_w = std::min(theta_o + theta_e, Pi);
        return 2 * Pi * (1 - cos_theta_o) +
               PiOver2 * (2 * theta_w * sin_theta_o - std::cos(theta_o - 2 * theta_w) - 2 * theta_o * sin_theta_o + cos_theta_o);
    }
    static LightCone Union(const LightCone &cone0, const LightCone &cone1)
    {
        LightCone const *c0 = &cone0;
        LightCone const *c1 = &cone1;
        if (c0->axis == Vector3f(0.f, 0.f, 0.f))
            return *c1;
        if (c1->axis == Vector3f(0.f, 0.f, 0.f))
            return *c0;

        if (c1->theta_o > c0->theta_o)
            std::swap(c0, c1);

        Float theta_e = std::max(c0->theta_e, c1->theta_e);
        Float theta_d = std::acos(Clamp(Dot(c0->axis, c1->axis), -1, 1));
        if (std::min(theta_d + c1->theta_o, Pi) <= c0->theta_o)
            return LightCone(c0->axis, c0->theta_o, theta_e);

        Float theta_o = (c0->theta_o + c1->theta_o + theta_d) / 2;
        if (Pi <= theta_o)
            return LightCone(c0->axis, Pi, theta_e);

        Vector3f pivot = Cross(c0->axis, c1->axis);
        if (pivot.LengthSquared() > MachineEpsilon)
            return LightCone(Rotate(Degrees(theta_o - c0->theta_o), pivot)(c0->axis), theta_o, theta_e);
        else
        {
            if (theta_d < PiOver2)
                return LightCone(c0->axis, c0->theta_o, theta_e);
            else
                return LightCone(c0->axis, Pi, theta_e);
        }
    }

    Vector3f axis;
    Float theta_o, theta_e;
};

struct BoundingSphere
{
    Point3f c;
    Float r;
};

struct LinearLightNode
{
    bool isLeaf() const { return axis == 3; }

    BoundingSphere boundingSphere;      // 16 bytes
    LightCone cone;                     // 20 bytes
    Float power;                        // 4 bytes
    Float power2;                       // 4 bytes              for compute var of power
    union
    {
        int lightsOffset;       // leaf
        int secondChildOffset;  // interior
    };                                  // 4 bytes
    int nLights;                        // 4 bytes
    int axis;                           // 4 bytes, leaf -> 3
    union
    {
        uint64_t bit;           // leaf
        uint64_t depth;         // interior
    };                                  // 8 bytes
};

// LightTree Declarations
class LightTree
{
  public:
    // LightTree Public Types
    enum class SplitMethod { SAOH, Middle, EqualCounts };

    // LightTree Public Methods
    LightTree(std::vector<std::shared_ptr<Light>> &l, int maxLightsInNode = 1, SplitMethod splitMethod = SplitMethod::SAOH);
    ~LightTree();

    // LightTree Public Data
    LinearLightNode *nodes = nullptr;

  private:
    // LightTree Private Methods
    LightNode *RecursiveBuild(MemoryArena &arena, std::vector<LightInfo> &lightInfo, int start, int end, int *totalNodes,
                              std::vector<std::shared_ptr<Light>> &orderdLights);
    int FlattenLightTree(LightNode *node, int *offset, uint64_t bit, int depth, std::vector<std::shared_ptr<Light>> &orderdLights);

    // LightTree Private Data
    const int maxLightsInNode;
    const SplitMethod splitMethod;
    std::vector<std::shared_ptr<Light>> &lights;
    size_t notConsideredOffset = 0;
};

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_LIGHTTREE_H