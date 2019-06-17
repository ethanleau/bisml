
// accelerators/lighttree.cpp*
#include "lighttree.h"
#include "stats.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Light tree", treeBytes);
STAT_RATIO("LightTree/Lights per leaf node", totalLights, totalLeafNodes);
STAT_COUNTER("LightTree/Interior nodes", interiorNodes);
STAT_COUNTER("LightTree/Leaf nodes", leafNodes);

// LightTree Local Declarations
struct LightInfo
{
    LightInfo() {}
    LightInfo(size_t lightNumber, Float power, const Bounds3f &bounds, Vector3f axis, Float theta_o, Float theta_e)
        : lightNumber(lightNumber), power(power), power2(power * power), bounds(bounds), centroid(.5f * bounds.pMin + .5f * bounds.pMax),
          cone(axis, theta_o, theta_e)
    {
    }
    size_t lightNumber;
    Float power, power2;
    Bounds3f bounds;
    Point3f centroid;
    LightCone cone;
};

struct LightNode
{
    // LightNode Public Mehthods
    void InitLeaf(int first, int n, Float p, Float p2, const Bounds3f &b, const LightCone &c)
    {
        firstLightOffset = first;
        nLights = n;
        power = p;
        power2 = p2;
        bounds = b;
        cone = c;
        children[0] = children[1] = nullptr;
        splitAxis = 3;
        ++leafNodes;
        ++totalLeafNodes;
        totalLights += n;
    }
    void InitInterior(int spltAx, LightNode *c0, LightNode *c1)
    {
        children[0] = c0;
        children[1] = c1;
        power = c0->power + c1->power;
        power2 = c0->power2 + c1->power2;
        nLights = c0->nLights + c1->nLights;
        bounds = Union(c0->bounds, c1->bounds);
        cone = LightCone::Union(c0->cone, c1->cone);
        splitAxis = spltAx;
        ++interiorNodes;
    }
    Float power, power2;
    Bounds3f bounds;
    LightCone cone;
    LightNode *children[2];
    int splitAxis, firstLightOffset, nLights;
};

LightTree::LightTree(std::vector<std::shared_ptr<Light>> &l, int maxLightsInNode /* = 1 */,
                     SplitMethod splitMethod /* = SplitMethod::SAOH */)
    : maxLightsInNode(std::min(255, maxLightsInNode)), splitMethod(splitMethod), lights(l)
{
    ProfilePhase _(Prof::LightTreeConstruction);
    if (lights.empty())
        return;

    // Initialize _lightInfo_ array for lights
    std::vector<LightInfo> lightInfo;
    std::vector<size_t> notConsiderdLightNumbers;
    for (size_t i = 0; i < lights.size(); ++i)
    {
        // Only consider triangle light, sphere light, disk light, point light and spot light right now
        Vector3f axis;
        Float theta_o, theta_e;
        if (lights[i]->GetOrientationAttributes(axis, theta_o, theta_e))
            lightInfo.push_back(LightInfo(i, lights[i]->LeWithArea(), lights[i]->WorldBound(), axis, theta_o, theta_e));
        else
            notConsiderdLightNumbers.push_back(i);
    }

    if (lightInfo.empty())
        return;

    // Build light tree for lights using _lightInfo_
    MemoryArena arena(1024 * 1024);
    int totalNodes = 0;
    std::vector<std::shared_ptr<Light>> orderedLights;
    orderedLights.reserve(lights.size());
    LightNode *root = RecursiveBuild(arena, lightInfo, 0, lightInfo.size(), &totalNodes, orderedLights);

    // Compute representation of depth-first traversal of light tree
    nodes = AllocAligned<LinearLightNode>(totalNodes);
    int offset = 0;
    FlattenLightTree(root, &offset, 0, 0, orderedLights);
    CHECK_EQ(totalNodes, offset);
    treeBytes += totalNodes * sizeof(LinearLightNode) + totalNodes * sizeof(uint32_t) + sizeof(*this);

    // push back not considered lights
    notConsideredOffset = orderedLights.size();
    for (auto number : notConsiderdLightNumbers)
        orderedLights.push_back(lights[number]);
    CHECK_EQ(lights.size(), orderedLights.size());
    lights.swap(orderedLights);  // need to test
    lightInfo.resize(0);
    LOG(INFO) << StringPrintf(
        "LightTree created with %d nodes for %d "
        "lights (%.2f MB), arena allocated %.2f MB",
        totalNodes, (int)lights.size(), float(totalNodes * sizeof(LightNode)) / (1024.f * 1024.f),
        float(arena.TotalAllocated()) / (1024.f * 1024.f));
}

LightTree::~LightTree()
{
    FreeAligned(nodes);
}

struct LightBucketInfo
{
    Float power = 0.0f;
    Bounds3f bounds;
    LightCone cone;
};

LightNode *LightTree::RecursiveBuild(MemoryArena &arena, std::vector<LightInfo> &lightInfo, int start, int end, int *totalNodes,
                                     std::vector<std::shared_ptr<Light>> &orderdLights)
{
    CHECK_NE(start, end);
    LightNode *node = arena.Alloc<LightNode>();
    (*totalNodes)++;
    // Compute bounds, cone and power of all lights in node
    Float power = 0, power2 = 0;
    Bounds3f bounds;
    LightCone cone;
    for (int i = start; i < end; ++i)
    {
        power += lightInfo[i].power;
        power2 += lightInfo[i].power2;
        bounds = Union(bounds, lightInfo[i].bounds);
        cone = LightCone::Union(cone, lightInfo[i].cone);
    }
    int nLights = end - start;
    if (nLights == 1)
    {
        // Create leaf _lightNode_
        int firstLightOffset = orderdLights.size();
        for (int i = start; i < end; ++i)
        {
            int lightNum = lightInfo[i].lightNumber;
            orderdLights.push_back(lights[lightNum]);
        }
        node->InitLeaf(firstLightOffset, nLights, power, power2, bounds, cone);
        return node;
    }
    else
    {
        // Compute bound of light centroids, choose split dimension _dim_
        Bounds3f centroidBounds;
        for (int i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, lightInfo[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition lights into two sets and build children
        int mid = (start + end) / 2;
        //if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim])
        //{
        //    // Create leaf _LightNode_
        //    int firstLightOffset = orderdLights.size();
        //    for (int i = start; i < end; ++i)
        //    {
        //        int lightNum = lightInfo[i].lightNumber;
        //        orderdLights.push_back(lights[lightNum]);
        //    }
        //    node->InitLeaf(firstLightOffset, nLights, power, power2, bounds, cone);
        //    return node;
        //}
        //else
        {
            // Partition lights based on _splitMethod_
            switch (splitMethod)
            {
            case SplitMethod::Middle:
            {
                // Partition lights through node's midpoint
                Float lmid = (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                LightInfo *midLight = std::partition(&lightInfo[start], &lightInfo[end - 1] + 1,
                                                     [dim, lmid](const LightInfo &li) { return li.centroid[dim] < lmid; });
                mid = midLight - &lightInfo[0];
                // For lots of lights with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall
                // through
                // to EqualCounts.
                if (mid != start && mid != end)
                    break;
            }
            case SplitMethod::EqualCounts:
            {
                // Partition lights into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&lightInfo[start], &lightInfo[mid], &lightInfo[end - 1] + 1,
                                 [dim](const LightInfo &a, const LightInfo &b) { return a.centroid[dim] < b.centroid[dim]; });
                break;
            }
            case SplitMethod::SAOH:
            default:
            {
                // Partition lights using approximate SAOH
                if (nLights <= 2)
                {
                    // Partition primitives into equally-sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&lightInfo[start], &lightInfo[mid], &lightInfo[end - 1] + 1,
                                     [dim](const LightInfo &a, const LightInfo &b) { return a.centroid[dim] < b.centroid[dim]; });
                }
                else
                {
                    // Allocate _BucketInfo_ for SAOH partition buckets
                    PBRT_CONSTEXPR int nBuckets = 12;
                    Float minCost = std::numeric_limits<Float>::max();
                    int minCostSplitBucket = 0;

                    // Explore all axes
                    Float maxLength = bounds.Diagonal()[bounds.MaximumExtent()];
                    for (int dim_i = 0; dim_i < 3; ++dim_i)
                    {
                        // Initiallize _BucketInfo_ for SAOH partition buckets
                        LightBucketInfo buckets[nBuckets];
                        for (int i = start; i < end; ++i)
                        {
                            int b = nBuckets * centroidBounds.Offset(lightInfo[i].centroid)[dim_i];
                            if (b == nBuckets)
                                b = nBuckets - 1;
                            CHECK_GE(b, 0);
                            CHECK_LT(b, nBuckets);
                            buckets[b].power += lightInfo[i].power;
                            buckets[b].bounds = Union(buckets[b].bounds, lightInfo[i].bounds);
                            buckets[b].cone = LightCone::Union(buckets[b].cone, lightInfo[i].cone);
                        }

                        // Compute costs for splitting after each bucket
                        Float k = maxLength / bounds.Diagonal()[dim_i];
                        for (int i = 0; i < nBuckets - 1; ++i)
                        {
                            Bounds3f b0, b1;
                            LightCone c0, c1;
                            Float power0 = 0.0f, power1 = 0.0f;
                            for (int j = 0; j <= i; ++j)
                            {
                                b0 = Union(b0, buckets[j].bounds);
                                c0 = LightCone::Union(c0, buckets[j].cone);
                                power0 += buckets[j].power;
                            }
                            for (int j = i + 1; j < nBuckets; ++j)
                            {
                                b1 = Union(b1, buckets[j].bounds);
                                c1 = LightCone::Union(c1, buckets[j].cone);
                                power1 += buckets[j].power;
                            }
                            Float cost = k * (power0 * b0.SurfaceArea() * c0.Measure() + power1 * b1.SurfaceArea() * c1.Measure()) /
                                         (bounds.SurfaceArea() * cone.Measure());
                            if (cost < minCost)
                            {
                                minCost = cost;
                                minCostSplitBucket = i;
                                dim = dim_i;
                            }
                        }
                    }
                    // Either create leaf or split lights at selected SAOH bucket
                    Float leafCost = power;
                    if (nLights > maxLightsInNode || minCost > leafCost)
                    {
                        // Split
                        LightInfo *lmid = std::partition(&lightInfo[start], &lightInfo[end - 1] + 1, [=](const LightInfo &li) {
                            int b = nBuckets * centroidBounds.Offset(li.centroid)[dim];
                            if (b == nBuckets)
                                b = nBuckets - 1;
                            CHECK_GE(b, 0);
                            CHECK_LT(b, nBuckets);
                            return b <= minCostSplitBucket;
                        });
                        mid = lmid - &lightInfo[0];
                    }
                    else
                    {
                        // Create leaf _LightNode_
                        int firstLightOffset = orderdLights.size();
                        for (int i = start; i < end; ++i)
                        {
                            int lightNum = lightInfo[i].lightNumber;
                            orderdLights.push_back(lights[lightNum]);
                        }
                        node->InitLeaf(firstLightOffset, nLights, power, power2, bounds, cone);
                        return node;
                    }
                }
                break;
            }
            }
            LightNode *leftNode = RecursiveBuild(arena, lightInfo, start, mid, totalNodes, orderdLights);
            LightNode *rightNode = RecursiveBuild(arena, lightInfo, mid, end, totalNodes, orderdLights);
            node->InitInterior(dim, leftNode, rightNode);
        }
    }
    return node;
}

int LightTree::FlattenLightTree(LightNode *node, int *offset, uint64_t bit, int depth, std::vector<std::shared_ptr<Light>> &orderdLights)
{
    LinearLightNode *linearNode = &nodes[*offset];
    node->bounds.BoundingSphere(&linearNode->boundingSphere.c, &linearNode->boundingSphere.r);
    linearNode->cone = node->cone;
    linearNode->power = node->power;
    linearNode->power2 = node->power2;
    int myOffset = (*offset)++;
    if (node->splitAxis == 3)
    {
        // Create leaf flattened light tree node
        CHECK(!node->children[0] && !node->children[1]);
        linearNode->lightsOffset = node->firstLightOffset;
        linearNode->nLights = node->nLights;
        linearNode->axis = 3;
        for (int i = 0; i < node->nLights; ++i)
        {
            orderdLights[node->firstLightOffset + i]->nodeID = myOffset;
            orderdLights[node->firstLightOffset + i]->arrayID = node->firstLightOffset + i;
        }
        linearNode->bit = bit;
    }
    else
    {
        // Create interior flattened light tree node
        linearNode->axis = node->splitAxis;
        linearNode->nLights = node->nLights;
        linearNode->depth = depth;
        FlattenLightTree(node->children[0], offset, bit, depth + 1, orderdLights);
        linearNode->secondChildOffset = FlattenLightTree(node->children[1], offset, (1ULL << depth) | bit, depth + 1, orderdLights);
    }
    return myOffset;
}

}  // namespace pbrt