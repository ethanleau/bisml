
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIAL_LTC_H
#define PBRT_MATERIAL_LTC_H

// material/ltc.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"

namespace pbrt {

struct SphericalCap
{
    SphericalCap(Vector3f dir, Float cosTheta) : dir(dir), cosTheta(cosTheta) {}

    Vector3f dir;
    Float cosTheta;
};

struct BRDFRecord
{
    Transform world2Local;
    Vector3f n, ng;
    Vector3f localWo;
    Float Kd = 0, Ks = 0, alpha = 1.f;
    Float specBrdfScale;
    Vector3f specPivot;
    SphericalCap specVisCap{Vector3f(0.f, 0.f, 1.f), 0.f};
};

struct BoundingSphere;
class LTC
{
  public:
    static bool GetBRDFRecord(const SurfaceInteraction &isect, BRDFRecord &record);
    static Float EvaluateLTCIntegral(const SurfaceInteraction &isect, const BRDFRecord &record, Point3f P[3]);
    static Float EvaluateLTC(const SurfaceInteraction &isect, const BRDFRecord &record, Point3f P[3]);
    static Float EvaluatePivot(const SurfaceInteraction &isect, const BRDFRecord &record, const BoundingSphere &boundingSphere);
    static Float EvaluatePivotIntegral(const SurfaceInteraction &isect, const BRDFRecord &record, const BoundingSphere &boundingSphere);

  private:
    static Vector3f IntegrateEdgeVec(const Vector3f &p1, const Vector3f &p2);
    static Vector3f IntegrateEdgeVectTest(const Vector3f &p1, const Vector3f &p2);
    static Vector3f ExtractPivot(const Vector3f &wo, Float alpha, Float &brdfScale);
    static Vector3f MaxValueVec(const SphericalCap &cap, const Vector3f &pivot);
    static SphericalCap Cap2PCap(const SphericalCap &cap, const Vector3f &pivot);
    static Vector2f Vec2PVec(const Vector2f &r, Float rp);
    static Vector3f Vec2PVec(const Vector3f &r, const Vector3f &rp);
    static Vector3f U2Cap(const Point2f &u, const SphericalCap &cap);
    static Float PdfCap(const Vector3f &wk, const SphericalCap &cap);
    static Float PivotJacobian(const Vector3f &wk, const Vector3f &rp);
    static Float CapSolidAngle(const SphericalCap &cap);
    static Float CapSolidAngle(const SphericalCap &cap1, const SphericalCap &cap2);
    static int ClipQuadToHorizon(Vector3f L[5]);

};

}  // namespace pbrt

#endif  // PBRT_MATERIAL_LTC_H