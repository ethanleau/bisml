
// mateial/ltc.cpp*
#include "ltc.h"
#include "ltc_table.h"
#include "interaction.h"
#include "reflection.h"
#include "lowdiscrepancy.h"
#include "shapes/sphere.h"
#include "accelerators/lighttree.h"
#include "stats.h"

namespace pbrt {

static PBRT_CONSTEXPR Float LUT_SIZE = 64.0f;
static PBRT_CONSTEXPR Float LUT_SCALE = LUT_SIZE - 1.0f;
static PBRT_CONSTEXPR Float LUT_BIAS = 0.5f;

bool LTC::GetBRDFRecord(const SurfaceInteraction &isect, BRDFRecord &record)
{
    Vector3f T1, T2, N;
    if (Dot(isect.wo, isect.n) > 0)
    {
        N = Vector3f(isect.shading.n);
        record.ng = Vector3f(isect.n);
    }

    else
    {
        N = Vector3f(-isect.shading.n);
        record.ng = -Vector3f(isect.n);
    }
    Float cosTheta = std::min(AbsDot(isect.wo, isect.shading.n), 1.0f);
    if (cosTheta < 0.999f)
        T1 = Normalize(isect.wo - N * cosTheta);
    else
        T1 = Normalize(isect.shading.dpdu);
    T2 = Cross(N, T1);
    Matrix4x4 mat(T1[0], T1[1], T1[2], 0.f, T2[0], T2[1], T2[2], 0.f, N[0], N[1], N[2], 0.f, 0.f, 0.f, 0.f, 1.f);
    record.world2Local = Transform(mat);
    record.n = N;
    record.localWo = record.world2Local(isect.wo);

    Float alpha = 0, Kd = 0, Ks = 0;
    for (int i = 0; i < isect.bsdf->nBxDFs; ++i)
        isect.bsdf->bxdfs[i]->PrepareForLTC(alpha, Kd, Ks);

    record.Kd = Kd;

    if (Ks > 0)
    {
        record.Ks = Ks;
        record.alpha = alpha;
        record.specPivot = ExtractPivot(record.localWo, alpha, record.specBrdfScale);
        record.specVisCap = Cap2PCap(SphericalCap(Vector3f(0, 0, 1), 0.f), record.specPivot);
    }

    return true;
}

Float LTC::EvaluateLTCIntegral(const SurfaceInteraction &isect, const BRDFRecord &record, Point3f P[3])
{
    int s, t;
    s = (int)(std::sqrt(1.0f - CosTheta(record.localWo)) * LUT_SCALE + LUT_BIAS);
    t = (int)(std::sqrt(record.alpha) * LUT_SCALE + LUT_BIAS);
    auto &entry = ltcTable[s][t];
    Matrix4x4 Minv =
        Matrix4x4(entry[0], 0.0f, entry[2], 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, entry[1], 0.0f, entry[3], 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
    Float weight = entry[4] + (1 - record.Ks) * entry[5];

    Transform trans(Minv);
    trans = trans * record.world2Local;

    Vector3f L[3];
    for (int i = 0; i < 3; ++i)
    {
        L[i] = Normalize(trans(P[i] - isect.p));
    }
    Vector3f vsum;
    vsum += IntegrateEdgeVec(L[0], L[1]);
    vsum += IntegrateEdgeVec(L[1], L[2]);
    vsum += IntegrateEdgeVec(L[2], L[0]);
    vsum = -vsum;

    Float len = vsum.Length();
    if (len == 0)
        return 0;

    Float z = vsum.z / len;
    int u, v;
    v = (z * 0.5f + 0.5f) * LUT_SCALE + LUT_BIAS;
    u = len * LUT_SCALE + LUT_BIAS;
    Float scale = ltcTable[u][v][6];

    return weight * scale * len;
}

Float LTC::EvaluateLTC(const SurfaceInteraction &isect, const BRDFRecord &record, Point3f P[3])
{
    int s, t;
    s = (int)(std::sqrt(1.0f - CosTheta(record.localWo)) * LUT_SCALE + LUT_BIAS);
    t = (int)(std::sqrt(record.alpha) * LUT_SCALE + LUT_BIAS);
    auto &entry = ltcTable[s][t];
    Matrix4x4 Minv =
        Matrix4x4(entry[0], 0.0f, entry[2], 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, entry[1], 0.0f, entry[3], 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
    Float weight = record.Ks * entry[4] + (1 - record.Ks) * entry[5];

    Transform trans(Minv);
    trans = trans * record.world2Local;

    Vector3f L[3];
    for (int i = 0; i < 3; ++i)
    {
        L[i] = Normalize(trans(P[i] - isect.p));
    }
    Vector3f vsum;
    vsum += IntegrateEdgeVec(L[0], L[1]);
    vsum += IntegrateEdgeVec(L[1], L[2]);
    vsum += IntegrateEdgeVec(L[2], L[0]);
    vsum = -vsum;

    Float len = vsum.Length();
    if (len == 0)
        return 0;

    Float z = vsum.z / len;
    int u, v;
    v = (z * 0.5f + 0.5f) * LUT_SCALE + LUT_BIAS;
    u = len * LUT_SCALE + LUT_BIAS;
    Float scale = ltcTable[u][v][6];

    // Float nom = weight * scale * len;
    Float nom = weight * std::abs(z);
    // return weight * scale * len;

    Vector3f M[3];
    for (int i = 0; i < 3; ++i)
    {
        M[i] = Normalize(record.world2Local(P[i] - isect.p));
    }
    Vector3f usum;
    usum += IntegrateEdgeVec(M[0], M[1]);
    usum += IntegrateEdgeVec(M[1], M[2]);
    usum += IntegrateEdgeVec(M[2], M[2]);
    usum = -usum;

    Float denom = usum.Length();
    if (denom == 0)
        return 0;

    return nom / denom;
}

Vector3f LTC::IntegrateEdgeVec(const Vector3f &p1, const Vector3f &p2)
{
    Float x = Clamp(Dot(p1, p2), -1, 1);
    Float y = std::abs(x);

    Float a = 0.8543985f + (0.4965155f + 0.0145206f * y) * y;
    Float b = 3.4175940f + (4.1616724f + y) * y;
    Float v = a / b;

    Float theta_sintheta = (x > 0) ? v : 0.5f / std::sqrt(std::max(1.0f - x * x, (Float)1e-7)) - v;

    return Cross(p1, p2) * theta_sintheta;
}

Vector3f LTC::IntegrateEdgeVectTest(const Vector3f &p1, const Vector3f &p2)
{
    Float cosTheta = Clamp(Dot(p1, p2), -1, 1);
    Float theta = acos(cosTheta);
    Vector3f res = Cross(p1, p2) * ((theta > 0.001) ? theta / sin(theta) : 1.0f);
    return res;
}

static PBRT_CONSTEXPR float pivotTable[64][64][4] = {
#include "pivot_data.inl"
};

Float LTC::EvaluatePivot(const SurfaceInteraction &isect, const BRDFRecord &record, const BoundingSphere &boundingSphere)
{
    ProfilePhase _(Prof::EvaluateBRDFIntegral);

    Vector3f center = record.world2Local(boundingSphere.c - isect.p);
    // Float tmp = radius * radius / center.LengthSquared();
    // SphericalCap cap(Normalize(center), 0);
    // if (tmp > 1.0f)
    //{
    //    cap.cosTheta = -0.9999f;
    //}
    // else
    //{
    //    cap.cosTheta = std::sqrt(1.0f - tmp);
    //}
    Float lenSqr = center.LengthSquared();
    Float tmp = Clamp(boundingSphere.r * boundingSphere.r / lenSqr, 0.f, 1.f);
    SphericalCap cap(center / std::sqrt(lenSqr), std::sqrt(1.0f - tmp));
    Float denom = CapSolidAngle(cap, SphericalCap(Vector3f(0, 0, 1), 0.f));
    if (denom == 0)
        return 0;

    SphericalCap pCap = Cap2PCap(cap, record.specPivot);
    Float solidAngle = CapSolidAngle(pCap, record.specVisCap) * 0.079577472f;
    solidAngle = Clamp(solidAngle, 0.f, 1.f);
    Float ret = record.specBrdfScale * record.Ks * solidAngle / denom;

    return ret;
}

Float LTC::EvaluatePivotIntegral(const SurfaceInteraction &isect, const BRDFRecord &record, const BoundingSphere &boundingSphere)
{
    ProfilePhase _(Prof::EvaluateBRDFIntegral);

    Vector3f center = record.world2Local(boundingSphere.c - isect.p);
    // Float tmp = radius * radius / center.LengthSquared();
    // SphericalCap cap(Normalize(center), 0);
    // if (tmp > 1.0f)
    //{
    //    cap.cosTheta = -0.9999f;
    //}
    // else
    //{
    //    cap.cosTheta = std::sqrt(1.0f - tmp);
    //}
    Float lenSqr = center.LengthSquared();
    Float tmp = Clamp(boundingSphere.r * boundingSphere.r / lenSqr, 0.f, 1.f);
    SphericalCap cap(center / std::sqrt(lenSqr), std::sqrt(1.0f - tmp));
    Float denom = CapSolidAngle(cap, SphericalCap(Vector3f(0, 0, 1), 0.f));
    if (denom == 0)
        return 0;

    SphericalCap pCap = Cap2PCap(cap, record.specPivot);
    Float solidAngle = CapSolidAngle(pCap, record.specVisCap) * 0.079577472f;
    solidAngle = Clamp(solidAngle, 0.f, 1.f);
    Float ret = record.specBrdfScale * solidAngle;
    return ret;
}

Vector3f LTC::ExtractPivot(const Vector3f &wo, Float alpha, Float &brdfScale)
{
    Float theta = std::acos(std::min(wo.z, 1.0f));
    int u, v;
    u = Clamp((int)(2 * theta * InvPi * LUT_SCALE + LUT_BIAS), 0, 63);
    v = Clamp((int)(std::sqrt(alpha) * LUT_SCALE + LUT_BIAS), 0, 63);
    const auto &params = pivotTable[u][v];
    Vector3f pivot = params[0] * Vector3f(std::sin(params[1]), 0.f, std::cos(params[1]));
    brdfScale = params[3];

    return pivot;
}

Vector3f LTC::MaxValueVec(const SphericalCap &cap, const Vector3f &pivot)
{
    Float pivotLength = pivot.Length();
    if (pivotLength < 0.001f)
        return cap.dir;

    Vector3f pivotDir = pivot / pivotLength;

    Float cosPhi = Clamp(Dot(cap.dir, pivotDir), -1, 1);
    if (cosPhi >= cap.cosTheta)
        return pivotDir;

    // Project to 2D
    Float sinPhi = std::sqrt(1.0f - cosPhi * cosPhi);
    Vector3f pivotOrthoDir;
    if (std::abs(cosPhi) < 0.9999)
        pivotOrthoDir = (cap.dir - cosPhi * pivotDir) / sinPhi;
    Float sinTheta = std::sqrt(1.0f - cap.cosTheta * cap.cosTheta);
    Float cosPhiCosTheta = cosPhi * cap.cosTheta;
    Float cosPhiSinTheta = cosPhi * sinTheta;
    Float sinPhiCosTheta = sinPhi * cap.cosTheta;
    Float sinPhiSinTheta = sinPhi * sinTheta;
    Vector2f dir1(cosPhiCosTheta + sinPhiSinTheta, sinPhiCosTheta - cosPhiSinTheta);
    Vector2f dir2(cosPhiCosTheta - sinPhiSinTheta, sinPhiCosTheta + cosPhiSinTheta);
    Vector3f v1 = dir1.x * pivotDir + dir1.y * pivotOrthoDir;
    Vector3f v2 = dir2.x * pivotDir + dir2.y * pivotOrthoDir;

    if (dir1.x > dir2.x && v1.z > 0)
        return v1;
    else if (dir2.x > dir1.x && v2.z > 0)
        return v2;
    else
    {
        Vector3f capAxis = Normalize(Vector3f(cap.dir.x, cap.dir.y, 0));
        Float sinAlpha = cap.dir.z;
        Float cosAlpha = std::sqrt(1.0f - sinAlpha * sinAlpha);
        Vector2f topVec2d(cosAlpha * cap.cosTheta - sinAlpha * sinTheta, sinAlpha * cap.cosTheta + cosAlpha * sinTheta);
        Vector3f topVec = topVec2d.x * capAxis + topVec2d.y * Vector3f(0, 0, 1);

        if (v1.z > 0)
        {
            if (Dot(v1, pivotDir) > Dot(topVec, pivotDir))
                return v1;
            else
                return topVec;
        }
        else if (v2.z > 0)
        {
            if (Dot(v2, pivotDir) > Dot(topVec, pivotDir))
                return v2;
            else
                return topVec;
        }
        else
            return topVec;
    }
}

SphericalCap LTC::Cap2PCap(const SphericalCap &cap, const Vector3f &pivot)
{
    Float pivotLength = pivot.Length();
    if (pivotLength < 0.001f)
        return SphericalCap(-cap.dir, cap.cosTheta);
    Vector3f pivotDir = pivot / pivotLength;

    // Project to 2D
    Float cosPhi = Clamp(Dot(cap.dir, pivotDir), -1, 1);
    Float sinPhi = std::sqrt(1.0f - cosPhi * cosPhi);
    Vector3f pivotOrthoDir;
    if (std::abs(cosPhi) < 0.9999)
        pivotOrthoDir = (cap.dir - cosPhi * pivotDir) / sinPhi;
    Float sinTheta = std::sqrt(1.0f - cap.cosTheta * cap.cosTheta);
    Float cosPhiCosTheta = cosPhi * cap.cosTheta;
    Float cosPhiSinTheta = cosPhi * sinTheta;
    Float sinPhiCosTheta = sinPhi * cap.cosTheta;
    Float sinPhiSinTheta = sinPhi * sinTheta;
    Vector2f dir1(cosPhiCosTheta + sinPhiSinTheta, sinPhiCosTheta - cosPhiSinTheta);
    Vector2f dir2(cosPhiCosTheta - sinPhiSinTheta, sinPhiCosTheta + cosPhiSinTheta);
    // Pivot transform
    Vector2f pDir1 = Vec2PVec(dir1, pivotLength);
    Vector2f pDir2 = Vec2PVec(dir2, pivotLength);
    // Transformed cap 2D direction
    Vector2f pDir = pDir1 + pDir2;
    Float lenSqr = pDir.LengthSquared();
    if (lenSqr < MachineEpsilon)
        pDir = Vector2f(-pDir1.y, pDir1.x);
    else
        pDir = pDir / std::sqrt(lenSqr);
    Float sinDelta = pDir2.y * pDir1.x - pDir2.x * pDir1.y;
    if (sinDelta < 0)
        pDir = -pDir;
    // Come back to 3D
    Vector3f capDir = pDir.x * pivotDir + pDir.y * pivotOrthoDir;
    Float capCos = Clamp(Dot(pDir, pDir1), -1, 1);

    return SphericalCap(capDir, capCos);
}

Vector2f LTC::Vec2PVec(const Vector2f &r, Float rp)
{
    Vector2f temp1 = Vector2f(r.x - rp, r.y);
    Vector2f temp2 = rp * r - Vector2f(1.f, 0.f);
    Float x = Dot(temp1, temp2);
    Float y = temp1.y * temp2.x - temp1.x * temp2.y;
    return Vector2f(x, y) / temp2.LengthSquared();
}

Vector3f LTC::Vec2PVec(const Vector3f &r, const Vector3f &rp)
{
    Vector3f tmp = r - rp;
    Vector3f cp1 = Cross(r, rp);
    Vector3f cp2 = Cross(tmp, cp1);
    Float dp = Dot(r, rp) - 1.f;
    Float qf = dp * dp + Dot(cp1, cp1);

    return ((dp * tmp - cp2) / qf);
}

Vector3f LTC::U2Cap(const Point2f &u, const SphericalCap &cap)
{
    // Generate the sample in the basis aligned with the cap
    Float z = (1.0f - cap.cosTheta) * u.x + cap.cosTheta;
    Float sinTheta = std::sqrt(1.0f - z * z);
    Float phi = 2 * Pi * u.y;
    Float x = sinTheta * std::cos(phi);
    Float y = sinTheta * std::sin(phi);

    // Compute basis vectors
    Vector3f t1, t2;
    if (cap.dir.z < -0.9999999f)
    {
        t1 = Vector3f(0.f, -1.0f, 0.f);
        t2 = Vector3f(-1.0f, 0.f, 0.f);
    }
    else
    {
        const Float a = 1.0f / (1.0f + cap.dir.z);
        const Float b = -cap.dir.x * cap.dir.y * a;
        t1 = Vector3f(1.0f - cap.dir.x * cap.dir.x * a, b, -cap.dir.x);
        t2 = Vector3f(b, 1.0f - cap.dir.y * cap.dir.y * a, -cap.dir.y);
    }

    Matrix4x4 mat(t1[0], t1[1], t1[2], 0.f, t2[0], t2[1], t2[2], 0.f, cap.dir[0], cap.dir[1], cap.dir[2], 0.f, 0.f, 0.f, 0.f, 1.f);
    Transform trans(mat);

    // Warp the sample in the proper basis
    return Normalize(trans(Vector3f(x, y, z)));
}

Float LTC::PdfCap(const Vector3f &wk, const SphericalCap &cap)
{
    // Make sure the sample lies in the cap
    if (Dot(wk, cap.dir) >= cap.cosTheta)
        return 1.0f / CapSolidAngle(cap);
    return 0.0f;
}

Float LTC::PivotJacobian(const Vector3f &wk, const Vector3f &rp)
{
    Float num = 1.0f - Dot(rp, rp);
    Vector3f tmp = wk - rp;
    Float den = Dot(tmp, tmp);

    return (num * num) / (den * den);
}

Float LTC::CapSolidAngle(const SphericalCap &cap)
{
    return 2 * Pi * (1.f - cap.cosTheta);
}

// Based on Oat and Sander's 2008 technique
Float LTC::CapSolidAngle(const SphericalCap &cap1, const SphericalCap &cap2)
{
    Float r1 = std::acos(cap1.cosTheta);
    Float r2 = std::acos(cap2.cosTheta);
    Float rd = std::acos(Clamp(Dot(cap1.dir, cap2.dir), -1, 1));

    if (rd <= std::max(r1, r2) - std::min(r1, r2))
    {
        // One cap in completely inside the other
        return 2 * Pi * (1.f - std::max(cap1.cosTheta, cap2.cosTheta));
    }
    else if (rd >= r1 + r2)
    {
        // No intersection
        return 0.f;
    }
    else
    {
        Float fDiff = std::abs(r1 - r2);
        Float den = r1 + r2 - fDiff;
        Float x = 1.0f - Clamp((rd - fDiff) / den, 0.f, 1.f);
        x = x * x * (3.f - 2.f * x);
        return x * 2 * Pi * (1.f - std::max(cap1.cosTheta, cap2.cosTheta));
    }
}

int LTC::ClipQuadToHorizon(Vector3f L[5])
{
    int n;
    // detect clipping config
    int config = 0;
    if (L[0].z > 0.0)
        config += 1;
    if (L[1].z > 0.0)
        config += 2;
    if (L[2].z > 0.0)
        config += 4;
    if (L[3].z > 0.0)
        config += 8;

    // clip
    n = 0;

    if (config == 0)
    {
        // clip all
    }
    else if (config == 1)  // V1 clip V2 V3 V4
    {
        n = 3;
        L[1] = -L[1].z * L[0] + L[0].z * L[1];
        L[2] = -L[3].z * L[0] + L[0].z * L[3];
    }
    else if (config == 2)  // V2 clip V1 V3 V4
    {
        n = 3;
        L[0] = -L[0].z * L[1] + L[1].z * L[0];
        L[2] = -L[2].z * L[1] + L[1].z * L[2];
    }
    else if (config == 3)  // V1 V2 clip V3 V4
    {
        n = 4;
        L[2] = -L[2].z * L[1] + L[1].z * L[2];
        L[3] = -L[3].z * L[0] + L[0].z * L[3];
    }
    else if (config == 4)  // V3 clip V1 V2 V4
    {
        n = 3;
        L[0] = -L[3].z * L[2] + L[2].z * L[3];
        L[1] = -L[1].z * L[2] + L[2].z * L[1];
    }
    else if (config == 5)  // V1 V3 clip V2 V4) impossible
    {
        n = 0;
    }
    else if (config == 6)  // V2 V3 clip V1 V4
    {
        n = 4;
        L[0] = -L[0].z * L[1] + L[1].z * L[0];
        L[3] = -L[3].z * L[2] + L[2].z * L[3];
    }
    else if (config == 7)  // V1 V2 V3 clip V4
    {
        n = 5;
        L[4] = -L[3].z * L[0] + L[0].z * L[3];
        L[3] = -L[3].z * L[2] + L[2].z * L[3];
    }
    else if (config == 8)  // V4 clip V1 V2 V3
    {
        n = 3;
        L[0] = -L[0].z * L[3] + L[3].z * L[0];
        L[1] = -L[2].z * L[3] + L[3].z * L[2];
        L[2] = L[3];
    }
    else if (config == 9)  // V1 V4 clip V2 V3
    {
        n = 4;
        L[1] = -L[1].z * L[0] + L[0].z * L[1];
        L[2] = -L[2].z * L[3] + L[3].z * L[2];
    }
    else if (config == 10)  // V2 V4 clip V1 V3) impossible
    {
        n = 0;
    }
    else if (config == 11)  // V1 V2 V4 clip V3
    {
        n = 5;
        L[4] = L[3];
        L[3] = -L[2].z * L[3] + L[3].z * L[2];
        L[2] = -L[2].z * L[1] + L[1].z * L[2];
    }
    else if (config == 12)  // V3 V4 clip V1 V2
    {
        n = 4;
        L[1] = -L[1].z * L[2] + L[2].z * L[1];
        L[0] = -L[0].z * L[3] + L[3].z * L[0];
    }
    else if (config == 13)  // V1 V3 V4 clip V2
    {
        n = 5;
        L[4] = L[3];
        L[3] = L[2];
        L[2] = -L[1].z * L[2] + L[2].z * L[1];
        L[1] = -L[1].z * L[0] + L[0].z * L[1];
    }
    else if (config == 14)  // V2 V3 V4 clip V1
    {
        n = 5;
        L[4] = -L[0].z * L[3] + L[3].z * L[0];
        L[0] = -L[0].z * L[1] + L[1].z * L[0];
    }
    else if (config == 15)  // V1 V2 V3 V4
    {
        n = 4;
    }

    if (n == 3)
        L[3] = L[0];
    if (n == 4)
        L[4] = L[0];

    return n;
}

}  // namespace pbrt