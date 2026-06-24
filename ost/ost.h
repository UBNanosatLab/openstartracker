#ifndef OST_H
#define OST_H

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef OST_RESTRICT
#ifdef __cplusplus
#define OST_RESTRICT __restrict__
#else
#define OST_RESTRICT restrict
#endif
#endif

#ifndef OST_UNUSED
#if defined(__GNUC__) || defined(__clang__)
#define OST_UNUSED __attribute__((unused))
#else
#define OST_UNUSED
#endif
#endif

#ifndef OST_DEF
#ifdef OST_EXPORT
#ifdef _WIN32
#define OST_DEF __declspec(dllexport)
#else
#define OST_DEF __attribute__((visibility("default")))
#endif
#else
#define OST_DEF
#endif
#endif

#define QMETHOD_ITER 0

#ifdef __cplusplus
extern "C" {
#endif

/* Single-header C core; all algorithm work buffers are caller-owned. */

#define OST_SWAP(type, a, b)              \
    do {                                  \
        type ost_swap_tmp = (a);          \
        (a)               = (b);          \
        (b)               = ost_swap_tmp; \
    } while (0)

#define OST_KD_COORD(base, elem_size, key_offset, idx, dim) \
    (*(float *)((char *)(base) + (size_t)(idx) * (elem_size) + \
                (key_offset) + (size_t)(dim) * sizeof(float)))

typedef int (*OSTKDVisit)(void *base, int idx, void *ctx);

static inline void ost_kdselect(void *base, size_t elem_size, size_t key_offset,
                                int l, int r, int k, int dim)
{
    r--;
    while (l < r) {
        float p = OST_KD_COORD(base, elem_size, key_offset, (l + r) >> 1, dim);
        int i = l, j = r;
        while (i <= j) {
            while (OST_KD_COORD(base, elem_size, key_offset, i, dim) < p) i++;
            while (OST_KD_COORD(base, elem_size, key_offset, j, dim) > p) j--;
            if (i <= j) {
                char *a = (char *)base + (size_t)i * elem_size;
                char *b = (char *)base + (size_t)j * elem_size;
                for (size_t n = 0; n < elem_size; n++) {
                    char t = a[n];
                    a[n] = b[n];
                    b[n] = t;
                }
                i++;
                j--;
            }
        }
        if (k <= j) r = j;
        else if (k >= i) l = i;
        else break;
    }
}

static inline void ost_kdbuild(void *base, size_t elem_size, size_t key_offset,
                               int min, int max, int bucket, int dim, int dims)
{
    int mid;
    if (max - min <= bucket || dims <= 0)
        return;
    mid = (min + max) / 2;
    ost_kdselect(base, elem_size, key_offset, min, max, mid, dim);
    ost_kdbuild(base, elem_size, key_offset, min, mid, bucket,
                (dim + 1) % dims, dims);
    ost_kdbuild(base, elem_size, key_offset, mid + 1, max, bucket,
                (dim + 1) % dims, dims);
}

static inline int ost_kd_in_box(void *base, size_t elem_size,
                                size_t key_offset, int idx, int dims,
                                const float *p, const float *r)
{
    const float *x = (const float *)((const unsigned char *)base +
                     (size_t)idx * elem_size + key_offset);
    for (int d = 0; d < dims; d++)
        if (fabsf(x[d] - p[d]) > r[d])
            return 0;
    return 1;
}

static inline int ost_kdsearch(void *base, size_t elem_size, size_t key_offset,
                               int min, int max, int bucket, int dim, int dims,
                               const float *p, const float *r,
                               OSTKDVisit visit, void *ctx)
{
    int mid, rc;
    if (max <= min || dims <= 0)
        return 0;
    if (max - min <= bucket) {
        for (int i = min; i < max; i++) {
            if (!ost_kd_in_box(base, elem_size, key_offset, i, dims, p, r))
                continue;
            rc = visit(base, i, ctx);
            if (rc)
                return rc;
        }
        return 0;
    }
    mid = (min + max) / 2;
    if (min < mid && p[dim] - r[dim] <=
        OST_KD_COORD(base, elem_size, key_offset, mid, dim)) {
        rc = ost_kdsearch(base, elem_size, key_offset, min, mid, bucket,
                          (dim + 1) % dims, dims, p, r, visit, ctx);
        if (rc < 0)
            return rc;
    }
    if (ost_kd_in_box(base, elem_size, key_offset, mid, dims, p, r)) {
        rc = visit(base, mid, ctx);
        if (rc)
            return rc;
    }
    if (mid + 1 < max && OST_KD_COORD(base, elem_size, key_offset, mid, dim)
        <= p[dim] + r[dim])
        return ost_kdsearch(base, elem_size, key_offset, mid + 1, max, bucket,
                            (dim + 1) % dims, dims, p, r, visit, ctx);
    return 0;
}

#ifndef OST_MIN_COMPONENT_AREA
#define OST_MIN_COMPONENT_AREA 4
#endif

/* Binary fields are area/sums; weighted fields are flux moments. */
typedef struct OSTCCComponent {
    int area;
    int sum_x;
    int sum_y;
    int min_x;
    int max_x;
    int min_y;
    int max_y;
    int signal;
    double wsum;
    double wx;
    double wy;
    double wxx;
    double wyy;
    double wxy;
    double eig_min;
} OSTCCComponent;

/* Workspace is width-bounded: labels are recycled when components close. */
typedef struct OSTCCBufferSizes {
    int max_labels;
    size_t components;
    size_t parent;
    size_t col_label;
    size_t active_count;
    size_t reuse_after_row;
    size_t total_bytes;
} OSTCCBufferSizes;

typedef struct OSTCCContext {
    int width;
    int max_labels;
    OSTCCComponent *components;
    int *parent;
    int *col_label;
    int *active_count;
    int *reuse_after_row;
} OSTCCContext;

int ost_cc_buffer_sizes(int width, OSTCCBufferSizes *sizes);
int ost_cc_init(OSTCCContext *ctx, int width,
                OSTCCComponent *components, int *parent,
                int *col_label, int *active_count,
                int *reuse_after_row);
int ost_cc_threshold_4(const unsigned char *image,
                       int width, int height, int stride,
                       unsigned char threshold,
                       OSTCCComponent *out, int out_max,
                       OSTCCContext *ctx);

/* Tiled background maps feed both local thresholds and photometric weights. */
typedef struct OSTBGConfig {
    int width;
    int height;
    int tile_size;
    int map_width;
    int map_height;
    int max_stars;
    int max_pixel_brightness;
    int sample_radius;
    double psf_sigma;
    double threshold_sigma;
    double detect_sigma;
} OSTBGConfig;

typedef struct OSTBGStats {
    const OSTBGConfig *cfg;
    double *mean;
    double *var;
    double *poisson;
    int *x_edge;
    int *y_edge;
} OSTBGStats;

/* Optional PSF-like refinement uses caller-owned least-squares workspace. */
typedef struct OSTBGFitStar {
    int xi;
    int yi;
} OSTBGFitStar;

typedef struct OSTBGFitWorkspace {
    OSTBGFitStar *stars1;
    OSTBGFitStar *stars2;
    double *params1;
    double *params2;
    double *normal;
    double *rhs;
    double *rhs_solve;
    double *cov_xy;
    double *dropped;
} OSTBGFitWorkspace;

int ost_png_dimensions(const char *filename, int *width, int *height);
int ost_png_read_rgba(const char *filename, unsigned char *dst,
                      int width, int height, int stride);
int ost_bg_config_init(OSTBGConfig *cfg, int width, int height);
void ost_bg_rgba_to_gray16(uint16_t *dst, const unsigned char *rgba,
                           int width, int height);
int ost_bg_compute_stats(const OSTBGConfig *cfg, const uint16_t *image,
                         int stride, double *mean, double *var,
                         double *poisson, int *x_edge, int *y_edge,
                         int *hist, int hist_len);
int ost_bg_extract_fused(const OSTBGConfig *cfg, const uint16_t *image,
                         int stride, OSTBGStats *stats, OSTCCContext *cc,
                         OSTCCComponent *stars, int stars_max);
int ost_bg_fit_stars(const OSTBGConfig *cfg, const uint16_t *image, int stride,
                     const OSTBGStats *stats, const OSTCCComponent *components,
                     int component_count, int num_iter,
                     OSTBGFitWorkspace *work, int max_stars,
                     double *params_out, double *cov_xy_out,
                     int *dropped_count_out);

static const double PI = 3.14159265358979323846;

typedef float Vec3[3];
typedef Vec3 Mat3[3];

/* Star identification and attitude primitives. */

/* File config plus derived field-of-view and tangent-plane calibration. */
typedef struct Config {
    int IMG_X, IMG_Y, MAX_FALSE_STARS, DB_REDUNDANCY, REQUIRED_STARS;
    int KDBUCKET_SIZE;
    float PIXSCALE, DOUBLE_STAR_PX, BASE_FLUX, IMAGE_VARIANCE;
    float THRESH_FACTOR, POS_VARIANCE, POS_ERR_SIGMA, PSF_SIGMA;
    float MAXFOV, MINFOV, MATCH_VALUE, PIXX_TANGENT, PIXY_TANGENT;
} Config;

/* Catalog and image stars share unit-vector storage; image stars also keep px/py. */
typedef struct Star {
    union {
        Vec3 v;
        struct { float x, y, z; };
    };
    float flux;
    float px, py, sigma_sq;
    int id, star_idx, unreliable;
} Star;

typedef struct StarDB {
    Star *v;
    int n, cap;
    float max_variance;
} StarDB;

/* KD-search map plus caller-owned result and mask arrays. */
typedef struct Query {
    Star *map;
    int n, kdsorted;
    int *kdresults;
    int kdresults_size, kdresults_maxsize;
    signed char *kdmask;
} Query;

/* Pairwise angular separation; sorted arrays enable candidate range lookup. */
typedef struct Constellation {
    float p;
    int s1, s2, idx;
} Constellation;

typedef struct CPair {
    float totalscore;
    int db_s1, db_s2, img_s1, img_s2;
} CPair;

typedef struct CDB {
    StarDB stars;
    Query results;
    Constellation *map;
    int map_size;
} CDB;

typedef struct StarFov {
    int *mask, *collision;
    int collision_size, collision_cap;
    float *s_px, *s_py, maxdist_sq, sigma_sq;
} StarFov;

typedef struct MatchResult {
    CPair match;
    Mat3 R;
    int *map;
    int map_size;
    CDB *db, *img;
    StarFov *img_mask;
} MatchResult;

typedef struct MatchWork {
    CPair *candidates;
    int candidate_cap;
    int *fov_mask, *collision;
    int collision_cap;
    float *fov_px, *fov_py, *scores;
    int *match_map, *work_map;
} MatchWork;

#define OST_MAX_CONSTELLATION_STARS 5
#define OST_MAX_CONSTELLATION_DIMS \
    (OST_MAX_CONSTELLATION_STARS * (OST_MAX_CONSTELLATION_STARS - 1) / 2)
#define OST_MAX_CONSTELLATION_PERMS 120

typedef struct OSTConstellationEdge {
    int star;
} OSTConstellationEdge;

typedef enum OSTConstellationDescriptorKind {
    OST_CONSTELLATION_PAIRDIST = 0,
    OST_CONSTELLATION_CROSSRATIO = 1
} OSTConstellationDescriptorKind;

typedef struct OSTConstellationIndex {
    CDB *pair;
    unsigned char *map;
    int map_size, cap, kd_ready, kd_bucket;
    int k, dims, descriptor_kind;
    size_t record_size;
} OSTConstellationIndex;

static inline float ost_vec_dot3(const float a[3], const float b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline void ost_vec_cross3(float r[3], const float a[3], const float b[3])
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

static inline int ost_vec_normalize3(float v[3])
{
    float n = sqrtf(ost_vec_dot3(v, v));
    if (n <= 1e-20f)
        return -1;
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;
    return 0;
}

static inline float ost_star_chord(const Star *a, const Star *b)
{
    float dx = a->v[0] - b->v[0], dy = a->v[1] - b->v[1];
    float dz = a->v[2] - b->v[2];
    return sqrtf(dx * dx + dy * dy + dz * dz);
}

static inline float ost_star_dist_arcsec(const Star *a, const Star *b)
{
    float x = a->v[0] * b->v[1] - b->v[0] * a->v[1];
    float y = a->v[0] * b->v[2] - b->v[0] * a->v[2];
    float z = a->v[1] * b->v[2] - b->v[1] * a->v[2];
    return (float)((3600.0 * 180.0 / PI) *
                   asinf(sqrtf(x * x + y * y + z * z)));
}

static inline float ost_wahba_det3(const float *a, const float *b,
                                          const float *c,
                                          int i, int j, int k)
{
    return a[i] * (b[j] * c[k] - b[k] * c[j]) -
           a[j] * (b[i] * c[k] - b[k] * c[i]) +
           a[k] * (b[i] * c[j] - b[j] * c[i]);
}

static inline int ost_wahba_null4(float q[4], const float a[4][4])
{
    float best[4] = {1.0f, 0.0f, 0.0f, 0.0f};
    float bestn = -1.0f;
    for (int skip = 0; skip < 4; skip++) {
        const float *r[3];
        float v[4], n;
        int m = 0;
        for (int i = 0; i < 4; i++)
            if (i != skip)
                r[m++] = a[i];
        v[0] =  ost_wahba_det3(r[0], r[1], r[2], 1, 2, 3);
        v[1] = -ost_wahba_det3(r[0], r[1], r[2], 0, 2, 3);
        v[2] =  ost_wahba_det3(r[0], r[1], r[2], 0, 1, 3);
        v[3] = -ost_wahba_det3(r[0], r[1], r[2], 0, 1, 2);
        n = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
        if (n > bestn) {
            bestn = n;
            memcpy(best, v, sizeof(best));
        }
    }
    if (bestn <= 1e-30f)
        return -1;
    bestn = 1.0f / sqrtf(bestn);
    for (int i = 0; i < 4; i++)
        q[i] = best[i] * bestn;
    return 0;
}

static inline void ost_wahba_quat_to_mat(Mat3 r, const float qin[4])
{
    float q0 = qin[0], q1 = qin[1], q2 = qin[2], q3 = qin[3];
    float n = 1.0f / sqrtf(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= n; q1 *= n; q2 *= n; q3 *= n;
    r[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    r[0][1] = 2.0f * (q1*q2 - q0*q3);
    r[0][2] = 2.0f * (q1*q3 + q0*q2);
    r[1][0] = 2.0f * (q1*q2 + q0*q3);
    r[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    r[1][2] = 2.0f * (q2*q3 - q0*q1);
    r[2][0] = 2.0f * (q1*q3 - q0*q2);
    r[2][1] = 2.0f * (q2*q3 + q0*q1);
    r[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}

static inline int weighted_wahba_vectors(
        Mat3 r, const Star *db_stars, const Star *img_stars,
        const int *db_idx, const int *img_idx, int n, int iter)
{
    float b[3][3] = {{0.0f}}, k[4][4] = {{0.0f}}, a[4][4], q[4];
    float sw = 0.0f, tr, lambda;
    if (!r || !db_stars || !img_stars || n < 2)
        return -1;
    if (iter < 0)
        iter = 0;
    if (iter > 8)
        iter = 8;
    for (int i = 0; i < n; i++) {
        const Star *ds = &db_stars[db_idx ? db_idx[i] : i];
        const Star *is = &img_stars[img_idx ? img_idx[i] : i];
        float w = 1.0f / (ds->sigma_sq + is->sigma_sq);
        sw += w;
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
                b[row][col] += w * ds->v[row] * is->v[col];
    }
    if (sw <= 0.0f)
        return -1;
    tr = b[0][0] + b[1][1] + b[2][2];
    k[0][0] = tr;
    k[0][1] = k[1][0] = b[1][2] - b[2][1];
    k[0][2] = k[2][0] = b[2][0] - b[0][2];
    k[0][3] = k[3][0] = b[0][1] - b[1][0];
    k[1][1] = b[0][0] - b[1][1] - b[2][2];
    k[1][2] = k[2][1] = b[0][1] + b[1][0];
    k[1][3] = k[3][1] = b[0][2] + b[2][0];
    k[2][2] = -b[0][0] + b[1][1] - b[2][2];
    k[2][3] = k[3][2] = b[1][2] + b[2][1];
    k[3][3] = -b[0][0] - b[1][1] + b[2][2];

    lambda = sw;
    for (int it = 0; it <= iter; it++) {
        for (int row = 0; row < 4; row++)
            for (int col = 0; col < 4; col++)
                a[row][col] = k[row][col] - (row == col ? lambda : 0.0f);
        if (ost_wahba_null4(q, a) < 0)
            return -1;
        lambda = 0.0f;
        for (int row = 0; row < 4; row++) {
            float kq = 0.0f;
            for (int col = 0; col < 4; col++)
                kq += k[row][col] * q[col];
            lambda += q[row] * kq;
        }
    }
    ost_wahba_quat_to_mat(r, q);
    return 0;
}

int ost_load_config(Config *c, const char *filename);
Star ost_make_db_star(const Config *c, float x, float y, float z, float flux, int id);
Star ost_make_img_star(const Config *c, float px, float py, float flux, int id);
void ost_star_db_init(StarDB *db, Star *storage, int cap);
int ost_db_add(StarDB *db, Star s);
int ost_load_catalog(const Config *cfg, StarDB *db, const char *filename,
                     float year, uint64_t *cat_keys, size_t cat_key_cap);
void ost_query_init(Query *q, StarDB *db, Star *map, int *res, signed char *mask);
void ost_query_reset_mask(Query *q);
void ost_query_clear_results(Query *q);
void ost_query_sort_flux(Query *q);
void ost_query_kdsort(Query *q, const Config *c);
void ost_query_search(Query *q, const Config *c, const float p[3],
                      float radius, float min_flux);
void ost_query_search_range(Query *q, const Config *c, const float p[3],
                            float r, float min_flux, int start, int end, int dim);
void ost_query_mask_filter(Query *q, const Config *c);
void ost_query_mask_uniform(Query *q, const Config *c, int min_stars, signed char *keep);
int ost_db_from_mask(StarDB *out, Query *q);
int ost_db_from_results(StarDB *out, Query *q);
int ost_db_from_image(CDB *cdb, StarDB *src, Star *star_storage, int star_cap,
                      Query *q, Star *query_map, int *query_results,
                      signed char *query_mask, Constellation *map, int map_cap,
                      int stars_per_fov);
int ost_db_from_catalog(CDB *cdb, StarDB *src, Star *star_storage, int star_cap,
                        Query *q, Star *query_map, int *query_results,
                        signed char *query_mask, Constellation *map, int map_cap,
                        int stars_per_fov, const Config *cfg, signed char *keep);
size_t ost_constellation_record_size(int k, int descriptor_kind);
int ost_constellation_index_init(OSTConstellationIndex *idx, CDB *pair,
                                 int k, int descriptor_kind,
                                 unsigned char *storage, int storage_cap);
int ost_constellation_count(CDB *pair, int k, uint64_t *count,
                            int *off, OSTConstellationEdge *edges,
                            int *tmp, int *common);
int ost_constellation_index_build(OSTConstellationIndex *idx, int *off,
                                  OSTConstellationEdge *edges, int *tmp,
                                  int *common);
void ost_constellation_index_kdsort(OSTConstellationIndex *idx);
int ost_db_match_constellations(OSTConstellationIndex *cat, CDB *img,
                                MatchResult *winner, const Config *cfg,
                                MatchWork *mw, float *p_match);
void ost_match_work_init(MatchWork *mw, CPair *candidates, int candidate_cap,
                         int *fov_mask, int *collision, int collision_cap,
                         float *fov_px, float *fov_py, float *scores,
                         int *match_map, int *work_map);
int fov_init(StarFov *fov, StarDB *stars, float db_max_variance,
             const Config *c, MatchWork *w);
float fov_score(StarFov *fov, int id, float px, float py);
int fov_resolve(StarFov *fov, int id, float px, float py);
int fov_get_id(StarFov *fov, const Config *c, float px, float py);
void mr_init(MatchResult *m, CDB *db, CDB *img, StarFov *mask, int *map);
void mr_copy(MatchResult *dst, MatchResult *src);
void compute_score(MatchResult *m, const Config *c, MatchWork *w);
int related(MatchResult *winner, CPair *p);
int ost_chol(const double *a, int n, double *l);
void ost_chol_solve(const double *l, int n, const double *b,
                    double *x, double *y);
void ost_chol_inv_diag(const double *l, int n, double *d,
                       double *e, double *x, double *y);
int ost_copy_n_brightest(StarDB *dst, StarDB *src, Star *tmp, int n);

#ifdef __cplusplus
}
#endif

#ifdef OST_IMPLEMENTATION

/* utilities */
typedef double OSTSym3[6]; /* packed lower: 00,10,11,20,21,22 */
typedef double OSTDVec3[3];

static inline int ost_tri_idx(int r, int c) { return r * (r + 1) / 2 + c; }
OST_DEF int ost_chol(const double *a, int n, double *l)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= i; j++) {
            double d = a[ost_tri_idx(i, j)];
            for (int k = 0; k < j; k++)
                d -= l[ost_tri_idx(i, k)] * l[ost_tri_idx(j, k)];
            if (i == j) {
                if (d <= 0.0) return -1;
                l[ost_tri_idx(i, j)] = sqrt(d);
            } else {
                l[ost_tri_idx(i, j)] = d / l[ost_tri_idx(j, j)];
            }
        }
    return 0;
}
OST_DEF void ost_chol_solve(const double *l, int n, const double *b,
                            double *x, double *y)
{
    for (int i = 0; i < n; i++) {
        double s = b[i];
        for (int k = 0; k < i; k++)
            s -= l[ost_tri_idx(i, k)] * y[k];
        y[i] = s / l[ost_tri_idx(i, i)];
    }
    for (int i = n; i-- > 0;) {
        double s = y[i];
        for (int k = i + 1; k < n; k++)
            s -= l[ost_tri_idx(k, i)] * x[k];
        x[i] = s / l[ost_tri_idx(i, i)];
    }
}
OST_DEF void ost_chol_inv_diag(const double *l, int n, double *d,
                               double *e, double *x, double *y)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            e[j] = 0.0;
        e[i] = 1.0;
        ost_chol_solve(l, n, e, x, y);
        d[i] = x[i];
    }
}

/* connected components */

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define OST_CC_RESTRICT restrict
#else
#define OST_CC_RESTRICT
#endif

typedef double (*OSTCCBackgroundFn)(void *opaque, int x, int y, double *var);

static void clear_components(OSTCCComponent *p, int n)
{
    memset(p, 0, (size_t)n * sizeof(*p));
}

static void insert_component(OSTCCComponent *OST_CC_RESTRICT out,
                             int *count, int out_max,
                             int area, int sum_x, int sum_y)
{
    int i;

    if (out_max <= 0 || area <= 0)
        return;
    if (*count >= out_max && area <= out[out_max - 1].area)
        return;

    i = (*count < out_max) ? (*count)++ : out_max - 1;
    while (i > 0 && area > out[i - 1].area) {
        out[i].area = out[i - 1].area;
        out[i].sum_x = out[i - 1].sum_x;
        out[i].sum_y = out[i - 1].sum_y;
        i--;
    }
    out[i].area = area;
    out[i].sum_x = sum_x;
    out[i].sum_y = sum_y;
}

static void component_clear(OSTCCComponent *c)
{
    memset(c, 0, sizeof(*c));
}

static void component_merge(OSTCCComponent *keep, const OSTCCComponent *merge)
{
    if (merge->area <= 0)
        return;
    if (keep->area <= 0) {
        *keep = *merge;
        return;
    }
    keep->area += merge->area;
    keep->sum_x += merge->sum_x;
    keep->sum_y += merge->sum_y;
    if (merge->min_x < keep->min_x) keep->min_x = merge->min_x;
    if (merge->max_x > keep->max_x) keep->max_x = merge->max_x;
    if (merge->min_y < keep->min_y) keep->min_y = merge->min_y;
    if (merge->max_y > keep->max_y) keep->max_y = merge->max_y;
    keep->signal |= merge->signal;
    keep->wsum += merge->wsum;
    keep->wx += merge->wx;
    keep->wy += merge->wy;
    keep->wxx += merge->wxx;
    keep->wyy += merge->wyy;
    keep->wxy += merge->wxy;
}

static void component_add_binary_pixel(OSTCCComponent *c, int x, int y)
{
    if (c->area <= 0) {
        c->min_x = x;
        c->max_x = x;
        c->min_y = y;
        c->max_y = y;
    } else {
        if (x < c->min_x) c->min_x = x;
        if (x > c->max_x) c->max_x = x;
        if (y < c->min_y) c->min_y = y;
        if (y > c->max_y) c->max_y = y;
    }
    c->area++;
    c->sum_x += x;
    c->sum_y += y;
}

static void component_add_weighted_pixel(OSTCCComponent *c,
                                         const unsigned short *image,
                                         int stride, int x, int y,
                                         OSTCCBackgroundFn background,
                                         void *background_opaque,
                                         double signal_sigma)
{
    double var;
    double mu;
    double v;

    component_add_binary_pixel(c, x, y);
    mu = background(background_opaque, x, y, &var);
    v = (double)image[(size_t)y * (size_t)stride + (size_t)x] - mu;
    if (v > signal_sigma * sqrt(var))
        c->signal = 1;
    if (v <= 0.0)
        return;
    c->wsum += v;
    c->wx += v * x;
    c->wy += v * y;
    c->wxx += v * x * x;
    c->wyy += v * y * y;
    c->wxy += v * x * y;
}

static int component_finish(OSTCCComponent *c)
{
    /* Convert raw weighted sums into centroid moments and a shape statistic. */
    double cx;
    double cy;
    double u20;
    double u02;
    double u11;
    double tr;
    double det;
    double d;

    if (c->area < OST_MIN_COMPONENT_AREA || c->wsum <= 0 || !c->signal ||
        c->min_x >= c->max_x || c->min_y >= c->max_y)
        return 0;

    cx = c->wx / c->wsum;
    cy = c->wy / c->wsum;
    u20 = c->wxx / c->wsum - cx * cx;
    u02 = c->wyy / c->wsum - cy * cy;
    u11 = c->wxy / c->wsum - cx * cy;
    tr = u20 + u02;
    det = u20 * u02 - u11 * u11;
    d = tr * tr - 4.0 * det;
    c->eig_min = (tr - sqrt(d > 0 ? d : 0)) * 0.5;
    return 1;
}

static void insert_weighted_component(OSTCCComponent *OST_CC_RESTRICT out,
                                      int *count, int out_max,
                                      OSTCCComponent c)
{
    int i;

    if (out_max <= 0 || !component_finish(&c))
        return;
    if (*count >= out_max && c.wsum <= out[out_max - 1].wsum)
        return;

    i = (*count < out_max) ? (*count)++ : out_max - 1;
    out[i] = c;
    while (i > 0 && out[i].wsum > out[i - 1].wsum) {
        OST_SWAP(OSTCCComponent, out[i], out[i - 1]);
        i--;
    }
}

static int root_compress(int *OST_CC_RESTRICT parent, int label)
{
    int root;

    root = label;
    while (parent[root] != root)
        root = parent[root];

    while (label != root) {
        int next;
        next = parent[label];
        parent[label] = root;
        label = next;
    }
    return root;
}

static int alloc_label(OSTCCContext *ctx, int row, int *next_reuse_label)
{
    int start;
    int label;

    start = *next_reuse_label;
    label = start;
    do {
        if (ctx->active_count[label] == 0 && ctx->reuse_after_row[label] < row) {
            ctx->parent[label] = label;
            component_clear(&ctx->components[label]);

            label++;
            if (label >= ctx->max_labels)
                label = 1;
            *next_reuse_label = label;
            return label == 1 ? ctx->max_labels - 1 : label - 1;
        }

        label++;
        if (label >= ctx->max_labels)
            label = 1;
    } while (label != start);

    return 0;
}

static int merge_roots(OSTCCContext *ctx, int a, int b, int row)
{
    int ra;
    int rb;
    int keep;
    int merge;

    ra = root_compress(ctx->parent, a);
    rb = root_compress(ctx->parent, b);
    if (ra == rb)
        return ra;

    keep = (ra < rb) ? ra : rb;
    merge = (ra < rb) ? rb : ra;

    component_merge(&ctx->components[keep], &ctx->components[merge]);
    ctx->parent[merge] = keep;
    ctx->active_count[keep] += ctx->active_count[merge];
    ctx->active_count[merge] = 0;
    component_clear(&ctx->components[merge]);
    ctx->reuse_after_row[merge] = row + 1;

    return keep;
}

static void close_binary_column(OSTCCContext *ctx, int x, int row,
                                OSTCCComponent *out, int *count,
                                int out_max)
{
    int label;
    int root;

    label = ctx->col_label[x];
    if (!label)
        return;

    ctx->col_label[x] = 0;
    root = root_compress(ctx->parent, label);
    if (--ctx->active_count[root] == 0) {
        insert_component(out, count, out_max,
                         ctx->components[root].area,
                         ctx->components[root].sum_x,
                         ctx->components[root].sum_y);
        component_clear(&ctx->components[root]);
        ctx->parent[root] = root;
        ctx->reuse_after_row[root] = row - 1;
    }
}

static void close_weighted_column(OSTCCContext *ctx, int x, int row,
                                  OSTCCComponent *out, int *count,
                                  int out_max)
{
    int label;
    int root;

    label = ctx->col_label[x];
    if (!label)
        return;

    ctx->col_label[x] = 0;
    root = root_compress(ctx->parent, label);
    if (--ctx->active_count[root] == 0) {
        insert_weighted_component(out, count, out_max, ctx->components[root]);
        component_clear(&ctx->components[root]);
        ctx->parent[root] = root;
        ctx->reuse_after_row[root] = row - 1;
    }
}

int ost_cc_buffer_sizes(int width, OSTCCBufferSizes *sizes)
{
    int max_labels;

    if (!sizes || width <= 0)
        return -1;

    max_labels = width + 1;
    sizes->max_labels = max_labels;
    sizes->components = (size_t)max_labels;
    sizes->parent = (size_t)max_labels;
    sizes->col_label = (size_t)width;
    sizes->active_count = (size_t)max_labels;
    sizes->reuse_after_row = (size_t)max_labels;
    sizes->total_bytes =
        sizes->components * sizeof(OSTCCComponent) +
        (sizes->parent + sizes->col_label + sizes->active_count +
         sizes->reuse_after_row) * sizeof(int);

    return 0;
}

int ost_cc_init(OSTCCContext *ctx, int width,
                OSTCCComponent *components,
                int *parent,
                int *col_label,
                int *active_count,
                int *reuse_after_row)
{
    OSTCCBufferSizes sizes;

    if (!ctx || !components || !parent || !col_label || !active_count ||
        !reuse_after_row)
        return -1;
    if (ost_cc_buffer_sizes(width, &sizes) < 0)
        return -1;

    ctx->width = width;
    ctx->max_labels = sizes.max_labels;
    ctx->components = components;
    ctx->parent = parent;
    ctx->col_label = col_label;
    ctx->active_count = active_count;
    ctx->reuse_after_row = reuse_after_row;
    return 0;
}

int ost_cc_threshold_4(const unsigned char *image,
                       int width, int height, int stride,
                       unsigned char threshold,
                       OSTCCComponent *out, int out_max,
                       OSTCCContext *ctx)
{
    int count;
    int next_reuse_label;

    if (!image || !ctx || !out || width <= 0 || height < 0 ||
        stride < width || out_max < 0 || ctx->width != width)
        return -1;

    count = 0;
    next_reuse_label = 1;
    clear_components(out, out_max);
    memset(ctx->col_label, 0, (size_t)width * sizeof(int));
    memset(ctx->active_count, 0, (size_t)ctx->max_labels * sizeof(int));
    memset(ctx->reuse_after_row, -1,
           (size_t)ctx->max_labels * sizeof(int));

    for (int y = 0; y < height; y++) {
        const unsigned char *row;
        int left;

        row = image + (size_t)y * (size_t)stride;
        left = 0;
        for (int x = 0; x < width; x++) {
            int top_label;
            int top;
            int label;

            top_label = ctx->col_label[x];
            if (row[x] <= threshold) {
                left = 0;
                close_binary_column(ctx, x, y, out, &count, out_max);
                continue;
            }

            top = top_label ? root_compress(ctx->parent, top_label) : 0;
            if (left)
                label = (top && top != left) ?
                    merge_roots(ctx, left, top, y) : left;
            else if (top)
                label = top;
            else {
                label = alloc_label(ctx, y, &next_reuse_label);
                if (!label)
                    return -2;
            }

            if (!top_label)
                ctx->active_count[label]++;
            ctx->col_label[x] = label;
            component_add_binary_pixel(&ctx->components[label], x, y);
            left = label;
        }
    }

    for (int x = 0; x < width; x++)
        close_binary_column(ctx, x, height, out, &count, out_max);
    return count;
}

/* background and fitting */
#include <math.h>
#include <string.h>

/* C port of FastSExtractorDecoupled.py's tiled histogram background pass. */

static int iround_even(double x)
{
    return (int)lrint(x);
}

int ost_bg_config_init(OSTBGConfig *cfg, int width, int height)
{
    if (!cfg || width <= 0 || height <= 0)
        return -1;
    cfg->width = width;
    cfg->height = height;
    cfg->tile_size = 64;
    cfg->map_width = (width + cfg->tile_size - 1) / cfg->tile_size;
    cfg->map_height = (height + cfg->tile_size - 1) / cfg->tile_size;
    cfg->max_stars = 256;
    cfg->max_pixel_brightness = 255 * 4;
    cfg->sample_radius = 2;
    cfg->psf_sigma = 0.0;
    cfg->threshold_sigma = 5.0;
    cfg->detect_sigma = 1.5;
    return 0;
}

void ost_bg_rgba_to_gray16(uint16_t *dst, const unsigned char *rgba,
                           int width, int height)
{
    int n;

    n = width * height;
    for (int i = 0; i < n; i++)
        dst[i] = (uint16_t)(rgba[4 * i] + 2 * rgba[4 * i + 1] + rgba[4 * i + 2]);
}

static void make_edges(int *edge, int n, int bins)
{
    if (bins <= 1) {
        edge[0] = 0;
        edge[1] = n;
        return;
    }
    for (int i = 0; i <= bins; i++)
        edge[i] = iround_even((double)i * (double)n / (double)bins);
}

/* Pick the largest low-intensity histogram interval that stays self-consistent. */
static void bg_from_hist(const int *hist, int hist_len, double sigma,
                         double *mean, double *var)
{
    double prev_c0 = 0.0;
    double prev_c1 = 0.0;
    double prev_c2 = 0.0;
    double prev_t = 0.0;
    int prev_bin = 0;
    double best_c0 = -1.0;
    double best_mean = 0.0;
    double best_var = 0.0;
    int have_prev = 0;

    for (int b = 0; b < hist_len; b++) {
        double c0;
        double c1;
        double c2;
        double m;
        double v;
        double t;

        if (hist[b] == 0)
            continue;
        c0 = hist[b];
        c1 = (double)hist[b] * b;
        c2 = (double)hist[b] * b * b;
        if (have_prev && prev_bin < prev_t) {
            c0 += prev_c0;
            c1 += prev_c1;
            c2 += prev_c2;
        }
        m = c1 / c0;
        /* +1/2 and +1/12 compensate for quantization inside integer bins. */
        v = c2 / c0 - m * m + 1.0 / 12.0;
        if (v < 0)
            v = 0;
        m += 0.5;
        t = m + sqrt(v) * sigma;
        if (c0 > best_c0) {
            best_c0 = c0;
            best_mean = m;
            best_var = v;
        }
        prev_c0 = c0;
        prev_c1 = c1;
        prev_c2 = c2;
        prev_t = t;
        prev_bin = b;
        have_prev = 1;
    }
    *mean = best_mean;
    *var = best_var;
}

int ost_bg_compute_stats(const OSTBGConfig *cfg, const uint16_t *image,
                         int stride, double *mean, double *var,
                         double *poisson, int *x_edge, int *y_edge,
                         int *hist, int hist_len)
{
    int mw;
    int mh;

    if (!cfg || !image || !mean || !var || !poisson || !x_edge || !y_edge ||
        !hist || hist_len <= 0 || stride < cfg->width)
        return -1;

    mw = cfg->map_width;
    mh = cfg->map_height;
    /* Tile edges are explicit, so image dimensions need not divide tile_size. */
    make_edges(x_edge, cfg->width, mw);
    make_edges(y_edge, cfg->height, mh);

    for (int ty = 0; ty < mh; ty++) {
        for (int tx = 0; tx < mw; tx++) {
            double m;
            double v;
            int x0 = x_edge[tx];
            int x1 = x_edge[tx + 1];
            int y0 = y_edge[ty];
            int y1 = y_edge[ty + 1];

            memset(hist, 0, (size_t)hist_len * sizeof(hist[0]));
            for (int y = y0; y < y1; y++) {
                const uint16_t *row = image + (size_t)y * (size_t)stride;
                for (int x = x0; x < x1; x++) {
                    int p = row[x];
                    if (p >= hist_len)
                        p = hist_len - 1;
                    hist[p]++;
                }
            }
            bg_from_hist(hist, hist_len, cfg->threshold_sigma, &m, &v);
            mean[ty * mw + tx] = m;
            var[ty * mw + tx] = v;
            poisson[ty * mw + tx] = m > 0 ? v / m : 0;
        }
    }
    return 0;
}

static double bilinear_grid(const double *map, int mw, int mh, double x, double y)
{
    int x0;
    int y0;
    int x1;
    int y1;
    double fx;
    double fy;
    double a;
    double b;

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > mw - 1) x = mw - 1;
    if (y > mh - 1) y = mh - 1;
    x0 = (int)floor(x);
    y0 = (int)floor(y);
    x1 = x0 + 1 < mw ? x0 + 1 : x0;
    y1 = y0 + 1 < mh ? y0 + 1 : y0;
    fx = x - x0;
    fy = y - y0;
    a = map[y0 * mw + x0] * (1.0 - fx) + map[y0 * mw + x1] * fx;
    b = map[y1 * mw + x0] * (1.0 - fx) + map[y1 * mw + x1] * fx;
    return a * (1.0 - fy) + b * fy;
}

static inline double tile_threshold_at(const OSTBGConfig *cfg,
                                       const double *mean, const double *var,
                                       int tx, int ty)
{
    int i = ty * cfg->map_width + tx;
    return mean[i] + cfg->detect_sigma * sqrt(var[i]);
}

static inline double threshold_grid(const OSTBGConfig *cfg, const double *mean,
                                    const double *var, double x, double y)
{
    int mw = cfg->map_width;
    int mh = cfg->map_height;
    int x0;
    int y0;
    int x1;
    int y1;
    double fx;
    double fy;
    double a;
    double b;

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > mw - 1) x = mw - 1;
    if (y > mh - 1) y = mh - 1;
    x0 = (int)floor(x);
    y0 = (int)floor(y);
    x1 = x0 + 1 < mw ? x0 + 1 : x0;
    y1 = y0 + 1 < mh ? y0 + 1 : y0;
    fx = x - x0;
    fy = y - y0;
    a = tile_threshold_at(cfg, mean, var, x0, y0) * (1.0 - fx) +
        tile_threshold_at(cfg, mean, var, x1, y0) * fx;
    b = tile_threshold_at(cfg, mean, var, x0, y1) * (1.0 - fx) +
        tile_threshold_at(cfg, mean, var, x1, y1) * fx;
    return a * (1.0 - fy) + b * fy;
}

static double center_coord(int p, const int *edge, int n)
{
    double c0;
    double c1;

    c0 = 0.5 * (edge[1] + edge[0]);
    if (p <= c0)
        return 0;
    c1 = 0.5 * (edge[n] + edge[n - 1]);
    if (p >= c1)
        return n - 1;
    for (int i = 0; i + 1 < n; i++) {
        double a = 0.5 * (edge[i + 1] + edge[i]);
        double b = 0.5 * (edge[i + 2] + edge[i + 1]);
        if (p <= b)
            return i + ((double)p - a) / (b - a);
    }
    return n - 1;
}

static double ost_bg_interpolate(void *opaque, int x, int y, double *var)
{
    OSTBGStats *s = (OSTBGStats *)opaque;
    const OSTBGConfig *cfg = s->cfg;
    double gx = center_coord(x, s->x_edge, cfg->map_width);
    double gy = center_coord(y, s->y_edge, cfg->map_height);

    if (var)
        *var = bilinear_grid(s->var, cfg->map_width, cfg->map_height, gx, gy);
    return bilinear_grid(s->mean, cfg->map_width, cfg->map_height, gx, gy);
}

typedef struct OSTBGThresholdTest {
    const OSTBGConfig *cfg;
    const double *mean;
    const double *var;
    double sx;
    double sy;
} OSTBGThresholdTest;

static double bg_threshold_y(const OSTBGThresholdTest *t, int tx,
                             int y0, int y1, double fy)
{
    const OSTBGConfig *cfg = t->cfg;
    double a = tile_threshold_at(cfg, t->mean, t->var, tx, y0);
    double b = tile_threshold_at(cfg, t->mean, t->var, tx, y1);

    return a * (1.0 - fy) + b * fy;
}

int ost_bg_extract_fused(const OSTBGConfig *cfg, const uint16_t *image,
                         int stride, OSTBGStats *stats, OSTCCContext *cc,
                         OSTCCComponent *stars, int stars_max)
{
    /* Thresholding, connected components, and first-pass photometry are fused. */
    OSTBGThresholdTest test;
    int count;
    int next_reuse_label;

    if (!cfg || !image || !stats || !cc || !stars ||
        stride < cfg->width || stars_max < 0 || cc->width != cfg->width)
        return -1;

    test.cfg = cfg;
    test.mean = stats->mean;
    test.var = stats->var;
    test.sx = cfg->width > 1 ?
        (double)(cfg->map_width - 1) / (cfg->width - 1) : 0.0;
    test.sy = cfg->height > 1 ?
        (double)(cfg->map_height - 1) / (cfg->height - 1) : 0.0;
    count = 0;
    next_reuse_label = 1;

    clear_components(stars, stars_max);
    memset(cc->col_label, 0, (size_t)cfg->width * sizeof(int));
    memset(cc->active_count, 0, (size_t)cc->max_labels * sizeof(int));
    memset(cc->reuse_after_row, -1, (size_t)cc->max_labels * sizeof(int));

    for (int y = 0; y < cfg->height; y++) {
        const uint16_t *row;
        double gy;
        int mh;
        int mw;
        int y0;
        int y1;
        double fy;
        int x;
        int left;

        row = image + (size_t)y * (size_t)stride;
        gy = y * test.sy;
        mh = cfg->map_height;
        mw = cfg->map_width;
        y0 = (int)floor(gy);
        if (y0 < 0)
            y0 = 0;
        if (y0 > mh - 1)
            y0 = mh - 1;
        y1 = y0 + 1 < mh ? y0 + 1 : y0;
        fy = gy - y0;
        x = 0;
        left = 0;

        while (x < cfg->width) {
            int tx;
            int tx1;
            int x_end;
            double gx;
            double a;
            double b;
            double th;
            double dth;

            if (test.sx <= 0.0) {
                tx = 0;
                tx1 = 0;
                x_end = cfg->width;
                gx = 0.0;
            } else {
                gx = x * test.sx;
                tx = (int)floor(gx);
                if (tx < 0)
                    tx = 0;
                if (tx > mw - 1)
                    tx = mw - 1;
                tx1 = tx + 1 < mw ? tx + 1 : tx;
                x_end = tx + 1 < mw ? (int)ceil((tx + 1) / test.sx) : cfg->width;
                if (x_end > cfg->width)
                    x_end = cfg->width;
            }
            a = bg_threshold_y(&test, tx, y0, y1, fy);
            b = bg_threshold_y(&test, tx1, y0, y1, fy);
            th = a + (b - a) * (gx - tx);
            dth = (b - a) * test.sx;

            while (x < x_end) {
                int top_label;

                top_label = cc->col_label[x];
                if (row[x] <= th) {
                    left = 0;
                    if (top_label)
                        close_weighted_column(cc, x, y, stars, &count, stars_max);
                } else {
                    int top;
                    int label;

                    top = top_label ? root_compress(cc->parent, top_label) : 0;
                    if (left)
                        label = (top && top != left) ?
                            merge_roots(cc, left, top, y) : left;
                    else if (top)
                        label = top;
                    else {
                        label = alloc_label(cc, y, &next_reuse_label);
                        if (!label)
                            return -2;
                    }

                    if (!top_label)
                        cc->active_count[label]++;
                    cc->col_label[x] = label;
                    component_add_weighted_pixel(&cc->components[label], image,
                                                 stride, x, y,
                                                 ost_bg_interpolate, stats,
                                                 cfg->threshold_sigma);
                    left = label;
                }
                x++;
                th += dth;
            }
        }
    }

    for (int x = 0; x < cfg->width; x++)
        close_weighted_column(cc, x, cfg->height, stars, &count, stars_max);
    return count;
}

static void bg_at(const OSTBGStats *s, int x, int y,
                  double *mu, double *var, double *poisson)
{
    const OSTBGConfig *cfg = s->cfg;
    double gx = center_coord(x, s->x_edge, cfg->map_width);
    double gy = center_coord(y, s->y_edge, cfg->map_height);

    *mu = bilinear_grid(s->mean, cfg->map_width, cfg->map_height, gx, gy);
    *var = bilinear_grid(s->var, cfg->map_width, cfg->map_height, gx, gy);
    *poisson = bilinear_grid(s->poisson, cfg->map_width, cfg->map_height, gx, gy);
}

/* Weighted moments seed x, y, and flux; calibration supplies PSF width. */
static int ost_bg_fit_init(const OSTBGConfig *cfg,
                           const OSTCCComponent *components,
                           int component_count, OSTBGFitStar *stars,
                           double *params, int max_stars)
{
    int n = 0;

    if (!cfg || !components || !stars || !params || max_stars <= 0)
        return -1;

    for (int i = 0; i < component_count && n < max_stars; i++) {
        double x;
        double y;
        int xi;
        int yi;

        if (components[i].wsum <= 0)
            continue;
        x = components[i].wx / components[i].wsum;
        y = components[i].wy / components[i].wsum;
        xi = (int)lrint(x);
        yi = (int)lrint(y);
        if (xi < cfg->sample_radius ||
            xi >= cfg->width - cfg->sample_radius ||
            yi < cfg->sample_radius ||
            yi >= cfg->height - cfg->sample_radius)
            continue;
        stars[n].xi = xi;
        stars[n].yi = yi;
        params[3 * n + 0] = x;
        params[3 * n + 1] = y;
        params[3 * n + 2] = components[i].wsum;
        n++;
    }
    params[3 * n] = cfg->psf_sigma;
    return n;
}

/* Pixel-integrated symmetric Gaussian model and Jacobian. */
static void psf_eval(double x0, double y0, double I, double sigma,
                     int px, int py, double *pred, double *J)
{
    const double sqrt2 = 1.4142135623730950488;
    const double sqrt_2pi = 2.5066282746310005024;
    double s = sigma * sqrt2;
    double dx = (double)px - x0;
    double dy = (double)py - y0;
    double x1 = (dx - 0.5) / s;
    double x2 = (dx + 0.5) / s;
    double y1 = (dy - 0.5) / s;
    double y2 = (dy + 0.5) / s;
    double ex = erf(x1) - erf(x2);
    double ey = erf(y1) - erf(y2);
    double e1 = exp(-x1 * x1);
    double e2 = exp(-x2 * x2);
    double e3 = exp(-y1 * y1);
    double e4 = exp(-y2 * y2);
    double norm = I / (sigma * 2.0 * sqrt_2pi);

    *pred = ex * ey * I * 0.25;
    J[0] = -norm * ey * (e1 - e2);
    J[1] = -norm * ex * (e3 - e4);
    J[2] = I != 0.0 ? *pred / I : 0.0;
}

/* Weighted least squares uses background plus Poisson-scaled image variance. */
static int build_fit_model(const OSTBGConfig *cfg, const uint16_t *image,
                           int stride, const OSTBGStats *stats,
                           const OSTBGFitStar *stars, const double *params,
                           int n, double sigma,
                           OSTBGFitStar *stars_out, double *params_out,
                           double *normal, double *rhs,
                           double *dropped, int *dropped_count,
                           int max_dropped)
{
    int m = 0;
    int r = cfg->sample_radius;

    for (int i = 0; i < n; i++) {
        OSTSym3 B = {0, 0, 0, 0, 0, 0};
        OSTDVec3 b = {0, 0, 0};
        double x0 = params[3 * i + 0];
        double y0 = params[3 * i + 1];
        double I = params[3 * i + 2];
        int n_valid = 0;
        int bad_var = 0;

        for (int x = stars[i].xi - r; x <= stars[i].xi + r; x++) {
            for (int y = stars[i].yi - r; y <= stars[i].yi + r; y++) {
                double mu;
                double var;
                double poisson;
                double obs;
                double pred;
                double J[3];

                bg_at(stats, x, y, &mu, &var, &poisson);
                obs = (double)image[(size_t)y * (size_t)stride + (size_t)x] - mu;
                var = fmax((double)image[(size_t)y * (size_t)stride + (size_t)x] *
                           poisson, 0.0) + var;
                psf_eval(x0, y0, I, sigma, x, y, &pred, J);
                if (pred >= cfg->max_pixel_brightness)
                    continue;
                n_valid++;
                if (var == 0.0) {
                    bad_var = 1;
                    continue;
                } else {
                    double w = 1.0 / var;
                    double e = obs - pred;
                    B[0] += J[0] * w * J[0];
                    B[1] += J[1] * w * J[0];
                    B[2] += J[1] * w * J[1];
                    B[3] += J[2] * w * J[0];
                    B[4] += J[2] * w * J[1];
                    B[5] += J[2] * w * J[2];
                    b[0] += J[0] * w * e;
                    b[1] += J[1] * w * e;
                    b[2] += J[2] * w * e;
                }
            }
        }
        if (n_valid < OST_MIN_COMPONENT_AREA || bad_var) {
            if (*dropped_count < max_dropped) {
                dropped[3 * *dropped_count + 0] = params[3 * i + 0];
                dropped[3 * *dropped_count + 1] = params[3 * i + 1];
                dropped[3 * *dropped_count + 2] = params[3 * i + 2];
            }
            (*dropped_count)++;
            continue;
        }
        stars_out[m] = stars[i];
        params_out[3 * m + 0] = params[3 * i + 0];
        params_out[3 * m + 1] = params[3 * i + 1];
        params_out[3 * m + 2] = params[3 * i + 2];
        if (normal) {
            memcpy(normal + 6 * m, B, sizeof(B));
            memcpy(rhs + 3 * m, b, sizeof(b));
        }
        m++;
    }
    params_out[3 * m] = sigma;
    return m;
}

static int solve_fit(double *params, int n, double *normal, double *rhs,
                     double *rhs_solve, double *cov_xy, int compute_cov)
{
    for (int i = 0; i < n; i++) {
        OSTSym3 l;
        OSTDVec3 y;

        if (ost_chol(normal + 6 * i, 3, l) < 0)
            return -1;
        ost_chol_solve(l, 3, rhs + 3 * i, rhs_solve + 3 * i, y);
        params[3 * i + 0] += 0.5 * rhs_solve[3 * i + 0];
        params[3 * i + 1] += 0.5 * rhs_solve[3 * i + 1];
        params[3 * i + 2] += 0.5 * rhs_solve[3 * i + 2];
        if (compute_cov) {
            OSTDVec3 d, e, x;
            ost_chol_inv_diag(l, 3, d, e, x, y);
            cov_xy[2 * i + 0] = d[0];
            cov_xy[2 * i + 1] = d[1];
        }
    }
    return 0;
}

int ost_bg_fit_stars(const OSTBGConfig *cfg, const uint16_t *image, int stride,
                     const OSTBGStats *stats, const OSTCCComponent *components,
                     int component_count, int num_iter,
                     OSTBGFitWorkspace *work, int max_stars,
                     double *params_out, double *cov_xy_out,
                     int *dropped_count_out)
{
    OSTBGFitStar *stars;
    OSTBGFitStar *stars_next;
    double *params;
    double *params_next;
    int n;
    int dropped_count = 0;

    if (!cfg || !image || !stats || !components || !work || !params_out ||
        !cov_xy_out || !dropped_count_out || stride < cfg->width ||
        max_stars <= 0 || cfg->psf_sigma <= 0.0)
        return -1;

    stars = work->stars1;
    stars_next = work->stars2;
    params = work->params1;
    params_next = work->params2;
    if (!stars || !stars_next || !params || !params_next || !work->normal ||
        !work->rhs || !work->rhs_solve || !work->cov_xy || !work->dropped)
        return -1;

    n = ost_bg_fit_init(cfg, components, component_count, stars, params, max_stars);
    if (n <= 0) {
        *dropped_count_out = 0;
        return n;
    }

    for (int it = 0; it < num_iter; it++) {
        double sigma = cfg->psf_sigma;
        int compute_cov = it == num_iter - 1;
        int m;

        m = build_fit_model(cfg, image, stride, stats, stars, params, n, sigma,
                            stars_next, params_next, work->normal, work->rhs,
                            work->dropped, &dropped_count, max_stars);
        if (m <= 0)
            break;
        if (solve_fit(params_next, m, work->normal, work->rhs,
                      work->rhs_solve, work->cov_xy, compute_cov) < 0)
            return -2;
        {
            OSTBGFitStar *ts = stars;
            double *tp = params;
            stars = stars_next;
            stars_next = ts;
            params = params_next;
            params_next = tp;
            n = m;
        }
    }

    {
        double sigma = cfg->psf_sigma;
        int m;

        m = build_fit_model(cfg, image, stride, stats, stars, params, n, sigma,
                            stars_next, params_next, NULL, NULL,
                            work->dropped, &dropped_count, max_stars);
        if (m > 0) {
            memcpy(params_out, params_next, (size_t)(3 * m + 1) * sizeof(*params_out));
            for (int i = 0; i < m; i++) {
                cov_xy_out[2 * i + 0] = work->cov_xy[2 * i + 0];
                cov_xy_out[2 * i + 1] = work->cov_xy[2 * i + 1];
            }
            n = m;
        }
    }
    *dropped_count_out = dropped_count;
    return n;
}

OST_DEF int ost_load_config(Config *c, const char *filename)
{
    FILE *f;
    char line[256], k[128], v[128];

    memset(c, 0, sizeof(*c));
    f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        return -1;
    }
    while (fgets(line, sizeof(line), f)) {
        if (sscanf(line, " %127[^=]=%127s", k, v) != 2)
            continue;
        if (!strcmp(k, "IMG_X")) c->IMG_X = atoi(v);
        else if (!strcmp(k, "IMG_Y")) c->IMG_Y = atoi(v);
        else if (!strcmp(k, "PIXSCALE")) c->PIXSCALE = (float)atof(v);
        else if (!strcmp(k, "POS_ERR_SIGMA")) c->POS_ERR_SIGMA = (float)atof(v);
        else if (!strcmp(k, "PSF_SIGMA")) c->PSF_SIGMA = (float)atof(v);
        else if (!strcmp(k, "POS_VARIANCE")) c->POS_VARIANCE = (float)atof(v);
        else if (!strcmp(k, "IMAGE_VARIANCE")) c->IMAGE_VARIANCE = (float)atof(v);
        else if (!strcmp(k, "THRESH_FACTOR")) c->THRESH_FACTOR = (float)atof(v);
        else if (!strcmp(k, "DOUBLE_STAR_PX")) c->DOUBLE_STAR_PX = (float)atof(v);
        else if (!strcmp(k, "MAX_FALSE_STARS")) c->MAX_FALSE_STARS = atoi(v);
        else if (!strcmp(k, "DB_REDUNDANCY")) c->DB_REDUNDANCY = atoi(v);
        else if (!strcmp(k, "REQUIRED_STARS")) c->REQUIRED_STARS = atoi(v);
        else if (!strcmp(k, "BASE_FLUX")) c->BASE_FLUX = (float)atof(v);
    }
    fclose(f);
    if (c->PSF_SIGMA <= 0.0f) {
        fprintf(stderr, "%s: missing or invalid PSF_SIGMA\n", filename);
        return -1;
    }
    c->MAXFOV = c->PIXSCALE * sqrt(c->IMG_X * c->IMG_X + c->IMG_Y * c->IMG_Y);
    c->MINFOV = c->PIXSCALE * c->IMG_Y;
    c->MATCH_VALUE = 4 * log(1.0 / (c->IMG_X * c->IMG_Y)) + log(2 * PI);
    c->PIXX_TANGENT = 2 * tan((c->IMG_X * c->PIXSCALE / 3600) * PI / (180 * 2)) / c->IMG_X;
    c->PIXY_TANGENT = 2 * tan((c->IMG_Y * c->PIXSCALE / 3600) * PI / (180 * 2)) / c->IMG_Y;
    c->KDBUCKET_SIZE = (c->IMG_X * c->PIXSCALE / 3600) *
                       (c->IMG_Y * c->PIXSCALE / 3600) * 3.5;
    return 0;
}

static uint64_t v3_key(const Vec3 v)
{
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < 3; i++) {
        uint32_t u;
        memcpy(&u, &v[i], sizeof(u));
        h = (h ^ u) * 1099511628211ULL;
    }
    return h ? h : 1;
}

OST_DEF Star ost_make_db_star(const Config *c, float x, float y, float z, float flux, int id)
{
    Star s;
    s.v[0] = x; s.v[1] = y; s.v[2] = z; s.flux = flux; s.id = id;
    s.px = y / (x * c->PIXX_TANGENT);
    s.py = z / (x * c->PIXY_TANGENT);
    s.star_idx = -1;
    s.unreliable = 0;
    s.sigma_sq = c->POS_VARIANCE;
    return s;
}

OST_DEF Star ost_make_img_star(const Config *c, float px, float py, float flux, int id)
{
    /* Image centroids live on the tangent plane and become camera-frame vectors. */
    Star s;
    float j, k;
    s.px = px; s.py = py; s.flux = flux; s.id = id;
    j = c->PIXX_TANGENT * px;
    k = c->PIXY_TANGENT * py;
    s.v[0] = 1. / sqrtf(j * j + k * k + 1);
    s.v[1] = j * s.v[0];
    s.v[2] = k * s.v[0];
    s.star_idx = -1;
    s.unreliable = 0;
    s.sigma_sq = c->IMAGE_VARIANCE / flux;
    return s;
}

OST_DEF void ost_star_db_init(StarDB *db, Star *storage, int cap)
{
    db->v = storage;
    db->n = 0;
    db->cap = cap;
    db->max_variance = 0.0f;
}

OST_DEF int ost_db_add(StarDB *db, Star s)
{
    if (db->n >= db->cap)
        return -1;
    if (db->max_variance < s.sigma_sq)
        db->max_variance = s.sigma_sq;
    s.star_idx = db->n;
    db->v[db->n++] = s;
    return 0;
}

static int db_copy(StarDB *dst, const StarDB *src)
{
    if (src->n > dst->cap)
        return -1;
    memcpy(dst->v, src->v, (size_t)src->n * sizeof(dst->v[0]));
    dst->n = src->n;
    dst->max_variance = src->max_variance;
    for (int i = 0; i < dst->n; i++)
        dst->v[i].star_idx = i;
    return 0;
}

static int key_seen(uint64_t *tab, size_t tab_len, uint64_t h)
{
    if (!tab || tab_len == 0)
        return -1;
    size_t i = (size_t)(h ^ (h >> 32)) % tab_len;
    for (size_t n = 0; n < tab_len; n++) {
        if (!tab[i]) {
            tab[i] = h;
            return 0;
        }
        if (tab[i] == h)
            return 1;
        i = (i + 1) % tab_len;
    }
    return -1;
}

OST_DEF int ost_load_catalog(const Config *cfg, StarDB *db,
                                  const char *filename, float year,
                                  uint64_t *cat_keys, size_t cat_key_cap)
{
    /* Caller sizes both the catalog DB and the de-duplication hash table. */
    FILE *f;
    char line[2048], *field[78];
    float yd = year - 1991.25f;

    if (!cat_keys || cat_key_cap == 0)
        return -1;
    f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        return -1;
    }
    memset(cat_keys, 0, cat_key_cap * sizeof(cat_keys[0]));
    db->max_variance = cfg->POS_VARIANCE;
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        int nf = 0;
        while (nf < 78) {
            field[nf++] = p;
            p = strchr(p, '|');
            if (!p)
                break;
            *p++ = 0;
        }
        if (nf > 29) {
            float mag = (float)atof(field[5]);
            float dec = yd * atof(field[13]) / 3600000.0 + atof(field[9]);
            float cosdec = cos(PI * dec / 180.0);
            float ra = yd * atof(field[12]) / (cosdec * 3600000.0) + atof(field[8]);
            Star s = ost_make_db_star(cfg,
                                  cos(PI * ra / 180.0) * cosdec,
                                  sin(PI * ra / 180.0) * cosdec,
                                  sin(PI * dec / 180.0),
                                  cfg->BASE_FLUX * powf(10.0f, -mag / 2.5f),
                                  atoi(field[1]));
            int seen = key_seen(cat_keys, cat_key_cap, v3_key(s.v));
            if (seen < 0 || (!seen && ost_db_add(db, s) < 0)) {
                fclose(f);
                return -1;
            }
        }
    }
    fclose(f);
    return 0;
}

static int cmp_flux_desc(const void *a, const void *b)
{
    float d = ((const Star *)b)->flux - ((const Star *)a)->flux;
    return (d < 0) ? -1 : (d > 0);
}

OST_DEF void ost_query_init(Query *q, StarDB *db, Star *map, int *res, signed char *mask)
{
    q->map = map;
    q->n = db->n;
    q->kdresults = res;
    q->kdmask = mask;
    q->kdresults_size = db->n;
    q->kdresults_maxsize = 0x7fffffff;
    q->kdsorted = 0;
    memcpy(q->map, db->v, (size_t)db->n * sizeof(q->map[0]));
    for (int i = 0; i < db->n; i++)
        q->kdresults[i] = i;
    memset(q->kdmask, 0, (size_t)db->n);
}

OST_DEF void ost_query_reset_mask(Query *q)
{
    memset(q->kdmask, 0, (size_t)(unsigned)q->n);
}

OST_DEF void ost_query_clear_results(Query *q)
{
    while (q->kdresults_size > 0)
        q->kdmask[q->kdresults[--q->kdresults_size]] = 0;
}

OST_DEF void ost_query_sort_flux(Query *q)
{
    qsort(q->map, (size_t)q->n, sizeof(q->map[0]), cmp_flux_desc);
    q->kdsorted = 0;
}

OST_DEF void ost_query_kdsort(Query *q, const Config *c)
{
    if (!q->kdsorted) {
        ost_kdbuild(q->map, sizeof(q->map[0]), 0, 0, q->n,
                    c->KDBUCKET_SIZE, 0, 3);
        q->kdsorted = 1;
    }
}

typedef struct OSTQuerySearchCtx {
    Query *q;
    const float *p;
    float r;
    float r2;
    float min_flux;
} OSTQuerySearchCtx;

static inline int kdcheck(void *base, int idx, void *vctx)
{
    OSTQuerySearchCtx *ctx = (OSTQuerySearchCtx *)vctx;
    Query *q = ctx->q;
    Star *s = &((Star *)base)[idx];
    float dx = ctx->p[0] - s->v[0], dy = ctx->p[1] - s->v[1], dz = ctx->p[2] - s->v[2];
    if (dx - ctx->r <= 0 && 0 <= dx + ctx->r &&
        dy - ctx->r <= 0 && 0 <= dy + ctx->r &&
        dz - ctx->r <= 0 && 0 <= dz + ctx->r &&
        ctx->min_flux <= s->flux && q->kdmask[idx] == 0 &&
        dx * dx + dy * dy + dz * dz <= ctx->r2) {
        int n = q->kdresults_size++;
        q->kdmask[idx] = 1;
        while (n > 0 && s->flux > q->map[q->kdresults[n - 1]].flux) {
            q->kdresults[n] = q->kdresults[n - 1];
            n--;
        }
        q->kdresults[n] = idx;
        if (q->kdresults_size > q->kdresults_maxsize) {
            q->kdresults_size = q->kdresults_maxsize;
            q->kdmask[q->kdresults[q->kdresults_size]] = 0;
        }
        if (q->kdresults_size == q->kdresults_maxsize)
            ctx->min_flux = q->map[q->kdresults[q->kdresults_size - 1]].flux;
    }
    return 0;
}

OST_DEF void ost_query_search(Query *q, const Config *c, const float p[3],
                         float arcsec, float min_flux)
{
    float a = arcsec / 3600.0f * (float)PI / 180.0f;
    float r = 2.0f * fabsf(sinf(a / 2.0f));
    float rv[3] = {r, r, r};
    OSTQuerySearchCtx ctx = {q, p, r, r * r, min_flux};
    ost_query_kdsort(q, c);
    ost_kdsearch(q->map, sizeof(q->map[0]), 0, 0, q->n,
                 c->KDBUCKET_SIZE, 0, 3, p, rv, kdcheck, &ctx);
}

OST_DEF void ost_query_search_range(Query *q, const Config *c, const float p[3],
                                  float arcsec, float min_flux, int min,
                                  int max, int dim)
{
    float a = arcsec / 3600.0f * (float)PI / 180.0f;
    float r = 2.0f * fabsf(sinf(a / 2.0f));
    float rv[3] = {r, r, r};
    OSTQuerySearchCtx ctx = {q, p, r, r * r, min_flux};
    ost_query_kdsort(q, c);
    ost_kdsearch(q->map, sizeof(q->map[0]), 0, min, max,
                 c->KDBUCKET_SIZE, dim, 3, p, rv, kdcheck, &ctx);
}

OST_DEF void ost_query_mask_filter(Query *q, const Config *c)
{
    /* Mask dim stars and likely double-stars before building the catalog DB. */
    ost_query_kdsort(q, c);
    for (int i = 0; i < q->n; i++) {
        int lastmask = q->kdmask[i];
        ost_query_search(q, c, q->map[i].v,
                     c->DOUBLE_STAR_PX * c->PIXSCALE,
                     c->THRESH_FACTOR * c->IMAGE_VARIANCE);
        if (q->kdresults_size > 1 || lastmask ||
            q->map[i].flux < c->THRESH_FACTOR * c->IMAGE_VARIANCE) {
            q->kdmask[i] = 1;
            q->kdresults_size = 0;
        } else {
            ost_query_clear_results(q);
        }
    }
}

OST_DEF void ost_query_mask_uniform(Query *q, const Config *c, int min_stars, signed char *keep)
{
    /* Keep enough stars per FOV without letting dense sky regions dominate. */
    int oldmax = q->kdresults_maxsize;
    memset(keep, 0, (size_t)q->n);
    q->kdresults_maxsize = min_stars;
    for (int i = 0; i < q->n; i++) if (!q->kdmask[i]) {
        ost_query_search(q, c, q->map[i].v,
                     c->MINFOV / 2, c->THRESH_FACTOR * c->IMAGE_VARIANCE);
        for (int j = 0; j < q->kdresults_size; j++)
            keep[q->kdresults[j]] = 1;
        ost_query_clear_results(q);
    }
    for (int i = 0; i < q->n; i++)
        q->kdmask[i] = keep[i] ? 0 : 1;
    q->kdresults_maxsize = oldmax;
}

OST_DEF int ost_db_from_mask(StarDB *out, Query *q)
{
    out->n = 0;
    for (int i = 0; i < q->n; i++)
        if (!q->kdmask[i] && ost_db_add(out, q->map[i]) < 0)
            return -1;
    return 0;
}

OST_DEF int ost_db_from_results(StarDB *out, Query *q)
{
    out->n = 0;
    for (int i = 0; i < q->kdresults_size; i++)
        if (ost_db_add(out, q->map[q->kdresults[i]]) < 0)
            return -1;
    return 0;
}

static int ost_constellation_before(const Star *stars, uint32_t a, uint32_t b)
{
    if (stars[a].flux > stars[b].flux)
        return 1;
    if (stars[a].flux < stars[b].flux)
        return 0;
    return a < b;
}

static void ost_constellation_canonicalize_ids(const Star *stars,
                                               uint32_t *ids, int k)
{
    for (int i = 1; i < k; i++) {
        uint32_t v = ids[i];
        int j = i;
        while (j > 0 && ost_constellation_before(stars, v, ids[j - 1])) {
            ids[j] = ids[j - 1];
            j--;
        }
        ids[j] = v;
    }
}

static int cmp_constellation(const void *pa, const void *pb)
{
    const Constellation *a = (const Constellation *)pa, *b = (const Constellation *)pb;
    if (a->p < b->p) return -1;
    if (a->p > b->p) return 1;
    if (a->s1 != b->s1) return a->s1 - b->s1;
    return a->s2 - b->s2;
}

OST_DEF int ost_db_from_image(CDB *cdb, StarDB *src,
                                  Star *star_storage, int star_cap,
                                  Query *q, Star *query_map,
                                  int *query_results, signed char *query_mask,
                                  Constellation *cmap, int cmap_cap,
                                  int stars_per_fov)
{
    /* Bright image stars generate all N*(N-1)/2 observed pair separations. */
    int ns, idx = 0;
    ost_star_db_init(&cdb->stars, star_storage, star_cap);
    if (db_copy(&cdb->stars, src) < 0)
        return -1;
    ost_query_init(q, &cdb->stars, query_map, query_results, query_mask);
    ost_query_sort_flux(q);
    cdb->results = *q;
    cdb->map = cmap;
    ns = cdb->stars.n;
    if (ns > stars_per_fov)
        ns = stars_per_fov;
    for (int j = 1; j < ns; j++) {
        for (int i = 0; i < j; i++) {
            if (idx >= cmap_cap)
                return -1;
            cdb->map[idx].p = ost_star_dist_arcsec(&q->map[i], &q->map[j]);
            cdb->map[idx].s1 = q->map[i].star_idx;
            cdb->map[idx].s2 = q->map[j].star_idx;
            idx++;
        }
    }
    qsort(cdb->map, (size_t)idx, sizeof(cdb->map[0]), cmp_constellation);
    for (int i = 0; i < idx; i++)
        cdb->map[i].idx = i;
    cdb->map_size = idx;
    return 0;
}

OST_DEF int ost_db_from_catalog(CDB *cdb, StarDB *src,
                               Star *star_storage, int star_cap,
                               Query *q, Star *query_map,
                               int *query_results, signed char *query_mask,
                               Constellation *cmap, int cmap_cap,
                               int stars_per_fov, const Config *cfg,
                               signed char *keep)
{
    /* Store only catalog pairs that could appear together inside the FOV. */
    int n = 0, out = 0;
    ost_star_db_init(&cdb->stars, star_storage, star_cap);
    if (db_copy(&cdb->stars, src) < 0)
        return -1;
    ost_query_init(q, &cdb->stars, query_map, query_results, query_mask);
    ost_query_mask_uniform(q, cfg, stars_per_fov, keep);
    for (int i = 0; i < q->n; i++) if (!q->kdmask[i]) {
        ost_query_search(q, cfg, q->map[i].v,
                     cfg->MAXFOV, cfg->THRESH_FACTOR * cfg->IMAGE_VARIANCE);
        for (int j = 0; j < q->kdresults_size; j++) {
            int k = q->kdresults[j];
            if (i != k && ost_constellation_before(cdb->stars.v,
                    (uint32_t)q->map[i].star_idx,
                    (uint32_t)q->map[k].star_idx)) {
                if (n >= cmap_cap)
                    return -1;
                cmap[n].p = ost_star_dist_arcsec(&q->map[i], &q->map[k]);
                cmap[n].s1 = q->map[i].star_idx;
                cmap[n].s2 = q->map[k].star_idx;
                n++;
            }
        }
        ost_query_clear_results(q);
    }
    ost_query_reset_mask(q);
    qsort(cmap, (size_t)n, sizeof(cmap[0]), cmp_constellation);
    for (int i = 0; i < n; i++) {
        if (out == 0 || cmp_constellation(&cmap[i], &cmap[out - 1])) {
            cmap[out] = cmap[i];
            cmap[out].idx = out;
            out++;
        }
    }
    cdb->results = *q;
    cdb->map = cmap;
    cdb->map_size = out;
    return 0;
}

OST_DEF int fov_init(StarFov *fov, StarDB *stars, float db_max_variance,
                     const Config *c, MatchWork *w)
{
    /* Rasterize image-star acceptance regions for fast projected-star scoring. */
    fov->mask = w->fov_mask;
    fov->collision = w->collision;
    fov->collision_size = 0;
    fov->collision_cap = w->collision_cap;
    fov->s_px = w->fov_px;
    fov->s_py = w->fov_py;
    float sigma_sq = stars->max_variance + db_max_variance;
    fov->sigma_sq = sigma_sq;
    fov->maxdist_sq = -sigma_sq * (log(sigma_sq) + c->MATCH_VALUE);
    for (int i = 0; i < c->IMG_X * c->IMG_Y; i++)
        fov->mask[i] = -1;
    float maxdist = sqrtf(fov->maxdist_sq);
    for (int id = 0; id < stars->n; id++) {
        int xmin, xmax, ymin, ymax;
        fov->s_px[id] = stars->v[id].px;
        fov->s_py[id] = stars->v[id].py;
        xmin = (int)(fov->s_px[id] - maxdist - 1);
        xmax = (int)(fov->s_px[id] + maxdist + 1);
        ymin = (int)(fov->s_py[id] - maxdist - 1);
        ymax = (int)(fov->s_py[id] + maxdist + 1);
        if (xmax > c->IMG_X / 2) xmax = c->IMG_X / 2;
        if (xmin < -c->IMG_X / 2) xmin = -c->IMG_X / 2;
        if (ymax > c->IMG_Y / 2) ymax = c->IMG_Y / 2;
        if (ymin < -c->IMG_Y / 2) ymin = -c->IMG_Y / 2;
        for (int x = xmin; x < xmax; x++) for (int y = ymin; y < ymax; y++) {
            float dx = x - fov->s_px[id], dy = y - fov->s_py[id];
            float score;
            int px, py, old;
            if (dx < -0.5f) dx += 1.0f;
            if (dy < -0.5f) dy += 1.0f;
            score = (fov->maxdist_sq - (dx * dx + dy * dy)) /
                    (2.0f * fov->sigma_sq);
            if (score <= 0)
                continue;
            px = x + c->IMG_X / 2;
            py = y + c->IMG_Y / 2;
            old = fov->mask[px + py * c->IMG_X];
            if (old != -1) {
                if (fov->collision_size + 2 > fov->collision_cap)
                    return -1;
                fov->collision_size += 2;
                fov->mask[px + py * c->IMG_X] = -fov->collision_size;
                fov->collision[fov->collision_size - 2] = id;
                fov->collision[fov->collision_size - 1] = old;
            } else {
                fov->mask[px + py * c->IMG_X] = id;
            }
        }
    }
    return 0;
}

OST_DEF float fov_score(StarFov *fov, int id, float px, float py)
{
    float dx = px - fov->s_px[id], dy = py - fov->s_py[id];
    if (dx < -0.5f) dx += 1.0f;
    if (dy < -0.5f) dy += 1.0f;
    return (fov->maxdist_sq - (dx * dx + dy * dy)) /
           (2.0f * fov->sigma_sq);
}

OST_DEF int fov_resolve(StarFov *fov, int id, float px, float py)
{
    /* Colliding acceptance regions are resolved by whichever star scores better. */
    int id1, id2;
    if (id >= -1)
        return id;
    id = -id;
    id1 = fov_resolve(fov, fov->collision[id - 2], px, py);
    id2 = fov_resolve(fov, fov->collision[id - 1], px, py);
    return fov_score(fov, id1, px, py) > fov_score(fov, id2, px, py) ? id1 : id2;
}

OST_DEF int fov_get_id(StarFov *fov, const Config *c, float px, float py)
{
    int nx = (int)(px + c->IMG_X / 2.0f);
    int ny = (int)(py + c->IMG_Y / 2.0f);
    int id = -1;
    if (nx == -1) nx++;
    else if (nx == c->IMG_X) nx--;
    if (ny == -1) ny++;
    else if (ny == c->IMG_Y) ny--;
    if (nx >= 0 && nx < c->IMG_X && ny >= 0 && ny < c->IMG_Y)
        id = fov->mask[nx + ny * c->IMG_X];
    return fov_resolve(fov, id, px, py);
}

OST_DEF void mr_init(MatchResult *m, CDB *db, CDB *img, StarFov *mask, int *map)
{
    memset(m, 0, sizeof(*m));
    m->db = db; m->img = img; m->img_mask = mask; m->map = map;
    m->map_size = img->stars.n;
    m->match.totalscore = -FLT_MAX;
}

OST_DEF void mr_copy(MatchResult *dst, MatchResult *src)
{
    CDB *db = dst->db, *img = dst->img;
    StarFov *mask = dst->img_mask;
    int *map = dst->map, map_size = dst->map_size;
    *dst = *src;
    dst->db = db; dst->img = img; dst->img_mask = mask;
    dst->map = map; dst->map_size = map_size;
    memcpy(dst->map, src->map, (size_t)map_size * sizeof(map[0]));
}

OST_DEF void compute_score(MatchResult *m, const Config *c, MatchWork *w)
{
    /* Start with a false-star penalty, then add the best projected match per star. */
    m->match.totalscore = log(1.0 / (c->IMG_X * c->IMG_Y)) * (2 * m->map_size);
    for (int i = 0; i < m->map_size; i++) {
        m->map[i] = -1;
        w->scores[i] = 0.0f;
    }
    for (int i = 0; i < m->db->results.kdresults_size; i++) {
        Star *s = &m->db->results.map[m->db->results.kdresults[i]];
        int o = s->star_idx;
        float x = ost_vec_dot3(s->v, m->R[0]);
        float y = ost_vec_dot3(s->v, m->R[1]);
        float z = ost_vec_dot3(s->v, m->R[2]);
        float px = y / (x * c->PIXX_TANGENT);
        float py = z / (x * c->PIXY_TANGENT);
        int n = fov_get_id(m->img_mask, c, px, py);
        if (n >= 0) {
            float score = fov_score(m->img_mask, n, px, py);
            if (score > w->scores[n]) {
                m->map[n] = o;
                w->scores[n] = score;
            }
        }
    }
    for (int i = 0; i < m->map_size; i++)
        m->match.totalscore += w->scores[i];
}

OST_DEF int related(MatchResult *winner, CPair *p)
{
    if (winner->match.totalscore == -FLT_MAX || p->totalscore == -FLT_MAX)
        return 0;
    return winner->map[p->img_s1] == p->db_s1 && winner->map[p->img_s2] == p->db_s2;
}

static void ost_sort_float_key(float *v, int n)
{
    for (int i = 1; i < n; i++) {
        float x = v[i];
        int j = i;
        while (j > 0 && x < v[j - 1]) {
            v[j] = v[j - 1];
            j--;
        }
        v[j] = x;
    }
}

static int ost_constellation_pairdist_key(const Star *stars,
                                          const uint32_t *ids, int k,
                                          float *out, int dims)
{
    int q = 0;
    if (!stars || !ids || !out || dims != k * (k - 1) / 2)
        return -1;
    for (int j = 1; j < k; j++)
        for (int i = 0; i < j; i++)
            out[q++] = ost_star_dist_arcsec(&stars[ids[i]], &stars[ids[j]]);
    return 0;
}

typedef struct OSTCpx {
    double r, i;
} OSTCpx;

static int ost_double_finite(double x)
{
    return x == x && x <= DBL_MAX && x >= -DBL_MAX;
}

static OSTCpx ost_cpx_make(double r, double i)
{
    OSTCpx z;
    z.r = r;
    z.i = i;
    return z;
}

static OSTCpx ost_cpx_add(OSTCpx a, OSTCpx b)
{
    return ost_cpx_make(a.r + b.r, a.i + b.i);
}

static OSTCpx ost_cpx_sub(OSTCpx a, OSTCpx b)
{
    return ost_cpx_make(a.r - b.r, a.i - b.i);
}

static OSTCpx ost_cpx_mul(OSTCpx a, OSTCpx b)
{
    return ost_cpx_make(a.r * b.r - a.i * b.i,
                        a.r * b.i + a.i * b.r);
}

static double ost_cpx_abs2(OSTCpx z)
{
    return z.r * z.r + z.i * z.i;
}

static int ost_cpx_div(OSTCpx a, OSTCpx b, OSTCpx *out)
{
    double den = ost_cpx_abs2(b);
    if (!out || den <= 1e-60 || !ost_double_finite(den))
        return -1;
    *out = ost_cpx_make((a.r * b.r + a.i * b.i) / den,
                        (a.i * b.r - a.r * b.i) / den);
    return ost_double_finite(out->r) && ost_double_finite(out->i) ? 0 : -1;
}

static int ost_constellation_crossratio_poly_complex(OSTCpx z, float *out)
{
    OSTCpx one = ost_cpx_make(1.0, 0.0);
    OSTCpx z2 = ost_cpx_mul(z, z);
    OSTCpx a = ost_cpx_add(ost_cpx_sub(z2, z), one);
    OSTCpx num = ost_cpx_mul(ost_cpx_mul(a, a), a);
    OSTCpx zm1 = ost_cpx_sub(z, one);
    OSTCpx den = ost_cpx_mul(z2, ost_cpx_mul(zm1, zm1));
    OSTCpx y;
    if (!out || ost_cpx_div(num, den, &y) < 0 ||
        fabs(y.r) > (double)FLT_MAX || fabs(y.i) > (double)FLT_MAX)
        return -1;
    out[0] = (float)y.r;
    out[1] = (float)y.i;
    return 0;
}

static int ost_constellation_crossratio_poly_real(double z, float *out)
{
    double a = z * z - z + 1.0;
    double den = z * z * (z - 1.0) * (z - 1.0);
    double y;
    if (!out || fabs(den) <= 1e-60 || !ost_double_finite(den))
        return -1;
    y = (a * a * a) / den;
    if (!ost_double_finite(y) || fabs(y) > (double)FLT_MAX)
        return -1;
    *out = (float)y;
    return 0;
}

static void ost_constellation_chart_components(const float v[3], int axis,
                                                double *u, double *a,
                                                double *b)
{
    if (axis == 1) {
        *u = v[1]; *a = v[2]; *b = v[0];
    } else if (axis == 2) {
        *u = v[2]; *a = v[0]; *b = v[1];
    } else {
        *u = v[0]; *a = v[1]; *b = v[2];
    }
}

static int ost_constellation_stereo_axis(const Star *stars,
                                         const uint32_t *ids, int k)
{
    int best_axis = -1;
    double best_min_den = -1.0;
    for (int axis = 0; axis < 3; axis++) {
        double min_den = DBL_MAX;
        for (int i = 0; i < k; i++) {
            double u, a, b, den;
            ost_constellation_chart_components(stars[ids[i]].v, axis,
                                               &u, &a, &b);
            den = 1.0 + u;
            if (den < min_den)
                min_den = den;
        }
        if (min_den > best_min_den) {
            best_min_den = min_den;
            best_axis = axis;
        }
    }
    return best_min_den > 1e-12 ? best_axis : -1;
}

static int ost_constellation_stereo_project(const Star *stars,
                                            const uint32_t *ids, int k,
                                            OSTCpx *z)
{
    int axis;
    if (!stars || !ids || !z)
        return -1;
    axis = ost_constellation_stereo_axis(stars, ids, k);
    if (axis < 0)
        return -1;
    for (int i = 0; i < k; i++) {
        double u, a, b, den;
        ost_constellation_chart_components(stars[ids[i]].v, axis,
                                           &u, &a, &b);
        den = 1.0 + u;
        if (den <= 1e-12)
            return -1;
        z[i] = ost_cpx_make(a / den, b / den);
        if (!ost_double_finite(z[i].r) || !ost_double_finite(z[i].i))
            return -1;
    }
    return 0;
}

static int ost_constellation_crossratio4_raw(const OSTCpx z[4], OSTCpx *cr)
{
    OSTCpx num, den;
    if (!z || !cr)
        return -1;
    num = ost_cpx_mul(ost_cpx_sub(z[0], z[2]), ost_cpx_sub(z[1], z[3]));
    den = ost_cpx_mul(ost_cpx_sub(z[1], z[2]), ost_cpx_sub(z[0], z[3]));
    return ost_cpx_div(num, den, cr);
}

static int ost_constellation_crossratio4_key(const Star *stars,
                                             const uint32_t *ids, int k,
                                             float *out, int dims)
{
    OSTCpx z[4], cr;
    if (!stars || !ids || !out || k != 4 || dims != 2 ||
        ost_constellation_stereo_project(stars, ids, k, z) < 0 ||
        ost_constellation_crossratio4_raw(z, &cr) < 0)
        return -1;
    return ost_constellation_crossratio_poly_complex(cr, out);
}

static double ost_det3_vec(const float a[3], const float b[3],
                           const float c[3])
{
    return (double)a[0] * ((double)b[1] * c[2] - (double)b[2] * c[1]) -
           (double)a[1] * ((double)b[0] * c[2] - (double)b[2] * c[0]) +
           (double)a[2] * ((double)b[0] * c[1] - (double)b[1] * c[0]);
}

static int ost_constellation_crossratio5_raw_ordered(const Star *stars,
        const uint32_t *ids, int i1, int i2, int i3, int i4, int i5,
        double *out, double *quality)
{
    const float *x1 = stars[ids[i1]].v;
    const float *x2 = stars[ids[i2]].v;
    const float *x3 = stars[ids[i3]].v;
    const float *x4 = stars[ids[i4]].v;
    const float *x5 = stars[ids[i5]].v;
    double d143 = ost_det3_vec(x1, x4, x3);
    double d125 = ost_det3_vec(x1, x2, x5);
    double num = ost_det3_vec(x1, x2, x4) * ost_det3_vec(x1, x5, x3);
    double den = d143 * d125;
    if (!out || fabs(den) <= 1e-60 || !ost_double_finite(num) ||
        !ost_double_finite(den))
        return -1;
    *out = num / den;
    if (quality)
        *quality = fabs(den);
    return ost_double_finite(*out) ? 0 : -1;
}

static int ost_constellation_crossratio5_anchor_key(const Star *stars,
        const uint32_t *ids, int anchor, float *out)
{
    int other[4], no = 0;
    double raw;
    double best_q = -1.0;
    float best = 0.0f;
    int have = 0;
    for (int i = 0; i < 5; i++)
        if (i != anchor)
            other[no++] = i;
    if (ost_constellation_crossratio5_raw_ordered(stars, ids, anchor,
            other[0], other[1], other[2], other[3], &raw, NULL) == 0 &&
        ost_constellation_crossratio_poly_real(raw, out) == 0)
        return 0;
    for (int i2 = 0; i2 < 5; i2++) if (i2 != anchor)
        for (int i3 = 0; i3 < 5; i3++) if (i3 != anchor && i3 != i2)
            for (int i4 = 0; i4 < 5; i4++)
                if (i4 != anchor && i4 != i2 && i4 != i3)
                    for (int i5 = 0; i5 < 5; i5++)
                        if (i5 != anchor && i5 != i2 && i5 != i3 &&
                            i5 != i4) {
                            double raw, q;
                            float val;
                            if (ost_constellation_crossratio5_raw_ordered(
                                    stars, ids, anchor, i2, i3, i4, i5,
                                    &raw, &q) < 0 ||
                                ost_constellation_crossratio_poly_real(raw,
                                    &val) < 0)
                                continue;
                            if (!have || q > best_q) {
                                have = 1;
                                best_q = q;
                                best = val;
                            }
                        }
    if (!have || !out)
        return -1;
    *out = best;
    return 0;
}

static int ost_constellation_crossratio5_key(const Star *stars,
                                             const uint32_t *ids, int k,
                                             float *out, int dims)
{
    if (!stars || !ids || !out || k != 5 || dims != 5)
        return -1;
    for (int a = 0; a < 5; a++)
        if (ost_constellation_crossratio5_anchor_key(stars, ids,
                a, &out[a]) < 0)
            return -1;
    ost_sort_float_key(out, 5);
    return 0;
}

static int ost_constellation_key(int descriptor_kind, const Star *stars,
                                 const uint32_t *ids, int k,
                                 float *out, int dims)
{
    switch (descriptor_kind) {
    case OST_CONSTELLATION_PAIRDIST:
        return ost_constellation_pairdist_key(stars, ids, k, out, dims);
    case OST_CONSTELLATION_CROSSRATIO:
        return k == 4 ? ost_constellation_crossratio4_key(stars, ids, k, out, dims) :
                        ost_constellation_crossratio5_key(stars, ids, k, out, dims);
    default:
        return -1;
    }
}

static int ost_constellation_dims(int k, int descriptor_kind)
{
    if (k < 2 || k > OST_MAX_CONSTELLATION_STARS)
        return 0;
    switch (descriptor_kind) {
    case OST_CONSTELLATION_PAIRDIST:
        return k * (k - 1) / 2;
    case OST_CONSTELLATION_CROSSRATIO:
        return k == 4 ? 2 : (k == 5 ? 5 : 0);
    default:
        return 0;
    }
}

static int ost_constellation_descriptor_symmetric(int k, int descriptor_kind)
{
    return (descriptor_kind == OST_CONSTELLATION_PAIRDIST && k == 2) ||
           (descriptor_kind == OST_CONSTELLATION_CROSSRATIO &&
            (k == 4 || k == 5));
}

OST_DEF size_t ost_constellation_record_size(int k, int descriptor_kind)
{
    int dims = ost_constellation_dims(k, descriptor_kind);
    if (dims <= 0)
        return 0;
    return (size_t)dims * sizeof(float) + (size_t)k * sizeof(uint32_t);
}

OST_DEF int ost_constellation_index_init(OSTConstellationIndex *idx,
        CDB *pair, int k, int descriptor_kind, unsigned char *storage,
        int storage_cap)
{
    int dims = ost_constellation_dims(k, descriptor_kind);
    size_t rec = ost_constellation_record_size(k, descriptor_kind);
    if (!idx || !pair || !rec || storage_cap < 0 || dims <= 0 ||
        dims > OST_MAX_CONSTELLATION_DIMS || (!storage && storage_cap > 0))
        return -1;
    memset(idx, 0, sizeof(*idx));
    idx->pair = pair;
    idx->map = storage;
    idx->cap = storage_cap;
    idx->k = k;
    idx->dims = dims;
    idx->descriptor_kind = descriptor_kind;
    idx->record_size = rec;
    idx->kd_bucket = 1 << (idx->k + 1);
    return 0;
}

static unsigned char *ost_constellation_index_ptr(OSTConstellationIndex *idx,
                                                  int row)
{
    return idx->map + (size_t)row * idx->record_size;
}

static float *ost_constellation_index_key(OSTConstellationIndex *idx, int row)
{
    return (float *)ost_constellation_index_ptr(idx, row);
}

static uint32_t *ost_constellation_index_stars(OSTConstellationIndex *idx,
                                               int row)
{
    return (uint32_t *)(ost_constellation_index_ptr(idx, row) +
                       (size_t)idx->dims * sizeof(float));
}

static int ost_constellation_edge_cmp(const void *pa, const void *pb)
{
    const OSTConstellationEdge *a = (const OSTConstellationEdge *)pa;
    const OSTConstellationEdge *b = (const OSTConstellationEdge *)pb;
    return (a->star > b->star) - (a->star < b->star);
}

static int ost_constellation_build_adj(CDB *pair, int *off,
                                       OSTConstellationEdge *edges, int *tmp)
{
    int nstars, npairs;
    if (!pair || !off || !edges || !tmp || pair->stars.n < 0 ||
        pair->map_size < 0)
        return -1;
    nstars = pair->stars.n;
    npairs = pair->map_size;
    for (int i = 0; i <= nstars; i++)
        off[i] = 0;
    for (int i = 0; i < npairs; i++) {
        int a = pair->map[i].s1, b = pair->map[i].s2;
        if ((unsigned)a >= (unsigned)nstars || (unsigned)b >= (unsigned)nstars)
            return -1;
        off[a + 1]++;
        off[b + 1]++;
    }
    for (int i = 0; i < nstars; i++) {
        off[i + 1] += off[i];
        tmp[i] = off[i];
    }
    for (int i = 0; i < npairs; i++) {
        int a = pair->map[i].s1, b = pair->map[i].s2;
        edges[tmp[a]++].star = b;
        edges[tmp[b]++].star = a;
    }
    for (int i = 0; i < nstars; i++)
        qsort(edges + off[i], (size_t)(off[i + 1] - off[i]),
              sizeof(*edges), ost_constellation_edge_cmp);
    return 0;
}

static int ost_constellation_adj_has(const int *off,
                                     const OSTConstellationEdge *edges,
                                     int a, int b)
{
    int l = off[a], r = off[a + 1];
    while (l < r) {
        int m = (l + r) >> 1;
        if (edges[m].star < b) l = m + 1;
        else r = m;
    }
    return l < off[a + 1] && edges[l].star == b;
}

static void ost_constellation_local_basis(const float v[3], float b0[3],
                                          float b1[3])
{
    float ref[3] = {0.0f, 0.0f, 1.0f};
    if (fabsf(v[2]) > 0.9f) {
        ref[1] = 1.0f;
        ref[2] = 0.0f;
    }
    ost_vec_cross3(b0, ref, v);
    if (ost_vec_normalize3(b0) < 0) {
        b0[0] = 1.0f;
        b0[1] = b0[2] = 0.0f;
    }
    ost_vec_cross3(b1, v, b0);
    ost_vec_normalize3(b1);
}

static void ost_constellation_perturb_catalog_star(const Config *cfg,
        const Star *src, const float axis[3], float theta, Star *dst)
{
    float c = cosf(theta), s = sinf(theta);
    *dst = *src;
    dst->v[0] = src->v[0] * c + axis[0] * s;
    dst->v[1] = src->v[1] * c + axis[1] * s;
    dst->v[2] = src->v[2] * c + axis[2] * s;
    ost_vec_normalize3(dst->v);
    if (fabsf(dst->v[0]) > 1e-20f) {
        dst->px = dst->v[1] / (dst->v[0] * cfg->PIXX_TANGENT);
        dst->py = dst->v[2] / (dst->v[0] * cfg->PIXY_TANGENT);
    }
}

static int ost_constellation_qerr_pairdist_fast(
        const Config *cfg, const Star *stars, const uint32_t *ids, int k,
        float db_max_variance, float *qkey, float *qerr)
{
    int q = 0;
    if (!cfg || !stars || !ids || !qkey || !qerr ||
        ost_constellation_pairdist_key(stars, ids, k, qkey,
            k * (k - 1) / 2) < 0)
        return -1;
    for (int j = 1; j < k; j++)
        for (int i = 0; i < j; i++) {
            float v = stars[ids[i]].sigma_sq + stars[ids[j]].sigma_sq +
                      2.0f * db_max_variance;
            qerr[q++] = cfg->POS_ERR_SIGMA * cfg->PIXSCALE *
                        sqrtf(fmaxf(0.0f, v));
        }
    return 0;
}

static int ost_constellation_ut_accum(int descriptor_kind, int k, int dims,
        const Star *stars, const uint32_t *ids, int sign,
        const float *center, float weight, float *sum, float *sumsq)
{
    float key[OST_MAX_CONSTELLATION_DIMS];
    if (ost_constellation_key(descriptor_kind, stars, ids, k, key, dims) < 0)
        return -1;
    for (int d = 0; d < dims; d++) {
        float delta = sign > 0 ? key[d] - center[d] : center[d] - key[d];
        sum[d] += weight * delta;
        sumsq[d] += weight * delta * delta;
    }
    return 0;
}

static int ost_constellation_qerr_ut(
        const Config *cfg, int descriptor_kind,
        const Star *stars, const uint32_t *ids, int k,
        float db_max_variance, float *qkey, float *qerr)
{
    Star base[OST_MAX_CONSTELLATION_STARS];
    Star pert[OST_MAX_CONSTELLATION_STARS];
    uint32_t local_ids[OST_MAX_CONSTELLATION_STARS] = {0};
    int dims = ost_constellation_dims(k, descriptor_kind);
    if (k < 2 || k > OST_MAX_CONSTELLATION_STARS ||
        dims <= 0 || dims > OST_MAX_CONSTELLATION_DIMS)
        return -1;
    if (descriptor_kind == OST_CONSTELLATION_PAIRDIST)
        return ost_constellation_qerr_pairdist_fast(cfg, stars, ids, k,
                                                    db_max_variance,
                                                    qkey, qerr);
    float sum[OST_MAX_CONSTELLATION_DIMS] = {0.0f};
    float sumsq[OST_MAX_CONSTELLATION_DIMS] = {0.0f};
    int nvar = 4 * k;
    float root_n = sqrtf((float)nvar);
    float weight = 1.0f / (float)(2 * nvar);
    float cat_sigma_rad = cfg->PIXSCALE * sqrtf(fmaxf(0.0f, db_max_variance)) /
                          ((float)(3600.0 * 180.0 / PI));

    for (int i = 0; i < k; i++) {
        base[i] = stars[ids[i]];
        local_ids[i] = (uint32_t)i;
    }
    if (ost_constellation_key(descriptor_kind, base, local_ids, k,
                              qkey, dims) < 0)
        return -1;

    for (int s = 0; s < k; s++) {
        float img_sigma = sqrtf(fmaxf(0.0f, base[s].sigma_sq));
        if (img_sigma > 0.0f) {
            float amp = root_n * img_sigma;
            for (int axis = 0; axis < 2; axis++) {
                memcpy(pert, base, (size_t)k * sizeof(base[0]));
                pert[s] = ost_make_img_star(cfg,
                    base[s].px + (axis == 0 ? amp : 0.0f),
                    base[s].py + (axis == 1 ? amp : 0.0f),
                    base[s].flux, base[s].id);
                pert[s].star_idx = base[s].star_idx;
                if (ost_constellation_ut_accum(descriptor_kind, k, dims, pert, local_ids, 1,
                        qkey, weight, sum, sumsq) < 0)
                    return -1;
                memcpy(pert, base, (size_t)k * sizeof(base[0]));
                pert[s] = ost_make_img_star(cfg,
                    base[s].px - (axis == 0 ? amp : 0.0f),
                    base[s].py - (axis == 1 ? amp : 0.0f),
                    base[s].flux, base[s].id);
                pert[s].star_idx = base[s].star_idx;
                if (ost_constellation_ut_accum(descriptor_kind, k, dims, pert, local_ids, 1,
                        qkey, weight, sum, sumsq) < 0)
                    return -1;
            }
        }
        if (cat_sigma_rad > 0.0f) {
            float b0[3], b1[3], amp = root_n * cat_sigma_rad;
            ost_constellation_local_basis(base[s].v, b0, b1);
            for (int axis = 0; axis < 2; axis++) {
                const float *b = axis == 0 ? b0 : b1;
                memcpy(pert, base, (size_t)k * sizeof(base[0]));
                ost_constellation_perturb_catalog_star(cfg, &base[s], b,
                                                       amp, &pert[s]);
                if (ost_constellation_ut_accum(descriptor_kind, k, dims, pert, local_ids, -1,
                        qkey, weight, sum, sumsq) < 0)
                    return -1;
                memcpy(pert, base, (size_t)k * sizeof(base[0]));
                ost_constellation_perturb_catalog_star(cfg, &base[s], b,
                                                       -amp, &pert[s]);
                if (ost_constellation_ut_accum(descriptor_kind, k, dims, pert, local_ids, -1,
                        qkey, weight, sum, sumsq) < 0)
                    return -1;
            }
        }
    }
    for (int d = 0; d < dims; d++)
        qerr[d] = cfg->POS_ERR_SIGMA *
            sqrtf(fmaxf(0.0f, sumsq[d] - sum[d] * sum[d]));
    return 0;
}

typedef struct OSTConstellationCliqueCtx {
    OSTConstellationIndex *out;
    const int *off;
    const OSTConstellationEdge *edges;
    int *common;
    int common_count, k;
    uint32_t chosen[OST_MAX_CONSTELLATION_STARS];
    uint64_t count;
} OSTConstellationCliqueCtx;

static int ost_constellation_add(OSTConstellationIndex *idx,
                                 const uint32_t *ids)
{
    uint32_t ordered[OST_MAX_CONSTELLATION_STARS];
    memcpy(ordered, ids, (size_t)idx->k * sizeof(ordered[0]));
    ost_constellation_canonicalize_ids(idx->pair->stars.v, ordered, idx->k);
    if (idx->map_size >= idx->cap)
        return -2;
    if (ost_constellation_key(idx->descriptor_kind, idx->pair->stars.v,
            ordered, idx->k, ost_constellation_index_key(idx, idx->map_size),
            idx->dims) < 0)
        return 0;
    memcpy(ost_constellation_index_stars(idx, idx->map_size), ordered,
           (size_t)idx->k * sizeof(ordered[0]));
    idx->map_size++;
    return 0;
}

static int ost_constellation_clique_rec(OSTConstellationCliqueCtx *c,
                                        int start, int depth)
{
    if (depth == c->k) {
        c->count++;
        return c->out ? ost_constellation_add(c->out, c->chosen) : 0;
    }
    if (depth >= OST_MAX_CONSTELLATION_STARS)
        return -1;
    for (int p = start; p <= c->common_count - (c->k - depth); p++) {
        int star = c->common[p], ok = 1;
        for (int i = 2; i < depth; i++)
            if (!ost_constellation_adj_has(c->off, c->edges, star,
                                           (int)c->chosen[i])) {
                ok = 0;
                break;
            }
        if (!ok)
            continue;
        c->chosen[depth] = (uint32_t)star;
        ok = ost_constellation_clique_rec(c, p + 1, depth + 1);
        if (ok < 0)
            return ok;
    }
    return 0;
}

static int ost_constellation_for_each_clique(CDB *pair, int k,
                                             OSTConstellationIndex *out,
                                             uint64_t *count, int *off,
                                             OSTConstellationEdge *edges,
                                             int *tmp, int *common)
{
    if (!pair || k < 2 || k > OST_MAX_CONSTELLATION_STARS ||
        pair->stars.n < k || pair->map_size < 0) {
        if (count) *count = 0;
        return pair && pair->stars.n < k ? 0 : -1;
    }
    if (count)
        *count = 0;
    if (k > 2 && (!off || !edges || !tmp || !common ||
        ost_constellation_build_adj(pair, off, edges, tmp) < 0))
        return -1;
    for (int i = 0; i < pair->map_size; i++) {
        int a = pair->map[i].s1, b = pair->map[i].s2;
        if (k > 2 && ost_constellation_before(pair->stars.v,
                (uint32_t)b, (uint32_t)a))
            OST_SWAP(int, a, b);
        OSTConstellationCliqueCtx c;
        memset(&c, 0, sizeof(c));
        c.out = out;
        c.off = off;
        c.edges = edges;
        c.common = common;
        c.k = k;
        c.chosen[0] = (uint32_t)a;
        c.chosen[1] = (uint32_t)b;
        if (k > 2) {
            int ia = off[a], ea = off[a + 1], ib = off[b], eb = off[b + 1];
            while (ia < ea && ib < eb) {
                int ca = edges[ia].star, cb = edges[ib].star;
                if (ca == cb) {
                    if (ost_constellation_before(pair->stars.v,
                            (uint32_t)b, (uint32_t)ca))
                        common[c.common_count++] = ca;
                    ia++;
                    ib++;
                } else if (ca < cb) {
                    ia++;
                } else {
                    ib++;
                }
            }
        }
        if (ost_constellation_clique_rec(&c, 0, 2) < 0)
            return -2;
        if (count)
            *count += c.count;
    }
    return 0;
}

OST_DEF int ost_constellation_count(CDB *pair, int k, uint64_t *count,
                                    int *off, OSTConstellationEdge *edges,
                                    int *tmp, int *common)
{
    return count ? ost_constellation_for_each_clique(pair, k, NULL, count,
                                                     off, edges, tmp, common)
                 : -1;
}

OST_DEF int ost_constellation_index_build(OSTConstellationIndex *idx, int *off,
                                          OSTConstellationEdge *edges, int *tmp,
                                          int *common)
{
    uint64_t unused = 0;
    if (!idx || !idx->pair || (!idx->map && idx->cap > 0))
        return -1;
    idx->map_size = 0;
    idx->kd_ready = 0;
    return ost_constellation_for_each_clique(idx->pair, idx->k, idx, &unused,
                                             off, edges, tmp, common);
}

OST_DEF void ost_constellation_index_kdsort(OSTConstellationIndex *idx)
{
    if (idx && !idx->kd_ready) {
        ost_kdbuild(idx->map, idx->record_size, 0, 0, idx->map_size,
                    idx->kd_bucket, 0, idx->dims);
        idx->kd_ready = 1;
    }
}

typedef struct OSTConstellationSearchCtx {
    OSTConstellationIndex *cat;
    CDB *img;
    MatchResult *m, *winner;
    const Config *cfg;
    MatchWork *w;
    CPair *candidates;
    int *nc;
    int candidate_cap;
    float qkey[OST_MAX_CONSTELLATION_DIMS];
    float qerr[OST_MAX_CONSTELLATION_DIMS];
    int img_perm_count, pair_dims;
    int img_perms[OST_MAX_CONSTELLATION_PERMS][OST_MAX_CONSTELLATION_STARS];
    float img_pair_key[OST_MAX_CONSTELLATION_PERMS][OST_MAX_CONSTELLATION_DIMS];
    float img_pair_err[OST_MAX_CONSTELLATION_PERMS][OST_MAX_CONSTELLATION_DIMS];
} OSTConstellationSearchCtx;

static int ost_constellation_pairwise_ok(OSTConstellationSearchCtx *ctx,
                                         const float *db_pair_key,
                                         int perm_slot)
{
    if (ctx->cat->descriptor_kind == OST_CONSTELLATION_PAIRDIST)
        return 1;
    if (!db_pair_key || perm_slot < 0 || perm_slot >= ctx->img_perm_count ||
        ctx->pair_dims <= 0 || ctx->pair_dims > OST_MAX_CONSTELLATION_DIMS)
        return 0;
    for (int d = 0; d < ctx->pair_dims; d++)
        if (fabsf(db_pair_key[d] - ctx->img_pair_key[perm_slot][d]) >
            ctx->img_pair_err[perm_slot][d])
            return 0;
    return 1;
}

static int ost_constellation_score_candidate(OSTConstellationSearchCtx *ctx,
                                             const int *db_idx,
                                             const int *img_idx)
{
    OSTConstellationIndex *cat = ctx->cat;
    ctx->m->match.db_s1 = db_idx[0];
    ctx->m->match.db_s2 = db_idx[1];
    ctx->m->match.img_s1 = img_idx[0];
    ctx->m->match.img_s2 = img_idx[1];
    if (weighted_wahba_vectors(ctx->m->R, cat->pair->stars.v,
                               ctx->img->stars.v, db_idx, img_idx,
                               cat->k, QMETHOD_ITER) < 0)
        return 0;
    if (cat->pair->results.kdsorted)
        ost_query_search(&cat->pair->results, ctx->cfg, ctx->m->R[0],
                         ctx->cfg->MAXFOV / 2.0f,
                         ctx->cfg->THRESH_FACTOR * ctx->cfg->IMAGE_VARIANCE);
    compute_score(ctx->m, ctx->cfg, ctx->w);
    if (ctx->m->match.totalscore > ctx->winner->match.totalscore) {
        if (ctx->winner->match.totalscore != -FLT_MAX) {
            if (*ctx->nc >= ctx->candidate_cap) {
                if (cat->pair->results.kdsorted)
                    ost_query_clear_results(&cat->pair->results);
                return -1;
            }
            ctx->candidates[(*ctx->nc)++] = ctx->winner->match;
        }
        mr_copy(ctx->winner, ctx->m);
    } else {
        if (*ctx->nc >= ctx->candidate_cap) {
            if (cat->pair->results.kdsorted)
                ost_query_clear_results(&cat->pair->results);
            return -1;
        }
        ctx->candidates[(*ctx->nc)++] = ctx->m->match;
    }
    if (cat->pair->results.kdsorted)
        ost_query_clear_results(&cat->pair->results);
    return 0;
}

static int ost_constellation_visit(void *base, int row, void *vctx)
{
    OSTConstellationSearchCtx *ctx = (OSTConstellationSearchCtx *)vctx;
    OSTConstellationIndex *cat = ctx->cat;
    float *key = (float *)((unsigned char *)base +
                 (size_t)row * cat->record_size);
    uint32_t *ids = (uint32_t *)((unsigned char *)key +
                    (size_t)cat->dims * sizeof(float));
    int db_idx[OST_MAX_CONSTELLATION_STARS];
    float db_pair_key[OST_MAX_CONSTELLATION_DIMS];
    for (int i = 0; i < cat->k; i++)
        db_idx[i] = (int)ids[i];
    if (cat->descriptor_kind != OST_CONSTELLATION_PAIRDIST &&
        ost_constellation_pairdist_key(cat->pair->stars.v, ids, cat->k,
            db_pair_key, ctx->pair_dims) < 0)
        return 0;
    for (int p = 0; p < ctx->img_perm_count; p++) {
        int rc;
        if (!ost_constellation_pairwise_ok(ctx, db_pair_key, p))
            continue;
        rc = ost_constellation_score_candidate(ctx, db_idx,
                                               ctx->img_perms[p]);
        if (rc < 0)
            return rc;
    }
    return 0;
}

typedef struct OSTConstellationQueryVariant {
    float key[OST_MAX_CONSTELLATION_DIMS];
    float qerr[OST_MAX_CONSTELLATION_DIMS];
    float pair_key[OST_MAX_CONSTELLATION_DIMS];
    float pair_err[OST_MAX_CONSTELLATION_DIMS];
    int img_idx[OST_MAX_CONSTELLATION_STARS];
} OSTConstellationQueryVariant;

static int ost_constellation_variant_same(
        const OSTConstellationIndex *cat,
        const OSTConstellationQueryVariant *a,
        const OSTConstellationQueryVariant *b)
{
    for (int d = 0; d < cat->dims; d++)
        if (fabsf(a->key[d] - b->key[d]) > 1e-5f ||
            fabsf(a->qerr[d] - b->qerr[d]) > 1e-5f)
            return 0;
    return 1;
}

static int ost_constellation_search_variants(OSTConstellationSearchCtx *s,
        OSTConstellationQueryVariant *variants, int variant_count)
{
    int used[OST_MAX_CONSTELLATION_PERMS] = {0};
    s->pair_dims = ost_constellation_dims(s->cat->k,
                                          OST_CONSTELLATION_PAIRDIST);
    if (variant_count > OST_MAX_CONSTELLATION_PERMS)
        return -1;
    for (int i = 0; i < variant_count; i++) if (!used[i]) {
        int rc;
        memcpy(s->qkey, variants[i].key,
               (size_t)s->cat->dims * sizeof(s->qkey[0]));
        memcpy(s->qerr, variants[i].qerr,
               (size_t)s->cat->dims * sizeof(s->qerr[0]));
        s->img_perm_count = 0;
        for (int j = i; j < variant_count; j++) {
            if (used[j] || !ost_constellation_variant_same(s->cat,
                                                           &variants[i],
                                                           &variants[j]))
                continue;
            if (s->img_perm_count >= OST_MAX_CONSTELLATION_PERMS)
                return -1;
            memcpy(s->img_perms[s->img_perm_count], variants[j].img_idx,
                   (size_t)s->cat->k * sizeof(s->img_perms[0][0]));
            if (s->cat->descriptor_kind != OST_CONSTELLATION_PAIRDIST) {
                memcpy(s->img_pair_key[s->img_perm_count], variants[j].pair_key,
                       (size_t)s->pair_dims * sizeof(s->img_pair_key[0][0]));
                memcpy(s->img_pair_err[s->img_perm_count], variants[j].pair_err,
                       (size_t)s->pair_dims * sizeof(s->img_pair_err[0][0]));
            }
            s->img_perm_count++;
            used[j] = 1;
        }
        rc = ost_kdsearch(s->cat->map, s->cat->record_size, 0,
                          0, s->cat->map_size, s->cat->kd_bucket,
                          0, s->cat->dims, s->qkey, s->qerr,
                          ost_constellation_visit, s);
        if (rc < 0)
            return rc;
    }
    s->img_perm_count = 0;
    return 0;
}

typedef struct OSTConstellationPermCtx {
    OSTConstellationSearchCtx *search;
    OSTConstellationQueryVariant *variants;
    int variant_count, variant_cap, fixed_descriptor;
    float fixed_key[OST_MAX_CONSTELLATION_DIMS];
    float fixed_qerr[OST_MAX_CONSTELLATION_DIMS];
    uint32_t ids[OST_MAX_CONSTELLATION_STARS];
    uint32_t perm_ids[OST_MAX_CONSTELLATION_STARS];
    int used[OST_MAX_CONSTELLATION_STARS];
} OSTConstellationPermCtx;

static int ost_constellation_query_perm_rec(OSTConstellationPermCtx *p,
                                            int depth)
{
    OSTConstellationSearchCtx *s = p->search;
    int k = s->cat->k;
    if (depth == k) {
        OSTConstellationQueryVariant *v;
        if (p->variant_count >= p->variant_cap)
            return -1;
        v = &p->variants[p->variant_count];
        if (p->fixed_descriptor) {
            memcpy(v->key, p->fixed_key,
                   (size_t)s->cat->dims * sizeof(v->key[0]));
            memcpy(v->qerr, p->fixed_qerr,
                   (size_t)s->cat->dims * sizeof(v->qerr[0]));
        } else if (ost_constellation_qerr_ut(s->cfg, s->cat->descriptor_kind,
                s->img->stars.v, p->perm_ids, k,
                s->cat->pair->stars.max_variance, v->key, v->qerr) < 0)
            return 0;
        if (s->cat->descriptor_kind != OST_CONSTELLATION_PAIRDIST &&
            ost_constellation_qerr_ut(s->cfg, OST_CONSTELLATION_PAIRDIST,
                s->img->stars.v, p->perm_ids, k,
                s->cat->pair->stars.max_variance,
                v->pair_key, v->pair_err) < 0)
            return 0;
        for (int i = 0; i < k; i++)
            v->img_idx[i] = (int)p->perm_ids[i];
        p->variant_count++;
        return 0;
    }
    for (int i = 0; i < k; i++) if (!p->used[i]) {
        int rc;
        p->used[i] = 1;
        p->perm_ids[depth] = p->ids[i];
        rc = ost_constellation_query_perm_rec(p, depth + 1);
        p->used[i] = 0;
        if (rc < 0)
            return rc;
    }
    return 0;
}

static int ost_constellation_query_feature(OSTConstellationSearchCtx *s,
                                           const uint32_t *ids)
{
    OSTConstellationQueryVariant variants[OST_MAX_CONSTELLATION_PERMS];
    OSTConstellationPermCtx p;
    memset(&p, 0, sizeof(p));
    p.search = s;
    p.variants = variants;
    p.variant_cap = OST_MAX_CONSTELLATION_PERMS;
    if (ost_constellation_descriptor_symmetric(s->cat->k,
                                               s->cat->descriptor_kind)) {
        if (ost_constellation_qerr_ut(s->cfg, s->cat->descriptor_kind,
                s->img->stars.v, ids, s->cat->k,
                s->cat->pair->stars.max_variance,
                p.fixed_key, p.fixed_qerr) < 0)
            return 0;
        p.fixed_descriptor = 1;
    }
    memcpy(p.ids, ids, (size_t)s->cat->k * sizeof(p.ids[0]));
    if (ost_constellation_query_perm_rec(&p, 0) < 0)
        return -1;
    if (p.variant_count <= 0)
        return 0;
    return ost_constellation_search_variants(s, variants, p.variant_count);
}

typedef struct OSTConstellationFeatureCtx {
    OSTConstellationSearchCtx *search;
    int anchor0, anchor1;
    uint32_t extra[OST_MAX_CONSTELLATION_STARS];
} OSTConstellationFeatureCtx;

static int ost_constellation_feature_rec(OSTConstellationFeatureCtx *f,
                                         int start, int depth)
{
    int need = f->search->cat->k - 2;
    if (depth == need) {
        uint32_t ids[OST_MAX_CONSTELLATION_STARS];
        ids[0] = (uint32_t)f->anchor0;
        ids[1] = (uint32_t)f->anchor1;
        memcpy(ids + 2, f->extra, (size_t)need * sizeof(ids[0]));
        ost_constellation_canonicalize_ids(f->search->img->stars.v, ids,
                                           f->search->cat->k);
        if (!((ids[0] == (uint32_t)f->anchor0 &&
               ids[1] == (uint32_t)f->anchor1) ||
              (ids[0] == (uint32_t)f->anchor1 &&
               ids[1] == (uint32_t)f->anchor0)))
            return 0;
        return ost_constellation_query_feature(f->search, ids);
    }
    for (int i = start; i <= f->search->img->stars.n - (need - depth); i++) {
        if (i == f->anchor0 || i == f->anchor1)
            continue;
        f->extra[depth] = (uint32_t)i;
        if (ost_constellation_feature_rec(f, i + 1, depth + 1) < 0)
            return -1;
    }
    return 0;
}

OST_DEF int ost_db_match_constellations(OSTConstellationIndex *cat, CDB *img,
                                        MatchResult *winner,
                                        const Config *cfg, MatchWork *w,
                                        float *p_match)
{
    StarFov fov;
    MatchResult m;
    CPair *candidates = w->candidates;
    int candidate_cap = w->candidate_cap, nc = 0;

    if (!cat || !cat->pair || !img || !winner || !cfg || !w || !p_match)
        return -1;
    *p_match = 0.0f;
    mr_init(winner, cat->pair, img, &fov, w->match_map);
    if (cat->pair->stars.n < cat->k || img->stars.n < cat->k ||
        cat->map_size <= 0)
        return 0;
    if (fov_init(&fov, &img->stars, cat->pair->stars.max_variance, cfg, w) < 0)
        return -1;
    ost_constellation_index_kdsort(cat);
    mr_init(&m, cat->pair, img, &fov, w->work_map);
    for (int n = 0; n < img->map_size; n++) {
        Constellation ic = img->map[n];
        OSTConstellationSearchCtx search;
        OSTConstellationFeatureCtx feat;
        if ((unsigned)ic.s1 >= (unsigned)img->stars.n ||
            (unsigned)ic.s2 >= (unsigned)img->stars.n)
            continue;
        memset(&search, 0, sizeof(search));
        search.cat = cat;
        search.img = img;
        search.m = &m;
        search.winner = winner;
        search.cfg = cfg;
        search.w = w;
        search.candidates = candidates;
        search.nc = &nc;
        search.candidate_cap = candidate_cap;
        memset(&feat, 0, sizeof(feat));
        feat.search = &search;
        feat.anchor0 = ic.s1;
        feat.anchor1 = ic.s2;
        if (ost_constellation_feature_rec(&feat, 0, 0) < 0)
            return -1;
    }
    if (winner->match.totalscore != -FLT_MAX) {
        double p = 1.0;
        for (int i = 0; i < nc; i++)
            if (!related(winner, &candidates[i]))
                p += exp((double)candidates[i].totalscore -
                         winner->match.totalscore);
        *p_match = (float)(1.0 / p);
    }
    return 0;
}

OST_DEF void ost_match_work_init(MatchWork *mw,
                                    CPair *candidates, int candidate_cap,
                                    int *fov_mask, int *collision,
                                    int collision_cap, float *fov_px,
                                    float *fov_py, float *scores,
                                    int *match_map, int *work_map)
{
    mw->candidates = candidates;
    mw->candidate_cap = candidate_cap;
    mw->fov_mask = fov_mask;
    mw->collision = collision;
    mw->collision_cap = collision_cap;
    mw->fov_px = fov_px;
    mw->fov_py = fov_py;
    mw->scores = scores;
    mw->match_map = match_map;
    mw->work_map = work_map;
}

OST_DEF int ost_copy_n_brightest(StarDB *dst, StarDB *src, Star *tmp, int n)
{
    memcpy(tmp, src->v, (size_t)src->n * sizeof(tmp[0]));
    qsort(tmp, (size_t)src->n, sizeof(tmp[0]), cmp_flux_desc);
    dst->n = 0;
    if (n > src->n)
        n = src->n;
    for (int i = 0; i < n; i++)
        if (ost_db_add(dst, tmp[i]) < 0)
            return -1;
    return 0;
}

#endif

#endif
