#ifndef BEAST_H
#define BEAST_H

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef BEAST_RESTRICT
#ifdef __cplusplus
#define BEAST_RESTRICT __restrict__
#else
#define BEAST_RESTRICT restrict
#endif
#endif

#ifndef BEAST_DEF
#ifdef BEAST_EXPORT
#ifdef _WIN32
#define BEAST_DEF __declspec(dllexport)
#else
#define BEAST_DEF __attribute__((visibility("default")))
#endif
#else
#define BEAST_DEF static
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define OST_TRACKER_MAX_STARS 1000
#define BEAST_SWAP(type, a, b) do { \
    type beast_swap_tmp = (a); \
    (a) = (b); \
    (b) = beast_swap_tmp; \
} while (0)

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

typedef struct OSTCCRun {
    int x0;
    int x1;
    int label;
} OSTCCRun;

typedef struct OSTCCBufferSizes {
    int max_labels;
    int max_runs;
    size_t components;
    size_t parent;
    size_t label_live;
    size_t free_after_row;
    size_t touched_stamp;
    size_t seen_stamp;
    size_t prev_runs;
    size_t curr_runs;
    size_t total_bytes;
} OSTCCBufferSizes;

typedef struct OSTCCContext {
    int width;
    int max_labels;
    int max_runs;
    OSTCCComponent *components;
    int *parent;
    int *label_live;
    int *free_after_row;
    int *touched_stamp;
    int *seen_stamp;
    OSTCCRun *prev_runs;
    OSTCCRun *curr_runs;
} OSTCCContext;

int ost_cc_buffer_sizes(int width, OSTCCBufferSizes *sizes);
int ost_cc_init(OSTCCContext *ctx, int width,
                OSTCCComponent *components, int *parent,
                int *label_live, int *free_after_row,
                int *touched_stamp, int *seen_stamp,
                OSTCCRun *prev_runs, OSTCCRun *curr_runs);
int ost_cc_threshold_4(const unsigned char *image,
                       int width, int height, int stride,
                       unsigned char threshold,
                       OSTCCComponent *out, int out_max,
                       OSTCCContext *ctx);

typedef struct OSTBGConfig {
    int width;
    int height;
    int tile_size;
    int map_width;
    int map_height;
    int max_stars;
    int max_pixel_brightness;
    int sample_radius;
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
    double *sigma_col;
    double *rhs;
    double *sigma_solve;
    double *rhs_solve;
    double *cov_xy;
    double *dropped;
} OSTBGFitWorkspace;

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

#ifdef __cplusplus
}
#endif

#ifdef BEAST_IMPLEMENTATION

/* connected components */
#include <math.h>
#include <string.h>

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define OST_CC_RESTRICT restrict
#else
#define OST_CC_RESTRICT
#endif

#if defined(__GNUC__)
#define OST_CC_NOINLINE __attribute__((noinline))
#else
#define OST_CC_NOINLINE
#endif

typedef struct OSTCCBinaryComponent {
    int area;
    int sum_x;
    int sum_y;
} OSTCCBinaryComponent;

typedef double (*OSTCCBackgroundFn)(void *opaque, int x, int y, double *var);
typedef int (*OSTCCRunFn)(void *opaque, const unsigned short *row,
                          int y, int width, OSTCCRun *runs);

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

static void component_add_weighted_run(OSTCCComponent *c,
                                       const unsigned short *image,
                                       int stride, int x0, int x1, int y,
                                       OSTCCBackgroundFn background,
                                       void *background_opaque,
                                       double signal_sigma)
{
    int len;

    len = x1 - x0 + 1;
    if (c->area <= 0) {
        c->min_x = x0;
        c->max_x = x1;
        c->min_y = y;
        c->max_y = y;
    } else {
        if (x0 < c->min_x) c->min_x = x0;
        if (x1 > c->max_x) c->max_x = x1;
        if (y < c->min_y) c->min_y = y;
        if (y > c->max_y) c->max_y = y;
    }
    c->area += len;
    c->sum_x += (x0 + x1) * len / 2;
    c->sum_y += y * len;

    for (int x = x0; x <= x1; x++) {
        double var;
        double mu;
        double v;

        mu = background(background_opaque, x, y, &var);
        v = (double)image[(size_t)y * (size_t)stride + (size_t)x] - mu;
        if (v > signal_sigma * sqrt(var))
            c->signal = 1;
        if (v <= 0)
            continue;
        c->wsum += v;
        c->wx += v * x;
        c->wy += v * y;
        c->wxx += v * x * x;
        c->wyy += v * y * y;
        c->wxy += v * x * y;
    }
}

static int component_finish(OSTCCComponent *c)
{
    double cx;
    double cy;
    double u20;
    double u02;
    double u11;
    double tr;
    double det;
    double d;

    if (c->area <= 0 || c->wsum <= 0 || !c->signal ||
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
        BEAST_SWAP(OSTCCComponent, out[i], out[i - 1]);
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

static int alloc_label(OSTCCContext *ctx, int row, int *next_free_label)
{
    int start;
    int label;

    start = *next_free_label;
    label = start;
    do {
        if (!ctx->label_live[label] && ctx->free_after_row[label] < row) {
            ctx->label_live[label] = 1;
            ctx->parent[label] = label;
            component_clear(&ctx->components[label]);

            label++;
            if (label >= ctx->max_labels)
                label = 1;
            *next_free_label = label;
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
    if (!ra)
        return rb;
    if (!rb)
        return ra;
    if (ra == rb)
        return ra;

    keep = (ra < rb) ? ra : rb;
    merge = (ra < rb) ? rb : ra;

    component_merge(&ctx->components[keep], &ctx->components[merge]);
    ctx->parent[merge] = keep;
    component_clear(&ctx->components[merge]);
    ctx->label_live[merge] = 0;
    ctx->free_after_row[merge] = row + 1;

    return keep;
}

static int alloc_binary_label(OSTCCBinaryComponent *OST_CC_RESTRICT components,
                              int *OST_CC_RESTRICT parent,
                              int *OST_CC_RESTRICT label_live,
                              int *OST_CC_RESTRICT free_after_row,
                              int max_labels, int row,
                              int *next_free_label)
{
    int start;
    int label;

    start = *next_free_label;
    label = start;
    do {
        if (!label_live[label] && free_after_row[label] < row) {
            label_live[label] = 1;
            parent[label] = label;
            components[label].area = 0;
            components[label].sum_x = 0;
            components[label].sum_y = 0;

            label++;
            if (label >= max_labels)
                label = 1;
            *next_free_label = label;
            return label == 1 ? max_labels - 1 : label - 1;
        }

        label++;
        if (label >= max_labels)
            label = 1;
    } while (label != start);

    return 0;
}

static int merge_binary_roots(OSTCCBinaryComponent *OST_CC_RESTRICT components,
                              int *OST_CC_RESTRICT parent,
                              int *OST_CC_RESTRICT label_live,
                              int *OST_CC_RESTRICT free_after_row,
                              int a, int b, int row)
{
    int ra;
    int rb;
    int keep;
    int merge;

    ra = root_compress(parent, a);
    rb = root_compress(parent, b);
    if (ra == rb)
        return ra;

    keep = (ra < rb) ? ra : rb;
    merge = (ra < rb) ? rb : ra;

    components[keep].area += components[merge].area;
    components[keep].sum_x += components[merge].sum_x;
    components[keep].sum_y += components[merge].sum_y;
    parent[merge] = keep;
    components[merge].area = 0;
    components[merge].sum_x = 0;
    components[merge].sum_y = 0;
    label_live[merge] = 0;
    free_after_row[merge] = row + 1;

    return keep;
}

static OST_CC_NOINLINE int extract_runs_thresh(const unsigned char *OST_CC_RESTRICT row,
                                               int width, unsigned char threshold,
                                               OSTCCRun *OST_CC_RESTRICT runs)
{
    int nr;
    int x;

    nr = 0;
    x = 0;
    while (x < width) {
        while (x < width && row[x] <= threshold)
            x++;
        if (x >= width)
            break;

        runs[nr].x0 = x;
        while (x + 1 < width && row[x + 1] > threshold)
            x++;
        runs[nr].x1 = x;
        runs[nr].label = 0;
        nr++;
        x++;
    }
    return nr;
}

int ost_cc_buffer_sizes(int width, OSTCCBufferSizes *sizes)
{
    int max_labels;
    int max_runs;

    if (!sizes || width <= 0)
        return -1;

    max_labels = width + 1;
    max_runs = (width + 1) / 2;

    sizes->max_labels = max_labels;
    sizes->max_runs = max_runs;
    sizes->components = (size_t)max_labels;
    sizes->parent = (size_t)max_labels;
    sizes->label_live = (size_t)max_labels;
    sizes->free_after_row = (size_t)max_labels;
    sizes->touched_stamp = (size_t)max_labels;
    sizes->seen_stamp = (size_t)max_labels;
    sizes->prev_runs = (size_t)max_runs;
    sizes->curr_runs = (size_t)max_runs;
    sizes->total_bytes =
        sizes->components * sizeof(OSTCCComponent) +
        (sizes->parent + sizes->label_live + sizes->free_after_row +
         sizes->touched_stamp + sizes->seen_stamp) * sizeof(int) +
        (sizes->prev_runs + sizes->curr_runs) * sizeof(OSTCCRun);

    return 0;
}

int ost_cc_init(OSTCCContext *ctx, int width,
                OSTCCComponent *components,
                int *parent,
                int *label_live,
                int *free_after_row,
                int *touched_stamp,
                int *seen_stamp,
                OSTCCRun *prev_runs,
                OSTCCRun *curr_runs)
{
    OSTCCBufferSizes sizes;

    if (!ctx || !components || !parent || !label_live || !free_after_row ||
        !touched_stamp || !seen_stamp || !prev_runs || !curr_runs)
        return -1;
    if (ost_cc_buffer_sizes(width, &sizes) < 0)
        return -1;

    ctx->width = width;
    ctx->max_labels = sizes.max_labels;
    ctx->max_runs = sizes.max_runs;
    ctx->components = components;
    ctx->parent = parent;
    ctx->label_live = label_live;
    ctx->free_after_row = free_after_row;
    ctx->touched_stamp = touched_stamp;
    ctx->seen_stamp = seen_stamp;
    ctx->prev_runs = prev_runs;
    ctx->curr_runs = curr_runs;
    return 0;
}

int ost_cc_threshold_4(const unsigned char *image,
                       int width, int height, int stride,
                       unsigned char threshold,
                       OSTCCComponent *out, int out_max,
                       OSTCCContext *ctx)
{
    OSTCCBinaryComponent *components;
    int *parent;
    int *label_live;
    int *free_after_row;
    int *touched_stamp;
    int *seen_stamp;
    OSTCCRun *prev_runs;
    OSTCCRun *curr_runs;
    int count;
    int next_free_label;
    int nr_prev;

    if (!image || !ctx || !out || width <= 0 || height < 0 ||
        stride < width || out_max < 0 || ctx->width != width)
        return -1;

    components = (OSTCCBinaryComponent *)ctx->components;
    parent = ctx->parent;
    label_live = ctx->label_live;
    free_after_row = ctx->free_after_row;
    touched_stamp = ctx->touched_stamp;
    seen_stamp = ctx->seen_stamp;
    prev_runs = ctx->prev_runs;
    curr_runs = ctx->curr_runs;
    count = 0;
    next_free_label = 1;
    nr_prev = 0;

    memset(label_live, 0, (size_t)ctx->max_labels * sizeof(int));
    memset(free_after_row, -1, (size_t)ctx->max_labels * sizeof(int));
    memset(touched_stamp, 0, (size_t)ctx->max_labels * sizeof(int));
    memset(seen_stamp, 0, (size_t)ctx->max_labels * sizeof(int));

    for (int y = 0; y < height; y++) {
        const unsigned char *row;
        int stamp;
        int nr_curr;
        int p;

        row = image + (size_t)y * (size_t)stride;
        stamp = y + 1;
        nr_curr = extract_runs_thresh(row, width, threshold, curr_runs);
        p = 0;

        for (int j = 0; j < nr_curr; j++) {
            OSTCCRun *cr;
            int assigned;
            int len;

            cr = &curr_runs[j];
            assigned = 0;
            len = cr->x1 - cr->x0 + 1;

            while (p < nr_prev && prev_runs[p].x1 < cr->x0)
                p++;

            for (int q = p; q < nr_prev && prev_runs[q].x0 <= cr->x1; q++) {
                int root;

                root = root_compress(parent, prev_runs[q].label);
                touched_stamp[root] = stamp;
                if (!assigned)
                    assigned = root;
                else
                    assigned = merge_binary_roots(components, parent,
                                                  label_live, free_after_row,
                                                  assigned, root, y);
            }

            if (!assigned) {
                assigned = alloc_binary_label(components, parent, label_live,
                                              free_after_row, ctx->max_labels,
                                              y, &next_free_label);
                if (!assigned) {
                    ctx->prev_runs = prev_runs;
                    ctx->curr_runs = curr_runs;
                    return -2;
                }
            }

            cr->label = assigned;
            label_live[assigned] = 1;
            components[assigned].area += len;
            components[assigned].sum_x += (cr->x0 + cr->x1) * len / 2;
            components[assigned].sum_y += y * len;
        }

        for (int i = 0; i < nr_prev; i++) {
            int root;

            root = root_compress(parent, prev_runs[i].label);
            if (seen_stamp[root] == stamp)
                continue;
            seen_stamp[root] = stamp;

            if (touched_stamp[root] != stamp) {
                insert_component(out, &count, out_max,
                                 components[root].area,
                                 components[root].sum_x,
                                 components[root].sum_y);
                components[root].area = 0;
                components[root].sum_x = 0;
                components[root].sum_y = 0;
                label_live[root] = 0;
                free_after_row[root] = y;
                parent[root] = root;
            }
        }

        BEAST_SWAP(OSTCCRun *, prev_runs, curr_runs);
        nr_prev = nr_curr;
    }

    {
        int stamp;

        stamp = height + 1;
        for (int i = 0; i < nr_prev; i++) {
            int root;

            root = root_compress(parent, prev_runs[i].label);
            if (seen_stamp[root] == stamp)
                continue;
            seen_stamp[root] = stamp;

            insert_component(out, &count, out_max,
                             components[root].area,
                             components[root].sum_x,
                             components[root].sum_y);
            components[root].area = 0;
            components[root].sum_x = 0;
            components[root].sum_y = 0;
            label_live[root] = 0;
            free_after_row[root] = height;
            parent[root] = root;
        }
    }

    ctx->prev_runs = prev_runs;
    ctx->curr_runs = curr_runs;
    return count;
}

static int cc_weighted_core(const unsigned short *image,
                            int width, int height, int stride,
                            OSTCCRunFn extract_runs, void *extract_opaque,
                            OSTCCBackgroundFn background,
                            void *background_opaque,
                            double signal_sigma,
                            OSTCCComponent *out, int out_max,
                            OSTCCContext *ctx)
{
    int count;
    int next_free_label;
    int nr_prev;
    int i;
    int y;

    if (!extract_runs || !image || !background || !ctx || !out || width <= 0 ||
        height < 0 || stride < width || out_max < 0 || ctx->width != width)
        return -1;

    count = 0;
    next_free_label = 1;
    nr_prev = 0;

    clear_components(out, out_max);
    clear_components(ctx->components, ctx->max_labels);
    memset(ctx->label_live, 0, (size_t)ctx->max_labels * sizeof(int));
    memset(ctx->free_after_row, -1, (size_t)ctx->max_labels * sizeof(int));
    memset(ctx->touched_stamp, 0, (size_t)ctx->max_labels * sizeof(int));
    memset(ctx->seen_stamp, 0, (size_t)ctx->max_labels * sizeof(int));

    for (i = 0; i < ctx->max_labels; i++)
        ctx->parent[i] = i;

    for (y = 0; y < height; y++) {
        int stamp;
        int nr_curr;
        int p;
        int j;

        stamp = y + 1;
        nr_curr = extract_runs(extract_opaque,
                               image + (size_t)y * (size_t)stride,
                               y, width, ctx->curr_runs);
        p = 0;

        for (j = 0; j < nr_curr; j++) {
            OSTCCRun *cr;
            int assigned;
            int q;

            cr = &ctx->curr_runs[j];
            assigned = 0;

            while (p < nr_prev && ctx->prev_runs[p].x1 < cr->x0)
                p++;

            for (q = p; q < nr_prev && ctx->prev_runs[q].x0 <= cr->x1; q++) {
                int root;

                root = root_compress(ctx->parent, ctx->prev_runs[q].label);
                ctx->touched_stamp[root] = stamp;

                if (!assigned)
                    assigned = root;
                else
                    assigned = merge_roots(ctx, assigned, root, y);
            }

            if (!assigned) {
                assigned = alloc_label(ctx, y, &next_free_label);
                if (!assigned)
                    return -2;
            }

            assigned = root_compress(ctx->parent, assigned);
            cr->label = assigned;
            ctx->label_live[assigned] = 1;
            component_add_weighted_run(&ctx->components[assigned], image, stride,
                                       cr->x0, cr->x1, y, background,
                                       background_opaque, signal_sigma);
        }

        for (i = 0; i < nr_prev; i++) {
            int root;

            root = root_compress(ctx->parent, ctx->prev_runs[i].label);
            if (ctx->seen_stamp[root] == stamp)
                continue;
            ctx->seen_stamp[root] = stamp;

            if (ctx->touched_stamp[root] != stamp) {
                insert_weighted_component(out, &count, out_max,
                                          ctx->components[root]);
                component_clear(&ctx->components[root]);
                ctx->label_live[root] = 0;
                ctx->free_after_row[root] = y;
                ctx->parent[root] = root;
            }
        }

        BEAST_SWAP(OSTCCRun *, ctx->prev_runs, ctx->curr_runs);
        nr_prev = nr_curr;
    }

    {
        int stamp;

        stamp = height + 1;
        for (i = 0; i < nr_prev; i++) {
            int root;

            root = root_compress(ctx->parent, ctx->prev_runs[i].label);
            if (ctx->seen_stamp[root] == stamp)
                continue;
            ctx->seen_stamp[root] = stamp;

            insert_weighted_component(out, &count, out_max,
                                      ctx->components[root]);
            component_clear(&ctx->components[root]);
            ctx->label_live[root] = 0;
            ctx->free_after_row[root] = height;
            ctx->parent[root] = root;
        }
    }

    return count;
}

static int ost_cc_weighted_runs_4(const unsigned short *image,
                                  int width, int height, int stride,
                                  OSTCCRunFn extract_runs,
                                  void *extract_opaque,
                                  OSTCCBackgroundFn background,
                                  void *background_opaque,
                                  double signal_sigma,
                                  OSTCCComponent *out, int out_max,
                                  OSTCCContext *ctx)
{
    return cc_weighted_core(image, width, height, stride, extract_runs,
                            extract_opaque, background,
                            background_opaque, signal_sigma,
                            out, out_max, ctx);
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
    cfg->max_stars = 1000;
    cfg->max_pixel_brightness = 255 * 4;
    cfg->sample_radius = 7;
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

static int bg_threshold_runs(void *opaque, const unsigned short *row,
                             int y, int width, OSTCCRun *runs)
{
    OSTBGThresholdTest *t = (OSTBGThresholdTest *)opaque;
    double gy = y * t->sy;
    int mh = t->cfg->map_height;
    int mw = t->cfg->map_width;
    int y0 = (int)floor(gy);
    int y1;
    double fy;
    int nr = 0;
    int x = 0;
    int in_run = 0;
    int x0 = 0;

    if (y0 < 0)
        y0 = 0;
    if (y0 > mh - 1)
        y0 = mh - 1;
    y1 = y0 + 1 < mh ? y0 + 1 : y0;
    fy = gy - y0;

    if (t->sx <= 0.0) {
        double th = bg_threshold_y(t, 0, y0, y1, fy);

        while (x < width) {
            while (x < width && row[x] <= th)
                x++;
            if (x >= width)
                break;
            runs[nr].x0 = x;
            while (x + 1 < width && row[x + 1] > th)
                x++;
            runs[nr].x1 = x;
            runs[nr].label = 0;
            nr++;
            x++;
        }
        return nr;
    }

    while (x < width) {
        double gx = x * t->sx;
        int tx = (int)floor(gx);
        int tx1;
        int x_end;
        double a;
        double b;
        double th;
        double dth;

        if (tx < 0)
            tx = 0;
        if (tx > mw - 1)
            tx = mw - 1;
        tx1 = tx + 1 < mw ? tx + 1 : tx;
        x_end = tx + 1 < mw ? (int)ceil((tx + 1) / t->sx) : width;
        if (x_end > width)
            x_end = width;
        a = bg_threshold_y(t, tx, y0, y1, fy);
        b = bg_threshold_y(t, tx1, y0, y1, fy);
        th = a + (b - a) * (gx - tx);
        dth = (b - a) * t->sx;

        while (x < x_end) {
            int above = row[x] > th;

            if (above) {
                if (!in_run) {
                    x0 = x;
                    in_run = 1;
                }
            } else if (in_run) {
                runs[nr].x0 = x0;
                runs[nr].x1 = x - 1;
                runs[nr].label = 0;
                nr++;
                in_run = 0;
            }
            x++;
            th += dth;
        }
    }
    if (in_run) {
        runs[nr].x0 = x0;
        runs[nr].x1 = width - 1;
        runs[nr].label = 0;
        nr++;
    }
    return nr;
}

int ost_bg_extract_fused(const OSTBGConfig *cfg, const uint16_t *image,
                         int stride, OSTBGStats *stats, OSTCCContext *cc,
                         OSTCCComponent *stars, int stars_max)
{
    OSTBGThresholdTest test;

    if (!cfg || !image || !stats || !cc || !stars)
        return -1;
    test.cfg = cfg;
    test.mean = stats->mean;
    test.var = stats->var;
    test.sx = cfg->width > 1 ?
        (double)(cfg->map_width - 1) / (cfg->width - 1) : 0.0;
    test.sy = cfg->height > 1 ?
        (double)(cfg->map_height - 1) / (cfg->height - 1) : 0.0;
    return ost_cc_weighted_runs_4(image, cfg->width, cfg->height, stride,
                                  bg_threshold_runs, &test,
                                  ost_bg_interpolate, stats,
                                  cfg->threshold_sigma, stars, stars_max, cc);
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

static int ost_bg_fit_init(const OSTBGConfig *cfg,
                           const OSTCCComponent *components,
                           int component_count, OSTBGFitStar *stars,
                           double *params, int max_stars)
{
    double eig_sum = 0.0;
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
    for (int i = 0; i < n; i++)
        eig_sum += components[i].eig_min;
    params[3 * n] = sqrt(fmax(eig_sum / (n ? n : 1), 1.0 / 12.0));
    return n;
}

static int chol3(const double *a, double *l)
{
    double d;

    d = a[0];
    if (d <= 0)
        return -1;
    l[0] = sqrt(d);
    l[1] = a[1] / l[0];
    l[3] = a[3] / l[0];
    d = a[2] - l[1] * l[1];
    if (d <= 0)
        return -1;
    l[2] = sqrt(d);
    l[4] = (a[4] - l[3] * l[1]) / l[2];
    d = a[5] - l[3] * l[3] - l[4] * l[4];
    if (d <= 0)
        return -1;
    l[5] = sqrt(d);
    return 0;
}

static void chol3_solve(const double *l, const double *b, double *x)
{
    double y0;
    double y1;
    double y2;

    y0 = b[0] / l[0];
    y1 = (b[1] - l[1] * y0) / l[2];
    y2 = (b[2] - l[3] * y0 - l[4] * y1) / l[5];

    x[2] = y2 / l[5];
    x[1] = (y1 - l[4] * x[2]) / l[2];
    x[0] = (y0 - l[1] * x[1] - l[3] * x[2]) / l[0];
}

static void chol3_inv_diag(const double *l, double *d)
{
    double e[3];
    double x[3];

    e[0] = 1; e[1] = 0; e[2] = 0;
    chol3_solve(l, e, x);
    d[0] = x[0];
    e[0] = 0; e[1] = 1;
    chol3_solve(l, e, x);
    d[1] = x[1];
    e[1] = 0; e[2] = 1;
    chol3_solve(l, e, x);
    d[2] = x[2];
}

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
    J[3] = -sqrt2 * norm *
        (ex * (e3 * y1 - e4 * y2) + ey * (e1 * x1 - e2 * x2));
}

static int build_fit_model(const OSTBGConfig *cfg, const uint16_t *image,
                           int stride, const OSTBGStats *stats,
                           const OSTBGFitStar *stars, const double *params,
                           int n, double sigma,
                           OSTBGFitStar *stars_out, double *params_out,
                           double *normal, double *sigma_col, double *rhs,
                           double *sigma_normal, double *sigma_rhs,
                           double *dropped, int *dropped_count,
                           int max_dropped)
{
    int m = 0;
    int r = cfg->sample_radius;

    *sigma_normal = 0.0;
    *sigma_rhs = 0.0;
    for (int i = 0; i < n; i++) {
        double B[6] = {0, 0, 0, 0, 0, 0};
        double g[3] = {0, 0, 0};
        double b[3] = {0, 0, 0};
        double gs = 0.0;
        double bs = 0.0;
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
                double J[4];

                bg_at(stats, x, y, &mu, &var, &poisson);
                obs = (double)image[(size_t)y * (size_t)stride + (size_t)x] - mu;
                var = fmax((double)image[(size_t)y * (size_t)stride + (size_t)x] *
                           poisson, 0.0) + var;
                psf_eval(x0, y0, I, sigma, x, y, &pred, J);
                if (!(pred < cfg->max_pixel_brightness && pred > sqrt(var)))
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
                    g[0] += J[0] * w * J[3];
                    g[1] += J[1] * w * J[3];
                    g[2] += J[2] * w * J[3];
                    b[0] += J[0] * w * e;
                    b[1] += J[1] * w * e;
                    b[2] += J[2] * w * e;
                    gs += J[3] * w * J[3];
                    bs += J[3] * w * e;
                }
            }
        }
        if (n_valid < 4 || bad_var) {
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
            memcpy(sigma_col + 3 * m, g, sizeof(g));
            memcpy(rhs + 3 * m, b, sizeof(b));
        }
        *sigma_normal += gs;
        *sigma_rhs += bs;
        m++;
    }
    params_out[3 * m] = sigma;
    return m;
}

static int solve_fit(double *params, int n, double *normal, double *sigma_col,
                     double *rhs, double sigma_normal, double sigma_rhs,
                     double *sigma_solve, double *rhs_solve,
                     double *cov_xy, int compute_cov)
{
    double S = sigma_normal;
    double bs = sigma_rhs;
    double cov_s;

    for (int i = 0; i < n; i++) {
        double l[6];

        if (chol3(normal + 6 * i, l) < 0)
            return -1;
        chol3_solve(l, sigma_col + 3 * i, sigma_solve + 3 * i);
        chol3_solve(l, rhs + 3 * i, rhs_solve + 3 * i);
        S -= sigma_col[3 * i + 0] * sigma_solve[3 * i + 0] +
             sigma_col[3 * i + 1] * sigma_solve[3 * i + 1] +
             sigma_col[3 * i + 2] * sigma_solve[3 * i + 2];
        bs -= sigma_col[3 * i + 0] * rhs_solve[3 * i + 0] +
              sigma_col[3 * i + 1] * rhs_solve[3 * i + 1] +
              sigma_col[3 * i + 2] * rhs_solve[3 * i + 2];
        if (compute_cov) {
            double d[3];
            chol3_inv_diag(l, d);
            cov_xy[2 * i + 0] = d[0];
            cov_xy[2 * i + 1] = d[1];
        }
    }
    if (S <= 0.0)
        return -1;
    cov_s = 1.0 / S;
    params[3 * n] += 0.5 * bs * cov_s;
    for (int i = 0; i < n; i++) {
        double ds = bs * cov_s;

        params[3 * i + 0] +=
            0.5 * (rhs_solve[3 * i + 0] - sigma_solve[3 * i + 0] * ds);
        params[3 * i + 1] +=
            0.5 * (rhs_solve[3 * i + 1] - sigma_solve[3 * i + 1] * ds);
        params[3 * i + 2] +=
            0.5 * (rhs_solve[3 * i + 2] - sigma_solve[3 * i + 2] * ds);
        if (compute_cov) {
            cov_xy[2 * i + 0] += sigma_solve[3 * i + 0] *
                                 sigma_solve[3 * i + 0] * cov_s;
            cov_xy[2 * i + 1] += sigma_solve[3 * i + 1] *
                                 sigma_solve[3 * i + 1] * cov_s;
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
        max_stars <= 0)
        return -1;

    stars = work->stars1;
    stars_next = work->stars2;
    params = work->params1;
    params_next = work->params2;
    if (!stars || !stars_next || !params || !params_next || !work->normal ||
        !work->sigma_col || !work->rhs || !work->sigma_solve ||
        !work->rhs_solve || !work->cov_xy || !work->dropped)
        return -1;

    n = ost_bg_fit_init(cfg, components, component_count, stars, params, max_stars);
    if (n <= 0) {
        *dropped_count_out = 0;
        return n;
    }

    for (int it = 0; it < num_iter; it++) {
        double sigma = fmax(params[3 * n], sqrt(1.0 / 12.0));
        double sigma_normal;
        double sigma_rhs;
        int compute_cov = it == num_iter - 1;
        int m;

        m = build_fit_model(cfg, image, stride, stats, stars, params, n, sigma,
                            stars_next, params_next, work->normal,
                            work->sigma_col, work->rhs,
                            &sigma_normal, &sigma_rhs, work->dropped,
                            &dropped_count, max_stars);
        if (m <= 0)
            break;
        if (solve_fit(params_next, m, work->normal, work->sigma_col, work->rhs,
                      sigma_normal, sigma_rhs, work->sigma_solve,
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
        double sigma = fmax(params[3 * n], sqrt(1.0 / 12.0));
        double sigma_normal;
        double sigma_rhs;
        int m;

        m = build_fit_model(cfg, image, stride, stats, stars, params, n, sigma,
                            stars_next, params_next, NULL, NULL, NULL,
                            &sigma_normal, &sigma_rhs, work->dropped,
                            &dropped_count, max_stars);
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

/* tracker */

#define PI 3.14159265358979323846
#define MAX_STARS 1000
#define MAX_CAT 120000
#define KEY_CAP 262144
#define MAX_FILTERED 30000
#define MAX_CDB 600000
#define MAX_NEAR 1024
#define MAX_LOCAL_CDB 8192
#define MAX_CANDIDATES 65536
#define MAX_COLLISION 16384

typedef float Vec3[3];
typedef float Mat3[9];
typedef const float *BEAST_RESTRICT Vec3In;
typedef float *BEAST_RESTRICT Vec3Out;
typedef const float *BEAST_RESTRICT Mat3In;
typedef float *BEAST_RESTRICT Mat3Out;

typedef struct BeastConfig {
    int IMG_X, IMG_Y, MAX_FALSE_STARS, DB_REDUNDANCY, REQUIRED_STARS;
    int KDBUCKET_SIZE;
    float PIXSCALE, DOUBLE_STAR_PX, BASE_FLUX, IMAGE_VARIANCE;
    float THRESH_FACTOR, POS_VARIANCE, POS_ERR_SIGMA;
    float MAXFOV, MINFOV, MATCH_VALUE, PIXX_TANGENT, PIXY_TANGENT;
} BeastConfig;
typedef BeastConfig Config;

typedef struct BeastStar {
    union {
        Vec3 v;
        struct { float x, y, z; };
    };
    float flux;
    float px, py, sigma_sq;
    int id, star_idx, unreliable;
} BeastStar;
typedef BeastStar Star;

typedef struct BeastStarDB {
    Star *v;
    int n, cap;
    float max_variance;
} BeastStarDB;
typedef BeastStarDB StarDB;

typedef struct BeastStarQuery {
    Star *map;
    int n, kdsorted;
    int *kdresults;
    int kdresults_size, kdresults_maxsize;
    signed char *kdmask;
} BeastStarQuery;
typedef BeastStarQuery Query;

typedef struct BeastConstellation {
    float p;
    int s1, s2, idx;
} BeastConstellation;
typedef BeastConstellation Constellation;

typedef struct BeastConstellationPair {
    float totalscore;
    int db_s1, db_s2, img_s1, img_s2;
} BeastConstellationPair;
typedef BeastConstellationPair CPair;

typedef struct BeastConstellationDB {
    StarDB stars;
    Query results;
    Constellation *map;
    int map_size;
} BeastConstellationDB;
typedef BeastConstellationDB CDB;

typedef struct BeastStarFov {
    int *mask, *collision;
    int collision_size, collision_cap;
    float *s_px, *s_py, maxdist_sq, sigma_sq;
} BeastStarFov;
typedef BeastStarFov StarFov;

typedef struct BeastMatchResult {
    CPair match;
    Vec3 R[3];
    int *map;
    int map_size;
    CDB *db, *img;
    StarFov *img_mask;
} BeastMatchResult;
typedef BeastMatchResult MatchResult;

typedef struct BeastMatchWork {
    CPair *candidates;
    int candidate_cap;
    int *fov_mask, *collision;
    int collision_cap;
    float *fov_px, *fov_py, *scores;
    int *match_map, *work_map;
} BeastMatchWork;

typedef struct Work {
    Config cfg;
    Star cat[MAX_CAT];
    Star filtered[MAX_FILTERED];
    Star near_stars[MAX_NEAR];
    Star img_stars[MAX_STARS], img_bright[MAX_STARS];
    Star tmp_stars[MAX_NEAR > MAX_STARS ? MAX_NEAR : MAX_STARS];
    Star q_full_map[MAX_CAT], q_filtered_map[MAX_FILTERED];
    Star q_near_map[MAX_NEAR], q_img_map[MAX_STARS], q_img2_map[MAX_STARS];
    int q_full_results[MAX_CAT + 1], q_filtered_results[MAX_FILTERED + 1];
    int q_near_results[MAX_NEAR + 1], q_img_results[MAX_STARS + 1];
    int q_img2_results[MAX_STARS + 1];
    signed char q_full_mask[MAX_CAT], q_filtered_mask[MAX_FILTERED];
    signed char q_near_mask[MAX_NEAR], q_img_mask[MAX_STARS], q_img2_mask[MAX_STARS];
    signed char keep_full[MAX_CAT], keep_filtered[MAX_FILTERED];
    Constellation cdb_map[MAX_CDB], local_cdb_map[MAX_LOCAL_CDB];
    Constellation img_cmap[16], img2_cmap[MAX_STARS * (MAX_STARS - 1) / 2];
    Constellation rel_cmap[MAX_STARS * (MAX_STARS - 1) / 2];
    CPair candidates[MAX_CANDIDATES];
    int *fov_mask;
    int collision[MAX_COLLISION];
    float fov_px[MAX_STARS], fov_py[MAX_STARS], scores[MAX_STARS];
    int match_map[MAX_STARS], work_map[MAX_STARS], results[MAX_STARS];
    uint64_t cat_keys[KEY_CAP];
    StarDB catalog_db, filtered_db;
    Query full_q, filtered_q;
    CDB global;
    int prepared;
} Work;

BEAST_DEF int beast_load_config(Config *c, const char *filename)
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

BEAST_DEF Star beast_make_db_star(const Config *c, float x, float y, float z, float flux, int id)
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

BEAST_DEF Star beast_make_img_star(const Config *c, float px, float py, float flux, int id)
{
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

BEAST_DEF void beast_star_db_init(StarDB *db, Star *storage, int cap)
{
    db->v = storage;
    db->n = 0;
    db->cap = cap;
    db->max_variance = 0.0f;
}

BEAST_DEF int beast_db_add(StarDB *db, Star s)
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

static int key_seen(uint64_t *tab, uint64_t h)
{
    uint32_t i = (uint32_t)(h ^ (h >> 32)) & (KEY_CAP - 1);
    while (tab[i]) {
        if (tab[i] == h)
            return 1;
        i = (i + 1) & (KEY_CAP - 1);
    }
    tab[i] = h;
    return 0;
}

BEAST_DEF int beast_load_catalog(const Config *cfg, StarDB *db,
                                  const char *filename, float year,
                                  uint64_t *cat_keys)
{
    FILE *f;
    char line[2048], *field[78];
    float yd = year - 1991.25f;

    f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        return -1;
    }
    memset(cat_keys, 0, KEY_CAP * sizeof(cat_keys[0]));
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
            Star s = beast_make_db_star(cfg,
                                  cos(PI * ra / 180.0) * cosdec,
                                  sin(PI * ra / 180.0) * cosdec,
                                  sin(PI * dec / 180.0),
                                  cfg->BASE_FLUX * powf(10.0f, -mag / 2.5f),
                                  atoi(field[1]));
            if (!key_seen(cat_keys, v3_key(s.v)) && beast_db_add(db, s) < 0) {
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

static void sort_flux_desc(Star *a, int n)
{
    for (int i = 1; i < n; i++) {
        Star s = a[i];
        int j = i;
        while (j > 0 && s.flux > a[j - 1].flux) {
            a[j] = a[j - 1];
            j--;
        }
        a[j] = s;
    }
}

static inline void kdselect(Star *a, int l, int r, int k, int dim)
{
    r--;
    while (l < r) {
        float p = a[(l + r) >> 1].v[dim];
        int i = l, j = r;
        while (i <= j) {
            while (a[i].v[dim] < p) i++;
            while (a[j].v[dim] > p) j--;
            if (i <= j) {
                BEAST_SWAP(Star, a[i], a[j]);
                i++;
                j--;
            }
        }
        if (k <= j) r = j;
        else if (k >= i) l = i;
        else break;
    }
}

static void kdbuild(Star *a, int min, int max, int bucket, int dim)
{
    int mid = (min + max) / 2;
    if (min + 1 < max) {
        int next = dim < 2 ? dim + 1 : 0;
        kdselect(a, min, max, mid, dim);
        if (mid - min > bucket) kdbuild(a, min, mid, bucket, next);
        else sort_flux_desc(a + min, mid - min);
        if (max - (mid + 1) > bucket) kdbuild(a, mid + 1, max, bucket, next);
        else sort_flux_desc(a + mid + 1, max - (mid + 1));
    }
}

BEAST_DEF void beast_query_init(Query *q, StarDB *db, Star *map, int *res, signed char *mask)
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

BEAST_DEF void beast_query_reset_mask(Query *q)
{
    memset(q->kdmask, 0, (size_t)(unsigned)q->n);
}

BEAST_DEF void beast_query_clear_results(Query *q)
{
    while (q->kdresults_size > 0)
        q->kdmask[q->kdresults[--q->kdresults_size]] = 0;
}

BEAST_DEF void beast_query_sort_flux(Query *q)
{
    qsort(q->map, (size_t)q->n, sizeof(q->map[0]), cmp_flux_desc);
    q->kdsorted = 0;
}

BEAST_DEF void beast_query_kdsort(Query *q, const Config *c)
{
    if (!q->kdsorted) {
        kdbuild(q->map, 0, q->n, c->KDBUCKET_SIZE, 0);
        q->kdsorted = 1;
    }
}

static inline void kdcheck(Query *q, int idx, const float p[3], float r, float min_flux)
{
    Star *s = &q->map[idx];
    float dx = p[0] - s->v[0], dy = p[1] - s->v[1], dz = p[2] - s->v[2];
    if (dx - r <= 0 && 0 <= dx + r &&
        dy - r <= 0 && 0 <= dy + r &&
        dz - r <= 0 && 0 <= dz + r &&
        min_flux <= s->flux && q->kdmask[idx] == 0 &&
        dx * dx + dy * dy + dz * dz <= r * r) {
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
    }
}

static void kdsearch(Query *q, const float p[3], float r, float min_flux,
                     int min, int max, int bucket, int dim)
{
    int mid = (min + max) / 2;
    int next = dim < 2 ? dim + 1 : 0;
    float t = p[dim];
    if (min < mid && t - r <= q->map[mid].v[dim]) {
        if (mid - min > bucket) {
            kdsearch(q, p, r, min_flux, min, mid, bucket, next);
        } else {
            for (int i = min; i < mid && min_flux <= q->map[i].flux; i++) {
                kdcheck(q, i, p, r, min_flux);
                if (q->kdresults_size == q->kdresults_maxsize)
                    min_flux = q->map[q->kdresults[q->kdresults_size - 1]].flux;
            }
        }
    }
    if (mid < max)
        kdcheck(q, mid, p, r, min_flux);
    if (q->kdresults_size == q->kdresults_maxsize)
        min_flux = q->map[q->kdresults[q->kdresults_size - 1]].flux;
    if (mid + 1 < max && q->map[mid].v[dim] <= t + r) {
        if (max - (mid + 1) > bucket) {
            kdsearch(q, p, r, min_flux, mid + 1, max, bucket, next);
        } else {
            if (q->kdresults_size == q->kdresults_maxsize)
                min_flux = q->map[q->kdresults[q->kdresults_size - 1]].flux;
            for (int i = mid + 1; i < max && min_flux <= q->map[i].flux; i++) {
                kdcheck(q, i, p, r, min_flux);
                if (q->kdresults_size == q->kdresults_maxsize)
                    min_flux = q->map[q->kdresults[q->kdresults_size - 1]].flux;
            }
        }
    }
}

BEAST_DEF void beast_query_search(Query *q, const Config *c, const float p[3],
                         float arcsec, float min_flux)
{
    float r = arcsec / 3600.0f * (float)PI / 180.0f;
    beast_query_kdsort(q, c);
    kdsearch(q, p, 2 * fabsf(sin(r / 2.0f)), min_flux, 0, q->n,
             c->KDBUCKET_SIZE, 0);
}

BEAST_DEF void beast_query_search_range(Query *q, const Config *c, const float p[3],
                                  float arcsec, float min_flux, int min,
                                  int max, int dim)
{
    float r = arcsec / 3600.0f * (float)PI / 180.0f;
    beast_query_kdsort(q, c);
    kdsearch(q, p, 2 * fabsf(sinf(r / 2.0f)), min_flux, min, max,
             c->KDBUCKET_SIZE, dim);
}

BEAST_DEF void beast_query_mask_filter(Query *q, const Config *c)
{
    beast_query_kdsort(q, c);
    for (int i = 0; i < q->n; i++) {
        int lastmask = q->kdmask[i];
        beast_query_search(q, c, q->map[i].v,
                     c->DOUBLE_STAR_PX * c->PIXSCALE,
                     c->THRESH_FACTOR * c->IMAGE_VARIANCE);
        if (q->kdresults_size > 1 || lastmask ||
            q->map[i].flux < c->THRESH_FACTOR * c->IMAGE_VARIANCE) {
            q->kdmask[i] = 1;
            q->kdresults_size = 0;
        } else {
            beast_query_clear_results(q);
        }
    }
}

BEAST_DEF void beast_query_mask_uniform(Query *q, const Config *c, int min_stars, signed char *keep)
{
    int oldmax = q->kdresults_maxsize;
    memset(keep, 0, (size_t)q->n);
    q->kdresults_maxsize = min_stars;
    for (int i = 0; i < q->n; i++) if (!q->kdmask[i]) {
        beast_query_search(q, c, q->map[i].v,
                     c->MINFOV / 2, c->THRESH_FACTOR * c->IMAGE_VARIANCE);
        for (int j = 0; j < q->kdresults_size; j++)
            keep[q->kdresults[j]] = 1;
        beast_query_clear_results(q);
    }
    for (int i = 0; i < q->n; i++)
        q->kdmask[i] = keep[i] ? 0 : 1;
    q->kdresults_maxsize = oldmax;
}

BEAST_DEF int beast_db_from_mask(StarDB *out, Query *q)
{
    out->n = 0;
    for (int i = 0; i < q->n; i++)
        if (!q->kdmask[i] && beast_db_add(out, q->map[i]) < 0)
            return -1;
    return 0;
}

BEAST_DEF int beast_db_from_results(StarDB *out, Query *q)
{
    out->n = 0;
    for (int i = 0; i < q->kdresults_size; i++)
        if (beast_db_add(out, q->map[q->kdresults[i]]) < 0)
            return -1;
    return 0;
}

static float star_dist_arcsec(const Star *a, const Star *b)
{
    float x = a->v[0] * b->v[1] - b->v[0] * a->v[1];
    float y = a->v[0] * b->v[2] - b->v[0] * a->v[2];
    float z = a->v[1] * b->v[2] - b->v[1] * a->v[2];
    return (float)((3600 * 180.0 / PI) * asinf(sqrtf(x * x + y * y + z * z)));
}

static int cmp_constellation(const void *pa, const void *pb)
{
    const Constellation *a = (const Constellation *)pa, *b = (const Constellation *)pb;
    if (a->p < b->p) return -1;
    if (a->p > b->p) return 1;
    if (a->s1 != b->s1) return a->s1 - b->s1;
    return a->s2 - b->s2;
}

BEAST_DEF int beast_db_from_image(CDB *cdb, StarDB *src,
                                  Star *star_storage, int star_cap,
                                  Query *q, Star *query_map,
                                  int *query_results, signed char *query_mask,
                                  Constellation *cmap, int cmap_cap,
                                  int stars_per_fov)
{
    int ns, idx = 0;
    beast_star_db_init(&cdb->stars, star_storage, star_cap);
    if (db_copy(&cdb->stars, src) < 0)
        return -1;
    beast_query_init(q, &cdb->stars, query_map, query_results, query_mask);
    beast_query_sort_flux(q);
    cdb->results = *q;
    cdb->map = cmap;
    ns = cdb->stars.n;
    if (ns > stars_per_fov)
        ns = stars_per_fov;
    for (int j = 1; j < ns; j++) {
        for (int i = 0; i < j; i++) {
            if (idx >= cmap_cap)
                return -1;
            cdb->map[idx].p = star_dist_arcsec(&q->map[i], &q->map[j]);
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

BEAST_DEF int beast_db_from_catalog(CDB *cdb, StarDB *src,
                               Star *star_storage, int star_cap,
                               Query *q, Star *query_map,
                               int *query_results, signed char *query_mask,
                               Constellation *cmap, int cmap_cap,
                               int stars_per_fov, const Config *cfg,
                               signed char *keep)
{
    int n = 0, out = 0;
    beast_star_db_init(&cdb->stars, star_storage, star_cap);
    if (db_copy(&cdb->stars, src) < 0)
        return -1;
    beast_query_init(q, &cdb->stars, query_map, query_results, query_mask);
    beast_query_mask_uniform(q, cfg, stars_per_fov, keep);
    for (int i = 0; i < q->n; i++) if (!q->kdmask[i]) {
        beast_query_search(q, cfg, q->map[i].v,
                     cfg->MAXFOV, cfg->THRESH_FACTOR * cfg->IMAGE_VARIANCE);
        for (int j = 0; j < q->kdresults_size; j++) {
            int k = q->kdresults[j];
            if (i != k && q->map[i].flux >= q->map[k].flux) {
                if (n >= cmap_cap)
                    return -1;
                cmap[n].p = star_dist_arcsec(&q->map[i], &q->map[k]);
                cmap[n].s1 = q->map[i].star_idx;
                cmap[n].s2 = q->map[k].star_idx;
                n++;
            }
        }
        beast_query_clear_results(q);
    }
    beast_query_reset_mask(q);
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

static void constellation_range(Constellation *a, int n, float p0, float p1, int *lo, int *hi)
{
    int l = 0, r = n;
    while (l < r) {
        int m = (l + r) >> 1;
        if (a[m].p < p0) l = m + 1;
        else r = m;
    }
    *lo = l;
    r = n;
    while (l < r) {
        int m = (l + r) >> 1;
        if (a[m].p <= p1) l = m + 1;
        else r = m;
    }
    *hi = l;
}

static int fov_init(StarFov *fov, StarDB *stars, float db_max_variance,
                    const Config *c, BeastMatchWork *w)
{
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

static inline float fov_score(StarFov *fov, int id, float px, float py)
{
    float dx = px - fov->s_px[id], dy = py - fov->s_py[id];
    if (dx < -0.5f) dx += 1.0f;
    if (dy < -0.5f) dy += 1.0f;
    return (fov->maxdist_sq - (dx * dx + dy * dy)) /
           (2.0f * fov->sigma_sq);
}

static int fov_resolve(StarFov *fov, int id, float px, float py)
{
    int id1, id2;
    if (id >= -1)
        return id;
    id = -id;
    id1 = fov_resolve(fov, fov->collision[id - 2], px, py);
    id2 = fov_resolve(fov, fov->collision[id - 1], px, py);
    return fov_score(fov, id1, px, py) > fov_score(fov, id2, px, py) ? id1 : id2;
}

static int fov_get_id(StarFov *fov, const Config *c, float px, float py)
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

static void mr_init(MatchResult *m, CDB *db, CDB *img, StarFov *mask, int *map)
{
    memset(m, 0, sizeof(*m));
    m->db = db; m->img = img; m->img_mask = mask; m->map = map;
    m->map_size = img->stars.n;
    m->match.totalscore = -FLT_MAX;
}

static void mr_set_pair(MatchResult *m, Constellation db, Constellation img)
{
    m->match.img_s1 = img.s1;
    m->match.img_s2 = img.s2;
    m->match.db_s1 = db.s1;
    m->match.db_s2 = db.s2;
}

static void mr_copy(MatchResult *dst, MatchResult *src)
{
    CDB *db = dst->db, *img = dst->img;
    StarFov *mask = dst->img_mask;
    int *map = dst->map, map_size = dst->map_size;
    *dst = *src;
    dst->db = db; dst->img = img; dst->img_mask = mask;
    dst->map = map; dst->map_size = map_size;
    memcpy(dst->map, src->map, (size_t)map_size * sizeof(map[0]));
}

#define M3(a, r, c) ((a)[(r) * 3 + (c)])

static inline float v3_dot(Vec3In a, Vec3In b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline void v3_cross(Vec3Out r, Vec3In a, Vec3In b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void v3_normalize(Vec3Out v)
{
    float n = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] /= n; v[1] /= n; v[2] /= n;
}

static inline void v3_neg(Vec3Out v)
{
    v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2];
}

static inline void m3_basis(Mat3Out m, Vec3In a, Vec3In b, Vec3In c)
{
    M3(m, 0, 0) = a[0]; M3(m, 0, 1) = b[0]; M3(m, 0, 2) = c[0];
    M3(m, 1, 0) = a[1]; M3(m, 1, 1) = b[1]; M3(m, 1, 2) = c[1];
    M3(m, 2, 0) = a[2]; M3(m, 2, 1) = b[2]; M3(m, 2, 2) = c[2];
}

static inline void m3_mul_bt(Mat3Out c, Mat3In a, Mat3In b)
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M3(c, i, j) = M3(a, i, 0) * M3(b, j, 0) +
                          M3(a, i, 1) * M3(b, j, 1) +
                          M3(a, i, 2) * M3(b, j, 2);
        }
    }
}

static void weighted_triad(MatchResult *m)
{
    Star *db_s1 = &m->db->stars.v[m->match.db_s1], *db_s2 = &m->db->stars.v[m->match.db_s2];
    Star *img_s1 = &m->img->stars.v[m->match.img_s1], *img_s2 = &m->img->stars.v[m->match.img_s2];
    const float *wa = db_s1->v, *wb = db_s2->v, *va = img_s1->v, *vb = img_s2->v;
    float weightA = 1.0f / (db_s1->sigma_sq + img_s1->sigma_sq);
    float weightB = 1.0f / (db_s2->sigma_sq + img_s2->sigma_sq);
    float sumAB = weightA + weightB;
    Vec3 wc, vc, waXwc, vaXvc, wbXwc, vbXvc;
    Mat3 wbase, vbase, A, B;
    float cz, sz, mz, cy, sy, my, cx, sx, mx;

    v3_cross(wc, wa, wb);
    v3_normalize(wc);
    v3_cross(vc, va, vb);
    v3_normalize(vc);
    v3_cross(waXwc, wa, wc);
    v3_cross(vaXvc, va, vc);
    m3_basis(wbase, wa, waXwc, wc);
    m3_basis(vbase, va, vaXvc, vc);
    m3_mul_bt(A, vbase, wbase);

    v3_neg(wc);
    v3_neg(vc);
    v3_cross(wbXwc, wb, wc);
    v3_cross(vbXvc, vb, vc);
    m3_basis(wbase, wb, wbXwc, wc);
    m3_basis(vbase, vb, vbXvc, vc);
    m3_mul_bt(B, vbase, wbase);

    weightA /= sumAB; weightB /= sumAB;
    cz = weightA * M3(A, 0, 0) + weightB * M3(B, 0, 0);
    sz = weightA * M3(A, 1, 0) + weightB * M3(B, 1, 0);
    mz = sqrtf(cz * cz + sz * sz); cz /= mz; sz /= mz;
    cy = weightA * sqrtf(M3(A, 2, 1) * M3(A, 2, 1) + M3(A, 2, 2) * M3(A, 2, 2)) +
         weightB * sqrtf(M3(B, 2, 1) * M3(B, 2, 1) + M3(B, 2, 2) * M3(B, 2, 2));
    sy = -weightA * M3(A, 2, 0) - weightB * M3(B, 2, 0);
    my = sqrtf(cy * cy + sy * sy); cy /= my; sy /= my;
    cx = weightA * M3(A, 2, 2) + weightB * M3(B, 2, 2);
    sx = weightA * M3(A, 2, 1) + weightB * M3(B, 2, 1);
    mx = sqrtf(cx * cx + sx * sx); cx /= mx; sx /= mx;
    m->R[0][0] = cy * cz;
    m->R[0][1] = cz * sx * sy - cx * sz;
    m->R[0][2] = sx * sz + cx * cz * sy;
    m->R[1][0] = cy * sz;
    m->R[1][1] = cx * cz + sx * sy * sz;
    m->R[1][2] = cx * sy * sz - cz * sx;
    m->R[2][0] = -sy;
    m->R[2][1] = cy * sx;
    m->R[2][2] = cx * cy;
}

#undef M3

static void compute_score(MatchResult *m, const Config *c, BeastMatchWork *w)
{
    m->match.totalscore = log(1.0 / (c->IMG_X * c->IMG_Y)) * (2 * m->map_size);
    for (int i = 0; i < m->map_size; i++) {
        m->map[i] = -1;
        w->scores[i] = 0.0f;
    }
    for (int i = 0; i < m->db->results.kdresults_size; i++) {
        Star *s = &m->db->results.map[m->db->results.kdresults[i]];
        int o = s->star_idx;
        float x = v3_dot(s->v, m->R[0]);
        float y = v3_dot(s->v, m->R[1]);
        float z = v3_dot(s->v, m->R[2]);
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

static int related(MatchResult *winner, CPair *p)
{
    if (winner->match.totalscore == -FLT_MAX || p->totalscore == -FLT_MAX)
        return 0;
    return winner->map[p->img_s1] == p->db_s1 && winner->map[p->img_s2] == p->db_s2;
}

BEAST_DEF int beast_db_match(CDB *db, CDB *img, MatchResult *winner,
                    const Config *cfg, BeastMatchWork *w, float *p_match)
{
    StarFov fov;
    MatchResult m;
    CPair *candidates = w->candidates;
    int candidate_cap = w->candidate_cap;
    int nc = 0;

    *p_match = 0.0f;
    mr_init(winner, db, img, &fov, w->match_map);
    if (db->stars.n < 3 || img->stars.n < 3)
        return 0;
    if (fov_init(&fov, &img->stars, db->stars.max_variance, cfg, w) < 0)
        return -1;
    mr_init(&m, db, img, &fov, w->work_map);
    for (int n = 0; n < img->map_size; n++) {
        Constellation ic = img->map[n];
        float err = cfg->POS_ERR_SIGMA * cfg->PIXSCALE *
            sqrtf(img->stars.v[ic.s1].sigma_sq + img->stars.v[ic.s2].sigma_sq +
                  2 * db->stars.max_variance);
        int lo, hi;
        constellation_range(db->map, db->map_size, ic.p - err, ic.p + err, &lo, &hi);
        if (lo >= hi)
            continue;
        for (int o = lo; o < hi; o++) {
            mr_set_pair(&m, db->map[o], ic);
            weighted_triad(&m);
            if (db->results.kdsorted)
                beast_query_search(&db->results, cfg, m.R[0], cfg->MAXFOV / 2,
                             cfg->THRESH_FACTOR * cfg->IMAGE_VARIANCE);
            for (int flip = 0; flip < 2; flip++) {
                compute_score(&m, cfg, w);
                if (m.match.totalscore > winner->match.totalscore) {
                    if (winner->match.totalscore != -FLT_MAX) {
                        if (nc >= candidate_cap)
                            return -1;
                        candidates[nc++] = winner->match;
                    }
                    mr_copy(winner, &m);
                } else {
                    if (nc >= candidate_cap)
                        return -1;
                    candidates[nc++] = m.match;
                }
                BEAST_SWAP(int, m.match.img_s1, m.match.img_s2);
                if (flip == 0)
                    weighted_triad(&m);
            }
            if (db->results.kdsorted)
                beast_query_clear_results(&db->results);
        }
    }
    if (winner->match.totalscore != -FLT_MAX) {
        double p = 1.0;
        for (int i = 0; i < nc; i++)
            if (!related(winner, &candidates[i]))
                p += exp((double)candidates[i].totalscore - winner->match.totalscore);
        *p_match = (float)(1.0 / p);
    }
    return 0;
}

BEAST_DEF void beast_match_work_init(BeastMatchWork *mw,
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

BEAST_DEF int beast_copy_n_brightest(StarDB *dst, StarDB *src, Star *tmp, int n)
{
    memcpy(tmp, src->v, (size_t)src->n * sizeof(tmp[0]));
    qsort(tmp, (size_t)src->n, sizeof(tmp[0]), cmp_flux_desc);
    dst->n = 0;
    if (n > src->n)
        n = src->n;
    for (int i = 0; i < n; i++)
        if (beast_db_add(dst, tmp[i]) < 0)
            return -1;
    return 0;
}

#endif

#endif
