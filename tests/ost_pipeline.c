#define BEAST_IMPLEMENTATION
#include "../beast/beast.h"

static void work_match_init(Work *w, BeastMatchWork *mw)
{
    beast_match_work_init(mw, w->candidates, MAX_CANDIDATES,
                           w->fov_mask, w->collision, MAX_COLLISION,
                           w->fov_px, w->fov_py, w->scores,
                           w->match_map, w->work_map);
}

static int star_measurements_to_img_db(Work *w, StarDB *db, Star *storage,
                                       const double *stars, int len,
                                       int ids_from_index);

static int match_catalog_stars(Work *w, Query *full_q, CDB *global,
                               const double *stars, int *result, int len)
{
    StarDB img, bright, near;
    CDB img_cdb, fov_cdb, img_full_cdb;
    Query q_img = {0}, q_near = {0}, q_img2 = {0};
    MatchResult winner, fov_winner;
    BeastMatchWork mw;
    float p_match;

    work_match_init(w, &mw);
    if (star_measurements_to_img_db(w, &img, w->img_stars, stars, len, 0) < 0)
        return -1;
    beast_star_db_init(&bright, w->img_bright, MAX_STARS);
    beast_star_db_init(&near, w->near_stars, MAX_NEAR);
    for (int i = 0; i < len; i++)
        result[i] = -1;
    if (beast_copy_n_brightest(&bright, &img, w->tmp_stars,
                         w->cfg.MAX_FALSE_STARS + w->cfg.REQUIRED_STARS) < 0)
        return -1;
    if (beast_db_from_image(&img_cdb, &bright, w->tmp_stars, MAX_STARS,
                               &q_img, w->q_img_map, w->q_img_results,
                               w->q_img_mask, w->img_cmap, 16,
                               w->cfg.MAX_FALSE_STARS + 2) < 0)
        return -1;
    if (beast_db_match(global, &img_cdb, &winner, &w->cfg, &mw, &p_match) < 0)
        return -1;
    if (p_match > 0.9f) {
        beast_query_search(full_q, &w->cfg, winner.R[0],
                     w->cfg.MAXFOV / 2, w->cfg.THRESH_FACTOR * w->cfg.IMAGE_VARIANCE);
        beast_query_search(&global->results, &w->cfg, winner.R[0],
                     w->cfg.MAXFOV / 2, w->cfg.THRESH_FACTOR * w->cfg.IMAGE_VARIANCE);
        if (beast_db_from_results(&near, full_q) < 0)
            return -1;
        if (beast_db_from_image(&fov_cdb, &near, w->near_stars, MAX_NEAR,
                                   &q_near, w->q_near_map, w->q_near_results,
                                   w->q_near_mask, w->local_cdb_map,
                                   MAX_LOCAL_CDB,
                                   global->results.kdresults_size) < 0)
            return -1;
        beast_query_clear_results(&global->results);
        beast_query_clear_results(full_q);

        if (beast_db_from_image(&img_full_cdb, &img, w->img_bright,
                                   MAX_STARS, &q_img2, w->q_img2_map,
                                   w->q_img2_results, w->q_img2_mask,
                                   w->img2_cmap,
                                   MAX_STARS * (MAX_STARS - 1) / 2,
                                   w->cfg.MAX_FALSE_STARS + 2) < 0)
            return -1;
        if (beast_db_match(&fov_cdb, &img_full_cdb, &fov_winner,
                     &w->cfg, &mw, &p_match) < 0)
            return -1;
        for (int i = 0; i < len; i++) {
            int dbi = fov_winner.map[i];
            result[i] = (dbi >= 0) ? fov_cdb.stars.v[dbi].id : -1;
        }
    }
    return 0;
}

static int star_measurements_to_img_db(Work *w, StarDB *db, Star *storage,
                                       const double *stars, int len,
                                       int ids_from_index)
{
    beast_star_db_init(db, storage, MAX_STARS);
    for (int i = 0; i < len; i++) {
        Star s = beast_make_img_star(&w->cfg,
                               (float)(stars[3 * i] - w->cfg.IMG_X / 2.0),
                               (float)(-(stars[3 * i + 1] - w->cfg.IMG_Y / 2.0)),
                               w->cfg.BASE_FLUX * powf(10.0f, (float)(-stars[3 * i + 2] / 2.5)),
                               ids_from_index ? i : -1);
        if (beast_db_add(db, s) < 0)
            return -1;
    }
    return 0;
}

static int prepare_catalog(Work *tw, int *fov_mask,
                           const char *catalog_path, float year)
{
    Work *w = (Work *)tw;
    const char *path = catalog_path ? catalog_path : "hip_main.dat";

    if (!w || !fov_mask || w->cfg.IMG_X <= 0 || w->cfg.IMG_Y <= 0)
        return -1;
    w->fov_mask = fov_mask;
    beast_star_db_init(&w->catalog_db, w->cat, MAX_CAT);
    beast_star_db_init(&w->filtered_db, w->filtered, MAX_FILTERED);
    if (beast_load_catalog(&w->cfg, &w->catalog_db, path, year,
                             w->cat_keys) < 0)
        return -1;
    beast_query_init(&w->full_q, &w->catalog_db,
               w->q_full_map, w->q_full_results, w->q_full_mask);
    beast_query_mask_filter(&w->full_q, &w->cfg);
    beast_query_mask_uniform(&w->full_q, &w->cfg, w->cfg.REQUIRED_STARS,
                       w->keep_full);
    if (beast_db_from_mask(&w->filtered_db, &w->full_q) < 0)
        return -1;
    beast_query_reset_mask(&w->full_q);
    if (beast_db_from_catalog(&w->global, &w->filtered_db, w->filtered,
                           MAX_FILTERED, &w->filtered_q,
                           w->q_filtered_map, w->q_filtered_results,
                           w->q_filtered_mask, w->cdb_map, MAX_CDB,
                           2 + w->cfg.DB_REDUNDANCY, &w->cfg,
                           w->keep_filtered) < 0)
        return -1;
    w->prepared = 1;
    return 0;
}


#include <errno.h>
#include <math.h>
#include <png.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Pipeline {
    Work *tracker;
    OSTBGConfig bg_cfg;
    OSTBGStats bg_stats;
    OSTCCContext cc;
    unsigned char *rgba;
    uint16_t *gray;
    double *bg_mean, *bg_var, *bg_poisson;
    int *x_edge, *y_edge, *hist;
    OSTCCComponent *cc_work, *components;
    int *parent, *label_live, *free_after_row, *touched_stamp, *seen_stamp;
    OSTCCRun *prev_runs, *curr_runs;
    OSTBGFitWorkspace fit_work;
    double *fit_params, *fit_cov, *stars;
    int *result;
} Pipeline;

static int read_png_rgba(const char *filename, int w, int h,
                         unsigned char *rgba)
{
    png_image image;

    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;
    if (!png_image_begin_read_from_file(&image, filename)) {
        fprintf(stderr, "%s: %s\n", filename, image.message);
        return -1;
    }
    if ((int)image.width != w || (int)image.height != h) {
        fprintf(stderr, "%s: got %ux%u, expected %dx%d\n",
                filename, image.width, image.height, w, h);
        png_image_free(&image);
        return -1;
    }
    image.format = PNG_FORMAT_RGBA;
    if (!png_image_finish_read(&image, NULL, rgba, 0, NULL)) {
        fprintf(stderr, "%s: %s\n", filename, image.message);
        png_image_free(&image);
        return -1;
    }
    png_image_free(&image);
    return 0;
}

static void write_stars(FILE *f, const double *stars, int n)
{
    for (int i = 0; i < 3 * n; i++)
        fprintf(f, "%.17g%c", stars[i], i == 3 * n - 1 ? '\n' : ',');
    if (n == 0)
        fputc('\n', f);
}

static void write_ids(const int *ids, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d%c", ids[i], i == n - 1 ? '\n' : ',');
    if (n == 0)
        putchar('\n');
}

static int process_image(Pipeline *p, const char *filename, FILE *stars_file)
{
    int n, fit_n, dropped, out_n = 0;
    float base_flux;

    if (read_png_rgba(filename, p->bg_cfg.width, p->bg_cfg.height,
                      p->rgba) < 0)
        return -1;
    ost_bg_rgba_to_gray16(p->gray, p->rgba,
                          p->bg_cfg.width, p->bg_cfg.height);
    if (ost_bg_compute_stats(&p->bg_cfg, p->gray, p->bg_cfg.width,
                             p->bg_mean, p->bg_var, p->bg_poisson,
                             p->x_edge, p->y_edge, p->hist, 65536) < 0)
        return -1;
    n = ost_bg_extract_fused(&p->bg_cfg, p->gray, p->bg_cfg.width,
                             &p->bg_stats, &p->cc, p->components,
                             p->bg_cfg.max_stars);
    if (n < 0)
        return -1;
    fit_n = ost_bg_fit_stars(&p->bg_cfg, p->gray, p->bg_cfg.width,
                             &p->bg_stats, p->components, n, 10,
                             &p->fit_work, p->bg_cfg.max_stars,
                             p->fit_params, p->fit_cov, &dropped);
    if (fit_n < 0)
        return -1;

    base_flux = p->tracker->cfg.BASE_FLUX;
    if (base_flux <= 0.0f)
        return -1;
    for (int i = 0; i < fit_n; i++) {
        double x = p->fit_params[3 * i + 0];
        double y = p->fit_params[3 * i + 1];
        double flux = p->fit_params[3 * i + 2];

        if (!isfinite(x) || !isfinite(y) || !isfinite(flux) || flux <= 0.0)
            continue;
        p->stars[3 * out_n + 0] = x;
        p->stars[3 * out_n + 1] = y;
        p->stars[3 * out_n + 2] = -2.5 * log10(flux / base_flux);
        out_n++;
    }
    if (stars_file)
        write_stars(stars_file, p->stars, out_n);
    if (match_catalog_stars(p->tracker, &p->tracker->full_q,
                            &p->tracker->global, p->stars,
                            p->result, out_n) < 0)
        return -1;
    write_ids(p->result, out_n);
    fprintf(stderr, "%s: components=%d fitted=%d used=%d dropped=%d sigma=%.9g\n",
            filename, n, fit_n, out_n, dropped,
            fit_n ? p->fit_params[3 * fit_n] : 0.0);
    return 0;
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "usage: %s [--catalog hip_main.dat] [--stars-out csv] "
            "calibration.txt year image.png [...]\n", prog);
}

int main(int argc, char **argv)
{
    Pipeline p;
    OSTCCBufferSizes cc_sizes;
    const char *catalog_path = "hip_main.dat";
    const char *stars_path = NULL;
    FILE *stars_file = NULL;
    int arg = 1, width, height, max_stars;
    int *fov_mask = NULL;
    size_t pixels, map_pixels;
    int rc = 1;

    memset(&p, 0, sizeof(p));
    while (arg < argc && !strncmp(argv[arg], "--", 2)) {
        if (!strcmp(argv[arg], "--catalog") && arg + 1 < argc) {
            catalog_path = argv[arg + 1];
            arg += 2;
        } else if (!strcmp(argv[arg], "--stars-out") && arg + 1 < argc) {
            stars_path = argv[arg + 1];
            arg += 2;
        } else {
            usage(argv[0]);
            return 1;
        }
    }
    if (argc - arg < 3) {
        usage(argv[0]);
        return 1;
    }

    p.tracker = (Work *)malloc(sizeof(*p.tracker));
    if (!p.tracker) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    memset(p.tracker, 0, sizeof(*p.tracker));
    if (beast_load_config(&p.tracker->cfg, argv[arg]) < 0)
        goto done;
    width = p.tracker->cfg.IMG_X;
    height = p.tracker->cfg.IMG_Y;
    if (ost_bg_config_init(&p.bg_cfg, width, height) < 0 ||
        ost_cc_buffer_sizes(width, &cc_sizes) < 0)
        goto done;
    max_stars = p.bg_cfg.max_stars;
    pixels = (size_t)width * (size_t)height;
    map_pixels = (size_t)p.bg_cfg.map_width * (size_t)p.bg_cfg.map_height;

    fov_mask = (int *)malloc(pixels * sizeof(*fov_mask));
    p.rgba = (unsigned char *)malloc(4 * pixels);
    p.gray = (uint16_t *)malloc(pixels * sizeof(*p.gray));
    p.bg_mean = (double *)malloc(map_pixels * sizeof(*p.bg_mean));
    p.bg_var = (double *)malloc(map_pixels * sizeof(*p.bg_var));
    p.bg_poisson = (double *)malloc(map_pixels * sizeof(*p.bg_poisson));
    p.x_edge = (int *)malloc((size_t)(p.bg_cfg.map_width + 1) *
                             sizeof(*p.x_edge));
    p.y_edge = (int *)malloc((size_t)(p.bg_cfg.map_height + 1) *
                             sizeof(*p.y_edge));
    p.hist = (int *)malloc(65536 * sizeof(*p.hist));
    p.cc_work = (OSTCCComponent *)malloc(cc_sizes.components *
                                         sizeof(*p.cc_work));
    p.parent = (int *)malloc(cc_sizes.parent * sizeof(*p.parent));
    p.label_live = (int *)malloc(cc_sizes.label_live * sizeof(*p.label_live));
    p.free_after_row = (int *)malloc(cc_sizes.free_after_row *
                                     sizeof(*p.free_after_row));
    p.touched_stamp = (int *)malloc(cc_sizes.touched_stamp *
                                    sizeof(*p.touched_stamp));
    p.seen_stamp = (int *)malloc(cc_sizes.seen_stamp * sizeof(*p.seen_stamp));
    p.prev_runs = (OSTCCRun *)malloc(cc_sizes.prev_runs *
                                     sizeof(*p.prev_runs));
    p.curr_runs = (OSTCCRun *)malloc(cc_sizes.curr_runs *
                                     sizeof(*p.curr_runs));
    p.components = (OSTCCComponent *)malloc((size_t)max_stars *
                                           sizeof(*p.components));
    p.fit_work.stars1 = (OSTBGFitStar *)malloc((size_t)max_stars *
                                               sizeof(*p.fit_work.stars1));
    p.fit_work.stars2 = (OSTBGFitStar *)malloc((size_t)max_stars *
                                               sizeof(*p.fit_work.stars2));
    p.fit_work.params1 = (double *)malloc((size_t)(3 * max_stars + 1) *
                                          sizeof(*p.fit_work.params1));
    p.fit_work.params2 = (double *)malloc((size_t)(3 * max_stars + 1) *
                                          sizeof(*p.fit_work.params2));
    p.fit_work.normal = (double *)malloc((size_t)(6 * max_stars) *
                                         sizeof(*p.fit_work.normal));
    p.fit_work.sigma_col = (double *)malloc((size_t)(3 * max_stars) *
                                            sizeof(*p.fit_work.sigma_col));
    p.fit_work.rhs = (double *)malloc((size_t)(3 * max_stars) *
                                      sizeof(*p.fit_work.rhs));
    p.fit_work.sigma_solve = (double *)malloc((size_t)(3 * max_stars) *
                                              sizeof(*p.fit_work.sigma_solve));
    p.fit_work.rhs_solve = (double *)malloc((size_t)(3 * max_stars) *
                                            sizeof(*p.fit_work.rhs_solve));
    p.fit_work.cov_xy = (double *)malloc((size_t)(2 * max_stars) *
                                         sizeof(*p.fit_work.cov_xy));
    p.fit_work.dropped = (double *)malloc((size_t)(3 * max_stars) *
                                          sizeof(*p.fit_work.dropped));
    p.fit_params = (double *)malloc((size_t)(3 * max_stars + 1) *
                                    sizeof(*p.fit_params));
    p.fit_cov = (double *)malloc((size_t)(2 * max_stars) *
                                 sizeof(*p.fit_cov));
    p.stars = (double *)malloc((size_t)(3 * max_stars) *
                               sizeof(*p.stars));
    p.result = (int *)malloc((size_t)max_stars * sizeof(*p.result));

    if (!fov_mask || !p.rgba || !p.gray || !p.bg_mean ||
        !p.bg_var || !p.bg_poisson || !p.x_edge || !p.y_edge || !p.hist ||
        !p.cc_work || !p.parent || !p.label_live || !p.free_after_row ||
        !p.touched_stamp || !p.seen_stamp || !p.prev_runs || !p.curr_runs ||
        !p.components || !p.fit_work.stars1 || !p.fit_work.stars2 ||
        !p.fit_work.params1 || !p.fit_work.params2 || !p.fit_work.normal ||
        !p.fit_work.sigma_col || !p.fit_work.rhs ||
        !p.fit_work.sigma_solve || !p.fit_work.rhs_solve ||
        !p.fit_work.cov_xy || !p.fit_work.dropped || !p.fit_params ||
        !p.fit_cov || !p.stars || !p.result) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    p.bg_stats.cfg = &p.bg_cfg;
    p.bg_stats.mean = p.bg_mean;
    p.bg_stats.var = p.bg_var;
    p.bg_stats.poisson = p.bg_poisson;
    p.bg_stats.x_edge = p.x_edge;
    p.bg_stats.y_edge = p.y_edge;
    if (ost_cc_init(&p.cc, width, p.cc_work, p.parent, p.label_live,
                    p.free_after_row, p.touched_stamp, p.seen_stamp,
                    p.prev_runs, p.curr_runs) < 0)
        goto done;
    if (prepare_catalog(p.tracker, fov_mask, catalog_path,
                        (float)atof(argv[arg + 1])) < 0)
        goto done;
    if (stars_path) {
        stars_file = fopen(stars_path, "w");
        if (!stars_file) {
            fprintf(stderr, "%s: %s\n", stars_path, strerror(errno));
            goto done;
        }
    }
    for (int i = arg + 2; i < argc; i++)
        if (process_image(&p, argv[i], stars_file) < 0)
            goto done;
    rc = 0;

done:
    if (stars_file)
        fclose(stars_file);
    free(p.result);
    free(p.stars);
    free(p.fit_cov);
    free(p.fit_params);
    free(p.fit_work.dropped);
    free(p.fit_work.cov_xy);
    free(p.fit_work.rhs_solve);
    free(p.fit_work.sigma_solve);
    free(p.fit_work.rhs);
    free(p.fit_work.sigma_col);
    free(p.fit_work.normal);
    free(p.fit_work.params2);
    free(p.fit_work.params1);
    free(p.fit_work.stars2);
    free(p.fit_work.stars1);
    free(p.components);
    free(p.curr_runs);
    free(p.prev_runs);
    free(p.seen_stamp);
    free(p.touched_stamp);
    free(p.free_after_row);
    free(p.label_live);
    free(p.parent);
    free(p.cc_work);
    free(p.hist);
    free(p.y_edge);
    free(p.x_edge);
    free(p.bg_poisson);
    free(p.bg_var);
    free(p.bg_mean);
    free(p.gray);
    free(p.rgba);
    free(fov_mask);
    free(p.tracker);
    return rc;
}
