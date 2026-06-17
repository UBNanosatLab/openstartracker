#include "../ost/ost.h"

#include <errno.h>
#include <math.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define IMG_MAX_STARS 2048

typedef struct {
    int width, height;
    double image_variance, thresh_factor;
} ImgConfig;

static int load_img_config(ImgConfig *cfg, const char *filename)
{
    FILE *f;
    char line[256], key[128], val[128];

    memset(cfg, 0, sizeof(*cfg));
    f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        return -1;
    }

    while (fgets(line, sizeof(line), f)) {
        if (sscanf(line, " %127[^=]=%127s", key, val) != 2)
            continue;
        if (!strcmp(key, "IMG_X"))
            cfg->width = atoi(val);
        else if (!strcmp(key, "IMG_Y"))
            cfg->height = atoi(val);
        else if (!strcmp(key, "IMAGE_VARIANCE"))
            cfg->image_variance = atof(val);
        else if (!strcmp(key, "THRESH_FACTOR"))
            cfg->thresh_factor = atof(val);
    }
    fclose(f);

    if (cfg->width <= 0 || cfg->height <= 0 ||
        cfg->image_variance <= 0 || cfg->thresh_factor <= 0) {
        fprintf(stderr, "%s: missing required calibration keys\n", filename);
        return -1;
    }
    return 0;
}

static int read_png_rgba(const char *filename, int w, int h, unsigned char *rgba)
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

static void make_threshold_image(unsigned char *dst,
                                 const unsigned char *img,
                                 const unsigned char *med,
                                 int w, int h)
{
    int n = w * h;

    for (int i = 0; i < n; i++) {
        int r = img[4 * i + 0] - med[4 * i + 0];
        int g = img[4 * i + 1] - med[4 * i + 1];
        int b = img[4 * i + 2] - med[4 * i + 2];
        if (r < 0) r = 0;
        if (g < 0) g = 0;
        if (b < 0) b = 0;
        dst[i] = (unsigned char)((77 * r + 150 * g + 29 * b) >> 8);
    }
}

static unsigned char clamp_threshold(double x)
{
    if (x < 0)
        return 0;
    if (x > 255)
        return 255;
    return (unsigned char)floor(x + 0.5);
}

static int cmp_component_flux(const void *pa, const void *pb)
{
    const OSTCCComponent *a = (const OSTCCComponent *)pa;
    const OSTCCComponent *b = (const OSTCCComponent *)pb;

    if (a->area != b->area)
        return b->area - a->area;
    if (a->sum_y != b->sum_y)
        return a->sum_y - b->sum_y;
    return a->sum_x - b->sum_x;
}

static double now_sec(void)
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

int main(int argc, char **argv)
{
    ImgConfig cfg;
    OSTBGConfig bg_cfg;
    OSTBGStats bg_stats;
    OSTCCBufferSizes sizes;
    OSTCCContext cc;
    unsigned char *median_rgba, *image_rgba, *gray;
    uint16_t *gray16;
    double *bg_mean, *bg_var, *bg_poisson;
    OSTCCComponent *work_comp, *stars;
    OSTBGFitWorkspace fit_work;
    OSTBGFitStar *fit_stars1, *fit_stars2;
    double *fit_params1, *fit_params2, *fit_params, *fit_cov;
    double *fit_normal, *fit_sigma_col, *fit_rhs;
    double *fit_sigma_solve, *fit_rhs_solve, *fit_cov_xy, *fit_dropped;
    int *parent, *col_label, *active_count, *free_after_row;
    int *x_edge, *y_edge, *hist;
    unsigned char threshold;
    size_t pixels;
    int bg_mode;
    int fit_mode;
    int image_mode;
    int calib_arg;
    int first_image_arg;
    int rc = 1;
    int profile = 0;
    int profile_images = 0;
    double profile_start = 0.0;
    double t_read = 0.0, t_gray = 0.0, t_stats = 0.0;
    double t_extract = 0.0, t_fit = 0.0;

    bg_mode = argc > 1 && !strcmp(argv[1], "--bg");
    fit_mode = argc > 1 && !strcmp(argv[1], "--fit");
    image_mode = bg_mode || fit_mode;
    calib_arg = image_mode ? 2 : 1;
    first_image_arg = image_mode ? 3 : 3;
    if ((!image_mode && argc < 4) || (image_mode && argc < 4)) {
        fprintf(stderr, "usage: %s calibration.txt median.png image.png [...]\n", argv[0]);
        fprintf(stderr, "       %s --bg calibration.txt image.png [...]\n", argv[0]);
        fprintf(stderr, "       %s --fit calibration.txt image.png [...]\n", argv[0]);
        return 1;
    }
    if (load_img_config(&cfg, argv[calib_arg]) < 0)
        return 1;
    if (ost_bg_config_init(&bg_cfg, cfg.width, cfg.height) < 0)
        return 1;
    if (ost_cc_buffer_sizes(cfg.width, &sizes) < 0)
        return 1;

    pixels = (size_t)cfg.width * (size_t)cfg.height;
    median_rgba = (unsigned char *)malloc(pixels * 4);
    image_rgba = (unsigned char *)malloc(pixels * 4);
    gray = (unsigned char *)malloc(pixels);
    gray16 = (uint16_t *)malloc(pixels * sizeof(*gray16));
    bg_mean = (double *)malloc((size_t)bg_cfg.map_width *
                               (size_t)bg_cfg.map_height * sizeof(*bg_mean));
    bg_var = (double *)malloc((size_t)bg_cfg.map_width *
                              (size_t)bg_cfg.map_height * sizeof(*bg_var));
    bg_poisson = (double *)malloc((size_t)bg_cfg.map_width *
                                  (size_t)bg_cfg.map_height * sizeof(*bg_poisson));
    x_edge = (int *)malloc((size_t)(bg_cfg.map_width + 1) * sizeof(*x_edge));
    y_edge = (int *)malloc((size_t)(bg_cfg.map_height + 1) * sizeof(*y_edge));
    hist = (int *)malloc(65536 * sizeof(*hist));
    work_comp = (OSTCCComponent *)malloc(sizes.components * sizeof(*work_comp));
    parent = (int *)malloc(sizes.parent * sizeof(*parent));
    col_label = (int *)malloc(sizes.col_label * sizeof(*col_label));
    active_count = (int *)malloc(sizes.active_count * sizeof(*active_count));
    free_after_row = (int *)malloc(sizes.free_after_row * sizeof(*free_after_row));
    stars = (OSTCCComponent *)malloc(IMG_MAX_STARS * sizeof(*stars));
    fit_stars1 = (OSTBGFitStar *)malloc(bg_cfg.max_stars * sizeof(*fit_stars1));
    fit_stars2 = (OSTBGFitStar *)malloc(bg_cfg.max_stars * sizeof(*fit_stars2));
    fit_params1 = (double *)malloc((3 * bg_cfg.max_stars + 1) * sizeof(*fit_params1));
    fit_params2 = (double *)malloc((3 * bg_cfg.max_stars + 1) * sizeof(*fit_params2));
    fit_params = (double *)malloc((3 * bg_cfg.max_stars + 1) * sizeof(*fit_params));
    fit_cov = (double *)malloc(2 * bg_cfg.max_stars * sizeof(*fit_cov));
    fit_normal = (double *)malloc(6 * bg_cfg.max_stars * sizeof(*fit_normal));
    fit_sigma_col = (double *)malloc(3 * bg_cfg.max_stars * sizeof(*fit_sigma_col));
    fit_rhs = (double *)malloc(3 * bg_cfg.max_stars * sizeof(*fit_rhs));
    fit_sigma_solve = (double *)malloc(3 * bg_cfg.max_stars * sizeof(*fit_sigma_solve));
    fit_rhs_solve = (double *)malloc(3 * bg_cfg.max_stars * sizeof(*fit_rhs_solve));
    fit_cov_xy = (double *)malloc(2 * bg_cfg.max_stars * sizeof(*fit_cov_xy));
    fit_dropped = (double *)malloc(3 * bg_cfg.max_stars * sizeof(*fit_dropped));

    if (!median_rgba || !image_rgba || !gray || !gray16 ||
        !bg_mean || !bg_var || !bg_poisson || !x_edge || !y_edge || !hist ||
        !work_comp ||
        !parent || !col_label || !active_count || !free_after_row || !stars ||
        !fit_stars1 || !fit_stars2 || !fit_params1 || !fit_params2 ||
        !fit_params || !fit_cov || !fit_normal || !fit_sigma_col ||
        !fit_rhs || !fit_sigma_solve || !fit_rhs_solve || !fit_cov_xy ||
        !fit_dropped) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    fit_work.stars1 = fit_stars1;
    fit_work.stars2 = fit_stars2;
    fit_work.params1 = fit_params1;
    fit_work.params2 = fit_params2;
    fit_work.normal = fit_normal;
    fit_work.sigma_col = fit_sigma_col;
    fit_work.rhs = fit_rhs;
    fit_work.sigma_solve = fit_sigma_solve;
    fit_work.rhs_solve = fit_rhs_solve;
    fit_work.cov_xy = fit_cov_xy;
    fit_work.dropped = fit_dropped;
    if (ost_cc_init(&cc, cfg.width, work_comp, parent, col_label,
                    active_count, free_after_row) < 0)
        goto done;
    if (!image_mode && read_png_rgba(argv[2], cfg.width, cfg.height, median_rgba) < 0)
        goto done;

    profile = getenv("OST_PROFILE") != NULL;
    profile_start = now_sec();
    threshold = clamp_threshold(cfg.thresh_factor * cfg.image_variance);
    for (int arg = first_image_arg; arg < argc; arg++) {
        int n;
        double t0;
        double t1;

        t0 = now_sec();
        if (read_png_rgba(argv[arg], cfg.width, cfg.height, image_rgba) < 0)
            goto done;
        t1 = now_sec();
        t_read += t1 - t0;
        if (image_mode) {
            int kept = 0;
            int fit_n;
            int dropped_count;
            double sigma = 0;

            t0 = now_sec();
            ost_bg_rgba_to_gray16(gray16, image_rgba, cfg.width, cfg.height);
            t1 = now_sec();
            t_gray += t1 - t0;
            t0 = now_sec();
            if (ost_bg_compute_stats(&bg_cfg, gray16, cfg.width,
                                     bg_mean, bg_var, bg_poisson,
                                     x_edge, y_edge, hist, 65536) < 0)
                goto done;
            t1 = now_sec();
            t_stats += t1 - t0;
            bg_stats.cfg = &bg_cfg;
            bg_stats.mean = bg_mean;
            bg_stats.var = bg_var;
            bg_stats.poisson = bg_poisson;
            bg_stats.x_edge = x_edge;
            bg_stats.y_edge = y_edge;
            t0 = now_sec();
            n = ost_bg_extract_fused(&bg_cfg, gray16, cfg.width,
                                     &bg_stats, &cc, stars,
                                     bg_cfg.max_stars);
            t1 = now_sec();
            t_extract += t1 - t0;
            if (n < 0) {
                fprintf(stderr, "%s: background extraction failed (%d)\n", argv[arg], n);
                goto done;
            }
            if (fit_mode) {
                t0 = now_sec();
                fit_n = ost_bg_fit_stars(&bg_cfg, gray16, cfg.width, &bg_stats,
                                         stars, n, 3, &fit_work,
                                         bg_cfg.max_stars, fit_params,
                                         fit_cov, &dropped_count);
                t1 = now_sec();
                t_fit += t1 - t0;
                if (fit_n < 0) {
                    fprintf(stderr, "%s: star fit failed (%d)\n", argv[arg], fit_n);
                    goto done;
                }
                profile_images++;
                printf("%s,%d,%d,%.12g", argv[arg], fit_n, dropped_count,
                       fit_n ? fit_params[3 * fit_n] : 0.0);
                for (int i = 0; i < fit_n; i++) {
                    printf(",%.9f,%.9f,%.9f,%.9f,%.9f",
                           fit_params[3 * i + 0],
                           fit_params[3 * i + 1],
                           fit_params[3 * i + 2],
                           sqrt(fmax(fit_cov[2 * i + 0], 0.0)),
                           sqrt(fmax(fit_cov[2 * i + 1], 0.0)));
                }
                putchar('\n');
                continue;
            }
            profile_images++;
            for (int i = 0; i < n; i++) {
                double x = stars[i].wx / stars[i].wsum;
                double y = stars[i].wy / stars[i].wsum;
                int xi = (int)lrint(x);
                int yi = (int)lrint(y);
                if (xi < bg_cfg.sample_radius ||
                    xi >= cfg.width - bg_cfg.sample_radius ||
                    yi < bg_cfg.sample_radius ||
                    yi >= cfg.height - bg_cfg.sample_radius)
                    continue;
                kept++;
            }
            for (int i = 0; i < kept; i++)
                sigma += stars[i].eig_min;
            sigma = sqrt(fmax(sigma / (kept ? kept : 1), 1.0 / 12.0));

            printf("%s,%d,%.12g", argv[arg], kept, sigma);
            for (int i = 0; i < n; i++) {
                double x = stars[i].wx / stars[i].wsum;
                double y = stars[i].wy / stars[i].wsum;
                int xi = (int)lrint(x);
                int yi = (int)lrint(y);
                if (xi < bg_cfg.sample_radius ||
                    xi >= cfg.width - bg_cfg.sample_radius ||
                    yi < bg_cfg.sample_radius ||
                    yi >= cfg.height - bg_cfg.sample_radius)
                    continue;
                printf(",%.9f,%.9f,%.9f,%.9f,%d",
                       x, y, stars[i].wsum, stars[i].eig_min, stars[i].area);
            }
            putchar('\n');
            continue;
        }
        make_threshold_image(gray, image_rgba, median_rgba, cfg.width, cfg.height);
        n = ost_cc_threshold_4(gray, cfg.width, cfg.height, cfg.width,
                               threshold, stars, IMG_MAX_STARS, &cc);
        if (n < 0) {
            fprintf(stderr, "%s: connected components failed (%d)\n", argv[arg], n);
            goto done;
        }
        qsort(stars, (size_t)n, sizeof(*stars), cmp_component_flux);

        printf("%s,%d", argv[arg], n);
        for (int i = 0; i < n; i++) {
            double cx = stars[i].sum_x / (double)stars[i].area;
            double cy = stars[i].sum_y / (double)stars[i].area;
            printf(",%.3f,%.3f,%d", cx, cy, stars[i].area);
        }
        putchar('\n');
    }
    rc = 0;

done:
    if (profile) {
        double total = now_sec() - profile_start;
        double measured = t_read + t_gray + t_stats + t_extract + t_fit;

        fprintf(stderr,
                "profile images=%d total=%.6f read=%.6f gray=%.6f "
                "stats=%.6f fused_extract=%.6f fit=%.6f other=%.6f\n",
                profile_images, total, t_read, t_gray, t_stats,
                t_extract, t_fit, total - measured);
    }
    free(fit_dropped);
    free(fit_cov_xy);
    free(fit_rhs_solve);
    free(fit_sigma_solve);
    free(fit_rhs);
    free(fit_sigma_col);
    free(fit_normal);
    free(fit_cov);
    free(fit_params);
    free(fit_params2);
    free(fit_params1);
    free(fit_stars2);
    free(fit_stars1);
    free(stars);
    free(free_after_row);
    free(active_count);
    free(col_label);
    free(parent);
    free(work_comp);
    free(hist);
    free(y_edge);
    free(x_edge);
    free(bg_poisson);
    free(bg_var);
    free(bg_mean);
    free(gray16);
    free(gray);
    free(image_rgba);
    free(median_rgba);
    return rc;
}
