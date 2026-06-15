#define OST_IMPLEMENTATION
#include "ost_test_work.h"

static void work_match_init(Work *w, MatchWork *mw)
{
    ost_match_work_init(mw, w->candidates, MAX_CANDIDATES,
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
    Query q_img, q_near, q_img2;
    MatchResult winner, fov_winner;
    MatchWork mw;
    float p_match;

    memset(&q_img, 0, sizeof(q_img));
    memset(&q_near, 0, sizeof(q_near));
    memset(&q_img2, 0, sizeof(q_img2));
    work_match_init(w, &mw);
    if (star_measurements_to_img_db(w, &img, w->img_stars, stars, len, 0) < 0)
        return -1;
    ost_star_db_init(&bright, w->img_bright, MAX_STARS);
    ost_star_db_init(&near, w->near_stars, MAX_NEAR);
    for (int i = 0; i < len; i++)
        result[i] = -1;
    if (ost_copy_n_brightest(&bright, &img, w->tmp_stars,
                         w->cfg.MAX_FALSE_STARS + w->cfg.REQUIRED_STARS) < 0)
        return -1;
    if (ost_db_from_image(&img_cdb, &bright, w->tmp_stars, MAX_STARS,
                               &q_img, w->q_img_map, w->q_img_results,
                               w->q_img_mask, w->img_cmap, 16,
                               w->cfg.MAX_FALSE_STARS + 2) < 0)
        return -1;
    if (ost_db_match(global, &img_cdb, &winner, &w->cfg, &mw, &p_match) < 0)
        return -1;
    if (p_match > 0.9f) {
        ost_query_search(full_q, &w->cfg, winner.R[0],
                     w->cfg.MAXFOV / 2, w->cfg.THRESH_FACTOR * w->cfg.IMAGE_VARIANCE);
        ost_query_search(&global->results, &w->cfg, winner.R[0],
                     w->cfg.MAXFOV / 2, w->cfg.THRESH_FACTOR * w->cfg.IMAGE_VARIANCE);
        if (ost_db_from_results(&near, full_q) < 0)
            return -1;
        if (ost_db_from_image(&fov_cdb, &near, w->near_stars, MAX_NEAR,
                                   &q_near, w->q_near_map, w->q_near_results,
                                   w->q_near_mask, w->local_cdb_map,
                                   MAX_LOCAL_CDB,
                                   global->results.kdresults_size) < 0)
            return -1;
        ost_query_clear_results(&global->results);
        ost_query_clear_results(full_q);

        if (ost_db_from_image(&img_full_cdb, &img, w->img_bright,
                                   MAX_STARS, &q_img2, w->q_img2_map,
                                   w->q_img2_results, w->q_img2_mask,
                                   w->img2_cmap,
                                   MAX_STARS * (MAX_STARS - 1) / 2,
                                   w->cfg.MAX_FALSE_STARS + 2) < 0)
            return -1;
        if (ost_db_match(&fov_cdb, &img_full_cdb, &fov_winner,
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
    ost_star_db_init(db, storage, MAX_STARS);
    for (int i = 0; i < len; i++) {
        Star s = ost_make_img_star(&w->cfg,
                               (float)(stars[3 * i] - w->cfg.IMG_X / 2.0),
                               (float)(-(stars[3 * i + 1] - w->cfg.IMG_Y / 2.0)),
                               w->cfg.BASE_FLUX * powf(10.0f, (float)(-stars[3 * i + 2] / 2.5)),
                               ids_from_index ? i : -1);
        if (ost_db_add(db, s) < 0)
            return -1;
    }
    return 0;
}

static int match_relative_stars(Work *w,
                                const double *reference_stars,
                                int reference_len,
                                const double *current_stars,
                                int *reference_index_result,
                                int current_len,
                                float *p_match)
{
    StarDB reference, current;
    CDB reference_cdb, current_cdb;
    Query q_reference, q_current;
    MatchResult winner;
    MatchWork mw;

    memset(&q_reference, 0, sizeof(q_reference));
    memset(&q_current, 0, sizeof(q_current));
    work_match_init(w, &mw);
    *p_match = 0.0f;
    for (int i = 0; i < current_len; i++)
        reference_index_result[i] = -1;
    if (star_measurements_to_img_db(w, &reference, w->near_stars,
                                    reference_stars, reference_len, 1) < 0 ||
        star_measurements_to_img_db(w, &current, w->img_stars,
                                    current_stars, current_len, 0) < 0)
        return -1;

    if (ost_db_from_image(&reference_cdb, &reference, w->near_stars,
                               MAX_STARS, &q_reference, w->q_near_map,
                               w->q_near_results, w->q_near_mask,
                               w->rel_cmap,
                               MAX_STARS * (MAX_STARS - 1) / 2,
                               w->cfg.MAX_FALSE_STARS + 2) < 0)
        return -1;

    if (ost_db_from_image(&current_cdb, &current, w->img_bright,
                               MAX_STARS, &q_current, w->q_img2_map,
                               w->q_img2_results, w->q_img2_mask,
                               w->img2_cmap,
                               MAX_STARS * (MAX_STARS - 1) / 2,
                               w->cfg.MAX_FALSE_STARS + 2) < 0)
        return -1;

    if (ost_db_match(&reference_cdb, &current_cdb, &winner,
                 &w->cfg, &mw, p_match) < 0)
        return -1;
    if (*p_match > 0.0f) {
        for (int i = 0; i < current_len; i++) {
            int dbi = winner.map[i];
            reference_index_result[i] = (dbi >= 0) ? reference_cdb.stars.v[dbi].id : -1;
        }
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
    ost_star_db_init(&w->catalog_db, w->cat, MAX_CAT);
    ost_star_db_init(&w->filtered_db, w->filtered, MAX_FILTERED);
    if (ost_load_catalog(&w->cfg, &w->catalog_db, path, year,
                             w->cat_keys, KEY_CAP) < 0)
        return -1;
    ost_query_init(&w->full_q, &w->catalog_db,
               w->q_full_map, w->q_full_results, w->q_full_mask);
    ost_query_mask_filter(&w->full_q, &w->cfg);
    ost_query_mask_uniform(&w->full_q, &w->cfg, w->cfg.REQUIRED_STARS,
                       w->keep_full);
    if (ost_db_from_mask(&w->filtered_db, &w->full_q) < 0)
        return -1;
    ost_query_reset_mask(&w->full_q);
    if (ost_db_from_catalog(&w->global, &w->filtered_db, w->filtered,
                           MAX_FILTERED, &w->filtered_q,
                           w->q_filtered_map, w->q_filtered_results,
                           w->q_filtered_mask, w->cdb_map, MAX_CDB,
                           2 + w->cfg.DB_REDUNDANCY, &w->cfg,
                           w->keep_filtered) < 0)
        return -1;
    w->prepared = 1;
    return 0;
}


int main(int argc, char **argv)
{
    Work *tracker;
    FILE *file;
    char line[131072];
    double data[3 * MAX_STARS];
    int result[MAX_STARS];
    int *fov_mask;
    clock_t time_sum = 0;
    int run_times = 0;
    int rc = 1;
    int relative_self = 0;
    float p_match;
    int arg = 1;

    if (argc > 1 && !strcmp(argv[1], "--relative-self")) {
        relative_self = 1;
        arg = 2;
    }
    if (argc - arg < 3) {
        printf("./test_ost_c [--relative-self] input.csv calibration.txt year\n");
        return 0;
    }

    tracker = (Work *)calloc(1, sizeof(*tracker));
    if (!tracker) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    if (ost_load_config(&tracker->cfg, argv[arg + 1]) < 0)
        goto done_tracker;
    fov_mask = (int *)malloc((size_t)tracker->cfg.IMG_X *
                             (size_t)tracker->cfg.IMG_Y *
                             sizeof(*fov_mask));
    if (!fov_mask) {
        fprintf(stderr, "out of memory\n");
        goto done_tracker;
    }
    if (prepare_catalog(tracker, fov_mask, "hip_main.dat",
                        (float)atof(argv[arg + 2])) < 0)
        goto done_fov;

    file = fopen(argv[arg], "r");
    if (!file) {
        fprintf(stderr, "%s: %s\n", argv[arg], strerror(errno));
        goto done_fov;
    }
    if (relative_self) {
        int i = 0, len;
        char *tok;

        if (!fgets(line, sizeof(line), file)) {
            fprintf(stderr, "empty input\n");
            fclose(file);
            goto done_fov;
        }
        tok = strtok(line, ",");
        while (tok && i < 3 * MAX_STARS) {
            data[i++] = atof(tok);
            tok = strtok(NULL, ",");
        }
        len = i / 3;
        if (match_relative_stars(tracker, data, len, data, result, len,
                                 &p_match) < 0) {
            fprintf(stderr, "relative star match failed\n");
            fclose(file);
            goto done_fov;
        }
        if (p_match <= 0.99f) {
            fprintf(stderr, "relative star self-test confidence too low: %g\n", p_match);
            fclose(file);
            goto done_fov;
        }
        for (i = 0; i < len; i++) {
            if (result[i] != i) {
                fprintf(stderr, "relative star self-test failed at %d: got %d\n",
                        i, result[i]);
                fclose(file);
                goto done_fov;
            }
        }
        fclose(file);
        rc = 0;
        goto done_fov;
    }

    while (fgets(line, sizeof(line), file)) {
        int i = 0;
        int len;
        char *tok = strtok(line, ",");
        clock_t start;
        clock_t end;

        while (tok && i < 3 * MAX_STARS) {
            data[i++] = atof(tok);
            tok = strtok(NULL, ",");
        }
        len = i / 3;
        start = clock();
        if (match_catalog_stars(tracker, &tracker->full_q, &tracker->global,
                                data, result, len) < 0) {
            fprintf(stderr, "star match failed\n");
            fclose(file);
            goto done_fov;
        }
        end = clock();
        time_sum += end - start;
        for (i = 0; i < len; i++)
            printf("%d%c", result[i], i == len - 1 ? '\n' : ',');
        run_times++;
        fprintf(stderr, "Time: %f\n",
                (float)time_sum / (CLOCKS_PER_SEC * run_times));
    }
    fclose(file);
    rc = 0;

done_fov:
    free(fov_mask);
done_tracker:
    free(tracker);
    return rc;
}
