#include "ost_test_tracker.h"

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
