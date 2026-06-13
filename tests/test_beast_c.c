#define C_BEAST_IMPLEMENTATION
#include "../c_beast/c_beast.h"

int main(int argc, char **argv)
{
    OSTTracker *tracker;
    FILE *file;
    char line[131072];
    double data[3 * OST_TRACKER_MAX_STARS];
    int result[OST_TRACKER_MAX_STARS];
    int *fov_mask;
    clock_t time_sum = 0;
    int run_times = 0;
    int rc = 1;

    if (argc < 4) {
        printf("./test_beast_c input.csv calibration.txt year\n");
        return 0;
    }

    tracker = (OSTTracker *)calloc(1, ost_tracker_work_size());
    if (!tracker) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    if (ost_tracker_configure(tracker, argv[2]) < 0)
        goto done_tracker;
    fov_mask = (int *)malloc((size_t)ost_tracker_width(tracker) *
                             (size_t)ost_tracker_height(tracker) *
                             sizeof(*fov_mask));
    if (!fov_mask) {
        fprintf(stderr, "out of memory\n");
        goto done_tracker;
    }
    if (ost_tracker_prepare(tracker, fov_mask, "hip_main.dat",
                            (float)atof(argv[3])) < 0)
        goto done_fov;

    file = fopen(argv[1], "r");
    if (!file) {
        fprintf(stderr, "%s: %s\n", argv[1], strerror(errno));
        goto done_fov;
    }
    while (fgets(line, sizeof(line), file)) {
        int i = 0;
        int len;
        char *tok = strtok(line, ",");
        clock_t start;
        clock_t end;

        while (tok && i < 3 * OST_TRACKER_MAX_STARS) {
            data[i++] = atof(tok);
            tok = strtok(NULL, ",");
        }
        len = i / 3;
        start = clock();
        if (ost_tracker_solve_spikes(tracker, data, result, len) < 0) {
            fprintf(stderr, "star_id failed\n");
            fclose(file);
            goto done_fov;
        }
        end = clock();
        time_sum += end - start;
        for (i = 0; i < len; i++)
            printf("%d%c", result[i], i == len - 1 ? '\n' : ',');
        run_times++;
        fprintf(stderr, "Time on edison: %f\n",
                ((float)time_sum * C_BEAST_EDISON_SPEED_FACTOR) /
                (CLOCKS_PER_SEC * run_times));
    }
    fclose(file);
    rc = 0;

done_fov:
    free(fov_mask);
done_tracker:
    free(tracker);
    return rc;
}
