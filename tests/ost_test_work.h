#ifndef OST_TEST_WORK_H
#define OST_TEST_WORK_H

#include <stdint.h>
#include "../ost/ost.h"

enum {
    MAX_STARS = 1000,
    MAX_CAT = 120000,
    KEY_CAP = 262144,
    MAX_FILTERED = 30000,
    MAX_CDB = 600000,
    MAX_NEAR = 1024,
    MAX_LOCAL_CDB = 8192,
    MAX_CANDIDATES = 65536,
    MAX_COLLISION = 16384
};

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

#endif
