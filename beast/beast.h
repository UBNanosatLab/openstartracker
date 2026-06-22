#ifndef OST_COMPAT_BEAST_H
#define OST_COMPAT_BEAST_H

#include "../ost/ost.h"

#if defined(__GNUC__) || defined(__clang__)
#warning "beast/beast.h is deprecated; include ost/ost.h and use ost_* C symbols instead."
#endif

#ifdef __cplusplus
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace beast_compat_detail {
static inline Config& cfg() {
    static Config c;
    return c;
}

static inline float flux_from_mag(float mag) {
    return cfg().BASE_FLUX * std::pow(10.0f, -mag / 2.5f);
}

}

#define IMG_X (beast_compat_detail::cfg().IMG_X)
#define IMG_Y (beast_compat_detail::cfg().IMG_Y)
#define MAX_FALSE_STARS (beast_compat_detail::cfg().MAX_FALSE_STARS)
#define DB_REDUNDANCY (beast_compat_detail::cfg().DB_REDUNDANCY)
#define REQUIRED_STARS (beast_compat_detail::cfg().REQUIRED_STARS)
#define PIXSCALE (beast_compat_detail::cfg().PIXSCALE)
#define DOUBLE_STAR_PX (beast_compat_detail::cfg().DOUBLE_STAR_PX)
#define BASE_FLUX (beast_compat_detail::cfg().BASE_FLUX)
#define IMAGE_VARIANCE (beast_compat_detail::cfg().IMAGE_VARIANCE)
#define THRESH_FACTOR (beast_compat_detail::cfg().THRESH_FACTOR)
#define POS_VARIANCE (beast_compat_detail::cfg().POS_VARIANCE)
#define POS_ERR_SIGMA (beast_compat_detail::cfg().POS_ERR_SIGMA)
#define PSF_SIGMA (beast_compat_detail::cfg().PSF_SIGMA)
#define MAXFOV (beast_compat_detail::cfg().MAXFOV)
#define MINFOV (beast_compat_detail::cfg().MINFOV)
#define MATCH_VALUE (beast_compat_detail::cfg().MATCH_VALUE)
#define PIXX_TANGENT (beast_compat_detail::cfg().PIXX_TANGENT)
#define PIXY_TANGENT (beast_compat_detail::cfg().PIXY_TANGENT)

struct star {
    float x, y, z, photons;
    int32_t star_idx, id, unreliable;
    float sigma_sq, px, py;

    star()
        : x(0)
        , y(0)
        , z(0)
        , photons(0)
        , star_idx(-1)
        , id(-1)
        , unreliable(0)
        , sigma_sq(0)
        , px(0)
        , py(0) { }

    star(float px_, float py_, float photons_, int id_) {
        Star s  = ost_make_img_star(&beast_compat_detail::cfg(), px_, py_, photons_, id_);
        x       = s.x;
        y       = s.y;
        z       = s.z;
        photons = s.flux;
        star_idx   = s.star_idx;
        id         = s.id;
        unreliable = s.unreliable;
        sigma_sq   = s.sigma_sq;
        px         = s.px;
        py         = s.py;
    }

    star(float x_, float y_, float z_, float photons_, int id_) {
        Star s = ost_make_db_star(&beast_compat_detail::cfg(), x_, y_, z_, photons_, id_);
        x      = s.x;
        y      = s.y;
        z      = s.z;
        photons    = s.flux;
        star_idx   = s.star_idx;
        id         = s.id;
        unreliable = s.unreliable;
        sigma_sq   = s.sigma_sq;
        px         = s.px;
        py         = s.py;
    }

    float dist_arcsec(const star& s) const {
        float a = x * s.y - s.x * y;
        float b = x * s.z - s.x * z;
        float c = y * s.z - s.y * z;
        return (3600 * 180.0f / (float)PI)
            * std::asin(std::sqrt(a * a + b * b + c * c));
    }
};

typedef Constellation constellation;
typedef CPair constellation_pair;

namespace beast_compat_detail {
static inline void mirror_star(const Star& s, ::star& d) {
    d.x          = s.x;
    d.y          = s.y;
    d.z          = s.z;
    d.photons    = s.flux;
    d.star_idx   = s.star_idx;
    d.id         = s.id;
    d.unreliable = s.unreliable;
    d.sigma_sq   = s.sigma_sq;
    d.px         = s.px;
    d.py         = s.py;
}
}

struct star_db {
    StarDB db;
    Star* map;
    int map_size;
    float max_variance;
    int kdsorted;
    std::vector<Star> storage;
    mutable star tmp_star;

    star_db()
        : map(NULL)
        , map_size(0)
        , max_variance(0.0f)
        , kdsorted(0) {
        sync();
    }

    void sync() {
        map             = storage.empty() ? NULL : &storage[0];
        map_size        = (int)storage.size();
        db.v            = map;
        db.n            = map_size;
        db.cap          = (int)storage.capacity();
        db.max_variance = max_variance;
    }

    void add(const Star& s) {
        storage.push_back(s);
        if (max_variance < s.sigma_sq) max_variance = s.sigma_sq;
        sync();
    }

    star_db& operator+=(const star& s) {
        Star cs = ost_make_db_star(&beast_compat_detail::cfg(), s.x, s.y, s.z, s.photons,
                                   s.id);
        cs.px   = s.px;
        cs.py   = s.py;
        cs.sigma_sq   = s.sigma_sq;
        cs.unreliable = s.unreliable;
        cs.star_idx   = (int)storage.size();
        add(cs);
        return *this;
    }

    void add_star(float x, float y, float z, float mag, int id) {
        Star s     = ost_make_db_star(&beast_compat_detail::cfg(), x, y, z,
                                      beast_compat_detail::flux_from_mag(mag), id);
        s.star_idx = (int)storage.size();
        add(s);
    }

    void add_star(float px, float py, float mag, int id) {
        Star s     = ost_make_img_star(&beast_compat_detail::cfg(), px, py,
                                       beast_compat_detail::flux_from_mag(mag), id);
        s.star_idx = (int)storage.size();
        add(s);
    }

    void load_catalog(const char* path = "hip_main.dat", float year = 1991.25f) {
        FILE* f = std::fopen(path, "r");
        if (!f) std::abort();
        std::fclose(f);

        size_t star_cap = std::max<size_t>(storage.capacity(), 120000);
        size_t key_cap  = 262144;
        for (;;) {
            storage.resize(star_cap);
            ost_star_db_init(&db, &storage[0], (int)storage.size());
            std::vector<uint64_t> keys(key_cap);
            if (ost_load_catalog(&beast_compat_detail::cfg(), &db, path, year,
                                 &keys[0], keys.size()) == 0)
                break;
            star_cap *= 2;
            key_cap *= 2;
        }
        storage.resize(db.n);
        max_variance = db.max_variance;
        sync();
    }

    star* get_star(int i) const {
        beast_compat_detail::mirror_star(storage[(size_t)i], tmp_star);
        return &tmp_star;
    }

    star_db* copy_n_brightest(int n) {
        star_db* out = new star_db;
        out->storage.resize(std::max(1, std::min(n, map_size)));
        out->sync();
        std::vector<Star> tmp(std::max(1, map_size));
        if (ost_copy_n_brightest(&out->db, &db, &tmp[0], n) < 0) std::abort();
        out->storage.resize(out->db.n);
        out->max_variance = out->db.max_variance;
        out->sync();
        return out;
    }
};

struct star_query {
    star_db* stars;
    Query q;
    Star* map;
    int n, kdsorted;
    int* kdresults;
    int kdresults_size, kdresults_maxsize;
    signed char* kdmask;
    std::vector<Star> qmap;
    std::vector<int> res_storage;
    std::vector<signed char> mask;

    star_query(star_db* s)
        : stars(s) {
        qmap.resize(std::max(1, s->map_size));
        res_storage.resize((size_t)s->map_size + 1);
        mask.resize((size_t)s->map_size + 1);
        ost_query_init(&q, &s->db, &qmap[0], &res_storage[0], &mask[0]);
        sync();
    }

    void sync() {
        map               = q.map;
        n                 = q.n;
        kdsorted          = q.kdsorted;
        kdresults         = q.kdresults;
        kdresults_size    = q.kdresults_size;
        kdresults_maxsize = q.kdresults_maxsize;
        kdmask            = q.kdmask;
    }

    void kdsort() {
        ost_query_kdsort(&q, &beast_compat_detail::cfg());
        sync();
    }
    void reset_kdresults() {
        ost_query_clear_results(&q);
        sync();
    }
    void undo_kdsearch() {
        ost_query_clear_results(&q);
        sync();
    }
    void reset_kdmask() {
        ost_query_reset_mask(&q);
        sync();
    }
    void kdmask_filter_catalog() {
        ost_query_mask_filter(&q, &beast_compat_detail::cfg());
        sync();
    }
    void kdmask_uniform_density(int min_stars) {
        std::vector<signed char> keep(std::max(1, stars->map_size));
        ost_query_mask_uniform(&q, &beast_compat_detail::cfg(), min_stars, &keep[0]);
        sync();
    }
    void kdsearch(float x, float y, float z, float arcsec, float min_flux) {
        float p[3] = { x, y, z };
        ost_query_search(&q, &beast_compat_detail::cfg(), p, arcsec, min_flux);
        sync();
    }

    int r_size() const { return kdresults_size; }
    void clear_kdresults() { reset_kdresults(); }

    star_db* from_kdmask() {
        star_db* out = new star_db;
        out->storage.resize(std::max(1, stars->map_size));
        out->sync();
        if (ost_db_from_mask(&out->db, &q) < 0) std::abort();
        out->storage.resize(out->db.n);
        out->max_variance = out->db.max_variance;
        out->sync();
        return out;
    }

    star_db* from_kdresults() {
        star_db* out = new star_db;
        out->storage.resize(std::max(1, kdresults_size));
        out->sync();
        if (ost_db_from_results(&out->db, &q) < 0) std::abort();
        out->storage.resize(out->db.n);
        out->max_variance = out->db.max_variance;
        out->sync();
        return out;
    }
};

struct constellation_db {
    constellation* map;
    star_db* stars;
    int map_size;
    CDB cdb;
    Query q;
    std::vector<Star> star_storage, qmap;
    star_query* results;
    std::vector<int> res_storage;
    std::vector<signed char> mask, keep;
    std::vector<Constellation> cmap;

    constellation_db()
        : map(NULL)
        , stars(NULL)
        , map_size(0)
        , results(NULL) { }

    constellation_db(star_db* s, int stars_per_fov, bool from_image)
        : map(NULL)
        , stars(s)
        , map_size(0)
        , results(NULL) {
        build(s, stars_per_fov, from_image);
    }

    void build(star_db* s, int stars_per_fov, bool from_image) {
        stars = s;
        star_storage.resize(std::max(1, s->map_size));
        qmap.resize(std::max(1, s->map_size));
        res_storage.resize((size_t)s->map_size + 1);
        mask.resize((size_t)s->map_size + 1);
        keep.resize(std::max(1, s->map_size));
        size_t cap = from_image ? (size_t)std::max(1, stars_per_fov * (stars_per_fov - 1) / 2)
                                : std::max<size_t>(1024, (size_t)s->map_size * stars_per_fov);
        int rc;
        for (;;) {
            cmap.resize(cap);
            rc = from_image
                ? ost_db_from_image(&cdb, &s->db, &star_storage[0],
                                    (int)star_storage.size(), &q, &qmap[0],
                                    &res_storage[0], &mask[0], &cmap[0], (int)cmap.size(),
                                    stars_per_fov)
                : ost_db_from_catalog(&cdb, &s->db, &star_storage[0],
                                      (int)star_storage.size(), &q, &qmap[0],
                                      &res_storage[0], &mask[0], &cmap[0],
                                      (int)cmap.size(), stars_per_fov,
                                      &beast_compat_detail::cfg(), &keep[0]);
            if (rc == 0) break;
            if (from_image) std::abort();
            cap *= 2;
        }
        map      = cdb.map;
        map_size = cdb.map_size;
        if (results) delete results;
        results    = new star_query(stars);
        results->q = cdb.results;
        results->sync();
    }

    ~constellation_db() { delete results; }
};

struct beast_db {
    star_db* stars;
    star_query* results;
    constellation_db* constellations;
    bool own_stars;

    beast_db()
        : stars(new star_db)
        , results(NULL)
        , constellations(NULL)
        , own_stars(true) {
        stars->load_catalog();
        results = new star_query(stars);
        results->kdmask_filter_catalog();
        results->kdmask_uniform_density(REQUIRED_STARS);
        constellations = new constellation_db(stars, 2 + DB_REDUNDANCY, false);
        results->reset_kdmask();
        results->reset_kdresults();
    }

    beast_db(star_db* s)
        : stars(s)
        , results(NULL)
        , constellations(NULL)
        , own_stars(false) {
        constellations = new constellation_db(stars, MAX_FALSE_STARS + 2, true);
    }

    ~beast_db() {
        delete results;
        delete constellations;
        if (own_stars) delete stars;
    }
};

struct match_result {
    CPair match;
    int32_t* map;
    float R11, R12, R13, R21, R22, R23, R31, R32, R33;
    MatchResult c;
    std::vector<int> map_storage;

    match_result()
        : map(NULL) {
        std::memset(&c, 0, sizeof(c));
    }

    void sync_from_c() {
        match = c.match;
        map_storage.assign(c.map, c.map + c.map_size);
        map = map_storage.empty() ? NULL : &map_storage[0];
        R11 = c.R[0][0];
        R21 = c.R[0][1];
        R31 = c.R[0][2];
        R12 = c.R[1][0];
        R22 = c.R[1][1];
        R32 = c.R[1][2];
        R13 = c.R[2][0];
        R23 = c.R[2][1];
        R33 = c.R[2][2];
    }

    star_db* from_match() {
        star_db* out = new star_db;
        int n        = c.img ? c.img->stars.n : 0;
        out->storage.resize(std::max(1, n));
        for (int i = 0; i < n; i++) {
            out->storage[(size_t)i] = c.img->stars.v[i];
            int dbi                 = (i < c.map_size) ? c.map[i] : -1;
            if (dbi >= 0)
                out->storage[(size_t)i] = c.db->stars.v[dbi];
            else
                out->storage[(size_t)i].id = -1;
            out->storage[(size_t)i].star_idx = i;
        }
        out->max_variance = c.db ? c.db->stars.max_variance : 0.0f;
        out->sync();
        return out;
    }
};

namespace beast_compat_detail {
struct match_context {
    std::vector<CPair> candidates;
    std::vector<int> fov_mask, collision, match_map, work_map;
    std::vector<float> fov_px, fov_py, scores;

    void prepare(int n, size_t candidate_cap, size_t collision_cap) {
        candidates.resize(std::max<size_t>(1, candidate_cap));
        fov_mask.resize((size_t)IMG_X * IMG_Y);
        collision.resize(std::max<size_t>(1, collision_cap));
        fov_px.resize(std::max(1, n));
        fov_py.resize(std::max(1, n));
        scores.resize(std::max(1, n));
        match_map.resize(std::max(1, n));
        work_map.resize(std::max(1, n));
    }
};

static inline match_context& match_scratch() {
    static match_context c;
    return c;
}
}

struct db_match {
    float p_match;
    int32_t* map;
    match_result* winner;
    std::vector<int32_t> id_map;

    db_match(beast_db* db, beast_db* img)
        : p_match(0.0f)
        , map(NULL)
        , winner(new match_result) {
        run(beast_compat_detail::match_scratch(), db->constellations, img->constellations,
            db->stars, img->stars);
    }

    db_match(constellation_db* db, constellation_db* img)
        : p_match(0.0f)
        , map(NULL)
        , winner(new match_result) {
        run(beast_compat_detail::match_scratch(), db, img, db->stars, img->stars);
    }

private:
    void run(beast_compat_detail::match_context& ctx, constellation_db* db,
             constellation_db* img, star_db*, star_db* img_stars) {
        int n = img_stars->map_size;
        size_t candidate_cap = std::max(ctx.candidates.size(),
                                        std::max<size_t>(65536, (size_t)img->map_size * 16));
        size_t collision_cap = std::max(ctx.collision.size(),
                                        std::max<size_t>(16384, (size_t)n * 8));
        MatchWork work;
        for (;;) {
            ctx.prepare(n, candidate_cap, collision_cap);
            ost_match_work_init(&work, &ctx.candidates[0], (int)ctx.candidates.size(),
                                &ctx.fov_mask[0], &ctx.collision[0],
                                (int)ctx.collision.size(), &ctx.fov_px[0],
                                &ctx.fov_py[0], &ctx.scores[0], &ctx.match_map[0],
                                &ctx.work_map[0]);
            if (ost_db_match(&db->cdb, &img->cdb, &winner->c,
                             &beast_compat_detail::cfg(), &work, &p_match) == 0)
                break;
            if (db->cdb.results.kdsorted)
                ost_query_clear_results(&db->cdb.results);
            candidate_cap *= 2;
            collision_cap *= 2;
        }
        winner->sync_from_c();
        int map_len = n;
        for (int i = 0; i < n; i++)
            if (map_len <= img_stars->storage[(size_t)i].star_idx)
                map_len = img_stars->storage[(size_t)i].star_idx + 1;
        id_map.assign((size_t)map_len, -1);
        for (int i = 0; i < n && i < winner->c.map_size; i++) {
            int o = winner->c.map[i];
            if (o >= 0)
                id_map[(size_t)img_stars->storage[(size_t)i].star_idx]
                    = winner->c.db->stars.v[o].id;
        }
        map = id_map.empty() ? NULL : &id_map[0];
    }

public:
    ~db_match() { delete winner; }
};

static inline int load_config(const char* filename) {
    return ost_load_config(&beast_compat_detail::cfg(), filename);
}

#endif /* __cplusplus */
#endif
