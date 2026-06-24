"""Native Python interface for OpenStartracker's OST core.

The C tests use ``ost.h`` by explicitly owning a config, star databases,
queries, constellation databases, and match workspaces.  This module exposes
that same model with Python-owned storage and a small ``Tracker`` convenience
wrapper for the test flow.
"""

import ctypes as _ct
import math as _math
import os as _os

_here = _os.path.dirname(__file__)
_lib = _ct.CDLL(_os.path.join(_here, "_ost.so"))

__all__ = [
    "Config",
    "Star",
    "StarDB",
    "Query",
    "ConstellationDB",
    "ConstellationIndex",
    "OST_CONSTELLATION_PAIRDIST",
    "OST_CONSTELLATION_CROSSRATIO",
    "MatchResult",
    "ImagePipeline",
    "Tracker",
    "load_config",
    "match",
    "match_constellations",
    "read_png_rgba",
]

MAX_STARS = 1000
MAX_CAT = 120000
KEY_CAP = 262144
MAX_CDB = 600000
MAX_CANDIDATES = 65536
MAX_COLLISION = 16384

class Config(_ct.Structure):
    """Camera/tracker configuration loaded from a calibration file."""
    _fields_ = [
        ("IMG_X", _ct.c_int), ("IMG_Y", _ct.c_int),
        ("MAX_FALSE_STARS", _ct.c_int), ("DB_REDUNDANCY", _ct.c_int),
        ("REQUIRED_STARS", _ct.c_int), ("KDBUCKET_SIZE", _ct.c_int),
        ("PIXSCALE", _ct.c_float), ("DOUBLE_STAR_PX", _ct.c_float),
        ("BASE_FLUX", _ct.c_float), ("IMAGE_VARIANCE", _ct.c_float),
        ("THRESH_FACTOR", _ct.c_float), ("POS_VARIANCE", _ct.c_float),
        ("POS_ERR_SIGMA", _ct.c_float), ("PSF_SIGMA", _ct.c_float),
        ("MAXFOV", _ct.c_float),
        ("MINFOV", _ct.c_float), ("MATCH_VALUE", _ct.c_float),
        ("PIXX_TANGENT", _ct.c_float), ("PIXY_TANGENT", _ct.c_float),
    ]

    @classmethod
    def load(cls, filename):
        """Load and return a ``Config`` from ``filename``."""
        cfg = cls()
        if _lib.ost_load_config(_ct.byref(cfg), _b(filename)) < 0:
            raise OSError(filename)
        return cfg

class Star(_ct.Structure):
    """Catalog or image star in the OST C layout."""
    _fields_ = [
        ("x", _ct.c_float), ("y", _ct.c_float), ("z", _ct.c_float),
        ("flux", _ct.c_float), ("px", _ct.c_float), ("py", _ct.c_float),
        ("sigma_sq", _ct.c_float), ("id", _ct.c_int),
        ("star_idx", _ct.c_int), ("unreliable", _ct.c_int),
    ]

    @classmethod
    def catalog(cls, cfg, x, y, z, flux, id):
        """Create a catalog star from unit-vector coordinates."""
        return _lib.ost_make_db_star(_ct.byref(cfg), x, y, z, flux, id)

    @classmethod
    def image(cls, cfg, px, py, flux, id=-1):
        """Create an image star from center-relative pixel coordinates."""
        return _lib.ost_make_img_star(_ct.byref(cfg), px, py, flux, id)

class _StarDB(_ct.Structure):
    _fields_ = [("v", _ct.POINTER(Star)), ("n", _ct.c_int),
                ("cap", _ct.c_int), ("max_variance", _ct.c_float)]

class _Query(_ct.Structure):
    _fields_ = [("map", _ct.POINTER(Star)), ("n", _ct.c_int),
                ("kdsorted", _ct.c_int),
                ("kdresults", _ct.POINTER(_ct.c_int)),
                ("kdresults_size", _ct.c_int),
                ("kdresults_maxsize", _ct.c_int),
                ("kdmask", _ct.POINTER(_ct.c_byte))]

class _Constellation(_ct.Structure):
    _fields_ = [("p", _ct.c_float), ("s1", _ct.c_int),
                ("s2", _ct.c_int), ("idx", _ct.c_int)]

class _CPair(_ct.Structure):
    _fields_ = [("totalscore", _ct.c_float), ("db_s1", _ct.c_int),
                ("db_s2", _ct.c_int), ("img_s1", _ct.c_int),
                ("img_s2", _ct.c_int)]

class _CDB(_ct.Structure):
    _fields_ = [("stars", _StarDB), ("results", _Query),
                ("map", _ct.POINTER(_Constellation)), ("map_size", _ct.c_int)]

class _ConstellationEdge(_ct.Structure):
    _fields_ = [("star", _ct.c_int)]

OST_CONSTELLATION_PAIRDIST = 0
OST_CONSTELLATION_CROSSRATIO = 1

_DESCRIPTOR_KINDS = {
    None: OST_CONSTELLATION_PAIRDIST,
    "pairdist": OST_CONSTELLATION_PAIRDIST,
    b"pairdist": OST_CONSTELLATION_PAIRDIST,
    "crossratio": OST_CONSTELLATION_CROSSRATIO,
    b"crossratio": OST_CONSTELLATION_CROSSRATIO,
}

class _ConstellationIndex(_ct.Structure):
    _fields_ = [("pair", _ct.POINTER(_CDB)),
                ("map", _ct.POINTER(_ct.c_ubyte)),
                ("map_size", _ct.c_int), ("cap", _ct.c_int),
                ("kd_ready", _ct.c_int), ("kd_bucket", _ct.c_int),
                ("k", _ct.c_int), ("dims", _ct.c_int),
                ("descriptor_kind", _ct.c_int),
                ("record_size", _ct.c_size_t)]

class _StarFov(_ct.Structure):
    pass

_Vec3 = _ct.c_float * 3
class _CMatchResult(_ct.Structure):
    _fields_ = [("match", _CPair), ("R", _Vec3 * 3),
                ("map", _ct.POINTER(_ct.c_int)), ("map_size", _ct.c_int),
                ("db", _ct.POINTER(_CDB)), ("img", _ct.POINTER(_CDB)),
                ("img_mask", _ct.POINTER(_StarFov))]

class _MatchWork(_ct.Structure):
    _fields_ = [("candidates", _ct.POINTER(_CPair)), ("candidate_cap", _ct.c_int),
                ("fov_mask", _ct.POINTER(_ct.c_int)),
                ("collision", _ct.POINTER(_ct.c_int)), ("collision_cap", _ct.c_int),
                ("fov_px", _ct.POINTER(_ct.c_float)),
                ("fov_py", _ct.POINTER(_ct.c_float)),
                ("scores", _ct.POINTER(_ct.c_float)),
                ("match_map", _ct.POINTER(_ct.c_int)),
                ("work_map", _ct.POINTER(_ct.c_int))]

_PConfig = _ct.POINTER(Config)
_PStar = _ct.POINTER(Star)
_PStarDB = _ct.POINTER(_StarDB)
_PQuery = _ct.POINTER(_Query)
_PCDB = _ct.POINTER(_CDB)
_PConstellationIndex = _ct.POINTER(_ConstellationIndex)
_PConstellationEdge = _ct.POINTER(_ConstellationEdge)

_lib.ost_load_config.argtypes = [_PConfig, _ct.c_char_p]
_lib.ost_load_config.restype = _ct.c_int
_lib.ost_make_db_star.argtypes = [_PConfig, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_int]
_lib.ost_make_db_star.restype = Star
_lib.ost_make_img_star.argtypes = [_PConfig, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_int]
_lib.ost_make_img_star.restype = Star
_lib.ost_star_db_init.argtypes = [_PStarDB, _PStar, _ct.c_int]
_lib.ost_db_add.argtypes = [_PStarDB, Star]
_lib.ost_db_add.restype = _ct.c_int
_lib.ost_copy_n_brightest.argtypes = [_PStarDB, _PStarDB, _PStar, _ct.c_int]
_lib.ost_copy_n_brightest.restype = _ct.c_int
_lib.ost_load_catalog.argtypes = [_PConfig, _PStarDB, _ct.c_char_p, _ct.c_float, _ct.POINTER(_ct.c_uint64), _ct.c_size_t]
_lib.ost_load_catalog.restype = _ct.c_int
_lib.ost_query_init.argtypes = [_PQuery, _PStarDB, _PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte)]
_lib.ost_query_reset_mask.argtypes = [_PQuery]
_lib.ost_query_clear_results.argtypes = [_PQuery]
_lib.ost_query_sort_flux.argtypes = [_PQuery]
_lib.ost_query_kdsort.argtypes = [_PQuery, _PConfig]
_lib.ost_query_search.argtypes = [_PQuery, _PConfig, _ct.POINTER(_ct.c_float), _ct.c_float, _ct.c_float]
_lib.ost_query_mask_filter.argtypes = [_PQuery, _PConfig]
_lib.ost_query_mask_uniform.argtypes = [_PQuery, _PConfig, _ct.c_int, _ct.POINTER(_ct.c_byte)]
_lib.ost_db_from_mask.argtypes = [_PStarDB, _PQuery]
_lib.ost_db_from_mask.restype = _ct.c_int
_lib.ost_db_from_results.argtypes = [_PStarDB, _PQuery]
_lib.ost_db_from_results.restype = _ct.c_int
_lib.ost_db_from_image.argtypes = [_PCDB, _PStarDB, _PStar, _ct.c_int, _PQuery, _PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte), _ct.POINTER(_Constellation), _ct.c_int, _ct.c_int]
_lib.ost_db_from_image.restype = _ct.c_int
_lib.ost_db_from_catalog.argtypes = [_PCDB, _PStarDB, _PStar, _ct.c_int, _PQuery, _PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte), _ct.POINTER(_Constellation), _ct.c_int, _ct.c_int, _PConfig, _ct.POINTER(_ct.c_byte)]
_lib.ost_db_from_catalog.restype = _ct.c_int
_lib.ost_match_work_init.argtypes = [_ct.POINTER(_MatchWork), _ct.POINTER(_CPair), _ct.c_int, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int), _ct.c_int, _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int)]
_lib.ost_constellation_record_size.argtypes = [_ct.c_int, _ct.c_int]
_lib.ost_constellation_record_size.restype = _ct.c_size_t
_lib.ost_constellation_index_init.argtypes = [_PConstellationIndex, _PCDB, _ct.c_int, _ct.c_int, _ct.POINTER(_ct.c_ubyte), _ct.c_int]
_lib.ost_constellation_index_init.restype = _ct.c_int
_lib.ost_constellation_count.argtypes = [_PCDB, _ct.c_int, _ct.POINTER(_ct.c_uint64), _ct.POINTER(_ct.c_int), _PConstellationEdge, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int)]
_lib.ost_constellation_count.restype = _ct.c_int
_lib.ost_constellation_index_build.argtypes = [_PConstellationIndex, _ct.POINTER(_ct.c_int), _PConstellationEdge, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int)]
_lib.ost_constellation_index_build.restype = _ct.c_int
_lib.ost_constellation_index_kdsort.argtypes = [_PConstellationIndex]
_lib.ost_db_match_constellations.argtypes = [_PConstellationIndex, _PCDB, _ct.POINTER(_CMatchResult), _PConfig, _ct.POINTER(_MatchWork), _ct.POINTER(_ct.c_float)]
_lib.ost_db_match_constellations.restype = _ct.c_int
class _CCComponent(_ct.Structure):
    _fields_ = [("area", _ct.c_int), ("sum_x", _ct.c_int), ("sum_y", _ct.c_int),
                ("min_x", _ct.c_int), ("max_x", _ct.c_int),
                ("min_y", _ct.c_int), ("max_y", _ct.c_int),
                ("signal", _ct.c_int), ("wsum", _ct.c_double),
                ("wx", _ct.c_double), ("wy", _ct.c_double),
                ("wxx", _ct.c_double), ("wyy", _ct.c_double),
                ("wxy", _ct.c_double), ("eig_min", _ct.c_double)]

class _CCBufferSizes(_ct.Structure):
    _fields_ = [("max_labels", _ct.c_int),
                ("components", _ct.c_size_t), ("parent", _ct.c_size_t),
                ("col_label", _ct.c_size_t), ("active_count", _ct.c_size_t),
                ("reuse_after_row", _ct.c_size_t),
                ("total_bytes", _ct.c_size_t)]

class _CCContext(_ct.Structure):
    _fields_ = [("width", _ct.c_int), ("max_labels", _ct.c_int),
                ("components", _ct.POINTER(_CCComponent)),
                ("parent", _ct.POINTER(_ct.c_int)),
                ("col_label", _ct.POINTER(_ct.c_int)),
                ("active_count", _ct.POINTER(_ct.c_int)),
                ("reuse_after_row", _ct.POINTER(_ct.c_int))]

class _BGConfig(_ct.Structure):
    _fields_ = [("width", _ct.c_int), ("height", _ct.c_int),
                ("tile_size", _ct.c_int), ("map_width", _ct.c_int),
                ("map_height", _ct.c_int), ("max_stars", _ct.c_int),
                ("max_pixel_brightness", _ct.c_int), ("sample_radius", _ct.c_int),
                ("psf_sigma", _ct.c_double),
                ("threshold_sigma", _ct.c_double), ("detect_sigma", _ct.c_double)]

class _BGStats(_ct.Structure):
    _fields_ = [("cfg", _ct.POINTER(_BGConfig)),
                ("mean", _ct.POINTER(_ct.c_double)),
                ("var", _ct.POINTER(_ct.c_double)),
                ("poisson", _ct.POINTER(_ct.c_double)),
                ("x_edge", _ct.POINTER(_ct.c_int)),
                ("y_edge", _ct.POINTER(_ct.c_int))]

class _BGFitStar(_ct.Structure):
    _fields_ = [("xi", _ct.c_int), ("yi", _ct.c_int)]

class _BGFitWorkspace(_ct.Structure):
    _fields_ = [("stars1", _ct.POINTER(_BGFitStar)),
                ("stars2", _ct.POINTER(_BGFitStar)),
                ("params1", _ct.POINTER(_ct.c_double)),
                ("params2", _ct.POINTER(_ct.c_double)),
                ("normal", _ct.POINTER(_ct.c_double)),
                ("rhs", _ct.POINTER(_ct.c_double)),
                ("rhs_solve", _ct.POINTER(_ct.c_double)),
                ("cov_xy", _ct.POINTER(_ct.c_double)),
                ("dropped", _ct.POINTER(_ct.c_double))]

_lib.ost_cc_buffer_sizes.argtypes = [_ct.c_int, _ct.POINTER(_CCBufferSizes)]
_lib.ost_cc_buffer_sizes.restype = _ct.c_int
_lib.ost_cc_init.argtypes = [_ct.POINTER(_CCContext), _ct.c_int,
                             _ct.POINTER(_CCComponent), _ct.POINTER(_ct.c_int),
                             _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int),
                             _ct.POINTER(_ct.c_int)]
_lib.ost_cc_init.restype = _ct.c_int
_lib.ost_png_dimensions.argtypes = [_ct.c_char_p, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int)]
_lib.ost_png_dimensions.restype = _ct.c_int
_lib.ost_png_read_rgba.argtypes = [_ct.c_char_p, _ct.POINTER(_ct.c_ubyte), _ct.c_int, _ct.c_int, _ct.c_int]
_lib.ost_png_read_rgba.restype = _ct.c_int
_lib.ost_bg_config_init.argtypes = [_ct.POINTER(_BGConfig), _ct.c_int, _ct.c_int]
_lib.ost_bg_config_init.restype = _ct.c_int
_lib.ost_bg_rgba_to_gray16.argtypes = [_ct.POINTER(_ct.c_uint16), _ct.POINTER(_ct.c_ubyte), _ct.c_int, _ct.c_int]
_lib.ost_bg_compute_stats.argtypes = [_ct.POINTER(_BGConfig), _ct.POINTER(_ct.c_uint16), _ct.c_int,
                                      _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_double),
                                      _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int),
                                      _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int), _ct.c_int]
_lib.ost_bg_compute_stats.restype = _ct.c_int
_lib.ost_bg_extract_fused.argtypes = [_ct.POINTER(_BGConfig), _ct.POINTER(_ct.c_uint16), _ct.c_int,
                                      _ct.POINTER(_BGStats), _ct.POINTER(_CCContext),
                                      _ct.POINTER(_CCComponent), _ct.c_int]
_lib.ost_bg_extract_fused.restype = _ct.c_int
_lib.ost_bg_fit_stars.argtypes = [_ct.POINTER(_BGConfig), _ct.POINTER(_ct.c_uint16), _ct.c_int,
                                  _ct.POINTER(_BGStats), _ct.POINTER(_CCComponent), _ct.c_int,
                                  _ct.c_int, _ct.POINTER(_BGFitWorkspace), _ct.c_int,
                                  _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_double),
                                  _ct.POINTER(_ct.c_int)]
_lib.ost_bg_fit_stars.restype = _ct.c_int

def _b(path):
    return path if isinstance(path, bytes) else str(path).encode()

def load_config(filename):
    """Load and return a ``Config``."""
    return Config.load(filename)

def read_png_rgba(filename):
    """Read a PNG with libpng and return ``width, height, rgba_bytes``."""
    width = _ct.c_int()
    height = _ct.c_int()
    path = _b(filename)
    if _lib.ost_png_dimensions(path, _ct.byref(width), _ct.byref(height)) < 0:
        raise OSError(filename)
    n = 4 * width.value * height.value
    rgba = (_ct.c_ubyte * n)()
    if _lib.ost_png_read_rgba(path, rgba, width.value, height.value,
                              4 * width.value) < 0:
        raise OSError(filename)
    return width.value, height.value, bytes(rgba)

class StarDB:
    """Resizable Python-owned ``StarDB`` storage for the OST C API."""
    def __init__(self, cap=MAX_STARS):
        self._cap = max(1, int(cap))
        self._stars = (Star * self._cap)()
        self._c = _StarDB()
        _lib.ost_star_db_init(_ct.byref(self._c), self._stars, self._cap)

    def __len__(self):
        return self._c.n

    def __getitem__(self, idx):
        if idx < 0:
            idx += self._c.n
        if idx < 0 or idx >= self._c.n:
            raise IndexError(idx)
        return self._stars[idx]

    @property
    def max_variance(self):
        return self._c.max_variance

    @max_variance.setter
    def max_variance(self, value):
        self._c.max_variance = value

    def _ensure(self, n):
        if n <= self._cap:
            return
        cap = max(n, self._cap * 2)
        stars = (Star * cap)()
        _ct.memmove(stars, self._stars, _ct.sizeof(Star) * self._c.n)
        self._stars = stars
        self._cap = cap
        self._c.v = self._stars
        self._c.cap = cap

    def add(self, star):
        """Append a ``Star``."""
        self._ensure(self._c.n + 1)
        if _lib.ost_db_add(_ct.byref(self._c), star) < 0:
            raise MemoryError("StarDB capacity exhausted")
        return self

    def add_image(self, cfg, px, py, flux, id=-1):
        """Append an image star from center-relative pixels."""
        return self.add(Star.image(cfg, float(px), float(py), float(flux), int(id)))

    def add_catalog(self, cfg, x, y, z, flux, id):
        """Append a catalog star from a unit vector."""
        return self.add(Star.catalog(cfg, float(x), float(y), float(z), float(flux), int(id)))

    def copy(self):
        """Deep-copy this database."""
        out = StarDB(max(1, self._c.n))
        out._c.n = self._c.n
        out._c.max_variance = self._c.max_variance
        _ct.memmove(out._stars, self._stars, _ct.sizeof(Star) * self._c.n)
        return out

    def copy_n_brightest(self, n):
        """Return a database containing the ``n`` brightest stars."""
        out = StarDB(max(1, min(int(n), self._c.n)))
        tmp = (Star * max(1, self._c.n))()
        if _lib.ost_copy_n_brightest(_ct.byref(out._c), _ct.byref(self._c), tmp, int(n)) < 0:
            raise MemoryError("copy_n_brightest")
        return out

    def load_catalog(self, cfg, filename="hip_main.dat", year=1991.25, key_cap=KEY_CAP):
        """Load Hipparcos catalog stars into this database, growing as needed."""
        star_cap = max(self._cap, MAX_CAT)
        while True:
            self._ensure(star_cap)
            self._c.n = 0
            keys = (_ct.c_uint64 * int(key_cap))()
            rc = _lib.ost_load_catalog(_ct.byref(cfg), _ct.byref(self._c),
                                       _b(filename), float(year), keys, key_cap)
            if rc == 0:
                return self
            star_cap *= 2
            key_cap *= 2

    @classmethod
    def from_catalog(cls, cfg, filename="hip_main.dat", year=1991.25):
        db = cls(MAX_CAT)
        return db.load_catalog(cfg, filename, year)

    @classmethod
    def from_measurements(cls, cfg, measurements, ids_from_index=False):
        """Build image stars from ``(x, y, mag)`` image measurements."""
        db = cls(max(1, len(measurements)))
        for i, m in enumerate(measurements):
            x, y, mag = m
            flux = cfg.BASE_FLUX * _math.pow(10.0, -float(mag) / 2.5)
            # PNG/image measurements use the same y-down pixel convention as
            # legacy startracker.py; simulator CSV tests use a separate y-up path.
            db.add_image(cfg, float(x) - cfg.IMG_X / 2.0,
                         float(y) - cfg.IMG_Y / 2.0,
                         flux, i if ids_from_index else -1)
        return db

class Query:
    """KD-search map, mask, and result arrays for a ``StarDB``."""
    def __init__(self, stars):
        self.stars = stars
        n = len(stars)
        self.map = (Star * max(1, n))()
        self.results = (_ct.c_int * (n + 1))()
        self.mask = (_ct.c_byte * (n + 1))()
        self._c = _Query()
        _lib.ost_query_init(_ct.byref(self._c), _ct.byref(stars._c),
                            self.map, self.results, self.mask)

    @property
    def result_size(self):
        return self._c.kdresults_size

    def clear_results(self):
        _lib.ost_query_clear_results(_ct.byref(self._c))

    def reset_mask(self):
        _lib.ost_query_reset_mask(_ct.byref(self._c))

    def sort_by_flux(self):
        _lib.ost_query_sort_flux(_ct.byref(self._c))

    def kdsort(self, cfg):
        _lib.ost_query_kdsort(_ct.byref(self._c), _ct.byref(cfg))

    def search(self, cfg, vector, arcsec, min_flux=0.0):
        p = (_ct.c_float * 3)(float(vector[0]), float(vector[1]), float(vector[2]))
        _lib.ost_query_search(_ct.byref(self._c), _ct.byref(cfg), p,
                              float(arcsec), float(min_flux))

    def mask_filter_catalog(self, cfg):
        _lib.ost_query_mask_filter(_ct.byref(self._c), _ct.byref(cfg))

    def mask_uniform_density(self, cfg, min_stars):
        keep = (_ct.c_byte * max(1, len(self.stars)))()
        _lib.ost_query_mask_uniform(_ct.byref(self._c), _ct.byref(cfg),
                                    int(min_stars), keep)
        self._keep = keep

    def from_mask(self):
        out = StarDB(max(1, len(self.stars)))
        out.max_variance = self.stars.max_variance
        if _lib.ost_db_from_mask(_ct.byref(out._c), _ct.byref(self._c)) < 0:
            raise MemoryError("db_from_mask")
        return out

    def from_results(self):
        out = StarDB(max(1, self._c.kdresults_size))
        out.max_variance = self.stars.max_variance
        if _lib.ost_db_from_results(_ct.byref(out._c), _ct.byref(self._c)) < 0:
            raise MemoryError("db_from_results")
        return out

class ConstellationDB:
    """Pairwise constellation database used by the OST matcher."""
    def __init__(self, stars, stars_per_fov, from_image=False, cfg=None, map_cap=None):
        self.stars = stars.copy()
        self.results = Query(self.stars)
        self._c = _CDB()
        self._keep = None
        if map_cap is None:
            n = len(self.stars)
            map_cap = max(1, n * (n - 1) // 2) if from_image else max(MAX_CDB, n * int(stars_per_fov), 1)
        while True:
            self._map = (_Constellation * int(map_cap))()
            if from_image:
                rc = _lib.ost_db_from_image(_ct.byref(self._c), _ct.byref(stars._c),
                                            self.stars._stars, self.stars._cap,
                                            _ct.byref(self.results._c), self.results.map,
                                            self.results.results, self.results.mask,
                                            self._map, map_cap, int(stars_per_fov))
            else:
                if cfg is None:
                    raise TypeError("cfg is required for catalog constellation DBs")
                self._keep = (_ct.c_byte * max(1, len(self.stars)))()
                rc = _lib.ost_db_from_catalog(_ct.byref(self._c), _ct.byref(stars._c),
                                              self.stars._stars, self.stars._cap,
                                              _ct.byref(self.results._c), self.results.map,
                                              self.results.results, self.results.mask,
                                              self._map, map_cap, int(stars_per_fov),
                                              _ct.byref(cfg), self._keep)
            if rc == 0:
                break
            if from_image:
                raise MemoryError("image constellation map too small")
            map_cap *= 2
        self.stars._c = self._c.stars
        self.results._c = self._c.results

    @classmethod
    def from_image(cls, stars, stars_per_fov, map_cap=None):
        return cls(stars, stars_per_fov, True, None, map_cap)

    @classmethod
    def from_catalog(cls, stars, cfg, stars_per_fov, map_cap=None):
        return cls(stars, stars_per_fov, False, cfg, map_cap)

    def _sync_results(self):
        self._c.results = self.results._c

class ConstellationIndex:
    """C-backed K-star constellation index over an OST pair catalog."""
    def __init__(self, pair_db, k=3, descriptor="pairdist"):
        self.pair_db = pair_db
        self.k = int(k)
        try:
            self.descriptor_kind = _DESCRIPTOR_KINDS[descriptor]
        except KeyError:
            raise ValueError("unsupported constellation descriptor") from None
        nstars = len(pair_db.stars)
        npairs = int(pair_db._c.map_size)
        self.off = (_ct.c_int * max(1, nstars + 1))()
        self.edges = (_ConstellationEdge * max(1, 2 * npairs))()
        self.tmp = (_ct.c_int * max(1, nstars))()
        self.common = (_ct.c_int * max(1, nstars))()
        count = _ct.c_uint64(0)
        if _lib.ost_constellation_count(
                _ct.byref(pair_db._c), self.k, _ct.byref(count), self.off,
                self.edges, self.tmp, self.common) < 0:
            raise RuntimeError("constellation count failed")
        rec = int(_lib.ost_constellation_record_size(
            self.k, self.descriptor_kind))
        if count.value > (1 << 31) - 1 or rec <= 0:
            raise MemoryError("constellation index too large")
        self._storage = (_ct.c_ubyte * max(1, int(count.value) * rec))()
        self._c = _ConstellationIndex()
        if _lib.ost_constellation_index_init(
                _ct.byref(self._c), _ct.byref(pair_db._c), self.k,
                self.descriptor_kind, self._storage, int(count.value)) < 0:
            raise RuntimeError("constellation index init failed")
        if _lib.ost_constellation_index_build(
                _ct.byref(self._c), self.off, self.edges, self.tmp,
                self.common) < 0:
            raise RuntimeError("constellation index build failed")
        _lib.ost_constellation_index_kdsort(_ct.byref(self._c))

    @property
    def count(self):
        return int(self._c.map_size)

class MatchResult:
    """Result from matching two constellation databases."""
    def __init__(self, c_result, p_match):
        self.p_match = float(p_match)
        self.score = float(c_result.match.totalscore)
        self.map = [int(c_result.map[i]) for i in range(c_result.map_size)]
        self.rotation = tuple(tuple(float(c_result.R[r][c]) for c in range(3))
                              for r in range(3))

    def __len__(self):
        return len(self.map)

class _Workspace:
    def __init__(self, cfg, image_stars, image_pairs, candidate_cap=None, collision_cap=None):
        n = max(1, int(image_stars))
        candidate_cap = max(MAX_CANDIDATES, int(image_pairs) * 16, 1) if candidate_cap is None else candidate_cap
        collision_cap = max(MAX_COLLISION, n * 8, 1) if collision_cap is None else collision_cap
        self.candidates = (_CPair * candidate_cap)()
        self.fov_mask = (_ct.c_int * (cfg.IMG_X * cfg.IMG_Y))()
        self.collision = (_ct.c_int * collision_cap)()
        self.fov_px = (_ct.c_float * n)()
        self.fov_py = (_ct.c_float * n)()
        self.scores = (_ct.c_float * n)()
        self.match_map = (_ct.c_int * n)()
        self.work_map = (_ct.c_int * n)()
        self._c = _MatchWork()
        _lib.ost_match_work_init(_ct.byref(self._c), self.candidates, candidate_cap,
                                 self.fov_mask, self.collision, collision_cap,
                                 self.fov_px, self.fov_py, self.scores,
                                 self.match_map, self.work_map)

def _match_with(cfg, db, img, call, candidate_scale=16):
    db._sync_results()
    img._sync_results()
    candidate_cap = max(MAX_CANDIDATES, img._c.map_size * candidate_scale, 1)
    collision_cap = max(MAX_COLLISION, len(img.stars) * 8, 1)
    while True:
        work = _Workspace(cfg, len(img.stars), img._c.map_size,
                          candidate_cap, collision_cap)
        c_result = _CMatchResult()
        p_match = _ct.c_float(0.0)
        rc = call(c_result, work, p_match)
        if rc == 0:
            return MatchResult(c_result, p_match.value)
        if db._c.results.kdsorted:
            _lib.ost_query_clear_results(_ct.byref(db._c.results))
        candidate_cap *= 2
        collision_cap *= 2

def match(cfg, db, img):
    """Match catalog/FOV ``db`` against image ``img`` and return ``MatchResult``."""
    return match_constellations(cfg, ConstellationIndex(db, 2), img)

def match_constellations(cfg, constellation_index, img):
    """Match image pairs against a precomputed K-star constellation index."""
    db = constellation_index.pair_db
    return _match_with(cfg, db, img,
        lambda c_result, work, p_match: _lib.ost_db_match_constellations(
            _ct.byref(constellation_index._c), _ct.byref(img._c),
            _ct.byref(c_result), _ct.byref(cfg), _ct.byref(work._c),
            _ct.byref(p_match)))

class ImagePipeline:
    """Image-to-measurement pipeline backed by ``ost_bg_*`` functions."""
    def __init__(self, cfg):
        self.cfg = cfg
        self.bg_cfg = _BGConfig()
        if _lib.ost_bg_config_init(_ct.byref(self.bg_cfg), cfg.IMG_X, cfg.IMG_Y) < 0:
            raise ValueError("invalid image dimensions")
        self.bg_cfg.psf_sigma = cfg.PSF_SIGMA
        sizes = _CCBufferSizes()
        if _lib.ost_cc_buffer_sizes(cfg.IMG_X, _ct.byref(sizes)) < 0:
            raise ValueError("invalid connected-component width")
        pixels = cfg.IMG_X * cfg.IMG_Y
        map_pixels = self.bg_cfg.map_width * self.bg_cfg.map_height
        max_stars = self.bg_cfg.max_stars
        self.rgba = (_ct.c_ubyte * (4 * pixels))()
        self.gray = (_ct.c_uint16 * pixels)()
        self.mean = (_ct.c_double * map_pixels)()
        self.var = (_ct.c_double * map_pixels)()
        self.poisson = (_ct.c_double * map_pixels)()
        self.x_edge = (_ct.c_int * (self.bg_cfg.map_width + 1))()
        self.y_edge = (_ct.c_int * (self.bg_cfg.map_height + 1))()
        self.hist = (_ct.c_int * 65536)()
        self.cc_components = (_CCComponent * sizes.components)()
        self.parent = (_ct.c_int * sizes.parent)()
        self.col_label = (_ct.c_int * sizes.col_label)()
        self.active_count = (_ct.c_int * sizes.active_count)()
        self.reuse_after_row = (_ct.c_int * sizes.reuse_after_row)()
        self.components = (_CCComponent * max_stars)()
        self.fit_stars1 = (_BGFitStar * max_stars)()
        self.fit_stars2 = (_BGFitStar * max_stars)()
        self.params1 = (_ct.c_double * (3 * max_stars + 1))()
        self.params2 = (_ct.c_double * (3 * max_stars + 1))()
        self.normal = (_ct.c_double * (6 * max_stars))()
        self.rhs = (_ct.c_double * (3 * max_stars))()
        self.rhs_solve = (_ct.c_double * (3 * max_stars))()
        self.cov_xy = (_ct.c_double * (2 * max_stars))()
        self.dropped_work = (_ct.c_double * (3 * max_stars))()
        self.fit_params = (_ct.c_double * (3 * max_stars + 1))()
        self.fit_cov = (_ct.c_double * (2 * max_stars))()
        self.bg_stats = _BGStats(_ct.pointer(self.bg_cfg), self.mean, self.var,
                                 self.poisson, self.x_edge, self.y_edge)
        self.cc = _CCContext()
        self.fit_work = _BGFitWorkspace(self.fit_stars1, self.fit_stars2,
                                        self.params1, self.params2,
                                        self.normal, self.rhs,
                                        self.rhs_solve, self.cov_xy,
                                        self.dropped_work)
        if _lib.ost_cc_init(_ct.byref(self.cc), cfg.IMG_X,
                            self.cc_components, self.parent,
                            self.col_label, self.active_count,
                            self.reuse_after_row) < 0:
            raise RuntimeError("connected-component init failed")

    def measure_rgba(self, rgba):
        """Return ``((x, y, mag), ...), info`` for an RGBA image buffer."""
        expected = 4 * self.cfg.IMG_X * self.cfg.IMG_Y
        if len(rgba) != expected:
            raise ValueError("RGBA buffer has %d bytes, expected %d" % (len(rgba), expected))
        _ct.memmove(self.rgba, bytes(rgba), expected)
        _lib.ost_bg_rgba_to_gray16(self.gray, self.rgba,
                                   self.cfg.IMG_X, self.cfg.IMG_Y)
        if _lib.ost_bg_compute_stats(_ct.byref(self.bg_cfg), self.gray,
                                     self.cfg.IMG_X, self.mean, self.var,
                                     self.poisson, self.x_edge, self.y_edge,
                                     self.hist, 65536) < 0:
            raise RuntimeError("background stats failed")
        n = _lib.ost_bg_extract_fused(_ct.byref(self.bg_cfg), self.gray,
                                      self.cfg.IMG_X, _ct.byref(self.bg_stats),
                                      _ct.byref(self.cc), self.components,
                                      self.bg_cfg.max_stars)
        if n < 0:
            raise RuntimeError("star extraction failed")
        dropped = _ct.c_int(0)
        fit_n = _lib.ost_bg_fit_stars(_ct.byref(self.bg_cfg), self.gray,
                                      self.cfg.IMG_X, _ct.byref(self.bg_stats),
                                      self.components, n, 3,
                                      _ct.byref(self.fit_work),
                                      self.bg_cfg.max_stars,
                                      self.fit_params, self.fit_cov,
                                      _ct.byref(dropped))
        if fit_n < 0:
            raise RuntimeError("star fitting failed")
        stars = []
        for i in range(fit_n):
            x = self.fit_params[3 * i]
            y = self.fit_params[3 * i + 1]
            flux = self.fit_params[3 * i + 2]
            if not _math.isfinite(x) or not _math.isfinite(y) or not _math.isfinite(flux):
                continue
            if flux <= 0.0 or self.cfg.BASE_FLUX <= 0.0:
                continue
            stars.append((x, y, -2.5 * _math.log10(flux / self.cfg.BASE_FLUX)))
        info = {
            "components": int(n),
            "fitted": int(fit_n),
            "used": len(stars),
            "dropped": int(dropped.value),
            "sigma": float(self.fit_params[3 * fit_n]) if fit_n else 0.0,
        }
        return stars, info

class Tracker:
    """Convenience wrapper matching the C test helper flow."""
    def __init__(self, cfg, k=2, descriptor="pairdist"):
        self.cfg = cfg
        self.k = int(k)
        self.descriptor = descriptor
        self.catalog = None
        self.full_query = None
        self.filtered = None
        self.global_db = None
        self.constellation_index = None

    @classmethod
    def from_config_file(cls, filename, k=2, descriptor="pairdist"):
        return cls(load_config(filename), k, descriptor)

    def prepare_catalog(self, catalog="hip_main.dat", year=1991.25):
        """Load, filter, and pair a catalog like ``prepare_catalog`` in tests."""
        self.catalog = StarDB.from_catalog(self.cfg, catalog, year)
        self.full_query = Query(self.catalog)
        self.full_query.mask_filter_catalog(self.cfg)
        self.full_query.mask_uniform_density(self.cfg, self.cfg.REQUIRED_STARS)
        self.filtered = self.full_query.from_mask()
        self.full_query.reset_mask()
        self.global_db = ConstellationDB.from_catalog(
            self.filtered, self.cfg, 2 + self.cfg.DB_REDUNDANCY)
        self.constellation_index = ConstellationIndex(
            self.global_db, self.k, self.descriptor)
        return self

    def measurements_to_db(self, measurements, ids_from_index=False):
        return StarDB.from_measurements(self.cfg, measurements, ids_from_index)

    def match_catalog_stars(self, measurements, max_false_p=0.1):
        """Return catalog ids for ``(x, y, mag)`` image measurements."""
        if self.global_db is None:
            raise RuntimeError("call prepare_catalog first")
        img = self.measurements_to_db(measurements)
        bright = img.copy_n_brightest(self.cfg.MAX_FALSE_STARS + self.cfg.REQUIRED_STARS)
        img_cdb = ConstellationDB.from_image(bright, self.cfg.MAX_FALSE_STARS + 2)
        winner = match_constellations(self.cfg, self.constellation_index, img_cdb)
        ids = [-1] * len(img)
        if 1.0 - winner.p_match > max_false_p:
            return ids

        min_flux = self.cfg.THRESH_FACTOR * self.cfg.IMAGE_VARIANCE
        self.full_query.search(self.cfg, winner.rotation[0], self.cfg.MAXFOV / 2, min_flux)
        self.global_db.results.search(self.cfg, winner.rotation[0], self.cfg.MAXFOV / 2, min_flux)
        near = self.full_query.from_results()
        fov_db = ConstellationDB.from_image(near, self.global_db.results.result_size)
        self.global_db.results.clear_results()
        self.full_query.clear_results()

        img_full = ConstellationDB.from_image(img, self.cfg.MAX_FALSE_STARS + 2)
        fov_winner = match_constellations(
            self.cfg, ConstellationIndex(fov_db, self.k, self.descriptor),
            img_full)
        for i, dbi in enumerate(fov_winner.map[:len(ids)]):
            ids[i] = int(fov_db.stars[dbi].id) if dbi >= 0 else -1
        return ids

    def match_relative_stars(self, reference, current):
        """Match current measurements to reference indices, like the C self-test."""
        ref = self.measurements_to_db(reference, True)
        cur = self.measurements_to_db(current, False)
        ref_cdb = ConstellationDB.from_image(ref, self.cfg.MAX_FALSE_STARS + 2)
        cur_cdb = ConstellationDB.from_image(cur, self.cfg.MAX_FALSE_STARS + 2)
        winner = match(self.cfg, ref_cdb, cur_cdb)
        ids = [-1] * len(cur)
        if winner.p_match > 0.0:
            for i, dbi in enumerate(winner.map[:len(ids)]):
                ids[i] = int(ref_cdb.stars[dbi].id) if dbi >= 0 else -1
        return ids, winner.p_match
