import ctypes as _ct
import math as _math
import os as _os

_here = _os.path.dirname(__file__)
_lib = _ct.CDLL(_os.path.join(_here, "_beast_py.so"))

class Config(_ct.Structure):
    _fields_ = [
        ("IMG_X", _ct.c_int), ("IMG_Y", _ct.c_int),
        ("MAX_FALSE_STARS", _ct.c_int), ("DB_REDUNDANCY", _ct.c_int),
        ("REQUIRED_STARS", _ct.c_int), ("KDBUCKET_SIZE", _ct.c_int),
        ("PIXSCALE", _ct.c_float), ("DOUBLE_STAR_PX", _ct.c_float),
        ("BASE_FLUX", _ct.c_float), ("IMAGE_VARIANCE", _ct.c_float),
        ("THRESH_FACTOR", _ct.c_float), ("POS_VARIANCE", _ct.c_float),
        ("POS_ERR_SIGMA", _ct.c_float), ("MAXFOV", _ct.c_float),
        ("MINFOV", _ct.c_float), ("MATCH_VALUE", _ct.c_float),
        ("PIXX_TANGENT", _ct.c_float), ("PIXY_TANGENT", _ct.c_float),
    ]

class BeastStar(_ct.Structure):
    _fields_ = [
        ("x", _ct.c_float), ("y", _ct.c_float), ("z", _ct.c_float),
        ("flux", _ct.c_float), ("px", _ct.c_float), ("py", _ct.c_float),
        ("sigma_sq", _ct.c_float), ("id", _ct.c_int),
        ("star_idx", _ct.c_int), ("unreliable", _ct.c_int),
    ]

class StarDB(_ct.Structure):
    _fields_ = [("v", _ct.POINTER(BeastStar)), ("n", _ct.c_int),
                ("cap", _ct.c_int), ("max_variance", _ct.c_float)]

class Query(_ct.Structure):
    _fields_ = [("map", _ct.POINTER(BeastStar)), ("n", _ct.c_int),
                ("kdsorted", _ct.c_int), ("kdresults", _ct.POINTER(_ct.c_int)),
                ("kdresults_size", _ct.c_int), ("kdresults_maxsize", _ct.c_int),
                ("kdmask", _ct.POINTER(_ct.c_byte))]

class Constellation(_ct.Structure):
    _fields_ = [("p", _ct.c_float), ("s1", _ct.c_int),
                ("s2", _ct.c_int), ("idx", _ct.c_int)]

class CPair(_ct.Structure):
    _fields_ = [("totalscore", _ct.c_float), ("db_s1", _ct.c_int),
                ("db_s2", _ct.c_int), ("img_s1", _ct.c_int),
                ("img_s2", _ct.c_int)]

class CDB(_ct.Structure):
    _fields_ = [("stars", StarDB), ("results", Query),
                ("map", _ct.POINTER(Constellation)), ("map_size", _ct.c_int)]

class StarFov(_ct.Structure):
    pass

Vec3 = _ct.c_float * 3
class MatchResultC(_ct.Structure):
    _fields_ = [("match", CPair), ("R", Vec3 * 3),
                ("map", _ct.POINTER(_ct.c_int)), ("map_size", _ct.c_int),
                ("db", _ct.POINTER(CDB)), ("img", _ct.POINTER(CDB)),
                ("img_mask", _ct.POINTER(StarFov))]

class MatchWork(_ct.Structure):
    _fields_ = [("candidates", _ct.POINTER(CPair)), ("candidate_cap", _ct.c_int),
                ("fov_mask", _ct.POINTER(_ct.c_int)),
                ("collision", _ct.POINTER(_ct.c_int)), ("collision_cap", _ct.c_int),
                ("fov_px", _ct.POINTER(_ct.c_float)),
                ("fov_py", _ct.POINTER(_ct.c_float)),
                ("scores", _ct.POINTER(_ct.c_float)),
                ("match_map", _ct.POINTER(_ct.c_int)),
                ("work_map", _ct.POINTER(_ct.c_int))]

PConfig = _ct.POINTER(Config)
PStar = _ct.POINTER(BeastStar)
PStarDB = _ct.POINTER(StarDB)
PQuery = _ct.POINTER(Query)
PCDB = _ct.POINTER(CDB)

MAX_STARS = 1000
MAX_CAT = 120000
MAX_CDB = 600000
MAX_CANDIDATES = 65536
MAX_COLLISION = 16384
KEY_CAP = 262144

_lib.beast_load_config.argtypes = [PConfig, _ct.c_char_p]
_lib.beast_load_config.restype = _ct.c_int
_lib.beast_make_db_star.argtypes = [PConfig, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_int]
_lib.beast_make_db_star.restype = BeastStar
_lib.beast_make_img_star.argtypes = [PConfig, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_int]
_lib.beast_make_img_star.restype = BeastStar
_lib.beast_star_db_init.argtypes = [PStarDB, PStar, _ct.c_int]
_lib.beast_db_add.argtypes = [PStarDB, BeastStar]
_lib.beast_db_add.restype = _ct.c_int
_lib.beast_copy_n_brightest.argtypes = [PStarDB, PStarDB, PStar, _ct.c_int]
_lib.beast_copy_n_brightest.restype = _ct.c_int
_lib.beast_load_catalog.argtypes = [PConfig, PStarDB, _ct.c_char_p, _ct.c_float, _ct.POINTER(_ct.c_uint64)]
_lib.beast_load_catalog.restype = _ct.c_int
_lib.beast_query_init.argtypes = [PQuery, PStarDB, PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte)]
_lib.beast_query_sort_flux.argtypes = [PQuery]
_lib.beast_query_kdsort.argtypes = [PQuery, PConfig]
_lib.beast_query_reset_mask.argtypes = [PQuery]
_lib.beast_query_clear_results.argtypes = [PQuery]
_lib.beast_query_search.argtypes = [PQuery, PConfig, _ct.POINTER(_ct.c_float), _ct.c_float, _ct.c_float]
_lib.beast_query_search_range.argtypes = [PQuery, PConfig, _ct.POINTER(_ct.c_float), _ct.c_float, _ct.c_float, _ct.c_int, _ct.c_int, _ct.c_int]
_lib.beast_query_mask_filter.argtypes = [PQuery, PConfig]
_lib.beast_query_mask_uniform.argtypes = [PQuery, PConfig, _ct.c_int, _ct.POINTER(_ct.c_byte)]
_lib.beast_db_from_mask.argtypes = [PStarDB, PQuery]
_lib.beast_db_from_mask.restype = _ct.c_int
_lib.beast_db_from_results.argtypes = [PStarDB, PQuery]
_lib.beast_db_from_results.restype = _ct.c_int
_lib.beast_db_from_image.argtypes = [PCDB, PStarDB, PStar, _ct.c_int, PQuery, PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte), _ct.POINTER(Constellation), _ct.c_int, _ct.c_int]
_lib.beast_db_from_image.restype = _ct.c_int
_lib.beast_db_from_catalog.argtypes = [PCDB, PStarDB, PStar, _ct.c_int, PQuery, PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte), _ct.POINTER(Constellation), _ct.c_int, _ct.c_int, PConfig, _ct.POINTER(_ct.c_byte)]
_lib.beast_db_from_catalog.restype = _ct.c_int
_lib.beast_match_work_init.argtypes = [_ct.POINTER(MatchWork), _ct.POINTER(CPair), _ct.c_int, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int), _ct.c_int, _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int)]
_lib.beast_db_match.argtypes = [PCDB, PCDB, _ct.POINTER(MatchResultC), PConfig, _ct.POINTER(MatchWork), _ct.POINTER(_ct.c_float)]
_lib.beast_db_match.restype = _ct.c_int

_cfg = Config()
class _CVar:
    pass
cvar = _CVar()

def _b(path):
    return path if isinstance(path, bytes) else str(path).encode()

def _publish_config():
    for name, _ in Config._fields_:
        setattr(cvar, name, getattr(_cfg, name))

def load_config(filename):
    if _lib.beast_load_config(_ct.byref(_cfg), _b(filename)) < 0:
        raise OSError(filename)
    _publish_config()

class star:
    __slots__ = ("_db", "_idx", "_own")
    def __init__(self, a=0.0, b=0.0, c=0.0, flux=None, id=-1, _db=None, _idx=None, _own=None):
        self._db = _db
        self._idx = _idx
        if _own is not None:
            self._own = _own
        elif flux is None:
            self._own = BeastStar()
        elif id == -1:
            self._own = _lib.beast_make_img_star(
                _ct.byref(_cfg),
                float(a),
                float(b),
                float(c),
                int(flux),
            )
        else:
            self._own = _lib.beast_make_db_star(
                _ct.byref(_cfg),
                float(a),
                float(b),
                float(c),
                float(flux),
                int(id),
            )
    @property
    def _c(self):
        return self._db._stars[self._idx] if self._db is not None else self._own
    def _set(self, v):
        if self._db is not None:
            self._db._stars[self._idx] = v
        else:
            self._own = v

    def _get(name):
        return property(
            lambda self: getattr(self._c, name),
            lambda self, v: setattr(self._c, name, v),
        )
    x = _get("x")
    y = _get("y")
    z = _get("z")
    flux = _get("flux")
    px = _get("px")
    py = _get("py")
    sigma_sq = _get("sigma_sq")
    id = _get("id")
    star_idx = _get("star_idx")
    unreliable = _get("unreliable")
    def dist_arcsec(self, s):
        a = self.x * s.y - s.x * self.y
        b = self.x * s.z - s.x * self.z
        c = self.y * s.z - s.y * self.z
        return (3600 * 180.0 / _math.pi) * _math.asin(_math.sqrt(a*a+b*b+c*c))

class star_db:
    def __init__(self, cap=MAX_STARS):
        self._cap = max(1, int(cap))
        self._stars = (BeastStar * self._cap)()
        self._db = StarDB()
        _lib.beast_star_db_init(_ct.byref(self._db), self._stars, self._cap)
    @property
    def max_variance(self):
        return self._db.max_variance

    @max_variance.setter
    def max_variance(self, v):
        self._db.max_variance = v

    def _ensure(self, n):
        if n <= self._cap:
            return
        cap = max(n, self._cap * 2)
        new = (BeastStar * cap)()
        _ct.memmove(new, self._stars, _ct.sizeof(BeastStar) * self._db.n)
        self._stars = new
        self._cap = cap
        self._db.v = self._stars
        self._db.cap = cap
    def size(self):
        return self._db.n
    def __iadd__(self, s):
        self._ensure(self._db.n + 1)
        cs = s._c if isinstance(s, star) else s
        if _lib.beast_db_add(_ct.byref(self._db), cs) < 0:
            raise MemoryError("star_db full")
        return self
    def get_star(self, idx):
        if self._db.n <= 0:
            return None
        return star(_db=self, _idx=int(idx))
    def copy(self):
        out = star_db(max(self._db.n, 1))
        out._db.n = self._db.n
        out._db.max_variance = self._db.max_variance
        _ct.memmove(out._stars, self._stars, _ct.sizeof(BeastStar) * self._db.n)
        return out
    def copy_n_brightest(self, n):
        out = star_db(max(1, min(int(n), self._db.n)))
        tmp = (BeastStar * max(1, self._db.n))()
        rc = _lib.beast_copy_n_brightest(
            _ct.byref(out._db),
            _ct.byref(self._db),
            tmp,
            int(n),
        )
        if rc < 0:
            raise MemoryError("copy_n_brightest")
        return out
    def load_catalog(self, catalog, year):
        self._ensure(MAX_CAT)
        self._db.n = 0
        keys = (_ct.c_uint64 * KEY_CAP)()
        rc = _lib.beast_load_catalog(
            _ct.byref(_cfg),
            _ct.byref(self._db),
            _b(catalog),
            float(year),
            keys,
        )
        if rc < 0:
            raise OSError(catalog)
    def count(self, s):
        sid = s.id
        return sum(1 for i in range(self._db.n) if self._stars[i].id == sid) if sid >= 0 else 0

class star_query:
    def __init__(self, db):
        self.stars = db
        n = db.size()
        self.map = (BeastStar * max(1, n))()
        self.kdresults = (_ct.c_int * (n + 1))()
        self._kdmask = (_ct.c_byte * (n + 1))()
        self._q = Query()
        _lib.beast_query_init(_ct.byref(self._q), _ct.byref(db._db), self.map, self.kdresults, self._kdmask)
        self.map_size = n
    def is_kdsorted(self):
        return self._q.kdsorted

    def sort(self):
        _lib.beast_query_sort_flux(_ct.byref(self._q))

    def kdsort(self):
        _lib.beast_query_kdsort(_ct.byref(self._q), _ct.byref(_cfg))

    def r_size(self):
        return self._q.kdresults_size

    def get_kdmask(self, i):
        return self._kdmask[i]

    def reset_kdmask(self):
        _lib.beast_query_reset_mask(_ct.byref(self._q))

    def clear_kdresults(self):
        _lib.beast_query_clear_results(_ct.byref(self._q))
    def kdsearch(self, x, y, z, r, min_flux, min=0, max=None, dim=0):
        if max is None:
            p = (_ct.c_float * 3)(x, y, z)
            _lib.beast_query_search(
                _ct.byref(self._q),
                _ct.byref(_cfg),
                p,
                r,
                min_flux,
            )
        else:
            p = (_ct.c_float * 3)(x, y, z)
            _lib.beast_query_search_range(
                _ct.byref(self._q),
                _ct.byref(_cfg),
                p,
                r,
                min_flux,
                min,
                max,
                dim,
            )

    def kdmask_filter_catalog(self):
        _lib.beast_query_mask_filter(_ct.byref(self._q), _ct.byref(_cfg))
    def kdmask_uniform_density(self, min_stars_per_fov):
        keep = (_ct.c_byte * max(1, self.map_size))()
        _lib.beast_query_mask_uniform(_ct.byref(self._q), _ct.byref(_cfg), int(min_stars_per_fov), keep)
    def from_kdmask(self):
        out = star_db(max(1, self.map_size))
        out.max_variance = self.stars.max_variance
        if _lib.beast_db_from_mask(_ct.byref(out._db), _ct.byref(self._q)) < 0:
            raise MemoryError("from_kdmask")
        return out
    def from_kdresults(self):
        out = star_db(max(1, self._q.kdresults_size))
        out.max_variance = self.stars.max_variance
        if _lib.beast_db_from_results(_ct.byref(out._db), _ct.byref(self._q)) < 0:
            raise MemoryError("from_kdresults")
        return out

class constellation_db:
    def __init__(self, s, stars_per_fov, from_image):
        self.stars = s.copy()
        self.results = star_query(self.stars)
        self._cdb = CDB()
        if from_image:
            ns = min(self.stars.size(), int(stars_per_fov))
            cap = ns * (ns - 1) // 2
        else:
            cap = MAX_CDB
        self._map_cap = max(1, cap)
        self._map = (Constellation * self._map_cap)()
        if from_image:
            rc = _lib.beast_db_from_image(_ct.byref(self._cdb), _ct.byref(s._db), self.stars._stars, self.stars._cap, _ct.byref(self.results._q), self.results.map, self.results.kdresults, self.results._kdmask, self._map, self._map_cap, int(stars_per_fov))
        else:
            keep = (_ct.c_byte * max(1, self.stars.size()))()
            rc = _lib.beast_db_from_catalog(_ct.byref(self._cdb), _ct.byref(s._db), self.stars._stars, self.stars._cap, _ct.byref(self.results._q), self.results.map, self.results.kdresults, self.results._kdmask, self._map, self._map_cap, int(stars_per_fov), _ct.byref(_cfg), keep)
        if rc < 0:
            raise MemoryError("constellation_db")
        self.stars._db = self._cdb.stars
        self.results._q = self._cdb.results
        self.map_size = self._cdb.map_size
        self.map = self._map
    def _sync_results(self):
        self._cdb.results = self.results._q

_workspace = {}
def _arr(key, typ, n):
    a = _workspace.get(key)
    if a is None or len(a) < n:
        a = (typ * max(1, n))()
        _workspace[key] = a
    return a

class match_result:
    def __init__(self, c, db, img):
        self.match = c.match
        self._db = db
        self._img = img
        self._map = [c.map[i] for i in range(c.map_size)]
        self._size = c.map_size
        self.R11, self.R21, self.R31 = c.R[0][0], c.R[0][1], c.R[0][2]
        self.R12, self.R22, self.R32 = c.R[1][0], c.R[1][1], c.R[1][2]
        self.R13, self.R23, self.R33 = c.R[2][0], c.R[2][1], c.R[2][2]
    def size(self):
        return self._size
    def from_match(self):
        if self.match.totalscore < -3e38:
            return None
        out = self._img.stars.copy()
        out.max_variance = self._db.stars.max_variance
        for n, dbi in enumerate(self._map):
            dst = self._img.stars._stars[n].star_idx
            if dbi >= 0:
                out._stars[dst] = self._db.stars._stars[dbi]
            else:
                out._stars[dst].id = -1
        return out

    def print_ori(self):
        print("DEC=%f" % ((360 + _math.asin(self.R31) * 180 / _math.pi) % 360))
        print("RA=%f" % ((360 + _math.atan2(self.R21, self.R11) * 180 / _math.pi) % 360))
        print("ORIENTATION=%f" % (-_math.atan2(self.R32, self.R33) * 180 / _math.pi))

class db_match:
    def __init__(self, db, img):
        self.p_match = 0.0
        self.winner = None
        db._sync_results()
        img._sync_results()
        n = img.stars.size()
        candidates = _arr('candidates', CPair, MAX_CANDIDATES)
        fov_mask = _arr('fov_mask', _ct.c_int, cvar.IMG_X * cvar.IMG_Y)
        collision = _arr('collision', _ct.c_int, MAX_COLLISION)
        fov_px = _arr('fov_px', _ct.c_float, n)
        fov_py = _arr('fov_py', _ct.c_float, n)
        scores = _arr('scores', _ct.c_float, n)
        match_map = _arr('match_map', _ct.c_int, n)
        work_map = _arr('work_map', _ct.c_int, n)
        work = MatchWork()
        _lib.beast_match_work_init(_ct.byref(work), candidates, len(candidates), fov_mask, collision, len(collision), fov_px, fov_py, scores, match_map, work_map)
        cwin = MatchResultC()
        p = _ct.c_float(0)
        rc = _lib.beast_db_match(
            _ct.byref(db._cdb),
            _ct.byref(img._cdb),
            _ct.byref(cwin),
            _ct.byref(_cfg),
            _ct.byref(work),
            _ct.byref(p),
        )
        if rc < 0:
            raise RuntimeError("db_match")
        self.p_match = p.value
        self.winner = match_result(cwin, db, img)
