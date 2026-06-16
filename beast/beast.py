"""Python bindings for OpenStartracker's legacy beast API.

New code should use the native ``ost`` Python module.  This module keeps the
older lower-case beast classes used by ``startracker.py``.

Typical legacy use::

    import beast

    beast.load_config("calibration.txt")
    catalog = beast.star_db()
    catalog.load_catalog("hip_main.dat", 1991.25)
    query = beast.star_query(catalog)
    query.kdmask_filter_catalog()
    query.kdmask_uniform_density(beast.cvar.REQUIRED_STARS)
    filtered = query.from_kdmask()
    catalog_constellations = beast.constellation_db(
        filtered, 2 + beast.cvar.DB_REDUNDANCY, 0)

Image stars are collected in another :class:`star_db`, converted to a
``constellation_db(..., from_image=1)``, and matched with :class:`db_match`.
"""

import ctypes as _ct
import math as _math
import os as _os

_here = _os.path.dirname(__file__)
_lib = _ct.CDLL(_os.path.join(_here, "..", "ost", "_ost.so"))

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

class Star(_ct.Structure):
    _fields_ = [
        ("x", _ct.c_float), ("y", _ct.c_float), ("z", _ct.c_float),
        ("flux", _ct.c_float), ("px", _ct.c_float), ("py", _ct.c_float),
        ("sigma_sq", _ct.c_float), ("id", _ct.c_int),
        ("star_idx", _ct.c_int), ("unreliable", _ct.c_int),
    ]

class StarDB(_ct.Structure):
    _fields_ = [("v", _ct.POINTER(Star)), ("n", _ct.c_int),
                ("cap", _ct.c_int), ("max_variance", _ct.c_float)]

class Query(_ct.Structure):
    _fields_ = [("map", _ct.POINTER(Star)), ("n", _ct.c_int),
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
PStar = _ct.POINTER(Star)
PStarDB = _ct.POINTER(StarDB)
PQuery = _ct.POINTER(Query)
PCDB = _ct.POINTER(CDB)

MAX_STARS = 1000
MAX_CAT = 120000
INITIAL_CDB = 600000
INITIAL_CANDIDATES = 65536
INITIAL_COLLISION = 16384
KEY_CAP = 262144

__all__ = [
    "cvar",
    "load_config",
    "star",
    "star_db",
    "star_query",
    "constellation_db",
    "db_match",
    "match_result",
]

_lib.ost_load_config.argtypes = [PConfig, _ct.c_char_p]
_lib.ost_load_config.restype = _ct.c_int
_lib.ost_make_db_star.argtypes = [PConfig, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_int]
_lib.ost_make_db_star.restype = Star
_lib.ost_make_img_star.argtypes = [PConfig, _ct.c_float, _ct.c_float, _ct.c_float, _ct.c_int]
_lib.ost_make_img_star.restype = Star
_lib.ost_star_db_init.argtypes = [PStarDB, PStar, _ct.c_int]
_lib.ost_db_add.argtypes = [PStarDB, Star]
_lib.ost_db_add.restype = _ct.c_int
_lib.ost_copy_n_brightest.argtypes = [PStarDB, PStarDB, PStar, _ct.c_int]
_lib.ost_copy_n_brightest.restype = _ct.c_int
_lib.ost_load_catalog.argtypes = [PConfig, PStarDB, _ct.c_char_p, _ct.c_float, _ct.POINTER(_ct.c_uint64), _ct.c_size_t]
_lib.ost_load_catalog.restype = _ct.c_int
_lib.ost_query_init.argtypes = [PQuery, PStarDB, PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte)]
_lib.ost_query_sort_flux.argtypes = [PQuery]
_lib.ost_query_kdsort.argtypes = [PQuery, PConfig]
_lib.ost_query_reset_mask.argtypes = [PQuery]
_lib.ost_query_clear_results.argtypes = [PQuery]
_lib.ost_query_search.argtypes = [PQuery, PConfig, _ct.POINTER(_ct.c_float), _ct.c_float, _ct.c_float]
_lib.ost_query_search_range.argtypes = [PQuery, PConfig, _ct.POINTER(_ct.c_float), _ct.c_float, _ct.c_float, _ct.c_int, _ct.c_int, _ct.c_int]
_lib.ost_query_mask_filter.argtypes = [PQuery, PConfig]
_lib.ost_query_mask_uniform.argtypes = [PQuery, PConfig, _ct.c_int, _ct.POINTER(_ct.c_byte)]
_lib.ost_db_from_mask.argtypes = [PStarDB, PQuery]
_lib.ost_db_from_mask.restype = _ct.c_int
_lib.ost_db_from_results.argtypes = [PStarDB, PQuery]
_lib.ost_db_from_results.restype = _ct.c_int
_lib.ost_db_from_image.argtypes = [PCDB, PStarDB, PStar, _ct.c_int, PQuery, PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte), _ct.POINTER(Constellation), _ct.c_int, _ct.c_int]
_lib.ost_db_from_image.restype = _ct.c_int
_lib.ost_db_from_catalog.argtypes = [PCDB, PStarDB, PStar, _ct.c_int, PQuery, PStar, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_byte), _ct.POINTER(Constellation), _ct.c_int, _ct.c_int, PConfig, _ct.POINTER(_ct.c_byte)]
_lib.ost_db_from_catalog.restype = _ct.c_int
_lib.ost_match_work_init.argtypes = [_ct.POINTER(MatchWork), _ct.POINTER(CPair), _ct.c_int, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int), _ct.c_int, _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_float), _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int)]
_lib.ost_db_match.argtypes = [PCDB, PCDB, _ct.POINTER(MatchResultC), PConfig, _ct.POINTER(MatchWork), _ct.POINTER(_ct.c_float)]
_lib.ost_db_match.restype = _ct.c_int

_cfg = Config()
class _CVar:
    """Configuration values loaded by :func:`load_config`.

    Attributes mirror the calibration/configuration fields used by the C core,
    such as ``IMG_X``, ``IMG_Y``, ``PIXSCALE``, ``REQUIRED_STARS``,
    ``MAX_FALSE_STARS``, ``DB_REDUNDANCY``, and matching thresholds.
    """
    def __dir__(self):
        """Return the config field names exposed after loading a config."""
        return [name for name, _ in Config._fields_]

cvar = _CVar()

def _b(path):
    return path if isinstance(path, bytes) else str(path).encode()

def _publish_config():
    for name, _ in Config._fields_:
        setattr(cvar, name, getattr(_cfg, name))

def load_config(filename):
    """Load a camera calibration/configuration file.

    After this succeeds, the parsed values are available as attributes on
    :data:`cvar`.  Load a config before constructing image stars, loading a
    catalog, creating constellation databases, or matching.

    Args:
        filename: Path to a calibration/configuration text file.

    Raises:
        OSError: If the file cannot be read or parsed by the OST backend.
    """
    if _lib.ost_load_config(_ct.byref(_cfg), _b(filename)) < 0:
        raise OSError(filename)
    _publish_config()

class star:
    """A single catalog or image star.

    ``star(x, y, z, flux, id)`` creates a catalog star when ``id >= 0``.
    ``star(px, py, flux, -1)`` creates an image star from pixel offsets relative
    to the image center; the loaded config supplies the camera model.  Stars
    returned by :meth:`star_db.get_star` are lightweight views into that
    database, so assigning properties updates the database entry.
    """
    __slots__ = ("_db", "_idx", "_own")
    def __init__(self, a=0.0, b=0.0, c=0.0, flux=None, id=-1, _db=None, _idx=None, _own=None):
        """Create a catalog star, image star, or internal database view."""
        self._db = _db
        self._idx = _idx
        if _own is not None:
            self._own = _own
        elif flux is None:
            self._own = Star()
        elif id == -1:
            self._own = _lib.ost_make_img_star(
                _ct.byref(_cfg),
                float(a),
                float(b),
                float(c),
                int(flux),
            )
        else:
            self._own = _lib.ost_make_db_star(
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

    def _get(name, doc):
        return property(
            lambda self: getattr(self._c, name),
            lambda self, v: setattr(self._c, name, v),
            doc=doc,
        )
    x = _get("x", "Unit-vector x coordinate.")
    y = _get("y", "Unit-vector y coordinate.")
    z = _get("z", "Unit-vector z coordinate.")
    flux = _get("flux", "Star brightness/flux value.")
    px = _get("px", "Image-space x coordinate in pixels relative to image center.")
    py = _get("py", "Image-space y coordinate in pixels relative to image center.")
    sigma_sq = _get("sigma_sq", "Estimated position variance for this star.")
    id = _get("id", "Catalog identifier, or -1 when unmatched/unknown.")
    star_idx = _get("star_idx", "Original index of this star in its source database.")
    unreliable = _get("unreliable", "Nonzero when the backend marked this star unreliable.")
    def dist_arcsec(self, s):
        """Return angular distance to another star in arcseconds."""
        a = self.x * s.y - s.x * self.y
        b = self.x * s.z - s.x * self.z
        c = self.y * s.z - s.y * self.z
        return (3600 * 180.0 / _math.pi) * _math.asin(_math.sqrt(a*a+b*b+c*c))

class star_db:
    """Resizable database of :class:`star` entries.

    A ``star_db`` owns the contiguous C array used by the OST backend.  It is
    used for both catalog stars and image detections.
    """
    def __init__(self, cap=MAX_STARS):
        """Create an empty database with at least ``cap`` star slots."""
        self._cap = max(1, int(cap))
        self._stars = (Star * self._cap)()
        self._db = StarDB()
        _lib.ost_star_db_init(_ct.byref(self._db), self._stars, self._cap)
    @property
    def max_variance(self):
        """Maximum position variance among stars in this database."""
        return self._db.max_variance

    @max_variance.setter
    def max_variance(self, v):
        self._db.max_variance = v

    def _ensure(self, n):
        if n <= self._cap:
            return
        cap = max(n, self._cap * 2)
        new = (Star * cap)()
        _ct.memmove(new, self._stars, _ct.sizeof(Star) * self._db.n)
        self._stars = new
        self._cap = cap
        self._db.v = self._stars
        self._db.cap = cap
    def size(self):
        """Return the number of stars currently stored."""
        return self._db.n
    def __iadd__(self, s):
        """Append a :class:`star` and return this database."""
        self._ensure(self._db.n + 1)
        cs = s._c if isinstance(s, star) else s
        if _lib.ost_db_add(_ct.byref(self._db), cs) < 0:
            raise MemoryError("star_db full")
        return self
    def get_star(self, idx):
        """Return a mutable view of the star at ``idx``, or ``None`` if empty."""
        if self._db.n <= 0:
            return None
        return star(_db=self, _idx=int(idx))
    def copy(self):
        """Return a deep copy of this database."""
        out = star_db(max(self._db.n, 1))
        out._db.n = self._db.n
        out._db.max_variance = self._db.max_variance
        _ct.memmove(out._stars, self._stars, _ct.sizeof(Star) * self._db.n)
        return out
    def copy_n_brightest(self, n):
        """Return a new database containing the ``n`` brightest stars."""
        out = star_db(max(1, min(int(n), self._db.n)))
        tmp = (Star * max(1, self._db.n))()
        rc = _lib.ost_copy_n_brightest(
            _ct.byref(out._db),
            _ct.byref(self._db),
            tmp,
            int(n),
        )
        if rc < 0:
            raise MemoryError("copy_n_brightest")
        return out
    def load_catalog(self, catalog, year):
        """Load Hipparcos catalog stars for ``year`` into this database.

        The database grows as needed.  Existing entries are replaced.
        ``catalog`` is usually ``"hip_main.dat"``.
        """
        with open(catalog, 'rb'):
            pass
        star_cap = max(self._cap, MAX_CAT)
        key_cap = KEY_CAP
        while True:
            self._ensure(star_cap)
            self._db.n = 0
            keys = (_ct.c_uint64 * key_cap)()
            rc = _lib.ost_load_catalog(
                _ct.byref(_cfg),
                _ct.byref(self._db),
                _b(catalog),
                float(year),
                keys,
                key_cap,
            )
            if rc == 0:
                break
            star_cap *= 2
            key_cap *= 2
    def count(self, s):
        """Return how many stars in this database have the same catalog id."""
        sid = s.id
        return sum(1 for i in range(self._db.n) if self._stars[i].id == sid) if sid >= 0 else 0

class star_query:
    """Search and filtering helper for a :class:`star_db`.

    Queries can sort a database, perform angular kd-tree searches, maintain a
    mask of selected stars, and create new databases from search results.
    """
    def __init__(self, db):
        """Create query workspace for ``db``."""
        self.stars = db
        n = db.size()
        self.map = (Star * max(1, n))()
        self.kdresults = (_ct.c_int * (n + 1))()
        self._kdmask = (_ct.c_byte * (n + 1))()
        self._q = Query()
        _lib.ost_query_init(_ct.byref(self._q), _ct.byref(db._db), self.map, self.kdresults, self._kdmask)
        self.map_size = n
    def is_kdsorted(self):
        """Return whether this query has been kd-tree sorted."""
        return self._q.kdsorted

    def sort(self):
        """Sort the query map by decreasing flux."""
        _lib.ost_query_sort_flux(_ct.byref(self._q))

    def kdsort(self):
        """Sort the query map for kd-tree angular searches."""
        _lib.ost_query_kdsort(_ct.byref(self._q), _ct.byref(_cfg))

    def r_size(self):
        """Return the number of stars in the current kd-search result set."""
        return self._q.kdresults_size

    def get_kdmask(self, i):
        """Return the current mask value for query-map entry ``i``."""
        return self._kdmask[i]

    def reset_kdmask(self):
        """Mark all query-map entries as available."""
        _lib.ost_query_reset_mask(_ct.byref(self._q))

    def clear_kdresults(self):
        """Clear the current kd-search result set."""
        _lib.ost_query_clear_results(_ct.byref(self._q))
    def kdsearch(self, x, y, z, r, min_flux, min=0, max=None, dim=0):
        """Search near unit vector ``(x, y, z)`` within ``r`` arcseconds.

        Results are available through :meth:`from_kdresults` and
        :meth:`r_size`.  ``min_flux`` rejects stars dimmer than the threshold.
        The optional ``min``, ``max``, and ``dim`` arguments search a kd-tree
        subrange and are mainly for backend/internal use.
        """
        if max is None:
            p = (_ct.c_float * 3)(x, y, z)
            _lib.ost_query_search(
                _ct.byref(self._q),
                _ct.byref(_cfg),
                p,
                r,
                min_flux,
            )
        else:
            p = (_ct.c_float * 3)(x, y, z)
            _lib.ost_query_search_range(
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
        """Mask catalog stars that are poor candidates for matching."""
        _lib.ost_query_mask_filter(_ct.byref(self._q), _ct.byref(_cfg))
    def kdmask_uniform_density(self, min_stars_per_fov):
        """Mask stars to keep roughly uniform sky density."""
        keep = (_ct.c_byte * max(1, self.map_size))()
        _lib.ost_query_mask_uniform(_ct.byref(self._q), _ct.byref(_cfg), int(min_stars_per_fov), keep)
    def from_kdmask(self):
        """Return a new database containing entries kept by the current mask."""
        out = star_db(max(1, self.map_size))
        out.max_variance = self.stars.max_variance
        if _lib.ost_db_from_mask(_ct.byref(out._db), _ct.byref(self._q)) < 0:
            raise MemoryError("from_kdmask")
        return out
    def from_kdresults(self):
        """Return a new database containing the current kd-search results."""
        out = star_db(max(1, self._q.kdresults_size))
        out.max_variance = self.stars.max_variance
        if _lib.ost_db_from_results(_ct.byref(out._db), _ct.byref(self._q)) < 0:
            raise MemoryError("from_kdresults")
        return out

class constellation_db:
    """Constellation-pair database used for star matching.

    Catalog databases are normally built with ``from_image=0`` after catalog
    filtering.  Image/FOV databases are normally built with ``from_image=1``
    from detected image stars or nearby catalog stars.  This mirrors the flow
    used by ``startracker.py``: filtered catalog -> catalog constellation DB,
    image detections -> image constellation DB, then :class:`db_match`.
    """
    def __init__(self, s, stars_per_fov, from_image):
        """Build constellation pairs from ``s``.

        Args:
            s: Source :class:`star_db`.
            stars_per_fov: Number of neighbors/pairs to keep per field of view.
            from_image: Nonzero for image/FOV databases, zero for catalogs.
        """
        self.stars = s.copy()
        self.results = star_query(self.stars)
        self._cdb = CDB()
        if from_image:
            ns = min(self.stars.size(), int(stars_per_fov))
            cap = max(1, ns * (ns - 1) // 2)
        else:
            cap = max(INITIAL_CDB, self.stars.size() * int(stars_per_fov), 1)
        while True:
            self._map_cap = cap
            self._map = (Constellation * self._map_cap)()
            if from_image:
                rc = _lib.ost_db_from_image(_ct.byref(self._cdb), _ct.byref(s._db), self.stars._stars, self.stars._cap, _ct.byref(self.results._q), self.results.map, self.results.kdresults, self.results._kdmask, self._map, self._map_cap, int(stars_per_fov))
            else:
                keep = (_ct.c_byte * max(1, self.stars.size()))()
                rc = _lib.ost_db_from_catalog(_ct.byref(self._cdb), _ct.byref(s._db), self.stars._stars, self.stars._cap, _ct.byref(self.results._q), self.results.map, self.results.kdresults, self.results._kdmask, self._map, self._map_cap, int(stars_per_fov), _ct.byref(_cfg), keep)
            if rc == 0:
                break
            if from_image:
                raise MemoryError("constellation_db")
            cap *= 2
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
    """Best match returned by :class:`db_match`.

    Matrix entries ``R11`` .. ``R33`` describe the attitude solution.  The
    result also stores the mapping from image stars to catalog stars and can
    convert that mapping back into a :class:`star_db` with :meth:`from_match`.
    """
    def __init__(self, c, db, img):
        """Wrap a backend match result."""
        self.match = c.match
        self._db = db
        self._img = img
        self._map = [c.map[i] for i in range(c.map_size)]
        self._size = c.map_size
        self.R11, self.R21, self.R31 = c.R[0][0], c.R[0][1], c.R[0][2]
        self.R12, self.R22, self.R32 = c.R[1][0], c.R[1][1], c.R[1][2]
        self.R13, self.R23, self.R33 = c.R[2][0], c.R[2][1], c.R[2][2]
    def size(self):
        """Return the number of image stars considered in this match."""
        return self._size
    def from_match(self):
        """Return image stars with matched catalog entries filled in.

        Unmatched stars have ``id == -1``.  Returns ``None`` if the backend did
        not produce a valid winner.
        """
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
        """Print DEC, RA, and ORIENTATION angles for debugging."""
        print("DEC=%f" % ((360 + _math.asin(self.R31) * 180 / _math.pi) % 360))
        print("RA=%f" % ((360 + _math.atan2(self.R21, self.R11) * 180 / _math.pi) % 360))
        print("ORIENTATION=%f" % (-_math.atan2(self.R32, self.R33) * 180 / _math.pi))

class db_match:
    """Match an image/FOV constellation database against a catalog database.

    After construction, :attr:`p_match` is the estimated match probability and
    :attr:`winner` is a :class:`match_result`.  ``startracker.py`` accepts a
    match when ``p_match`` exceeds its threshold and enough stars were matched.
    """
    def __init__(self, db, img):
        """Run the matcher between catalog/FOV ``db`` and image ``img``."""
        self.p_match = 0.0
        self.winner = None
        db._sync_results()
        img._sync_results()
        n = img.stars.size()
        candidate_cap = max(INITIAL_CANDIDATES, img._cdb.map_size * 16, 1)
        collision_cap = max(INITIAL_COLLISION, n * 8, 1)
        while True:
            candidates = _arr('candidates', CPair, candidate_cap)
            fov_mask = _arr('fov_mask', _ct.c_int, cvar.IMG_X * cvar.IMG_Y)
            collision = _arr('collision', _ct.c_int, collision_cap)
            fov_px = _arr('fov_px', _ct.c_float, n)
            fov_py = _arr('fov_py', _ct.c_float, n)
            scores = _arr('scores', _ct.c_float, n)
            match_map = _arr('match_map', _ct.c_int, n)
            work_map = _arr('work_map', _ct.c_int, n)
            work = MatchWork()
            _lib.ost_match_work_init(_ct.byref(work), candidates, len(candidates), fov_mask, collision, len(collision), fov_px, fov_py, scores, match_map, work_map)
            cwin = MatchResultC()
            p = _ct.c_float(0)
            rc = _lib.ost_db_match(
                _ct.byref(db._cdb),
                _ct.byref(img._cdb),
                _ct.byref(cwin),
                _ct.byref(_cfg),
                _ct.byref(work),
                _ct.byref(p),
            )
            if rc == 0:
                break
            if db._cdb.results.kdsorted:
                _lib.ost_query_clear_results(_ct.byref(db._cdb.results))
            candidate_cap *= 2
            collision_cap *= 2
        self.p_match = p.value
        self.winner = match_result(cwin, db, img)
