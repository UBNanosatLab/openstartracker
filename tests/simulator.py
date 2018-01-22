#!/usr/bin/python
# coding: utf-8

# Copyright (c) 2016 Joerg H. Mueller <nexyon@gmail.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, distribute with modifications, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE ABOVE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# Except as contained in this notice, the name(s) of the above copyright
# holders shall not be used in advertising or otherwise to promote the
# sale, use or other dealings in this Software without prior written
# authorization.

# Imports
# -------

import pandas as pd
import numpy as np
import math
import sys

BASE_PHOTONS = 19100 # photoelectrons per mm^2 and second of a magnitude 0 G2 star
            

# Geometric Transformations
# -------------------------

def angles_to_vector(azimuth, altitude):
    """Transforms azimuth altitude representation to unit-length vectors."""
    caz = np.cos(azimuth)
    saz = np.sin(azimuth)

    cal = np.cos(altitude)
    sal = np.sin(altitude)

    x = caz * cal
    y = saz * cal
    z = sal

    return np.array([x, y, z]).transpose()

def vector_to_angles(vectors):
    """Transforms unit-length vectors to the azimuth altitude representation."""
    x, y, z = split_vectors(vectors)

    az = np.arctan2(y, x)
    alt = np.arcsin(z)

    return np.array([az, alt])

def split_vectors(vectors):
    """Splits vectors into their x, y and z components."""
    return vectors[..., 0], vectors[..., 1], vectors[..., 2]

def randomVectors(num):
    """Generates `num` random three dimensional unit-length vectors."""
    rands = np.random.rand(num, 2)

    ca = np.cos(2. * np.pi * rands[:, 0])
    sa = np.sin(2. * np.pi * rands[:, 0])

    z = rands[:, 1] * 2 - 1
    r = np.sqrt(1 - z * z)

    return np.vstack([r * ca, r * sa, z]).transpose()

def random_matrix():
    """Generates a random 3x3 orientation matrix.
    Based on James Arvo's algorithm from Graphics Gems III, pages 117-120"""
    rands = np.random.rand(3)

    ca = np.cos(2. * np.pi * rands[0])
    sa = np.sin(2. * np.pi * rands[0])

    R = ([[ca, sa, 0],
        [-sa, ca, 0],
        [0, 0, 1.]])

    cb = np.cos(2 * np.pi * rands[1])
    sb = np.sin(2 * np.pi * rands[1])

    sc = np.sqrt(rands[2])
    mc = np.sqrt(1 - rands[2])

    v = np.array([[cb * sc, sb * sc, mc]])

    H = np.eye(3) - 2 * np.dot(v.transpose(), v)

    matrix = -np.dot(H, R)

    return matrix

def add_vector_noise(base_vectors, stddev):
    """Adds Gaussian noise to a list of vectors."""
    num = len(base_vectors)

    v_r = randomVectors(num)

    v_r = np.cross(base_vectors, v_r)

    v_r /= np.linalg.norm(v_r, axis=1).reshape((-1, 1))

    v_r = base_vectors + stddev * (np.random.randn(num, 1) * v_r)

    return v_r / np.linalg.norm(v_r, axis=1).reshape((-1, 1))

# Camera
# ------

class Camera(object):
    def __init__(self, resolution, pixel_ar = 1, principal_point = (0.5, 0.5)):
        self.resolution = resolution
        self.pixel_ar = pixel_ar
        self.principal_point = principal_point

    def from_angles(self, azimuth, altitude):
        """Transforms camera-relative azimuth altitude information to pixel coordinates."""
        res_x, res_y = self.resolution
        pp_x, pp_y = self.principal_point

        ar = self.pixel_ar * res_x / res_y

        theta = np.pi / 2 - altitude

        # x-axis flip: note that we look in the positive z-axis direction,
        # so in a right-handed coordinate system the x-axis goes to the left,
        # but we want it to go to the right, so we flip the x-axis
        # the y-axis still goes up
        # this equals transforming the azimuth as follows:
        alpha = np.pi - azimuth

        r = self.project(theta)
        x = (     np.cos(alpha) * r + pp_x) * res_x
        y = (ar * np.sin(alpha) * r + pp_y) * res_y

        return np.array((y, x)).transpose()

    def to_angles(self, pos):
        """Transforms pixel coordinates to camera-relative azimuth altitude information."""
        res_x, res_y = self.resolution
        pp_x, pp_y = self.principal_point

        ar = self.pixel_ar * res_x / res_y

        x = pos[:, 1]
        y = pos[:, 0]

        x = (x / res_x - pp_x)
        y = (y / res_y - pp_y) / ar

        az = np.pi - np.arctan2(y, x)

        r = np.sqrt(x ** 2 + y ** 2)
        theta = self.unproject(r)

        alt = np.pi / 2 - theta

        return az, alt

class EquidistantCamera(Camera):
    def __init__(self, f, resolution, pixel_ar = 1., principal_point = (0.5, 0.5)):
        self.f = f
        super(EquidistantCamera, self).__init__(resolution, pixel_ar, principal_point)

    def project(self, theta):
        return theta * self.f

    def unproject(self, r):
        return r / self.f

class RectilinearCamera(Camera):
    def __init__(self, f, resolution, pixel_ar = 1., principal_point = (0.5, 0.5)):
        self.f = f
        super(RectilinearCamera, self).__init__(resolution, pixel_ar, principal_point)

    def project(self, theta):
        result = np.tan(theta) * self.f
        result[theta > np.pi / 2.] = 1.e9
        return result

    def unproject(self, r):
        return np.arctan(r / self.f)

class OrthographicCamera(Camera):
    def __init__(self, f, resolution, pixel_ar = 1., principal_point = (0.5, 0.5)):
        self.f = f
        super(OrthographicCamera, self).__init__(resolution, pixel_ar, principal_point)

    def project(self, theta):
        result = np.sin(theta) * self.f
        result[theta > np.pi / 2.] = 1.e9
        return result

    def unproject(self, r):
        return np.arcsin(r / self.f)

class EquisolidAngleCamera(Camera):
    def __init__(self, f, resolution, pixel_ar = 1., principal_point = (0.5, 0.5)):
        self.f = f
        super(EquisolidAngleCamera, self).__init__(resolution, pixel_ar, principal_point)

    def project(self, theta):
        result = np.sin(theta / 2.) * 2. * self.f
        result[theta > np.pi] = 1.e9
        return result

    def unproject(self, r):
        return 2. * np.arcsin(r / (2. * self.f))

class StereographicCamera(Camera):
    def __init__(self, f, resolution, pixel_ar = 1., principal_point = (0.5, 0.5)):
        self.f = f
        super(StereographicCamera, self).__init__(resolution, pixel_ar, principal_point)

    def project(self, theta):
        result = np.tan(theta / 2.) * 2. * self.f
        result[theta > np.pi] = 1.e9
        return result

    def unproject(self, r):
        return 2. * np.arctan(r / (2. * self.f))

class CubicCamera(Camera):
    def __init__(self, k1, k2, resolution, pixel_ar = 1., principal_point = (0.5, 0.5)):
        self.k1 = k1
        self.k2 = k2
        super(CubicCamera, self).__init__(resolution, pixel_ar, principal_point)

    def project(self, theta):
        return self.k1 * theta + self.k2 * theta ** 3

    def unproject(self, r):
        # 0 = a * theta ^ 3 + c * theta + d
        a = self.k2
        c = self.k1
        d = -r

        if a == 0:
            theta = r/c
            return theta

        delta_0 = -3. * a * c
        delta_1 = 27. * a * a * d

        if a < 0:
            C = ((delta_1 + np.sqrt(np.complex64(delta_1 ** 2 - 4. * delta_0 ** 3))) / 2.) ** (1. / 3.)

            theta = np.real(-1. / (3. * a) * (C + delta_0 / (C)))
        else:
            C = ((delta_1 + np.sqrt(delta_1 ** 2 - 4. * delta_0 ** 3)) / 2) ** (1. / 3.)

            theta = -1. / (3. * a) * (C + delta_0 / (C))


        return theta


# Star Catalog
# ------------

class StarCatalog:
    def __init__(self, filename='hip_main.dat'):
        if filename is not None:
            self.read(filename)
            self.preprocess()

    def preprocess(self):
        filter_index = np.logical_not(np.logical_or(np.isnan(self.catalog['RAdeg']), np.isnan(self.catalog['Vmag'])))

        self.catalog = self.catalog[filter_index]

        self.star_vectors = angles_to_vector(np.deg2rad(self.catalog['RAdeg']), np.deg2rad(self.catalog['DEdeg']))

        self.magnitudes = self.catalog['Vmag'].values

        # alternative magnitude
        #VT = self.catalog['VTmag']
        #BT = self.catalog['BTmag']
        ## http://www.aerith.net/astro/color_conversion.html
        ## http://ads.nao.ac.jp/cgi-bin/nph-bib_query?bibcode=2002AJ....124.1670M&db_key=AST&high=3d1846678a19297
        #self.magnitudes = VT + 0.00097 - 0.1334 * (BT - VT) + 0.05486 * (BT - VT) ** 2 - 0.01998 * (BT - VT) ** 3

        # randomly fill out missing magnitudes
        #magnitudes = self.catalog['Vmag'].values
        #nans = np.isnan(magnitudes)
        #magnitudes[nans] = np.random.choice(magnitudes[~nans], np.sum(nans))
        #self.catalog['Vmag'] = magnitudes

    def lookup_indices(self, indices):
        result = indices.copy()
        result[indices >= 0] = self.catalog['HIP'].iloc[indices[indices >= 0]]
        return result

    def read(self, filename='hip_main.dat'):
        """Loads the Hipparchos star catalog."""

        columns = [
            "Catalog",
            "HIP",
            "Proxy",
            "RAhms",
            "DEdms",
            "Vmag",
            "VarFlag",
            "r_Vmag",
            "RAdeg",
            "DEdeg",
            "AstroRef",
            "Plx",
            "pmRA",
            "pmDE",
            "e_RAdeg",
            "e_DEdeg",
            "e_Plx",
            "e_pmRA",
            "e_pmDE",
            "DE:RA",
            "Plx:RA",
            "Plx:DE",
            "pmRA:RA",
            "pmRA:DE",
            "pmRA:Plx",
            "pmDE:RA",
            "pmDE:DE",
            "pmDE:Plx",
            "pmDE:pmRA",
            "F1",
            "F2",
            "---",
            "BTmag",
            "e_BTmag",
            "VTmag",
            "e_VTmag",
            "m_BTmag",
            "B-V",
            "e_B-V",
            "r_B-V",
            "V-I",
            "e_V-I",
            "r_V-I",
            "CombMag",
            "Hpmag",
            "e_Hpmag",
            "Hpscat",
            "o_Hpmag",
            "m_Hpmag",
            "Hpmax",
            "HPmin",
            "Period",
            "HvarType",
            "moreVar",
            "morePhoto",
            "CCDM",
            "n_CCDM",
            "Nsys",
            "Ncomp",
            "MultFlag",
            "Source",
            "Qual",
            "m_HIP",
            "theta",
            "rho",
            "e_rho",
            "dHp",
            "e_dHp",
            "Survey",
            "Chart",
            "Notes",
            "HD",
            "BD",
            "CoD",
            "CPD",
            "(V-I)red",
            "SpType",
            "r_SpType",
        ]

        self.catalog = pd.read_csv(filename, sep='|', names=columns, skipinitialspace=True)

# Star Detector
# -----

class StarDetector:
    def __init__(self, sigma_psf, t_exp, aperture, base_flux):
        self.sigma_psf = sigma_psf
        self.t_exp = t_exp
        self.aperture = aperture
        self.base_flux = base_flux
        
    @staticmethod
    def norm_gaussian(sigma):
        return math.erf((2*sigma)**-0.5)**2/4.
        
    def compute_flux(self, magnitude, add_noise=True):
        flux = self.base_flux * (10 ** (-magnitude / 2.5))

        if add_noise:
            flux = flux + np.random.normal(0, IMAGE_VARIANCE, len(flux))
            
            flux_per_photon = self.base_flux* self.t_exp * self.aperture ** 2 * np.pi/BASE_PHOTONS
            flux += flux_per_photon*np.random.normal(0, self.norm_gaussian(self.sigma_psf)*np.sqrt(flux/flux_per_photon), len(flux))

        return flux

    def compute_magnitude(self, flux):
        return -2.5 * np.log10(flux / self.base_flux)

    def compute_magnitude_threshold(self):
        threshold = THRESH_FACTOR*IMAGE_VARIANCE
        return self.compute_magnitude(threshold)

    def add_noise(self, magnitude):
        return self.compute_magnitude(self.compute_flux(magnitude))

# Scene
# -----

class Scene:
    def __init__(self, catalog, camera, detector, gaussian_noise_sigma=None, quantization_noise=None, magnitude_gaussian=None):
        self.catalog = catalog
        self.camera = camera
        self.detector = detector
        self.gaussian_noise_sigma = gaussian_noise_sigma
        self.quantization_noise = quantization_noise
        self.magnitude_gaussian = magnitude_gaussian
        self.orientation = None
        self.pos = None
        self.ids = None
        self.magnitude_threshold = detector.compute_magnitude_threshold()
    def compute(self, orientation=None):
        """Generates a scene for the star tracker.
        If not orientation is given a random one is generated.
        Gaussian noise is applied to star positions if enabled."""

        res_x, res_y = self.camera.resolution

        if orientation is None:
            orientation = random_matrix()

        self.orientation = orientation

        star_ids = np.arange(len(self.catalog.star_vectors))
        pos = np.dot(self.catalog.star_vectors, orientation.transpose())


        # instead vector noise
        if self.gaussian_noise_sigma:
            pos = add_vector_noise(pos, self.gaussian_noise_sigma)

        az, alt = vector_to_angles(pos)

        scene = self.camera.from_angles(az, alt)

        if self.quantization_noise:
            scene = np.int32(scene)

        selection = np.logical_and(np.logical_and(scene[:, 0] >= 0, scene[:, 0] < res_y), np.logical_and(scene[:, 1] >= 0, scene[:, 1] < res_x))

        scene = scene[selection]
        scene_ids = star_ids[selection]

        self.pos = scene
        self.magnitudes = self.catalog.magnitudes[scene_ids]
        self.ids = self.catalog.lookup_indices(scene_ids)

    def add_false_stars(self, false_stars):
        """Adds randomly generated false stars to a scene."""

        if isinstance(false_stars, int):
            res_x, res_y = self.camera.resolution

            false_star_pos = np.random.rand(false_stars, 2) * (res_y, res_x)
        else:
            false_star_pos = false_stars
            false_stars = len(false_stars)

        if self.quantization_noise:
            false_star_pos = np.int32(false_star_pos)

        self.pos = np.concatenate([self.pos, false_star_pos])
        self.magnitudes = np.concatenate([self.magnitudes, np.random.choice(self.catalog.magnitudes[self.catalog.magnitudes < self.magnitude_threshold], size=false_stars)])
        self.ids = np.concatenate([self.ids, -np.ones(false_stars, np.int32)])

    def copy(self, camera=None, orientation=None, copy_false_stars=True):

        false_stars = self.pos[self.ids == -1] if copy_false_stars else []

        if camera is None:
            camera = self.camera
        else:
            false_stars = camera.from_angles(*self.camera.to_angles(false_stars))

        if orientation is None:
            orientation = self.orientation

        scene = Scene(self.catalog, camera, self.detector, self.gaussian_noise_sigma, self.quantization_noise, self.magnitude_gaussian)
        scene.compute(orientation)

        scene.add_false_stars(false_stars)

        scene.scramble()

        scene.add_magnitude_noise(self.magnitude_gaussian)

        scene.filter_magnitudes()

        return scene

    def scramble(self):
        """Scrambles the order of stars in a scene."""

        scramble_index = np.random.permutation(range(len(self.ids)))

        self.pos = self.pos[scramble_index, ...]
        self.magnitudes = self.magnitudes[scramble_index]
        self.ids = self.ids[scramble_index]

    def add_magnitude_noise(self, gaussian=None):
        #if catalog:
        #    self.magnitudes[self.scene_ids >= 0] += np.random.normal(0, self.catalog.catalog['e_VTmag'])

        self.magnitudes = self.detector.add_noise(self.magnitudes)

        if gaussian is not None:
            self.magnitudes += np.random.normal(0, gaussian, size=len(self.magnitudes))

    def filter_magnitudes(self):
        filter_index = self.magnitudes <= self.magnitude_threshold
        self.pos = self.pos[filter_index]
        self.magnitudes = self.magnitudes[filter_index]
        self.ids = self.ids[filter_index]

    @staticmethod
    def random(catalog, camera, detector, min_true, max_true, min_false, max_false, min_stars=1, max_tries = 1000, gaussian_noise_sigma=None, quantization_noise=None, magnitude_gaussian=None):
        scene = Scene(catalog, camera, detector, gaussian_noise_sigma, quantization_noise, magnitude_gaussian)

        num_stars = 0
        tries = 0

        ok = False

        while not ok and tries < max_tries:
            scene.compute()

            scene.add_false_stars(np.random.randint(min_false, max_false + 1))

            scene.scramble()

            scene.add_magnitude_noise(magnitude_gaussian)

            scene.filter_magnitudes()

            num_stars = np.sum(scene.ids >= 0)

            ok = num_stars >= min_true and num_stars <= max_true and len(scene.ids) > min_stars

            tries += 1

        if tries == max_tries:
            return None

        return scene

#only do this part if we were run as a python script
if __name__ == '__main__':
	if len(sys.argv)!=4:
		sys.stdout.write("Usage: python simulator.py calibration.txt input.csv results.csv\n")
		exit()

	exec(open(sys.argv[1]).read())

	res_x = float(IMG_X) # pixels
	res_y = float(IMG_Y) # pixels

	# normalized focal length
	f = 0.5 / np.tan(np.deg2rad((IMG_X*PIXSCALE/3600)) / 2.)
	# pixel aspect ratio
	pixel_ar = 1.

	# normalized principal point
	ppx = 0.5
	ppy = 0.5

	gaussian_noise_sigma = (np.pi*PIXSCALE)*np.sqrt(POS_VARIANCE)/(180*3600.0) # rad

	cam = 0

	sigma_psf = 0.5*float(DOUBLE_STAR_PX) # pixel


	t_exp = float(EXPOSURE_TIME) # s
	aperture = float(APERTURE) # mm
	base_flux = float(BASE_FLUX)

	magnitude_gaussian = 0.02 # mag

	min_true = 0
	max_true = 1000
	min_false = 0
	max_false = MAX_FALSE_STARS


	catalog = StarCatalog()
	cameras = [
		RectilinearCamera,
		EquidistantCamera,
		EquisolidAngleCamera,
		StereographicCamera,
		OrthographicCamera,
	]

	camera = cameras[cam](f, (res_x, res_y), pixel_ar, (ppx, ppy))
	detector = StarDetector(sigma_psf, t_exp, aperture, base_flux)
	num_scenes = 100
	inputs = []
	outputs = []

	for i in range(num_scenes):
		scene = Scene.random(catalog, camera, detector, min_true, max_true, min_false, max_false, gaussian_noise_sigma=gaussian_noise_sigma, magnitude_gaussian=magnitude_gaussian)
		inputs.append(np.hstack((scene.pos[::, ::-1], scene.magnitudes.reshape(-1, 1))).flatten())
		outputs.append(scene.ids)

	def write_csv(filename, lines):
		with open(filename, 'w') as f:
			for line in lines:
				f.write(','.join(str(value) for value in line) + '\n')

	write_csv(sys.argv[2], inputs)
	write_csv(sys.argv[3], outputs)
