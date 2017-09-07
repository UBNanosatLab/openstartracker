from simulator import *
#from example.py
res_x = 1392 # pixels
res_y = 1080 # pixels

# normalized focal length
f = 0.5 / np.tan(np.deg2rad(6) / 2)
# pixel aspect ratio
pixel_ar = 1

# normalized principal point
ppx = 0.5
ppy = 0.5

gaussian_noise_sigma = 20e-6 # rad

cam = 0
# magnitude parameters

A_pixel = 525 # photonelectrons/s mm
sigma_pixel = 525 # photonelectrons/s mm

sigma_psf = 0.5 # pixel
t_exp = 0.05 # s
aperture = 60.7 # mm

base_photons = 19100 # photoelectrons per mm^2 and second of a magnitude 0 G2 star

magnitude_gaussian = 0.01 # mag
# star count

min_true = 0
max_true = 100
min_false = 0
max_false = 5


catalog = StarCatalog()
cameras = [
    RectilinearCamera,
    EquidistantCamera,
    EquisolidAngleCamera,
    StereographicCamera,
    OrthographicCamera,
]

camera = cameras[cam](f, (res_x, res_y), pixel_ar, (ppx, ppy))
detector = StarDetector(A_pixel, sigma_pixel, sigma_psf, t_exp, aperture, base_photons)
num_scenes = 300
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

write_csv('input.csv', inputs)
write_csv('result.csv', outputs)

#_ = scene.render(False)
