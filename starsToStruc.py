import ctypes
import socket 

from io import StringIO, BytesIO

class stars_to_struct():
    fields = [("stars", ctypes.c_int),
              ("declination", ctypes.c_float),
              ("right_ascension", ctypes.c_float),
              ("body_x", ctypes.c_float*"stars"),
              ("body_y", ctypes.c_float*"stars"),
              ("body_z", ctypes.c_float*"stars"),
              ("ref_x", ctypes.c_float*"stars"),
              ("ref_y", ctypes.c_float*"stars"),
              ("ref_z", ctypes.c_float*"stars"),
              ("uncertainty", ctypes.c_float*"stars")]
    
    def _init_(self):
        super(stars_to_struct, self)._init_()
        self.stars = 0
        self.declination = 0.0
        self.right_ascension = 0.0
        self.body_x = [0.0]*self.stars
        self.body_y = [0.0]*self.stars
        self.body_z = [0.0]*self.stars
        self.ref_x = [0.0]*self.stars
        self.ref_y = [0.0]*self.stars
        self.ref_z = [0.0]*self.stars
        self.uncertainty = [0.0]*self.stars
    
    def load_image(self, stars, declination, right_ascension, body_x, body_y, body_z, ref_x, ref_y, ref_z, uncertainty): 
        self.stars = stars
        self.declination = declination
        self.right_ascension = right_ascension
        self.body_x = body_x
        self.body_y = body_y
        self.body_z = body_z
        self.ref_x = ref_x
        self.ref_y = ref_y
        self.ref_z = ref_z
        self.uncertainty = uncertainty
