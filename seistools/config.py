from obspy.clients.fdsn.mass_downloader import RectangularDomain
import configparser
import sys

para = configparser.ConfigParser()
para.read('para.ini')

domain = RectangularDomain(
    minlatitude = para["Map Info"].getfloat("minlatitude"),
    maxlatitude = para["Map Info"].getfloat("maxlatitude"),
    minlongitude = para["Map Info"].getfloat("minlongitude"),
    maxlongitude = para["Map Info"].getfloat("maxlongitude"))