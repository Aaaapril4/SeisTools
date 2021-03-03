import obspy
from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn.mass_downloader import RectangularDomain,Restrictions,MassDownloader

import numpy as np
import configparser
import os
import multiprocessing as mp
from tqdm import tqdm
import pyasdf

client = Client("IRIS")



def get_station():
    '''
    Get the info of all available networks and stations satisfying requirements
    Return:
        station inventory
    '''

    Netinv = client.get_stations(
        network = para['Station Info'].get('network', '*'),
        station = para['Station Info'].get('station', '*'),
        channel = para['Station Info'].get('channel', '*'),
        starttime = para['Station Info'].get('starttime', '19900101T00:00:00'),
        endtime = para['Station Info'].get('endtime', '20201231T23:59:59'),
        maxlatitude = para['Map Info'].getfloat('maxlatitude'),
        minlatitude = para['Map Info'].getfloat('minlatitude'),
        maxlongitude = para['Map Info'].getfloat('maxlongitude'),
        minlongitude = para['Map Info'].getfloat('minlongitude')
    )

    return Netinv



def get_event():
    start = UTCDateTime('19940101T000000')
    end = UTCDateTime('20201231T000000')
    file = '/mnt/home/jieyaqi/Documents/wavecat.xml'

    client = Client("IRIS")
    cat = client.get_events(starttime=start,
                endtime=end,
                latitude=-4.5,
                longitude=35,
                minradius=30,
                maxradius=90,
                minmagnitude=6,
                filename=file)

    cat = obspy.read_events(file)

    for event in cat:
        extra = {'downloaded': {'value': 'False',
                'namespace':'http://test.org/xmlns/1.0'}}
        event.extra = extra
    cat.write(file,format='QUAKEML')
    
    return



def _get_nettime(nw):
    '''
    Calculate the time range of each network for downloading
    Return:
        starttime, endtime
    '''

    nwbegin = min([sta.start_date for sta in nw.stations])
    if None in [sta.end_date for sta in nw.stations]:
        nwend = UTCDateTime.now()
    else:
        nwend = max([sta.end_date for sta in nw.stations])

    return max(nwbegin, para['Station Info'].get('starttime')), min(nwend, para['Station Info'].get('endtime'))



def get_mseed_storage(network, station, location, channel, starttime, endtime):
    '''
    Check if mseed file exists, get name if not
    '''

    mseedname = os.path.join(para['DEFAULT'].get('projdir'), "data/waveform/%s.%s.%s.%s_%s_%s.mseed" % (network, station, location, channel, starttime, endtime))
    if os.path.exists(mseedname):
        True
    else:
        return mseedname



def get_station_storage(network, station):
    '''
    Check if mseed file exists, get name if not
    '''

    stationname = os.path.join(para['DEFAULT'].get('projdir'), "data/station/%s.%s.mseed" % (network, station))
    if os.path.exists(stationname):
        True
    else:
        return stationname



def _download_cont(nwinv):
    '''
    Download data for each network
    Parameter:
        Network inventory
    Return:
        None
    '''

    domain = RectangularDomain(
        minlatitude = para['Map Info'].getfloat('minlatitude'),
        maxlatitude = para['Map Info'].getfloat('maxlatitude'),
        minlongitude = para['Map Info'].getfloat('minlongitude'),
        maxlongitude = para['Map Info'].getfloat('maxlongitude')) #para

    starttime, endtime = _get_nettime(nwinv)

    restrictions = Restrictions(
        starttime = starttime,
        endtime = endtime,
        chunklength_in_sec = 86400,
        network = nwinv.code,
        station = para['Station Info'].get('station'),
        channel_priorities = ["L??", "B??", "H??"],
        reject_channels_with_gaps = True,
        minimum_length = 0.1,
        minimum_interstation_distance_in_m = 100.0)

    mdl = MassDownloader(providers=["IRIS"])
    mdl.download(domain, restrictions, mseed_storage = get_mseed_storage, stationxml_storage = get_station_storage)

    return



def saveasdf(filename, mseeddir, stationdir, eventdir = None):
    '''
    Create asdf file for data
    Parameters:
        filename of asdf, path waveform directory, station inventory or path of stationxml dictionary, event catalog or path of eventxml
    Return:
        None
    '''

    ds = pyasdf.ASDFDataSet(f'{para["DEFAULT"].get("projdir")}/{filename}.h5', compression="gzip-3")

    for mseed in os.listdir(mseeddir):
        if ".mseed" not in mseed:
            pass
        else:
            ds.add_waveforms(mseed, tag = 'raw')
    
    if isinstance(stationdir, obspy.core.inventory):
        for nwinv in stationdir:
            for stainv in stationdir:
                ds.add_stationxml(stationdir)
    elif isinstance(stationdir, str):
        for xml in os.listdir(stationdir):
            if ".xml" not in xml:
                pass
            else:
                ds.add_stationxml(xml)
    
    if isinstance(eventdir, obspy.core.event.catalog):
        ds.add_quakeml(eventdir)
    elif isinstance(eventdir, str):
        for xml in os.listdir(eventdir):
            if ".xml" not in xml:
                pass
            else:
                ds.add_quakeml(xml)
    
    return 



def download_cont():
    '''
    Download continuous waveform for all networks and save as asdf
    Return:
        None
    '''

    Netinv = get_station()
    with mp.pool(para['DEFAULT'].get('ncpu')) as p:
        list(tqdm(p.imap(_download_cont, [nw for nw in Netinv]),total=len(Netinv)))
    
    saveasdf('amnoise', f'{para["DEFAULT"].get("projdir")}/data/waveform', Netinv)

    return
    

if __name__ == '__main__':
    para = configparser.ConfigParser()
    para.read('para.ini')
    download_cont()
