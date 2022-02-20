import obspy
from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn.mass_downloader import Restrictions,MassDownloader

import numpy as np
import os
import multiprocessing as mp
from tqdm import tqdm
import pyasdf

from .config import domain, para

client = Client("IRIS")



def get_station():
    '''
    Get the info of all available networks and stations satisfying requirements
    Return:
        station inventory
    '''

    Netinv = client.get_stations(
        network = para["Station Info"].get("network", "*"),
        station = para["Station Info"].get("station", "*"),
        channel = para["Station Info"].get("channelpri", "*"),
        starttime = para["Station Info"].get("starttime", "19900101T00:00:00"),
        endtime = para["Station Info"].get("endtime", UTCDateTime.now()),
        maxlatitude = para["Map Info"].getfloat("maxlatitude"),
        minlatitude = para["Map Info"].getfloat("minlatitude"),
        maxlongitude = para["Map Info"].getfloat("maxlongitude"),
        minlongitude = para["Map Info"].getfloat("minlongitude")
    )

    filter_net = para["Station Info"].get("network_filter").split(",")
    for net in filter_net:
        Netinv = Netinv.remove(network=net)

    if para["Save Data"].getboolean("stationinfo"):
        Netinv.write(f'{para["DEFAULT"].get("projdir")}/station.txt', format="STATIONTXT", level='station')

    return Netinv



def get_event_radius(minradius, maxradius, minmag):
    '''
    Get all the events in time range
    Return:
        Event catalog
    '''
    
    maxlatitude = para["Map Info"].getfloat("maxlatitude")
    minlatitude = para["Map Info"].getfloat("minlatitude")
    maxlongitude = para["Map Info"].getfloat("maxlongitude")
    minlongitude = para["Map Info"].getfloat("minlongitude")
    
    cat = client.get_events(
        starttime=para["Station Info"].get("starttime", "19900101T00:00:00"),
        endtime=para["Station Info"].get("endtime", UTCDateTime.now()),
        latitude=(maxlatitude+minlatitude) / 2,
        longitude=(maxlongitude+minlongitude) / 2,
        minradius=minradius,
        maxradius=maxradius,
        minmagnitude=minmag)

    for event in cat:
        extra = {"downloaded": {"value": "False",
                "namespace":"http://test.org/xmlns/1.0"}}
        event.extra = extra
    
    if para["Save Data"].getboolean("eventcatlog"):
        cat.write(f'{para["DEFAULT"].get("projdir")}/evcatalog.xml', format="QUAKEML")
    
    return cat



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

    return max(nwbegin, UTCDateTime(para["Station Info"].get("starttime"))), min(nwend, UTCDateTime(para["Station Info"].get("endtime")))



def get_nettime(Netinv):
    nettime = []
    for nw in Netinv:
        starttime, endtime = _get_nettime(nw)
        nettime.append((nw.code, starttime, endtime))
    return nettime



def _get_downloadlist(Netinv):
    '''
    Split the network time into several slots for parallal downloading
    Return:
        list of tunnel (networkname, begintime, endtime)
    '''

    downloadlist = []
    for nw in Netinv:
        starttime, endtime = _get_nettime(nw)
        stationday = len(nw) * (endtime - starttime) / 60 / 60 / 24
        
        if stationday <= 50 * 365:
            downloadlist.append((nw.code, starttime, endtime))
        elif stationday > 50 * 365:
            num = np.ceil(stationday / 50 / 365)
            interval = np.ceil((endtime - starttime) / num)
            timeslot = np.arange(starttime, endtime+1, interval)
            if endtime not in timeslot:
                timeslot = np.append(timeslot, endtime)

            for i in range(len(timeslot)-1):
                downloadlist.append((nw.code, timeslot[i], timeslot[i+1]))
    
    return downloadlist



def _download_cont(nw, starttime, endtime):
    '''
    Download data for each network
    Return:
        None
    '''
    print(f'======Download data for network {nw}: {starttime} - {endtime}======')
    restrictions = Restrictions(
        starttime = starttime,
        endtime = endtime,
        chunklength_in_sec = para["Station Info"].getint("chunksize", 1) * 24 * 60 * 60,
        network = nw,
        station = para["Station Info"].get("station", "*"),
        channel_priorities = para["Station Info"].get("channelpri", "*").split(","),
        reject_channels_with_gaps = False,
        minimum_length = 0.0,
        minimum_interstation_distance_in_m = 100.0)

    mdl = MassDownloader(providers=["IRIS"])
    mdl.download(
        domain, 
        restrictions, 
        mseed_storage = para["DEFAULT"].get("projdir") + "/waveform/{station}/{network}.{station}.{location}.{channel}__{starttime}__{endtime}.mseed", 
        stationxml_storage = para["DEFAULT"].get("projdir") + "/station/{network}.{station}.xml")

    return



def download_cont():
    '''
    Download continuous waveform for all networks and save as asdf
    Return:
        None
    '''

    Netinv = get_station()
    downloadlist = _get_downloadlist(Netinv)
    if para["DEFAULT"].getint("ncpu") == 1:
        for i in range(len(downloadlist)):
            _download_cont(*downloadlist[i])
    else:
        with mp.Pool(para["DEFAULT"].getint("ncpu")) as p:
            list(tqdm(p.starmap(_download_cont, downloadlist),total=len(downloadlist)))
    
    if para["Save Data"].getboolean("asdffile"):
        saveasdf("amnoise", f'{para["DEFAULT"].get("projdir")}/waveform', Netinv)

    return



def _download_event(eventtime, starttime, endtime):
    '''
    Download each event for all stations
    Return:
        None
    '''

    restrictions = Restrictions(
                        starttime = eventtime - starttime * 60,
                        endtime = eventtime + endtime * 60,
                        network = para["Station Info"].get("network", "*"),
                        station = para["Station Info"].get("station", "*"),
                        reject_channels_with_gaps = False,
                        minimum_length = 0.0,
                        channel_priorities = para["Station Info"].get("channelpri", "*").split(","))

    mdl = MassDownloader(providers=["IRIS"])
    mdl.download(
        domain, 
        restrictions, 
        mseed_storage = f'{para["DEFAULT"].get("projdir")}/waveform', 
        stationxml_storage = f'{para["DEFAULT"].get("projdir")}/station')
    return



def _download_each_event(event):
    '''
    Download one event teleseismic P wave for all stations
    Output:
        None
    '''

    if event.extra['downloaded']['value'] == "False":
            eventtime = event.origins[0].time
            try:
                _download_event(eventtime, 5, 60)
            except:
                return
            event.extra['downloaded']['value'] = "True"
    
    return event



def download_event():
    '''
    Download events for teleseismic P wave
    Output:
        None
    '''

    cat = get_event_radius(30,90,6)
    with mp.Pool(para["DEFAULT"].getint("ncpu")) as p:
        cat_list = list(tqdm(p.imap(_download_each_event, [event for event in cat]),total=len(cat)))
    
    cat_new = cat.copy()
    cat_new.clear()
    for event in cat_list:
        cat_new.append(event)
    if para["Save Data"].getboolean("asdffile"):
        saveasdf("teleseismic", f'{para["DEFAULT"].get("projdir")}/waveform', f'{para["DEFAULT"].get("projdir")}/station', cat_new)

    return



def saveasdf(filename, mseeddir, stationdir, eventdir = None):
    '''
    Create asdf file for data
    Return:
        None
    '''

    ds = pyasdf.ASDFDataSet(f'{para["DEFAULT"].get("projdir")}/{filename}.h5', compression="gzip-3")

    for mseed in os.listdir(mseeddir):
        if ".mseed" not in mseed:
            pass
        else:
            ds.add_waveforms(f'{mseeddir}/{mseed}', tag = "raw")
    
    if isinstance(stationdir, obspy.core.inventory):
        for nwinv in stationdir:
            for stainv in nwinv:
                ds.add_stationxml(f'{stationdir}/{stainv}')
    elif isinstance(stationdir, str):
        for xml in os.listdir(stationdir):
            if ".xml" not in xml:
                pass
            else:
                ds.add_stationxml(f'{stationdir}/{xml}')
    
    if isinstance(eventdir, obspy.core.event.catalog):
        ds.add_quakeml(eventdir)
    elif isinstance(eventdir, str):
        for xml in os.listdir(eventdir):
            if ".xml" not in xml:
                pass
            else:
                ds.add_quakeml(f'{eventdir}/{xml}')
    
    return 



if __name__ == '__main__':
    # download_event()
    
    download_cont()
