#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
A simple script to plot the visibility of an asteroid
from a selected observatory in a given time. 
It is based on this astropy tutorial:
https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_obs-planning.html

The get_indexfile function is from this project:
https://github.com/rszakats/astname

Usage: 
bokeh serve --show asteroid_visibility.py
"""
import getpass
import os
import platform
import warnings
from datetime import datetime

import astropy.units as u
import numpy as np
import pytz
from astropy.coordinates import (AltAz, EarthLocation, SkyCoord, get_moon,
                                 get_sun)
from astropy.table import Table
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (AutocompleteInput, ColorBar, ColumnDataSource,
                          DatePicker, HoverTool, LinearColorMapper, Select,
                          Span, TextInput)
from bokeh.palettes import Viridis256, d3
from bokeh.plotting import Figure, figure
from pytz import timezone
from timezonefinder import TimezoneFinder

from astname import get_indexfile

warnings.filterwarnings("ignore")

def download_asteroid_data():
    """
    Downloads the NASA JPL/Horizons index file, DASTCOM.IDX and saves it to
    the cache folder.

    Parameters
    ----------
    None

    Returns
    -------
    inputfile: str
        The name of the input file.
    """
    if platform.system() == 'Linux':
        cachedir = '/home/'+str(getpass.getuser())+'/.cache/'
    elif platform.system() == 'Windows':
        cachedir = "c:\\Users\\"+str(getpass.getuser())+"\\.cache\\"
    if os.path.exists(cachedir) is False:
        os.mkdir(cachedir)
    datafile = 'DASTCOM.IDX'
    url = 'ftp://ssd.jpl.nasa.gov/pub/xfr/DASTCOM.IDX'
    # If the index file is missing, the script will try to download it.
    if os.path.exists(cachedir+datafile) is False:
        print("No index file was found! Downloading...")
        get_indexfile(url, cachedir+datafile)
    inputfile = f"{cachedir}{datafile}"
    return inputfile

#-----------------------------------------------------------------------------------------------------------
def calc_visibility(name, otime, observing_location, observatoryid):
    """
    Calculates the object's azimuth and altitude from the given
    location at the specified time. Updates the sources for the 
    Bokeh plots

    Parameters
    ----------
    name: str
        Name(number) of the target. 

    otime: str
        Observing time. Format: YYYY-MM-DD HH:mm:ss
    
    observing_location: EarthLocation object
        Location of the observatory.
    
    observatoryid: str
        Observatory ID of the observing_location. 
        See the available IDs here:
        https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html
    
    Returns
    -------
    None
    """
    otime = str(otime).split()[0]
    utc = pytz.utc
    loc = TimezoneFinder()
    timez = loc.timezone_at(lng=observing_location.lon.value, lat=observing_location.lat.value)
    loc_timez = timezone(timez)
    y = int(str(otime).split("-")[0])
    m = int(str(otime).split("-")[1])
    d = int(str(otime).split("-")[2])
    loc_otime = datetime(y,m,d,0,0,0,tzinfo=utc)
    loc_dt = loc_otime.astimezone(loc_timez)
    utcoffset = loc_dt.utcoffset().total_seconds()/3600
    otime = Time(Time(otime), format='iso', scale='utc',
         location=(observing_location.lon, observing_location.lat))
    
    jd = Time(otime - 2*utcoffset*u.hour).jd
    print(f"Processing: {name} ")
    obj = Horizons(id=f"{name}:", location=observatoryid, epochs=jd)
    target = obj.ephemerides()
    center = SkyCoord(target['RA'][0]*u.deg, target['DEC'][0]*u.deg,
                        frame='icrs')
    center = center.transform_to(AltAz(obstime=otime,
                                    location=observing_location))

    midnight = otime - 2*utcoffset*u.hour
    delta_midnight = np.linspace(-10, 10, 1000)*u.hour
    times = midnight + delta_midnight + utcoffset*u.hour
    frame2 = AltAz(obstime=times, location=observing_location)
    sunaltazs = get_sun(times).transform_to(frame2)
    moon = get_moon(times)
    moonaltazs = moon.transform_to(frame2)

    center2 = center.transform_to(frame2)
    sep = moon.separation(center)
    print(f"Moon separation for {name}: {np.mean(sep)}")

    label1 = []
    label2 = []
    label3 = []
    color3 = []
    for _ in range(len(delta_midnight)):
        label1.append("Sun")
        label2.append("Moon")
        label3.append(name)
        color3.append("green")
    
    source1.data = {'xs':delta_midnight.value,
                   'ys': sunaltazs.alt.value,
                   'labels': label1,
                   'time': times.value}
    source2.data = {'xs':delta_midnight.value,
                   'ys': moonaltazs.alt.value,
                   'labels': label2,
                   'time': times.value}
    source3.data = {'xs':delta_midnight.value.tolist(),
                   'ys': center2.alt.value.tolist(),
                   'labels': label3,
                   'color': color3,
                   'az': center2.az.value.tolist(),
                   'time': times.value.tolist()}

    interval4 = source1.data['xs'][sunaltazs.alt < -0*u.deg]
    center4 = interval4[int(len(interval4)/2)]
    interval5 = source1.data['xs'][sunaltazs.alt < -12*u.deg]
    center5 = interval5[int(len(interval5)/2)]
    source4.data = {'center': [center4], 'width': [max(source1.data['xs'][sunaltazs.alt < -0*u.deg])-min(source1.data['xs'][sunaltazs.alt < -0*u.deg])]}
    source5.data = {'center': [center5], 'width': [max(source1.data['xs'][sunaltazs.alt < -12*u.deg])-min(source1.data['xs'][sunaltazs.alt < -12*u.deg])]}

#-----------------------------------------------------------------------------------------------------------
def callback(attr, old, new):
    """
    Callback function for the input field widgets.
    If nothing is set, uses default values for input.
    Calls the calc_visibility() function to update
    the plot sources.

    Parameters
    ----------
    attr: str
        String passed to the callback function.

    old: str
        Old value of the input field.
    
    new: str
        New value of the input field.
    
    Returns
    -------
    None
    """
    if sname.value == "":
        sname.value = '1 Ceres'
    if date.value == "":
        date.value = "2023-05-05"
    if autoc.value == "":
        autoc.value = "Piszkesteto Stn. (Konkoly)"
    sname.value = str(sname.value).split()[0]
    obsdata = get_obsdata(observatories, autoc.value)
    observing_location = EarthLocation(lat=obsdata['gclat'][0], lon=obsdata['long'][0],
                                   height=obsdata['height'][0]*u.m)
    calc_visibility(sname.value, date.value, observing_location, obsdata['code'][0])
    update_info()

#-----------------------------------------------------------------------------------------------------------
def search_callback2(attr, old, new):
    """
    Callback function for the object search field widget.
    It uses the data from the DASTCOM.IDX file to find
    the asteroid. It extracts all the matches in case 
    of a string search.

    Calls the update_info() function to update
    the plot sources.

    Parameters
    ----------
    attr: str
        String passed to the callback function.

    old: str
        Old value of the input field.
    
    new: str
        New value of the input field.
    
    Returns
    -------
    None
    """
    try:
        targ = int(new)
        if isinstance(targ, int):
            if int(targ) > 2000000 and int(targ) < 50000000:
                result = fulldata[np.where(fulldata['naifid'] == f'{new}')]
                sname.options = [i +" "+ j  for i, j in zip(result['astnum'].tolist(), result['astid'].tolist())]
                source0.data = {'ssos': sname.options}
            else:
                result = fulldata[np.where(fulldata['astnum'] == f'{new}')]
                sname.options = [i +" "+ j  for i, j in zip(result['astnum'].tolist(), result['astid'].tolist())]
                source0.data = {'ssos': sname.options}
    except ValueError:
        if isinstance(new, str):
            result = fulldata[np.core.defchararray.find(np.char.lower(fulldata['astid'].value), f"{new}".lower())!=-1]
            sname.options = [i +" "+ j  for i, j in zip(result['astnum'].tolist(), result['astid'].tolist())]
            source0.data = {'ssos': sname.options}
    if len(sname.options) > 0:
        sname.value = sname.options[0]
        update_info()

#-----------------------------------------------------------------------------------------------------------
def update_info():
    """
    Function to update the source data for the plot
    with the selected asteroid's data.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    data = fulldata[np.where(fulldata['astnum'] == f'{sname.value}')]
    info = f"Number: {data['astnum'][0]}\nName/ID: {data['astid'][0]}\n"\
           f"NAIFID: {data['naifid'][0]}\nOther ID(s): {str(data['altid'][0]).replace(',', ';')}"
    infos = []
    labels = []
    for _ in range(len(source3.data['xs'])):
        labels.append(f"{data['astnum'][0]} {data['astid'][0]}")
        infos.append(info)
    source3.data['info'] = infos
    source3.data['labels'] = labels
    p.title.text = f"Visibility plot for {data['astnum'][0]} {data['astid'][0]} at {date.value} from {autoc.value}"

    
#-----------------------------------------------------------------------------------------------------------
def read_obsdata():
    """
    Function to read the observatory data file.

    Parameters
    ----------
    None

    Returns
    -------
    t: astropy.table object
        The observatory data in an astropy Table.
        <Table length=2408>
        name   dtype 
        ------ -------
        name   str48    # observatory name
        code    str3    # Observatory Code assigned by the Minor Planet Center
        long float64    # Longitude of the observatory (deg)
        gclat float64   # Geocentric latitude of the observatory (deg)
        gdlat float64   # Geodetic latitude of the observatory WGS84 (deg)
        height float64  # Altitude of the observatory (m)
        map   str65     # link to the location on openstreet map

    """
    # https://dc.zah.uni-heidelberg.de/obscode/q/query/form
    t = Table.read('observatories.csv', format='ascii.csv')
    return t

#-----------------------------------------------------------------------------------------------------------
def get_obsdata(observatories, observatory):
    """
    Gets the observatory data for the selected observatory.


    Parameters
    ----------
    observatories: astropy.table object
        The return value of the read_obsdata() function.

    observatory: str
        Value of the autocomplete input widget, i.e.
        the name of the selected osbervatory.


    Returns
    -------
    obs: astropy.table object
        The observatory data in an astropy Table.
        <Table length=1>
        name   dtype 
        ------ -------
        name   str48    # observatory name
        code    str3    # Observatory Code assigned by the Minor Planet Center
        long float64    # Longitude of the observatory (deg)
        gclat float64   # Geocentric latitude of the observatory (deg)
        gdlat float64   # Geodetic latitude of the observatory WGS84 (deg)
        height float64  # Altitude of the observatory (m)
        map   str65     # link to the location on openstreet map
    """

    obs = observatories[observatories['name'] == observatory]
    return obs

#-----------------------------------------------------------------------------------------------------------
# def conv(s):
#     return s.strip()

#-----------------------------------------------------------------------------------------------------------
def read_astdata(infile):
    """
    Function to read the DASTCOM.IDX file to
    an astropy.table object.

    Parameters
    ----------
    infile: str
        Full path to the DASTCOM.IDX file.

    Returns
    t: astropy.table
        The astropy.table object containing the 
        whole data from the input file.
        <Table length=1286232>
        name  dtype
        ------ -----
        astnum  str8  # Asteroid IAU number, see: https://sci.esa.int/web/home/-/30244-asteroid-numbers-and-names
        astid str28   # Asteroid name/designation
        naifid  str8  # NAIF Integer ID code, see: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
        altid str70   # Other designations

    """
    astnums = []
    astids = []
    naifids = []
    altids = []
    i = 0
    with open(infile, "r") as file:
        for line in file:
            data = line.split(",")
            astnums.append(data[0].split()[0])
            if len(data[0].split()[1:]) > 1:
                ids = ""
                for ii in data[0].split()[1:]:
                    ids += f"{ii} "
                astids.append(ids)
            else:
                astids.append(data[0].split()[1:][0])
            naifids.append(data[1])
            alt = ""
            for j in range(2, 8):
                try:
                    alt += str(line.split(',')[j])+", "
                except IndexError:
                    alt = alt
            alt = alt.replace(", ,", "")
            alt = alt.rstrip(" ,\n")
            altids.append(alt)
            i += 1
    t = Table([astnums, astids, naifids, altids],
            names=("astnum", "astid", "naifid", "altid"),meta={'name': 'first table'})
    return t


#-----------------------------------------------------------------------------------------------------------
astdat_file = download_asteroid_data()
fulldata = read_astdata(astdat_file)
astdat = [i +" "+ j  for i, j in zip(fulldata['astnum'].tolist(), fulldata['astid'].tolist())]

data = {'ssos': astdat}
source = ColumnDataSource(data)

user = getpass.getuser()
if os.path.exists("/tmp/"+str(user)) is False:
    os.mkdir("/tmp/"+str(user))

observatories = read_obsdata()

source0 = ColumnDataSource(data={'ssos':['1 Ceres']})
source1 = ColumnDataSource(data={'xs':[], 'ys': [], 'labels': []})
source2 = ColumnDataSource(data={'xs':[], 'ys': [], 'labels': []})
source3 = ColumnDataSource(data={'xs':[], 'ys': [], 'labels': ['1'],
                                 'color': [], 'az': [1.],
                                 'time': [datetime.now(), datetime.now()],
                                 'info': []})
source4 = ColumnDataSource(data={'center': [0], 'width': [0]})
source5 = ColumnDataSource(data={'center': [0], 'width': [0]})

date = DatePicker(title='Choose date', value=f"{datetime.now()}")
name = TextInput(title='Search for target')
sname = Select(title='Choose a target:', options=source0.data['ssos'], value=f"{source0.data['ssos'][0]}")
autoc = AutocompleteInput(title='Search Observatory:', completions=observatories['name'].tolist(),
                          case_sensitive=False, placeholder="Search for an observatory.")

date.on_change('value', callback)
name.on_change('value', search_callback2)
name.on_change('value', callback)
sname.on_change('value', callback)
autoc.on_change('value', callback)

p = figure(height=1000, width=1200, y_range=(0, 90) , x_range=(-10, 10))#, tools=[hover]) #  , x_range=(-12, 12)
p.rect(x='center', y=45, width='width', height=90, source=source5,
       fill_alpha=0.8, line_color=d3['Category20c'][20][-1], fill_color=d3['Category20c'][20][-1])
p.rect(x='center', y=45, width='width', height=90, source=source4,
       fill_alpha=0.4, line_color=d3['Category20c'][20][16], fill_color=d3['Category20c'][20][16])

color_mapper = LinearColorMapper(palette=Viridis256, low=0, high=360)
sun = p.line(x='xs', y='ys', source=source1, legend_label="Sun", color="orange", line_width=1)
moon = p.line(x='xs', y='ys', source=source2, legend_label="Moon", color=d3['Category10'][3][0], line_width=2, line_dash='dashed')
targ = p.circle(x='xs', y='ys', source=source3, legend_field="labels", color={'field': 'az', 'transform': color_mapper})
hline = Span(location=20, dimension='width', line_color='red', line_width=2, line_dash='dashed', line_alpha=0.5)
p.renderers.extend([hline])

bar = ColorBar(color_mapper=color_mapper, location=(0,0), title="Azimuth (deg)")
p.add_layout(bar, "right")

p.add_tools(HoverTool(renderers=[sun, moon, targ],line_policy='interp',
                  tooltips=[("object", "@labels"), ("UTC", "@time"), ("delta Time from midnight", "@xs"),
                            ("altitude", "@ys"), ("azimuth", "@az"),
                            ("Info:" ,"@info")],
                  point_policy='snap_to_data'))

p.xaxis.axis_label = 'Hours from Local Midnight'
p.yaxis.axis_label = 'Altitude (deg)'
p.legend.click_policy="hide"
p.title.text = f"Visibility plot for {str(source3.data['labels'][0])} at {date.value} from {autoc.value}"

top_row = row(date, name, sname, autoc)
middle_row = row(p)

curdoc().add_root(column(top_row, middle_row))
curdoc().theme = 'caliber'
curdoc().title = "Visibility plot"
