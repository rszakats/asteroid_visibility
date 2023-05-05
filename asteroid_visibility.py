#!/usr/bin/env python3.8

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astroquery.jplhorizons import Horizons
# import datetime
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
from astropy.table import Table
# import time
import getpass
import os
import warnings
from bokeh.plotting import Figure, figure
from bokeh.models import ColumnDataSource, WheelZoomTool, DatePicker, TextInput, ResetTool, SaveTool, PanTool, Select, AutocompleteInput, MultiLine, Select
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.palettes import d3, Viridis256
from bokeh.models import LinearColorMapper, HoverTool, Span
from bokeh.models import ColorBar
from timezonefinder import TimezoneFinder
import pytz
from pytz import timezone
from datetime import datetime
# from bokeh.models.glyphs import VBar
from astname import get_indexfile
import platform
import copy

warnings.filterwarnings("ignore")

def download_asteroid_data():
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


def calc_visibility(name, otime, observing_location, observatoryid):
    fmt = '%Y-%m-%d %H:%M:%S'
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
    obj = Horizons(id=name, location=observatoryid, epochs=jd)
    target = obj.ephemerides()
    center = SkyCoord(target['RA'][0]*u.deg, target['DEC'][0]*u.deg,
                        frame='icrs')
    center = center.transform_to(AltAz(obstime=otime,
                                    location=observing_location))

    print(f"utcoffset (h): {utcoffset}")
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
    for i in range(len(delta_midnight)):
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
    source3.data = {'xs':delta_midnight.value,
                   'ys': center2.alt.value,
                   'labels': label3,
                   'color': color3,
                   'az': center2.az.value,
                   'time': times.value}
    interval4 = source1.data['xs'][sunaltazs.alt < -0*u.deg]
    center4 = interval4[int(len(interval4)/2)]
    interval5 = source1.data['xs'][sunaltazs.alt < -12*u.deg]
    center5 = interval5[int(len(interval5)/2)]
    source4.data = {'center': [center4], 'width': [max(source1.data['xs'][sunaltazs.alt < -0*u.deg])-min(source1.data['xs'][sunaltazs.alt < -0*u.deg])]}
    source5.data = {'center': [center5], 'width': [max(source1.data['xs'][sunaltazs.alt < -12*u.deg])-min(source1.data['xs'][sunaltazs.alt < -12*u.deg])]}


def callback(attr, old, new):
    if sname.value == "":
        sname.value = '1'
    if date.value == "":
        date.value = "2023-05-05"
    if autoc.value == "":
        autoc.value = "Piszkesteto Stn. (Konkoly)"
    print(f"sname.value: {sname.value}")
    orig = copy.copy(sname.value)
    sname.value = str(sname.value).split()[0]
    print(f"sname.value after split: {sname.value}")
    obsdata = get_obsdata(observatories, autoc.value)
    print(f"obsdata:\n{obsdata}")
    observing_location = EarthLocation(lat=obsdata['gclat'][0], lon=obsdata['long'][0],
                                   height=obsdata['height'][0]*u.m)
    calc_visibility(sname.value, date.value, observing_location, obsdata['code'][0])
    p.title.text = f"Visibility plot for {orig} at {date.value} from {autoc.value}"

def search_callback(attr, old, new):
    # get the search term from the input widget
    term = new.lower()
    print(f"term: {term}")

    # filter the data based on the search term
    filtered_data = {'ssos': []}
    for i in range(len(data['ssos'])):
        row = {'ssos': data['ssos'][i]}

        # loop over each value in the row
        for value in row.values():
            value_str = str(value).lower()

            # check if the value contains the search term
            if term in value_str:
                filtered_data['ssos'].append(row['ssos'])
                print(f"Match! {term} {value_str}")
                break

    # update the data source with the filtered data
    # source.data = filtered_data
    # print(f"source.data:\n{source.data}")
    # print(f"search.value: {search.value}")
    source0.data = {'ssos': filtered_data['ssos']}
    sname.options = filtered_data['ssos']
    # sname.value = filtered_data['ssos'][0]

# def create_renderers(colormap, name):
#     xs = [source3.data['xs']]
#     ys = [source3.data['ys']]
#     l = source3.data['az']

#     color_mapper = LinearColorMapper(palette=colormap, 
#         low=min(l), high=max(l))

#     glyph = MultiLine(xs="xs", ys="ys", 
#                       line_color={'field': 'az', 'transform': color_mapper}, 
#                       line_width=2, name=name)

#     color_bar = ColorBar(color_mapper=color_mapper, title=name)
#     return glyph, color_bar


def read_obsdata():
    # https://dc.zah.uni-heidelberg.de/obscode/q/query/form
    t = Table.read('observatories.csv', format='ascii.csv')
    return t

def get_obsdata(observatories, observatory):
    obs = observatories[observatories['name'] == observatory]
    return obs

def conv(s):
    return s.strip()

astdat_file = download_asteroid_data()

astdat = np.loadtxt(astdat_file, delimiter=",", dtype=str, usecols=0, converters={0: conv}).tolist()
# print(f"astdat: {astdat[0:10]}")
data = {'ssos': astdat}
source = ColumnDataSource(data)

user = getpass.getuser()
if os.path.exists("/tmp/"+str(user)) is False:
    os.mkdir("/tmp/"+str(user))
plt.style.use(astropy_mpl_style)
quantity_support()
hover = HoverTool(line_policy='interp',
                  tooltips=[("object", "@labels"), ("UTC", "@time"), ("delta Time from midnight", "@xs"), ("altitude", "@ys"), ("azimuth", "@az")],
                  point_policy='snap_to_data')

p = figure(height=1000, width=1200, y_range=(0, 90) , x_range=(-10, 10))#, tools=[hover]) #  , x_range=(-12, 12)

observatories = read_obsdata()

source0 = ColumnDataSource(data={'ssos':['1']})
source1 = ColumnDataSource(data={'xs':[], 'ys': [], 'labels': []})
source2 = ColumnDataSource(data={'xs':[], 'ys': [], 'labels': []})
source3 = ColumnDataSource(data={'xs':[], 'ys': [], 'labels': ['1'], 'color': [], 'az': [1.], 'time': [datetime(2023, 5, 5), datetime(2023, 5, 6)]})
source4 = ColumnDataSource(data={'center': [0], 'width': [0]})
source5 = ColumnDataSource(data={'center': [0], 'width': [0]})

date = DatePicker(title='Choose date', value="2023-05-05")
name = TextInput(title='Write target name')
# sname = AutocompleteInput(title='Search Target:', completions=astdat, case_sensitive=False, placeholder="Search for an asteroid.")
sname = Select(title='Choose a target:', options=source0.data['ssos'], value=f"{source0.data['ssos'][0]}")
autoc = AutocompleteInput(title='Search Observatory:', completions=observatories['name'].tolist(), case_sensitive=False)

date.on_change('value', callback)
name.on_change('value', search_callback)
sname.on_change('value', callback)
autoc.on_change('value', callback)

top_row = row(date, name, sname, autoc)
middle_row = row(p)

p.vbar(source=source5, x='center' ,width='width', bottom=0,top=90, color=d3['Category20c'][20][-1], fill_alpha=0.8, line_color=None)
p.vbar(source=source4, x='center' ,width='width', bottom=0,top=90, color=d3['Category20c'][20][16], fill_alpha=0.4, line_color=None)

color_mapper = LinearColorMapper(palette=Viridis256, low=0, high=360)
p.line(x='xs', y='ys', source=source1, legend_label="Sun", color="orange", line_width=1)
p.line(x='xs', y='ys', source=source2, legend_label="Moon", color=d3['Category20c'][4][1], line_width=2, line_dash='dashed')
p.circle(x='xs', y='ys', source=source3, legend_field="labels", color={'field': 'az', 'transform': color_mapper})
hline = Span(location=20, dimension='width', line_color='red', line_width=2, line_dash='dashed', line_alpha=0.5)
p.renderers.extend([hline])
p.add_tools(HoverTool(line_policy='interp',
                  tooltips=[("object", "@labels"), ("UTC", "@time"), ("delta Time from midnight", "@xs"), ("altitude", "@ys"), ("azimuth", "@az")],
                  point_policy='snap_to_data'))
p.xaxis.axis_label = 'Hours from Local Midnight'
p.yaxis.axis_label = 'Altitude (deg)'
p.legend.click_policy="hide"
bar = ColorBar(color_mapper=color_mapper, location=(0,0), title="Azimuth (deg)")
p.add_layout(bar, "right")
p.title.text = f"Visibility plot for {str(source3.data['labels'][0])} at {date.value} from {autoc.value}"
curdoc().add_root(column(top_row, middle_row))
curdoc().theme = 'caliber'
curdoc().title = "Visibility plot"
