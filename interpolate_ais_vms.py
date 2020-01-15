# -*- coding: utf-8 -*-
"""
Created on Wed May 29 11:26:59 2019

@author: WSBARATA01
"""

import os, glob, shutil
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import nearest_points
from datetime import datetime, timedelta
from math import radians, asin, sin, cos, atan2, degrees, isclose
import pyproj
from zipfile import ZipFile
from  tkinter import Tk, filedialog

base_path = "D:\\BARATA"

#input directory path
Tk().withdraw()
#data_folder = filedialog.askdirectory(initialdir=f'{base_path}\\2.seonse_outputs',title='Pick a directory')
data_folder = filedialog.askdirectory(initialdir=f'{base_path}\\2.seonse_outputs',title='Pick a directory')[:-4] + "*"

data_list = glob.glob(data_folder)

data_date = os.path.basename(data_list[-1])[-15:-7]

aiszip_path = glob.glob(f'{base_path}\\10.ais\\*\\*{data_date}*.zip')[0]
ais_path = f'{os.path.dirname(aiszip_path)}\\indo_{os.path.basename(aiszip_path)[6:14]}_ais.csv'
vms_path = glob.glob(f'{base_path}\\9.vms\\*\\*{data_date}*csv')[0]
vms_info_path = glob.glob(f'{base_path}\\9.vms\\vms_info_fix.csv')[0]

boundary_threshold = 0.5 #km
distance_threshold = 0.5 #km
length_threshold = 1 #ratio
time_threshold = 30 #minutes

def calculate_coordinates(departure_point, bearing, distance):
    R = 63781000
    d = distance

    lat1 = radians(departure_point[1]) #Current lat point converted to radians
    lon1 = radians(departure_point[0]) #Current long point converted to radians

    bearing = radians(bearing)

    lat2 = asin(sin(lat1)*cos(d/R) + cos(lat1)*sin(d/R)*cos(bearing))

    lon2 = lon1 + atan2(sin(bearing)*sin(d/R)*cos(lat1), cos(d/R)-sin(lat1)*sin(lat2))

    lat2 = degrees(lat2)
    lon2 = degrees(lon2)

    return (round(lon2, 4), round(lat2, 4))

def nearest(row, geom_union, df1, df2, geom1_col='geometry', geom2_col='geometry', src_column=None):
    """Find the nearest point and return the corresponding value from specified column."""
    # Find the geometry that is closest
    nearest = df2[geom2_col] == nearest_points(row[geom1_col], geom_union)[1]
    # Get the corresponding value from df2 (matching is based on the geometry)
    value = df2[nearest][src_column].array[0]
    
    return value

def distance(row, geom_union, geom1_col='geometry', geom2_col='geometry'):
    nearest = nearest_points(row[geom1_col], geom_union)
    #distance = sqrt((nearest[0].x - nearest[1].x)**2+(nearest[0].y - nearest[1].y)**2)
    geod = pyproj.Geod(ellps='WGS84')
    angle1,angle2,distance = geod.inv(nearest[0].x, nearest[0].y, nearest[1].x, nearest[1].y)
    
    return distance
    
def interpolate(gdf, lon_column, lat_column, time_column, heading_column, speed_column):
    time_list = []
    head_list = []
    speed_list = []
    intpol_list = []
    
    for i, row in gdf.iterrows():
        date = row[time_column]
        targetdate = datetime.strptime(date, '%Y-%m-%dT%H:%M:%S')
        
        if date <= shipdate.strftime('%Y-%m-%dT%H:%M:%S'):
            timedelta = shipdate - targetdate
        else:
            timedelta = targetdate - shipdate
            
        speed = (row[speed_column])*0.514444
        
        dist = speed*timedelta.seconds
        
        head = row[heading_column]
        
        geom = row.geometry
        geom_tuple = (geom.x, geom.y)
        
        next_geom_tuple = calculate_coordinates(geom_tuple, head, dist)
        
        next_geom = Point(next_geom_tuple)
        
        time_list.append(str(shipdate))
        head_list.append(row[heading_column])
        speed_list.append(row[speed_column])
        
        intpol_list.append(next_geom)
    
    gdf[lon_column] = [i.x for i in intpol_list]
    gdf[lat_column] = [i.y for i in intpol_list]
    gdf.geometry = intpol_list
    return gdf

def df_to_gdf(df, x_column, y_column):
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[x_column], df[y_column]))
    return gdf

def feat_extent(gdf):
    xmin = gdf.geometry.x.min() - boundary_threshold
    xmax = gdf.geometry.x.max() + boundary_threshold
    ymin = gdf.geometry.y.min() - boundary_threshold
    ymax = gdf.geometry.y.max() + boundary_threshold
    return xmin,xmax,ymin,ymax

def filter_time(gdf, time_column):
    ftimegdf = gdf.loc[(gdf[time_column] >= startdate) & (gdf[time_column] <= stopdate)]
    return ftimegdf

def vms_gt(beacon):
    try:
        gt = int(vms_info_df[vms_info_df['beacon'] == str(beacon)]['gt'].array[0])
    except:
        gt = None
    return gt
    

#perform script    
for data_path in data_list:
    print ('-------------------------')
    print (os.path.basename(data_path))
    ship_list = glob.glob(f'{data_path}\\*SHIP.shp')
    if len(ship_list) > 0:
        ship_path = ship_list[0]
        
        #mendefinisikan tanggal data
        shipdate = datetime.strptime(data_path[-15:], '%Y%m%d_%H%M%S')
        startdate = (shipdate - timedelta(minutes=time_threshold)).strftime('%Y-%m-%dT%H:%M:%S')
        stopdate = (shipdate + timedelta(minutes=time_threshold)).strftime('%Y-%m-%dT%H:%M:%S')
        
        #mengekstrak data csv di dalam zip
        if not os.path.exists(ais_path):
            with ZipFile(aiszip_path) as theZip:
                fileNames = theZip.namelist()
                for fileName in fileNames:
                    if fileName.endswith('csv'):
                        with theZip.open(fileName) as f:
                            with open(ais_path, 'wb') as outfile:
                                shutil.copyfileobj(f, outfile)
        else:
            pass
        
        #membaca data ais dan vms
        ais_df = pd.read_csv(ais_path)
        vms_df = pd.read_csv(vms_path)
        vms_info_df = pd.read_csv(vms_info_path)

        vms_df['Location date'] = [(datetime.strptime(i, '%m/%d/%Y %H:%M:%S')).strftime('%Y-%m-%dT%H:%M:%S') for i in vms_df['Location date']]
        
        #konvert ke geodataframe
        ship_gdf = gpd.read_file(ship_path)
        ais_gdf = df_to_gdf(ais_df, 'longitude', 'latitude') 
        vms_gdf = df_to_gdf(vms_df, 'Longitude', 'Latitude')
        
        #mendapatkan informasi batas koordinat kapal
        xmin,xmax,ymin,ymax = feat_extent(ship_gdf)
        
        #filter data berdasarkan waktu
        ais_filter1 = filter_time(ais_gdf, 'time')
        vms_filter1 = filter_time(vms_gdf, 'Location date')
        
        #filter data berdasarkan lokasi setelah filter waktu
        ais_filter2 = ais_filter1.cx[xmin:xmax, ymin:ymax]
        vms_filter2 = vms_filter1.cx[xmin:xmax, ymin:ymax]
        
        #filter data berdasarkan posisi terakhir
        ais_filter3 = ais_filter2.groupby('mmsi', group_keys=False).apply(lambda x: x[x['time'] == min(x['time'], key=lambda y: abs(pd.to_datetime(y) - shipdate))])
        vms_filter3 = vms_filter2.groupby('Beacon ID', group_keys=False).apply(lambda x: x[x['Location date'] == min(x['Location date'], key=lambda y: abs(pd.to_datetime(y) - shipdate))])
        
        #get gt information
        if not vms_filter3.empty:
            vms_filter3['GT'] = [vms_gt(beacon) for beacon in vms_filter3['Beacon ID']]
        else:
            vms_filter3 = vms_filter2
        
        #eksekusi
        if ais_filter3.empty == False or vms_filter3.empty == False:
            print ('AIS data is available around the area:', len(ais_filter3))
            print ('VMS data is available around the area:', len(vms_filter3))
            
            associated_gdf = ship_gdf.copy()
            
            if ais_filter3.empty == False and vms_filter3.empty == False:
                ais_intpol_gdf = interpolate(ais_filter3, 'longitude', 'latitude', 'time', 'cog', 'sog')
                vms_intpol_gdf = interpolate(vms_filter3, 'Longitude', 'Latitude', 'Location date', 'Heading', 'Speed')
                
                ais_union = ais_intpol_gdf.unary_union
                vms_union = vms_intpol_gdf.unary_union
                
                associated_gdf['NEAREST_MMSI'] = associated_gdf.apply(nearest, geom_union=ais_union, df1=associated_gdf, df2=ais_intpol_gdf, src_column='mmsi', axis=1)
                associated_gdf['NEAREST_BEACON'] = associated_gdf.apply(nearest, geom_union=vms_union, df1=associated_gdf, df2=vms_intpol_gdf, src_column='Beacon ID', axis=1)
                
                associated_gdf['AIS_LENGTH'] = associated_gdf.apply(nearest, geom_union=ais_union, df1=associated_gdf, df2=ais_intpol_gdf, src_column='length', axis=1)
                associated_gdf['VMS_GT'] = associated_gdf.apply(nearest, geom_union=vms_union, df1=associated_gdf, df2=vms_intpol_gdf, src_column='GT', axis=1)
                
                associated_gdf['AIS_DISTANCE'] = associated_gdf.apply(distance, geom_union=ais_union, axis=1)
                associated_gdf['AIS_DISTANCE'] = associated_gdf['AIS_DISTANCE']/1000
                
                associated_gdf['VMS_DISTANCE'] = associated_gdf.apply(distance, geom_union=vms_union, axis=1)
                associated_gdf['VMS_DISTANCE'] = associated_gdf['VMS_DISTANCE']/1000
                
                len_tol = [isclose(i,j,rel_tol=length_threshold) for i,j in zip(associated_gdf['AIS_LENGTH'], associated_gdf['LENGTH'])]
                
                associated_gdf['STATUS'] = ['AIS' if x <= distance_threshold and y is True else None for x,y in zip(associated_gdf['AIS_DISTANCE'],len_tol)]
                associated_gdf['STATUS'] = ['VMS' if x <= distance_threshold and y == None else None if x > distance_threshold and y == None else 'AIS' for x,y in zip(associated_gdf['VMS_DISTANCE'],associated_gdf['STATUS'])]
            
            elif ais_filter3.empty == False and vms_filter3.empty == True:
                ais_intpol_gdf = interpolate(ais_filter3, 'longitude', 'latitude', 'time', 'cog', 'sog')
                ais_union = ais_intpol_gdf.unary_union
                
                associated_gdf['NEAREST_MMSI'] = associated_gdf.apply(nearest, geom_union=ais_union, df1=associated_gdf, df2=ais_intpol_gdf, src_column='mmsi', axis=1)
                associated_gdf['AIS_LENGTH'] = associated_gdf.apply(nearest, geom_union=ais_union, df1=associated_gdf, df2=ais_intpol_gdf, src_column='length', axis=1)
                
                associated_gdf['AIS_DISTANCE'] = associated_gdf.apply(distance, geom_union=ais_union, axis=1)
                associated_gdf['AIS_DISTANCE'] = associated_gdf['AIS_DISTANCE']/1000
                
                len_tol = [isclose(i,j,rel_tol=length_threshold) for i,j in zip(associated_gdf['AIS_LENGTH'], associated_gdf['LENGTH'])]
                
                associated_gdf['STATUS'] = ['AIS' if x <= distance_threshold and y is True else None for x,y in zip(associated_gdf['AIS_DISTANCE'],len_tol)]
            
            else:
                vms_intpol_gdf = interpolate(vms_filter3, 'Longitude', 'Latitude', 'Location date', 'Heading', 'Speed')
                vms_union = vms_intpol_gdf.unary_union
                
                associated_gdf['NEAREST_BEACON'] = associated_gdf.apply(nearest, geom_union=vms_union, df1=associated_gdf, df2=vms_intpol_gdf, src_column='Beacon ID', axis=1)
                associated_gdf['VMS_GT'] = associated_gdf.apply(nearest, geom_union=vms_union, df1=associated_gdf, df2=vms_intpol_gdf, src_column='GT', axis=1)
                
                associated_gdf['VMS_DISTANCE'] = associated_gdf.apply(distance, geom_union=vms_union, axis=1)
                associated_gdf['VMS_DISTANCE'] = associated_gdf['VMS_DISTANCE']/1000
                
                associated_gdf['STATUS'] = ['VMS' if x <= distance_threshold else None for x in associated_gdf['VMS_DISTANCE']]

            associated_aisgdf = associated_gdf[associated_gdf['STATUS'] == 'AIS']
            associated_vmsgdf = associated_gdf[associated_gdf['STATUS'] == 'VMS']  
            
            #remove vms duplicate
            associated_vmsfiltergdf = associated_vmsgdf.groupby('NEAREST_BEACON', group_keys=False).apply(lambda x: x[x['VMS_DISTANCE'] == min(x['VMS_DISTANCE'], key=lambda y: y)])
            
            associated_ais = len(associated_aisgdf)
            associated_vms = len(associated_vmsfiltergdf)
            
            associated_aisvmsgdf = gpd.GeoDataFrame(pd.concat([associated_aisgdf, associated_vmsfiltergdf], ignore_index=True))
            
            associated_shipgdf = pd.merge(ship_gdf, associated_aisvmsgdf, how='left')
            
            #convert data type to object
            if set(['NEAREST_MMSI', 'NEAREST_BEACON']).issubset(associated_shipgdf.columns):
                try:
                    associated_shipgdf.iloc[:,16:18] = associated_shipgdf.iloc[:,16:-5].astype('Int64').astype(object)
                except:
                    pass
            else:
                try:
                    associated_shipgdf.iloc[:,16:17] = associated_shipgdf.iloc[:,16:-5].astype('Int64').astype(object)
                except:
                    pass


            if associated_ais > 0 or associated_vms > 0:
                print(associated_aisvmsgdf)
                print ('Total of ship associated with AIS:', associated_ais)
                print ('Total of ship associated with VMS:', associated_vms)
            else:
                print ('No ship associated with AIS and VMS')
                
            export = input('Eksport atau tidak (y/n): ')
            if export == 'y':
                correlated_basepath = os.path.dirname(ship_path).replace("2.seonse_outputs","12.correlated_ship")
                output_path = f'{correlated_basepath}\\interpolated\\{os.path.basename(ship_path)[:-4]}_CORRELATED.shp'
                ais_output_path = f'{os.path.dirname(output_path)}\\{os.path.basename(os.path.dirname(ship_path))}ais.csv'
                vms_output_path = f'{os.path.dirname(output_path)}\\{os.path.basename(os.path.dirname(ship_path))}_vms.csv'
                
                if not os.path.exists(os.path.dirname(output_path)):
                    os.makedirs(os.path.dirname(output_path))
                else:
                    pass
                associated_shipgdf.to_file(output_path)

                if ais_filter3.empty == False:
                    ais_intpol_gdf.to_csv(ais_output_path, index=False)
                else:
                    pass
                
                if vms_filter3.empty == False:
                    vms_intpol_gdf.to_csv(vms_output_path, index=False)
                else:
                    pass
                
                print ('Export successfully')
            else:
                pass
        else:
            print ('AIS and VMS data is unavailable around the area')
    else:
        print ('Ship data is unavailable')
