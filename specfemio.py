import os
import copy
import pickle

    
def write_stations_from_headers(fqp,stations):
    
    path = os.path.relpath(fqp, start=os.curdir)
    fqpname = os.path.join(path, 'STATIONS')
    header_fqpname = os.path.join(path, 'pyheader.stations')
    
    print('fqpname:', fqpname)

    str_stations = []
    for i in range(len(stations)):
        s = stations[i]
        sname = f'{s.name}_{s.trid}_{s.gid}_{s.sid}'
        snet  = s.network
        slat  = s.lat_yc # or Y coordinate
        slon  = s.lon_xc # or X coordinate
        selev = s.elevation
        sbur  = s.burial
        str_stations.append('%s %s %.2f %.2f %.2f %.2f\n' %(sname,snet,slat,slon,selev,sbur))

    f = open(fqpname, 'w')
    f.writelines(str_stations)
    f.close()
    
    f = open(header_fqpname, 'wb')
    pickle.dump(stations,f)
    f.close()
    

def read_stations_from_headers(fqp):
    
    path = os.path.relpath(fqp, start=os.curdir)
    header_fqpname = os.path.join(path, 'pyheader.stations')
    
    f = open(header_fqpname, 'rb')
    stations = pickle.load(f)
    f.close()
    
    return stations
    

def make_station_sym_cross_group(station=None,delta=None):
    
    station_group = []
    
    cpy_station = copy.deepcopy(station)
    cpy_station.gid = station.gid = 0
    station_group.append(cpy_station)
    
    # add x - delta
    station_xm = copy.deepcopy(station)
    station_xm.lon_xc -= delta
    station_xm.gid  = 1
    station_group.append(station_xm)
    
    # add x + delta
    station_xp = copy.deepcopy(station)
    station_xp.lon_xc += delta
    station_xp.gid  = 2
    station_group.append(station_xp)
    
    # add y - delta
    station_ym = copy.deepcopy(station)
    station_ym.lat_yc -= delta
    station_ym.gid  = 3
    station_group.append(station_ym)
    
    # add y + delta
    station_yp = copy.deepcopy(station)
    station_yp.lat_yc += delta
    station_yp.gid  = 4
    station_group.append(station_yp)
    
    # add z - delta
    station_zm = copy.deepcopy(station)
    station_zm.burial -= delta
    station_zm.gid  = 5
    station_group.append(station_zm)
    
    # add y + delta
    station_zp = copy.deepcopy(station)
    station_zp.burial += delta
    station_zp.gid  = 6
    station_group.append(station_zp)
    
    return station_group
    
    
def make_station_sym_cross_group_list(stations,delta):
    
    
    group_station_list = []
    for i in range(len(stations)):
        group_station = make_station_group(station=stations[i],delta=delta)
        #group_station_list = list(set(group_station_list + group_station))
        group_station_list += group_station
    
    return list(set(group_station_list))
