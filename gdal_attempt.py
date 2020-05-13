def execute_gdal(target_folder, fpath, res, fmt, epsg):
    from gdal import Grid, GridOptions
    with open(target_folder + "GDAL-IDWexposed-options.txt", 'r') as file_in:
        params = file_in.readlines()
    if fmt == "ASC":    
        base_params = ['format = "XYZ"',
                       'width = ' + str(res[0]),
                       'height = ' + str(res[1]),
                       'noData = -9999',
                       'layers = "idw"']
    if fmt == "GeoTIFF":    
        base_params = ['format = "GTiff"',
                       'width = ' + str(res[0]),
                       'height = ' + str(res[1]),
                       'noData = -9999',
                       'outputSRS = "EPSG:' + epsg + '"',
                       'layers = "idw"']
    for param in base_params: params.insert(0, param.strip("\n").strip("\t"))
    #def printer(p1, p2, p3, p4, p5, p6, p7):
    #    print(p1, p1, p3, p4, p5, p6, p7)
    #printer(*params)
    #options = GridOptions(*params)
    
    options = GridOptions(format='GTiff', width = 100, height = 100,
                          algorithm = 'invdist:power=2',
                          layers=[fpath[:-4] + '_gf.shp'])
    Grid(destName = fpath[:-4] + '_IDWexposed.tif',
              srcDS = fpath[:-4] + '_gf.shp',
              options = options)