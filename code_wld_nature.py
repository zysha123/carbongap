import ee
from functools import partial

# All datasets are avaiable at https://code.earthengine.google.com/?asset=users/zongyaosha/imgs2018
ee.Initialize()
 
# Asset location
asset_location="users/username/this_dir"
#Local file store for statistics files
file_location="/usr/temp"

#MOD17 time series
day8_psnnet = ee.ImageCollection("MODIS/006/MOD17A2H")
#Recoded Harmonized World Soil Dataset from HWSD.
hwsd = ee.Image(asset_location+"/soil/recoded_hwsd")
#Land cover types from MCD12Q1 time series
landcoverMOD12Q1 = ee.ImageCollection("MODIS/006/MCD12Q1")
#World land cover area. For filtering processing region
roi = ee.FeatureCollection(asset_location+"/world_land_area")
#Landforms layer, from https://esdac.jrc.ec.europa.eu
landforms = ee.Image(asset_location+"/Iwahashi")
#Monthly temperture and precipitation, for processing potential NPP (PNPP)
climate = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
#This feature collection contains countries with continent attribute, from which statistics on continent scale are possible
countries2 = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
#Uniform spatial scale
scaleSize = 500 #pixel size

landforms_clip = landforms.clip(roi).reproject('EPSG:4326', None, scaleSize)
soil_clip = hwsd.clip(roi).reproject('EPSG:4326', None, scaleSize)

#Optimal window distance ~20km estimated from function estimateDistance
neiborghhoodSize=40 # 40pixels=20km
neighborhoodKernel=ee.Kernel.circle(neiborghhoodSize)
#percentile 90, for computing the target maximum NPP within each landform+soil+vegetation cover neighborhood
percentile=ee.Reducer.percentile([90])



#-----------------------------Prepare (preprocess) datasets--------------------------------------------
# Export land (vegetation) cover from landcoverMOD12Q1 and save to local asset for convenience
def saveLandCover():
    resultmap=[None]*18
    task=[None]*18
    for i in range(18): #2001-2018
       year=str(2001+i)
       resultmap[i] =  landcoverMOD12Q1.filterDate(year+'-1-1', year+'-12-31').select('LC_Type1').first().clip(roi).reproject('EPSG:4326', None, scaleSize)
       print("begin saveLandCover task:"+year)
       task[i] = ee.batch.Export.image.toAsset(
           image=resultmap[i],
           description='landcoverMOD12Q1_'+year,
           assetId=asset_location+'/landcoverMOD12Q1_'+year,
           scale=scaleSize,
           maxPixels=1e13,
           crs='EPSG:4326',
           region=roi.geometry()
       )
       task[i].start()
       print("end saveLandCover task: "+year)
#saveLandCover()

#annual total NPP
def getNPP(year):
    npp = (day8_psnnet.filterDate(year + '-1-1', year + '-12-31').select('PsnNet')).sum().multiply(0.1).int16()
    return npp.reproject('EPSG:4326', None, scaleSize)

#Calculate NPP and save to user asset
def saveNPP():
    resultmap=[None]*18
    task=[None]*18
    for i in range(18): #2001-2018
       year=str(2001+i)
       resultmap[i] =getNPP(year)
       print("begin npp task:"+year)
       task[i] = ee.batch.Export.image.toAsset(
           image=resultmap[i],
           description='NPP'+year+"new",
           assetId=asset_location+'/NPP'+year,
           scale=scaleSize,
           maxPixels=1e13,
           crs='EPSG:4326',
           region=roi.geometry()
       )
       task[i].start()
       print("end npp task: "+year)
#saveNPP()

# Segmentation LVS patches based on landform (L),vegType(V) and soiltype(S)
# Region segmentation of homogeneous natural environments based on LVS
def saveUniqueUnits():
    resultmap=[None]*18
    task=[None]*18
    for i in range(18): #2001-2018
       year=str(2001+i)
       resultmap[i] =createUniqueUnits(year) # call function
       print("begin createUniqueUnits task:"+year)
       task[i] = ee.batch.Export.image.toAsset(
           image=resultmap[i],
           description='UniqueUnits'+year,
           assetId=asset_location+'/UniqueUnits'+year,
           scale=scaleSize,
           maxPixels=1e13,
           crs='EPSG:4326',
           region=roi.geometry()
       )
       task[i].start()
       print("end createUniqueUnits task: "+year)
#saveUniqueUnits()

# Each DN value in the image is uniquely identified by the combination of landform (L),vegType(V) and soiltype(S)
def createUniqueUnits(year):
    # make soil ids 3 digits, i.e., from 1-450 to 100-550(3 digits)
    x4 = soil_clip.add(ee.Image.constant(99))
    # make landform 2 digits, i.e., original land form codes 1-16 were reassigned to be 10-25 (2 digits)
    x6 = landforms_clip.add(9)
    # land cover
    landcover = ee.Image(asset_location+"/landcoverMOD12Q1_"+str(year))
    # make landcover id into 2 digits, i.e., land cover code 1-17 converted to id: 10-26 (2 digits)
    x8 = landcover.add(9)
    # Then we will be able to have unique value for each combination of land forms, soil types and vegetation types: landforms*100000+soil*100+landcover, range: 1010010~2969929
    output_landforms_soil = x6.multiply(1000).add(x4)
    # get unique value for the combination of land forms and soil types
    output = output_landforms_soil.multiply(100).add(x8)  # combine output_landforms_soil and landcover
    # unique id=land forms *100000 + soil types*100 + land cover types
    return output.reproject('EPSG:4326', None, scaleSize).rename(['uid'])
#createUniqueUnits(2001)

#average monthly temperature for processing Potential NPP (PNPP)
#note that deriving average temperture using tmin and tmax might be biased
# but we eventually evaluate the relative contribution (RC) to carbon sequestration from climate factors between locations
# the bias will have very small impact on the RC
def averageMonthlyTemp(year,i):
   climate_tmmx = climate.filterDate(year + '-' + str(i + 1) + '-1', year + '-' + str(i + 1) + '-31').select("tmmx").first()
   climate_tmmn = climate.filterDate(year + '-' + str(i + 1) + '-1', year + '-' + str(i + 1) + '-31').select("tmmn").first()
   r = (climate_tmmx.add(climate_tmmn).divide(2)).rename(['tm_avg'])
   return r.reproject('EPSG:4326', None, scaleSize)

# biomass to carbon conversion coefficiet based on land (vegetation) cover type.0.50 for woody ecosystems and 0.45 for herbaceous
def getConversionCoefficient(landcover):
    return landcover.remap([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45])

#Potential NPP using Miami model
def computePNPP(year):
    #Lopp through each month of the year
    for i in range(12):
        #monthly total precipitation
        climate_pr = climate.filterDate(year + '-'+str(i)+'-1', year + '-'+str(i)+'-31').select("pr").sum()
        x1 = ((climate_pr.multiply(-0.000664)).exp().add(-1)).multiply(-3000)
        x22 = (averageMonthlyTemp(year,i).multiply(-0.119).multiply(0.1).add(1.315)).exp().add(1)
        x2 = ee.Image.constant(3000).divide(x22)
        landcover = ee.Image(asset_location+"/landcoverMOD12Q1_"+str(year))
        coef=getConversionCoefficient(landcover)
        if pnpp is None: # first month
            pnpp = x1.min(x2).multiply(coef).rename(['pnpp']).int16()
        else: # monthly sum
            pnpp = pnpp.add(x1.min(x2).multiply(coef).rename(['pnpp']).int16())
    return pnpp.reproject('EPSG:4326', None, scaleSize)

def savePNPP():
    resultmap=[None]*18
    task=[None]*18
    for i in range(18): #2001-2018
       year=str(2001+i)
       resultmap[i] =computePNPP(year)
       print("begin pnpp task:"+year)
       task[i] = ee.batch.Export.image.toAsset(
           image=resultmap[i],
           description='PNPP'+year,
           assetId=asset_location+'/PNPP'+year,
           scale=scaleSize,
           maxPixels=1e13,
           crs='EPSG:4326',
           region=roi.geometry()
       )
       task[i].start()
       print("end pnpp task: "+year)
#savePNPP()

#Relative NPP (i.e., relative contribution from climate impact), which records the relative observed NPP rectified by PNPP
#we could process RC for each LVS zone, but it is more efficient to process the whole region all together
#then RC is taken to process each LVS zone
def getRelativeNPP(year):
    npp = getNPP(year)
    pnpp = computePNPP(year)
    # remove climate effect. The result is the relative contribution effect from climate impact between locations
    bias = npp.subtract(pnpp).int16()
    return bias.reproject('EPSG:4326', None, scaleSize)

#save RelativeNpp (Relative contribution from climate impact)
def saveRelativeNpp():
    resultmap=[None]*18
    task=[None]*18
    for i in range(18): #2001-2018
       year=str(2001+i)
       resultmap[i] =getRelativeNPP(year)
       print("begin task:"+year)
       task[i] = ee.batch.Export.image.toAsset(
           image=resultmap[i],
           description='relativeNpp'+year,
           assetId=asset_location+'/relativeNpp'+year,
           scale=scaleSize,
           maxPixels=1e13,
           crs='EPSG:4326',
           region=roi.geometry()
       )
       task[i].start()
       print("end task: "+year)
#saveRelativeNpp()
#---------------------------------------End Prepare (preprocess) datasets-------------------------------


#------------------------------------------Compute carbon gap-------------------------------------------
#Compute carbon gap from preprocessed images, including relativeNpp-yyyy and UniqueUnits-yyyy
#theory behind: in each unique LVS zone constrained by neighborhoodKernel, the relative NPP should have no internal variation
# if there is no land management difference; conversely, if variations show up, then the higher relative NPP indicates optimized land mangement
# Those locations showing low relative NPP in the LVS zone can also adopt optimized land management practices so NPP there will be improved
def computeResult(year):
    relativeNpp = ee.Image(asset_location+"/relativeNpp"+str(year))
    uniqueunit = ee.Image(asset_location+"/UniqueUnits"+str(year))

    #convert neighborhood of uniqueunit to bands to facilitate the analysis
    nbhBands_unitClass=uniqueunit.neighborhoodToBands(neighborhoodKernel)

    #function for filtering particular image layers from the bands
    def getUnitClassBands(bandName):
        return nbhBands_unitClass.select([bandName])

    # map
    bandList_unitClass = list(map(getUnitClassBands, (nbhBands_unitClass.bandNames().getInfo())))

    imgCollection_unitClass=ee.ImageCollection.fromImages(bandList_unitClass)

    def maskUnitClassImg(img):
      #using unique_units to mask neighborhood pixel layers
      return img.neq(uniqueunit) #Create a binary mask

    imgCollection_unitClass_masked=imgCollection_unitClass.map(maskUnitClassImg)

    def maskNPPImg(img):
      return img.mask(imgCollection_unitClass_masked.filterMetadata('system:index','equals',img.get('system:index')).first()).rename(['b'])

    nbhBands_npp=relativeNpp.neighborhoodToBands(neighborhoodKernel)
    def getNppBands(bandName):
        return nbhBands_npp.select([bandName])

    bandList_npp = list(map(partial(getNppBands), nbhBands_npp.bandNames().getInfo()))

    imgCollection_npp=ee.ImageCollection.fromImages(bandList_npp)
    imgCollection_npp_masked=imgCollection_npp.map(maskNPPImg)

    focal_pct90 = imgCollection_npp_masked.reduce(reducer=percentile).rename('npp')
    gap=focal_pct90.subtract(relativeNpp)
    mask=gap.gt(ee.Image.constant(0)) # mask out the pixel values that are less than 0
    # pixel value>0 will be kept and value=0 will be masked out; pixel values that are less than 0 will be replaced with 0
    gap=gap.mask(mask).unmask(0) # replace with 0

    print("begin task: "+str(year))
    # export carbon gap images to my asset
    task = ee.batch.Export.image.toAsset(
        image=gap,
        description='Gap'+str(year),
        assetId=asset_location+'/Gap'+str(year),
        scale=scaleSize,
        maxPixels=1e13,
        crs='EPSG:4326',
        region=roi.geometry()
    )
    task.start()
    print("end task： "+str(year))

#compute every year
def computeResultAll():
    for i in range(18):
        computeResult(2001+i)
#computeResultAll()

# compute the ration of carbon gap to NPP
def saveGapRatio():
    resultmap=[None]*18
    task=[None]*18
    for i in range(18): #2000-2018
       year=str(2001+i)
       npp=ee.Image(asset_location+"/NPP"+str(year))
       gap=ee.Image(asset_location+"/Carbongap"+str(year))
       resultmap[i] =  gap.clip(roi).divide(npp).multiply(100)
       print("begin GapRatio task:"+year)
       task[i] = ee.batch.Export.image.toAsset(
           image=resultmap[i].round().int8(),
           description='GapRatio'+year,
           assetId=asset_location+'/GapRatio'+year,
           scale=scaleSize,
           maxPixels=1e13,
           region=roi.geometry()
       )
       task[i].start()
       print("end GapRatio task: "+year)
#saveGapRatio()

#-------------------------------------End Compute carbon gap-------------------------------------------


#----------Optimal window size (~20km). Data is output from the function, which can them be analyzed--------------------
# estimate the constrained distance
def estimateDistance(year):
    #save the statistics to a file
    file_name=file_location+"/distance.csv"
    landcover_full= ee.Image(asset_location+"/landcoverMOD12Q1_2018")
    randomsamples = ee.FeatureCollection(asset_location+"/randomsamples2")
    # 5000 samples (features)
    for id in range(5000):
        sample_pnt=randomsamples.filter(ee.Filter.eq('id', id)).first()
        if sample_pnt is None:
            break
        sample_pnt_Info=sample_pnt.getInfo()
        land_id=sample_pnt_Info['properties']["LC_Type1"]

        #if land_id==11 or land_id==13 or land_id>=15: # land cover id=11, 13, 15,16,17 are not vegetated
        #largest distance, here=200pixels=100km
        max_distance_pixels=200 # set a miximum distance as 100km. Disadvantage of too long distance is discussed in the m.s.
        coords=sample_pnt_Info["geometry"]["coordinates"]
        sample_pnt=ee.Geometry.Point(coords)

        buffer_Area=sample_pnt.buffer(max_distance_pixels*500) #pixel size=500 m
        landcover = landcover_full.clip(buffer_Area).reproject('EPSG:4326', None, scaleSize)
        #get the area around the current point
        gap= ee.Image(asset_location+"/relativeNpp"+str(year)).clip(buffer_Area)
        #zonal analysis needs to select pixels that have the same land cover id as that of the current location

        # one point, one line showing PCT90 changes with the buffer distance
        gap_with_buffer=""
        for buff in range(1,max_distance_pixels,1):
            geom_buffer=sample_pnt.buffer((buff+1)*500-1) #pixel to meters by multiplying by 500m
            data = gap.clip(geom_buffer).mask(landcover.eq(ee.Image.constant(land_id))).reduceRegion(ee.Reducer.percentile([90]),geom_buffer,500)
            #data = gap.clip(geom_buffer).mask(landcover.eq(ee.Image.constant(land_id))).reduceRegion(ee.Reducer.toList(),geom_buffer,500)
            dataN = ee.Number(data)
            x=(dataN.getInfo())
            gap_with_buffer+=str(x["PsnNet"])+","
        f = open(file_name, "a")
        f.write(str(id)+","+str(land_id)+","+str(coords[0])+","+str(coords[1])+","+gap_with_buffer+"\n")
        f.close()
        #Analyze statistics offline. PCT90 changes with the distance.
#estimateDistance(2018)
#------End Optimal window size (~20km). Data is output from the function, which can them be analyzed--------------------

#---------------------------------------Compute Statistics--------------------------------------------------------------
#Pixel area of each location
def createPixelArea():
    pa=ee.Image.pixelArea().clip(roi).reproject(crs='EPSG:4326',scale=scaleSize)
    task = ee.batch.Export.image.toAsset(
               image=pa,
               description='pixelArea',
               assetId=asset_location+'/pixelArea',
               scale=scaleSize,
               maxPixels=1e13,
               crs='EPSG:4326',
               region=roi.geometry()
     )
    task.start()
#createPixelArea()

# union land (vegetation) cover area for all the years (without need to be common in land cover types).
# there might be a small area "loss" during the union if compared to the area of each year. DN=1 denotes MOD vegetated land cover area
def SaveVegetatedArea():
    #exclude land cover classes 11,13,15-17
    #mask=gap.gt(ee.Number(previous_break)).And(gap.lte(v[t]))
    temp_mask=ee.Image(asset_location+"/landcoverMOD12Q1_"+str(2001)).clip(roi)
    mask= temp_mask.lt(ee.Number(11)).Or(temp_mask.eq(ee.Number(12))).Or(temp_mask.eq(ee.Number(14)))
    for i in range(17):
        temp_mask=ee.Image(asset_location+"/landcoverMOD12Q1_"+str(2002+i)).clip(roi)
        temp_mask= temp_mask.lt(ee.Number(11)).Or(temp_mask.eq(ee.Number(12))).Or(temp_mask.eq(ee.Number(14)))
        mask= mask.And(temp_mask)
    task = ee.batch.Export.image.toAsset(
           image=mask,
           description='Vegetated_AllYears',
           assetId=asset_location+'/Vegetated_AllYears',
           scale=scaleSize,
           maxPixels=1e13,
           crs='EPSG:4326',
           region=roi.geometry()
       )
    task.start()
#SaveVegetatedArea()

#Compute temporally mean value
def processCarbongap_and_NPP_mean():
    npp=ee.Image(asset_location+"/NPP"+str(2001))
    gap=ee.Image(asset_location+"/Carbongap"+str(2001))
    for i in range(17): #2000-2018
       year=str(2002+i)
       npp=npp.add(ee.Image(asset_location+"/NPP"+str(year)))
       gap=gap.add(ee.Image(asset_location+"/Carbongap"+str(year)))

    task1 = ee.batch.Export.image.toAsset(
                  image=npp.divide(ee.Number(18)),
                  description='CarbongapMean',
                  assetId=asset_location+'/CarbongapMean',
                  scale=scaleSize,
                  maxPixels=1e13,
                  crs='EPSG:4326',
                  region=roi.geometry()
    )
    task1.start()
    task2 = ee.batch.Export.image.toAsset(
                  image=gap.divide(ee.Number(18)),
                  description='CarbongapMean',
                  assetId=asset_location+'/NPPMean',
                  scale=scaleSize,
                  maxPixels=1e13,
                  crs='EPSG:4326',
                  region=roi.geometry()
    )
    task2.start()
#processCarbongap_and_NPP_mean()

#compute statistis for vegetated area, mean carbon gap,  mean npp, total carbon gap, and total npp for each continent/region for each year
def processStatistics(year):
    region_list=[]
    region_list.append("North America")
    region_list.append("Europe")
    region_list.append("Africa")
    region_list.append("Australia")
    region_list.append("E Asia")
    region_list.append("N Asia")
    region_list.append("S Asia")
    region_list.append("SE Asia")
    region_list.append("SW Asia")
    region_list.append("Central Asia")
    region_list.append("Central America")
    region_list.append("South America")

    for j in range(len(region_list)):
        a=countries2.filterMetadata("wld_rgn","contains",region_list[j])
        if region_list[j]=="North America":
            a=a.filterMetadata("country_na","not_equals","Greenland")

        f1 = open(file_location+"/continent_landcoverarea.csv", "a")
        f_gap = open(file_location+"/continent_meangap.csv", "a")
        f_npp = open(file_location+"/continent_meannpp.csv", "a")
        f_npptotal = open(file_location+"/continent_totalnpp.csv", "a")
        f_gaptotal = open(file_location+"/continent_totalgap.csv", "a")

        # make statistics for each land (veg) cover types within each continent/region (a)
        mask=ee.Image(asset_location+"/landcoverMOD12Q1_"+str(year)).clip(a)

        area=ee.Image(asset_location+"/pixelArea").clip(a)
        gap= ee.Image(asset_location+"/Carbongap"+str(year)).clip(a)
        npp= ee.Image(asset_location+"/NPP"+str(year)).clip(a)

        gap_total = gap.mask(mask).addBands(mask).reduceRegion(
              reducer=ee.Reducer.sum().group(
                groupField=1,
                groupName='code',
              ),
              geometry=a,
              maxPixels=1e13,
              bestEffort=True
        )

        npp_total = npp.mask(mask).addBands(mask).reduceRegion(
              reducer=ee.Reducer.sum().group(
                groupField=1,
                groupName='code',
              ),
              geometry=a, #rectangle,
              maxPixels=1e13,
              bestEffort=True
        )
        #gap_total
        label=region_list[j]+","+str(year)
        b=gap_total.getInfo()["groups"]
        code_sum=[]
        for m in range(17):
            code_sum.append(0) # set initial value 0 for all

        for k in range(len(b)):
            code_sum[b[k]["code"]-1]=b[k]["sum"]

        for i in range(17):
            label+=","+str(code_sum[i]) #set label
            i += 1

        f_gaptotal.write(label+"\n")
        f_gaptotal.close()

       #npp_total
        label=region_list[j]+","+str(year)
        b=npp_total.getInfo()["groups"]
        code_sum=[]
        for m in range(17):
            code_sum.append(0)

        for k in range(len(b)):
            code_sum[b[k]["code"]-1]=b[k]["sum"]

        for i in range(17):
            label+=","+str(code_sum[i])
            i += 1

        f_npptotal.write(label+"\n")
        f_npptotal.close()

        gap_mean = gap.mask(mask).addBands(mask).reduceRegion(
              reducer=ee.Reducer.mean().group(
                groupField=1,
                groupName='code',
              ),
              geometry=a, #rectangle,
              maxPixels=1e13,
              bestEffort=True
        )
        npp_mean = npp.mask(mask).addBands(mask).reduceRegion(
              reducer=ee.Reducer.mean().group(
                groupField=1,
                groupName='code',
              ),
              geometry=a, #rectangle,
              maxPixels=1e13,
              bestEffort=True
        )

        area_f_veg = area.mask(mask).addBands(mask).reduceRegion(
              reducer=ee.Reducer.sum().group(
                groupField=1,
                groupName='code',
              ),
              geometry=a, #rectangle,
              maxPixels=1e13,
              bestEffort=True
        )

        #area
        label=region_list[j]+","+str(year)
        b=area_f_veg.getInfo()["groups"]

        code_sum=[]
        for m in range(17):
            code_sum.append(0)

        for k in range(len(b)):
            code_sum[b[k]["code"]-1]=b[k]["sum"]

        for i in range(17):
            label+=","+str(code_sum[i])
            i += 1

        f1.write(label+"\n")
        f1.close()

        #gap mean
        label=region_list[j]+","+str(year)
        b=gap_mean.getInfo()["groups"]
        code_mean=[]
        for m in range(17):
            code_mean.append(0)

        for k in range(len(b)):
            code_mean[b[k]["code"]-1]=b[k]["mean"]

        for i in range(17):
            label+=","+str(code_mean[i])
            i += 1

        f_gap.write(label+"\n")
        f_gap.close()

        #npp mean
        label=region_list[j]+","+str(year)
        b=npp_mean.getInfo()["groups"]
        code_mean=[]
        for m in range(17):
            code_mean.append(0)

        for k in range(len(b)):
            code_mean[b[k]["code"]-1]=b[k]["mean"]

        for i in range(17):
            label+=","+str(code_mean[i])
            i += 1

        f_npp.write(label+"\n")
        f_npp.close()
    print('finished')

def processStatisticsBatch():
    for i in range(18):
        processStatistics(str(2018-i))
#processStatisticsBatch()


# statistics on carbon gap and population density
# geo: contiment feature, label: contiment/region name
def computeStatistics_Carbongap_NPP_Population_Percentile(geo,label):
    landcover= ee.Image(asset_location+"/Vegetated_AllYears").clip(geo)
    gap= ee.Image(asset_location+"/CarbongapMean").clip(geo)
    npp= ee.Image(asset_location+"/NPPMean").clip(geo)
    area=ee.Image(asset_location+"/pixelArea").clip(geo)
    popu = ee.ImageCollection("CIESIN/GPWv411/GPW_UNWPP-Adjusted_Population_Count").first().clip(geo)
    maskmap=gap.gt(ee.Number(0)).clip(geo).mask(landcover)

    f = open(file_location+"/percentilestat.csv", "a")
    p_list=[5,10, 15,20, 25,30, 35,40, 45,50, 55,60,65, 70, 75,80,85, 90,95,100] # percentile divisions
    p=ee.Reducer.percentile(p_list)
    string_list_gap=label+",Carbongap_Mean"
    string_list_gap_pct_threshold=label+",Carbongap_Max"
    string_list_area=label+",Area_this_PCT"
    string_list_popu=label+",Popu_Mean"
    string_list_npp=label+",NPP_Mean"
    threshold = gap.mask(maskmap).reduceRegion(
          reducer=p,
          geometry=geo,
          maxPixels=1e13,
          bestEffort=True
    )
    v=threshold.getInfo()
    #list.append(v)
    previous_break=0
    for indx in range(len(p_list)):
        t="mean_p"+str(p_list[indx])
        mask=gap.gt(ee.Number(previous_break)).And(gap.lte(v[t]))
        previous_break=v[t]
        string_list_gap_pct_threshold+=","+str(v[t])
        area_f = area.mask(mask).reduceRegion(
              reducer=ee.Reducer.sum(),
              geometry=geo,
              maxPixels=1e13,
              bestEffort=True
        )
        popu_f = popu.mask(mask).reduceRegion(
              reducer=ee.Reducer.mean(),
              geometry=geo,
              maxPixels=1e13,
              bestEffort=True
        )
        gap_mean= gap.mask(mask).reduceRegion(
              reducer=ee.Reducer.mean(),
              geometry=geo,
              maxPixels=1e13,
              bestEffort=True
        )
        npp_mean= npp.mask(mask).reduceRegion(
              reducer=ee.Reducer.mean(),
              geometry=geo,
              maxPixels=1e13,
              bestEffort=True
        )
        string_list_area+=","+str((area_f.getInfo())["area"])
        string_list_popu+=","+str((popu_f.getInfo())["unwpp-adjusted_population_count"])
        string_list_gap+=","+str((gap_mean.getInfo())["mean"])
        string_list_npp+=","+str((npp_mean.getInfo())["mean"])

    f.write(string_list_gap_pct_threshold+"\n")
    f.write(string_list_area+"\n")
    f.write(string_list_popu+"\n")
    f.write(string_list_gap+"\n")
    f.write(string_list_npp+"\n")
    f.close()
    print('finished')

# statistics on carbon gap and population density for all continents/regions
def processAllStatistics_Carbongap_NPP_Population_Percentile():
    region_list=[]
    region_list.append("North America")
    region_list.append("Europe")
    region_list.append("Africa")
    region_list.append("Australia")
    region_list.append("E Asia")
    region_list.append("N Asia")
    region_list.append("S Asia")
    region_list.append("SE Asia")
    region_list.append("SW Asia")
    region_list.append("Central Asia")

    region_list.append("Central America")
    region_list.append("South America")

    for j in range(len(region_list)):
        a=countries2.filterMetadata("wld_rgn","contains",region_list[j])
        if region_list[j]=="North America":
            a=a.filterMetadata("country_na","not_equals","Greenland")
        computeStatistics_Carbongap_NPP_Population_Percentile(a,region_list[j])
#processAllStatistics_Carbongap_NPP_Population_Percentile()

#--------------------------------------End computing statistics------------------------------------------
