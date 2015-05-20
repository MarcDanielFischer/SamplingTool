
package com.arbonaut.sampling;

import java.io.File;
import java.util.Collection;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;
import java.util.ArrayList;

import com.arbonaut.sampling.design.SamplingDesign;
import com.arbonaut.sampling.SamplePlot;
import com.arbonaut.sampling.io.SampleWriter;
import com.arbonaut.sampling.util.UTMFinder;

import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.filter.Filter;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.geometry.jts.JTS;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.filter.text.cql2.CQL;
import org.geotools.referencing.CRS;
import org.opengis.referencing.operation.MathTransform;


import com.arbonaut.sampling.io.CsvSampleWriter; 


/**
 * Provides a set of handy methods for carrying out spatial sampling using this
 * library. 
 * 
 * Please note that this is supposed to serve as an example rather than a fully
 * functional production code. Different users of this library may have 
 * different needs (such as the CRS in which sampling is done etc...) so it is
 * advised to roll out your own implementation if this one is not enough.            
 */
public final class Sampling {
	
    private Sampling() { }
	
    
    /**
     * Finds all unique stratum names in the shapefile.
     * @param shapeFile         Shapefile to read
     * @param stratumColumn     Name of the column that holds the stratum name       
     */         
    public static Set<String> findUniqueStrata(File shapeFile, 
                                               String stratumColumn) throws Exception
    {
        if (!shapeFile.exists())
            throw new Exception(shapeFile.toString() + " does not exist");

        FileDataStore dataStore = FileDataStoreFinder.getDataStore(shapeFile);
        if (dataStore == null)
            throw new Exception("Cannot find data store for " + shapeFile.toString());
        
        SimpleFeatureSource featureSource = dataStore.getFeatureSource();
        return findUniqueStrata(featureSource.getFeatures(), stratumColumn);
    }
    
    
    /**
     * Finds all unique stratum names in the feature collection.
     * @param features          Feature collection to scan 
     * @param stratumColumn     Name of the column that holds the stratum name       
     */         
    public static Set<String> findUniqueStrata(SimpleFeatureCollection features, 
                                               String stratumColumn) throws Exception
    {
        // Is there a better way than just going through all?
        Set<String> strata = new HashSet<String>();
        SimpleFeatureIterator iterator = features.features();
		try {
			while( iterator.hasNext()) {
				SimpleFeature f = iterator.next();
                strata.add(f.getAttribute(stratumColumn).toString());
			}
		}
		finally {
			iterator.close(); 
		}
        return strata;        
    }
    
    
    /**
     * Generates spatial sample over a census defined by features in the 
     * shapefile without stratification.
     * 
     * @param shapeFile Shapefile defining the census to be sampled
     * @param sampleSize     The amount of plots to place
     * @param forceUTM       If true, sampling is done in a suitable UTM zone
     *                       instead of the native CRS of the dataset                                                                
     * @param sampler        Sample design to use
     * @param writer         Sample writer to use     
     */              
    public static void generateSample(File shapeFile, int sampleSize,
                                      boolean forceUTM,
                                      SamplingDesign sampler, 
                                      SampleWriter writer) throws Exception 
    {
        Map<String, Integer> N = new HashMap<String, Integer>();
        N.put("", sampleSize);
        generateSample(shapeFile, null, N, forceUTM, sampler, writer);
    }
    

    /**
     * Generates spatial sample over a census defined by features in the 
     * collection without stratification.
     * 
     * @param features       Feature collection defining the census to be sampled
     * @param sampleSize     The amount of plots to place                                                      
     * @param forceUTM       If true, sampling is done in a suitable UTM zone
     *                       instead of the native CRS of the dataset                                                                
     * @param sampler        Sample design to use
     * @param writer         Sample writer to use     
     */              
    public static void generateSample(SimpleFeatureSource features, int sampleSize,
                                      boolean forceUTM,
                                      SamplingDesign sampler, SampleWriter writer) throws Exception 
    {
        Map<String, Integer> N = new HashMap<String, Integer>();
        N.put("", sampleSize);
        generateSample(features, null, N, forceUTM, sampler, writer);
    }
    
    
    /**
     * Generates spatial sample over a census defined by features in the 
     * shapefile. The sample can be drawn per-stratum or without taking strata
     * into account.
     * 
     * @param shapeFile Shapefile defining the census to be sampled
     * @param stratumColumn  Name of a column holding stratum names. Can be null 
     *                       if stratification should not be used
     * @param sampleSize     The amount of plots to place per stratum. If 
     *                       stratification is not used, the map should contain 
     *                       a single item (key does not matter) with the number
     *                       of plots to be placed. Otherwise the key should be
     *                       the stratum name                                                                  
     * @param forceUTM       If true, sampling is done in a suitable UTM zone
     *                       instead of the native CRS of the dataset                                                                
     * @param sampler        Sample design to use
     * @param writer         Sample writer to use     
     */              
    public static void generateSample(File shapeFile, String stratumColumn, 
                                      Map<String, Integer> sampleSize, boolean forceUTM,
                                      SamplingDesign sampler, SampleWriter writer) throws Exception 
    {
        if (!shapeFile.exists())
            throw new Exception(shapeFile.toString() + " does not exist");
        
        FileDataStore dataStore = FileDataStoreFinder.getDataStore(shapeFile);
        if (dataStore == null)
            throw new Exception("Cannot find data store for " + shapeFile.toString());
        
		SimpleFeatureSource featureSource = dataStore.getFeatureSource();
        generateSample(featureSource, stratumColumn, sampleSize, forceUTM, sampler, writer);
    }
    

    /**
     * Generates spatial sample over a census defined by features in the 
     * collection. The sample can be drawn per-stratum or without taking strata
     * into account.
     * 
     * @param features       Feature collection defining the census to be sampled
     * @param stratumColumn  Name of a column holding stratum names. Can be null 
     *                       if stratification should not be used
     * @param sampleSize     The amount of plots to place per stratum. If 
     *                       stratification is not used, the map should contain 
     *                       a single item (key does not matter) with the number
     *                       of plots to be placed. Otherwise the key should be
     *                       the stratum name                                                                  
     * @param forceUTM       If true, sampling is done in a suitable UTM zone
     *                       instead of the native CRS of the dataset                                                                
     * @param sampler        Sample design to use
     * @param writer         Sample writer to use     
     */              
    public static void generateSample(SimpleFeatureSource features, String stratumColumn, 
                                      Map<String, Integer> sampleSize, boolean forceUTM,
                                      SamplingDesign sampler, SampleWriter writer) throws Exception 
    {
        // Don't have stratum column, just use all features.
        if (stratumColumn == null)
        {
            int size = sampleSize.values().iterator().next();
            generateSampleForStratum(features.getFeatures(), size, null, forceUTM, sampler, writer);            
        }
        // Have stratum columns so generate distinct sample for each stratum
        else
        {
            Set<String> strata = findUniqueStrata(features.getFeatures(), stratumColumn);
            for (String stratum : strata)
            {
                Integer size = sampleSize.get(stratum);
                if (size == null)
                    throw new Exception("No size given for stratum " + stratum);    
            
                // Get features that belong to this stratum from the feature 
                // collection 
                
                Filter filter = CQL.toFilter(
                    String.format("%s='%s'", stratumColumn, stratum));
                      
                generateSampleForStratum(
                    features.getFeatures(filter), 
                    size, stratum, forceUTM, sampler, writer);   
            }
        }    
    
    
    }
  
  
    /*
     * Generates sample taking all features in the collection as belonging to 
     * the same stratum (or no stratum if stratumName is null).   
     */     
    private static void generateSampleForStratum(SimpleFeatureCollection stratumFeatures,
                                                 int sampleSize, String stratumName, boolean forceUTM,
                                                 SamplingDesign sampler, SampleWriter writer) throws Exception 
    {
        // We'll go over all the features and union them together. This is best 
        // done by building a GeometryCollection and the calling its union() 
        // method, which unions all geoms at once (requires JTS 1.9 or newer).
         
        Collection<Geometry> allGeoms = new ArrayList<Geometry>();
        SimpleFeatureIterator iterator = stratumFeatures.features();
        
        CoordinateReferenceSystem crs = null;
		try 
        {
			while( iterator.hasNext()) 
            {
				SimpleFeature f = iterator.next();
                Geometry geom = (Geometry)f.getDefaultGeometry();
                
                // Assume all features have the same CRS
                if (crs == null) 
                {
                    crs = f.getType().getCoordinateReferenceSystem();
                    if (crs == null)
                        throw new BadCRSException("Feature has no coordinate system");
                }
                allGeoms.add(geom); 
			}
		}
		finally {
			iterator.close(); // prevents memory leaks and data loss
		}
        
        // Since source CRS may be missing some parameters (if coming from Shapefiles), 
        // we try to retrieve the full definition if we can find a EPSG code for it
        
        Integer epsg = CRS.lookupEpsgCode(crs, true);
        if (epsg != null)
        {
            System.out.println("epsg: " + epsg.toString()); 
            crs = CRS.decode("EPSG:"+epsg.toString());
        } 
        
        
        GeometryFactory gf = JTSFactoryFinder.getGeometryFactory( null );
        GeometryCollection geometryCollection =
            (GeometryCollection) gf.buildGeometry( allGeoms );  

        Geometry combinedGeom = geometryCollection.union();
        
        
        // Now that we have the combined geometry, check if the features are in
        // a geographic system while projected system is required by the design.
        // If yes, reproject census to a suitable UTM CRS.
        
        if (forceUTM)
        {
            CoordinateReferenceSystem utmCRS = UTMFinder.findBestUTM(combinedGeom, crs);
            MathTransform xform = CRS.findMathTransform(crs, utmCRS);
    
            combinedGeom = JTS.transform(combinedGeom, xform);
            crs = utmCRS;
            
            // Just for testing the reprojection
            // CsvSampleWriter csv = new CsvSampleWriter(new File("c:\\temp\\censusInUTM.csv"), sampler, false);
            // SamplePlot p = new SamplePlot(combinedGeom, crs, "");
            // csv.write(p);
            // csv.close(); 
            //System.out.println("Using CRS: " + utmCRS.toString()); 
        } 
        
        
        // Now that we have the combined geometry for this stratum, run sampling    
        sampler.generate(combinedGeom, crs, sampleSize, stratumName, writer);
    }

}
