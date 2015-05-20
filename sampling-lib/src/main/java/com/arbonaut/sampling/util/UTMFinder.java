package com.arbonaut.sampling.util;

import java.util.ArrayList;

import org.geotools.feature.simple.SimpleFeatureImpl;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.geotools.referencing.CRS;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.BoundingBox;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.Polygonal;


/**
 * This class provides static methods to locate suitable UTM coordinate system
 * for given geometry.   
 * @author daniel
 */
public final class UTMFinder {

    private UTMFinder() {}

	/**
	 * Finds the best-fitting UTM zone for the input geometry (and its CRS) and 
	 * returns the CoordinateReferenceSystem that describes said UTM zone.
	 * 
	 * If the geometry does not already come in EPSG:4326 CRS (WGS84), this 
	 * method first reprojects the geometry into this particular CRS in order to 
	 * obtain a bounding box that uses lon/lat values which are in turn required 
	 * for the method to properly derive an adequate UTM zone (finding a 
	 * matching UTM zone depends mainly on Longitude values). 
	 * 
	 * If a feature spans more than one UTM zone, the average of all zones is 
	 * used as the target UTM zone. If the average is not an integer number, 
	 * it is rounded down to the next integer value (by convention, it might as 
	 * well be rounded up).
	 * 
	 * Behaviour for cross-equator-geometries is probably not tested sufficiently.
	 * 
	 *           
	 * @param geometry  Geometry for which a fitting UTM zone is to be retrieved
	 * @param crs       Coordinate system of the geometry
	 */
	public static CoordinateReferenceSystem 
    findBestUTM(Geometry geom, CoordinateReferenceSystem crs) throws Exception 
    {
        // forces x/y order, otherwise axis order might be messed up by the JTS.transform() method
        //System.setProperty("org.geotools.referencing.forceXY", "true"); 
        
		CoordinateReferenceSystem standardGeographicCRS = CRS.decode("EPSG:4326");
        
		// check if feature CRS is different from EPSG:4326
		boolean needsReproject = !CRS.equalsIgnoreMetadata(standardGeographicCRS, crs);
		
		// values we need to find a matching UTM zone for the input feature
		
        
        // if feature CRS is different from EPSG:4326
		if (needsReproject) {
        
            System.out.println("towgs "); 
         
			MathTransform transform = CRS.findMathTransform(crs, standardGeographicCRS);
			geom = JTS.transform(geom, transform); // transform feature Geometry to EPSG:4326
			
			// feature Geometry might either be a Polygon or a MultiPolygon object depending on input data
			if (!(geom instanceof Polygonal))
            	throw new Exception("Geometry is neither a Polygon nor a MultiPolygon object");
        }
		
        Envelope envelope = geom.getEnvelopeInternal();
        double lonMin = envelope.getMinX();
		double lonMax = envelope.getMaxX();
	    double latMin = envelope.getMinY();
		double latMax = envelope.getMaxY();
		
        
		int zoneMin = longitude2UTM(lonMin);
		int zoneMax = longitude2UTM(lonMax);
		/*
		 * If a feature spans more than one UTM zone,
		 * the average of all zones is used as the target UTM zone.
		 * If the average is not an integer number, it is rounded down
		 * to the next integer value (by convention, it might as well be rounded up).
		 */
		int utmZone = (int)Math.floor((zoneMin + zoneMax) / 2); 

		/*
		 * derive EPSG code for UTM Zone:
		 * UTM                    EPSG
		 * 01 N                   32601
		 * 02 N                   32602
		 * 60 N                   32660
		 * EPSG codes for UTM zones 1-60 N : 32601 - 32660 (all for WGS84 Datum, there are others, too)
		 * EPSG codes for UTM zones 1-60 S : 32701 - 32760 (all for WGS84 Datum, there are others, too)
		 * --> Formula: if(northernHemisphere) EPSG = 32600 + utmZone
		 * if(southernHemisphere) EPSG = 32700 + utmZone
		 * bei Hemisphärenüberschreitenden features: Lat-Mittelwert ((latMin + latMax) / 2) versuchen,
		 * könnte unerwartetes Verhalten verursachen bei evtl. negativen Hochwerten 
		 * (wenn latMean auf Nordhalbkugel liegt und sich ein feature auch auf die Südhalbkugel erstreckt)
		 * --> Testen
		 * (dann evtl mit UTM-Südzonen probieren, die haben false northing, da gibts keine negativen Werte)
		 */
		
		double latMean = (latMin + latMax)/2;

		/*
		 * The EPSG Code ist determined depending on whether
		 * the mean latitude value is on the northern or the southern hemisphere
		 * (mean latitude value: just the mean between min and max latitude values of the feature 
		 * bounding box; the fact which hemisphere the biggest area proportion 
		 * of the feature is located on is not taken into account here)
		 */
		int epsgCode;
		if(latMean > 0){ 
			epsgCode = 32600 + utmZone; // northern hemisphere
		}else{ 
			epsgCode = 32700 + utmZone; // southern hemisphere
		}

		// create target CRS from EPSG code
		String epsg = "EPSG:" + epsgCode;
		return CRS.decode(epsg);
	}

	
	
    
    /**
	 * Derive UTM Zone for a given longitude value.
	 * @param longitude
	 * @return UTM Zone
	 */
	public static int longitude2UTM(double longitude) {
		// TODO Ausnahmen: Norwegen etc. --> also auch abhängig von Latitude --> evtl behandeln --> start: welche Zonen sind überhaupt irregulär?
		int utmZone = (int)Math.floor(((longitude + 180) / 6) +1);
		if(utmZone > 60) utmZone = utmZone % 60; // if input longitude is > 180 for some reason (and output UTM zone is > 60 then)
		return utmZone;
	}

	

}
