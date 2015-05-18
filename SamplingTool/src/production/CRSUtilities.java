package production;

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
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;


/**
 * This class provides static methods to perform CRS operations such as conversion 
 * between UTM CRS and geographic CRS.
 * @author daniel
 *
 */
public class CRSUtilities {

	/**
	 * Finds the best-fitting UTM zone for the input SimpleFeature and returns the CoordinateReferenceSystem
	 * that describes said UTM zone.
	 * If the SimpleFeature does not already come in EPSG:4326 CRS (standard geographic projection using lon/lat values),
	 * this method first reprojects the SimpleFeature into this particular CRS in order to obtain a
	 * bounding box that uses lon/lat values which are in turn required for the method to properly derive an 
	 * adequate UTM zone (finding a matching UTM zone depends mainly on Longitude values). 
	 * 
	 *
	 * If a feature spans more than one UTM zone,
	 * the average of all zones is used as the target UTM zone.
	 * If the average is not an integer number, it is rounded down
	 * to the next integer value (by convention, it might as well be rounded up).
	 * 
	 * 
	 * Behaviour for cross-equator-features is probably not tested sufficiently. 
	 * @param feature the feature to find a fitting UTM zone for.
	 * @param featureCRS required additional information as the feature itself does not contain the CRS it uses.
	 * @return
	 */
	public static CoordinateReferenceSystem getTargetUTM_CRS(SimpleFeature feature, CoordinateReferenceSystem featureCRS) throws Exception{
		System.setProperty("org.geotools.referencing.forceXY", "true"); // forces x/y order, otherwise axis order might be messed up by the JTS.transform() method
		CoordinateReferenceSystem standardGeographicCRS = CRS.decode("EPSG:4326");
		// check if feature CRS is different from EPSG:4326
		boolean needsReproject = !CRS.equalsIgnoreMetadata(standardGeographicCRS, featureCRS);
		
		// values we need to find a matching UTM zone for the input feature
		double lonMin;
		double lonMax;
		double latMin;
		double latMax;
		
		if(needsReproject){ // if feature CRS is different from EPSG:4326
			MathTransform transform = CRS.findMathTransform(featureCRS, standardGeographicCRS, true);
			Geometry geom = (Geometry)feature.getAttribute( "the_geom" );
			geom = JTS.transform(geom, transform); // transform feature Geometry to EPSG:4326
			Envelope envelope = null;
			// feature Geometry might either be a Polygon or a MultiPolygon object depending on input data
			if(geom instanceof Polygon){
				geom = (Polygon) geom;
				envelope = geom.getEnvelopeInternal(); // getEnvelopeInternal() returns the bounding box of the Geometry
			}else if(geom instanceof MultiPolygon){
				geom = (MultiPolygon) geom;
				envelope = geom.getEnvelopeInternal();
			}else{
				throw new Exception("Feature Geometry is neither a Polygon nor a MultiPolygon object");
			}
			lonMin = envelope.getMinX();
			lonMax = envelope.getMaxX();
			latMin = envelope.getMinY();
			latMax = envelope.getMaxY();
		}else{ // if there is no need for reprojection (feature CRS == EPSG:4326)
			BoundingBox boundingBox = feature.getBounds(); // runtime class: ReferencedEnvelope
			
			lonMin = boundingBox.getMinX();
			lonMax = boundingBox.getMaxX();
			latMin = boundingBox.getMinY();
			latMax = boundingBox.getMaxY();
		}

		int zoneMin = long2UTM(lonMin);
		int zoneMax = long2UTM(lonMax);
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
		String crs = "EPSG:" + epsgCode;
		CoordinateReferenceSystem targetCRS = CRS.decode(crs);
		return targetCRS;
	}

	
	
	
	
	/**
	 * TODO MEthode von oben umschreiben für Stratum
	 * Soll ich dann die beiden Varianten drinlassen? Man weiß ja nie...
	 * @return
	 */
	public static CoordinateReferenceSystem getTargetUTM_CRS(Stratum stratum)throws Exception{

		// reproject Geometry to EPSG:4326 if needed
		System.setProperty("org.geotools.referencing.forceXY", "true"); // forces x/y order, otherwise axis order might be messed up by the JTS.transform() method
		CoordinateReferenceSystem stratumCRS = stratum.getCRS();
		CoordinateReferenceSystem standardGeographicCRS = CRS.decode("EPSG:4326");
		// check if feature CRS is different from EPSG:4326
		boolean needsReproject = !CRS.equalsIgnoreMetadata(standardGeographicCRS, stratumCRS);
		Geometry geom = stratum.getGeometry();
		if(needsReproject){ // if feature CRS is different from EPSG:4326
			MathTransform transform = CRS.findMathTransform(stratumCRS, standardGeographicCRS, true);
			geom = JTS.transform(geom, transform); // transform feature Geometry to EPSG:4326
		}
		
		// get Bounding Box and extract Min/Max values
		Envelope envelope = geom.getEnvelopeInternal(); // getEnvelopeInternal() returns the bounding box of the Geometry
		double lonMin = envelope.getMinX();
		double lonMax = envelope.getMaxX();
		double latMin = envelope.getMinY();
		double latMax = envelope.getMaxY();

		// find matching UTM zone 
		int zoneMin = long2UTM(lonMin);
		int zoneMax = long2UTM(lonMax);
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
		String crs = "EPSG:" + epsgCode;
		CoordinateReferenceSystem targetCRS = CRS.decode(crs);
		return targetCRS;
	}
	
	
	
	
	
	/**
	 * Derive UTM Zone for a given longitude value.
	 * @param longitude
	 * @return UTM Zone
	 */
	public static int long2UTM(double longitude){
		// TODO Ausnahmen: Norwegen etc. --> also auch abhängig von Latitude --> evtl behandeln --> start: welche Zonen sind überhaupt irregulär?
		int utmZone = (int)Math.floor(((longitude + 180) / 6) +1);
		if(utmZone > 60) utmZone = utmZone % 60; // if input longitude is > 180 for some reason (and output UTM zone is > 60 then)
		return utmZone;
	}

	/**
	 * Convenience method to convert processed data from UTM projection to geographic LonLat projection  
	 */
	public static void UTM2LonLat(ArrayList<Plot> plots) throws Exception{

		// iterate over plots and convert their Point property to geographic LonLat
		for(Plot plot : plots){

			// get the plot's current CRS 
			CoordinateReferenceSystem sourceCRS = plot.getCRS();

			// get Point property (UTM projection at this stage)
			Point plotPoint = plot.getPoint();

			// convert Point property to LonLat
			Point plotLatLon = (Point)JTS.toGeographic(plotPoint, sourceCRS);

			// save changed Point property
			plot.setPoint(plotLatLon);

			// change plot's CRS property (just to work properly, this property should not be needed any more at this stage of the process)
			plot.setCRS(org.geotools.referencing.crs.DefaultGeographicCRS.WGS84);

		}

	}
	
	
	
	/**
	 * Reproject Strata CRS to UTM.
	 * @param inputStrata
	 * @return
	 */
	public static Stratum[] convert2UTM(Stratum[] inputStrata)throws Exception{
		Stratum[] outputStrata = new Stratum[inputStrata.length];
		for (int i = 0; i < inputStrata.length; i++) {
			outputStrata[i] = convert2UTM(inputStrata[i]);
		}
		return outputStrata;
	}

	/**
	 * Reproject Stratum CRS to UTM.
	 * @param inputStratum
	 * @return
	 */
	public static Stratum convert2UTM(Stratum inputStratum)throws Exception{
		CoordinateReferenceSystem sourceCRS = inputStratum.getCRS();
		CoordinateReferenceSystem targetUTM_CRS = CRSUtilities.getTargetUTM_CRS(inputStratum); // targetCRS:in UTM projection
		// transform source CRS directly (without looking up the corresponding EPSG code first and use that CRS instead)
		MathTransform transform = CRS.findMathTransform(sourceCRS, targetUTM_CRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
		// get feature Geometry
		Geometry sourceGeometry = inputStratum.getGeometry();
		// convert Geometry to target CRS
		Geometry utmGeometry = JTS.transform( sourceGeometry, transform);

		Stratum outputStratum = new Stratum(utmGeometry, targetUTM_CRS, inputStratum.getName());
		return outputStratum;
	}

	

}
