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
	 * Gets a matching UTM CoordinateReferenceSystem for the input SimpleFeature.
	 * If the feature spans more than one UTM zone, the formula (zoneMin + zoneMax) / 2 is used.
	 * Behaviour for cross-equator-features is not tested yet. 
	 * @param feature
	 * @return
	 */
	public static CoordinateReferenceSystem getTargetUTM_CRS(SimpleFeature feature, CoordinateReferenceSystem featureCRS) throws Exception{
		// find matching UTM zone for feature (UTM zones depend mainly on Longitude values)
		

		
		// boundingBox must have latLong values in order to derive an UTM zone
		// that makes sense
		// --> convert BBOX to WGS84 CRS if needed
		
		// how to find Default Geographic CRS in Geotools? --> see tutorial
		// CRS.decode("EPSG:4326");
		// DefaultGeographicCRS.WGS84;`
		
		// check if CRS from Plot and Raster are different and reproject if needed
		System.setProperty("org.geotools.referencing.forceXY", "true"); // sonst werden evtl X und Y vertauscht
		CoordinateReferenceSystem standardGeographicCRS = CRS.decode("EPSG:4326");
		boolean needsReproject = !CRS.equalsIgnoreMetadata(standardGeographicCRS, featureCRS);
		
		double lonMin;
		double lonMax;
		double latMin;
		double latMax;
		
		if(needsReproject){
			MathTransform transform = CRS.findMathTransform(featureCRS, standardGeographicCRS, true);
			Geometry geom1 = (Geometry)feature.getAttribute( "the_geom" );
			Geometry geom2 = JTS.transform(geom1, transform);
			// feature.setDefaultGeometry(geom2); // schlecht --> permanent change, auch bei copied feature
			// wie komme ich an minX etc. über Geometry?
			Envelope envelope = null;
			if(geom2 instanceof Polygon){
				geom2 = (Polygon) geom2;
				envelope = geom2.getEnvelopeInternal();
			}else if(geom2 instanceof MultiPolygon){
				geom2 = (MultiPolygon) geom2;
				envelope = geom2.getEnvelopeInternal();
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
		int utmZone = (int)Math.floor((zoneMin + zoneMax) / 2); // for when the feature spans more than one UTM zone
		/*
		 *  von Katja übernommen, ist Konventionssache: 
		 *  wenn ein stratum sich über mehrere UTM-Zonen erstreckt, 
		 *  wird die mittlere UTM-Zone für das gesamt Stratum ausgewählt (bei ungerader Anzahl von UTM-Zonen) 
		 *  bzw. abgerundet bei gerade Anzahl von UTM-Zonen (man könnte zb auch aufrunden)
		 *  --> gibts auch metrische CRS für größere Ausdehnungen (UTM-Zonen sind immer nur 6 Grad breit)??
		 */


		// derive EPSG code for UTM Zone:
		//  UTM                    EPSG
		//  01 N                   32601
		//  02 N                   32602
		//  60 N                   32660
		// EPSG codes for UTM zones 1-60 N : 32601 - 32660 (all for WGS84 Datum, there are others, too)
		// EPSG codes for UTM zones 1-60 S : 32701 - 32760 (all for WGS84 Datum, there are others, too)
		// --> Formel: if(northernHemisphere) EPSG = 32600 + utmZone
		// if(southernHemisphere) EPSG = 32700 + utmZone
		// bei Hemisphärenüberschreitenden features: Lat-Mittelwert ((latMin + latMax) / 2) versuchen,
		// könnte unerwartetes Verhalten verursachen bei evtl. negativen Hochwerten 
		// (wenn latMean auf Nordhalbkugel liegt und sich ein feature auch auf die Südhalbkugel erstreckt)
		// --> Testen
		// (dann evtl mit UTM-Südzonen probieren, die haben false northing, da gibts keine negativen Werte)

		
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
}
