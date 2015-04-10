/*
 * Hier versuche ich jetzt erstmals, auf ein Polygon(Feature) zuzugreifen,
 * das CRS auszulesen und dann zu transformieren
 * 
 * Wofür brauche ich das? um Cluster zu bauen, muss ich mit Meterangaben rechnen ("und jetzt 200m nach Westen.."),
 * und dafür brauche ich ein metrisches CRS (wohl idR UTM, es gibt aber anscheinend auch andere)
 */


package test;

import java.io.File;
import java.util.ArrayList;

import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.BoundingBox;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Geometry;

public class LonLat2UTM_Test {
	
	
	public static int long2UTM(double longitude){
		// TODO Ausnahmen: Norwegen etc. --> also auch abhängig von Latitude --> evtl behandeln --> start: welche Zonen sind überhaupt irregulär?
		int utmZone = (int)Math.floor(((longitude + 180) / 6) +1);
		if(utmZone > 60) utmZone = utmZone % 60; // if input longitude is > 180 for some reason (and output UTM zone is > 60 then)
		return utmZone;
	}
	
	
	public static void main(String [] args) throws Exception{
		
		// Access first SHP Feature (GeoTools logic): Input File --> DataStore --> FeatureSource --> FeatureCollection --> Iterator
		File inputFile = new File("D:\\Daniel\\TestData\\Deutschland\\Deutschland.shp");
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(inputFile);
		SimpleFeatureSource featureSource = dataStore.getFeatureSource(); // wird benötigt, um an einzelne Features ranzukommen (Feature = Zeile in SHP Attribute Table)
		final SimpleFeatureCollection featureCollection = featureSource.getFeatures(); // wieso muss die final sein?		
		SimpleFeatureIterator iterator = featureCollection.features(); // an die Features(=Zeilen in SHP attribute table) kommt man, indem man über die FeatureCollection iteriert 
		// unsecured access to  first feature
		 SimpleFeature firstFeature = iterator.next();
		iterator.close(); // prevents memory leaks or data loss or whatever

			
		// read SHP CRS (sourceCRS)
		CoordinateReferenceSystem sourceCRS = dataStore.getSchema().getCoordinateReferenceSystem();
		
		// find matching UTM zone for feature
		BoundingBox boundingBox = firstFeature.getBounds();
		double lonMin = boundingBox.getMinX();
		double lonMax = boundingBox.getMaxX();
		
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
		
		double latMin = boundingBox.getMinY();
		double latMax = boundingBox.getMaxY();
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
		
	
		// find MathTransform
		// brauche ich dafür den EPSG-Code oder gehts auch direkt?? 
		// --> geht auch direkt, kommen aber andere Sachen raus 
		// transform source CRS directly (without looking up the corresponding EPSG code first and use that CRS instead)
		MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
		
		// transform source CRS after having looked up the corresponding EPSG code
		String lookedUpEPSGcodeForSourceCRS = CRS.lookupIdentifier( sourceCRS, true );
		CoordinateReferenceSystem formalizedSourceCRS = CRS.decode(lookedUpEPSGcodeForSourceCRS);
		MathTransform transform2 = CRS.findMathTransform(formalizedSourceCRS, targetCRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
		
		// get feature Geometry
		Geometry sourceGeometry = (Geometry) firstFeature.getAttribute( "the_geom" );
//		Geometry geom2 = (Geometry)firstFeature.getDefaultGeometry(); // so gehts auch, 2 Wege

		
		// convert Geometry to target CRS
		Geometry targetGeometry = JTS.transform( sourceGeometry, transform);
		
		// sample points inside targetGeometry
		
		
		/*
		 *  und jetzt? 
		 *  
		 *  So wie ich es verstanden habe, muss man von jedem Feature einzeln die Geometry transformieren,
		 *  also muss ich diesen Transform-Schritt in einer Schleife machen; und wie speichert man dann?
		 *  
		 *  ich habe jetzt nur die Geometrien in der Hand, losgelöst vom restlichen Feature-Inhalt (dh, auch keine Namen)
		 *  wenn ich jetzt nur mit den Geometrien weiterbastele, sind die Namen weg und so
		 *  --> ich brauche ja auch gar keine Namen, ich will nur Punkte Sampeln innerhalb der jeweils ausgewählen Geometry
		 *  
		 *  
		 */
		
		
		
		
	}
}
