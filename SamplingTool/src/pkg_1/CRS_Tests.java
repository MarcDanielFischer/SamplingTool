/*
 *  Coordinate Reference Systems CRS
 *  
 *  GeoTools-Kapitel "Referencing" ist hier interessant
 *  
 *  
 *  Was will ich hier?
 *  Um die Cluster fürs Sampling bauen zu können, muss ich das CRS der Inputdaten 
 *  (bis jetzt immer GeoLatLong mit WGS84, kann aber auch was anderes sein) transformieren nach UTM.
 *  Besondere Schwierigkeit: UTM ist nicht weltweit gleich (dh, ich kann nicht die eine Standardtransformation nehmen)
 *  sondern die UTM-Zone ist abhängig von Longitude.
 *  Was muss ich also machen? 
 *  1) find out source data CRS --> dataStore.getSchema().getCoordinateReferenceSystem()
 *  2) find out target UTM zone depending on lon values of input data
 *  3) transform source CRS to target CRS 
 *  
 *  - Hier mache ich ein paar Beispiele durch von GeoTools und schaue, was mir hilft (http://docs.geotools.org/stable/userguide/library/referencing/crs.html)
 *  
 *  
 *  
 *  - Frage: wie findet man den zugehörigen EPSG-Code zu einem CRS, wenn es einen gibt? --> Abschnitt "Matching a CoordinateReferenceSystem" könnte die Lösung sein
 *  - Frage: wie liest man das CRS von einem SHP aus (aus dem zugehörigen .prj file)?
 */

package pkg_1;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Map;
import java.util.Random;

import javax.swing.*;

import org.geotools.metadata.iso.citation.Citations;
import org.geotools.referencing.CRS;
import org.geotools.referencing.ReferencingFactoryFinder;
import org.geotools.referencing.operation.DefiningConversion;
import org.geotools.referencing.wkt.Formattable;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.crs.CRSFactory;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.crs.GeographicCRS;
import org.opengis.referencing.crs.ProjectedCRS;
import org.opengis.referencing.cs.CartesianCS;
import org.opengis.referencing.operation.Conversion;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.MathTransformFactory;

public class CRS_Tests {
	
	
	public static void main(String[] args){
		/*
		 * http://docs.geotools.org/stable/userguide/library/referencing/crs.html:
		 * 
		 *   - CoordinateReferenceSystem can also be defined by a text format 
		 *    ((called “Well Known Text” or WKT). 
		 *    This is a standard provided by the OGC and shows up in inside a 
		 *    shapefile “prj” file, or in a databases such as PostGIS and Oracle.
		 *    To parse WKT please use the CRS.parseWKT( txt ) method
		 *  
		 *  - So sieht das .prj file von den Ghana-Testdaten aus:
		 *    GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]
		 *  
		 */
			 
		
		
		try {
			
			
			
//			// vorhandenes CRS anzeigen
//			CoordinateReferenceSystem sourceCRS = CRS.decode("EPSG:4326");
//			System.out.println(sourceCRS.toString());
			
			
			// eigenes CRS aus WKT bauen
//			String wkt = "GEOGCS[" + "\"WGS 84\"," + "  DATUM[" + "    \"WGS_1984\","
//			        + "    SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],"
//			        + "    TOWGS84[0,0,0,0,0,0,0]," + "    AUTHORITY[\"EPSG\",\"6326\"]],"
//			        + "  PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],"
//			        + "  UNIT[\"DMSH\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9108\"]],"
//			        + "  AXIS[\"Lat\",NORTH]," + "  AXIS[\"Long\",EAST],"
//			        + "  AUTHORITY[\"EPSG\",\"4326\"]]";
//			CoordinateReferenceSystem crs = CRS.parseWKT(wkt);
//			System.out.println(crs.toString());
			
			
//			// ich kann aus den .prj files von meinen Test-SHPS auch eigene CRS erzeugen:
//			String crsGhanaSHP = "GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]";
//			String crsBcBorderSHP = "GEOGCS[\"GCS_North_American_1983\",DATUM[\"North_American_Datum_1983\",SPHEROID[\"GRS_1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]]";
//			CoordinateReferenceSystem ghanaCRS = CRS.parseWKT(crsGhanaSHP);
//			CoordinateReferenceSystem bcBorderCRS = CRS.parseWKT(crsBcBorderSHP);
//			System.out.println(ghanaCRS.toString());
//			System.out.println(bcBorderCRS.toString());
			
			
			/*
			 * Matching a CoordinateReferenceSystem
			 * You can actually search based on any metadata, not just name, the way you do it is 
			 * you construct an example of what you are looking for - and than ask for the best match.
			 * This functionality is especially useful when you have produced a CoordinateReferenceSystem 
			 * by parsing WKT and you would like to find the “official” code for it.
			 * 
			 * --> so bekommt man den EPSG-Code zu einem gegebenen CRS raus: 
			 */
//			String wkt =
//					"GEOGCS[\"ED50\",\n" +
//							"  DATUM[\"European Datum 1950\",\n" +
//							"  SPHEROID[\"International 1924\", 6378388.0, 297.0]],\n" +
//							"PRIMEM[\"Greenwich\", 0.0],\n" +
//							"UNIT[\"degree\", 0.017453292519943295]]";
//			CoordinateReferenceSystem example = CRS.parseWKT(wkt);
//
//			String code = CRS.lookupIdentifier( example, true ); // should be "EPSG:4230"
//			CoordinateReferenceSystem crs = CRS.decode( code ); // "gescheit" umwandeln in ursprünglich gegebenes CRS (?)
//			System.out.println( crs);
			
			
			/*
			 * Finding a Math Transform
			 * 
			 * - When using a CoordinateReferenceSystem that has been parsed from WKT you will often need 
			 *   to “relax” the accuracy by setting the lenient parameter to true when searching with 
			 *   findMathTransform. The official CoordinateReferenceSystem definitions 
			 *   provided by the EPSG database have extra metadata (describing how to do Datum 
			 *   shifts for example), beyond what can be provided using WKT.
			 * 
			 */
//			CoordinateReferenceSystem sourceCRS = CRS.decode("EPSG:4326");
//			CoordinateReferenceSystem targetCRS = CRS.decode("EPSG:23032");
//
//			MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
//			System.out.println(transform.toString());
			
			
			/*
			 * Transforming a Geometry
			 * A MathTransform, as generated above, can be used by bashing away at the interface 
			 * and feeding it DirectPosition objects one at a time.
			 * Or you could break out the JTS utility class where this work has been done for you:
			 */
//			Geometry targetGeometry = JTS.transform( sourceGeometry, transform);//an die sourceGeometry muss man ersma rankommen
//			// Transforming an ISO Geometry is more straight forward:
//			Geometry target = geometry.transform( targetCRS );
			
			
			/*
			 * Axis Order
			 * One thing that often comes up is the question of axis order.
			 * The EPSG database often defines axis in an order that is inconvenient for display; 
			 * we have a method to quickly check what is going on
			 */
//			if( CRS.getAxisOrder( coordinateReferenceSystem ) == CRS.AxisOrder.LAT_LON){
//				// lat lon
//			}
			
			
			
			/*
			 * Several often-used values are available as predefined constants.
			 * Here is an example of accessing several of the predefined constants:
			 */
//			GeographicCRS geoCRS = org.geotools.referencing.crs.DefaultGeographicCRS.WGS84;
//		    GeodeticDatum wgs84Datum = org.geotools.referencing.datum.DefaultGeodeticDatum.WGS84;
//		    PrimeMeridian greenwichMeridian = org.geotools.referencing.datum.DefaultPrimeMeridian.GREENWICH;
//		    CartesianCS cartCS = org.geotools.referencing.cs.DefaultCartesianCS.GENERIC_2D;
//		    CoordinateSystemAxis latAxis = org.geotools.referencing.cs.DefaultCoordinateSystemAxis.GEODETIC_LATITUDE;
			
			
			/*
			 * Creating a CoordinateReferenceSystem
			 * You can use factories defined by the referencing system to create things by 
			 * hand using java code.
			 * This example shows the creation of a WGS84 / UTM 10N CoordinateReferenceSystem:
			 */
//			MathTransformFactory mtFactory = ReferencingFactoryFinder.getMathTransformFactory(null);
//		    CRSFactory crsFactory = ReferencingFactoryFinder.getCRSFactory(null);
//		    
//		    GeographicCRS geoCRS = org.geotools.referencing.crs.DefaultGeographicCRS.WGS84;
//		    CartesianCS cartCS = org.geotools.referencing.cs.DefaultCartesianCS.GENERIC_2D;
//		    
//		    ParameterValueGroup parameters = mtFactory.getDefaultParameters("Transverse_Mercator");
//		    parameters.parameter("central_meridian").setValue(-111.0);
//		    parameters.parameter("latitude_of_origin").setValue(0.0);
//		    parameters.parameter("scale_factor").setValue(0.9996);
//		    parameters.parameter("false_easting").setValue(500000.0);
//		    parameters.parameter("false_northing").setValue(0.0);
//		    Conversion conversion = new DefiningConversion("Transverse_Mercator", parameters);
//		    
//		    Map<String, ?> properties = Collections.singletonMap("name", "WGS 84 / UTM Zone 12N");
//		    ProjectedCRS projCRS = crsFactory.createProjectedCRS(properties, geoCRS, conversion, cartCS);
			
			
			
			
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		
		

	}
}
