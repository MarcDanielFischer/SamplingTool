package production;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.AttributeDescriptor;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.vividsolutions.jts.geom.Geometry;

/**
 * This class provides static methods to access data contained in Shapefiles.
 * @author daniel
 *
 */
public class ShapeFileUtilities {


	/**
	 * Calls getStratum() for the length of the strataNames array.
	 * @param shapefile
	 * @param column
	 * @param strataNames
	 * @return
	 * @throws Exception
	 */
	public static Stratum[] getStrata(File shapefile, String column, String[] strataNames) throws Exception{
		Stratum[] strataArray = new Stratum[strataNames.length];
		for(int i = 0; i < strataNames.length; i++){
			strataArray[i] = getStratum(shapefile, column, strataNames[i]);
		}
		return strataArray;
	}

	
	
	/**
	 * Builds a Stratum from Shapefile Features that have the same name,
	 * @param shapeFile
	 * @param column
	 * @param stratumName
	 * @return
	 * @throws Exception
	 */
	public static Stratum getStratum(File shapeFile, String column, String stratumName) throws Exception{ 
		// 1 Stratum may consist of one or more Features
		ArrayList<SimpleFeature> features = getFeatures(shapeFile, column, stratumName);
		
		/*
		 *  Get Stratum geometry.
		 *  If there is more than one feature, merge geometries.
		 */
		Geometry geom = null;
		if(features.size() > 1){
			 geom = mergeFeatureGeometries(features); 
		}else{
			geom = (Geometry) features.get(0).getAttribute("the_geom");
		}
		
		// build Stratum object 
		CoordinateReferenceSystem crs = ShapeFileUtilities.getSHPCRS(shapeFile);
		Stratum returnStratum = new Stratum(geom, crs, stratumName);
		return returnStratum;
	}
	
	

	/**
	 * Gets the Features that constitute a Stratum.
	 * A Stratum may consist of a single Feature or
	 * multiple Features with the same name.
	 * As this is unknown beforehand, this methods returns an ArrayList object containing one
	 * or more elements, according to the input data.
	 * @param shapeFile
	 * @param column
	 * @param stratumName
	 * @return an ArrayList containing the Stratum Features.
	 * @throws Exception
	 */
	private static ArrayList<SimpleFeature> getFeatures(File shapeFile, String column, String stratumName) throws Exception{
		// Access Shapefile
		// GeoTools logic: Input File --> DataStore --> FeatureSource --> FeatureCollection --> Iterator
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(shapeFile);
		SimpleFeatureSource featureSource = dataStore.getFeatureSource(); // wird benötigt, um an einzelne Features ranzukommen (Feature = Zeile in SHP Attribute Table)
		final SimpleFeatureCollection featureCollection = featureSource.getFeatures(); // wieso muss die final sein?		
		SimpleFeatureIterator iterator = featureCollection.features(); // an die Features(=Zeilen in SHP attribute table) kommt man, indem man über die FeatureCollection iteriert 

		// iterate over features
		// it may happen that a feature name is not unique, so we write all features with the specified name to an ArrayList
		ArrayList<SimpleFeature> features = new ArrayList<SimpleFeature>(); //
		try {
			while( iterator.hasNext()){
				SimpleFeature f = iterator.next();
				//features.add(feature);
				if(f.getAttribute(column).toString().equals(stratumName)){
					features.add(f);
				}

			}
		}
		finally {
			iterator.close(); // prevents memory leaks and data loss
		}

		return features;
	}



	/**
	 * Multi-Feature-Strata problem: 
	 * 
	 * Strata may consist of more than a single feature. In order to treat them as a single entity, 
	 * they must somehow be merged into a single object. 
	 * 
	 * The way we do this here is to combine multi-feature-strata geometries into 
	 * one single geometry (if the geometries are adjacent, they can be combined into a Polygon object,
	 * otherwise they will be a MultiPolygon) before reprojecting them to UTM (otherwise they would not be combinable any more 
	 * as they would possibly be reprojected to different UTM zones)
	 * @param features The features to be merged. Must all use the same CRS
	 * @return the merged Geometries of all input features  
	 */
	private static Geometry mergeFeatureGeometries(ArrayList<SimpleFeature> features){
		Geometry geomMerge = (Geometry) features.get(0).getAttribute("the_geom");
		for(int i = 1; i < features.size(); i++){
			Geometry g = (Geometry) features.get(i).getAttribute("the_geom");
			geomMerge = geomMerge.union(g);
		}

		return geomMerge;
	}

	
	
	
	
	
	/**
	 * This method takes a Shapefile as a File object as input parameter and returns
	 * the columns of the Shapefile as a String[] array. 
	 * @param inputShapeFile File object. Must be a Shapefile
	 * @return the Shapefile column names as String[] array
	 * @throws Exception
	 */
	public static String[] getSHPColNames(File inputShapeFile) throws Exception {
		// Access SHP Column Names (GeoTools logic): Input File --> DataStore --> FeatureType (Anzahl und Art der Spalten in Attribute Table) --> AttributeDescriptors --> getLocalName()
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(inputShapeFile);
		SimpleFeatureType simpleFeatureType = dataStore.getSchema();
		List<AttributeDescriptor> attributeDescriptors = simpleFeatureType.getAttributeDescriptors(); // gibt noch andere Möglichkeiten, an die Spaltennamen zu kommen, zb direkt über die features selbst
		String[] columnNames = new String[simpleFeatureType.getAttributeCount()];
		for(int i = 0; i < attributeDescriptors.size(); i++){
			columnNames[i] = attributeDescriptors.get(i).getLocalName();
		}
		return columnNames;
		
	}
	
	/**
	 * Returns the values (in the intended use case of this method: feature names, eg strata names) contained 
	 * in the specified column of the input Shapefile.
	 * @param inputShapeFile
	 * @param strataColumn
	 * @return Values contained in a specified column of the input Shapefile.
	 * @throws Exception
	 */
	public static ArrayList<String> getSHPColumnValues (File inputShapeFile, String strataColumn) throws Exception{
		// initialize output ArrayList
		ArrayList<String> values = new ArrayList<String>();
		// get all features contained in the Shapefile
		ArrayList<SimpleFeature> features = getAllSHPFeatures(inputShapeFile);
		// iterate over features and extract the specified column´s attribute values 
		for(SimpleFeature feature : features){
			values.add((String)feature.getAttribute(strataColumn));
		}
		// sort values alphabetically 
		Collections.sort(values);
		// remove duplicate values
		for(int i = 0; i < values.size() -1 ; i++){
			// wenn es einen nächsten gibt, anschauen
			while(true){
				// if there is a next element in the ArrayList...
				if(i+1 <= values.size() -1){ // i+1 is the index for the value following the current value. We must make sure that it is not ouf of bounds (ie, pointing to an index not contained in the ArrayList)
					if(values.get(i).equals(values.get(i+1))){
						values.remove(i+1);
					}else{
						break;
					}
				}else{
					break;

				}
			}
		}
		return values;
	}
	
	/**
	 * Convenience Method to get all features of a Shapefile nicely accessible in an ArrayList object.
	 * @param inputShapeFile
	 * @return The features contained in the input Shapefile.
	 * @throws Exception
	 */
	public static ArrayList<SimpleFeature> getAllSHPFeatures(File inputShapeFile) throws Exception{ 
		// Access Features in an SHP Column (GeoTools logic): Input File --> DataStore --> FeatureSource --> FeatureCollection --> Iterator
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(inputShapeFile);
		SimpleFeatureSource featureSource = dataStore.getFeatureSource(); // wird benötigt, um an einzelne Features ranzukommen (Feature = Zeile in SHP Attribute Table)
		final SimpleFeatureCollection featureCollection = featureSource.getFeatures(); // wieso muss die final sein?		
		SimpleFeatureIterator iterator = featureCollection.features(); // an die Features(=Zeilen in SHP attribute table) kommt man, indem man über die FeatureCollection iteriert 
		ArrayList<SimpleFeature> features = new ArrayList<SimpleFeature>(); //
		try {
			while( iterator.hasNext()  ){
				SimpleFeature feature = iterator.next();
				features.add(feature);
				
			}
		}
		finally {
			iterator.close(); // prevents memory leaks and data loss
		}
		return features;
	}
	
	/**
	 * Convenience method to get the CoordinateReferenceSystem of an input SHP.
	 * @param inputFile
	 * @return the CRS of the input file
	 * @throws Exception
	 */
	public static CoordinateReferenceSystem getSHPCRS(File inputFile) throws Exception{
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(inputFile);
		CoordinateReferenceSystem crs = dataStore.getSchema().getCoordinateReferenceSystem();
		return crs;
	}
	
	
}
