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

public class ShapeFileUtilities {

	
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
		ArrayList<SimpleFeature> features = getSHPFeatures(inputShapeFile);
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
	public static ArrayList<SimpleFeature> getSHPFeatures(File inputShapeFile) throws Exception{ 
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
