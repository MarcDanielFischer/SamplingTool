/**
 * Hier kommen die Methoden rein, die die eigentliche Funktionalität des 
 * Sampling-Tools ausmachen (Trennung GUI / Funktionalität)
 * 
 * - TODO: Methoden von Katjas Draft hier reinpacken
 */

package pkg_1;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.swing.data.JFileDataStoreChooser;
import org.geotools.referencing.CRS;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.AttributeDescriptor;
import org.opengis.geometry.BoundingBox;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;

public class SamplingFunctionalityMethods {
	

	/**
	 * This method takes a Shapefile as a File object as input parameter and returns
	 * the columns of the Shapefile as a String[] array. 
	 * @param inputFile File object. Must be a Shapefile
	 * @return the Shapefile column names as String[] array
	 * @throws Exception
	 */
	public static String[] getSHPColNames(File inputFile) throws Exception {
		// Access SHP Column Names (GeoTools logic): Input File --> DataStore --> FeatureType (Anzahl und Art der Spalten in Attribute Table) --> AttributeDescriptors --> getLocalName()
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(inputFile);
		SimpleFeatureType simpleFeatureType = dataStore.getSchema();
		List<AttributeDescriptor> attributeDescriptors = simpleFeatureType.getAttributeDescriptors(); // gibt noch andere Möglichkeiten, an die Spaltennamen zu kommen, zb direkt über die features selbst
		String[] columnNames = new String[simpleFeatureType.getAttributeCount()];
		for(int i = 0; i < attributeDescriptors.size(); i++){
			columnNames[i] = attributeDescriptors.get(i).getLocalName();
		}
		return columnNames;
		
	}
	
	/**
	 * Returns the values (in the intended use case of this method: feature names, eg strata names) contained in the specified column
	 * @param inputFile
	 * @param strataColumn
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<String> getStrataColumnValues (File inputFile, String strataColumn) throws Exception{
		ArrayList<String> values = new ArrayList<String>();
		ArrayList<SimpleFeature> features = getFeatures(inputFile);
		for(SimpleFeature feature : features){
			values.add((String)feature.getAttribute(strataColumn));
		}
		return values;
	}
	
	
	
	/**
	 * Convenience Method to get all features of a Shapefile nicely accessible in an ArrayList object
	 * --> check for possible performance losses, as FeatureCollection does not load everything into memory
	 * @param inputFile
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<SimpleFeature> getFeatures(File inputFile) throws Exception{ 
		// TODO klären: wie ist das, wenn mehrere Polygone gleich heißen? (EIN Stratum besteht aus mehreren Polygonen, alle sollen besampelt werden)
		// Access Features in an SHP Column (GeoTools logic): Input File --> DataStore --> FeatureSource --> FeatureCollection --> Iterator
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(inputFile);
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
			iterator.close(); // prevents memory leaks or data loss or whatever
		}
		return features;
	}
	
	
	
	/**
	 * Convenience method to convert input data in geographic LonLat projection to UTM projection
	 * Method: - long2UTM #function to derive UTM zone from longitude coordinate (hat sie selbst geschrieben)
     * verwendete Formel: (floor((longitude + 180)/6) + 1) %% 60 --> UTM-System verstehen --> nachlesen
     * - funktioniert wahrscheinlich nicht für alle UTM-Zonen (es gibt diese irregulär geformten Zonen im Norden
     * --> evtl für Skandinavien relevant)
	 */
	public static void lonLat2UTM(){
		/*
		 * evtl hilfreiche Code Snippets
		 * (source: http://gis.stackexchange.com/questions/110249/coordinate-conversion-epsg3857-to-epsg4326-using-opengis-jts-too-slow)
		 * - org.geotools.referencing.GeodeticCalculator geht anscheinend auch
		 */
//		CoordinateReferenceSystem targetCRS = CRS.decode("EPSG:4326");
//        CoordinateReferenceSystem sourceCRS = CRS.decode("EPSG:3857");
//        MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS, false);
//        GeometryFactory geometryFactory = new GeometryFactory(new PrecisionModel(), 3857);
//        Point minPoint = geometryFactory.createPoint(new Coordinate(xMin, yMin));
//        Point minTargetPoint = (Point) JTS.transform(minPoint, transform);
//        Point maxPoint = geometryFactory.createPoint(new Coordinate(xMax, yMax));
//        PointmaxTargetPoint = (Point) JTS.transform(maxPoint, transform);
	}
	
	/**
	 * Convenience method to convert processed data in UTM projection to geographic LonLat projection  
	 */
	public static void UTM2LonLat(){

	}
	
	
	/**
	 * entry point: von hier wird weiterverzweigt zu 
	 * simpleRandomSampling() bzw. systematicSampling()
	 * --> runSampling() hat ALLE Parameter
	 */
	public static boolean runSampling(File inputFile, String sampleColumn, int numStrata,  String[] selectedStrata, int samplingDesign,  int[] numPlotsToBeSampled, int gridDistX, int gridDistY, int startingPoint, int startX, int startY, int clusterSampling, int clusterShape, int numSubPlotsinHVerticalLine, int numClusterSubPlots, int distBetweenSubPlots ){
		// TODO Methode umstrukturieren: alle output-Zeilen sammeln und dann auf einen Schlag schreiben 
		// TODO Sampling Methods so umschreiben, dass sie selber intern über mehrere Straten iterieren (for-loop dort und nicht hier)
		/*
		 *  gewünschte Struktur hier:
		 *  1) calls to individual sampling methods (simple random, cluster..)
		 *  2) save output lines in List or so
		 *  3) write entire output at once 
		 */
		
		
		// get all features in the shapefile
		try {
			ArrayList<SimpleFeature> allFeatures = getFeatures(inputFile); // dieser SChritt könnte performancegefährdend sein, weil ALLE features eingelsen werden --> evtl ändern
			// filter selected features
			ArrayList<SimpleFeature> selectedFeatures = new ArrayList<SimpleFeature>();
			for(int i = 0; i < selectedStrata.length; i++){ // evtl Zeit sparen, indem ich hier optimiert vorgehe
				for(SimpleFeature feature : allFeatures){
					if(feature.getAttribute(sampleColumn).toString().equals(selectedStrata[i])){
						selectedFeatures.add(feature);
					}
				}
			}
			/*
			 *  convert features to UTM (only necessary if SHP not in metric system --> check)
			 *  - this is only needed for cluster and systematic sampling, not for SRS
			 */
			// if(SHP not already in metric projection so that we can calculate distances directly):
			for(SimpleFeature feature : selectedFeatures){
				// Katja: 
				// determine which UTM Zone a stratum lies in 
				   //  BBOX --> lonMin, lonMax // --> wie geht das mit GeoTools?
				   //  derive UTM zones for lonMin, lonMax (only lon values relevant for UTM zones (except exceptions))

				// reproject Stratum to specified UTM Zone
				   // spTransform(feature, targetCRS)
				
				
				//-----------------------------
				// --> wie geht das mit GeoTools?
				// UTM Zone 32N (Deutschland): EPSG code 32632
				// GeometryCRS lab
				// class Mathtransform...
				
				
//              // Bsp aus GeoTools Documentation:  
//				CoordinateReferenceSystem sourceCRS = CRS.decode("EPSG:4326");
//				CoordinateReferenceSystem targetCRS = CRS.decode("EPSG:23032");
//				MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS, true);
				
//				// So kommt man an das aktuelle CRS der verwendeten Daten heran:
//				SimpleFeatureType schema = featureSource.getSchema();
//				CoordinateReferenceSystem dataCRS = schema.getCoordinateReferenceSystem();
				
				
				
			}
			
			
			
			
			
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		//----------------------------------------------------------------------------------------
		// ich übernehme jetzt Katjas Struktur, alles andere ist mir gerade zu schwierig
//		for(strata_utm in strata_utm_list){
//			
//		}
		
		
		//----------------------------------------------------------------------------------------
		
		
		
		
		
		
		
		
		// SRS hard-wired --> funktioniert
		try{
			File outputFile = getFile(); // creates output file and writes header line into it
			// loop over polygons and call sampling method for each of them
			for(int i = 0; i < numStrata; i++){
				Point[] outputPlots = simpleRandomSampling(inputFile, sampleColumn, selectedStrata[i], numPlotsToBeSampled[i]);
				// write plots to file in each iteration:
				writeOutput(outputFile, outputPlots, selectedStrata[i]);
			}
			return true;
		}catch(Exception e){
			JOptionPane.showMessageDialog(null, e.toString());
		}
		
		
		
		
		
//      // eigener Versuch
//		// suboptimale Struktur
//		try{
//			// first create output file and writes header line into it --> move to one single location insted of "distributed writing"
//			File outputFile = getFile(); 
//			Point[] outputPlots;
//			// delegate to appropriate sampling Method
//			
//			// Simple Random Sampling
//			if(samplingDesign == GUI_Designer.SIMPLE_RANDOM_SAMPLING){
//				// call simpleRandomSampling() multiple times (one for each selected stratum) 
//				// TODO for loop von hier wegbewegen: simpleRandomSampling() soll selbstständig über mehrere Straten iterieren und jeweils das Stratum mit dazu schreiben
//				for(int i = 0; i < numStrata; i++){
//					outputPlots = simpleRandomSampling(inputFile, sampleColumn, selectedPolygons[i], numPlotsToBeSampled[i]);
//					// write plots to file in each iteration:
//					writeOutput(outputFile, outputPlots, selectedPolygons[i]);
//				}
//			}
//			
//			if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
//				clusterSampling(File inputFile, String sampleColumn, int numStrata,  String[] selectedPolygons, int samplingDesign,  int[] numPlotsToBeSampled, int gridDistX, int gridDistY, int startingPoint, int startX, int startY, int clusterSampling, int clusterShape, int numSubPlotsinHVerticalLine, int numClusterSubPlots, int distBetweenSubPlots);
//			}
//
//
//			return true;
//		}catch(Exception e){
//			JOptionPane.showMessageDialog(null, e.toString());
//		}
		
		return false;
	}
	
	
	public static void writeOutput(File file, Point[] samplePlots, String selectedPolygon) throws Exception{
		// Problem hier: für jede iteration wird das file ganz überschrieben
		// --> ich will append, nicht overwrite
		// TODO vl doch BufferedWriter nehmen --> ist der irgendwie sicherer?

		FileWriter fileWriter = new FileWriter(file, true);
		for(int j = 0; j <= samplePlots.length -1; j++){
			double x = samplePlots[j].getCoordinate().x;
			double y = samplePlots[j].getCoordinate().y;
			fileWriter.write(Double.toString(x) + ",");
			fileWriter.write(Double.toString(y) + ",");
			fileWriter.write((j+1) + ","); // "Plot_nr" derived from samplePlots-array index position 
			fileWriter.write(selectedPolygon + "\n");
		}
		fileWriter.close();
}
	

	/**
	 * Simple Random Sampling umsetzen
	 * Was brauche ich hier?
	 * - Zugriff auf einzelne Polygone
	 * - BoundingBox  --> getBounds()
	 * - Random Coordinates innerhalb dieser BoundingBox
	 * - Punkte erzeugen --> mit GeometryFactory, siehe Bsp unten
	 * - schauen, wie es die Katja gemacht hat
	 * 
	 * So hat Katja es gemacht:
	 * - Simple Random Sample Method (Params: shape, number of points to be sampled)
	 * - dafür wird die fertige Funktion "spsample" aus package "sp" verwendet
	 *   generiert Punkte in UTM-Projektion, numeriert sie und benennt X- und Y-Spalte (in Java wahrsch. nicht notwendig)
	 * - shape = SpatialPolygonsDataFrame created through import of shapefile -> inventory area (UTM projection) 
	 */
	public static Point[] simpleRandomSampling(File file, String column,  String selectedPolygon, int numPlots) throws Exception{ // TODO param list erweitern
		// TODO klären: wie ist das, wenn mehrere Polygone gleich heißen? (EIN Stratum besteht aus mehreren Polygonen, alle sollen besampelt werden)
		
		SimpleFeature samplePolygon = null; // Initialize to "null" and then search it in while loop
		ArrayList<SimpleFeature> features = getFeatures(file);
		for(SimpleFeature feature : features){
			if(feature.getAttribute(column).toString().equals(selectedPolygon)) samplePolygon= feature;
		}
		
		// Array der spezifizierten Länge für Outputpunkte anlegen
		// (darf nicht im if-Block stehen, da return value)
		Point[] samplePlots = new Point[numPlots];
		
		//-------------------------------------------------------------------
		// Start Sampling Loop
		/*
		 * Sampling Loop
		 * 
		 *  so lange Punkte sampeln, bis
		 *  gewünschte Anzahl erreicht ist, und dann auf output csv schreiben
		 */
		if(samplePolygon != null){ // wenn das gesuchte Polygon (Stratum) gefunden ist...
			// Punkte ins gewählte Polygon hineinsampeln: ersma BoundingBox holen
			BoundingBox boundingBox = samplePolygon.getBounds();

			
			int i = 0; // Index für Array-Befüllung (unschöne Lösung hier, geht vl besser)
			//Flag hier notwendig, weil gesampelte Punkte auch außerhalb des Polygons liegen können, daher keine bestimmte Anzahl Iterationen spezifizierbar
			boolean samplePointsArrayFullFlag = false;
			// Array mit Punkten befüllen, bis er voll ist 
			while(!samplePointsArrayFullFlag){
				// jetzt einen random Point innerhalb der BoundingBox-Grenzen erzeugen
				double minX = boundingBox.getMinX();
				double maxX = boundingBox.getMaxX();
				double minY = boundingBox.getMinY();
				double maxY = boundingBox.getMaxY();
				// X(lon) and Y(lat): random double values between min and max values of BoundingBox
				double lon = 1000; // nur damit unten nicht kommt "variable may not have been initialized" --> unschön, evtl, verbessern
				boolean lonWithinBBox = false;
				/*
				 * generate random number in a given interval:
				 * double random = (Math.random() * (upper - lower)) + lower; // --> das ist die Lösung , so gehts!
				 */
				while(!lonWithinBBox){
					lon = (Math.random() * (maxX - minX)) + minX;
					// TODO prüfen: Überprüfung (ob lonWithinBBox) überhaupt noch notwendig?
					if (lon >= minX && lon <= maxX) lonWithinBBox = true;
				}
				double lat = 1000; // nur damit unten nicht kommt "variable may not have been initialized" --> unschön, evtl, verbessern
				boolean latWithinBBox = false;
				while(!latWithinBBox){
					lat = (Math.random() * (maxY - minY)) + minY;
					// TODO prüfen: Überprüfung (ob latWithinBBox) überhaupt noch notwendig?
					if (lat >= minY && lat <= maxY) latWithinBBox = true;
				}
				// generate a point from lon,lat values created above
				GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
				Coordinate coord = new Coordinate( lon, lat );
				Point point = geometryFactory.createPoint( coord );
				/*
				 *  check, ob Punkte innerhalb des POlygons liegen:
				 *  Operation "Within": One geometry is completely within another (no touching edges):
				 *  wenn die Punkte genau AUF der Außenbegrenzung liegen dürfen, muss ich hier
				 *  eine andere Geometrieoperation nehmen
				 */
				MultiPolygon geomSamplepolygon = (MultiPolygon)samplePolygon.getAttribute("the_geom");
				// Add point to output points array, but only if it is inside polygon
				if(point.within(geomSamplepolygon)){
					samplePlots[i] = point; 
					i++;
					if(i >= samplePlots.length) samplePointsArrayFullFlag = true; 
				}
			}
			// End Sampling Loop
			//-------------------------------------------------------------------
			
			
			
//			//-------------------------------------------------------------------
//			// Start Save Output stuff
//			// --> try and move this block out of the sampling method for multiple polygon sampling
//			int confirmDialogResult = JOptionPane.showConfirmDialog(null, "samplePlots[] - Array successfully generated. Do you want to save your output?");
//			if(confirmDialogResult == JOptionPane.YES_OPTION){
//				writeCsvOutput(samplePlots, selectedPolygon);
//			}else{
//				return;
//			}
//			// End Save Output stuff
//			//-------------------------------------------------------------------
			
		}
		return samplePlots;
	}

	public static void systematicSampling(){ // (Params: shape, x, y, dx, dy)
		// systematic sample with given starting point 
		// if no starting point given: call simpleRandomSampling() to create one

	}
	
	public static void clusterSampling(){
		// Weiterverzweigung zu den verschiedenen Clusterformen (I,L,H,Square)
	}
	
	public static void create_I_clusters(){

	}
	
	public static void create_L_clusters(){

	}
	
	public static void create_H_clusters(){

	}
	
	public static void create_Square_clusters(){

	}
	
	/**
	 * Get output file. This method shows a saveFileDialog etc.
	 * This method also writes the file header to prevent writing it multiple times, as the subsequent calls to writeOutput()
	 * are being executed inside a loop. 
	 * @return
	 */
	private static File getFile() throws Exception{
		// Problem hier: Methode muss File-Objekt zurückgeben. Wenn User den SaveFileDialog abbricht, ohne eine Datei anzugeben,
		// muss ich hier nach der Programmierlogik eine Exception werfen
		// TODO automatically add ".csv" file extension
		JFileChooser saveFileDialog = new JFileChooser(); // JFileChooser kommt aus javax.swing und hat nichts mit GeoTools zu tun 
		//    FileNameExtensionFilter filter = new FileNameExtensionFilter(
		//        "CSV files", "csv");
		//    chooser.setFileFilter(filter);
		saveFileDialog.setDialogType(JFileChooser.SAVE_DIALOG);
		if(saveFileDialog.showOpenDialog(null) != JFileChooser.APPROVE_OPTION) { // return from method if users cancels the dialog (ie, user does not click on "Save"-Button)
		throw new Exception("SaveFileDialog cancelled without specifying an output file");
		}
		// create output file 
		File saveFile = saveFileDialog.getSelectedFile();
		if(saveFile.exists()){
			int fileOverwriteChoice = JOptionPane.showConfirmDialog(null, "File already exists. Do you want to overwrite the existing file?");
			if(fileOverwriteChoice == JOptionPane.YES_OPTION){
				saveFile.createNewFile();
			}else{
				throw new Exception("you have chosen not to overwrite an existing file");
			}

		}else{
			saveFile.createNewFile();
		}
		
		//-----------------------------
		//write file header here, because  writeOutput() is called multiple times in a loop and the header is only needed once
		FileWriter fileWriter = new FileWriter(saveFile, true);
		fileWriter.write("\"X\",\"Y\",\"Plot_nr\",\"Stratum\"\n");
		fileWriter.close();
		//-----------------------------
		
		
		return saveFile;
	}
	
	
	
	/**
	 * TODO OBSOLET (?)
	 * Writes the generated sample plots to output CSV file 
	 */
	public static void writeCsvOutput(Point[] samplePlots, String selectedPolygon) throws Exception{ // Exception muss immer dazu, sonst beschwert sich Swing oder GeoTools
		// 1) get File to write to
		// TODO automatically add ".csv" file extension
		// saveFileDialog
		JFileChooser saveFileDialog = new JFileChooser(); // JFileChooser kommt aus javax.swing und hat nichts mit GeoTools zu tun 
		//    FileNameExtensionFilter filter = new FileNameExtensionFilter(
		//        "CSV files", "csv");
		//    chooser.setFileFilter(filter);
		saveFileDialog.setDialogType(JFileChooser.SAVE_DIALOG);
		if(saveFileDialog.showOpenDialog(null) != JFileChooser.APPROVE_OPTION) { // return from method if users cancels the dialog (ie, user does not click on "Save"-Button)
			return; // beendet nur die aktuelle Methode(File wird nicht gespeichert, Programm läuft weiter) --> System.exit(0) beendet das ganze Programm
		}
		// create output file 
		File saveFile = saveFileDialog.getSelectedFile();
		if(saveFile.exists()){
			int fileOverwriteChoice = JOptionPane.showConfirmDialog(null, "File already exists. Do you want to overwrite the existing file?");
			if(fileOverwriteChoice == JOptionPane.YES_OPTION){
				saveFile.createNewFile();
			}else{
				JOptionPane.showMessageDialog(null, "you have chosen not to overwrite an existing file");
				return; // beendet nur die aktuelle Methode(File wird nicht gespeichert, Programm läuft weiter) --> System.exit(0) beendet das ganze Programm
			}

		}else{
			saveFile.createNewFile();
		}
		//-------------------------------------------------------------------
		// 2) write data to output file using FileWriter
		// TODO vl doch BufferedWriter nehmen --> ist der irgendwie sicherer?

		FileWriter fileWriter = new FileWriter(saveFile);
		fileWriter.write("\"X\",\"Y\",\"Plot_nr\",\"Stratum\"\n");
		for(int j = 0; j <= samplePlots.length -1; j++){
			double x = samplePlots[j].getCoordinate().x;
			double y = samplePlots[j].getCoordinate().y;
			fileWriter.write(Double.toString(x) + ",");
			fileWriter.write(Double.toString(y) + ",");
			fileWriter.write((j+1) + ","); // "Plot_nr" derived from samplePlots-array index position 
			fileWriter.write(selectedPolygon + "\n");
		}
		fileWriter.close();

		//-------------------------------------------------------------------
		// 3) show success status message window
		JOptionPane.showMessageDialog(null, "File data successfully written");
	}
	
}
