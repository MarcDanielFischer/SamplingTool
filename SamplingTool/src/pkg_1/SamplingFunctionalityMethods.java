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
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.swing.data.JFileDataStoreChooser;
import org.geotools.referencing.CRS;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.AttributeDescriptor;
import org.opengis.geometry.BoundingBox;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
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
	 * Derive UTM Zone for a given longitude value.
	 * @param longitude
	 * @return
	 */
	public static int long2UTM(double longitude){
		// TODO Ausnahmen: Norwegen etc. --> also auch abhängig von Latitude --> evtl behandeln --> start: welche Zonen sind überhaupt irregulär?
		int utmZone = (int)Math.floor(((longitude + 180) / 6) +1);
		if(utmZone > 60) utmZone = utmZone % 60; // if input longitude is > 180 for some reason (and output UTM zone is > 60 then)
		return utmZone;
	}
	
	
	/**
	 * Convenience method to convert processed data in UTM projection to geographic LonLat projection  
	 */
	public static void UTM2LonLat(ArrayList<Point> plots) throws Exception{
		
		// iterate over plots and convert them
		for(int i = 0; i < plots.size(); i++){
			Point currentPlot = plots.get(i);
			int srid = currentPlot.getSRID(); // TODO check for srid to be set (must not be 0)
			CoordinateReferenceSystem sourceCRS = CRS.decode("EPSG:" + srid);
			Point plotLatLon = (Point)JTS.toGeographic(currentPlot, sourceCRS);
			plots.set(i, plotLatLon); // replace elements in ArrayList instead of writing them to a new one in order to save memory
			
		}

	}
	
	
	/**
	 * This method is responsible for the entire Sampling process.
	 * It takes the sampling params from the GUI, implements the sampling logic
	 * and delegates to specific methods (eg, Cluster Building) 
	 */
	public static boolean runSampling(File inputFile, String sampleColumn, int numStrata,  String[] selectedStrata, int samplingDesign,  int[] numPlotsToBeSampled, int gridDistX, int gridDistY, int startingPoint, int startX, int startY, int clusterSampling, int clusterShape, int numSubPlotsinHVerticalLine, int numClusterSubPlots, int distBetweenSubPlots ){
		
		
		// get all features in the shapefile
		try {
			ArrayList<SimpleFeature> allFeatures = getFeatures(inputFile); // dieser SChritt könnte performancegefährdend sein, weil ALLE features eingelsen werden --> evtl ändern
		
		// filter selected features
			ArrayList<SimpleFeature> selectedFeatures = new ArrayList<SimpleFeature>();
			for(int i = 0; i < selectedStrata.length; i++){ // evtl Zeit sparen, indem ich hier optimiert vorgehe
				for(SimpleFeature feature : allFeatures){ // TODO hier mit while-loop, damit der Schleifendurchlauf abgekürzt wird
					if(feature.getAttribute(sampleColumn).toString().equals(selectedStrata[i])){
						selectedFeatures.add(feature);
					}
				}
			}
			
			
		// convert to UTM
			// read SHP CRS (sourceCRS)
			CoordinateReferenceSystem sourceCRS = getCRS(inputFile); // note: there is only ONE source CRS, as a shapefile has only one associated .prj file
			
			
			// iterate over selected features, convert them to UTM(more specifically, only their geometries, as other attribs are irrelevant here)
			// and add converted Geometries to ArrayList
			// --> this CRS transformation has to be done in a loop because each feature possibly needs to be transformed to its own CRS
			ArrayList<Geometry> strataUTM= new ArrayList<Geometry>();
			
			for(SimpleFeature feature : selectedFeatures){
				CoordinateReferenceSystem targetCRS = getTargetCRS(feature); // targetCRS:in UTM projection
				// transform source CRS directly (without looking up the corresponding EPSG code first and use that CRS instead)
				MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
				// get feature Geometry
				Geometry sourceGeometry = (Geometry) feature.getAttribute( "the_geom" );
				// convert Geometry to target CRS
				Geometry targetGeometry = JTS.transform( sourceGeometry, transform);
				
				// add CRS info to targetGeometry so that the points can have it, too, and then they can be 
				// reprojected from UTM to LatLon which they could not i
				int srid = CRS.lookupEpsgCode(targetCRS, true);
				targetGeometry.setSRID(srid);
				// add stratum name to targetGeometry so that we can retrieve that name easily when it comes to writing output data
				targetGeometry.setUserData(feature.getAttribute(sampleColumn));
				strataUTM.add(targetGeometry);
				
			}
			
			
			// create ArrayList to store and append output plots (ArrayList is a more convenient type to append data to than a regular array)
			ArrayList<Point> outputPlots = new ArrayList<Point>();
			
			// iterate over converted Geometries and sample plots according to input params
			for(int i = 0; i < strataUTM.size(); i++){
				// sample plots


				// muss man diese if-Verzweigungen IN den for loop reinschreiben oder kann man das irgendwie trennen?
				// --> ich könnte erst die Samplingart bestimmen und dann
				// den jeweiligen Fach-Sample-Methoden einfach die ganze ArrayList übergeben und die loopen dann 
				// --> aber wie ist das, wenn die Straten unterschiedliche CRS haben? dann müsste eine MEthode über mehrere Straten mit versch.
				// CRS iterieren und in jedes auf unterschiedlich Weise Punkte hineinsampeln


				// hier die Samplinglogik
				// simple random sample
				if(samplingDesign == GUI_Designer.SIMPLE_RANDOM_SAMPLING){
					if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_NO){
						// 1)simple random sample, no cluster plots
						ArrayList<Point> simpleRandomPlots = simpleRandomSampling(strataUTM.get(i), numPlotsToBeSampled[i]);
						outputPlots.addAll(simpleRandomPlots);

					}
					// 2)simple random sample with cluster plots
					if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
						// first create cluster center points using simple random Sampling // cluster_coord <- coordinates(SRS(strata_utm, n = number_of_plots[l]))

						// Second, make a cluster out of each cluster center point
						if(clusterShape == GUI_Designer.I_SHAPE){
							// output_plots <- I_plots(cluster_coord[], plot_dist[l], cluster_nrplots[l])
						}
						if(clusterShape == GUI_Designer.L_SHAPE){

						}
						if(clusterShape == GUI_Designer.SQUARE_SHAPE){

						}
						if(clusterShape == GUI_Designer.H_SHAPE){

						}
					}
				}

				// systematic sampling
				if(samplingDesign == GUI_Designer.SYSTEMATIC_SAMPLING){

					// random or specified starting point?
					if(startingPoint == GUI_Designer.STARTING_POINT_RANDOM){

					}
					if(startingPoint == GUI_Designer.STARTING_POINT_SPECIFIED){

					}

					// 3) systematic sample , no cluster plots
					if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_NO){

					}
					// 4) systematic sample with cluster plots
					if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
						if(clusterShape == GUI_Designer.I_SHAPE){
							// output_plots <- I_plots(cluster_coord, plot_dist[l], cluster_nrplots[l])
						}
						if(clusterShape == GUI_Designer.L_SHAPE){

						}
						if(clusterShape == GUI_Designer.SQUARE_SHAPE){

						}
						if(clusterShape == GUI_Designer.H_SHAPE){

						}
					}
				}

			}
			
			
			// reproject output plots to LatLon
			UTM2LonLat(outputPlots);
			
			// write output to file
			File outputFile = getFile();
			writeOutput(outputFile, outputPlots);
			
			// if everything went well...
			return true;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		// if things went not so well (this statement is never reached if everything works fine)
		return false;
	}
	
	
	
	
	
	
	
	/**
	 * Gets a matching UTM CoordinateReferenceSystem for the input SimpleFeature.
	 * If the feature spans more than one UTM zone, the formula (zoneMin + zoneMax) / 2 is used.
	 * Behaviour for cross-equator-features is not tested yet. 
	 * @param feature
	 * @return
	 */
	private static CoordinateReferenceSystem getTargetCRS(SimpleFeature feature) throws Exception{
		// find matching UTM zone for feature (UTM zones depend mainly on Longitude values)
		BoundingBox boundingBox = feature.getBounds();
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
		return targetCRS;
	}
	
	

	/**
	 * Convenience method to get the CoordinateReferenceSystem of an input SHP.
	 * @param inputFile
	 * @return the CRS of the input file
	 * @throws Exception
	 */
	private static CoordinateReferenceSystem getCRS(File inputFile) throws Exception{
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(inputFile);
		CoordinateReferenceSystem crs = dataStore.getSchema().getCoordinateReferenceSystem();
		return crs;
	}


	/**
	 * new Version, use this one
	 * @param stratum
	 * @param numPlots
	 * @return
	 */
	public static ArrayList<Point> simpleRandomSampling(Geometry stratum, int numPlots){
		// initialize output ArrayList
		ArrayList<Point> output = new ArrayList<Point>();
		
		// BBOX, min/max values --> aber das ist doch doof, dann benutze ich an 2 Stellen BBOX (umständlich), einmal bei CRS-Konversion und hier --> kann man das vermeiden?
		Envelope envelope = stratum.getEnvelopeInternal();
		double minX = envelope.getMinX();
		double maxX = envelope.getMaxX();
		double minY = envelope.getMinY();
		double maxY = envelope.getMaxY();
		
		// sample Plots within bbox
		while(output.size() < numPlots){
			// generate random X
			double X = (Math.random() * (maxX - minX)) + minX;
			// generate random Y
			double Y = (Math.random() * (maxY - minY)) + minY;
			// generate a point from X, Y values created above
			GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
			Coordinate coord = new Coordinate( X, Y );
			Point point = geometryFactory.createPoint( coord );
			// add two values to the point that are not automatically copied from stratum
			point.setSRID(stratum.getSRID()); // point needs srid property so it can be reprojected to LatLon later on 
			point.setUserData(stratum.getUserData()); // we want to add the stratum name to the point so that we can write it to output file later on

			// check if Point inside stratum (not just within bbox) and add to output ArrayList
			if(point.within(stratum)){
				output.add(point);
			}
		}
		
		// return ArrayList
		return output;
	}
	
	
	/**
	 * DEPRECATED, use overload method instead
	 * 
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
	
	
	public static void writeOutput(File file, ArrayList<Point> samplePlots) throws Exception{
		// TODO vl doch BufferedWriter nehmen --> ist der irgendwie sicherer?

		FileWriter fileWriter = new FileWriter(file, true); // second param is for appending to file (Yes/No)
		
		// the following two lines are needed to derive appropriate index numbers for the plots
		String plotStratumName = null;
		int j = 0;
		
		// iterate over ArrayList and write Point values to output file
		for(int i = 0; i < samplePlots.size(); i++){
			Point currentPlot = samplePlots.get(i);
			double x = currentPlot.getCoordinate().x;
			double y = currentPlot.getCoordinate().y;
			fileWriter.write(Double.toString(x) + ",");
			fileWriter.write(Double.toString(y) + ",");
			// the following if-else block is to derive the plot index number for the output file (index starts at one for each stratum)
			if (plotStratumName == (String)currentPlot.getUserData()){ // if plotStratumName is same as for the plot before, then continue increasing index j
				j++;				
			}else{
				plotStratumName = (String)currentPlot.getUserData(); // if plotStratumName has changed, then reset index j to 1
				j = 1;
			}
			fileWriter.write((j) + ","); // "Plot_nr" derived from samplePlots-ArrayList index position 
			fileWriter.write((String)currentPlot.getUserData() + "\n"); // der Polygonname fehlt jetzt noch
		}
		fileWriter.close();
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
