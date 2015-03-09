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
import org.opengis.referencing.crs.GeographicCRS;
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
	 * This method is responsible for the entire Sampling process.
	 * It takes the sampling params from the GUI, implements the sampling logic
	 * and delegates to specific methods (eg, Cluster Building) 
	 */
	public static boolean runSampling(File inputFile, String sampleColumn, int numStrata,  String[] selectedStrata, int samplingDesign,  int[] numPlotsToBeSampled, int gridDistX, int gridDistY, int startingPoint, int startX, int startY, int clusterSampling, int clusterShape, int numSubPlotsinHVerticalLine, int numClusterSubPlots, int distBetweenSubPlots ){
	// TODO check: param "numStrata" obsolet?
		
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
			
			
		// reproject selected features to UTM
			// read input SHP CRS (sourceCRS)
			CoordinateReferenceSystem sourceCRS = getCRS(inputFile); // note: there is only ONE source CRS, as a shapefile has only one associated .prj file
			
			
			// iterate over selected features, convert them to UTM(more specifically, only their geometries, as other attribs are irrelevant here)
			// and add converted Geometries to ArrayList
			// --> this CRS transformation has to be done in a loop because each feature possibly needs to be transformed to its own CRS (different UTM zones) according to its position
			ArrayList<Stratum> strataUTM= new ArrayList<Stratum>();
			
			for(SimpleFeature feature : selectedFeatures){
				CoordinateReferenceSystem targetCRS = getTargetCRS(feature); // targetCRS:in UTM projection
				// transform source CRS directly (without looking up the corresponding EPSG code first and use that CRS instead)
				MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
				// get feature Geometry
				Geometry sourceGeometry = (Geometry) feature.getAttribute( "the_geom" );
				// convert Geometry to target CRS
				Geometry targetGeometry = JTS.transform( sourceGeometry, transform);
				
				
				// make a Stratum object out of the targetGeometry that holds additional information (CRS, stratumName )
				Stratum stratum = new Stratum(targetGeometry, targetCRS, (String)feature.getAttribute(sampleColumn));
				strataUTM.add(stratum);
				

				
			}
			
			
			// create ArrayList to store and append output plots (ArrayList is a more convenient type to append data to than a regular array)
			ArrayList<Plot> outputPlots = new ArrayList<Plot>();
			
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
						ArrayList<Plot> simpleRandomPlots = simpleRandomSampling(strataUTM.get(i), numPlotsToBeSampled[i]);
						outputPlots.addAll(simpleRandomPlots);

					}
					
					// 2)simple random sample with cluster plots
					if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
						
						// first create cluster center points using simple random Sampling // cluster_coord <- coordinates(SRS(strata_utm, n = number_of_plots[l]))
						ArrayList<Plot> clusterSeedPoints = simpleRandomSampling(strataUTM.get(i), numPlotsToBeSampled[i]);
						
						// Second, make a cluster out of each cluster center point
						if(clusterShape == GUI_Designer.I_SHAPE){
							// output_plots <- I_plots(cluster_coord[], plot_dist[l], cluster_nrplots[l])
							ArrayList<Plot> i_Plots = create_I_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, strataUTM.get(i));
							outputPlots.addAll(i_Plots);
						}
						if(clusterShape == GUI_Designer.L_SHAPE){
							ArrayList<Plot> l_Plots = create_L_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, strataUTM.get(i));
							outputPlots.addAll(l_Plots);
						}
						if(clusterShape == GUI_Designer.SQUARE_SHAPE){
							ArrayList<Plot> square_Plots = create_Square_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, strataUTM.get(i));
							outputPlots.addAll(square_Plots);
						}
						if(clusterShape == GUI_Designer.H_SHAPE){
							ArrayList<Plot> h_Plots = create_H_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, numSubPlotsinHVerticalLine,  strataUTM.get(i));
							outputPlots.addAll(h_Plots);
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
			writeOutput(outputFile, outputPlots, clusterSampling);
			
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
	public static ArrayList<Plot> simpleRandomSampling(Stratum stratum, int numPlots){
		// initialize output ArrayList
		ArrayList<Plot> output = new ArrayList<Plot>();
		
		// BBOX, min/max values 
		// TODO  BBOX wird an 2 Stellen verwendet (umständlich), einmal bei CRS-Konversion und hier --> kann man das vermeiden?
		Geometry stratumGeometry = stratum.getGeometry();
		Envelope envelope = stratumGeometry.getEnvelopeInternal(); // Envelope is Bounding Box
		double minX = envelope.getMinX();
		double maxX = envelope.getMaxX();
		double minY = envelope.getMinY();
		double maxY = envelope.getMaxY();
		
		int plotNr = 1; // initialize outside while loop and increase only if plot added
		
		// randomly create Point objects within bbox
		while(output.size() < numPlots){
			// generate random X
			double X = (Math.random() * (maxX - minX)) + minX;
			// generate random Y
			double Y = (Math.random() * (maxY - minY)) + minY;
			
			// generate a point from X, Y values created above using GeometryFactory
			GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
			Coordinate coord = new Coordinate( X, Y );
			Point point = geometryFactory.createPoint( coord );
			

			// check if Point inside stratum (not just within bbox), create Plot and add Plot to output ArrayList
			if(point.within(stratumGeometry)){
				// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
				Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), plotNr);
				output.add(plot);
				plotNr ++; // increase plotNr only after successfully adding new Plot to output (eg, when randomly generated point falls inside BBOX)
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
	
	
	/**
	 * Grow i-shaped clusters from given seed points with a total of numClusterSubPlots Plots in each cluster and distBetweenSubPlots distance between them (distance in Meters).
	 * As for now, the clusters are linearly grown from south to north (the cluster seed point is the southernmost plot). 
	 * TODO check if cluster growing mode must be changed (seed point in the center of the line)
	 * @param clusterSeedPoints
	 * @param distBetweenSubPlots
	 * @param numClusterSubPlots
	 * @param stratum this param is needed in order to check whether all generated cluster plots are inside the stratum 
	 */
	public static ArrayList<Plot> create_I_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, Stratum stratum){
		
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();
		
		// iterate over seed points
		for(Plot seedPoint : clusterSeedPoints){
			
			
			// use seedPoint.plotNr as clusterNr and set new seedPoint.plotNr to 1
			int clusterNr = seedPoint.getPlotNr();
			int subPlotNr = 1;
			seedPoint.setClusterNr(clusterNr);
			seedPoint.setPlotNr(subPlotNr);
			
			// add each seed point to output after changing its numbering
			outputPlots.add(seedPoint);
			
			double x = seedPoint.getPoint().getCoordinate().x;
			double y = seedPoint.getPoint().getCoordinate().y;
			
			// create additional points until reaching the specified total number of plots per cluster
			for(int i = 0; i < numClusterSubPlots - 1; i++){
				
				y += distBetweenSubPlots; // build i-clusters along North-South axis (change only y coordinate)				
				
				// create Points using GeometryFactory
				GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
				Coordinate coord = new Coordinate( x, y );
				Point point = geometryFactory.createPoint( coord );

				// check if Point inside stratum (not just within bbox), create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					outputPlots.add(plot);
					
				}
			}
		}
		return outputPlots;
	}
	
	
	
	/**
	 * Grow L-shaped clusters from given seed points with a total of numClusterSubPlots Plots in each cluster and distBetweenSubPlots distance between them (distance in Meters).
	 * In case of an even total number of Plots shaping the cluster, the vertical axis will be one Plot longer than the horizontal axis (just like the real character "L").
	 * @param clusterSeedPoints
	 * @param distBetweenSubPlots
	 * @param numClusterSubPlots
	 * @param stratum this param is needed in order to check whether all generated cluster plots are inside the stratum 
	 */
	public static ArrayList<Plot> create_L_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, Stratum stratum){
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();

		// iterate over seed points
		for(Plot seedPoint : clusterSeedPoints){

			// use seedPoint.plotNr as clusterNr and set new seedPoint.plotNr to 1
			int clusterNr = seedPoint.getPlotNr();
			int subPlotNr = 1;
			seedPoint.setClusterNr(clusterNr);
			seedPoint.setPlotNr(subPlotNr);

			// add each seed point to output after changing its numbering
			outputPlots.add(seedPoint);

			if(numClusterSubPlots > 1){ // if numClusterSubPlots == 1, the seed point alone will be its own cluster
				
				// Determine number of plots to be created along each axis of the L
				int numSubPlotsVerticalAxis = 0;
				int numSubPlotsHorizontalAxis = 0;
				
				if(numClusterSubPlots % 2 != 0){ // odd total number of Plots shaping the cluster
					// in case of an odd total Plot number, both axes are of the same length
					numSubPlotsVerticalAxis = (int)Math.floor(numClusterSubPlots / 2);
					numSubPlotsHorizontalAxis = (int)Math.floor(numClusterSubPlots / 2);
				} else{ // even total number of Plots shaping the cluster
					// in case of an even total Plot number, the vertical axis will be one Plot longer than the horizontal axis
					numSubPlotsVerticalAxis = numClusterSubPlots / 2;
					numSubPlotsHorizontalAxis = (numClusterSubPlots / 2) -1;
				}


				// extract seed point coords and use them to build the other Plots
				double x = seedPoint.getPoint().getCoordinate().x;
				double y = seedPoint.getPoint().getCoordinate().y;


				// build Plots along vertical axis
				for(int i = 1; i <= numSubPlotsVerticalAxis; i++){
					y += distBetweenSubPlots; // as we build along the vertical axis, only y coords are affected

					// create Points using GeometryFactory
					GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
					Coordinate coord = new Coordinate( x, y );
					Point point = geometryFactory.createPoint( coord );

					// check if Point inside stratum (not just within bbox), create Plot and add Plot to output ArrayList
					if(point.within(stratum.getGeometry())){
						subPlotNr = i * 2; // vertical axis Plots have  even numbers
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						outputPlots.add(plot);

					}
				}


				// reset y coord (which has been increased during vertical axis Plot construction)  so that we can start building plots along the horizontal axis
				y = seedPoint.getPoint().getCoordinate().y;

				// build plots along horizontal axis
				for(int i = 1; i <= numSubPlotsHorizontalAxis; i++){
					x += distBetweenSubPlots; // as we build along the vertical axis, only x coords are affected

					// create Points using GeometryFactory
					GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
					Coordinate coord = new Coordinate( x, y );
					Point point = geometryFactory.createPoint( coord );

					// check if Point inside stratum (not just within bbox), create Plot and add Plot to output ArrayList
					if(point.within(stratum.getGeometry())){
						subPlotNr = (i * 2)  + 1  ; // horizontal axis Plots have odd numbers
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						outputPlots.add(plot);

					}
				}

			}
		}
		return outputPlots;

	}
	
	
	
	
	public static ArrayList<Plot> create_H_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, int numSubPlotsinHVerticalLine, Stratum stratum){

//		ArrayList<Plot> outputPlots = new ArrayList<Plot>();
//
//		// iterate over seed points
//		for(Plot seedPoint : clusterSeedPoints){
//
//
//			// use seedPoint.plotNr as clusterNr and set new seedPoint.plotNr to 1
//			int clusterNr = seedPoint.getPlotNr();
//			int subPlotNr = 1;
//			seedPoint.setClusterNr(clusterNr);
//			seedPoint.setPlotNr(subPlotNr);
//			
//			// add each seed point to output after changing its numbering
//			outputPlots.add(seedPoint);
//			
//			double x = seedPoint.getPoint().getCoordinate().x;
//			double y = seedPoint.getPoint().getCoordinate().y;
//			
//			// create additional points until reaching the specified total number of plots per cluster
//			for(int i = 0; i < numClusterSubPlots - 1; i++){
//				
//				y += distBetweenSubPlots; // build i-clusters along North-South axis (change only y coordinate)				
//				
//				// create Points using GeometryFactory
//				GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
//				Coordinate coord = new Coordinate( x, y );
//				Point point = geometryFactory.createPoint( coord );
//
//				// check if Point inside stratum (not just within bbox), create Plot and add Plot to output ArrayList
//				if(point.within(stratum.getGeometry())){
//					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
//					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
//					Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
//					outputPlots.add(plot);
//					
//				}
//			}
//		}
//		return outputPlots;
	}
	
	
	
	
	
	public static ArrayList<Plot> create_Square_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, Stratum stratum){
		return null;
	}
	
	/**
	 * Get output file. This method shows a saveFileDialog etc.
	 * 
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
		
			
		
		return saveFile;
	}
	
	
	/**
	 * 
	 * @param file
	 * @param samplePlots
	 * @param clusterSampling needed in order to adapt file header (determine whether or not to add column "cluster_nr")
	 * @throws Exception
	 */
	public static void writeOutput(File file, ArrayList<Plot> samplePlots, int clusterSampling) throws Exception{
		// TODO vl doch BufferedWriter nehmen --> ist der irgendwie sicherer?

		FileWriter fileWriter = new FileWriter(file, true); // second param is for appending to file (Yes/No)
		
		// write file header
		// adapt header for CLUSTER_SAMPLING_YES option
		if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){ 
			fileWriter.write("\"X\",\"Y\",\"Cluster_nr\",\"Plot_nr\",\"Stratum\"\n");
			
		} else{
			fileWriter.write("\"X\",\"Y\",\"Plot_nr\",\"Stratum\"\n"); // header without Clusters
		}
		
		
		for(Plot plot : samplePlots){
			
			
			double x = plot.getPoint().getCoordinate().x;
			double y = plot.getPoint().getCoordinate().y;
			fileWriter.write(Double.toString(x) + ",");
			fileWriter.write(Double.toString(y) + ",");
			// write clusterNr only to output if Cluster Sampling is the chosen sampling option
			if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
				fileWriter.write(Integer.toString(plot.getClusterNr())+ ",");
			}
			fileWriter.write(Integer.toString(plot.getPlotNr()) + ",");
			fileWriter.write(plot.getStratumName() + "\n");
			
			
		}
		
		fileWriter.close();
		
		
//		// the following two lines are needed to derive appropriate index numbers for the plots
//		String plotStratumName = null;
//		int j = 0;
//		
//		// iterate over ArrayList and write Point values to output file
//		for(int i = 0; i < samplePlots.size(); i++){
//			Point currentPlot = samplePlots.get(i);
//			double x = currentPlot.getCoordinate().x;
//			double y = currentPlot.getCoordinate().y;
//			fileWriter.write(Double.toString(x) + ",");
//			fileWriter.write(Double.toString(y) + ",");
//			// the following if-else block is to derive the plot index number for the output file (index starts at one for each stratum)
//			if (plotStratumName == (String)currentPlot.getUserData()){ // if plotStratumName is same as for the plot before, then continue increasing index j
//				j++;				
//			}else{
//				plotStratumName = (String)currentPlot.getUserData(); // if plotStratumName has changed, then reset index j to 1
//				j = 1;
//			}
//			fileWriter.write((j) + ","); // "Plot_nr" derived from samplePlots-ArrayList index position 
//			fileWriter.write((String)currentPlot.getUserData() + "\n"); // der Polygonname fehlt jetzt noch
//		}
//		fileWriter.close();
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
