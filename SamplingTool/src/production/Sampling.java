
package production;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.geotools.coverage.grid.GridCoverage2D;
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
import com.vividsolutions.jts.geom.Polygon;


/**
 * This class contains methods needed to execute the sampling process.
 * Its core part is the runSampling() method that performs all the steps of the actual sampling operation
 * and calls the other methods of this class (and those contained in RasterProcessing).
 * @author daniel
 *
 */
public class Sampling {
	
	// TODO die ganzen input params evtl als ParameterValueGroup (Key-Value_Pair) übergeben
	/**
	 * This method is the entry point for the entire Sampling process.
	 * It takes the sampling parameters from the GUI, implements the sampling logic
	 * and delegates to specific methods (eg, Cluster Building). 
	 * 
	 * @param inputShapeFile
	 * @param sampleColumn
	 * @param selectedStrata
	 * @param samplingDesign
	 * @param numPlotsToBeSampled
	 * @param gridDistX
	 * @param gridDistY
	 * @param startingPoint
	 * @param startX
	 * @param startY
	 * @param clusterSampling
	 * @param clusterShape
	 * @param numSubPlotsinHVerticalLine
	 * @param numSubPlotsinHhorizontalLine
	 * @param numClusterSubPlots
	 * @param distBetweenSubPlots
	 * @param bufferSize
	 * @param weightedSampling
	 * @param inputRasterFile
	 * @throws Exception
	 */
	public static void runSampling(File inputShapeFile, String sampleColumn, String[] selectedStrata, int samplingDesign,  int[] numPlotsToBeSampled, int gridDistX, int gridDistY, int startingPoint, double startX, double startY, int clusterSampling, int clusterShape, int numSubPlotsinHVerticalLine, int numSubPlotsinHhorizontalLine, int numClusterSubPlots, int distBetweenSubPlots, double bufferSize, boolean weightedSampling, File inputRasterFile ) throws Exception {

		// TODO getFeatures() probably not necessary --> just iterate is enough 
		// get all features in the shapefile
		ArrayList<SimpleFeature> allFeatures = ShapeFileUtilities.getSHPFeatures(inputShapeFile); // dieser Schritt könnte performancegefährdend sein, weil ALLE features eingelsen werden --> evtl ändern

		// filter selected features
		ArrayList<SimpleFeature> selectedFeatures = new ArrayList<SimpleFeature>();

		for(int i = 0; i < selectedStrata.length; i++){ // evtl Zeit sparen, indem ich hier optimiert vorgehe
			for(SimpleFeature currentFeature : allFeatures){ 
				if(currentFeature.getAttribute(sampleColumn).toString().equals(selectedStrata[i])){
					selectedFeatures.add(currentFeature);
				}
			}
		}


		/*
		 * Multi-Feature-Strata problem: 
		 * 
		 * Strata may consist of more than a single feature. In order to treat them as a single entity, 
		 * they must somehow be merged into a single object. 
		 * 
		 * The way we do this here is to combine multi-feature-strata geometries into 
		 * one single geometry (if the geometries are adjacent, they can be combined into a Polygon object,
		 * otherwise they will be a MultiPolygon) before reprojecting them to UTM (otherwise they would not be combinable any more 
		 * as they would possibly be reprojected to different UTM zones)
		 */

		// iterate over selectedFeatures and merge Geometries using union() if they have the same name (eg, belong to the same stratum)
		if(selectedFeatures.size() > selectedStrata.length){ // perform this step only if there actually is at least one Multi-Feature stratum
			for(int i = 0; i < selectedFeatures.size()-1; i++){ // iterate until penultimate element so that there is always a following element (avoid ArrayIndexOutOfBoundsException)
				// while the next feature belongs to the same stratum as the current feature --> merge Geometries
				String nameCurrentFeature = selectedFeatures.get(i).getAttribute(sampleColumn).toString();
				while(true){
					if(i+1 <=  selectedFeatures.size()-1){ // avoid IndexOutOfBoundsException while looking at the next element 
						String nameNextFeature = selectedFeatures.get(i+1).getAttribute(sampleColumn).toString();

						if(nameCurrentFeature.equals(nameNextFeature)){
							// extract Geometries from both Features
							Geometry geom1 = (Geometry) selectedFeatures.get(i).getAttribute("the_geom");
							Geometry geom2 = (Geometry) selectedFeatures.get(i+1).getAttribute("the_geom");

							// merge Geometries using union()
							// geomMerge may be of type Polygon or MultiPolygon, depending on whether the merged Geometries are adjacent to each other or not
							Geometry geomMerge = geom1.union(geom2);

							// set as Geometry for feature i
							// works also for Polygon objects although column type is specified as MultiPolygon 
							selectedFeatures.get(i).setAttribute("the_geom", geomMerge);

							// remove feature i+1 from ArrayList
							selectedFeatures.remove(i+1);
						}else{
							break;
						}
					}else{
						break;
					}
				}
			}
		}




		// reproject selected features to UTM
		// read input SHP CRS (sourceCRS)
		CoordinateReferenceSystem sourceCRS = ShapeFileUtilities.getSHPCRS(inputShapeFile); // note: there is only ONE source CRS, as a shapefile has only one associated .prj file


		// iterate over selected features, convert them to UTM(more specifically, only their geometries, as other attribs are irrelevant here),
		// apply a negative buffer to the Geometries 
		// and add converted Geometries to ArrayList
		// --> this CRS transformation has to be done in a loop because each feature possibly needs to be transformed to its own CRS (different UTM zones) according to its position
		ArrayList<Stratum> strataUTM= new ArrayList<Stratum>();

		for(SimpleFeature currentFeature : selectedFeatures){
			CoordinateReferenceSystem targetUTM_CRS = CRSUtilities.getTargetUTM_CRS(currentFeature, sourceCRS); // targetCRS:in UTM projection
			// transform source CRS directly (without looking up the corresponding EPSG code first and use that CRS instead)
			MathTransform transform = CRS.findMathTransform(sourceCRS, targetUTM_CRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
			// get feature Geometry
			Geometry sourceGeometry = (Geometry) currentFeature.getAttribute( "the_geom" );
			// convert Geometry to target CRS
			Geometry utmGeometry = JTS.transform( sourceGeometry, transform);

			// apply negative buffer to targetGeometry before creating Stratum object with it and starting the actual sampling process 
			Geometry bufferedUTMGeometry = utmGeometry.buffer(-bufferSize);

			// make a Stratum object out of the targetGeometry that holds additional information (CRS, stratumName )
			Stratum stratum = new Stratum(bufferedUTMGeometry, targetUTM_CRS, (String)currentFeature.getAttribute(sampleColumn));
			strataUTM.add(stratum);
		}



		// create ArrayList to store and append output plots (ArrayList is a more convenient type to append data to than a regular array)
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();

		// iterate over converted Geometries and sample plots according to input params
		for(int i = 0; i < strataUTM.size(); i++){

			// muss man diese if-Verzweigungen IN den for loop reinschreiben oder kann man das irgendwie trennen?
			// --> ich könnte erst die Samplingart bestimmen und dann
			// den jeweiligen Fach-Sample-Methoden einfach die ganze ArrayList übergeben und die loopen dann 
			// --> aber wie ist das, wenn die Straten unterschiedliche CRS haben? dann müsste eine MEthode über mehrere Straten mit versch.
			// CRS iterieren und in jedes auf unterschiedlich Weise Punkte hineinsampeln


			//-------------------------------------------------------------------------------------------------------------
			// Sampling logic
			
			ArrayList<Plot> clusterSeedPoints = null; // must initialize variable outside if block so that it is visible in other if blocks
			
			// create simple Plots first (random or systematic)
			
			// Simple random sample, no weighted Sampling 
			if(samplingDesign == GUI_Designer.SIMPLE_RANDOM_SAMPLING && weightedSampling == false){

				// Plots are generated regardless of Clustering choice
				ArrayList<Plot> simpleRandomPlots = simpleRandomSampling(strataUTM.get(i), numPlotsToBeSampled[i]);
				
				// in case of no clustering write generated Plots to output directly
				if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_NO){
					outputPlots.addAll(simpleRandomPlots);
				}

				// in case of clustering use generated Plots as cluster seed points
				if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
					clusterSeedPoints = simpleRandomPlots;
				}
			}
			
			
			///////////////////////////////////////////////////////////////////////////
			// weighted Sampling 
			///////////////////////////////////////////////////////////////////////////
			// Simple random sample, weighted Sampling 
			if(samplingDesign == GUI_Designer.SIMPLE_RANDOM_SAMPLING && weightedSampling == true){
				// TODO: weighted Sampling: how can all this go into a method on its own?
				GridCoverage2D coverage = RasterProcessing.readGeoTiff(inputRasterFile);
				
				
				// die Geometries sind hier schon alle zu UTM konvertiert (verschieden Zonen, je nach Lage)
				// --> ich könnte sie wieder zurückkonvertieren, aber das wäre irgendwie blöd
				// gibts die auch noch im Original? und was ist dann mit den MultiPolyStrata?
				// RasterProcessing.getClipGeometries() holt die Geometries nochmal frisch aus dem
				// SHP raus und reprojected sie zum rasterCRS und kommt auch mit MultiPolyStrata gut klar
				ArrayList<Geometry> clipGeoms = RasterProcessing.getClipGeometries(coverage.getCoordinateReferenceSystem(), inputShapeFile, sampleColumn, strataUTM.get(i).getName());
				GridCoverage2D clippedCoverage = RasterProcessing.getClippedCoverage(clipGeoms, coverage);
				
				double maxValue = RasterProcessing.getCoverageMaxValue(clippedCoverage, 0);
				double noDataValue = RasterProcessing.getNoDataValue(coverage, 0);
				
				// rejection testing: call  simpleRandomSampling() one at a time until desired number of plots is reached
				/*
				 * look at simpleRandomSampling() implementation more closely 
				 * --> reduce processing costs by not calculating a stratum BBOX for every single plot etc.
				 */
				
				ArrayList<Plot> weightedSimpleRandomPlots = new ArrayList<Plot>();
				int numSampledPlots = 0;
				do{
					// create a Plot (one at a time, has to be rejection tested)
					/**
					 * TODO simpleRandomSampling() berechnet hier für jeden einzelnen 
					 * plot neu die BBOX vom gleichen Stratum --> könnte ich 
					 * auf einmal verkürzen, spart vl Rechenzeit
					 */
					Plot plot = simpleRandomSampling(strataUTM.get(i), 1).get(0);
					// get plot value 
					// TODO warum geht double[] hier auch für ganzzahlige Rasters? --> debug step through
					double plotValue = RasterProcessing.getValueAtPosition(coverage, plot, (double[]) null)[0];
					// throw exception if plot has nodata value
					if(plotValue == noDataValue){
						throw new Exception("Plot does not have a corresponding raster pixel value.\n"
								+ " Make sure the weight raster completely covers the stratum area.");
					}
					// get plot weight 
					double plotWeight = plotValue / maxValue;
					// rejection test the plot
					boolean keepPlot = RasterProcessing.rejectionTesting(plotWeight);
					if(keepPlot == false){ //if plot is rejected, sample a new one
						continue;
					}else{
						numSampledPlots++;
						plot.setPlotNr(numSampledPlots);
						plot.setWeight(plotWeight);
						weightedSimpleRandomPlots.add(plot);
						
					}
					
				}
				while(numSampledPlots < numPlotsToBeSampled[i]);
				
				// Plots are generated regardless of Clustering choice
//				ArrayList<Plot> weightedSimpleRandomPlots = simpleRandomSampling(strataUTM.get(i), numPlotsToBeSampled[i]);

				// in case of no clustering write generated Plots to output directly
				if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_NO){
					outputPlots.addAll(weightedSimpleRandomPlots);
				}

				// in case of clustering use generated Plots as cluster seed points
				if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
					clusterSeedPoints = weightedSimpleRandomPlots;
				}
			}

			// Systematic sampling
			if(samplingDesign == GUI_Designer.SYSTEMATIC_SAMPLING){

				Point startPointUTM = null; // // must initialize variable outside if block so that it is visible in other if blocks

				// get a starting point both for random and specified starting point options
				if(startingPoint == GUI_Designer.STARTING_POINT_RANDOM){
					// get a random starting point using simpleRandomSampling()
					startPointUTM = simpleRandomSampling(strataUTM.get(i), 1).get(0).getPoint(); // dieser Punkt hat bereits UTM-Koords, weil input stratumUTM schon in UTM ist

				}
				if(startingPoint == GUI_Designer.STARTING_POINT_SPECIFIED){
					// make a Point out of given start point coords (LonLat)
					GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
					Coordinate coord = new Coordinate( startX, startY );
					Point startPointLonLat = geometryFactory.createPoint( coord ); // dieser Punkt ist in LonLat, weil startX und startY LonLat sind und hat kein associated CRS

					// reproject startPointLonLat using a MathTransform that transforms between sourceCRS and targetCRS
					// source CRS: we use the already defined input file's CRS here as the grid start point refers to a location within the input file area
					// target CRS: use the current stratum's CRS (as grid points are to be located inside this stratum, even if the starting point is located oustide the stratum)
					CoordinateReferenceSystem targetCRS = strataUTM.get(i).getCRS();
					// transform source CRS directly (without looking up the corresponding EPSG code first and use that CRS instead)
					MathTransform transform = CRS.findMathTransform(sourceCRS, targetCRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
					// convert Point to target CRS
					startPointUTM = (Point) JTS.transform( startPointLonLat, transform);
				}

				// call systematicSampling() using the newly created startPoint
				// Plots are generated regardless of Clustering choice
				ArrayList<Plot> systematicPlots = systematicSampling(strataUTM.get(i),startPointUTM, gridDistX, gridDistY );

				// in case of no clustering write generated Plots to output directly
				if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_NO){
					outputPlots.addAll(systematicPlots);
				}
				
				// in case of clustering use generated Plots as cluster seed points
				if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){
					clusterSeedPoints = systematicPlots;
				}
			}
			
			// ----------------------------------------------------------
			// make Clusters out of simple Plots if needed
			
			if(clusterSampling == GUI_Designer.CLUSTER_SAMPLING_YES){

				// Second, make a cluster out of each cluster center point
				if(clusterShape == GUI_Designer.I_SHAPE){
					// output_plots <- I_plots(cluster_coord[], plot_dist[l], cluster_nrplots[l])
					ArrayList<Plot> i_Plots = Clusters.create_I_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, strataUTM.get(i));
					outputPlots.addAll(i_Plots);
				}
				if(clusterShape == GUI_Designer.L_SHAPE){
					ArrayList<Plot> l_Plots = Clusters.create_L_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, 1, strataUTM.get(i));
					outputPlots.addAll(l_Plots);
				}
				if(clusterShape == GUI_Designer.L_SHAPE_UPSIDE_DOWN){
					ArrayList<Plot> l_Plots = Clusters.create_L_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, -1, strataUTM.get(i));
					outputPlots.addAll(l_Plots);
				}
				if(clusterShape == GUI_Designer.SQUARE_SHAPE){
					ArrayList<Plot> square_Plots = Clusters.create_Square_clusters(clusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, strataUTM.get(i));
					outputPlots.addAll(square_Plots);
				}
				if(clusterShape == GUI_Designer.SQUARE_SHAPE_ROTATED){
					ArrayList<Plot> square_Plots = Clusters.create_rotated_Square_clusters(clusterSeedPoints, distBetweenSubPlots, strataUTM.get(i));
					outputPlots.addAll(square_Plots);
				}
				if(clusterShape == GUI_Designer.H_SHAPE){
					ArrayList<Plot> h_Plots = Clusters.create_H_clusters(clusterSeedPoints, distBetweenSubPlots, numSubPlotsinHVerticalLine,  numSubPlotsinHhorizontalLine, strataUTM.get(i));
					outputPlots.addAll(h_Plots);
				}
			}
			// End Cluster options
			// ----------------------------------------------------------
			
			// End Sampling logic
			//-------------------------------------------------------------------------------------------------------------

		}

		if(outputPlots.size()== 0){
			throw new Exception("No output plots have been generated");
		}else{
			// reproject output plots to LatLon
			CRSUtilities.UTM2LonLat(outputPlots);

		// write output to file
		File outputFile = OutputUtilities.getFile();
		OutputUtilities.writeCSVoutput(outputFile, outputPlots, clusterSampling, weightedSampling);
		}
		
	}

	

	/**
	 * Generates a specified number of output Plots in a given Stratum using simple random sampling. 
	 * @param stratum
	 * @param numPlots
	 * @return output plots
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
	 * Generates plots in a given Stratum using systematic sampling.
	 * 
	 * @param stratum The stratum to generate the plots in.
	 * @param startPoint grid starting point
	 * @param gridDistX distance between grid points in x direction
	 * @param gridDistY distance between grid points in y direction
	 * @return output plots
	 */
	public static ArrayList<Plot> systematicSampling(Stratum stratum, Point startPoint, int gridDistX, int  gridDistY ){
		
		// initialize output ArrayList
		ArrayList<Plot> output = new ArrayList<Plot>();
				
		// BBOX
		Envelope boundingBbox = stratum.getGeometry().getEnvelopeInternal();
		
		// BBBOX min / max values
		double minX = boundingBbox.getMinX();
		double maxX = boundingBbox.getMaxX();
		double minY = boundingBbox.getMinY();
		double maxY = boundingBbox.getMaxY();
		
		// determine numPoints in all directions (seen from the starting point, how many points are there in each direction until reaching the edges of the boundingBox)
		double xStart = startPoint.getCoordinate().x;
		double yStart = startPoint.getCoordinate().y;
		int numPointsLeft = (int)Math.floor((xStart - minX) / gridDistX); 
		int numPointsRight = (int)Math.floor((maxX - xStart) / gridDistX); 
		int numPointsUp = (int)Math.floor((maxY -yStart) / gridDistY); 
		int numPointsDown = (int)Math.floor((yStart - minY) / gridDistY); 
		
		// calculate totals: numPointsHorizontal (number of points in one grid row covering the whole width of the bounding box) and numPointsVertical
		// (number of points in one grid column covering the whole length of the bounding box)
		int numPointsHorizontal = numPointsLeft + numPointsRight + 1; // numPoints to the right and left of the starting point + starting point itself
		int numPointsVertical = numPointsUp + numPointsDown + 1; // numPoints up and down + starting point itself
		
		// create first grid point in the upper left corner
		// coords
		double x = xStart - (gridDistX * numPointsLeft);
		double y = yStart + (gridDistY * numPointsUp);
		// generate a point from X, Y values created above using GeometryFactory
		GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
		Coordinate coord;
		Point point;
		
		// create the grid (must use Point objects instead of just x,y coords to be able to check whether they are within the stratum later on)
		// initialize ArrayList object to store all grid Points
		ArrayList<Point> gridPoints = new ArrayList<Point>(); 
		// nested for loop: create grid Points starting from upper left corner of the grid, from top to bottom (same order as when reading a page of text) 
		for(int i = 0; i < numPointsVertical; i++ ){
			 
			 for(int k = 0; k < numPointsHorizontal; k++){
				 // create a Point 
				 coord = new Coordinate( x, y );
				 point = geometryFactory.createPoint( coord );
				 gridPoints.add(point);
				 // increase x coord
				 x += gridDistX;
				  
			 }
			 // reset x to its lefternmost value ("carriage return")
			 x = xStart - (gridDistX * numPointsLeft);
			 
			 // shift y coord one row down
			 y -= gridDistY;
		 }
		
		// check whether Points are within the stratum and only create Plots from Points within stratum
		int plotNr = 1; // index for the Plots to be created
		for(Point p : gridPoints){
			if(p.within(stratum.getGeometry())){
				Plot plot = new Plot(p, stratum.getName(), stratum.getCRS(), plotNr);
				output.add(plot);
				plotNr ++;
			}
		}
		return output;
	}
	

}
