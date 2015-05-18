
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
	 * @param strataNames
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
	public static ArrayList<Plot> runSampling(File inputShapeFile, String sampleColumn, String[] strataNames, int samplingDesign,  int[] numPlotsToBeSampled, int gridDistX, int gridDistY, boolean startPointSpecified, double startX, double startY, boolean clusterSampling, int clusterShape, int numSubPlotsinHVerticalLine, int numSubPlotsinHhorizontalLine, int numClusterSubPlots, int distBetweenSubPlots, double bufferSize, boolean weightedSampling, File inputRasterFile ) throws Exception {

		// get strata out of SHP features 
		Stratum[] strata = ShapeFileUtilities.getStrata(inputShapeFile, sampleColumn, strataNames);		

		// reproject strata to UTM 
		Stratum[] strataUTM = CRSUtilities.convert2UTM(strata);
		
		// buffer Stratum geometries --> for Buffering, they have to be in UTM already 
		Stratum[] bufferedStrataUTM = applyBuffer(strataUTM, bufferSize);
		
		
		
		ArrayList<Plot> plots = null;
		// SRS
		if(samplingDesign == GUI_Designer.SIMPLE_RANDOM_SAMPLING && weightedSampling == false){
			plots = Sampling.simpleRandomSampling(bufferedStrataUTM, numPlotsToBeSampled);
		}
		
		// Systematic Sampling
		if(samplingDesign == GUI_Designer.SYSTEMATIC_SAMPLING){
			plots = Sampling.systematicSampling(bufferedStrataUTM, gridDistX, gridDistY, startPointSpecified, ShapeFileUtilities.getSHPCRS(inputShapeFile), startX, startY);
		}
		
		
		
		// Weighted Sampling (only makes sense for random Sampling)
		if(samplingDesign == GUI_Designer.SIMPLE_RANDOM_SAMPLING && weightedSampling == true){
			// call weightedSampling() using unmodified strata array!
			plots = weightedSampling(strata, inputRasterFile, numPlotsToBeSampled, bufferSize);
		}
		

		// Cluster Sampling --> call after Weighted Sampling so that weighted plots can be clustered as well
		if(clusterSampling == true){
			plots = clusterSampling(bufferedStrataUTM, plots, clusterShape, numClusterSubPlots, distBetweenSubPlots, numSubPlotsinHVerticalLine, numSubPlotsinHhorizontalLine);
		}
		
		return plots;
	}

	
	public static ArrayList<Plot> simpleRandomSampling(Stratum[] strata, int numPlots[]){
		ArrayList<Plot> plots = new ArrayList<Plot>();
		for (int i = 0; i < strata.length; i++) {
			plots.addAll(simpleRandomSampling(strata[i], numPlots[i]));
		}
		return plots;
	}
	
	
	
	/**
	 * Generates a specified number of output Plots in a given Stratum using simple random sampling. 
	 * @param stratum
	 * @param numPlots
	 * @return output plots
	 */
	public static ArrayList<Plot> simpleRandomSampling(Stratum stratum, int numPlots){
		// initialize output ArrayList
		ArrayList<Plot> plots = new ArrayList<Plot>();
		
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
		while(plots.size() < numPlots){
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
				plots.add(plot);
				plotNr ++; // increase plotNr only after successfully adding new Plot to output (eg, when randomly generated point falls inside BBOX)
			}
		}
		
		// return ArrayList
		return plots;
	}
	


	/**
	 * Creates a Point from X and Y values.
	 * @param x
	 * @param y
	 * @return
	 */
	public static Point createPoint(double x, double y){
		GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
		Coordinate coord = new Coordinate( x, y );
		Point point = geometryFactory.createPoint( coord );
		return point;
	}
	

	/**
	 * Delegates to systematicSampling() for each stratum.
	 * Startpoint is the same for each Stratum (if specified).
	 * @param strata
	 * @param gridDistX
	 * @param gridDistY
	 * @param startPointSpecified
	 * @param startPointCRS
	 * @param startX
	 * @param startY
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<Plot>  systematicSampling(Stratum[] strata, int gridDistX, int gridDistY, boolean startPointSpecified, CoordinateReferenceSystem startPointCRS, double startX, double startY)throws Exception{
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();
		for(int i = 0; i < strata.length; i++){
			ArrayList<Plot> plots = systematicSampling(strata[i],  gridDistX,  gridDistY,  startPointSpecified,  startPointCRS,  startX,  startY);
			outputPlots.addAll(plots);
		}
		return outputPlots;
	}
	
	
	/**
	 * startX and startX (if specified) come in the input SHP´s CRS whereas the strata 
	 * should already be reprojected to an UTM CRS.
	 * @param stratum
	 * @param gridDistX
	 * @param gridDistY
	 * @param startPointSpecified
	 * @param startPointCRS the CRS the startPoint Coordinates come in (should be the input shapefile´s CRS)
	 * @param startX
	 * @param startY
	 * @return
	 */
	public static ArrayList<Plot>  systematicSampling(Stratum stratum, int gridDistX, int gridDistY, boolean startPointSpecified, CoordinateReferenceSystem startPointCRS, double startX, double startY)throws Exception{
		Point startPoint = null;
		// either startPoint is already specified...
		if(startPointSpecified == true){
			// make a Point out of given start point coords
			startPoint = createPoint(startX, startY);
			// reproject startPoint using a MathTransform that transforms between sourceCRS and targetCRS
			// source CRS: we use the already defined input file's CRS here as the grid start point refers to a location within the input file area
			// target CRS: use the current stratum's CRS (as grid points are to be located inside this stratum, even if the starting point is located oustide the stratum)
			CoordinateReferenceSystem targetCRS = stratum.getCRS();
			// transform source CRS directly (without looking up the corresponding EPSG code first and use that CRS instead)
			MathTransform transform = CRS.findMathTransform(startPointCRS, targetCRS, true); // last param is the "lenient" param which can be important when there is not much transform info (?)
			// convert Point to target CRS
			startPoint = (Point) JTS.transform( startPoint, transform);
		}
		// or we create one using simple random sampling
		else{
			// call SRS to create startPoint
			ArrayList<Plot> plot = simpleRandomSampling(stratum, 1);
			startPoint =  plot.get(0).getPoint();
		}
		return systematicSampling(stratum, startPoint, gridDistX, gridDistY); 
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
	
	
	/**
	 * Applies an inwardly (i.e. negative) buffer to the Stratum.
	 * The input stratum has to be in a metric CRS (e.g. UTM), otherwise the method will likely
	 * produce an empty Geometry.
	 * 
	 * This method calls applyBuffer() for each Stratum in the array.
	 * @param stratum
	 * @param bufferSize in Meters
	 * @return
	 */
	public static Stratum[] applyBuffer(Stratum[] strata, double bufferSize){
		for(int i = 0; i < strata.length; i++){
			strata[i] = applyBuffer(strata[i], bufferSize);
		}
		return strata;
	}
	
	
	/**
	 * Applies an inwardly (i.e. negative) buffer to the Stratum.
	 * The input stratum has to be in a metric CRS (e.g. UTM), otherwise the method will likely
	 * produce an empty Geometry.
	 * @param stratum
	 * @param bufferSize in Meters
	 * @return
	 */
	public static Stratum applyBuffer(Stratum stratum, double bufferSize){
		Geometry geom = stratum.getGeometry();
		geom= geom.buffer(-bufferSize);
		stratum.setGeometry(geom);
		return stratum;
	}
	
	public static ArrayList<Plot> clusterSampling(Stratum[] strata, ArrayList<Plot> plots, int clusterShape, int numClusterSubPlots, int distBetweenSubPlots, int numSubPlotsinHVerticalLine, int numSubPlotsinHhorizontalLine){
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();
		if(clusterShape == GUI_Designer.I_SHAPE){
			outputPlots = Clusters.create_I_clusters(plots, distBetweenSubPlots, numClusterSubPlots, strata);
		}
		if(clusterShape == GUI_Designer.L_SHAPE){
			outputPlots = Clusters.create_L_clusters(plots, distBetweenSubPlots, numClusterSubPlots, 1, strata);
		}
		if(clusterShape == GUI_Designer.L_SHAPE_UPSIDE_DOWN){
			outputPlots = Clusters.create_L_clusters(plots, distBetweenSubPlots, numClusterSubPlots, -1, strata);
		}
		if(clusterShape == GUI_Designer.SQUARE_SHAPE){
			outputPlots = Clusters.create_Square_clusters(plots, distBetweenSubPlots, numClusterSubPlots, strata);
		}
		if(clusterShape == GUI_Designer.SQUARE_SHAPE_ROTATED){
			outputPlots = Clusters.create_rotated_Square_clusters(plots, distBetweenSubPlots, strata);
		}
		if(clusterShape == GUI_Designer.H_SHAPE){
			outputPlots = Clusters.create_H_clusters(plots, distBetweenSubPlots, numSubPlotsinHVerticalLine,  numSubPlotsinHhorizontalLine, strata);
		}
		return outputPlots;
	}
	
	/**
	 * Creates weighted sample plots. 
	 * This method can only be used for simple random sampling (not for Systematic Sampling). 
	 * The method will buffer any input stratum before generating the actual plots, so
	 * that buffered input strata will undergo the buffering process again.
	 * 
	 * This method calls weightedSampling() for every Stratum in the input array.
	 * 
	 * @param strata The strata to sample plots inside. 
	 * Must be unbuffered in order to find the correct zonalMax value.
	 * @param inputRasterFile
	 * @param numPlotsToBeSampled This array must be of the same length as the strata array. 
	 * @param bufferSize in Meters. This parameter is needed so that plots are generated 
	 * at a sufficient distance to the Stratum boundary (so that circular plots will not be cut off).
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<Plot> weightedSampling(Stratum[] strata, File inputRasterFile, int[] numPlotsToBeSampled, double bufferSize) throws Exception{
		ArrayList<Plot> plots = new ArrayList<Plot>();
		for(int i = 0; i < strata.length; i++){
			plots.addAll(weightedSampling(strata[i], inputRasterFile, numPlotsToBeSampled[i], bufferSize));
		}
		return plots;
	}
	
	/**
	 * Creates weighted sample plots. 
	 * This method can only be used for simple random sampling (not for Systematic Sampling). 
	 * The method will buffer any input stratum before generating the actual plots, so
	 * that buffered input strata will undergo the buffering process again.
	 * @param stratum The stratum to sample plots inside. Must be unbuffered in order to find the correct zonalMax value.
	 * @param inputRasterFile
	 * @param numPlotsToBeSampled
	 * @param bufferSize in Meters. This parameter is needed so that plots are generated 
	 * at a sufficient distance to the Stratum boundary (so that circular plots will not be cut off).
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<Plot> weightedSampling(Stratum stratum, File inputRasterFile, int numPlotsToBeSampled, double bufferSize) throws Exception{
		GridCoverage2D coverage = RasterProcessing.readGeoTiff(inputRasterFile);
		Geometry clipGeom = stratum.getGeometry();

		// equalize crsStratum and crsRaster
		// check if CRS are different
		CoordinateReferenceSystem crsStratum = stratum.getCRS();
		CoordinateReferenceSystem crsRaster = coverage.getCoordinateReferenceSystem();
		boolean needsReproject = !CRS.equalsIgnoreMetadata(crsStratum, crsRaster);
		if (needsReproject) {
			MathTransform transform = CRS.findMathTransform(crsStratum, crsRaster, true);
			clipGeom = JTS.transform(clipGeom, transform);
		}
		
		GridCoverage2D clippedCoverage = RasterProcessing.getClippedCoverage(clipGeom, coverage);
		double maxValue = RasterProcessing.getCoverageMaxValue(clippedCoverage, 0);
		double noDataValue = RasterProcessing.getNoDataValue(coverage, 0);

		// rejection testing: call  simpleRandomSampling() one at a time until desired number of plots is reached
		/*
		 * look at simpleRandomSampling() implementation more closely 
		 * --> reduce processing costs by not calculating a stratum BBOX for every single plot etc.
		 */

		// buffer input stratum before sampling plots --> reproject stratum to UTM before, otherwise the buffer cannot be applied
		Stratum stratumUTM = CRSUtilities.convert2UTM(stratum);
		Stratum bufferedStratum = Sampling.applyBuffer(stratumUTM, bufferSize);
		ArrayList<Plot> weightedSimpleRandomPlots = new ArrayList<Plot>();
		int numSampledPlots = 0;
		do{
			// create a Plot (one at a time, has to be rejection tested)
			Plot plot = Sampling.simpleRandomSampling(bufferedStratum, 1).get(0);
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
		while(numSampledPlots < numPlotsToBeSampled);
		
		return weightedSimpleRandomPlots;
	}


	

}
