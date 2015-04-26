/*
 * Ziel: Clip Raster by Polygon
 * exakter: Clip Geotiff (Raster) by SimpleFeature (Polygon) aus Shapefile
 * 
 * was geht: Shapefile einlesen, auf Polygon zugreifen, auf Raster zugreifen
 * 
 * 
 * Ansätze: 
 *   1) CoverageProcessor --> CoverageCrop Operation --> ROI param: clip Polygon
 *
*
 * Hier sind die Möglichkeiten, die ich für die aktuelle Problemstellung
 * "Clip Raster by Polygon"
 * herausgefunden habe. Vl gibts noch mehr, ich suche aber jetzt nicht mehr weiter.
 * 
 * org.geotools.process.raster.CropCoverage // --> ist das nicht das gleiche wie die CoverageCrop operation über CoverageProcessor?
 * org.geotools.process.raster.RasterZonalStatistics
 * org.geotools.coverage.processing.operation.ZonalStats
 * org.geotools.process.vector.ClipProcess //--> geht glaube ich nur für Polygone gegenseitig, nicht mit Raster
 * 
 * 
 * Usage: 
 * - Process Tutorial --> nicht so besonders informativ
 * - zu jeder Klasse gibts eine Testklasse mit Anwendungsbeispielen, da kann ich mir was abschauen
 * - sonst: Mailing Lists
 */


package test;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.SampleModel;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.geotools.coverage.Category;
import org.geotools.coverage.GridSampleDimension;
import org.geotools.coverage.TypeMap;
import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridCoverageFactory;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.grid.io.AbstractGridCoverage2DReader;
import org.geotools.coverage.grid.io.AbstractGridFormat;
import org.geotools.coverage.grid.io.GridCoverage2DReader;
import org.geotools.coverage.grid.io.GridFormatFinder;
import org.geotools.coverage.processing.AbstractOperation;
import org.geotools.coverage.processing.AbstractProcessor;
import org.geotools.coverage.processing.CoverageProcessor;
import org.geotools.coverage.processing.OperationJAI;
import org.geotools.coverage.processing.operation.Extrema;
import org.geotools.data.FeatureSource;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.gce.geotiff.GeoTiffFormat;
import org.geotools.gce.geotiff.GeoTiffReader;
import org.geotools.gce.geotiff.GeoTiffWriter;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.opengis.coverage.Coverage;
import org.opengis.coverage.SampleDimensionType;
import org.opengis.coverage.grid.GridCoverageReader;
import org.opengis.coverage.grid.GridCoverageWriter;
import org.opengis.coverage.grid.GridEnvelope;
import org.opengis.coverage.processing.Operation;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.geometry.DirectPosition;
import org.opengis.geometry.Envelope;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.datum.PixelInCell;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;
import org.opengis.util.InternationalString;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.geometry.GeneralEnvelope;
import org.geotools.process.ProcessFactory;
import org.geotools.process.Processors;
import org.geotools.process.raster.CoverageUtilities;
import org.geotools.process.raster.RasterProcessFactory;
import org.geotools.process.raster.RasterZonalStatistics;
import org.geotools.referencing.CRS;
import org.geotools.resources.coverage.IntersectUtils;
import org.geotools.util.NumberRange;
import org.jaitools.media.jai.zonalstats.ZonalStatsDescriptor;
import org.jaitools.numeric.SampleStats;

import production.SamplingFunctionalityMethods;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;

import javax.media.jai.ROI;
import javax.media.jai.ROIShape;
import javax.media.jai.iterator.RectIterFactory;
import javax.media.jai.iterator.RectIter;



public class RasterTests {


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	 * http://gis.stackexchange.com/questions/73210/how-to-crop-an-image-based-on-a-shapefile-using-geotools
	 * 
	 * it turns out that GeoTools does 
	 * indeed support the ability to crop a raster based on an arbitrary (....) geometry. 
	 * There are no good examples out there but a colleague passed me a tip suggesting the 
	 * CoverageProcessor and its "ROI" argument. // ROI: region of interest
	 * 
	 * Idee: CoverageProcessor --> CoverageCrop Operation --> ROI param: clip Polygon
	 */
	// TODO diees Methode gebrauchstüchtig machen --> siehe Anpassungen in main()	
	private static Coverage clipImageToFeatureSource(RenderedImage image, ReferencedEnvelope bounds,
			FeatureSource<SimpleFeatureType, SimpleFeature> featureSource)
					throws IOException, FactoryException, MismatchedDimensionException, TransformException {
		FeatureCollection<SimpleFeatureType, SimpleFeature> collection = featureSource
				.getFeatures();

		CoordinateReferenceSystem crsFeatures = featureSource.getSchema().getCoordinateReferenceSystem();
		CoordinateReferenceSystem crsRaster = bounds.getCoordinateReferenceSystem();
		boolean needsReproject = !CRS.equalsIgnoreMetadata(crsFeatures, crsRaster);
		MathTransform transform = CRS.findMathTransform(crsFeatures, crsRaster, true);

		FeatureIterator<SimpleFeature> iterator = collection.features();
		List<Geometry> all = new ArrayList<Geometry>();
		try {
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
				Geometry geometry = (Geometry) feature.getDefaultGeometry();
				if (geometry == null)
					continue;
				if (!geometry.isSimple())
					continue;
				if (needsReproject) {
					geometry = JTS.transform(geometry, transform);
					System.out.println("Reprojected a geometry.  Result is " + geometry.toString());
				}
				Geometry intersection = geometry.intersection(JTS.toGeometry(bounds)); // Intersection between feature Geometry and Raster BBox --> the part of the geometry that is inside the raster
				if (intersection.isEmpty()) {
					continue;
				}
				if(intersection instanceof MultiPolygon) {
					MultiPolygon mp = (MultiPolygon)intersection;
					for (int i = 0; i < mp.getNumGeometries(); i++) {
						com.vividsolutions.jts.geom.Polygon g = (com.vividsolutions.jts.geom.Polygon)mp.getGeometryN(i);
						Geometry gIntersection = IntersectUtils.intersection(g, JTS.toGeometry(bounds));
						if (gIntersection.isEmpty()) {
							continue;
						}
						all.add(g);
					}
				}
				else if (intersection instanceof Polygon)
					all.add(intersection);
				else
					continue;
			}
		} finally {
			if (iterator != null) {
				iterator.close();
			}
		}
		GridCoverageFactory gridCoverageFactory = new GridCoverageFactory();
		Coverage coverage = gridCoverageFactory.create("Raster", image, bounds); // hier Fehler: Illegal transform of type "LinearTransform1D".
		Coverage clippedCoverage = null;
		if (all.size() > 0) {
			CoverageProcessor processor = new CoverageProcessor();
			ParameterValueGroup params = processor.getOperation("CoverageCrop")
					.getParameters();
			params.parameter("Source").setValue(coverage);
			GeometryFactory factory = JTSFactoryFinder.getGeometryFactory(null);
			Geometry[] a = all.toArray(new Geometry[0]);
			GeometryCollection c = new GeometryCollection(a, factory);
			//params.parameter("ENVELOPE").setValue(bounds);
			params.parameter("ROI").setValue(c);
			params.parameter("ForceMosaic").setValue(true);
			clippedCoverage = processor.doOperation(params);
		}
		if (all.size() == 0){
			//logger.info("Crop by shapefile requested but no simple features matched extent!");
			System.out.println("Crop by shapefile requested but no simple features matched extent!");
		}
		return clippedCoverage;
	}


	/*
	 * http://gis.stackexchange.com/questions/73210/how-to-crop-an-image-based-on-a-shapefile-using-geotools
	 * 
	 * This clips the image but may also reduce the original extent. 
	 * If you want to preserve that, then you'll need to "matt" the clippedCoverage like this:
	 */
	// TODO rausfinden: wofür ist diese Methode gut???	
	private BufferedImage mattCroppedImage(final BufferedImage source, GridCoverage2D cropped) 
	{
		RenderedImage raster = cropped.getRenderedImage();
		int height = source.getHeight();
		int width = source.getWidth();
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D gr = image.createGraphics();
		gr.setPaint(Color.green);
		gr.fill(new Rectangle2D.Double(0,0, image.getWidth(), image.getHeight()));
		AffineTransform at = AffineTransform.getTranslateInstance(0, 0);
		gr.drawRenderedImage(cropped.getRenderedImage(), at);
		return image;
	}
	
	


	/*
	 * Find raster pixel value at given location (x, y)
	 * 
	 * wenn dieser Teil nicht als eigene Methode ausgelagert wird,
	 * erkennt er nicht mehr das richtige Format (WorldImageFormat statt GeoTiff)
	 */
	public static Object getValueAtPosition(GridCoverage2D coverage, double x, double y){
		/* You can also evaluate the coverage at a specific point in order to determine 
		 * what values are present at that location.:
		 * http://docs.geotools.org/latest/userguide/library/coverage/grid.html
		 */
		// direct access
		CoordinateReferenceSystem crs = coverage.getCoordinateReferenceSystem();
		DirectPosition position = new DirectPosition2D( crs, x, y);

		double[] value = (double[]) coverage.evaluate( position ); // assume double --> return type kann sich ändern, je nach Werten im Kanal
		//		byte[] value = (byte[]) coverage.evaluate( position ); // für clouds.jpg gibts einene Byte-Array zurück
		return value;
	}



	
	/*
	 * (manipulated) Code snippet from
	 * https://svn.osgeo.org/geotools/trunk/modules/plugin/geotiff/src/test/java/org/geotools/gce/geotiff/GeoTiffWriterTest.java 
	 */
	public static void writeGeoTiff(String outputFilePath, GridCoverage2D coverage)throws IOException{
		File outputFile = new File(outputFilePath);
		final GeoTiffFormat format = new GeoTiffFormat();
		final GridCoverageWriter writer = format.getWriter(outputFile);
		writer.write(coverage, null); // kein sidecar file wird geschrieben, korrektes CRS ist anscheinend trotzdem im Output-GeoTiff drin (sagt QGIS)
	}



	public static void main(String[] args)throws Exception{

		// use fictive Ghana weight Raster (GeoTiff format)
		File rasterFile = new File("D:\\_HCU\\_Masterarbeit\\TestData\\Globcover2009_V2.3_Global\\GLOBCOVER_L4_200901_200912_V2.3.tif");
		// File rasterFile = new File("D:\\_HCU\\_Masterarbeit\\TestData\\uDig_Sample_Data\\clouds.jpg");
		System.out.println("input raster: " + rasterFile.toString());
		
		// automatische Rasterformaterkennung --> ist mir zu störanfällig, geht theoretisch auch 
		//AbstractGridFormat format = GridFormatFinder.findFormat( rasterFile ); // fehleranfällig !!
		//AbstractGridCoverage2DReader reader = format.getReader(rasterFile);
		//GridCoverage2DReader reader = format.getReader(rasterFile); // macht das einen Unterschied? --> hier nicht

		// create Reader directly: scheint weniger fehleranfällig zu sein,  
		GeoTiffReader reader = new GeoTiffReader(rasterFile);

		GridCoverage2D coverage = null;
		// GridCoverage2D coverage = reader.read(null); // direkt-ohne Exception handling
		try {
			coverage = reader.read(null);
		} catch (IOException giveUp) {
			throw new RuntimeException(giveUp);
		}
		// int numBands = coverage.getNumSampleDimensions();

		// GridSampleDimension dim = coverage.getSampleDimension(0); // SampleDimension: Band; should always be index position 0 as there is only one band in a weight raster

		// String sampleDimensionName = dim.getDescription().toString();

		// GridGeometry2D gridGeom = coverage.getGridGeometry();

		
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Clip Polygon from Raster 
		
		// Specify Input Shapefile and stratum
		File shapeFile = new File("D:\\_HCU\\_Masterarbeit\\TestData\\TestSmallPolygons\\TestSmallPolygons.shp");
		String sampleColumn = "VEGZONE";
		String clipPolygon = "Test4";
		
		System.out.println("shapeFile: " + shapeFile.toString());
		System.out.println("Stratum: " + clipPolygon.toString());
		
		
		
		
		// connect to Shapefile
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(shapeFile);
		SimpleFeatureSource featureSource = dataStore.getFeatureSource(); // wird benötigt, um an einzelne Features ranzukommen (Feature = Zeile in SHP Attribute Table)
		SimpleFeatureCollection collection = featureSource.getFeatures();

		// check if CRS from SHP and Raster are different
		CoordinateReferenceSystem crsFeatures = featureSource.getSchema().getCoordinateReferenceSystem();
		CoordinateReferenceSystem crsRaster = coverage.getCoordinateReferenceSystem();
		boolean needsReproject = !CRS.equalsIgnoreMetadata(crsFeatures, crsRaster);
		MathTransform transform = CRS.findMathTransform(crsFeatures, crsRaster, true);

		//TODO weiter vereinfachen
		// iterate over features inside SHP
		SimpleFeatureIterator iterator = collection.features();
		List<Geometry> clipGeometries = new ArrayList<Geometry>();
		try {
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
				// check for desired clipPolygon and do the rest of stuff only if it is found
				if(feature.getAttribute(sampleColumn).equals(clipPolygon))
				{
					// extract and reproject Geometry
					Geometry geometry = (Geometry) feature.getDefaultGeometry();
					if (needsReproject) {
						geometry = JTS.transform(geometry, transform);
					}
					
					clipGeometries.add(geometry);
					
//					//TODO intersection-part rausschmeißen
//					Envelope envelope = coverage.getEnvelope();
//					ReferencedEnvelope refEnv = new ReferencedEnvelope(envelope);
//					
//					Geometry intersection = geometry.intersection(JTS.toGeometry(refEnv)); // intersection Polygonfläche (feature) mit Raster BBox
//					if (intersection.isEmpty()) {
//						continue;
//					}
//					
//					if(intersection instanceof MultiPolygon) {
//						MultiPolygon mp = (MultiPolygon)intersection;
//						for (int i = 0; i < mp.getNumGeometries(); i++) {
//							com.vividsolutions.jts.geom.Polygon g = (com.vividsolutions.jts.geom.Polygon)mp.getGeometryN(i);
//							Geometry gIntersection = IntersectUtils.intersection(g, JTS.toGeometry(refEnv));
//							if (gIntersection.isEmpty()) {
//								continue;
//							}
//							clipGeometries.add(g);
//						}
//					}
//					else if (intersection instanceof Polygon)
//						clipGeometries.add(intersection);
				}
				
				else
					continue;
			}
		} finally {
			if (iterator != null) {
				iterator.close();
			}
		}
		
		
		GridCoverage2D clippedCoverage = null;
		if (clipGeometries.size() > 0) {
			CoverageProcessor processor = new CoverageProcessor();
			Operation operation = processor.getOperation("CoverageCrop");
			ParameterValueGroup params = operation.getParameters();
			params.parameter("Source").setValue(coverage);
			GeometryFactory factory = JTSFactoryFinder.getGeometryFactory(null);
			Geometry[] a = clipGeometries.toArray(new Geometry[0]);
			GeometryCollection c = new GeometryCollection(a, factory);
			params.parameter("ROI").setValue(c);
			params.parameter("ForceMosaic").setValue(true); // was macht ForceMosaic?
			clippedCoverage = (GridCoverage2D)processor.doOperation(params);
		}

		
		
		// write clippedCoverage to file and ceck in QGIS
		//writeGeoTiff("D:\\_HCU\\_Masterarbeit\\TestData\\TestClip\\TestSmallPolygon4.tif", clippedCoverage);
		
		
		
		
		/////////////////////////////////////////////////////////////////////////////////////
		// Iterieren über Rasterbild:                                                      //
		// Vergleich von 3 Varianten                                                       //
        /////////////////////////////////////////////////////////////////////////////////////
		
		
		// ignore null value when screening raster for max and sums
		// find null value in Categories
		// TODO find null value code must be generalized to fit other input files
		GridSampleDimension band = clippedCoverage.getSampleDimension(0); // weight raster is supposed to only have 1 dimension (band)
		List<Category> categories = band.getCategories();
		System.out.println("categories: " + categories.toString());
		Category noData = categories.iterator().next();
		NumberRange<? extends Number> range= noData.getRange();
		double noDataValue = range.getMaximum(); // geht auch mit getMinimum(), da nur ein Wert
		System.out.println(noData.toString());
		
        ////////////////////////////////////////////////////////
		//            Variante 1: using RectIter              //
        ////////////////////////////////////////////////////////
		// maxValue von clippedCoverage
		// Lösung: org.geotools.coverage.grid.RenderedSampleDimension ermittelt min/max --> Code snippet kopiert und angepasst
        /*
         * Code copied from: org.geotools.coverage.grid.RenderedSampleDimension, Zeile 283 ff
         * Computes the minimal and maximal values, if not explicitely provided.
         * This information is required for determining the range of geophysics
         * values.
         * 
         * iterator: javax.media.jai.iterator.RectIter
         * runtime class: com.sun.media.jai.iterator.RectIterFallback
         * 
         * min/ max sind double[]
         */
		
		// Zugriff auf Raster-Bild über RenderedImage
		RenderedImage renderedImage = clippedCoverage.getRenderedImage(); // runtime class of renderedImage: javax.media.jai.RenderedOp
		double[] min = null;
		double[] max = null;
		double suma = 0;
		int numBands = clippedCoverage.getNumSampleDimensions();
		RectIter coverageIterator = RectIterFactory.create(renderedImage, null);
		min = new double[numBands];
		Arrays.fill(min, Double.POSITIVE_INFINITY);
		max = new double[numBands];
		Arrays.fill(max, Double.NEGATIVE_INFINITY);
		int b = 0;
		coverageIterator.startBands(); // Sets the iterator to the first band of the image. The pixel column and line are unchanged. 
		if (!coverageIterator.finishedBands()) do { // Returns true if the max band in the image has been exceeded. 
			coverageIterator.startLines(); // Sets the iterator to the first line of its bounding rectangle. The pixel and band offsets are unchanged. 
			if (!coverageIterator.finishedLines()) do { // Returns true if the bottom row of the bounding rectangle has been passed. 
				coverageIterator.startPixels(); // Sets the iterator to the leftmost pixel of its bounding rectangle. The line and band offsets are unchanged. 
				if (!coverageIterator.finishedPixels()) do { // Returns true if the right edge of the bounding rectangle has been passed. 
					final double currentValue = coverageIterator.getSampleDouble();// Returns the current sample as a double. 
					if(currentValue != noDataValue){ // ignore NULL values (might be different according to input file)
						if (currentValue<min[b]) min[b]=currentValue;
						if (currentValue>max[b]) max[b]=currentValue;
						suma = suma + currentValue;
					}

				} while (!coverageIterator.nextPixelDone());
			} while (!coverageIterator.nextLineDone());
			if (!(min[b] < max[b])) {
				min[b] = 0;
				max[b] = 1;
			}
			b++;
		} while (!coverageIterator.nextBandDone());

		System.out.println("min mit RectIter: " + min[0]);
		System.out.println("max mit RectIter: " + max[0]);
		System.out.println("Summe mit RectIter: " + suma);


		
        ////////////////////////////////////////////////////////
		//            Variante 2: Eigenbau                    //
        ////////////////////////////////////////////////////////
		// Bsp: http://de.slideshare.net/moovida/opensource-gis-development-part-4, slide 7

		// könnte evtl. Probleme geben, weil das Raster nicht bei 0,0 losgeht, sondern bei 8,118 oder so
		// difference between image space and world space
		// beware of NODATA values!
		
		// über renderedImage
		//numColumns
		int numColumns = renderedImage.getWidth();
		//numRows
		int numRows = renderedImage.getHeight();
		int minX = renderedImage.getMinX();
		int minY = renderedImage.getMinY();
		System.out.println("renderedImage width: " + numColumns);
		System.out.println("renderedImage height: " + numRows);
		System.out.println("renderedImage minX: " + minX);
		System.out.println("renderedImage minY: " + minY);
		
		

		
		// how to access pixel values
		// raster.getSample()
		// coverage.evaluate()
		
		double sum = 0;
		double maxValue = 0;
		double[] values = new double[1];
		double currentValue = 0; 
		for(int r = minY; r < minY + numRows; r++){ // loop over rows
			for(int c = minX; c < minX + numColumns; c++){ // loo over columns
				clippedCoverage.evaluate(new GridCoordinates2D(c, r), values);
				currentValue = values[0];
				if(currentValue != noDataValue){ // ignore NULL values (might be different according to input file)
					sum = sum + currentValue;
					if(currentValue > maxValue) maxValue = currentValue;
				}

			}
		}
		System.out.println("max mit Eigenbau: " + maxValue);
		System.out.println("sum mit Eigenbau: " + sum);

		
		
		
        ////////////////////////////////////////////////////////
		//            Variante 3: Internetvorlage             //
        ////////////////////////////////////////////////////////
		// mal ausprobieren als Alternative
		// so gehts anscheinend, aber unschöne Lösung
		Raster raster = clippedCoverage.getRenderedImage().getData();
        double[] data = new double[raster.getHeight()*raster.getWidth()];        
        raster.getSamples(raster.getMinX(),
                raster.getMinY(),
                raster.getWidth(), 
                raster.getHeight(), 0, data);
        // convenience method to find max in double[]-array
        // braucht import org.apache.commons.lang3.ArrayUtils;
        //System.out.println(SampleStats.max(ArrayUtils.toObject(data), true));
		// ich machs hier selbst for speed: loop over array and find max
        double maximum = 0;
        double summe = 0;
        double currentVal;
        for(int i = 0; i < data.length; i++){
        	currentVal = data[i];
        	if(currentVal != noDataValue){ // ignore NULL values (might be different according to input file)
        		if(currentVal > maximum) maximum = currentVal;
        		summe = summe + currentVal;
        	}

		}
		System.out.println("max mit Internetvorlage: " + maximum);
		System.out.println("sum mit Internetvorlage: " + summe);
		
		
		
		// is there a NULL value in the clipped region? then throw exception
		// --> schwierig rauszufinden, weil das clippedraster rechteckig ist und deswegen ganz viele Nullwerte hat
		
		// Summe aller Pixelwerte im Clip --> RectIter einsetzbar?
		
		
		
		

		
		
		
		
		
		
		
		
		
		
		
		
		
//		// kucken, was der CoverageProcessor alles kann
//		CoverageProcessor cProcess = new CoverageProcessor();
//		Collection<Operation> operations = cProcess.getOperations();
//		System.out.println("All available operations forCoverageCrop operation:  ");
//		Iterator<Operation> i = operations.iterator();
//		while(i.hasNext()){
//			System.out.println(i.next().toString());
//		}
//		Operation operation = cProcess.getOperation("CoverageCrop");
//		System.out.println("Description: " + operation.getDescription());
//		System.out.println("Name: " + operation.getName());
//		System.out.println("Vendor: " + operation.getVendor());
//		System.out.println("Parameters: " + operation.getParameters());
//		System.out.println("runtime Class: " + operation.getClass());
//		// org.geotools.coverage.processing.operation.Crop // CoverageCrop class
		
		
		
		// test this:  	
//		Set<ProcessFactory> processFactories = org.geotools.process.Processors.getProcessFactories();
//		Iterator<ProcessFactory> i = processFactories.iterator();
//		while(i.hasNext()){
//			ProcessFactory pf = i.next();
//			System.out.println(pf.toString());
//			System.out.println(pf.getNames());
//			
//		}
		//org.geotools.process.Processors.getParameterInfo()
		//org.geotools.process.Processors.getResultInfo()


		////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//------------------------------------------------------------------------------------------------------
		// Find raster pixel value at given location (x, y)
		//double[] result = (double[]) getValueAtPosition(coverage,653672.0, 738800.0 ); // je nach inputfile können unterschiedliche Arrayklassen rauskommen
		//System.out.println("Raster cell Value at specified position: " + result[0]);
		//------------------------------------------------------------------------------------------------------


//		//------------------------------------------------------------------------------------------------------
//		/*
//		 * crop --> brauche ich vl, um Stratumflächen aus dem Raster auszuschneiden
//		 * 
//		 * The Crop operation provides a way to crop (cut) a GridCoverage in order to obtain 
//		 * a sub-area defined by an Envelope.:
//		 * http://docs.geotools.org/latest/userguide/library/coverage/grid.html
//		 */
//		//final AbstractProcessor processor = new DefaultProcessor(null);
//		final CoverageProcessor processor = new CoverageProcessor(null);
//
////		// list operations available through CoverageProcessor  
////		Collection<Operation> operations = processor.getOperations();
////		Iterator iterator = operations.iterator();
////		while(iterator.hasNext()){
////			Operation operation = (Operation)iterator.next();
////			System.out.println(operation.toString());
////		}
//		
//		// Test operations Warp and ZonalStats
//
////		// get some specific operation and describe it
////		Operation op = processor.getOperation("Warp");
////		System.out.println(op.getDescription());
//		
//		final ParameterValueGroup param = processor.getOperation("CoverageCrop").getParameters();
//		
//		//GridCoverage2D coverage = ...{get a coverage from somewhere}...;
//		final GeneralEnvelope crop = new GeneralEnvelope(new double[]{731048.0, 582189.0}, new double[]{800440.0, 630292.0});
//		param.parameter("Source").setValue( coverage );
//		param.parameter("Envelope").setValue( crop );
//
//		System.out.println("Bis hierher alles gut");
//
//		GridCoverage2D cropped = (GridCoverage2D) processor.doOperation(param);
//		//GridCoverageWriter writer = format.getWriter(new File("D:\\_HCU\\_Masterarbeit\\TestData\\goldengod.tif"));
//		//writer.write(cropped, null);
//		//------------------------------------------------------------------------------------------------------











		//------------------------------------------------------------------------------------------------------ 
		//        // Code Snippet
		//        // https://github.com/geotools/geotools/blob/master/modules/library/coverage/src/test/java/org/geotools/coverage/processing/ZonalStasTest.java
		//        /*
		//         * crop on region of interest
		//         */
		//        final CoverageProcessor processor = CoverageProcessor.getInstance();
		//        final ParameterValueGroup param = processor.getOperation("CoverageCrop")
		//                .getParameters();
		//        param.parameter("Source").setValue(gridCoverage2D);
		//        param.parameter("Envelope").setValue(new GeneralEnvelope(bbox));
		//        final GridCoverage2D cropped = (GridCoverage2D) processor.doOperation(param);
		//------------------------------------------------------------------------------------------------------

	}

	//	//------------------------------------------------------------------------------------------------------
	//	/*
	//	 * http://stackoverflow.com/questions/27761590/geotools-total-area-where-gridcoverage-has-value-x
	//	 * 
	//	 * I'm working in Geotools with Java. So far, I have a GridCoverage2D and a List of Geometries. 
	//	 * The GridCoverage2D is a digital elevation model, originating from a geotiff. Everything 
	//	 * works fine till here.
	//	 * Now I want to get the area for each polygon where the elevation has a certain value. 
	//	 * For example for this Geometry, I want to know the total area where the elevation is 27 m.
	//	 * How can I achieve this?
	//	 * 
	//	 *   
	//	 *   
	//	 */
	//	private double selectCells(GridCoverage2D cov, int value) {
	//		GridGeometry2D geom = cov.getGridGeometry();
	//		// cov.show();
	//		final OperationJAI op = new OperationJAI("Histogram");
	//		ParameterValueGroup params = op.getParameters();
	//		GridCoverage2D coverage;
	//		params.parameter("Source").setValue(cov);
	//		coverage = (GridCoverage2D) op.doOperation(params, null);
	//		javax.media.jai.Histogram hist = (Histogram) coverage
	//				.getProperty("histogram");
	//
	//		int total = hist.getSubTotal(0, value, value);
	//		double area = calcAreaOfCell(geom);
	//		Unit<?> unit = cov.getCoordinateReferenceSystem().getCoordinateSystem()
	//				.getAxis(0).getUnit();
	//		System.out.println("which gives " + (area * total) + " " + unit
	//				+ "^2 area with value " + value);
	//		return area * total;
	//	}
	//
	//	private double calcAreaOfCell(GridGeometry2D geom) {
	//		GridEnvelope gridRange = geom.getGridRange();
	//		int width = gridRange.getHigh(0) - gridRange.getLow(0) + 1; // allow for the
	//		int height = gridRange.getHigh(1) - gridRange.getLow(1) + 1;// 0th row/col
	//		Envelope envelope = geom.getEnvelope();
	//		double dWidth = envelope.getMaximum(0) - envelope.getMinimum(0);
	//		double dHeight = envelope.getMaximum(1) - envelope.getMinimum(1);
	//		double cellWidth = dWidth / width;
	//		double cellHeight = dHeight / height;
	//
	//		return cellWidth * cellHeight;
	//	}
	//------------------------------------------------------------------------------------------------------
}
