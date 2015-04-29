package test;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
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
import org.opengis.coverage.CannotEvaluateException;
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

import production.Plot;
import production.Sampling;

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



public class RasterTests_2 {



	/**
	 * Find raster pixel value at given Plot location.
	 * This method is overloaded so that it fits both raster files with integer as well as floating point values.
	 * - works with multiband rasters too
	 * http://docs.geotools.org/latest/userguide/library/coverage/grid.html
	 * @param coverage
	 * @param plot
	 * @param dest
	 * @return
	 * @throws Exception
	 */
	public static byte[] getValueAtPosition(GridCoverage2D coverage, Plot plot, byte[] dest) throws Exception{
		Point point = plot.getPoint();
		// check if CRS from Plot and Raster are different and reproject if needed
		CoordinateReferenceSystem crsPlot = plot.getCRS();
		CoordinateReferenceSystem crsRaster = coverage.getCoordinateReferenceSystem();
		boolean needsReproject = !CRS.equalsIgnoreMetadata(crsPlot, crsRaster);
		MathTransform transform = CRS.findMathTransform(crsPlot, crsRaster, true);

		if(needsReproject){
			point = (Point)JTS.transform(point, transform);
		}
		
		DirectPosition position = new DirectPosition2D( crsRaster, point.getX(), point.getY());
		
		byte[] value =coverage.evaluate(position, (byte  []) null);
		return value;
	}


	/**
	 * Find raster pixel value at given Plot location.
	 * This method is overloaded so that it fits both raster files with integer as well as floating point values.
	 * - works with multiband rasters too
	 * http://docs.geotools.org/latest/userguide/library/coverage/grid.html
	 * @param coverage
	 * @param plot
	 * @param dest
	 * @return
	 * @throws Exception
	 */
	public static double[] getValueAtPosition(GridCoverage2D coverage, Plot plot, double[] dest) throws Exception{
		Point point = plot.getPoint();
		// check if CRS from Plot and Raster are different and reproject if needed
		CoordinateReferenceSystem crsPlot = plot.getCRS();
		CoordinateReferenceSystem crsRaster = coverage.getCoordinateReferenceSystem();
		boolean needsReproject = !CRS.equalsIgnoreMetadata(crsPlot, crsRaster);
		MathTransform transform = CRS.findMathTransform(crsPlot, crsRaster, true);

		if(needsReproject){
			point = (Point)JTS.transform(point, transform);
		}
		
		DirectPosition position = new DirectPosition2D( crsRaster, point.getX(), point.getY());
		
		double[] value =coverage.evaluate(position, (double[]  ) null);
		return value;
	}

	
	/**
	 * Writes a GridCoverage2D to the specified output file path.
	 * https://svn.osgeo.org/geotools/trunk/modules/plugin/geotiff/src/test/java/org/geotools/gce/geotiff/GeoTiffWriterTest.java 
	 * @param outputFilePath
	 * @param coverage
	 * @throws IOException
	 */
	public static void writeGeoTiff(String outputFilePath, GridCoverage2D coverage)throws IOException{
		File outputFile = new File(outputFilePath);
		final GeoTiffFormat format = new GeoTiffFormat();
		final GridCoverageWriter writer = format.getWriter(outputFile);
		writer.write(coverage, null); // no sidecar file; CRS is written inside .tif file
	}


	/**
	 * Gets the Geometries to be clipped from a Raster as an ArrayList object. 
	 * If the input SHP does not have the same CRS as the input raster,
	 * the Geometries are converted to the raster CRS.
	 * 
	 * @param coverage
	 * @param shapeFile
	 * @param sampleColumn
	 * @param clipStratum
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<Geometry> getClipGeometries (GridCoverage2D coverage, File shapeFile, String sampleColumn, String clipStratum) throws Exception{
		// connect to Shapefile
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(shapeFile);
		SimpleFeatureSource featureSource = dataStore.getFeatureSource(); // wird benötigt, um an einzelne Features ranzukommen (Feature = Zeile in SHP Attribute Table)
		SimpleFeatureCollection collection = featureSource.getFeatures();

		// check if CRS from SHP and Raster are different
		CoordinateReferenceSystem crsFeatures = featureSource.getSchema().getCoordinateReferenceSystem();
		CoordinateReferenceSystem crsRaster = coverage.getCoordinateReferenceSystem();
		boolean needsReproject = !CRS.equalsIgnoreMetadata(crsFeatures, crsRaster);
		MathTransform transform = CRS.findMathTransform(crsFeatures, crsRaster, true);

		// iterate over features inside SHP
		SimpleFeatureIterator iterator = collection.features();
		ArrayList<Geometry> clipGeometries = new ArrayList<Geometry>();
		try {
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
				// check for desired clipPolygon perform following operations only if it is found
				if(feature.getAttribute(sampleColumn).equals(clipStratum))
				{
					// extract Geometry
					Geometry geometry = (Geometry) feature.getDefaultGeometry();
					// reproject Geometry if needed
					if (needsReproject) {
						geometry = JTS.transform(geometry, transform);
					}
					clipGeometries.add(geometry);
				}
			}
		} finally {
			if (iterator != null) {
				iterator.close();
			}
		}
		return clipGeometries;
	}

	
	
	/**
	 * Clips an input GridCoverage2D using ArrayList<Geometry> clipGeometries.
	 * @param clipGeometries must have the same CRS as coverage. This is provided when 
	 * clipGeometries have been extracted using the getClipGeometries() method.
	 * @param coverage
	 * @return
	 */
	public static GridCoverage2D getClippedCoverage(ArrayList<Geometry> clipGeometries, GridCoverage2D coverage){
		CoverageProcessor processor = new CoverageProcessor();
		Operation operation = processor.getOperation("CoverageCrop");
		ParameterValueGroup params = operation.getParameters();
		params.parameter("Source").setValue(coverage);
		GeometryFactory factory = JTSFactoryFinder.getGeometryFactory(null);
		Geometry[] a = clipGeometries.toArray(new Geometry[0]);
		GeometryCollection c = new GeometryCollection(a, factory);
		params.parameter("ROI").setValue(c);
		params.parameter("ForceMosaic").setValue(true); 
		System.setProperty("com.sun.media.jai.disableMediaLib", "true"); // gets rid of the annoying Exception in the following line
		GridCoverage2D clippedCoverage = (GridCoverage2D)processor.doOperation(params); // this line throws the following annoying Exception if JAI MediaLib is not disabled: Error: Could not find mediaLib accelerator wrapper classes. Continuing in pure Java mode.
		return clippedCoverage;
	}

	
	/**
	 * Convenience method for reading GeoTIFF raster files.
	 * @param rasterFile
	 * @return
	 * @throws Exception
	 */
	public static GridCoverage2D readGeoTiff(File rasterFile) throws Exception{
		// create Reader  
		GeoTiffReader reader = new GeoTiffReader(rasterFile);
		// read Raster data
		GridCoverage2D coverage = reader.read(null);
		return coverage;
	}
	
	
	/**
	 * Gets the NODATA value of the specified band of a GridCoverage2D object. For use with weight rasters: bandIndex should always be 0
	 * as weight rasters are supposed to only have 1 band. If there is no category named "No data", return value will be
	 * Double.NEGATIVE_INFINITY.
	 * @param coverage input raster 
	 * @param bandIndex the band to get the NODATA value from
	 * @return noDataValue
	 * @throws Exception
	 */
	public static double getNoDataValue(GridCoverage2D coverage, int bandIndex) throws Exception{
		double noDataValue = Double.NEGATIVE_INFINITY; // is that a wise thing to do here? --> primitive variables cannot be set to NULL, so i need a dummy value that at the same time is unlikely to be used as a valid raster data value

		GridSampleDimension band = coverage.getSampleDimension(bandIndex); 
		List<Category> categories = band.getCategories();
		Iterator<Category> iterator = categories.iterator();
		while(iterator.hasNext()){
			Category category = iterator.next();
			InternationalString name = category.getName();

			if(name.toString().equals("No data")){
				NumberRange<? extends Number> range = category.getRange();

				if((double)range.getMinValue() != (double)range.getMaxValue()){ //"No data" category must only contain one value
					throw new Exception("Category \"No data\" contains more than 1 value");
				}else{
					noDataValue = (double)range.getMinValue(); //might as well use geMaxValue() as min must be equal to max
				}
				break;
			}
		}
		return noDataValue;
	}
	
	
	/**
	 * Gets the maximum value out of all not-NULL values in the specified band of the input GridCoverage2D object (NULL values will be ignored).
	 * If the operation somehow does not work correctly, the return value will be Double.NEGATIVE_INFINITY
	 * @param coverage
	 * @param band index of the band of the raster file that is to be used (should always be 0 in the case of weight rasters)
	 * @return maximum value
	 * @throws Exception 
	 */
	public static double getCoverageMaxValue(GridCoverage2D coverage, int band) throws Exception{
		// get renderedImage
		RenderedImage renderedImage = coverage.getRenderedImage();
		// get image params
		int numColumns = renderedImage.getWidth();
		int numRows = renderedImage.getHeight();
		int minX = renderedImage.getMinX();
		int minY = renderedImage.getMinY();
		double noDataValue = getNoDataValue(coverage, band);// raster noDataValue might differ from 0 and must therefore be specifically queried and ignored
		
		double maxValue = Double.NEGATIVE_INFINITY;
		double[] values = new double[1];
		double currentValue = 0; 
		
		// iterate over renderedImage using image params
		for(int r = minY; r < minY + numRows; r++){ // loop over rows
			for(int c = minX; c < minX + numColumns; c++){ // loo over columns
				coverage.evaluate(new GridCoordinates2D(c, r), values);
				currentValue = values[0];
				if(currentValue != noDataValue){ // ignore NULL values (might be different according to input file)
					if(currentValue > maxValue) maxValue = currentValue;
				}

			}
		}
		return maxValue;
	}	
	

	/**
	 * Calculates the sum of all not-NULL values in the specified band of the input GridCoverage2D object(NULL values will be ignored).
	 * @param coverage
	 * @param band index of the band of the raster file that is to be used (should always be 0 in the case of weight rasters)
	 * @return sum
	 * @throws Exception 
	 */
	public static double getCoverageSum(GridCoverage2D coverage, int band) throws Exception{
		// get renderedImage
		RenderedImage renderedImage = coverage.getRenderedImage();
		// get image params
		int numColumns = renderedImage.getWidth();
		int numRows = renderedImage.getHeight();
		int minX = renderedImage.getMinX();
		int minY = renderedImage.getMinY();
		double noDataValue = getNoDataValue(coverage, band); // raster noDataValue might differ from 0 and must therefore be specifically queried and ignored

		double sum = 0;
		double[] values = new double[1];
		double currentValue = 0; 
		
		// iterate over renderedImage using image params
		for(int r = minY; r < minY + numRows; r++){ // loop over rows
			for(int c = minX; c < minX + numColumns; c++){ // loo over columns
				coverage.evaluate(new GridCoordinates2D(c, r), values);
				currentValue = values[0];
				if(currentValue != noDataValue){ // ignore NULL values (might be different according to input file)
					sum = sum + currentValue;
				}
			}
		}
		return sum;
	}

	
	/**
	 * Applies a rejection testing operation to a normalizedPlotWeight value in order to 
	 * decide on whether to keep a sampled plot in the sample.
	 * The normalizedPlotWeight Parameter is compared to a uniformly distributed random number.
	 * If normalizedPlotWeight < random number, return value will be false, otherwise true.
	 * For the uniformly distributed random number, Math.random() is used; the method Javadoc says that 
	 * "values are chosen pseudorandomly with (approximately) uniform distribution"
	 * @param normalizedPlotWeight: must be in range 0-1
	 * @return
	 */
	public static boolean rejectionTesting(double normalizedPlotWeight){
		boolean keepPlot = true;
		// use uniformly distributed numbers
		//Math.random() documentation: values are chosen pseudorandomly with (approximately) uniform distribution 
		double randomNumber = Math.random();
		if (normalizedPlotWeight < randomNumber){
			keepPlot = false;
		}
		return keepPlot;
	}


	public static void main(String[] args)throws Exception{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// input --> GUI-Erfassung, hier nur hilfsweise für self-contained example
		// input stuff
		//File rasterFile = new File("D:\\_HCU\\_Masterarbeit\\_TestData\\Globcover2009_V2.3_Global\\GLOBCOVER_L4_200901_200912_V2.3.tif"); // GUI-Erfassung
		File rasterFile = new File("D:\\_HCU\\_Masterarbeit\\_TestData\\weightraster\\fictive_weightraster.tif"); // GUI-Erfassung
		File shapeFile = new File("D:\\_HCU\\_Masterarbeit\\_TestData\\TestSmallPolygons\\TestSmallPolygons.shp"); // GUI-Erfassung
		String sampleColumn = "VEGZONE"; // GUI-Erfassung
		String clipStratum = "irregular"; // GUI-Erfassung
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		
		
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//clip(crop) raster:
// read raster
GridCoverage2D coverage = readGeoTiff(rasterFile);
// get Geometries to be clipped from the raster
ArrayList<Geometry> clipGeometries = getClipGeometries(coverage, shapeFile, sampleColumn, clipStratum);
// clip raster using clipGeometries
GridCoverage2D clippedCoverage = getClippedCoverage(clipGeometries, coverage);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		


		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//create Plot -> hier nur hilfsweise für self-contained example, kommt aus Sampling Tool
		FileDataStore dataStore = FileDataStoreFinder.getDataStore(shapeFile);
		SimpleFeatureSource featureSource = dataStore.getFeatureSource(); // wird benötigt, um an einzelne Features ranzukommen (Feature = Zeile in SHP Attribute Table)
		SimpleFeatureCollection collection = featureSource.getFeatures();

		
		CoordinateReferenceSystem crsRaster = coverage.getCoordinateReferenceSystem();

		GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
		Coordinate coord = new Coordinate( 795989.242986, 609624.337442 );
		Point point = geometryFactory.createPoint( coord );
		Plot plot = new Plot(point, "",  crsRaster); // wie komme ich an WGS84 CRS? --> vl über SHP
////////////////////////////////////////////////////////////////////////
		
		
		
		


		

		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this is the central code line!!!
		
		//test mit floating point und int rasters
		
		
		double pixelValue = 0;
		// find out raster pixel value number type before calling getValueAtPosition() method
		final int dataType = coverage.getRenderedImage().getSampleModel().getDataType();
		switch (dataType) {
		case DataBuffer.TYPE_BYTE:   byte[] pixelValueByte = getValueAtPosition(coverage, plot, (byte  []) null);
		pixelValue = (double)pixelValueByte[0];
		break;
		case DataBuffer.TYPE_DOUBLE: double[] pixelValueDouble = getValueAtPosition(coverage, plot, (double  []) null);
		pixelValue = (double)pixelValueDouble[0];
		break;
		default: throw new CannotEvaluateException("Raster with unknown pixel value number format");
		}

//		//Vorlage: 
//		/**
//	     * Returns the value vector for a given location (world coordinates).
//	     * A value for each sample dimension is included in the vector.
//	     */
//	    public Object evaluate(final DirectPosition point) throws CannotEvaluateException {
//	        final int dataType = image.getSampleModel().getDataType();
//	        switch (dataType) {
//	            case DataBuffer.TYPE_BYTE:   return evaluate(point, (byte  []) null);
//	            case DataBuffer.TYPE_SHORT:  // Fall through
//	            case DataBuffer.TYPE_USHORT: // Fall through
//	            case DataBuffer.TYPE_INT:    return evaluate(point, (int   []) null);
//	            case DataBuffer.TYPE_FLOAT:  return evaluate(point, (float []) null);
//	            case DataBuffer.TYPE_DOUBLE: return evaluate(point, (double[]) null);
//	            default: throw new CannotEvaluateException();
//	        }
//	    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
	

		
		
		
		
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double noDataValue = getNoDataValue(coverage, 0); // second param is the band index (weight raster should only have 1 band, so index is 0)
		double maxValue = getCoverageMaxValue(clippedCoverage, 0); // second param is the band index (weight raster should only have 1 band, so index is 0)
		double sum = getCoverageSum(clippedCoverage, 0); // only for pdf sampling report (?) // second param is the band index
		// normalize pixel value using maxValue
		double normalizedPixelValue =  pixelValue / maxValue; // make pixelValue be in the range 0-1
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		
		
		//writeGeoTiff("D:\\_HCU\\_Masterarbeit\\_TestData\\TestClip\\TestClipIrregular.tif", clippedCoverage);
		
		// rejection test
		boolean keepPlot = rejectionTesting(normalizedPixelValue);
		
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		System.out.println("input raster: " + rasterFile.toString());
		System.out.println("shapeFile: " + shapeFile.toString());
		System.out.println("Stratum: " + clipStratum.toString());
		System.out.println("Plot X: " + point.getX());
		System.out.println("Plot Y: " + point.getY());
		if (pixelValue == noDataValue){
			System.out.println("Error: pixelValue == noDataValue! ");
		}
		System.out.println("Raster pixel Value at Plot location: " + pixelValue);
		System.out.println("noDataValue: " + noDataValue);
		System.out.println("maxValue for stratum area: " + maxValue);
		System.out.println("sum of all pixel values for stratum area: " + sum);	
		System.out.println("normalizedPixelValue: " + normalizedPixelValue);
		System.out.println("rejection test: " + keepPlot);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		
	}
}
