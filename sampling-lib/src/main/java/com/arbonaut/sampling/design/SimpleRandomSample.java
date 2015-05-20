
package com.arbonaut.sampling.design;

import com.arbonaut.sampling.SamplePlot;
import com.arbonaut.sampling.io.SampleWriter;
import com.arbonaut.sampling.design.plot.PlotGeometry;
import com.arbonaut.sampling.design.cluster.PlotCluster;
import com.arbonaut.sampling.util.PlotInPolygonChecker;

import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Polygonal;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Location;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.crs.ProjectedCRS;
import org.geotools.geometry.jts.JTSFactoryFinder;


/**
 * Implements simple random sampling strategy.
 * 
 */
public class SimpleRandomSample extends SamplingDesignBase {
	
    int maxNumAttempts;
    
    
    /**
     * Constructor.
     * 
     * @param pgeom   Plot geometry factory defining the shape of each plot     
     * @param cluster Plot cluster instance or null not to generate clustered  
     *                sample                   
     */         
    public SimpleRandomSample(PlotGeometry plotGeometry, PlotCluster plotClustering) 
    {
        super(plotGeometry, plotClustering);
        this.maxNumAttempts = 100;  // Try 100x more than sample size, then quit
    }
    
    
    /**
     * Constructor.
     * 
     * @param pgeom   Plot geometry factory defining the shape of each plot     
     * @param cluster Plot cluster instance or null not to generate clustered  
     *                sample                   
     * @param attempts Maximum number of attempts for placing the desired number 
     *                 of plots before chickening out, relative to the sample 
     *                 size. This is still multiplied by sampleSize to get the 
     *                 actual limit. 
     */         
    public SimpleRandomSample(PlotGeometry plotGeometry, PlotCluster plotClustering,
                              int maxNumAttempts) 
    {
        super(plotGeometry, plotClustering);
        this.maxNumAttempts = maxNumAttempts;
    }
    
    
    @Override
    protected void doGenerate(Geometry censusGeometry, 
                              CoordinateReferenceSystem censusCRS,
                              int sampleSize, String stratumName,
                              SampleWriter writer) throws Exception
    {
        Envelope envelope = censusGeometry.getEnvelopeInternal(); 
		double minX = envelope.getMinX();
		double maxX = envelope.getMaxX();
		double minY = envelope.getMinY();
		double maxY = envelope.getMaxY();    

        PlotCluster cluster = this.getPlotClustering();

        PlotInPolygonChecker plotChecker = new PlotInPolygonChecker(
            (Polygonal)censusGeometry, this.getPlotGeometry());
            
        Coordinate c = new Coordinate();   
        Coordinate c2 = new Coordinate();
        
        int numPlots = 0;
        int clusterNr = 0;
        int attemptsRemaining = this.maxNumAttempts * sampleSize;
        
        while (numPlots < sampleSize && attemptsRemaining > 0)
        {
            --attemptsRemaining;
            
            // generate random X
			c.x = (Math.random() * (maxX - minX)) + minX;
			// generate random Y
			c.y = (Math.random() * (maxY - minY)) + minY;
            
            
            // Clustered plot
            if (cluster != null)
            {
                CoordinateSequence subplots = cluster.create(c);
                int clusterPlotNr = 0;
                for (int i = 0; i != subplots.size(); i++)
                {
                    subplots.getCoordinate(i, c2);
                    boolean added = sampleSinglePoint(
                        c2, censusCRS, clusterPlotNr, clusterNr, stratumName, 
                        writer, plotChecker);
                         
                    if (added) { ++numPlots; ++clusterPlotNr; }                    
                }
                
                // Placed at least 1 plot in te cluster -> increase cluster number
                if (clusterPlotNr > 0) ++clusterNr;
            }
            // Single plot
            else
            {
                boolean added = sampleSinglePoint(
                    c, censusCRS, numPlots, SamplePlot.NO_CLUSTER, stratumName, 
                    writer, plotChecker); 
                if (added) ++numPlots;                                  
            }
		}
      
    }
    
    
    public SampleInfo getSampleInfo()
    {
        AreaSampleInfo info = new AreaSampleInfo();
        
        // Todo: fill this
        
        return info;
    }    
    
    
    
    private boolean sampleSinglePoint(Coordinate coord, CoordinateReferenceSystem crs,
                                      int plotNr, int clusterNr,
                                      String stratumName, SampleWriter writer,
                                      PlotInPolygonChecker plotChecker) throws Exception
    {
        // Quickly check if this point will do just based on the centre so
        // we could possibly skip creation of the geometry completely 
        int location = plotChecker.centreTest(coord);
        
        // Absolutely outside
        if (location == Location.EXTERIOR) return false;
        
        // Need geometry test
        Geometry plotGeom = null;
        if (location == Location.BOUNDARY)
        {
            plotGeom = this.getPlotGeometry().create(coord);
            // Outside
            if (!plotChecker.contains(plotGeom)) return false;
        }   
        
        // Got so far so we're inside                
        if (plotGeom == null)
            plotGeom = this.getPlotGeometry().create(coord);
        
        SamplePlot plot = 
            new SamplePlot(plotGeom, crs, stratumName, plotNr, clusterNr);  
        writer.write(plot);
        
        return true;
    }
    

	
}
