
package com.arbonaut.sampling.design;

import com.arbonaut.sampling.io.SampleWriter;
import com.arbonaut.sampling.design.plot.PlotGeometry;
import com.arbonaut.sampling.design.cluster.PlotCluster;
                                                       
import com.vividsolutions.jts.geom.Geometry;
import org.opengis.referencing.crs.CoordinateReferenceSystem;


/**
 * Sampling design interface.
 * 
 * Sampling design encapsulates the data and logic needed to draw a spatial 
 * sample. The typical usage is such that once set up, the generate() method
 * is called one or more times (once per stratum if strata are in use). Once 
 * finished, the getSampleInfo() method can be called to collect information 
 * needed for inferring statistical properties of the sample.     
 */
public interface SamplingDesign {
	
    /**
     * Returns the plot geometry.
     */         
    public PlotGeometry getPlotGeometry();


    /**
     * Returns the plot cluster shape or null if clustering is not used.
     */         
    public PlotCluster getPlotClustering();
    
    
    /**
     * Generates new sample according to this sampling design. 
     * 
     * The sample is drawn from withing the supplied census, which must be a  
     * polygon or multipolygon geometry. The CRS must be a projected coordinate 
     * system, such as UTM. Sampling is always carried out in the coordinate 
     * system of the census; to output plots in different spatial reference, 
     * use the ReprojectingSampleWriter class. 
     *
     * @param censusGeometry  Geometry defining the census to be sampled
     * @param censusCRS       Spatial reference of the census geometry       
     * @param sampleSize      Number of plots to place       
     * @param writer          SampleWriter instance that stores the generated 
     *                        plots
     *                        
     * This overload assumes no stratum.          
     */         
    public void generate(Geometry censusGeometry,
                         CoordinateReferenceSystem censusCRS,
                         int sampleSize, 
                         SampleWriter writer) throws Exception;
                         
                         
    /**
     * Generates new sample according to this sampling design. 
     * 
     * The sample is drawn from withing the supplied census, which must be a  
     * polygon or multipolygon geometry. The CRS must be a projected coordinate 
     * system, such as UTM. Sampling is always carried out in the coordinate 
     * system of the census; to output plots in different spatial reference, 
     * use the ReprojectingSampleWriter class. 
     *
     * @param censusGeometry  Geometry defining the census to be sampled
     * @param censusCRS       Spatial reference of the census geometry       
     * @param sampleSize      Number of plots to place       
     * @param stratumName     Name of the stratum we're sampling or null if 
     *                        there is no stratum
     * @param writer          SampleWriter instance that stores the generated 
     *                        plots
     */         
    public void generate(Geometry censusGeometry,
                         CoordinateReferenceSystem censusCRS,
                         int sampleSize, 
                         String stratumName,
                         SampleWriter writer) throws Exception;                         
  
 	
    /**
     * Returns information about the sample generated under this design. 
     */              
    public SampleInfo getSampleInfo(); 
    
}
