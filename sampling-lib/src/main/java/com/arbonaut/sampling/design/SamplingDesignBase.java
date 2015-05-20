
package com.arbonaut.sampling.design;

import com.arbonaut.sampling.BadCRSException;
import com.arbonaut.sampling.io.SampleWriter;
import com.arbonaut.sampling.design.plot.PlotGeometry;
import com.arbonaut.sampling.design.cluster.PlotCluster;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Lineal;
import com.vividsolutions.jts.geom.Polygonal;
import com.vividsolutions.jts.operation.valid.IsValidOp; 
import com.vividsolutions.jts.operation.valid.TopologyValidationError; 

import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.crs.ProjectedCRS;
import org.geotools.referencing.CRS;


/**
 * Abstract base for sampling design implementations. 
 *  
 * At the moment does not do so much, namely validates the geometry of the 
 * census and its CRS, which must be a valid and projected reference system.        
 */
public abstract class SamplingDesignBase implements SamplingDesign {
	
    private PlotGeometry plotGeometry;
    private PlotCluster plotClustering;    
    
    public SamplingDesignBase(PlotGeometry plotGeometry, PlotCluster plotClustering) {
        this.plotGeometry = plotGeometry;
        this.plotClustering = plotClustering;
    }
    
    public PlotGeometry getPlotGeometry() { return this.plotGeometry; }
    public PlotCluster getPlotClustering() { return this.plotClustering; }
    
    
    
    public void generate(Geometry censusGeometry, 
                         CoordinateReferenceSystem censusCRS,
                         int sampleSize,
                         SampleWriter writer) throws Exception
    {
        generate(censusGeometry, censusCRS, sampleSize, null, writer);    
    }
    
    
    public void generate(Geometry censusGeometry, 
                         CoordinateReferenceSystem censusCRS,
                         int sampleSize,
                         String stratumName,
                         SampleWriter writer) throws Exception
    {
        boolean geomOk = 
            (censusGeometry instanceof Polygonal && canSampleArea()) ||
            (censusGeometry instanceof Lineal && canSampleAlongCurve());
    
        if (!geomOk)             
            throw new Exception("Census geometry not supported by this sampling design");
            
        // Check the CRS
            
        if (!(censusCRS instanceof ProjectedCRS))
            throw new BadCRSException("Sampling design requires projected coordinate system");             
            
        // Validate topology of the census geometry since incorrect geometries 
        // may cause problems later             
            
        //IsValidOp validOp = new IsValidOp(censusGeometry);
        //TopologyValidationError err = validOp.getValidationError();
        
        //if (err != null)
        //    throw new Exception("Census geometry is invalid: " + err.toString());            
            
        // All ok 
        doGenerate(censusGeometry, censusCRS, sampleSize, stratumName, writer);                              
    }                     
  
  
    /**
     * Subclasses can override this to specify whether the sampling design can 
     * be applied to a draw sample along a curve. Default returns false.     
     */         
    protected boolean canSampleAlongCurve() { return false; }
    
    /**
     * Subclasses can override this to specify whether the sampling design can 
     * be applied to a draw sample from within an area. Default returns true.     
     */         
    protected boolean canSampleArea() { return true; }
   

    /**
     * To be implemented by the sampling design to do the actual work.
     */         
    protected abstract void doGenerate(Geometry censusGeometry,
                                       CoordinateReferenceSystem censusCRS,
                                       int sampleSize, String stratumName, 
                                       SampleWriter writer) throws Exception;

 	
}
