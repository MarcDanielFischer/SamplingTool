
package com.arbonaut.sampling.io.filter;

import java.io.IOException;
import com.arbonaut.sampling.SamplePlot;
import com.arbonaut.sampling.BadCRSException;

import com.vividsolutions.jts.geom.Geometry;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.geotools.referencing.CRS;
import org.geotools.geometry.jts.JTS;


/**
 * Sample writer that reprojects plot geometry to different CRS and then 
 * delegates the work to the underlying writer.
 * 
 * ReprojectingSampleWriter is used for writing the sample plots in a different
 * CRS than in which they were sampled.   
 */
public class ReprojectingSampleFilter implements SampleFilter {

    CoordinateReferenceSystem targetCRS;
    
    /**
     * ReprojectingSampleFilter constructor.
     * 
     * @param targetEPSGCode  EPSG code of the target coordinate system     
     * @throws NoSuchAuthorityCodeException if epsg code is invalid                    
     */         
    public ReprojectingSampleFilter(int targetEPSGCode) 
        throws NoSuchAuthorityCodeException, FactoryException
    {
        this.targetCRS = CRS.decode(String.format("EPSG:%d", targetEPSGCode));
    }
    
    
    /**
     * ReprojectingSampleFilter constructor.
     * @param targetCRS  Target coordinate system               
     */         
    public ReprojectingSampleFilter(CoordinateReferenceSystem targetCRS)
    {
        this.targetCRS = targetCRS;
    }
    
    
    public SamplePlot filter(SamplePlot plot) throws Exception 
    {
        if (plot.getCRS() == null) 
            throw new BadCRSException("Sample plot has no CRS"); 
            
        MathTransform xform = CRS.findMathTransform(
            plot.getCRS(), this.targetCRS);
    
        Geometry projected = JTS.transform(plot.getGeometry(), xform);

        return new SamplePlot(plot, projected, this.targetCRS);
    }	

}
