
package com.arbonaut.sampling.io.filter;

import java.io.IOException;
import com.arbonaut.sampling.SamplePlot;


/**
 * Applies transformation to a sample plot.  
 */
public interface SampleFilter {

    /**
     * Applies transformation to a sample plot and returns the transformed 
     * result.
     */         
    public SamplePlot filter(SamplePlot plot) throws Exception;

}
